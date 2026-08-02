## A01 Buoy Covariates for 2025-present
## Appends preliminary (uncorrected, not MWRA-validated) buoy sensor covariates
## for site 2 (A01) to a01_drivers.csv from the same realtime NetCDF
## data used for the chlora_pre target. Historical (MWRA-corrected coverage ends 2024-12-14)

library(aws.s3)
library(ncdf4)
library(dplyr)
library(lubridate)
library(readr)
library(yaml)

config <- yaml::read_yaml("challenge_configuration.yaml")

base_url <- gsub("https://", "", config$endpoint)
Sys.setenv(AWS_ACCESS_KEY_ID     = Sys.getenv("OSN_KEY"),
           AWS_SECRET_ACCESS_KEY = Sys.getenv("OSN_SECRET"),
           AWS_DEFAULT_REGION    = "")

drivers_key <- config$target_groups$Coastal$a01_drivers_filepath
origin      <- as.POSIXct("1858-11-17 00:00:00", tz = "UTC")

doppler_depths <- c(10, 14, 18, 22, 26, 30, 34, 38, 42, 46, 50, 54, 58)

## Existing drivers file 

old_data <- tryCatch({
  tmp <- tempfile(fileext = ".csv")
  aws.s3::save_object(object = drivers_key, bucket = config$s3_bucket_read,
                      file = tmp, base_url = base_url, use_https = TRUE, region = "")
  readr::read_csv(tmp, show_col_types = FALSE)
}, error = function(e) {
  message("No existing a01_drivers.csv found: ", conditionMessage(e))
  NULL
})

last_dt <- if (!is.null(old_data) && nrow(old_data) > 0) {
  max(old_data$datetime_hour, na.rm = TRUE)
} else {
  as.POSIXct("2005-10-22", tz = "UTC")
}
message("Last existing datetime_hour: ", last_dt)

## Download realtime NetCDF for one sensor group

fetch_nc <- function(url) {
  tmp <- tempfile(fileext = ".nc")
  ok <- tryCatch({
    utils::download.file(url, destfile = tmp, mode = "wb", quiet = TRUE)
    TRUE
  }, error = function(e) { message("download failed: ", conditionMessage(e)); FALSE })
  if (!ok || !file.exists(tmp) || file.info(tmp)$size == 0) return(NULL)
  nc <- tryCatch(ncdf4::nc_open(tmp), error = function(e) NULL)
  if (is.null(nc)) { unlink(tmp); return(NULL) }
  nc
}

to_hourly <- function(df) {
  df %>%
    dplyr::mutate(datetime_hour = lubridate::floor_date(datetime, "hour")) %>%
    dplyr::select(-datetime) %>%
    dplyr::group_by(datetime_hour) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
}

## Met: air_temperature (single depth, -3m)

message("Downloading A01 met (air_temperature)...")
met_hourly <- NULL
nc <- fetch_nc(config$target_groups$Coastal$a01_met_url)
if (!is.null(nc)) {
  t <- ncvar_get(nc, "time")
  datetime <- origin + t * 86400
  vals <- as.numeric(ncvar_get(nc, "air_temperature"))
  qc   <- as.integer(ncvar_get(nc, "air_temperature_qc"))
  vals[is.na(qc) | qc != 0] <- NA
  df <- data.frame(datetime = datetime, x = vals)
  names(df)[2] <- "air_temperature[C](-3m)"
  met_hourly <- to_hourly(df)
  nc_close(nc)
}

## SBE37 @ 1m and 20m: temperature, salinity, conductivity

read_sbe37 <- function(url, depth_label) {
  nc <- fetch_nc(url)
  if (is.null(nc)) return(NULL)
  t <- ncvar_get(nc, "time")
  datetime <- origin + t * 86400

  extract <- function(var) {
    vals <- as.numeric(ncvar_get(nc, var))
    qc   <- as.integer(ncvar_get(nc, paste0(var, "_qc")))
    vals[is.na(qc) | qc != 0] <- NA
    vals
  }

  df <- data.frame(
    datetime = datetime,
    temperature  = extract("temperature"),
    salinity     = extract("salinity"),
    conductivity = extract("conductivity")
  )
  names(df) <- c("datetime",
                 paste0("temperature[C](", depth_label, ")"),
                 paste0("salinity[psu](", depth_label, ")"),
                 paste0("conductivity[msiemens/cm](", depth_label, ")"))
  nc_close(nc)
  to_hourly(df)
}

message("Downloading A01 SBE37 @ 1m...")
sbe1m_hourly  <- read_sbe37(config$target_groups$Coastal$a01_sbe37_1m_url,  "1m")
message("Downloading A01 SBE37 @ 20m...")
sbe20m_hourly <- read_sbe37(config$target_groups$Coastal$a01_sbe37_20m_url, "20m")

## Doppler: current_u/current_v at 10-58m 

message("Downloading A01 doppler currents...")
doppler_hourly <- NULL
nc <- fetch_nc(config$target_groups$Coastal$a01_doppler_url)
if (!is.null(nc)) {
  t <- ncvar_get(nc, "time")
  datetime <- origin + t * 86400
  depths <- ncvar_get(nc, "depth")
  keep_idx <- which(depths %in% doppler_depths)

  cu <- ncvar_get(nc, "current_u")
  cv <- ncvar_get(nc, "current_v")
  qu <- ncvar_get(nc, "current_u_qc")
  qv <- ncvar_get(nc, "current_v_qc")

  df <- data.frame(datetime = datetime)
  for (di in keep_idx) {
    d <- depths[di]
    u <- cu[di, ]; u[is.na(qu[di, ]) | qu[di, ] != 0] <- NA
    v <- cv[di, ]; v[is.na(qv[di, ]) | qv[di, ] != 0] <- NA
    df[[paste0("current_u[cm/s](", d, "m)")]] <- u
    df[[paste0("current_v[cm/s](", d, "m)")]] <- v
  }
  doppler_hourly <- to_hourly(df)
  nc_close(nc)
}

## Combine all sensor groups

pieces <- Filter(Negate(is.null), list(met_hourly, sbe1m_hourly, sbe20m_hourly, doppler_hourly))

if (length(pieces) == 0) {
  message("No sensor data retrieved this run; nothing to append.")
} else {
  combined <- Reduce(function(a, b) dplyr::full_join(a, b, by = "datetime_hour"), pieces) %>%
    dplyr::arrange(datetime_hour) %>%
    dplyr::filter(datetime_hour > last_dt)

  message("New rows to append: ", nrow(combined))

  if (nrow(combined) == 0) {
    message("Already up-to-date; no new rows to write.")
  } else {
    new_data <- if (!is.null(old_data) && nrow(old_data) > 0) {
      dplyr::bind_rows(old_data, combined)
    } else {
      combined
    }
    new_data <- new_data %>% dplyr::arrange(datetime_hour)

    tmp_out <- tempfile(fileext = ".csv")
    readr::write_csv(new_data, tmp_out)
    aws.s3::put_object(file = tmp_out, object = drivers_key, bucket = config$s3_bucket_read,
                       base_url = base_url, use_https = TRUE, region = "")
    message("Uploaded updated a01_drivers.csv (", nrow(new_data), " total rows) to s3://",
            config$s3_bucket_read, "/", drivers_key)
  }
}

message("A01 covariates script complete!")
