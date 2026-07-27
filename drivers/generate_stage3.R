library(minioclient)
library(arrow)
library(dplyr)
library(yaml)
source("drivers/to_hourly.R")

minioclient::install_mc()
mc_alias_set("osn", "uri.osn.mghpcc.org", Sys.getenv("OSN_KEY"), Sys.getenv("OSN_SECRET"))

config <- yaml::read_yaml("challenge_configuration.yaml")

s3 <- arrow::s3_bucket(
  config$s3_bucket_read,
  endpoint_override = config$endpoint,
  access_key = Sys.getenv("OSN_KEY"),
  secret_key = Sys.getenv("OSN_SECRET"),
  scheme = "https"
)

metadata_path <- gsub(paste0("^", config$s3_bucket_read, "/"), "", config$target_metadata_bucket)
drivers_path  <- gsub(paste0("^", config$s3_bucket_read, "/"), "", config$drivers_bucket)

site_coords <- arrow::read_csv_arrow(
  s3$path(config$field_sites_path)
) %>%
  as.data.frame() %>%
  transmute(
    site_id   = as.character(field_site_id),
    latitude  = as.numeric(latitude),
    longitude = as.numeric(longitude)
  )

site_list <- site_coords$site_id
message("Sites loaded: ", length(site_list))

mc_alias_set("osn", "uri.osn.mghpcc.org", Sys.getenv("OSN_KEY"), Sys.getenv("OSN_SECRET"))
mc_mirror("osn/bu4cast-ci/read/challenges/project_id=bu4cast/drivers/pseudo", "pseudo")

future::plan("future::multisession", workers = 4)
furrr::future_walk(site_list, function(curr_site_id) {
  source("drivers/to_hourly.R")
  library(arrow)
  library(dplyr)
  library(yaml)

  config       <- yaml::read_yaml("challenge_configuration.yaml")
  drivers_path <- gsub(paste0("^", config$s3_bucket_read, "/"), "", config$drivers_bucket)

  s3_w <- arrow::s3_bucket(
    config$s3_bucket_read,
    endpoint_override = config$endpoint,
    access_key = Sys.getenv("OSN_KEY"),
    secret_key = Sys.getenv("OSN_SECRET"),
    scheme = "https"
  )

  # skip if stage3 already exists for this site
  existing <- tryCatch(
    arrow::open_dataset(s3_w$path(paste0(drivers_path, "/stage3"))) %>%
      dplyr::filter(site_id == curr_site_id) %>%
      dplyr::collect(),
    error = function(e) data.frame()
  )
  if (nrow(existing) > 0) {
    message("Skipping ", curr_site_id, " - stage3 already exists")
    return(NULL)
  }

  df <- arrow::open_dataset("pseudo") |>
    dplyr::filter(variable %in% c("PRES", "TMP", "RH", "UGRD", "VGRD", "APCP", "DSWRF", "DLWRF")) |>
    dplyr::filter(site_id == curr_site_id) |>
    dplyr::collect() |>
    dplyr::mutate(date         = lubridate::as_date(reference_datetime),
                  new_datetime = date + lubridate::hours(as.numeric(horizon)) + lubridate::hours(as.numeric(cycle)),
                  datetime     = ifelse(datetime != new_datetime, new_datetime, datetime),
                  datetime     = lubridate::as_datetime(datetime)) |>
    dplyr::select(-date, -new_datetime)

    print(curr_site_id)
    
      hourly_df <- df |>
        to_hourly(use_solar_geom = TRUE, psuedo = TRUE) |>
        dplyr::mutate(ensemble = as.numeric(stringr::str_sub(ensemble, start = 4, end = 5))) |>
        dplyr::rename(parameter = ensemble)
    
      arrow::write_dataset(hourly_df,
                           path         = s3_w$path(paste0(drivers_path, "/stage3")),
                           partitioning = "site_id")
    },
    .options = furrr::furrr_options(globals = "site_coords"))
