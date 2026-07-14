library(tidyverse)
library(lubridate)
library(aws.s3)
library(yaml)
library(httr)

source("baseline_models/R/randomWalkDailyFunction.R")
source("baseline_models/models/coastal_climatology.R")
source("baseline_models/models/coastal_random_walk.R")

Sys.setenv("AWS_DEFAULT_REGION" = "")

config <- yaml::read_yaml("challenge_configuration.yaml")
null_start_date <- as_date(config$target_groups$Coastal$null_start_date)
base_url <- gsub("https://", "", config$endpoint)

# Read corrected targets (chlora_cci_corrected for both sites)
corrected_url <- paste0(config$endpoint, "/", config$s3_bucket_read, "/",
                        config$target_groups$Coastal$targets_corrected_filepath)
targets_all <- readr::read_csv(corrected_url, guess_max = 10000) %>%
  mutate(datetime = as_date(datetime), site_id = as.character(site_id))

# Read raw targets for chlora_mrwa (site 2 in-situ buoy)
raw_url <- paste0(config$endpoint, "/", config$s3_bucket_read, "/",
                  config$target_groups$Coastal$targets_filepath)
raw_targets_mrwa <- readr::read_csv(raw_url, guess_max = 10000) %>%
  mutate(datetime = as_date(datetime), site_id = as.character(site_id)) %>%
  filter(site_id == "2", variable == "chlora_mrwa")

# Site 1: chlora_cci_corrected only; forecasts from site 1 null start date to yesterday
targets_site1 <- targets_all %>% filter(site_id == "1", variable == "chlora_cci_corrected")
site1_dates   <- seq(null_start_date, as_date(config$target_groups$Coastal$site_2_forecast_end), by = "day")

# Site 2: both chlora_cci_corrected and chlora_mrwa; forecasts from site 2 null start date to forecast end
SITE2_START <- as_date(config$target_groups$Coastal$site_2_null_start_date)
SITE2_END   <- as_date(config$target_groups$Coastal$site_2_forecast_end)
targets_site2 <- bind_rows(
  targets_all %>% filter(site_id == "2", variable == "chlora_cci_corrected"),
  raw_targets_mrwa
)
site2_dates <- seq(SITE2_START, SITE2_END, by = "day")

# Get reference dates already uploaded for a given model
get_existing_dates <- function(model_name) {
  get_dates_from_prefix <- function(prefix) {
    tryCatch({
      files <- aws.s3::get_bucket_df(
        bucket = config$s3_bucket_write,
        prefix = paste0(prefix, "/coastal-"),
        base_url = base_url,
        use_https = TRUE,
        region = "",
        max = Inf
      )
      if (nrow(files) == 0) return(as_date(character(0)))
      files %>%
        pull(Key) %>%
        .[stringr::str_detect(., model_name)] %>%
        stringr::str_extract("\\d{4}-\\d{2}-\\d{2}") %>%
        as_date() %>%
        na.omit()
    }, error = function(e) as_date(character(0)))
  }
  unique(c(
    get_dates_from_prefix(config$forecasts_bucket),
    get_dates_from_prefix("challenges/project_id=bu4cast/raw-submissions")
  ))
}

run_site_dates <- function(run_fn, dates, targets, label) {
  existing <- get_existing_dates(label)
  missing  <- dates[!dates %in% existing]
  message(length(missing), " ", label, " dates to run")
  for (ref_date in as.list(missing)) {
    run_fn(as_date(ref_date), config, targets)
  }
}

# Climatology
run_site_dates(run_coastal_climatology, site1_dates, targets_site1, "climatology")
run_site_dates(run_coastal_climatology, site2_dates, targets_site2, "climatology")
httr::GET(config$target_groups$Coastal$health_checks$climatology_null)

# Random walk
run_site_dates(run_coastal_random_walk, site1_dates, targets_site1, "randomWalk")
run_site_dates(run_coastal_random_walk, site2_dates, targets_site2, "randomWalk")
httr::GET(config$target_groups$Coastal$health_checks$random_walk_null)

