library(tidyverse)
library(lubridate)
library(aws.s3)
library(yaml)
library(httr)

source("baseline_models/R/randomWalkDailyFunction.R")
source("baseline_models/models/coastal_climatology.R")
source("baseline_models/models/coastal_random_walk.R")
source("targets/target_helper_functions.R") # read_s3_csv_weblink()

Sys.setenv("AWS_DEFAULT_REGION" = "")

config <- yaml::read_yaml("challenge_configuration.yaml")
null_start_date <- as_date(config$target_groups$Coastal$null_start_date)
base_url <- gsub("https://", "", config$endpoint)

# Two coverage windows: the 2021 exploratory backtest, then a gap, then daily
# coverage resuming from daily_start_date through yesterday (extends each run).
daily_start_date <- as_date(config$target_groups$Coastal$daily_start_date)
daily_end_date    <- Sys.Date() - 1
daily_dates <- if (daily_start_date <= daily_end_date) {
  seq(daily_start_date, daily_end_date, by = "day")
} else {
  as_date(character(0))
}

# Read corrected targets (chlora_cci_corrected for both sites)
corrected_url <- paste0(config$endpoint, "/", config$s3_bucket_read, "/",
                        config$target_groups$Coastal$targets_corrected_filepath)
targets_all <- read_s3_csv_weblink(corrected_url, config, guess_max = 10000) %>%
  mutate(datetime = as_date(datetime), site_id = as.character(site_id))

# Read raw targets for chlora_mrwa (site 2 in-situ buoy)
raw_url <- paste0(config$endpoint, "/", config$s3_bucket_read, "/",
                  config$target_groups$Coastal$targets_filepath)
raw_targets_mrwa <- read_s3_csv_weblink(raw_url, config, guess_max = 10000) %>%
  mutate(datetime = as_date(datetime), site_id = as.character(site_id)) %>%
  filter(site_id == "2", variable == "chlora_mrwa")

# Site 1: chlora_cci_corrected only; 2021 exploratory backtest + daily from daily_start_date
targets_site1 <- targets_all %>% filter(site_id == "1", variable == "chlora_cci_corrected")
site1_dates   <- c(
  seq(null_start_date, as_date(config$target_groups$Coastal$site_2_forecast_end), by = "day"),
  daily_dates
)

# Site 2: both chlora_cci_corrected and chlora_mrwa; same two-window coverage
SITE2_START <- as_date(config$target_groups$Coastal$site_2_null_start_date)
SITE2_END   <- as_date(config$target_groups$Coastal$site_2_forecast_end)
targets_site2 <- bind_rows(
  targets_all %>% filter(site_id == "2", variable == "chlora_cci_corrected"),
  raw_targets_mrwa
)
site2_dates <- c(seq(SITE2_START, SITE2_END, by = "day"), daily_dates)

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
    get_dates_from_prefix("write/challenges/project_id=bu4cast/raw-submissions")
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

targets_combined <- bind_rows(targets_site1, targets_site2)
all_dates        <- site1_dates  # site2_dates are all within site1_dates

# Climatology
run_site_dates(run_coastal_climatology, all_dates, targets_combined, "climatology")
httr::GET(config$target_groups$Coastal$health_checks$climatology_null)

# Random walk
run_site_dates(run_coastal_random_walk, all_dates, targets_combined, "randomWalk")
httr::GET(config$target_groups$Coastal$health_checks$random_walk_null)
