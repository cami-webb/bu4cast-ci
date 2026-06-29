#devtools::install_version("duckdb", "1.2.2")
#remotes::install_github('cboettig/duckdbfs', upgrade = 'never')

library(dplyr)
library(duckdbfs)
library(progress)
library(bench)
library(yaml)
library(stringr)
library(minioclient)
library(DBI)

con <- duckdbfs::cached_connection(tempfile())
DBI::dbExecute(con, "SET THREADS=64;")

setup_s3 <- function(con) {
  DBI::dbExecute(con, paste0("SET s3_access_key_id='",     Sys.getenv("OSN_KEY"),    "';"))
  DBI::dbExecute(con, paste0("SET s3_secret_access_key='", Sys.getenv("OSN_SECRET"), "';"))
  DBI::dbExecute(con, "SET s3_endpoint='minio-s3.apps.shift.nerc.mghpcc.org';")
  DBI::dbExecute(con, "SET s3_url_style='path';")
  DBI::dbExecute(con, "SET s3_use_ssl=true;")
  invisible(con)
}
setup_s3(con)

install_mc()
mc_alias_set("minio", "minio-s3.apps.shift.nerc.mghpcc.org", Sys.getenv("OSN_KEY"), Sys.getenv("OSN_SECRET"))

remove_dir <- function(path) {
  tryCatch(
    {
      minioclient::mc_rm(path, recursive = TRUE)
      message('directory successfully removed...')
    },
    error = function(cond) {
      message("The removal directory could not be found...")
      message(conditionMessage(cond))
      NA
    },
    warning = function(cond) {
      message('Deleting the directory caused a warning...')
      message(conditionMessage(cond))
      NULL
    },
    finally = {
      message("Finished the delete portion...")
    }
  )
}

remove_dir("minio/bu4cast-ci-write/challenges/project_id=bu4cast/tmp/score_me")
remove_dir("minio/bu4cast-ci-write/challenges/project_id=bu4cast/tmp/forecasts")
remove_dir("minio/bu4cast-ci-write/challenges/project_id=bu4cast/tmp/targets")
remove_dir("minio/bu4cast-ci-write/challenges/project_id=bu4cast/tmp/scores")

config <- read_yaml("challenge_configuration.yaml")

forecast_bundled_parquet_bucket <- paste0(config$s3_bucket_read, "/", config$forecasts_bucket, "/bundled-parquet/")
scores_bundled_parquet_bucket   <- paste0(config$scores_bucket, "/bundled-parquet/")

project        <- config$project_id
cut_off_date   <- Sys.Date() - lubridate::dmonths(6)
rescore        <- FALSE
obs_key_cols   <- c("project_id", "site_id", "datetime", "duration", "variable")
score_key_cols <- c(obs_key_cols, "model_id", "family", "reference_datetime")

# Build target file list, skipping any that don't exist in S3
target_files <- NULL
for (i in 1:length(config$target_groups)) {
  path <- paste0("s3://", config$s3_bucket_read, "/", config$target_groups[[i]]$targets_filepath)
  exists <- tryCatch({
    duckdbfs::open_dataset(path, format = "csv") |> dplyr::count() |> dplyr::collect()
    TRUE
  }, error = function(e) {
    message("Skipping ", basename(path), ": not found in S3")
    FALSE
  })
  if (exists) target_files <- c(target_files, path)
}

targets <-
  open_dataset(target_files,
               recursive = FALSE,
               format = "csv",
               parser_options = list(nullstr = "NA")
               ) |>
  filter(project_id == {project},
         datetime > {cut_off_date},
         !is.na(observation))

last_observed_date <- targets |> select(datetime) |> distinct() |>
  filter(datetime == max(datetime)) |> pull(datetime)

forecasts <-
  open_dataset(paste0("s3://", forecast_bundled_parquet_bucket)) |>
  filter(project_id == {project},
         datetime > {cut_off_date},
         datetime <= {last_observed_date},
         !is.na(model_id),
         !is.na(parameter),
         !is.na(prediction)
  ) |>
  mutate(family = ifelse(family == 'ensemble', "sample", family)) |>
  mutate(horizon = date_diff('day', as.POSIXct(reference_datetime), as.POSIXct(datetime))) |>
  filter(!(duration == "P1D" & horizon > 35))

scores <-
  open_dataset(paste0("s3://", scores_bundled_parquet_bucket)) |>
  filter(project_id == {project},
         datetime > {cut_off_date},
         !is.na(observation))

tol <- 1e-2
if (rescore) {
  print("rescoring changed observations")
  scores <- scores |>
    inner_join(targets, by = obs_key_cols) |>
    filter(abs(observation.x - observation.y) / observation.x < {tol})
}

print("Caching forecasts, scores, targets...")

bench::bench_time({
  forecasts |> group_by(variable) |> write_dataset("s3://bu4cast-ci-write/challenges/project_id=bu4cast/tmp/forecasts")
})
bench::bench_time({
  scores |> group_by(variable) |> write_dataset("s3://bu4cast-ci-write/challenges/project_id=bu4cast/tmp/scores")
})
bench::bench_time({
  targets |> group_by(variable) |> write_dataset("s3://bu4cast-ci-write/challenges/project_id=bu4cast/tmp/targets")
})
bench::bench_time({
  forecasts <- open_dataset("s3://bu4cast-ci-write/challenges/project_id=bu4cast/tmp/forecasts/**")
  scores    <- open_dataset("s3://bu4cast-ci-write/challenges/project_id=bu4cast/tmp/scores/**")
  targets   <- open_dataset("s3://bu4cast-ci-write/challenges/project_id=bu4cast/tmp/targets/**")
})

print("Compute who needs to be scored...")
bench::bench_time({
  forecasts |>
    anti_join(select(scores, all_of(score_key_cols))) |>
    inner_join(targets) |>
    group_by(variable) |>
    write_dataset("s3://bu4cast-ci-write/challenges/project_id=bu4cast/tmp/score_me")
})
