## Urban Targets file for the BU Forecasting Challenge
## Created: 10/20/2025

library(dplyr)
library(tidyr)
library(arrow)
library(httr)
library(jsonlite)
library(lubridate)
library(RCurl)
library(aws.s3)
source("targets/target_helper_functions.R")
print(paste0("Running urban_targets.R at ", Sys.time()))

# Step 0: Reload data for appending ---------------------------------------

# Get Configurations
config <- yaml::read_yaml("challenge_configuration.yaml")
challenge_name = config$target_groups$Urban$target_name
filename = config$target_groups$Urban$targets_filepath
read_bucket = config$s3_bucket_read
endpoint = config$endpoint

# Set Up Connection
s3_read <- arrow::s3_bucket(read_bucket,
                            endpoint_override = endpoint,
                            access_key = Sys.getenv("OSN_KEY"),
                            secret_key = Sys.getenv("OSN_SECRET"),
                            scheme = "https")

# Read in old data

# Write out file and challenge name for debugging
print(paste("challenge_name:", challenge_name))
print(paste("filename:", filename))

# Read in old data
urban_data_url = config$target_groups$Urban$targets_weblink
old_data = read_s3_csv_weblink(urban_data_url, config,
                          col_types = cols(project_id = col_character(),
                                           site_id = col_character(),
                                           datetime = col_character(),
                                           duration = col_character(),
                                           variable = col_character(),
                                           observation = col_double()))

old_data$datetime <- as.POSIXct(old_data$datetime, format = "%Y-%m-%d %H:%M", tz = "GMT") 

# Step 1: Download Data (last year and this year) -------------------------

# Set variables for date, counties, and parameters

# Dates
last_year = as.numeric(format(Sys.Date(), '%Y')) - 1
today = format(Sys.Date(), '%Y%m%d')

# Parameter codes
pollutant_codes = c('88101', '81102', '44201', '42602')
pollutant_names = c('PM2.5', 'PM10', 'O3', 'NO2')

# County Codes
county_codes = c('025', '009', '021')
county_names = c('Suffolk', 'Essex', 'Norfolk')

# Write out parameter codes for debugging
print(paste0('Pollutant Codes: ', paste(pollutant_codes, collapse = ', ')))
print(paste0('Pollutant Names: ', paste(pollutant_names, collapse = ', ')))

# Create an empty df
updated_data = data.frame()

# Loop through counties
for (i in seq_along(county_codes)){
  county_code = county_codes[[i]]
  county_name = county_names[[i]]
  
  print(paste0('Downloading ', county_name, ' data...'))
  
  # Clear temp_df so we can check if data was downloaded
  temp_df1 = data.frame()
  temp_df2 = data.frame()
  
  # Download last year's data
  url = paste0(
    'https://aqs.epa.gov/data/api/sampleData/byCounty?',
    'email=', Sys.getenv("EPA_EMAIL"),
    '&key=', Sys.getenv("EPA_KEY"),
    '&param=', paste(pollutant_codes, collapse = ','),
    '&bdate=', last_year, '0101',
    '&edate=', last_year, '1231',
    '&state=25',
    '&county=', county_code
  )
  
  flag = tryCatch({
    # Download/Convert the data
    response = GET(url)
    data = fromJSON(content(response, 'text'))
    temp_df1 = data$Data
    
    # Append and overwrite data in old_data
    updated_data = rbind(updated_data, temp_df1)

    print(paste0('Processing ', last_year,' data'))
    
    # Set false to not flag
    FALSE
  }, error = function(e) {
    message(paste0('Error for ', last_year, ': ', e$message))
    
    # Set True for error flag
    TRUE
  })
  
  # Download this year's data
  url = paste0(
    'https://aqs.epa.gov/data/api/sampleData/byCounty?',
    'email=', Sys.getenv("EPA_EMAIL"),
    '&key=', Sys.getenv("EPA_KEY"),
    '&param=', paste(pollutant_codes, collapse = ','),
    '&bdate=', format(Sys.Date(), '%Y'), '0101',
    '&edate=', today,
    '&state=25',
    '&county=', county_code
  )
  
  flag = tryCatch({
    # Download/Convert the data
    response = GET(url)
    data = fromJSON(content(response, 'text'))
    temp_df2 = data$Data
    
    # Append and overwrite data in old_data
    updated_data = rbind(updated_data, temp_df2)

    print(paste0('Processing ', format(Sys.Date(), '%Y'),' data'))
    
    # Set false to not flag
    FALSE
  }, error = function(e) {
    message(paste0('Error for ', format(Sys.Date(), '%Y'), ': ', e$message))
    
    # Set True for error flag
    TRUE
  })
  
  if (flag) next

  # Let sleep after each county to prevent high influx of requests
  Sys.sleep(10)
}

# Stop if no new data downloaded
if (nrow(updated_data) == 0) {
  message("No new data downloaded. Stopping script.")
  stop("Exiting script because updated_data is empty.")
}


# Step 2: Organize -----------------------------------------------------

copy_updated_data = updated_data

# Set date as datetime
copy_updated_data$datetime <- as.POSIXct(
  paste(copy_updated_data$date_gmt, copy_updated_data$time_gmt),
  format = "%Y-%m-%d %H:%M",
  tz = "GMT"
)

# Update duration colum to ISO 8601 format
copy_updated_data$sample_duration = gsub("1 HOUR", "PT1H", copy_updated_data$sample_duration)
copy_updated_data$sample_duration = gsub("24 HOUR", "P1D", copy_updated_data$sample_duration)

# Update variable name
copy_updated_data$parameter = gsub("PM2.5 - Local Conditions", 'PM2.5', copy_updated_data$parameter)
copy_updated_data$parameter = gsub("PM10 Total 0-10um STP", 'PM10', copy_updated_data$parameter)
copy_updated_data$parameter = gsub("Ozone", 'O3', copy_updated_data$parameter)
copy_updated_data$parameter = gsub("Nitrogen dioxide \\(NO2\\)", "NO2", copy_updated_data$parameter)
copy_updated_data <- copy_updated_data %>%
  mutate(parameter = case_when(
    parameter == "PM2.5" & sample_duration == "PT1H" ~ "PM2.5_P1H",
    parameter == "PM2.5" & sample_duration != "PT1H" ~ "PM2.5_P1D",
    parameter == "PM10"  & sample_duration == "PT1H" ~ "PM10_P1H",
    parameter == "PM10"  & sample_duration != "PT1H" ~ "PM10_P1D",
    parameter == "NO2"   & sample_duration == "PT1H" ~ "NO2_P1H",
    parameter == "NO2"   & sample_duration != "PT1H" ~ "NO2_P1D",
    TRUE ~ parameter
  ))

copy_updated_data$site_id = paste(copy_updated_data$state_code,
                                            copy_updated_data$county_code,
                                            copy_updated_data$site_number,
                                      sep = '-')

# Remove duplicates
copy_updated_data <- copy_updated_data %>%
  mutate(
    date_gmt = as.POSIXct(date_gmt, tz = "GMT"),
    date_of_last_change = as.POSIXct(date_of_last_change, tz = "GMT")
  ) %>%
  drop_na(sample_measurement) %>% # remove NA's
  # Most recently updated for datetime, variable, site_id, duration, and poc
  group_by(datetime, parameter, site_id, sample_duration, poc) %>%
  slice_max(date_of_last_change, n = 1, with_ties = FALSE) %>%  
  ungroup() %>%
  # Longer duration for datetime, variable, site_id, and duration
  group_by(site_id, parameter, poc) %>%
  mutate(
    sensor_duration = as.numeric(difftime(
      max(date_gmt, na.rm = TRUE),
      min(date_gmt, na.rm = TRUE),
      units = "days"
    ))
  ) %>%
  ungroup() %>%
  group_by(datetime, parameter, site_id, sample_duration) %>%
  slice_max(sensor_duration, n = 1, with_ties = FALSE) %>%
  ungroup()

# Select data
data = copy_updated_data[, c('site_id', 'datetime', 'sample_duration',
                       'parameter', 'sample_measurement')]

# Rename columns
colnames(data) = c('site_id', 'datetime', 'duration', 'variable', 'observation')

# Set project id
data$project_id = 'bu4cast'

# Reorder columns
data = data[, c('project_id', 'site_id', 'datetime', 'duration', 'variable',
                'observation')]

# Switch to dataframes
old_data <- as.data.frame(old_data)
data <- as.data.frame(data)

# Clean date type
old_data <- old_data %>%
  mutate(datetime_str = format(datetime, "%Y-%m-%d %H:%M:%S"))
data <- data %>%
  mutate(datetime_str = format(datetime, "%Y-%m-%d %H:%M:%S"))

# Round observations
old_data$observation <- round(old_data$observation, 6)
data$observation <- round(data$observation, 6)

# Create primary key
old_data$key <- paste(old_data$project_id, old_data$site_id, 
                      old_data$datetime_str, old_data$duration, 
                      old_data$variable, sep = "|")

data$key <- paste(data$project_id, data$site_id, 
                  data$datetime_str, data$duration, 
                  data$variable, sep = "|")

# Merge old_data and data to new_data
# Remove duplicates (keeping the data version)
new_data <- bind_rows(data, old_data) %>%
  distinct(key, .keep_all = TRUE) %>%
  select(-key)


# Organize by date
new_data = new_data[order(new_data$datetime), ]

# Print Row Counts for QC
cat("QC - Row Counts:\n")
cat("  Old Data: ", nrow(old_data), "\n")
cat("  New Downloaded Data:     ", nrow(data),     "\n")
cat("  New Combined Data:     ", nrow(new_data),     "\n")

# Force datetime column as string
new_data$datetime <- format(new_data$datetime, format = "%Y-%m-%d %H:%M")

# Select columns
new_data = new_data[, c('project_id', 'site_id', 'datetime', 'duration', 
                        'variable','observation')]

# Create metadata
site_metadata_df = urban_metadata_sites(copy_updated_data)
pollutant_metadata_df = urban_metadata_pollutant(copy_updated_data)

# Step 3: Write to S3 Bucket ----------------------------------------------

# Set Up S3 Connection
s3_read <- arrow::s3_bucket('bu4cast-ci',
                            endpoint_override = 'https://uri.osn.mghpcc.org',
                            access_key = Sys.getenv("OSN_KEY"),
                            secret_key = Sys.getenv("OSN_SECRET"),
                            scheme = "https")

# Write to S3 bucket
arrow::write_csv_arrow(new_data, sink = s3_read$path(filename))

# Write metadata to bucket
site_metadata_df_filename = config$target_groups$Urban$site_metadata_filepath
pollutant_metadata_df_filename = config$target_groups$Urban$units_metadata_filepath

arrow::write_csv_arrow(site_metadata_df, sink = s3_read$path(site_metadata_df_filename))
arrow::write_csv_arrow(pollutant_metadata_df, sink = s3_read$path(pollutant_metadata_df_filename))

# Step 4: Clean Up and Health Check ---------------------------------------

# Remove file from working directory
csv_filename = paste(challenge_name, "-targets.csv")
unlink(csv_filename)

# Health Check
# Created at www.healthchecks.io
# Currently set to urban-targets
RCurl::getURL("https://hc-ping.com/79b757b6-fd76-4844-aa88-ee24344e0ab7")
