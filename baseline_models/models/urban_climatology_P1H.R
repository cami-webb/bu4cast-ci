## Climatology Model for Urban Targets at Hourly Rate

# Hourly Variables to Work with: O3, PM2.5, PM10, NO2

# Daily Variables to Work with: PM10, PM2.5

# Set Up ------------------------------------------------------------------

# Library Imports
library(readr)
library(dplyr)
library(lubridate)

# Define constants
ne = 31

# Read in data
urban_data_url = 'https://uri.osn.mghpcc.org/bu4cast-ci/read/challenges/targets/project_id=bu4cast/urban-targets.csv'
big_urban_data = read_csv(urban_data_url, 
                          col_types = cols(project_id = col_character(),
                                           site_id = col_character(),
                                           datetime = col_date(),
                                           duration = col_character(),
                                           variable = col_character(),
                                           observation = col_double()))
# Read in all sites
sites_metadata_url = 'https://uri.osn.mghpcc.org/bu4cast-ci/read/challenges/targets/project_id=bu4cast/urban-targets-sites.csv'

sites_metadata <- read_csv(
  sites_metadata_url,
  col_types = cols(
    # ID and Location: Character
    site_id = col_character(),
    site_lat = col_character(),
    site_long = col_character(),
    
    # Start/End Date: Date
    PM2.5_P1D_StartDate = col_date(),
    PM2.5_P1D_EndDate = col_date(),
    PM2.5_P1H_StartDate = col_date(),
    PM2.5_P1H_EndDate = col_date(),
    PM10_P1D_StartDate = col_date(),
    PM10_P1D_EndDate = col_date(),
    PM10_P1H_StartDate = col_date(),
    PM10_P1H_EndDate = col_date(),
    O3_StartDate = col_date(),
    O3_EndDate = col_date(),
    NO2_P1D_StartDate = col_date(),
    NO2_P1D_EndDate = col_date(),
    NO2_P1H_StartDate = col_date(),
    NO2_P1H_EndDate = col_date(),
    
    # Active: Logical
    PM2.5_P1D_Active = col_logical(),
    PM2.5_P1H_Active = col_logical(),
    PM10_P1D_Active = col_logical(),
    PM10_P1H_Active = col_logical(),
    O3_Active = col_logical(),
    NO2_P1D_Active = col_logical(),
    NO2_P1H_Active = col_logical(),
    
    # Defaults
    .default = col_character()
  )
)


# Find Active Sites -------------------------------------------------------

# Calculate today but last year
current_doy = yday(Sys.Date())
today = Sys.Date()
last_year = year(today) - 1


# Create a list of all active sites and if there is data for last year

# PM2.5 Hourly
PM2.5_P1H_active_sites = sites_metadata %>%
  filter(
    PM2.5_P1H_Active == TRUE &
    (yday(PM2.5_P1H_StartDate) < yday(Sys.Date()) | 
    year(PM2.5_P1H_StartDate) < last_year)
  )  %>%
  pull(site_id)

# PM10 Hourly
PM10_P1H_active_sites = sites_metadata %>%
  filter(
    PM10_P1H_Active == TRUE &
      (yday(PM10_P1H_StartDate) < yday(Sys.Date()) | 
         year(PM10_P1H_StartDate) < last_year)
  )  %>%
  pull(site_id)

# O3 Hourly
O3_active_sites = sites_metadata %>%
  filter(
    O3_Active == TRUE &
      (yday(O3_StartDate) < yday(Sys.Date()) | 
         year(O3_StartDate) < last_year)
  )  %>%
  pull(site_id)

# NO2 Hourly
NO2_P1H_active_sites = sites_metadata %>%
  filter(
    NO2_P1H_Active == TRUE &
      (yday(NO2_P1H_StartDate) < yday(Sys.Date()) | 
         year(NO2_P1H_StartDate) < last_year)
  )  %>%
  pull(site_id)

# Compile all active hourly sites
active_P1H_sites = unique(c(PM2.5_P1H_active_sites,
                            PM10_P1H_active_sites,
                            O3_active_sites,
                            NO2_P1H_active_sites))

# PM2.5 Daily
PM2.5_P1D_active_sites = sites_metadata %>%
  filter(
    PM2.5_P1D_Active == TRUE &
      (yday(PM2.5_P1D_StartDate) < yday(Sys.Date()) | 
         year(PM2.5_P1D_StartDate) < last_year)
  )  %>%
  pull(site_id)

# PM2.5 Daily
PM10_P1D_active_sites = sites_metadata %>%
  filter(
    PM10_P1D_Active == TRUE &
      (yday(PM10_P1D_StartDate) < yday(Sys.Date()) | 
         year(PM10_P1D_StartDate) < last_year)
  )  %>%
  pull(site_id)

# Compile all active daily sites
active_P1D_sites = unique(c(PM2.5_P1D_active_sites,
                            PM10_P1D_active_sites))

# Climatology Forecast ----------------------------------------------------

# Which variables to use
active_P1D_variables = c('PM2.5 - Daily', 'PM10 - Daily')
active_P1H_variables = c('PM2.5 - Hourly',
                         'PM10 - Hourly',
                         'NO2 - Hourly',
                         'O3')

# Find the most recent datetime and predict from there
start_date = max(big_urban_data$datetime)
end_date = as_date(start_date + days(35))

forecast_daily_dates = seq(as.POSIXct(start_date), as.POSIXct(end_date), by = "1 day")
forecast_doy = yday(forecast_daily_dates)
forecast_hourly_times = seq(as.POSIXct(start_date), as.POSIXct(end_date), by = "hour")
forecast_doy_hourly = yday(forecast_hourly_times) + hour(forecast_hourly_times)/24.0

# Add yod to data (both daily and hourly)
big_urban_data$doy = yday(big_urban_data$datetime)
big_urban_data$doy_hourly = yday(big_urban_data$datetime) + 
  hour(big_urban_data$datetime)/24.0

# Find mean and sd for each of the sites
climatology_daily_stats <- big_urban_data %>%
  filter(
    duration == 'P1D',
    site_id %in% active_P1D_sites,
    doy %in% forecast_doy,
    variable %in% active_P1D_variables
  ) %>%
  group_by(site_id, doy, variable) %>%
  summarise(
    mean_value = mean(observation, na.rm = TRUE),
    sd_value = sd(observation, na.rm = TRUE)
  )

climatology_hourly_stats <- big_urban_data %>%
  filter(
    duration == 'PT1H',
    site_id %in% active_P1H_sites,
    doy_hourly %in% forecast_doy_hourly,
    variable %in% active_P1H_variables
  ) %>%
  group_by(site_id, doy_hourly, variable) %>%
  summarise(
    mean_value = mean(observation, na.rm = TRUE),
    sd_value = sd(observation, na.rm = TRUE)
  )

# Build 35 days out Climatology Prediction


# Columns: datetime, duration, site_id, parameter, variable, prediction

# Formatting for Submission -----------------------------------------------

# Columns: project_id, model_id, datetime, reference_datetime, duration, site_id,
# family, parameter, variable, prediction

project_id = 'bu4cast'
model_id = 'bu4cast-test' # same as team name
reference_datetime = start_date # what is the datetime from which you are starting your prediction
family = 'normal'
