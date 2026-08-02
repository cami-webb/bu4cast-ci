## Example Targets file for the BU Forecasting Challenges
## Created: 9/26/2025

# This template will walk you through each step that needs to be done in a generic format
# The target files should only contain the variable(s) you want the students to predict
# If you have other data to help predict those variable(s), that will go in the driver files


## Step 0: Reload data for appending
# First you must decide if it will be easier (code wise and computationally) to either
# download all of the data again each time or just append any new data
# If you are appending the data, keep the following steps to read the old data from
# the S3 bucket
# Most of this code is identical to the writing data portion in Step 4
# You only need to set up the connection and file name once, so feel free to delete
# the duplicate lines

# Set Up Connection
s3_read <- arrow::s3_bucket('bu4cast-ci/read',
                            endpoint_override = 'https://uri.osn.mghpcc.org',
                            access_key = Sys.getenv("OSN_KEY"),
                            secret_key = Sys.getenv("OSN_SECRET"),
                            scheme = "https")

# Create file name/folder
challenge_name = '' # Fill in string with name of challenge (i.e Disease, Urban, Coastal)
filename = paste("challenges/targets/project_id=bu4cast/", challenge_name,
                 "-targets.csv", sep = "")

# Read in old data
old_data <- arrow::read_csv_arrow(s3_read$path(filename))

## Step 1: Download Data
# If you need a key to access the data, create a secret key
# Github Repo -> Settings -> Secrets and Variables -> Actions -> New Repository Secret
# Access the variable (i.e MY_VAR) in your script using Sys.getenv("MY_VAR")

# Download the list of sites you want to work with


# Download data for target variable(s)


## Step 2: Combine data together for each of your sites

data = data.frame(matrix(ncol = 6, nrow = 0))

# This is the order and names you must use for the data when writing to the S3 Bucket
colnames(data) = c('project_id', 'site_id', 'datetime', 'duration', 'variable', 'observation')



## Column Explanations
# project_id = bu4cast
# site_id
# datetime
# observation = actual value of variable
# duration = how often the data is collected
  # "P1M" = monthly
  # "P1D" = daily
  # "PT30M" = 30 minute
  # General Formula found here: https://courses.ischool.berkeley.edu/i290-14/s05/lecture-15/slide31.html
# variable = which type of variable we're using (since we don't have a different column for each variable)

## Step 3: Prep data for storage
# Depending on your process for step 2, this may not be necessary
# If appending to old data: append and sort


## Step 4: Write to S3 Bucket

# Set Up S3 Connection
s3_read <- arrow::s3_bucket('bu4cast-ci/read',
                            endpoint_override = 'https://uri.osn.mghpcc.org',
                            access_key = Sys.getenv("OSN_KEY"),
                            secret_key = Sys.getenv("OSN_SECRET"),
                            scheme = "https")
# Create file name/folder
challenge_name = '' # Fill in string with name of challenge (i.e Disease, Urban, Coastal)
filename = paste("challenges/targets/project_id=bu4cast/", challenge_name,
                 "-targets.csv", sep = "")

# Write to S3 bucket
# Change data if you changed the name of the cleaned data in the script
arrow::write_csv_arrow(data, sink = s3_read$path(filename))

## Step 5: Clean Up and Health Check

# Remove file from working directory
csv_filename = paste(challenge_name, "-targets.csv.gz")
unlink(csv_filename)

# Health Check
# Created at www.healthchecks.io
# Currently set to bu4cast-ci-example
RCurl::getURL("https://hc-ping.com/9836911f-f0e5-485b-91a0-f1d98f11dd7a")
