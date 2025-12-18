# R/01_data_cleaning.R

# 01_data_cleaning.R
# - read extraction Excel
# - normalize registry names
# - create .study_id

library(readxl)
library(dplyr)
library(stringr)
source('R/utils.R')

raw_path <- 'data/raw/HEPATIC-T2DM_SYNTHETIC_DATA.xlsx' # place file here
if(!file.exists(raw_path)) stop('Put your raw extraction Excel at: ', raw_path)

df <- readxl::read_excel(raw_path, sheet = 1)

names(df) <- str_trim(names(df))

# normalize registry column name
registry_col <- 'IF DATA FROM LARGE CLAIMS DATASET, WHICH ONE?'
if(!(registry_col %in% names(df))) stop('Registry column not found: ', registry_col)

# standardize registries 
df <- df %>% mutate(
  registry_raw = str_trim(.data[[registry_col]]),
  registry = case_when(
    str_detect(registry_raw, regex('TriNetX', ignore_case = TRUE)) ~ 'TriNetX',
    str_detect(registry_raw, regex('MarketScan|Merative', ignore_case = TRUE)) ~ 'MarketScan',
    str_detect(registry_raw, regex('CDARS|Clinical Data Analysis and Reporting System', ignore_case = TRUE)) ~ 'CDARS',
    str_detect(registry_raw, regex('NHIRD|National Health Insurance', ignore_case = TRUE)) ~ 'NHIRD Taiwan',
    str_detect(registry_raw, regex('NHIS', ignore_case = TRUE)) ~ 'NHIS Korea',
    str_detect(registry_raw, regex('VA', ignore_case = TRUE)) ~ 'VA Corporate Data Warehouse',
    str_detect(registry_raw, regex('CPRD', ignore_case = TRUE)) ~ 'CPRD UK',
    TRUE ~ NA_character_
  ),
  # If registry is missing, treat each study as its own data source
  registry = if_else(is.na(registry) | registry == '' , as.character(AUTHOR), registry)
)

# create study id
df$.study_id <- make_study_id(df, author_col = 'AUTHOR', year_col = 'YEAR')

if(!dir.exists('data')) dir.create('data')
write.csv(df, 'data/cleaned_extraction.csv', row.names = FALSE)
message('Wrote data/cleaned_extraction.csv')
