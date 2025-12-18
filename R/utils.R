# R/utils.R

library(dplyr); library(stringr); library(tibble)

safe_get_char <- function(df_row, col) {
  if(!(col %in% names(df_row))) return(NA_character_)
  val <- df_row[[col]]
  if(length(val) == 0) return(NA_character_)
  as.character(val)[1]
}

safe_get_int <- function(df_row, col) {
  if(!(col %in% names(df_row))) return(NA_integer_)
  val <- df_row[[col]]
  if(length(val) == 0) return(NA_integer_)
  suppressWarnings(as.integer(gsub("[^0-9]", "", as.character(val))))
}

extract_hr_ci_robust <- function(text) {
  if(is.na(text) || length(text)==0) return(c(hr=NA_real_, low=NA_real_, high=NA_real_))
  txt <- as.character(text)
  txt <- gsub("–", "-", txt)
  txt <- gsub('\u2013', '-', txt) # en dash
  txt <- gsub('[^0-9.() ,:-]', ' ', txt)
  txt <- gsub("\\s+", " ", trimws(txt))

  m <- regmatches(
    txt,
    regexec(
      "([0-9]+(?:\\.[0-9]+)?)\\s*\\(\\s*([0-9]+(?:\\.[0-9]+)?)\\s*[-–—]\\s*([0-9]+(?:\\.[0-9]+)?)\\s*\\)",
      txt
    )
  )[[1]]
  if(length(m) >= 4) {
    return(c(hr = as.numeric(m[2]), low = as.numeric(m[3]), high = as.numeric(m[4])))
  }
  nums <- unlist(regmatches(txt, gregexpr('[0-9]+(?:\\.[0-9]+)?', txt)))
  if(length(nums) >= 3) return(c(hr = as.numeric(nums[1]), low = as.numeric(nums[2]), high = as.numeric(nums[3])))
  return(c(hr=NA_real_, low=NA_real_, high=NA_real_))
}

# safe numeric extractor from text (like "61.1 (11.6)")
extract_first_numeric <- function(x) {
  if(is.na(x) || x == "") return(NA_real_)
  m <- regmatches(as.character(x), regexec('([+-]?[0-9]+\\.?[0-9]*)', as.character(x)))[[1]]
  if(length(m) >= 2) return(as.numeric(m[2]))
  return(NA_real_)
}

# compose study id
make_study_id <- function(df, author_col = 'AUTHOR', year_col = 'YEAR') {
  stopifnot(author_col %in% names(df), year_col %in% names(df))
  paste0(df[[author_col]], '.', df[[year_col]], '.')
}
