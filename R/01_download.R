# =============================================================================
# 01_download.R
# MEPS High-Cost Patient Prediction Project
# Step 1: Data Ingestion
#
# Downloads MEPS 2022 Full-Year Consolidated File (HC-243),
# filters to adults 18+, and saves output.
#
# Outputs:
#   data/processed/meps_adults_raw.rds
#   data/raw/download_log.csv
# =============================================================================

library(haven)

# -----------------------------------------------------------------------------
# 0. Directory setup
# -----------------------------------------------------------------------------
dir.create("data/raw",       recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. Download and unzip
# -----------------------------------------------------------------------------
meps_url <- "https://meps.ahrq.gov/mepsweb/data_files/pufs/h243/h243dta.zip"
zip_file <- "data/raw/h243dta.zip"
dta_file <- "data/raw/h243.dta"

if (!file.exists(dta_file)) {
  cat("Downloading MEPS 2022 FYC (HC-243)...\n")
  download.file(meps_url, destfile = zip_file, mode = "wb")
  unzip(zip_file, exdir = "data/raw")
  file.remove(zip_file)
} else {
  cat("MEPS 2022 FYC file already exists locally.\n")
}

if (!file.exists(dta_file)) {
  stop("Download/unzip failed: h243.dta was not found after extraction.")
}

# -----------------------------------------------------------------------------
# 2. Load raw data
# -----------------------------------------------------------------------------
cat("Loading raw data...\n")
df_raw <- read_dta(dta_file)
n_raw  <- nrow(df_raw)
cat(sprintf("Raw sample size (all ages): %d\n", n_raw))

# -----------------------------------------------------------------------------
# 3. Verify age variable
# -----------------------------------------------------------------------------
if (!"AGE22X" %in% names(df_raw)) {
  stop("AGE22X not found. Check variable names with: names(df_raw)")
}

n_age_na <- sum(is.na(df_raw$AGE22X))
if (n_age_na > 0) cat(sprintf("Warning: %d NA values in AGE22X — excluded.\n", n_age_na))

# -----------------------------------------------------------------------------
# 4. Filter adults 18+
# -----------------------------------------------------------------------------
df_adults <- df_raw[!is.na(df_raw$AGE22X) & df_raw$AGE22X >= 18, ]
n_adults  <- nrow(df_adults)

cat(sprintf("Adult sample (18+):         %d\n", n_adults))
cat(sprintf("Excluded (minors + NA age): %d\n", n_raw - n_adults))

# -----------------------------------------------------------------------------
# 5. Save
# -----------------------------------------------------------------------------
saveRDS(df_adults, "data/processed/meps_adults_raw.rds")
cat("Saved: data/processed/meps_adults_raw.rds\n")

write.csv(
  data.frame(source_url    = meps_url,
             download_date = as.character(Sys.Date()),
             raw_rows      = n_raw,
             adult_rows    = n_adults),
  "data/raw/download_log.csv",
  row.names = FALSE
)
cat("Saved: data/raw/download_log.csv\n")

cat("\nStep 1 complete.\n")
