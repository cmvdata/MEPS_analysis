# 02_clean.R
# MEPS High-Cost Patient Prediction Project
# Step 2: Variable Construction and Survey Design

# -----------------------------------------------------------------------------
# 0. Setup
# -----------------------------------------------------------------------------

if (!require("haven"))  install.packages("haven")
if (!require("survey")) install.packages("survey")
if (!require("dplyr"))  install.packages("dplyr")

library(haven)
library(survey)
library(dplyr)

options(survey.lonely.psu = "adjust")

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables",  recursive = TRUE, showWarnings = FALSE)

input_path  <- "data/processed/meps_adults_raw.rds"
output_path <- "data/processed/meps_analytical.rds"
design_path <- "data/processed/meps_design.rds"

cat("Loading filtered dataset...\n")
df <- readRDS(input_path)
cat(sprintf("Loaded: %d rows, %d columns\n", nrow(df), ncol(df)))

# -----------------------------------------------------------------------------
# 1. Variable checks (CRITICAL)
# -----------------------------------------------------------------------------

required_vars <- c(
  "TOTEXP22", "AGE22X", "SEX", "POVCAT22", "INSURC22",
  "VARSTR", "VARPSU", "PERWT22F"
)

missing_vars <- setdiff(required_vars, names(df))
if (length(missing_vars) > 0) {
  stop("Missing required variables: ", paste(missing_vars, collapse = ", "))
}

cat("All required variables found.\n")

# -----------------------------------------------------------------------------
# 2. Clean negative codes (MEPS convention)
#    -1 = inapplicable, -7 = refused, -8 = don't know, -9 = not ascertained
# -----------------------------------------------------------------------------

clean_negatives <- function(x) {
  if (is.numeric(x)) x[x < 0] <- NA
  return(x)
}

df <- df %>% mutate(across(everything(), clean_negatives))
cat("Negative MEPS codes recoded to NA.\n")

# -----------------------------------------------------------------------------
# 3. Core variables
# -----------------------------------------------------------------------------

df$total_exp <- df$TOTEXP22

df$age_cat <- cut(
  df$AGE22X,
  breaks = c(18, 35, 55, 65, 75, Inf),
  labels = c("18-34", "35-54", "55-64", "65-74", "75+"),
  right  = FALSE
)

df$sex <- factor(df$SEX, levels = c(1, 2), labels = c("Male", "Female"))

df$income_pct_fpl <- factor(
  df$POVCAT22,
  levels = 1:5,
  labels = c("<100% FPL", "100-125% FPL", "125-200% FPL", "200-400% FPL", ">=400% FPL")
)
cat("income_pct_fpl distribution (FPL labels):\n"); print(table(df$income_pct_fpl, useNA = "ifany"))

# -----------------------------------------------------------------------------
# 4. Insurance status
#
# INSURC22 codebook HC-243:
#   1 = <65 Any Private
#   2 = <65 Public Only
#   3 = <65 Uninsured
#   4 = 65+ Medicare Only
#   5 = 65+ Medicare + Private
#   6 = 65+ Medicare + Other Public
#   7 = 65+ No Medicare + Any Public/Private  <- has SOME private/public, not uninsured
#   8 = 65+ Uninsured
#
# Conservative collapse (codebook HC-243):
#   Private              : 1, 5     (any private, <65 and 65+ Medicare+Private)
#   Medicare/Other Public: 4, 6, 7  (65+ Medicare only, Medicare+Other Public,
#                                    65+ No Medicare + Any Public/Private —
#                                    code 7 is ambiguous; kept here as public
#                                    to avoid overstating private coverage)
#   Medicaid/Public Only : 2        (<65 public only)
#   Uninsured            : 3, 8
# -----------------------------------------------------------------------------

cat("Inspecting INSURC22 distribution...\n")
print(table(df$INSURC22, useNA = "ifany"))

df$insurance_status <- case_when(
  df$INSURC22 %in% c(1, 5)    ~ "Private",
  df$INSURC22 %in% c(4, 6, 7) ~ "Medicare/Other Public",
  df$INSURC22 == 2             ~ "Medicaid/Public Only",
  df$INSURC22 %in% c(3, 8)    ~ "Uninsured",
  TRUE                         ~ NA_character_
)

df$insurance_status <- factor(
  df$insurance_status,
  levels = c("Private", "Medicare/Other Public", "Medicaid/Public Only", "Uninsured")
)

cat("insurance_status distribution:\n")
print(table(df$insurance_status, useNA = "ifany"))

# -----------------------------------------------------------------------------
# 5. Education (robust handling)
# -----------------------------------------------------------------------------

if ("EDUCYR" %in% names(df)) {
  cat("EDUCYR found. Distribution:\n")
  print(table(df$EDUCYR, useNA = "ifany"))
  df$education <- case_when(
    df$EDUCYR < 12                   ~ "Less than HS",
    df$EDUCYR == 12                  ~ "HS",
    df$EDUCYR > 12 & df$EDUCYR < 16 ~ "Some college",
    df$EDUCYR >= 16                  ~ "College+",
    TRUE                             ~ NA_character_
  )
  df$education <- factor(df$education,
                         levels = c("Less than HS", "HS", "Some college", "College+"))
} else {
  # Education missing entirely: create explicit "Not observed" level so
  # downstream models do not silently drop all rows on complete-case.
  warning("EDUCYR not found. education set to 'Not observed' for all rows.")
  df$education <- factor("Not observed",
                         levels = c("Less than HS", "HS", "Some college", "College+", "Not observed"))
}

# -----------------------------------------------------------------------------
# 6. Race / Ethnicity
# -----------------------------------------------------------------------------

df$race_ethnicity <- case_when(
  df$HISPANX == 1     ~ "Hispanic",
  df$RACEV2X == 1     ~ "White Non-Hispanic",
  df$RACEV2X == 2     ~ "Black Non-Hispanic",
  df$RACEV2X == 4     ~ "Asian Non-Hispanic",
  !is.na(df$RACEV2X) ~ "Other/Multiple Non-Hispanic",
  TRUE                ~ NA_character_
)

df$race_ethnicity <- factor(
  df$race_ethnicity,
  levels = c("White Non-Hispanic", "Black Non-Hispanic", "Hispanic",
             "Asian Non-Hispanic", "Other/Multiple Non-Hispanic")
)

# -----------------------------------------------------------------------------
# 7. Morbidity — resolved variable by variable
#
# In HC-243, conditions diagnosed after age 18 use suffix _M18.
# Each condition is resolved independently: _M18 first, plain name as fallback.
# This avoids the block-replacement bug where finding one missing _M18
# would discard all valid _M18 variables.
# Coding: 1 = Yes, 2 = No (negatives already cleaned above)
# -----------------------------------------------------------------------------

resolve_chronic <- function(df, name_m18, name_plain) {
  if (name_m18   %in% names(df)) return(df[[name_m18]])
  if (name_plain %in% names(df)) return(df[[name_plain]])
  warning(sprintf("Neither %s nor %s found. Setting to NA.", name_m18, name_plain))
  rep(NA_real_, nrow(df))
}

df$diab_dx <- resolve_chronic(df, "DIABDX_M18", "DIABDX")
df$hibp_dx <- resolve_chronic(df, "HIBPDX_M18", "HIBPDX")
df$chd_dx  <- resolve_chronic(df, "CHDDX_M18",  "CHDDX")
df$strk_dx <- resolve_chronic(df, "STRKDX_M18", "STRKDX")
df$asth_dx <- resolve_chronic(df, "ASTHDX",     "ASTHDX")
df$arth_dx <- resolve_chronic(df, "ARTHDX",     "ARTHDX")

chronic_cols <- c("diab_dx", "hibp_dx", "chd_dx", "strk_dx", "asth_dx", "arth_dx")

# Convert: 1 = Yes → 1, 2 = No → 0, NA → NA
chronic_binary <- lapply(chronic_cols, function(v) {
  x <- df[[v]]
  ifelse(is.na(x), NA_real_, ifelse(x == 1, 1, 0))
})
names(chronic_binary) <- chronic_cols

# Count observed conditions and positive diagnoses separately.
# If at least one condition is observed, count positives over observed ones.
# Only set NA if ALL six conditions are missing (extremely rare in MEPS).
chronic_df     <- as.data.frame(chronic_binary)
n_obs_chronic  <- rowSums(!is.na(chronic_df))
n_yes_chronic  <- rowSums(chronic_df == 1, na.rm = TRUE)
df$n_chronic_conditions <- ifelse(n_obs_chronic == 0, NA_real_, n_yes_chronic)

df$n_chronic_conditions_cat <- case_when(
  df$n_chronic_conditions == 0          ~ "0",
  df$n_chronic_conditions == 1          ~ "1",
  df$n_chronic_conditions %in% c(2, 3) ~ "2-3",
  df$n_chronic_conditions >= 4         ~ "4+",
  TRUE                                  ~ NA_character_
)
df$n_chronic_conditions_cat <- factor(
  df$n_chronic_conditions_cat,
  levels = c("0", "1", "2-3", "4+")
)

# Cancer
if ("CANCERDX" %in% names(df)) {
  df$cancer_ever <- case_when(
    df$CANCERDX == 1 ~ 1L,
    df$CANCERDX == 2 ~ 0L,
    TRUE             ~ NA_integer_
  )
} else {
  df$cancer_ever <- NA_integer_
}

# Mental health — K6SUM22 (Kessler-6, serious psychological distress >= 13)
# If K6SUM22 is absent in the FYC, download and merge SAQ file (HC-245)
get_k6 <- function(df) {
  # Check all known K6 variable names across MEPS years
  k6_candidates <- c("K6SUM22", "K6SUM42", "K6SUM21", "K6SUM20")
  for (varname in k6_candidates) {
    if (varname %in% names(df) && sum(!is.na(df[[varname]])) > 100) {
      cat(sprintf("%s found in FYC: %d non-missing values.\n",
                  varname, sum(!is.na(df[[varname]]))))
      return(df[[varname]])
    }
  }
  cat("No K6 variable found in FYC. Downloading SAQ file (HC-245)...\n")
  saq_url <- "https://meps.ahrq.gov/mepsweb/data_files/pufs/h245/h245dta.zip"
  saq_zip <- "data/raw/h245dta.zip"
  saq_dta <- "data/raw/h245.dta"
  if (!file.exists(saq_dta)) {
    if (!file.exists(saq_zip)) {
      download.file(saq_url, destfile = saq_zip, mode = "wb")
    }
    unzip(saq_zip, exdir = "data/raw")
    if (!file.exists(saq_dta)) stop("SAQ file not found after unzip: ", saq_dta)
  }
  saq <- read_dta(saq_dta)
  if (!"K6SUM22" %in% names(saq)) {
    warning("K6SUM22 not found in SAQ either. mental_health = NA.")
    return(rep(NA_real_, nrow(df)))
  }
  saq <- saq[, c("DUPERSID", "K6SUM22")]
  saq$K6SUM22[saq$K6SUM22 < 0] <- NA
  merged <- merge(data.frame(DUPERSID = df$DUPERSID), saq,
                  by = "DUPERSID", all.x = TRUE)
  merged <- merged[match(df$DUPERSID, merged$DUPERSID), ]
  cat(sprintf("K6SUM22 merged: %d non-missing (%.1f%%).\n",
              sum(!is.na(merged$K6SUM22)),
              mean(!is.na(merged$K6SUM22)) * 100))
  return(merged$K6SUM22)
}

df$K6SUM22_clean <- get_k6(df)

# IMPORTANT: NA stays NA — do NOT recode missing to 0
df$mental_health <- case_when(
  df$K6SUM22_clean >= 13 ~ 1L,
  df$K6SUM22_clean <  13 ~ 0L,
  TRUE                   ~ NA_integer_
)

cat(sprintf("mental_health: %d SPD, %d non-SPD, %d missing\n",
            sum(df$mental_health == 1, na.rm = TRUE),
            sum(df$mental_health == 0, na.rm = TRUE),
            sum(is.na(df$mental_health))))

# Categorical version — keeps NA as explicit level for complete-case safety
df$mental_health_cat <- factor(
  ifelse(is.na(df$mental_health), "Not observed",
         ifelse(df$mental_health == 1, "SPD", "No SPD")),
  levels = c("No SPD", "SPD", "Not observed")
)

# -----------------------------------------------------------------------------
# 7B. Self-rated health status & functional limitation
# -----------------------------------------------------------------------------

# RTHLTH31: self-rated general health, round 1 (entry into 2022).
# Values: 1=Excellent, 2=Very good, 3=Good, 4=Fair, 5=Poor.
df$health_status_3 <- factor(
  case_when(
    df$RTHLTH31 %in% c(1, 2) ~ "Excellent/Very good",
    df$RTHLTH31 == 3         ~ "Good",
    df$RTHLTH31 %in% c(4, 5) ~ "Fair/Poor",
    TRUE                      ~ NA_character_
  ),
  levels = c("Excellent/Very good", "Good", "Fair/Poor")
)

# Functional limitation: needs help with ADL or IADL.
# ADLHLP31/IADLHP31 are round-1 versions; HC-243 lacks round-4
# equivalents (only desaggregated items exist).
df$any_functional_limitation <- case_when(
  df$ADLHLP31 == 1 | df$IADLHP31 == 1 ~ 1L,
  df$ADLHLP31 == 2 & df$IADLHP31 == 2 ~ 0L,
  TRUE                                ~ NA_integer_
)

cat("\nhealth_status_3 distribution:\n"); print(table(df$health_status_3, useNA = "ifany"))
cat(sprintf("any_functional_limitation: yes=%d  no=%d  NA=%d\n",
            sum(df$any_functional_limitation == 1, na.rm = TRUE),
            sum(df$any_functional_limitation == 0, na.rm = TRUE),
            sum(is.na(df$any_functional_limitation))))

# -----------------------------------------------------------------------------
# 7C. Concurrent utilization (calendar year 2022)
#
# These describe contemporaneous resource use; admissible because the
# analysis is cross-sectional within 2022 (see README "Why cross-sectional").
# HC-243 does not ship pre-computed ANYOP22 / ANYER22 / RXUSE22 flags;
# we construct them from the count variables that ARE shipped.
# -----------------------------------------------------------------------------

df$n_inpatient_2022   <- df$IPDIS22
df$any_inpatient_2022 <- as.integer(df$IPDIS22 > 0)
df$any_outpatient     <- as.integer(df$OBTOTV22 > 0 | df$OPTOTV22 > 0)
df$any_er             <- as.integer(df$ERTOT22  > 0)
df$any_rx_use         <- as.integer(df$RXTOT22  > 0)

cat("\nUtilization (any_X among adults 18+):\n")
cat(sprintf("  inpatient  : %d (%.1f%%)\n",
            sum(df$any_inpatient_2022, na.rm = TRUE),
            mean(df$any_inpatient_2022, na.rm = TRUE) * 100))
cat(sprintf("  outpatient : %d (%.1f%%)\n",
            sum(df$any_outpatient, na.rm = TRUE),
            mean(df$any_outpatient, na.rm = TRUE) * 100))
cat(sprintf("  ER         : %d (%.1f%%)\n",
            sum(df$any_er, na.rm = TRUE),
            mean(df$any_er, na.rm = TRUE) * 100))
cat(sprintf("  Rx         : %d (%.1f%%)\n",
            sum(df$any_rx_use, na.rm = TRUE),
            mean(df$any_rx_use, na.rm = TRUE) * 100))

# -----------------------------------------------------------------------------
# 7D. Cancer active vs remission — proxy via concurrent inpatient utilization
#
# CANCER ACTIVE PROXY (Producto 2.6)
# ---------------------------------------------------------
# Approximation rationale (no HC-242B condition-level data):
#
# We classify cancer_ever cases as 'active' if they had any
# inpatient stay in 2022, 'remission' otherwise.
#
# KNOWN LIMITATIONS:
#  - Patients on outpatient chemotherapy without admission
#    are misclassified as remission (false negatives for active).
#  - Patients hospitalized in 2022 for non-cancer reasons
#    are misclassified as active (false positives).
#  - The split is therefore noisy. Use as descriptor only,
#    not as causal subgroup.
#
# REFINEMENT PATH:
#  - Download HC-242B Medical Conditions File 2022 (~30 MB).
#  - Filter conditions with CCS codes in oncology range.
#  - Match to inpatient/outpatient events at encounter level.
# -----------------------------------------------------------------------------

df$cancer_active <- case_when(
  is.na(df$cancer_ever) | is.na(df$any_inpatient_2022) ~ NA_integer_,
  df$cancer_ever == 1 & df$any_inpatient_2022 == 1     ~ 1L,
  TRUE                                                  ~ 0L
)
df$cancer_remission <- case_when(
  is.na(df$cancer_ever) | is.na(df$any_inpatient_2022) ~ NA_integer_,
  df$cancer_ever == 1 & df$any_inpatient_2022 == 0     ~ 1L,
  TRUE                                                  ~ 0L
)

n_cancer_ever      <- sum(df$cancer_ever == 1,      na.rm = TRUE)
n_cancer_active    <- sum(df$cancer_active == 1,    na.rm = TRUE)
n_cancer_remission <- sum(df$cancer_remission == 1, na.rm = TRUE)
pct_active <- if (n_cancer_ever > 0) 100 * n_cancer_active / n_cancer_ever else NA_real_

cat(sprintf("\nCancer split among cancer_ever (n=%d):\n", n_cancer_ever))
cat(sprintf("  active    : %d (%.1f%% of cancer_ever)\n", n_cancer_active, pct_active))
cat(sprintf("  remission : %d (%.1f%% of cancer_ever)\n", n_cancer_remission, 100 - pct_active))

if (!is.na(pct_active) && pct_active < 5) {
  warning(sprintf(
    "cancer_active = %.1f%% of cancer_ever (<5%%). Likely under-detection: outpatient chemo without admission is misclassified as remission. Consider adding a high-Rx or high-cost criterion as secondary classifier.",
    pct_active
  ))
}

# -----------------------------------------------------------------------------
# 8. Survey design
# -----------------------------------------------------------------------------

cat("Creating survey design...\n")

design <- svydesign(
  ids     = ~VARPSU,
  strata  = ~VARSTR,
  weights = ~PERWT22F,
  data    = df,
  nest    = TRUE
)

# -----------------------------------------------------------------------------
# 9. Thresholds (survey-weighted)
# -----------------------------------------------------------------------------

cat("Computing weighted thresholds...\n")

q95 <- svyquantile(~total_exp, design, quantiles = 0.95, na.rm = TRUE)
q99 <- svyquantile(~total_exp, design, quantiles = 0.99, na.rm = TRUE)

thresh_95 <- as.numeric(coef(q95))
thresh_99 <- as.numeric(coef(q99))

cat(sprintf("Top 5%% threshold: $%.0f\n", thresh_95))
cat(sprintf("Top 1%% threshold: $%.0f\n", thresh_99))

write.csv(
  data.frame(top5_threshold = thresh_95, top1_threshold = thresh_99),
  "output/tables/thresholds.csv",
  row.names = FALSE
)

df$top5 <- as.integer(!is.na(df$total_exp) & df$total_exp >= thresh_95)
df$top1 <- as.integer(!is.na(df$total_exp) & df$total_exp >= thresh_99)

design <- update(design, top5 = df$top5, top1 = df$top1)

# -----------------------------------------------------------------------------
# 10. Zero expenditure report
# -----------------------------------------------------------------------------

n_zero <- sum(df$total_exp == 0, na.rm = TRUE)
p_zero <- mean(df$total_exp == 0, na.rm = TRUE)

cat(sprintf("Zero expenditure: %d obs (%.1f%%)\n", n_zero, p_zero * 100))

write.csv(
  data.frame(n_total  = nrow(df),
             n_zero   = n_zero,
             pct_zero = round(p_zero * 100, 1),
             n_pos    = nrow(df) - n_zero),
  "output/tables/zero_expenditure_stats.csv",
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 11. Save outputs
# -----------------------------------------------------------------------------

saveRDS(df,     output_path)
saveRDS(design, design_path)

cat("Saved:", output_path, "\n")
cat("Saved:", design_path, "\n")

cat(sprintf("\nFinal sample: %d rows | top5: %d (%.1f%%) | top1: %d (%.1f%%)\n",
            nrow(df),
            sum(df$top5), mean(df$top5) * 100,
            sum(df$top1), mean(df$top1) * 100))

cat("Step 2 complete.\n")
