# 03_descriptive.R
# MEPS High-Cost Patient Prediction Project
# Step 3: Descriptive Statistics
#
# Produces:
#   output/tables/table1_overall_vs_top5.html   — Overall vs Top 5% vs Non-top-5%
#   output/tables/table1_overall_vs_top1.html   — Overall vs Top 1% vs Non-top-1%
#   output/tables/table1_descriptives.csv       — Machine-readable version
#   output/tables/zero_expenditure_stats.txt

if (!require("survey"))    install.packages("survey")
if (!require("dplyr"))     install.packages("dplyr")
if (!require("gtsummary")) install.packages("gtsummary")
if (!require("gt"))        install.packages("gt")

library(survey)
library(dplyr)
library(gtsummary)
library(gt)

dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)

cat("Loading survey design object...\n")
design <- readRDS("data/processed/meps_design.rds")
df     <- design$variables

# -----------------------------------------------------------------------------
# 1. Readable factor labels for binary variables
# -----------------------------------------------------------------------------
df$cancer_ever_f <- factor(df$cancer_ever,
                           levels = c(0, 1),
                           labels = c("No", "Yes"))

# Cancer split (P2): active vs remission proxy via concurrent inpatient utilization
df$cancer_active_f    <- factor(df$cancer_active,
                                levels = c(0, 1),
                                labels = c("No", "Yes"))
df$cancer_remission_f <- factor(df$cancer_remission,
                                levels = c(0, 1),
                                labels = c("No", "Yes"))

# Functional limitation and concurrent inpatient utilization (P2)
df$any_functional_limitation_f <- factor(df$any_functional_limitation,
                                         levels = c(0, 1),
                                         labels = c("No", "Yes"))
df$any_inpatient_2022_f <- factor(df$any_inpatient_2022,
                                  levels = c(0, 1),
                                  labels = c("No", "Yes"))

# Use mental_health_cat if available (preferred), else build from mental_health
if ("mental_health_cat" %in% names(df)) {
  df$mental_health_f <- df$mental_health_cat
} else {
  df$mental_health_f <- factor(
    ifelse(is.na(df$mental_health), "Not observed",
           ifelse(df$mental_health == 1, "SPD", "No SPD")),
    levels = c("No SPD", "SPD", "Not observed")
  )
}

# Top-5% and Top-1% comparison groups (binary: in tail vs not)
df$top5_group <- factor(df$top5, levels = c(0, 1),
                        labels = c("Non-top 5%", "Top 5%"))
df$top1_group <- factor(df$top1, levels = c(0, 1),
                        labels = c("Non-top 1%", "Top 1%"))

# Update design with new variables
design <- update(design,
  cancer_ever_f               = df$cancer_ever_f,
  cancer_active_f             = df$cancer_active_f,
  cancer_remission_f          = df$cancer_remission_f,
  any_functional_limitation_f = df$any_functional_limitation_f,
  any_inpatient_2022_f        = df$any_inpatient_2022_f,
  mental_health_f             = df$mental_health_f,
  top5_group                  = df$top5_group,
  top1_group                  = df$top1_group
)

# -----------------------------------------------------------------------------
# 2. Variables to summarize
# -----------------------------------------------------------------------------
vars_to_summarize <- c(
  "age_cat", "sex", "income_pct_fpl", "insurance_status",
  "education", "race_ethnicity",
  # Primary morbidity metric (P3 AME-weighted, standardized)
  "morbidity_score_std",
  # Secondary morbidity (kept for continuity with legacy analysis)
  "n_chronic_conditions_cat",
  # P2: self-rated health and functional limitation
  "health_status_3",
  "any_functional_limitation_f",
  # P2: concurrent utilization
  "any_inpatient_2022_f",
  # Cancer (ever / active / remission)
  "cancer_ever_f",
  "cancer_active_f",
  "cancer_remission_f",
  # Mental health (SAQ subsample)
  "mental_health_f",
  # Spending
  "total_exp"
)

# Drop any that don't exist
missing_v <- setdiff(vars_to_summarize, names(design$variables))
if (length(missing_v) > 0) {
  warning("Dropping missing variables: ", paste(missing_v, collapse = ", "))
  vars_to_summarize <- setdiff(vars_to_summarize, missing_v)
}

stat_spec <- list(
  all_categorical() ~ "{p}%",
  total_exp         ~ "{median} ({p25}, {p75})"
)
dig_spec <- list(
  all_categorical() ~ 1,
  all_continuous()  ~ 0
)

# -----------------------------------------------------------------------------
# 3a. Table: Overall vs Top 5%
#     Aligned with README: compare high-cost tail (top5==1) to the rest
# -----------------------------------------------------------------------------
cat("Generating Table 1a: Overall vs Top 5%...\n")

tbl_overall_5 <- design %>%
  tbl_svysummary(
    include   = all_of(vars_to_summarize),
    statistic = stat_spec,
    digits    = dig_spec,
    missing   = "no"
  ) %>%
  modify_header(label = "**Variable**") %>%
  bold_labels()

tbl_by_top5 <- design %>%
  tbl_svysummary(
    by        = top5_group,
    include   = all_of(vars_to_summarize),
    statistic = stat_spec,
    digits    = dig_spec,
    missing   = "no"
  ) %>%
  modify_header(label = "**Variable**") %>%
  bold_labels()

tbl_top5_final <- tbl_merge(
  tbls        = list(tbl_overall_5, tbl_by_top5),
  tab_spanner = c("**Overall**", "**Top 5% vs Non-top 5%**")
)

gtsummary::as_gt(tbl_top5_final) %>%
  gt::gtsave(filename = "output/tables/table1_overall_vs_top5.html")

cat("Saved: output/tables/table1_overall_vs_top5.html\n")

# -----------------------------------------------------------------------------
# 3b. Table: Overall vs Top 1%
# -----------------------------------------------------------------------------
cat("Generating Table 1b: Overall vs Top 1%...\n")

tbl_by_top1 <- design %>%
  tbl_svysummary(
    by        = top1_group,
    include   = all_of(vars_to_summarize),
    statistic = stat_spec,
    digits    = dig_spec,
    missing   = "no"
  ) %>%
  modify_header(label = "**Variable**") %>%
  bold_labels()

tbl_top1_final <- tbl_merge(
  tbls        = list(tbl_overall_5, tbl_by_top1),
  tab_spanner = c("**Overall**", "**Top 1% vs Non-top 1%**")
)

gtsummary::as_gt(tbl_top1_final) %>%
  gt::gtsave(filename = "output/tables/table1_overall_vs_top1.html")

cat("Saved: output/tables/table1_overall_vs_top1.html\n")

# -----------------------------------------------------------------------------
# 4. Machine-readable CSV (basic weighted stats)
# -----------------------------------------------------------------------------
cat("Saving machine-readable descriptive stats...\n")

# Morbidity score quantiles (P3 AME-weighted, standardized)
ms_q <- quantile(df$morbidity_score_std,
                 probs = c(0.25, 0.5, 0.75),
                 na.rm = TRUE)

basic_stats <- data.frame(
  n_analytical_sample = nrow(df),
  n_top5              = sum(df$top5, na.rm = TRUE),
  pct_top5            = round(mean(df$top5, na.rm = TRUE) * 100, 1),
  n_top1              = sum(df$top1, na.rm = TRUE),
  pct_top1            = round(mean(df$top1, na.rm = TRUE) * 100, 1),
  median_total_exp    = round(median(df$total_exp, na.rm = TRUE), 0),
  n_zero_exp          = sum(df$total_exp == 0, na.rm = TRUE),
  pct_zero_exp        = round(mean(df$total_exp == 0, na.rm = TRUE) * 100, 1),
  morbidity_score_p25 = round(ms_q[[1]], 3),
  morbidity_score_p50 = round(ms_q[[2]], 3),
  morbidity_score_p75 = round(ms_q[[3]], 3)
)

write.csv(basic_stats, "output/tables/table1_descriptives.csv", row.names = FALSE)
cat("Saved: output/tables/table1_descriptives.csv\n")

# -----------------------------------------------------------------------------
# 5. Zero-expenditure analysis
# -----------------------------------------------------------------------------
cat("Analyzing zero-expenditure observations...\n")

n_total  <- nrow(df)
n_zero   <- sum(df$total_exp == 0, na.rm = TRUE)
pct_zero <- n_zero / n_total * 100

cat(sprintf("Total analytical sample (n):      %d\n", n_total))
cat(sprintf("Zero-expenditure observations:    %d (%.1f%%)\n", n_zero, pct_zero))
cat(sprintf("Positive-expenditure (for GLM):   %d (%.1f%%)\n", n_total - n_zero, 100 - pct_zero))

sink("output/tables/zero_expenditure_stats.txt")
cat(sprintf("Total analytical sample (n):      %d\n", n_total))
cat(sprintf("Zero-expenditure observations:    %d (%.1f%%)\n", n_zero, pct_zero))
cat(sprintf("Positive-expenditure (for GLM):   %d (%.1f%%)\n", n_total - n_zero, 100 - pct_zero))
sink()

cat("\nStep 3 complete.\n")
