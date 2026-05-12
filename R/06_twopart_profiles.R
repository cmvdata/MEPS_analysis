# =============================================================================
# 06_twopart_profiles.R — Two-part expected spending by profile (post-P4 refactor)
# =============================================================================
#
# Methodology:
#   - INTENSIVE margin (Gamma log-link): reuses model_main from
#     output/models/p4_interaction_models.rds (Producto 4). No refit.
#   - EXTENSIVE margin (quasibinomial logit): new fit. Specification matches
#     model_main EXCEPT it OMITS any_inpatient_2022. Reason: an inpatient stay
#     in 2022 implies total_exp > 0 by construction, so including it as a
#     predictor of Pr(Y > 0) is circular endogeneity.
#   - Profiles: 4 insurance categories x 3 morbidity levels (z = -1, 0, +1).
#     Reference for other covariates: age 45, Female, Middle Income (POVCAT22=4),
#     South, Good self-rated health, no functional limitation, no inpatient 2022.
#
# Outputs:
#   output/tables/table_twopart_profiles.csv
#   output/tables/table_twopart_profiles.html
#   output/models/participation_model.rds
# =============================================================================

suppressPackageStartupMessages({
  library(survey)
  library(splines)
  library(marginaleffects)
  library(dplyr)
  library(readr)
  library(tibble)
  library(gt)
})
options(survey.lonely.psu = "adjust")

dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("output/models", recursive = TRUE, showWarnings = FALSE)

# ---- 0. Load --------------------------------------------------------------
cat("Loading data, design, and Producto 4 model...\n")
meps        <- readRDS("data/processed/meps_analytical.rds")
meps_design <- readRDS("data/processed/meps_design.rds")
p4          <- readRDS("output/models/p4_interaction_models.rds")
model_main  <- p4$model_main  # intensive margin: Gamma log-link

# Strip haven_labelled / vctrs class so predict / model.matrix do not trip
strip_attrs <- function(x) { attributes(x) <- NULL; x }
for (v in c("POVCAT22", "AGE22X", "total_exp",
            "any_inpatient_2022", "any_functional_limitation")) {
  meps[[v]] <- strip_attrs(meps[[v]])
}

# Standardize morbidity score (consistent with Producto 4)
meps$morbidity_score_std <- as.numeric(scale(meps$morbidity_score))

# Build participation outcome: any_exp = 1[total_exp > 0]
meps$any_exp <- as.integer(meps$total_exp > 0)

# Regenerate design with new column
meps_design <- svydesign(
  id = ~VARPSU, strata = ~VARSTR, weights = ~PERWT22F,
  data = meps, nest = TRUE
)

# ---- 1. Extensive margin: Pr(Y > 0) ---------------------------------------
# CRITICAL: any_inpatient_2022 is OMITTED from the participation logit to
# avoid circular endogeneity (inpatient stay implies total_exp > 0).
cat("\nFitting extensive-margin model (Pr(Y > 0))...\n")

formula_ext <- any_exp ~ insurance_status * morbidity_score_std +
               splines::ns(AGE22X, df = 4) + sex +
               POVCAT22 + REGION22 +
               health_status_3 + any_functional_limitation

model_ext <- svyglm(formula_ext, design = meps_design,
                    family = quasibinomial(link = "logit"))

n_ext <- length(model_ext$residuals)
n_full <- nrow(meps)
cat(sprintf("Extensive model: n_eff = %d (of %d total adults; %.1f%% retained)\n",
            n_ext, n_full, 100 * n_ext / n_full))

saveRDS(model_ext, "output/models/participation_model.rds")
cat("Saved: output/models/participation_model.rds\n")

# ---- 2. Profile grid -------------------------------------------------------
# 4 insurance x 3 morbidity levels = 12 rows
profiles <- expand.grid(
  insurance_status    = levels(meps$insurance_status),
  morbidity_score_std = c(-1, 0, 1),
  stringsAsFactors    = FALSE
)
profiles$AGE22X                    <- 45
profiles$sex                       <- "Female"
profiles$POVCAT22                  <- 4    # numeric: ~200-400% FPL Middle Income
profiles$REGION22                  <- "South"
profiles$health_status_3           <- "Good"
profiles$any_functional_limitation <- 0
profiles$any_inpatient_2022        <- 0    # only used by intensive margin

# Coerce factors to match model frames
profiles$insurance_status <- factor(profiles$insurance_status,
                                    levels = levels(meps$insurance_status))
profiles$sex              <- factor(profiles$sex,
                                    levels = levels(meps$sex))
profiles$REGION22         <- factor(profiles$REGION22,
                                    levels = levels(meps$REGION22))
profiles$health_status_3  <- factor(profiles$health_status_3,
                                    levels = levels(meps$health_status_3))

cat("\nProfile grid (12 rows):\n")
print(profiles[, c("insurance_status", "morbidity_score_std")])

# ---- 3. Predictions --------------------------------------------------------
# Workaround for svyglm not preserving predvars for splines::ns(): we stack
# the profiles onto meps before predict so ns() can recompute knots from the
# full AGE22X distribution (instead of the single value AGE22X=45 in profiles
# alone). After prediction we extract only the profile rows.
cat("\nGenerating predictions (stack-and-extract for spline knots)...\n")

cols_needed <- intersect(names(profiles), names(meps))
meps_pad <- meps[, cols_needed, drop = FALSE]
# Filter to complete cases on the columns the formulas need so model.matrix
# does not drop rows downstream (which would break profile_idx).
keep_meps <- stats::complete.cases(meps_pad)
meps_pad  <- meps_pad[keep_meps, , drop = FALSE]
combined  <- dplyr::bind_rows(meps_pad, profiles[, cols_needed, drop = FALSE])

n_meps <- nrow(meps_pad)
n_prof <- nrow(profiles)
profile_idx <- (n_meps + 1):(n_meps + n_prof)
cat(sprintf("Combined data: %d meps complete-case rows + %d profile rows = %d total\n",
            n_meps, n_prof, nrow(combined)))

# Build design matrices from the combined data so splines::ns gets enough
# variation in AGE22X to compute its knots. Then extract just the profile
# rows. This bypasses predict.svyglm / predictions() entirely and preserves
# row order with certainty.
mm_ext_combined  <- model.matrix(delete.response(terms(model_ext)),  data = combined)
mm_int_combined  <- model.matrix(delete.response(terms(model_main)), data = combined)

stopifnot(nrow(mm_ext_combined)  == nrow(combined),
          nrow(mm_int_combined)  == nrow(combined))

mm_ext_prof <- mm_ext_combined[profile_idx, , drop = FALSE]
mm_int_prof <- mm_int_combined[profile_idx, , drop = FALSE]

# Linear predictor on link scale, then inverse link
eta_ext <- as.numeric(mm_ext_prof %*% coef(model_ext))
eta_int <- as.numeric(mm_int_prof %*% coef(model_main))

profiles$Pr_Y_pos      <- 1 / (1 + exp(-eta_ext))    # logit inverse
profiles$E_Y_given_pos <- exp(eta_int)               # log inverse

# SE on link scale via vcov, then delta-method to response scale
SE_eta_ext <- sqrt(diag(mm_ext_prof %*% vcov(model_ext)  %*% t(mm_ext_prof)))
SE_eta_int <- sqrt(diag(mm_int_prof %*% vcov(model_main) %*% t(mm_int_prof)))

# Delta method:
#   Pr(Y>0):    d(plogis(eta))/d(eta)  = p * (1 - p)
#   E[Y|Y>0]:   d(exp(eta))/d(eta)     = exp(eta) = mu
profiles$Pr_Y_pos_SE      <- profiles$Pr_Y_pos * (1 - profiles$Pr_Y_pos) * SE_eta_ext
profiles$E_Y_given_pos_SE <- profiles$E_Y_given_pos * SE_eta_int

profiles$total_E_Y <- profiles$Pr_Y_pos * profiles$E_Y_given_pos

# Approximate delta-method SE for the product (assumes independence between
# the two model fits — covariance not propagated since they are estimated
# on different objectives).
profiles$total_E_Y_SE <- sqrt(
  (profiles$Pr_Y_pos_SE * profiles$E_Y_given_pos)^2 +
  (profiles$Pr_Y_pos    * profiles$E_Y_given_pos_SE)^2
)

# ---- 4. Output table -------------------------------------------------------
out <- profiles %>%
  mutate(morbidity_level = case_when(
    morbidity_score_std == -1 ~ "Low (z = -1)",
    morbidity_score_std ==  0 ~ "Medium (z = 0)",
    morbidity_score_std ==  1 ~ "High (z = +1)"
  )) %>%
  select(insurance_status, morbidity_level, morbidity_score_std,
         Pr_Y_pos,         Pr_Y_pos_SE,
         E_Y_given_pos,    E_Y_given_pos_SE,
         total_E_Y,        total_E_Y_SE) %>%
  arrange(insurance_status, morbidity_score_std)

write_csv(out, "output/tables/table_twopart_profiles.csv")
cat("\nSaved: output/tables/table_twopart_profiles.csv\n")

# ---- 5. HTML table ---------------------------------------------------------
fmt_pct <- function(x, se) sprintf("%.1f%% (%.1f)", 100 * x, 100 * se)
fmt_dol <- function(x, se) sprintf("$%s ($%s)",
                                   formatC(round(x),  format = "d", big.mark = ","),
                                   formatC(round(se), format = "d", big.mark = ","))

tab_html <- out %>%
  mutate(
    `Pr(Y > 0)`    = fmt_pct(Pr_Y_pos,      Pr_Y_pos_SE),
    `E[Y | Y > 0]` = fmt_dol(E_Y_given_pos, E_Y_given_pos_SE),
    `E[Y]`         = fmt_dol(total_E_Y,     total_E_Y_SE),
    Insurance      = as.character(insurance_status),
    Morbidity      = morbidity_level
  ) %>%
  select(Insurance, Morbidity, `Pr(Y > 0)`, `E[Y | Y > 0]`, `E[Y]`)

gt(tab_html) %>%
  tab_header(
    title = "Two-part expected healthcare spending by profile",
    subtitle = paste0(
      "Reference profile: age 45, Female, Middle Income (200-400% FPL), ",
      "South, Good self-rated health, no functional limitation, no inpatient stay 2022"
    )
  ) %>%
  tab_source_note(source_note = paste0(
    "Intensive margin (E[Y|Y>0]): Gamma GLM log-link, reuses model_main from Producto 4. ",
    "Extensive margin (Pr[Y>0]): survey-weighted quasibinomial logit, OMITS any_inpatient_2022 ",
    "to avoid circular endogeneity. ",
    "E[Y] = Pr(Y>0) x E[Y|Y>0]. SE via approximate delta method assuming independence between ",
    "the two stages."
  )) %>%
  gtsave("output/tables/table_twopart_profiles.html")
cat("Saved: output/tables/table_twopart_profiles.html\n")

# ---- 6. Console summary ----------------------------------------------------
cat("\n=============================================================\n")
cat("Two-part expected spending - 12 profiles\n")
cat("=============================================================\n")
print(
  out %>% transmute(
    Insurance  = insurance_status,
    Morbidity  = morbidity_level,
    `Pr(Y>0)`  = sprintf("%.3f", Pr_Y_pos),
    `E[Y|Y>0]` = sprintf("$%s", formatC(round(E_Y_given_pos), big.mark=",")),
    `E[Y]`     = sprintf("$%s", formatC(round(total_E_Y),     big.mark=","))
  ),
  row.names = FALSE
)

cat("\nStep 6 complete.\n")
