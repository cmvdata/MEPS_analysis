# =============================================================================
# 02b_morbidity_score.R — AME-weighted morbidity score (revision 2)
# =============================================================================
#
# Revision 2 changes (post chd/strk audit):
#   - any_inpatient_2022 enters the GLM as a CONTROL (not as a score
#     component). Its inclusion separates the cancer_active effect from
#     the inpatient-stay-per-se effect, since cancer_active = cancer_ever
#     AND any_inpatient_2022 by Producto 2 construction.
#   - AMEs are still extracted for the 8 morbidity conditions only;
#     any_inpatient_2022 is NOT in the score.
#   - Proportional weights scaled to [0.5, 3.0] for numerical stability:
#         w_scaled = 0.5 + (w_prop - min) / (max - min) * 2.5
#     The score uses w_scaled (not w_prop).
#   - A sensitivity score morbidity_score_no_cancer is computed by
#     dropping cancer_active from the sum of weighted indicators.
#
# CAVEAT: pesos calibrados sobre subset total_exp > 0. Aplicarlos a la
# población completa es proyección, no observación, para pacientes con
# gasto cero. Documentado en limitations.
#
# -- Variable mappings (verified vs HC-243 analytical dataset) --------------
#   age          -> AGE22X (numeric 18-85)
#   insurance    -> insurance_status (factor)
#   REGION       -> REGION22 (haven_labelled, 1-4)
#   POVCAT22     -> income_pct_fpl (factor with FPL labels — same data)
#   URBAN22      :: NOT in analytical dataset; OMITTED from model.
#
# -- Encoding fixes ---------------------------------------------------------
#   chronic_dx (1=Yes/2=No haven_labelled) -> 1L/0L integer.
#   cancer_active, mental_health: already 0/1, just strip class.
#   total_exp, AGE22X: strip haven_labelled to keep model.matrix happy.
#   REGION22: convert to factor with census-region labels.
#
# -- NA handling ------------------------------------------------------------
#   FIT (model_inc): NA preserved; svyglm complete-cases automatically.
#     mental_health (~40% NA) drops rows from the fit, but the weight is
#     estimated where data exists.
#   SCORE: cond_mat NA -> 0. "Not diagnosed/observed" treated as
#     "condition absent" for purposes of summing weights.
#   The two populations differ — intentional.
# =============================================================================

suppressPackageStartupMessages({
  library(survey)
  library(splines)
  library(marginaleffects)
  library(dplyr)
  library(readr)
})

options(survey.lonely.psu = "adjust")

dir.create("output/tables",  recursive = TRUE, showWarnings = FALSE)
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

cat("Loading analytical dataset and survey design...\n")
meps        <- readRDS("data/processed/meps_analytical.rds")
meps_design <- readRDS("data/processed/meps_design.rds")
cat(sprintf("Loaded: %d rows.\n", nrow(meps)))

# ---- Strip haven_labelled and binarize chronic_dx ---------------------------
strip_attrs <- function(x) { attributes(x) <- NULL; x }

# 1=Yes / 2=No haven_labelled -> 1L / 0L integer
for (v in c("diab_dx","hibp_dx","chd_dx","strk_dx","asth_dx","arth_dx")) {
  x <- strip_attrs(meps[[v]])
  meps[[v]] <- ifelse(is.na(x), NA_integer_,
                      ifelse(x == 1, 1L, 0L))
}

# Already 0/1 -- just strip class/labels
for (v in c("cancer_active","mental_health","any_inpatient_2022")) {
  meps[[v]] <- as.integer(strip_attrs(meps[[v]]))
}

# Strip raw numerics
for (v in c("total_exp","AGE22X")) {
  meps[[v]] <- strip_attrs(meps[[v]])
}

# REGION22: 4-category census region as factor
meps$REGION22 <- factor(strip_attrs(meps$REGION22),
                        levels = 1:4,
                        labels = c("Northeast","Midwest","South","West"))

# Mirror cleaned columns onto design$variables
mirror_vars <- c("diab_dx","hibp_dx","chd_dx","strk_dx","asth_dx","arth_dx",
                 "cancer_active","mental_health","any_inpatient_2022",
                 "total_exp","AGE22X","REGION22")
for (v in mirror_vars) meps_design$variables[[v]] <- meps[[v]]

# ---- 1. Margen intensivo ----------------------------------------------------
design_pos <- subset(meps_design, total_exp > 0)

# ---- 2. Modelo Gamma log-link (any_inpatient_2022 as CONTROL) ---------------
model_inc <- svyglm(
  total_exp ~ diab_dx + hibp_dx + chd_dx + strk_dx + asth_dx +
              arth_dx + cancer_active + mental_health +
              any_inpatient_2022 +
              splines::ns(AGE22X, df = 4) + sex +
              income_pct_fpl + insurance_status + REGION22 +
              health_status_3 + any_functional_limitation,
  design = design_pos,
  family = Gamma(link = "log")
)

cat("\n=== Model coefficients (full table) ===\n")
print(round(summary(model_inc)$coefficients, 4))

# Coefficient highlight: any_inpatient_2022 (control)
inp_row <- summary(model_inc)$coefficients["any_inpatient_2022", , drop = FALSE]
cat(sprintf("\n[any_inpatient_2022]  beta=%.3f  SE=%.3f  p=%.4f\n",
            inp_row[1, "Estimate"], inp_row[1, "Std. Error"],
            inp_row[1, "Pr(>|t|)"]))

n_eff <- length(model_inc$residuals)
n_pos <- sum(meps$total_exp > 0, na.rm = TRUE)
cat(sprintf("\nEffective N in fit: %d (of %d positive-exp adults; %.1f%% retained)\n",
            n_eff, n_pos, 100 * n_eff / n_pos))
if (n_eff < 10000) {
  warning(sprintf("Effective N=%d < 10,000. Investigate dropped rows.", n_eff))
}

max_se <- max(summary(model_inc)$coefficients[, "Std. Error"], na.rm = TRUE)
cat(sprintf("Max coefficient SE: %.3f\n", max_se))
if (max_se > 5) {
  warning(sprintf("Max SE=%.3f > 5: severe multicollinearity probable.", max_se))
}

# ---- 3. AMEs en USD (only for the 8 score components) -----------------------
inc_eff <- avg_slopes(
  model_inc,
  variables = c("diab_dx","hibp_dx","chd_dx","strk_dx",
                "asth_dx","arth_dx","cancer_active","mental_health"),
  type = "response"
)

w_raw      <- as.numeric(inc_eff$estimate)
cond_names <- inc_eff$term

cat("\n=== Raw AMEs (USD) ===\n")
print(round(setNames(w_raw, cond_names), 2))

# 3b. AMEs no positivos
if (any(w_raw <= 0)) {
  warning(sprintf("AME <= 0 in: %s. Setting to 0.1 * min positive AME.",
                  paste(cond_names[w_raw <= 0], collapse = ", ")))
  min_pos <- min(w_raw[w_raw > 0])
  w_raw[w_raw <= 0] <- 0.1 * min_pos
}

w_prop <- w_raw / min(w_raw)
names(w_prop) <- cond_names

# 3c. Scale to [0.5, 3.0] for numerical stability
w_scaled <- 0.5 + (w_prop - min(w_prop)) /
                  (max(w_prop) - min(w_prop)) * 2.5
names(w_scaled) <- names(w_prop)

cat("\n=== Proportional weights (w_prop) and scaled to [0.5, 3.0] (w_scaled) ===\n")
print(data.frame(condition = names(w_prop),
                 w_prop    = round(w_prop, 2),
                 w_scaled  = round(w_scaled, 3),
                 row.names = NULL))

# ---- 4. Matriz de condiciones (NA -> 0 para SCORE) --------------------------
cond_cols <- c("diab_dx","hibp_dx","chd_dx","strk_dx",
               "asth_dx","arth_dx","cancer_active","mental_health")
cond_mat <- as.matrix(meps[, cond_cols])
cond_mat[is.na(cond_mat)] <- 0

# ---- 5. Backup + new score (using w_scaled) --------------------------------
meps$morbidity_score_flat <- meps$morbidity_score
meps$morbidity_score      <- as.numeric(cond_mat %*% w_scaled[cond_cols])
meps$morbidity_score_std  <- as.numeric(scale(meps$morbidity_score))

# ---- 5b. Sensitivity: score without cancer_active ---------------------------
cond_cols_no_cancer  <- setdiff(cond_cols, "cancer_active")
w_no_cancer          <- w_scaled[cond_cols_no_cancer]
cond_mat_no_cancer   <- as.matrix(meps[, cond_cols_no_cancer])
cond_mat_no_cancer[is.na(cond_mat_no_cancer)] <- 0
meps$morbidity_score_no_cancer     <- as.numeric(cond_mat_no_cancer %*% w_no_cancer)
meps$morbidity_score_no_cancer_std <- as.numeric(scale(meps$morbidity_score_no_cancer))

# ---- 6. Diagnóstico ---------------------------------------------------------
old_sd   <- sd(meps$morbidity_score_flat,        na.rm = TRUE)
new_sd   <- sd(meps$morbidity_score,             na.rm = TRUE)
nc_sd    <- sd(meps$morbidity_score_no_cancer,   na.rm = TRUE)
corr_val <- cor(meps$morbidity_score_flat, meps$morbidity_score,
                use = "complete.obs")
corr_nc  <- cor(meps$morbidity_score, meps$morbidity_score_no_cancer,
                use = "complete.obs")

cat("\n=== VALIDATION ===\n")
cat("SD score original (cost-weighted log1p) :", round(old_sd, 3), "\n")
cat("SD score nuevo  (AME w_scaled)          :", round(new_sd, 3), "\n")
cat("SD score sensitivity (no cancer)        :", round(nc_sd, 3), "\n")
cat("Correlación old vs new                  :", round(corr_val, 3), "\n")
cat("Correlación new vs no-cancer            :", round(corr_nc, 3), "\n")
cat("Rango pesos w_scaled                    :",
    paste(round(range(w_scaled), 2), collapse = " - "), "\n")
cat("Ratio max/min de w_scaled               :",
    round(max(w_scaled)/min(w_scaled), 2), "\n")
cat("Orden pesos w_scaled (desc)             :",
    paste(names(sort(w_scaled, decreasing = TRUE)), collapse = ", "), "\n")

# ---- 7. Persistencia --------------------------------------------------------
saveRDS(meps, "data/processed/meps_analytical.rds")

meps_design <- svydesign(
  id      = ~VARPSU,
  strata  = ~VARSTR,
  weights = ~PERWT22F,
  data    = meps,
  nest    = TRUE
)
saveRDS(meps_design, "data/processed/meps_design.rds")

# ---- 8. Tabla pesos ---------------------------------------------------------
weights_tbl <- tibble(
  condition   = cond_names,
  ame_usd     = round(w_raw, 2),
  w_prop      = round(w_prop, 2),
  w_scaled    = round(w_scaled, 3)
)
write_csv(weights_tbl, "output/tables/morbidity_cost_weights.csv")
cat("\nSaved: output/tables/morbidity_cost_weights.csv\n")

# ---- 9. Criterios de paso a Producto 4 -------------------------------------
cat("\n=== CRITERIOS DE PASO A PRODUCTO 4 ===\n")
cat(sprintf("  new_sd > 0.8                : %s (observed: %.3f)\n",
            new_sd > 0.8, new_sd))
cat(sprintf("  corr en [0.6, 0.95]         : %s (observed: %.3f)\n",
            corr_val >= 0.6 & corr_val <= 0.95, corr_val))
cat(sprintf("  ratio max/min > 3           : %s (observed: %.2f)\n",
            max(w_scaled)/min(w_scaled) > 3,
            max(w_scaled)/min(w_scaled)))
top3 <- names(sort(w_scaled, decreasing = TRUE))[1:3]
cat(sprintf("  chd_dx y strk_dx en top 3   : %s (top3: %s)\n",
            all(c("chd_dx","strk_dx") %in% top3),
            paste(top3, collapse = ", ")))

cat("\nStep 02b (AME-weighted morbidity score, rev 2) complete.\n")
