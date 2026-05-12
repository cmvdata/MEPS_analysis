# R/04_main_models.R
# ---------------------------------------------------------------
# Producto 4: Insurance × Morbidity interaction
#
# Tests whether insurance type modulates the morbidity-cost
# gradient. Two-part-style: this script handles intensive margin
# (Gamma log-link on positive expenditure). Extensive margin in
# Producto 5 (sensitivity).
#
# Primary specification: morbidity_score (AME-weighted, Producto 3 rev2).
# Sensitivity: morbidity_score_flat (cost-mean log1p, rev1).
# ---------------------------------------------------------------

suppressPackageStartupMessages({
  library(survey)
  library(splines)
  library(marginaleffects)
  library(car)
  library(emmeans)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tibble)
})

# ---------------------------------------------------------------
# 0. Load and verify
# ---------------------------------------------------------------
meps <- readRDS("data/processed/meps_analytical.rds")
meps_design <- readRDS("data/processed/meps_design.rds")

required_vars <- c("total_exp", "insurance_status", "morbidity_score",
                   "morbidity_score_flat", "AGE22X", "sex", "POVCAT22",
                   "REGION22", "health_status_3",
                   "any_functional_limitation", "any_inpatient_2022")

missing_vars <- setdiff(required_vars, names(meps))
if (length(missing_vars) > 0) {
  stop(sprintf("Variables faltantes: %s", paste(missing_vars, collapse = ", ")))
}

cat("=== PRODUCTO 4 - Pre-fit verification ===\n")
cat(sprintf("Total rows: %d\n", nrow(meps)))
cat(sprintf("Rows with total_exp > 0: %d\n",
            sum(meps$total_exp > 0, na.rm = TRUE)))

cat("\nInsurance levels (table):\n")
print(table(meps$insurance_status, useNA = "ifany"))

ins_counts <- table(meps$insurance_status)
if (any(ins_counts < 100)) {
  warning(sprintf("Categoria con N<100: %s",
                  paste(names(ins_counts)[ins_counts < 100],
                        collapse = ", ")))
}

# Standardize morbidity scores
meps$morbidity_score_std <- as.numeric(scale(meps$morbidity_score))
meps$morbidity_flat_std <- as.numeric(scale(meps$morbidity_score_flat))

cat("\nmorbidity_score_std summary:\n")
print(summary(meps$morbidity_score_std))

# Regenerate survey design with standardized vars
meps_design <- svydesign(
  id = ~VARPSU,
  strata = ~VARSTR,
  weights = ~PERWT22F,
  data = meps,
  nest = TRUE
)

# ---------------------------------------------------------------
# 1. Main model: intensive margin (Gamma log-link, total_exp > 0)
# ---------------------------------------------------------------
design_pos <- subset(meps_design, total_exp > 0)

cat("\n=== Fitting main model (Gamma, AME-weighted score) ===\n")

formula_main <- total_exp ~ insurance_status * morbidity_score_std +
                splines::ns(AGE22X, df = 4) + sex + POVCAT22 + REGION22 +
                health_status_3 + any_functional_limitation +
                any_inpatient_2022

model_main <- svyglm(formula_main, design = design_pos,
                     family = Gamma(link = "log"))

n_main <- nrow(model.frame(model_main))
cat(sprintf("N efectivo (main): %d\n", n_main))

# ---------------------------------------------------------------
# 2. Sensitivity model: flat (cost-mean log1p) score
# ---------------------------------------------------------------
cat("\n=== Fitting sensitivity model (Gamma, flat score) ===\n")

formula_flat <- total_exp ~ insurance_status * morbidity_flat_std +
                splines::ns(AGE22X, df = 4) + sex + POVCAT22 + REGION22 +
                health_status_3 + any_functional_limitation +
                any_inpatient_2022

model_flat <- svyglm(formula_flat, design = design_pos,
                     family = Gamma(link = "log"))

n_flat <- nrow(model.frame(model_flat))
cat(sprintf("N efectivo (sensitivity): %d\n", n_flat))

# ---------------------------------------------------------------
# 3. Marginal effects of morbidity by insurance
# ---------------------------------------------------------------
cat("\n=== Computing AMEs by insurance (population-weighted) ===\n")

# Population-weighted AMEs: pass survey weights so avg_slopes averages
# unit-level slopes weighted by PERWT22F (otherwise the average is over
# the sample, not the survey-target population).
# Use model$prior.weights (length = n_used after NA drops) instead of
# weights(design_pos) (length = full positive subset incl. dropped NAs).
wts_vec <- model_main$prior.weights
stopifnot(length(wts_vec) == n_main)

mfx_main <- avg_slopes(
  model_main,
  variables  = "morbidity_score_std",
  by         = "insurance_status",
  type       = "response",
  conf_level = 0.95,
  wts        = wts_vec
)

mfx_flat <- avg_slopes(
  model_flat,
  variables  = "morbidity_flat_std",
  by         = "insurance_status",
  type       = "response",
  conf_level = 0.95,
  wts        = wts_vec
)

# ---------------------------------------------------------------
# 3b. Pairwise tests across insurance AMEs (audit fix)
# ---------------------------------------------------------------
# Operates on the population-weighted slopes object (mfx_main) so the
# pairwise differences carry the same weighting as the headline AMEs.
# Wald p-values, no multiplicity adjustment (companion to the emmeans
# Tukey-adjusted pairs in section 5 below).
cat("\n=== Pairwise AME comparisons (Wald, audit fix) ===\n")

pairwise_tests <- hypotheses(
  mfx_main,
  hypothesis = ~pairwise
)

print(pairwise_tests)

# Relabel b1..b4 -> insurance categories (order from mfx_main rows)
ins_labels <- as.character(as.data.frame(mfx_main)$insurance_status)
pairwise_df <- as.data.frame(pairwise_tests)
pairwise_df$comparison <- vapply(pairwise_df$hypothesis, function(h) {
  m <- regmatches(h, regexec("\\(b(\\d+)\\) - \\(b(\\d+)\\)", h))[[1]]
  if (length(m) == 3L) {
    sprintf("%s - %s", ins_labels[as.integer(m[2])],
                       ins_labels[as.integer(m[3])])
  } else h
}, character(1))
pairwise_df <- pairwise_df[, c("comparison", "hypothesis", "estimate",
                                "std.error", "statistic", "p.value",
                                "conf.low", "conf.high")]
write.csv(
  pairwise_df,
  "output/tables/AME_pairwise_tests.csv",
  row.names = FALSE
)
cat("\nPairwise table (relabelled):\n")
print(pairwise_df)

# ---------------------------------------------------------------
# 4. Wald test: are the slopes different across insurance?
# ---------------------------------------------------------------
cat("\n=== Wald test (interaction terms) ===\n")

interaction_terms <- grep("insurance_status.*morbidity_score_std",
                          names(coef(model_main)), value = TRUE)
cat("Interaction terms tested:\n")
print(interaction_terms)

wald_main <- linearHypothesis(model_main,
                              paste(interaction_terms, "= 0"))
print(wald_main)

# ---------------------------------------------------------------
# 5. Pairwise contrasts (emmeans, all insurance pairs)
# ---------------------------------------------------------------
cat("\n=== Pairwise slope comparisons ===\n")

# Try emmeans first; fall back to marginaleffects::avg_slopes(hypothesis="pairwise")
emm_trends <- tryCatch({
  emtrends(model_main, ~ insurance_status, var = "morbidity_score_std")
}, error = function(e) {
  cat("emtrends failed:", e$message, "\n")
  NULL
})

if (!is.null(emm_trends)) {
  cat("\nSlopes per insurance (emmeans):\n")
  print(summary(emm_trends, infer = TRUE))
  cat("\nPairwise differences (Tukey):\n")
  print(pairs(emm_trends, adjust = "tukey"))
} else {
  cat("\nFallback: marginaleffects pairwise via formula syntax\n")
  pairwise_main <- tryCatch(
    avg_slopes(
      model_main,
      variables  = "morbidity_score_std",
      by         = "insurance_status",
      type       = "response",
      hypothesis = ~pairwise,
      wts        = wts_vec
    ),
    error = function(e) {
      cat("Pairwise fallback also failed:", e$message, "\n")
      cat("Reporting per-group AMEs only (already in mfx_main).\n")
      NULL
    }
  )
  if (!is.null(pairwise_main)) print(as.data.frame(pairwise_main))
}

# ---------------------------------------------------------------
# 6. Visualization
# ---------------------------------------------------------------
plot_data <- bind_rows(
  as.data.frame(mfx_main) %>% mutate(score = "AME-weighted (primary)"),
  as.data.frame(mfx_flat) %>% mutate(score = "Cost-mean log1p (sensitivity)")
)

p <- ggplot(plot_data,
            aes(x = insurance_status, y = estimate,
                ymin = conf.low, ymax = conf.high,
                color = score)) +
  geom_pointrange(position = position_dodge(width = 0.4), size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  scale_color_manual(values = c("AME-weighted (primary)" = "#2c3e50",
                                "Cost-mean log1p (sensitivity)" = "#e74c3c")) +
  labs(
    title = "Marginal effect of morbidity on annual healthcare spending",
    subtitle = "By insurance type, USA 2022 (1-SD increase in morbidity score)",
    x = "Insurance category",
    y = expression(paste(Delta, " Annual Spending (USD)")),
    color = "Score version",
    caption = sprintf("MEPS 2022 (HC-243). N effective = %d. Gamma GLM, log link, intensive margin.\nControls: age (spline df=4), sex, FPL, region, self-rated health, functional limitation, inpatient 2022.", n_main)
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        plot.caption = element_text(size = 8))

dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)
ggsave("output/figures/mfx_insurance_interaction.png", p,
       width = 9, height = 6, dpi = 300, bg = "white")
cat("\nFigura guardada: output/figures/mfx_insurance_interaction.png\n")

# ---------------------------------------------------------------
# 7. Save tables
# ---------------------------------------------------------------
dir.create("output/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("output/models", showWarnings = FALSE, recursive = TRUE)

mfx_main_tbl <- as.data.frame(mfx_main) %>%
  select(insurance_status, estimate, std.error, conf.low, conf.high, p.value)
write_csv(mfx_main_tbl, "output/tables/AME_morbidity_by_insurance_main.csv")

mfx_flat_tbl <- as.data.frame(mfx_flat) %>%
  select(insurance_status, estimate, std.error, conf.low, conf.high, p.value)
write_csv(mfx_flat_tbl, "output/tables/AME_morbidity_by_insurance_sensitivity.csv")

# linearHypothesis on svyglm returns Chisq / Pr(>Chisq) (not F)
wald_tbl <- tibble(
  test            = "H0: all interaction terms = 0",
  chisq_statistic = wald_main$Chisq[2],
  df              = wald_main$Df[2],
  res_df          = wald_main$Res.Df[2],
  p_value         = wald_main$`Pr(>Chisq)`[2]
)
write_csv(wald_tbl, "output/tables/wald_test_modulation.csv")

saveRDS(list(
  model_main = model_main,
  model_flat = model_flat,
  mfx_main = mfx_main,
  mfx_flat = mfx_flat,
  wald_main = wald_main,
  emm_trends = emm_trends,
  n_main = n_main,
  n_flat = n_flat
), "output/models/p4_interaction_models.rds")

# ---------------------------------------------------------------
# 8. Diagnostic summary
# ---------------------------------------------------------------
cat("\n", strrep("=", 70), "\n", sep = "")
cat("PRODUCTO 4 - DIAGNOSTIC SUMMARY\n")
cat(strrep("=", 70), "\n", sep = "")

cat(sprintf("\nN effective (main):        %d\n", n_main))
cat(sprintf("N effective (sensitivity): %d\n", n_flat))
cat(sprintf("Total positive sample:     %d\n",
            sum(meps$total_exp > 0, na.rm = TRUE)))

retention <- 100 * n_main / sum(meps$total_exp > 0, na.rm = TRUE)
cat(sprintf("Retention rate:            %.1f%%\n", retention))

if (retention < 90) {
  cat("\n[NOTE] N effective < 90% del positive sample.\n")
  cat("Probable cause: mental_health NA en SAQ subsample (39.6% por diseno MEPS).\n")
  cat("Pesos del score se calibraron con complete-case; aplicados a muestra completa.\n")
}

# AIC.svyglm returns c(eff.p, AIC, deltabar); extract the AIC scalar
aic_main <- as.numeric(AIC(model_main)["AIC"])
aic_flat <- as.numeric(AIC(model_flat)["AIC"])
cat(sprintf("\nAIC main:        %.2f\n", aic_main))
cat(sprintf("AIC sensitivity: %.2f\n", aic_flat))
cat(sprintf("AIC delta (main - sens): %.2f\n", aic_main - aic_flat))

cat("\nDispersion (main):", summary(model_main)$dispersion, "\n")
cat("Dispersion (sens):", summary(model_flat)$dispersion, "\n")

# Wald test (svyglm linearHypothesis -> Chisq, not F)
chisq_stat <- wald_main$Chisq[2]
chisq_p    <- wald_main$`Pr(>Chisq)`[2]
cat("\nWald test on interaction terms:\n")
cat(sprintf("  Chisq = %.3f, df=%.0f, res_df=%.0f, p = %.5f\n",
            chisq_stat, wald_main$Df[2], wald_main$Res.Df[2], chisq_p))
if (chisq_p < 0.05) {
  cat("  -> SIGNIFICANT modulation by insurance type.\n")
} else {
  cat("  -> No significant modulation by insurance type.\n")
}

cat("\n=== AME by insurance (PRIMARY: AME-weighted score) ===\n")
print(mfx_main_tbl)

cat("\n=== AME by insurance (SENSITIVITY: flat score) ===\n")
print(mfx_flat_tbl)

cat("\n=== Robustness check ===\n")
sign_match <- sign(mfx_main_tbl$estimate) == sign(mfx_flat_tbl$estimate)
cat(sprintf("Signs match across specifications: %d / %d\n",
            sum(sign_match), length(sign_match)))
if (all(sign_match)) {
  cat("-> Robustness: signs of AMEs CONSISTENT.\n")
} else {
  cat("-> ALERT: signs DIFFER. Investigate before publication.\n")
}

cat("\n=== Outputs ===\n")
cat("  output/tables/AME_morbidity_by_insurance_main.csv\n")
cat("  output/tables/AME_morbidity_by_insurance_sensitivity.csv\n")
cat("  output/tables/wald_test_modulation.csv\n")
cat("  output/models/p4_interaction_models.rds\n")
cat("  output/figures/mfx_insurance_interaction.png\n")

cat("\nProducto 4 completado. Espera luz verde antes de Producto 5.\n")
