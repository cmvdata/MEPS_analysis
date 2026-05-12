# R/05_sensitivity.R
# ---------------------------------------------------------------
# Producto 5: Sensitivity and robustness tests
#
# Four orthogonal sensitivity dimensions (scope acotado):
#   A) Outcome definition: top5/top10 + binary thresholds $15k/$50k
#   B) Morbidity score variant: AME-weighted (primary) / cost-mean
#      log1p (sensitivity) / AME-weighted no cancer
#   C) Subpopulation: age <65 / no functional limit / SAQ subsample
#   D) Influential observations: exclude top 1% spending (>p99)
#
# All tests reuse the SAME RHS as model_main from 04_main_models.R.
# AMEs population-weighted via weights(design, "sampling"). For svyglm,
# mod$prior.weights == weights(design, "sampling")[rows_used] (svyglm
# passes 1/design$prob as glm's weights arg). When NA drops occur in
# covariates, length(weights(design)) > length(mod$prior.weights), and
# avg_slopes requires a vector matching the model's effective N -- the
# helper falls back to mod$prior.weights in that case (still sampling
# weights, just for the rows actually used).
# ---------------------------------------------------------------

suppressPackageStartupMessages({
  library(survey)
  library(splines)
  library(marginaleffects)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tibble)
})
options(survey.lonely.psu = "adjust")

dir.create("output/tables",  recursive = TRUE, showWarnings = FALSE)
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

# ---- 0. Load + standardize scores + construct alt outcomes ------------------
meps        <- readRDS("data/processed/meps_analytical.rds")
meps_design <- readRDS("data/processed/meps_design.rds")

meps$morbidity_score_std           <- as.numeric(scale(meps$morbidity_score))
meps$morbidity_flat_std            <- as.numeric(scale(meps$morbidity_score_flat))
meps$morbidity_score_no_cancer_std <- as.numeric(scale(meps$morbidity_score_no_cancer))

# Alt outcomes for test A
q90 <- svyquantile(~total_exp, meps_design, quantiles = 0.90, na.rm = TRUE)
thresh_90 <- as.numeric(coef(q90))
meps$top10        <- as.integer(!is.na(meps$total_exp) & meps$total_exp >= thresh_90)
meps$exp_over_15k <- as.integer(!is.na(meps$total_exp) & meps$total_exp > 15000)
meps$exp_over_50k <- as.integer(!is.na(meps$total_exp) & meps$total_exp > 50000)

cat(sprintf("Top10 threshold (p90): $%.0f\n", thresh_90))
cat(sprintf("Top10 share         : %.2f%%\n", 100*mean(meps$top10)))
cat(sprintf("> $15k share        : %.2f%%\n", 100*mean(meps$exp_over_15k)))
cat(sprintf("> $50k share        : %.2f%%\n", 100*mean(meps$exp_over_50k)))

# Regenerate design with new columns
meps_design <- svydesign(
  id = ~VARPSU, strata = ~VARSTR, weights = ~PERWT22F,
  data = meps, nest = TRUE
)
design_pos <- subset(meps_design, total_exp > 0)

# Common controls (same as model_main in 04)
control_terms <- paste(
  "splines::ns(AGE22X, df = 4)", "sex", "POVCAT22", "REGION22",
  "health_status_3", "any_functional_limitation",
  "any_inpatient_2022",
  sep = " + "
)

# ---- Helpers ---------------------------------------------------------------
# Get sampling weights matching the model's effective rows.
# For svyglm: weights(design, "sampling") == mod$prior.weights for shared rows.
# When length(weights(design)) != length(prior.weights) (NA drops in covariates),
# fall back to prior.weights (correct length for avg_slopes).
get_wts_vec <- function(mod, design) {
  n_design <- nrow(design$variables)
  n_used   <- length(mod$prior.weights)
  if (n_design == n_used) {
    weights(design, "sampling")
  } else {
    cat(sprintf("    [note] %d NA drops -> using mod$prior.weights (n_used=%d)\n",
                n_design - n_used, n_used))
    mod$prior.weights
  }
}

ame_by_insurance_gamma <- function(score_var, design, label) {
  fml <- as.formula(sprintf(
    "total_exp ~ insurance_status * %s + %s",
    score_var, control_terms
  ))
  mod <- svyglm(fml, design = design, family = Gamma(link = "log"))
  wts_vec <- get_wts_vec(mod, design)
  mfx <- avg_slopes(mod, variables = score_var, by = "insurance_status",
                    type = "response", wts = wts_vec)
  as.data.frame(mfx) %>%
    mutate(spec = label,
           n_eff = length(mod$prior.weights)) %>%
    select(spec, n_eff, insurance_status, estimate, std.error,
           conf.low, conf.high, p.value)
}

ame_by_insurance_logit <- function(outcome, design, label) {
  fml <- as.formula(sprintf(
    "%s ~ insurance_status * morbidity_score_std + %s",
    outcome, control_terms
  ))
  mod <- svyglm(fml, design = design, family = quasibinomial(link = "logit"))
  wts_vec <- get_wts_vec(mod, design)
  mfx <- avg_slopes(mod, variables = "morbidity_score_std",
                    by = "insurance_status", type = "response", wts = wts_vec)
  as.data.frame(mfx) %>%
    mutate(spec = label,
           n_eff = length(mod$prior.weights)) %>%
    select(spec, n_eff, insurance_status, estimate, std.error,
           conf.low, conf.high, p.value)
}

# ---- A) Outcome robustness (logistic on full sample) -----------------------
cat("\n=== Test A: outcome robustness (logistic, prob scale) ===\n")
A_results <- bind_rows(
  ame_by_insurance_logit("top5",         meps_design, "top5 (p95)"),
  ame_by_insurance_logit("top10",        meps_design, "top10 (p90)"),
  ame_by_insurance_logit("exp_over_15k", meps_design, "> $15,000"),
  ame_by_insurance_logit("exp_over_50k", meps_design, "> $50,000")
)
write_csv(A_results, "output/tables/sensitivity_outcomes.csv")
print(A_results, row.names = FALSE)

# ---- B) Score robustness (Gamma on positive subset) ------------------------
cat("\n=== Test B: score robustness (Gamma USD scale) ===\n")
B_results <- bind_rows(
  ame_by_insurance_gamma("morbidity_score_std",           design_pos, "AME-weighted (primary)"),
  ame_by_insurance_gamma("morbidity_flat_std",            design_pos, "Cost-mean log1p"),
  ame_by_insurance_gamma("morbidity_score_no_cancer_std", design_pos, "AME-weighted no cancer")
)
write_csv(B_results, "output/tables/sensitivity_scores.csv")
print(B_results, row.names = FALSE)

# ---- C) Subpopulation robustness -------------------------------------------
cat("\n=== Test C: subpopulation robustness (Gamma USD scale) ===\n")
design_under65 <- subset(meps_design, AGE22X < 65 & total_exp > 0)
design_no_func <- subset(meps_design, any_functional_limitation == 0 & total_exp > 0)
design_saq     <- subset(meps_design, !is.na(mental_health) & total_exp > 0)

C_results <- bind_rows(
  ame_by_insurance_gamma("morbidity_score_std", design_under65, "Adults <65"),
  ame_by_insurance_gamma("morbidity_score_std", design_no_func, "No functional limit"),
  ame_by_insurance_gamma("morbidity_score_std", design_saq,     "SAQ subsample")
)
write_csv(C_results, "output/tables/sensitivity_subgroups.csv")
print(C_results, row.names = FALSE)

# ---- D) Influential observations (exclude p99+) ----------------------------
cat("\n=== Test D: influence -- exclude top 1% spending ===\n")
q99 <- svyquantile(~total_exp, meps_design, quantiles = 0.99, na.rm = TRUE)
thresh_99 <- as.numeric(coef(q99))
cat(sprintf("p99 threshold (excluded above): $%.0f\n", thresh_99))

design_no_extremes <- subset(meps_design, total_exp > 0 & total_exp < thresh_99)

D_results <- bind_rows(
  ame_by_insurance_gamma("morbidity_score_std", design_pos,         "Reference (full positive)"),
  ame_by_insurance_gamma("morbidity_score_std", design_no_extremes, "Excluding p99+")
)
write_csv(D_results, "output/tables/sensitivity_extremes.csv")
print(D_results, row.names = FALSE)

# Delta % vs reference (sanity for influence)
ref_d  <- D_results %>% filter(spec == "Reference (full positive)") %>%
          select(insurance_status, ref_est = estimate)
test_d <- D_results %>% filter(spec == "Excluding p99+") %>%
          select(insurance_status, test_est = estimate)
deltas <- left_join(ref_d, test_d, by = "insurance_status") %>%
  mutate(delta_pct = round(100 * (test_est - ref_est) / ref_est, 1))
cat("\nDelta % AMEs (excluding p99+ vs reference):\n")
print(deltas, row.names = FALSE)

# ---- Pattern check: Medicaid > Private > Uninsured -------------------------
cat("\n=== Pattern check: Medicaid > Private > Uninsured ===\n")
all_results <- bind_rows(
  A_results %>% mutate(test = "A: Outcome"),
  B_results %>% mutate(test = "B: Score"),
  C_results %>% mutate(test = "C: Subpop"),
  D_results %>% mutate(test = "D: Extremes")
)

pattern_check <- all_results %>%
  group_by(test, spec) %>%
  summarize(
    medicaid    = estimate[insurance_status == "Medicaid/Public Only"],
    private     = estimate[insurance_status == "Private"],
    uninsured   = estimate[insurance_status == "Uninsured"],
    pattern_ok  = (medicaid > private) & (private > uninsured),
    .groups = "drop"
  )
print(pattern_check, n = Inf)

n_match <- sum(pattern_check$pattern_ok, na.rm = TRUE)
n_total <- nrow(pattern_check)
cat(sprintf("\nPattern Medicaid>Private>Uninsured holds in %d / %d fits.\n",
            n_match, n_total))

# ---- Consolidated forest plot ----------------------------------------------
all_results <- all_results %>%
  mutate(panel = test)

p_summary <- ggplot(all_results,
                    aes(x = estimate, y = spec,
                        xmin = conf.low, xmax = conf.high,
                        color = insurance_status)) +
  geom_pointrange(position = position_dodge(width = 0.5), size = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ panel, scales = "free", ncol = 2) +
  scale_color_brewer(palette = "Set1") +
  labs(
    title    = "Sensitivity analyses: AME of morbidity by insurance",
    subtitle = "Test A: probability scale (logit). Tests B/C/D: USD (Gamma log-link).",
    x = "AME estimate (CI 95%)",
    y = NULL,
    color = "Insurance"
  ) +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold"))

ggsave("output/figures/sensitivity_summary.png", p_summary,
       width = 12, height = 9, dpi = 200, bg = "white")
cat("\nFigura guardada: output/figures/sensitivity_summary.png\n")

cat("\nProducto 5 completado.\n")
