# 05_models_glm.R
# MEPS High-Cost Patient Prediction Project
# Step 5: Gamma GLM on Positive Expenditures
#
# Rationale:
#   - Restricted to total_exp > 0 (zero-expenditure excluded by design;
#     zeros are structurally different and handled by the logistic models).
#   - log(total_exp + 1) is NOT used: it biases estimates and is indefensible
#     when zeros are already excluded.
#   - Gamma with log link is standard for right-skewed positive expenditure data.
#   - AMEs are reported in dollars (not log-scale), computed via marginaleffects.
#
# Outputs:
#   output/models/glm_models.rds
#   output/tables/glm_sample_size.csv
#   output/tables/ame_glm_GB_main.csv
#   output/tables/ame_glm_GC_interaction.csv
#   output/tables/deviance_r2_glm.csv

if (!require("survey"))          install.packages("survey")
if (!require("dplyr"))           install.packages("dplyr")
if (!require("marginaleffects")) install.packages("marginaleffects")

library(survey)
library(dplyr)
library(marginaleffects)

dir.create("output/models", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)

cat("Loading survey design object...\n")
design_full <- readRDS("data/processed/meps_design.rds")

# -----------------------------------------------------------------------------
# 1. Restrict to positive expenditures
#    This is the analytical sample for the GLM stage.
#    Zeros are excluded by design — not imputed, not log(x+1).
# -----------------------------------------------------------------------------
cat("Restricting to total_exp > 0...\n")

design_pos <- subset(design_full, total_exp > 0)

n_full <- nrow(design_full$variables)
n_pos  <- nrow(design_pos$variables)
cat(sprintf("Full sample:              %d\n", n_full))
cat(sprintf("Positive-expenditure:     %d (%.1f%%)\n",
            n_pos, n_pos / n_full * 100))

# -----------------------------------------------------------------------------
# 2. Ensure mental_health_cat and factor references are set
# -----------------------------------------------------------------------------
if (!"mental_health_cat" %in% names(design_pos$variables)) {
  design_pos <- update(design_pos,
    mental_health_cat = factor(
      ifelse(is.na(design_pos$variables$mental_health), "Not observed",
             ifelse(design_pos$variables$mental_health == 1, "SPD", "No SPD")),
      levels = c("No SPD", "SPD", "Not observed")
    )
  )
}

design_pos <- update(design_pos,
  age_cat                  = relevel(factor(design_pos$variables$age_cat),
                                     ref = "18-34"),
  insurance_status         = relevel(factor(design_pos$variables$insurance_status),
                                     ref = "Private"),
  income_pct_fpl           = relevel(factor(design_pos$variables$income_pct_fpl),
                                     ref = "High Income"),
  n_chronic_conditions_cat = relevel(factor(design_pos$variables$n_chronic_conditions_cat),
                                     ref = "0"),
  mental_health_cat        = relevel(factor(design_pos$variables$mental_health_cat),
                                     ref = "No SPD")
)

cat("Factor references set.\n")

# -----------------------------------------------------------------------------
# 3. Sparsity check for interaction model
# -----------------------------------------------------------------------------
cat("\nCross-tab: n_chronic_conditions_cat x insurance_status (positive-exp sample)\n")
xtab_pos <- table(
  design_pos$variables$n_chronic_conditions_cat,
  design_pos$variables$insurance_status,
  useNA = "ifany"
)
print(xtab_pos)

sparse_cells <- sum(xtab_pos < 30, na.rm = TRUE)
if (sparse_cells > 0) {
  cat(sprintf("WARNING: %d cells with n < 30. Interaction model may be unstable.\n",
              sparse_cells))
} else {
  cat("All cells >= 30. Interaction model should be stable.\n")
}

# -----------------------------------------------------------------------------
# 4. Model formulas (same structure as logistic, outcome = total_exp)
# -----------------------------------------------------------------------------

# Model G-A: Demographics only
form_GA <- total_exp ~ age_cat + sex + income_pct_fpl + insurance_status +
                       education + race_ethnicity

# Model G-B: Demographics + Morbidity (main results)
form_GB <- total_exp ~ age_cat + sex + income_pct_fpl + insurance_status +
                       education + race_ethnicity +
                       n_chronic_conditions_cat + cancer_ever + mental_health_cat

# Model G-C: Morbidity x Insurance interaction (extension)
form_GC <- total_exp ~ age_cat + sex + income_pct_fpl + education + race_ethnicity +
                       n_chronic_conditions_cat * insurance_status +
                       cancer_ever + mental_health_cat

# -----------------------------------------------------------------------------
# 5. Fit Gamma GLM with log link
# -----------------------------------------------------------------------------
cat("\nFitting Gamma GLM models (log link)...\n")

mod_GA <- svyglm(form_GA, design = design_pos, family = Gamma(link = "log"))
cat(sprintf("G-A done | n = %d\n", nrow(model.frame(mod_GA))))

mod_GB <- svyglm(form_GB, design = design_pos, family = Gamma(link = "log"))
cat(sprintf("G-B done | n = %d\n", nrow(model.frame(mod_GB))))

mod_GC <- svyglm(form_GC, design = design_pos, family = Gamma(link = "log"))
cat(sprintf("G-C done | n = %d\n", nrow(model.frame(mod_GC))))

# Save models
saveRDS(
  list(GA = mod_GA, GB = mod_GB, GC = mod_GC),
  "output/models/glm_models.rds"
)
cat("GLM models saved.\n")

# Sample sizes
n_table_glm <- data.frame(
  model = c("G-A", "G-B", "G-C"),
  spec  = c("Demographics", "Demo+Morbidity", "Demo+Morbidity*Insurance"),
  n     = c(nrow(model.frame(mod_GA)),
            nrow(model.frame(mod_GB)),
            nrow(model.frame(mod_GC)))
)
write.csv(n_table_glm, "output/tables/glm_sample_size.csv", row.names = FALSE)
print(n_table_glm)

# -----------------------------------------------------------------------------
# 6. Convergence check
#    svyglm does not always warn on non-convergence; check manually
# -----------------------------------------------------------------------------
check_convergence <- function(model, name) {
  if (!model$converged) {
    cat(sprintf("WARNING: %s did not converge!\n", name))
  } else {
    cat(sprintf("%s converged OK.\n", name))
  }
  # Flag implausibly large coefficients (|coef| > 10 on log scale)
  big_coefs <- names(which(abs(coef(model)) > 10))
  if (length(big_coefs) > 0) {
    cat(sprintf("  Large coefficients in %s: %s\n", name, paste(big_coefs, collapse = ", ")))
  }
}

cat("\nConvergence checks:\n")
check_convergence(mod_GA, "G-A")
check_convergence(mod_GB, "G-B")
check_convergence(mod_GC, "G-C")

# -----------------------------------------------------------------------------
# 7. AMEs in dollars — Model G-B (main results)
#    avg_slopes() on Gamma log-link returns marginal effects on the response
#    scale (dollars), not on the log scale. This is what we want.
# -----------------------------------------------------------------------------
cat("\nComputing AMEs in dollars for G-B...\n")

ame_GB <- avg_slopes(mod_GB)
ame_GB_df <- as.data.frame(ame_GB)[, c("term","contrast","estimate","std.error",
                                        "statistic","p.value","conf.low","conf.high")]

write.csv(ame_GB_df, "output/tables/ame_glm_GB_main.csv", row.names = FALSE)
cat("\nAMEs (dollars) — Model G-B:\n")
print(ame_GB_df)

# -----------------------------------------------------------------------------
# 8. Interaction AMEs — Model G-C
#    AME of morbidity on expenditure by insurance type (in dollars)
# -----------------------------------------------------------------------------
cat("\nComputing interaction AMEs for G-C...\n")

ame_GC <- avg_slopes(
  mod_GC,
  variables = "n_chronic_conditions_cat",
  by        = "insurance_status"
)

ame_GC_df <- as.data.frame(ame_GC)
write.csv(ame_GC_df, "output/tables/ame_glm_GC_interaction.csv", row.names = FALSE)

cat("\nInteraction AMEs (dollars) — Model G-C:\n")
print(ame_GC_df[, intersect(c("term","contrast","insurance_status",
                               "estimate","std.error","p.value"), names(ame_GC_df))])

# -----------------------------------------------------------------------------
# 9. Deviance-based pseudo-R2 (incremental comparison G-A -> G-B -> G-C)
# -----------------------------------------------------------------------------
cat("\nComputing deviance-based pseudo-R2 (GLM)...\n")

calc_dev_r2 <- function(model) round(1 - (model$deviance / model$null.deviance), 4)

r2_glm <- data.frame(
  Model = c("G-A: Demographics only",
            "G-B: Demo + Morbidity",
            "G-C: Demo + Morbidity x Insurance"),
  DevR2 = c(calc_dev_r2(mod_GA), calc_dev_r2(mod_GB), calc_dev_r2(mod_GC))
)

write.csv(r2_glm, "output/tables/deviance_r2_glm.csv", row.names = FALSE)
cat("\nDeviance-based pseudo-R2 (GLM, incremental interpretation only):\n")
print(r2_glm)

cat("\nStep 5 complete.\n")
