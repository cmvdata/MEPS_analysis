# 04_models_logistic.R
# MEPS High-Cost Patient Prediction Project
# Step 4: Logistic Regression Models (Top 5% and Top 1%)
#
# Outputs:
#   output/models/logistic_models.rds
#   output/models/ame_main.rds
#   output/models/ame_interaction.rds
#   output/tables/model_sample_sizes.csv
#   output/tables/ame_B1_top5.csv
#   output/tables/ame_B2_top1.csv
#   output/tables/ame_C1_interaction.csv
#   output/tables/ame_C2_interaction.csv
#   output/tables/deviance_r2.csv

if (!require("survey"))          install.packages("survey")
if (!require("marginaleffects")) install.packages("marginaleffects")

library(survey)
library(marginaleffects)

dir.create("output/models", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)

cat("Loading survey design object...\n")
design <- readRDS("data/processed/meps_design.rds")

# Ensure mental_health_cat exists
if (!"mental_health_cat" %in% names(design$variables)) {
  cat("Building mental_health_cat...\n")
  design <- update(design,
    mental_health_cat = factor(
      ifelse(is.na(design$variables$mental_health), "Not observed",
             ifelse(design$variables$mental_health == 1, "SPD", "No SPD")),
      levels = c("No SPD", "SPD", "Not observed")
    )
  )
}

# -----------------------------------------------------------------------------
# 1. Set factor reference levels explicitly
#    Choosing the largest / most natural reference for each predictor
# -----------------------------------------------------------------------------
design <- update(design,
  age_cat                = relevel(factor(design$variables$age_cat),
                                   ref = "18-34"),
  insurance_status       = relevel(factor(design$variables$insurance_status),
                                   ref = "Private"),
  income_pct_fpl         = relevel(factor(design$variables$income_pct_fpl),
                                   ref = "High Income"),
  n_chronic_conditions_cat = relevel(factor(design$variables$n_chronic_conditions_cat),
                                     ref = "0"),
  mental_health_cat      = relevel(factor(design$variables$mental_health_cat),
                                   ref = "No SPD")
)

cat("Factor reference levels set.\n")

# -----------------------------------------------------------------------------
# 2. Sparsity check for Set C interaction
#    If any cell < 30, flag it before running the interaction model
# -----------------------------------------------------------------------------
cat("\nCross-tab: n_chronic_conditions_cat x insurance_status\n")
xtab <- table(
  design$variables$n_chronic_conditions_cat,
  design$variables$insurance_status,
  useNA = "ifany"
)
print(xtab)

sparse_cells <- sum(xtab < 30, na.rm = TRUE)
if (sparse_cells > 0) {
  cat(sprintf("WARNING: %d cells with n < 30. Interaction model may be unstable.\n",
              sparse_cells))
} else {
  cat("All cells >= 30. Interaction model should be stable.\n")
}

# -----------------------------------------------------------------------------
# 3. Model formulas
# -----------------------------------------------------------------------------

# Set A: Demographics only
form_A1 <- top5 ~ age_cat + sex + income_pct_fpl + insurance_status + education + race_ethnicity
form_A2 <- top1 ~ age_cat + sex + income_pct_fpl + insurance_status + education + race_ethnicity

# Set B: Demographics + Morbidity (main results)
form_B1 <- top5 ~ age_cat + sex + income_pct_fpl + insurance_status + education + race_ethnicity +
                  n_chronic_conditions_cat + cancer_ever + mental_health_cat
form_B2 <- top1 ~ age_cat + sex + income_pct_fpl + insurance_status + education + race_ethnicity +
                  n_chronic_conditions_cat + cancer_ever + mental_health_cat

# Set C: Morbidity x Insurance interaction (extension)
form_C1 <- top5 ~ age_cat + sex + income_pct_fpl + education + race_ethnicity +
                  n_chronic_conditions_cat * insurance_status + cancer_ever + mental_health_cat
form_C2 <- top1 ~ age_cat + sex + income_pct_fpl + education + race_ethnicity +
                  n_chronic_conditions_cat * insurance_status + cancer_ever + mental_health_cat

# -----------------------------------------------------------------------------
# 4. Fit models
# -----------------------------------------------------------------------------
cat("\nFitting models (quasibinomial)...\n")

mod_A1 <- svyglm(form_A1, design = design, family = quasibinomial())
cat(sprintf("A1 done  | n = %d\n", nrow(model.frame(mod_A1))))

mod_A2 <- svyglm(form_A2, design = design, family = quasibinomial())
cat(sprintf("A2 done  | n = %d\n", nrow(model.frame(mod_A2))))

mod_B1 <- svyglm(form_B1, design = design, family = quasibinomial())
cat(sprintf("B1 done  | n = %d\n", nrow(model.frame(mod_B1))))

mod_B2 <- svyglm(form_B2, design = design, family = quasibinomial())
cat(sprintf("B2 done  | n = %d\n", nrow(model.frame(mod_B2))))

mod_C1 <- svyglm(form_C1, design = design, family = quasibinomial())
cat(sprintf("C1 done  | n = %d\n", nrow(model.frame(mod_C1))))

mod_C2 <- svyglm(form_C2, design = design, family = quasibinomial())
cat(sprintf("C2 done  | n = %d\n", nrow(model.frame(mod_C2))))

# Save models
saveRDS(
  list(A1 = mod_A1, A2 = mod_A2, B1 = mod_B1, B2 = mod_B2, C1 = mod_C1, C2 = mod_C2),
  "output/models/logistic_models.rds"
)

# Save sample sizes per model
n_table <- data.frame(
  model   = c("A1","A2","B1","B2","C1","C2"),
  outcome = c("top5","top1","top5","top1","top5","top1"),
  spec    = c("Demographics","Demographics","Demo+Morbidity","Demo+Morbidity",
              "Demo+Morbidity*Insurance","Demo+Morbidity*Insurance"),
  n       = c(nrow(model.frame(mod_A1)), nrow(model.frame(mod_A2)),
              nrow(model.frame(mod_B1)), nrow(model.frame(mod_B2)),
              nrow(model.frame(mod_C1)), nrow(model.frame(mod_C2)))
)
write.csv(n_table, "output/tables/model_sample_sizes.csv", row.names = FALSE)
cat("\nModel sample sizes saved.\n")
print(n_table)

# -----------------------------------------------------------------------------
# 5. Average Marginal Effects — Set B (main results)
#    avg_slopes() returns AME for each predictor level vs reference
# -----------------------------------------------------------------------------
cat("\nComputing AMEs for Set B...\n")

ame_B1 <- avg_slopes(mod_B1)
ame_B2 <- avg_slopes(mod_B2)

saveRDS(list(B1 = ame_B1, B2 = ame_B2), "output/models/ame_main.rds")

# Export as CSV
ame_B1_df <- as.data.frame(ame_B1)[, c("term","contrast","estimate","std.error","statistic","p.value","conf.low","conf.high")]
ame_B2_df <- as.data.frame(ame_B2)[, c("term","contrast","estimate","std.error","statistic","p.value","conf.low","conf.high")]

write.csv(ame_B1_df, "output/tables/ame_B1_top5.csv", row.names = FALSE)
write.csv(ame_B2_df, "output/tables/ame_B2_top1.csv", row.names = FALSE)

cat("\nAMEs — Model B1 (Top 5%):\n")
print(ame_B1_df)

# -----------------------------------------------------------------------------
# 6. Interaction AMEs — Set C
#    avg_slopes(..., by = "insurance_status") gives AME of morbidity
#    within each insurance stratum — the key extension result
# -----------------------------------------------------------------------------
cat("\nComputing interaction AMEs (morbidity by insurance)...\n")

ame_C1 <- avg_slopes(
  mod_C1,
  variables = "n_chronic_conditions_cat",
  by        = "insurance_status"
)

ame_C2 <- avg_slopes(
  mod_C2,
  variables = "n_chronic_conditions_cat",
  by        = "insurance_status"
)

saveRDS(list(C1 = ame_C1, C2 = ame_C2), "output/models/ame_interaction.rds")

ame_C1_df <- as.data.frame(ame_C1)
ame_C2_df <- as.data.frame(ame_C2)

write.csv(ame_C1_df, "output/tables/ame_C1_interaction.csv", row.names = FALSE)
write.csv(ame_C2_df, "output/tables/ame_C2_interaction.csv", row.names = FALSE)

cat("\nInteraction AMEs — Model C1 (Top 5%):\n")
print(ame_C1_df[, intersect(c("term","contrast","insurance_status","estimate","std.error","p.value"), names(ame_C1_df))])

# -----------------------------------------------------------------------------
# 7. Deviance-based pseudo-R2
#    Note: this is 1 - deviance/null.deviance, a deviance ratio.
#    Useful for incremental model comparison (A vs B vs C), not absolute fit.
#    Comparable only if models A/B/C share the same analytical sample.
# -----------------------------------------------------------------------------
cat("\nComputing deviance-based pseudo-R2...\n")

calc_dev_r2 <- function(model) round(1 - (model$deviance / model$null.deviance), 4)

r2_results <- data.frame(
  Model      = c("A: Demographics only",
                 "B: Demo + Morbidity",
                 "C: Demo + Morbidity x Insurance"),
  Top5_DevR2 = c(calc_dev_r2(mod_A1), calc_dev_r2(mod_B1), calc_dev_r2(mod_C1)),
  Top1_DevR2 = c(calc_dev_r2(mod_A2), calc_dev_r2(mod_B2), calc_dev_r2(mod_C2))
)

write.csv(r2_results, "output/tables/deviance_r2.csv", row.names = FALSE)
cat("\nDeviance-based pseudo-R2 (interpret incrementally, not as absolute fit):\n")
print(r2_results)

cat("\nStep 4 complete.\n")
