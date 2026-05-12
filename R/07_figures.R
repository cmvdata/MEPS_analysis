# =============================================================================
# 07_figures.R
# MEPS High-Cost Patient Prediction Project
# Step 7: Publication figures
#
# Produces 6 publication-quality figures:
#   Fig 1 — Expenditure distribution (log scale, annotated thresholds)
#   Fig 2 — Forest plot: AMEs from logistic Model B (top 5% and top 1%)
#   Fig 3 — Forest plot: AMEs from GLM Model G-B (dollars, WEIGHTED)
#   Fig 4 — Interaction plot: morbidity x insurance (GLM G-C, WEIGHTED)
#   Fig 5 — Incremental DevR2 across model sets (now includes participation)
#   Fig 6 — Two-part expected spending by profile (6 profiles)
#
# All figures saved as PNG (300 dpi) in output/figures/
# =============================================================================

if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr"))   install.packages("dplyr")
if (!require("scales"))  install.packages("scales")
if (!require("tidyr"))   install.packages("tidyr")

library(ggplot2)
library(dplyr)
library(scales)
library(tidyr)

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Shared theme
# -----------------------------------------------------------------------------
theme_meps <- function() {
  theme_bw(base_size = 12) +
    theme(
      panel.grid.minor  = element_blank(),
      panel.grid.major  = element_line(colour = "grey90"),
      plot.title        = element_text(face = "bold", size = 13),
      plot.subtitle     = element_text(colour = "grey40", size = 10),
      plot.caption      = element_text(colour = "grey50", size = 8, hjust = 0),
      axis.title        = element_text(size = 11),
      legend.position   = "bottom",
      legend.title      = element_text(face = "bold"),
      strip.background  = element_rect(fill = "grey95"),
      strip.text        = element_text(face = "bold")
    )
}

pal_ins <- c(
  "Private"               = "#2166ac",
  "Medicare/Other Public" = "#4dac26",
  "Medicaid/Public Only"  = "#d7191c",
  "Uninsured"             = "#fdae61"
)

# Palette for morbidity levels used in Fig 6 (3-level z-score scale post-P6 refactor)
pal_morb <- c(
  "Low (z = -1)"   = "#9ecae1",
  "Medium (z = 0)" = "#4292c6",
  "High (z = +1)"  = "#08519c"
)

# -----------------------------------------------------------------------------
# Load data and model outputs
# -----------------------------------------------------------------------------
cat("Loading data and model outputs...\n")

df         <- readRDS("data/processed/meps_analytical.rds")

thresholds <- tryCatch(read.csv("output/tables/thresholds.csv"),
                       error = function(e) NULL)
ame_B1_df  <- tryCatch(read.csv("output/tables/ame_B1_top5.csv"),
                       error = function(e) NULL)
ame_B2_df  <- tryCatch(read.csv("output/tables/ame_B2_top1.csv"),
                       error = function(e) NULL)

# Prefer weighted AMEs produced in Step 06; fall back to legacy unweighted.
ame_GB_df <- tryCatch(read.csv("output/tables/ame_GB_weighted.csv"),
                      error = function(e) NULL)
if (is.null(ame_GB_df)) {
  ame_GB_df <- tryCatch(read.csv("output/tables/ame_glm_GB_main.csv"),
                        error = function(e) NULL)
  if (!is.null(ame_GB_df)) cat("  Note: using legacy unweighted G-B AMEs.\n")
}

ame_GC_df <- tryCatch(read.csv("output/tables/ame_GC_weighted.csv"),
                      error = function(e) NULL)
if (is.null(ame_GC_df)) {
  ame_GC_df <- tryCatch(read.csv("output/tables/ame_glm_GC_interaction.csv"),
                        error = function(e) NULL)
  if (!is.null(ame_GC_df)) cat("  Note: using legacy unweighted G-C AMEs.\n")
}

ame_C1_df  <- tryCatch(read.csv("output/tables/ame_C1_interaction.csv"),
                       error = function(e) NULL)
r2_log     <- tryCatch(read.csv("output/tables/deviance_r2.csv"),
                       error = function(e) NULL)
r2_glm     <- tryCatch(read.csv("output/tables/deviance_r2_glm.csv"),
                       error = function(e) NULL)

# Step 06 outputs
twopart_df <- tryCatch(read.csv("output/tables/table_twopart_profiles.csv"),
                       error = function(e) NULL)

thresh_95 <- if (!is.null(thresholds)) thresholds$top5_threshold[1] else quantile(df$total_exp[df$total_exp > 0], 0.95, na.rm = TRUE)
thresh_99 <- if (!is.null(thresholds)) thresholds$top1_threshold[1] else quantile(df$total_exp[df$total_exp > 0], 0.99, na.rm = TRUE)

# -----------------------------------------------------------------------------
# Figure 1: Expenditure distribution
# -----------------------------------------------------------------------------
cat("Fig 1: Expenditure distribution...\n")

df_pos <- df %>% filter(total_exp > 0)

fig1 <- ggplot(df_pos, aes(x = total_exp)) +
  geom_histogram(bins = 80, fill = "#4393c3", colour = "white", alpha = 0.85) +
  geom_vline(xintercept = thresh_95, linetype = "dashed",
             colour = "#d6604d", linewidth = 0.9) +
  geom_vline(xintercept = thresh_99, linetype = "dotted",
             colour = "#762a83", linewidth = 0.9) +
  annotate("text", x = thresh_95 * 1.2, y = Inf, vjust = 1.8,
           label = sprintf("Top 5%%\n$%s",
                           format(round(thresh_95), big.mark = ",")),
           colour = "#d6604d", size = 3.2, hjust = 0) +
  annotate("text", x = thresh_99 * 1.2, y = Inf, vjust = 4.2,
           label = sprintf("Top 1%%\n$%s",
                           format(round(thresh_99), big.mark = ",")),
           colour = "#762a83", size = 3.2, hjust = 0) +
  scale_x_log10(
    labels = label_dollar(scale = 1e-3, suffix = "k"),
    breaks = c(100, 500, 1000, 5000, 10000, 50000, 200000)
  ) +
  labs(
    title    = "Distribution of Total Healthcare Expenditure (2022)",
    subtitle = "Adults 18+, positive expenditure only (log scale). Unweighted counts.",
    x        = "Total expenditure (USD, log scale)",
    y        = "Count",
    caption  = "Source: MEPS 2022 Full-Year Consolidated File (HC-243)."
  ) +
  theme_meps()

ggsave("output/figures/fig1_expenditure_distribution.png",
       fig1, width = 9, height = 5, dpi = 300)
cat("  Saved: fig1_expenditure_distribution.png\n")

# -----------------------------------------------------------------------------
# Figure 2: Logistic forest (B1 and B2)
# -----------------------------------------------------------------------------
cat("Fig 2: Forest plot logistic AMEs...\n")

if (!is.null(ame_B1_df) && !is.null(ame_B2_df)) {

  prep_ame <- function(df_in, outcome_label) {
    df_in %>%
      mutate(
        outcome = outcome_label,
        label   = if ("contrast" %in% names(df_in))
                    paste0(term, ": ", contrast) else term,
        sig     = p.value < 0.05
      ) %>%
      filter(!grepl("Not observed|Intercept|(Intercept)",
                    label, ignore.case = TRUE))
  }

  ame_combined <- bind_rows(
    prep_ame(ame_B1_df, "Top 5%"),
    prep_ame(ame_B2_df, "Top 1%")
  )

  term_order <- prep_ame(ame_B1_df, "Top 5%") %>%
    arrange(estimate) %>% pull(label)
  ame_combined$label <- factor(ame_combined$label, levels = term_order)

  fig2 <- ggplot(ame_combined,
                 aes(x = estimate, y = label,
                     colour = outcome, shape = sig)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                   height = 0.3, alpha = 0.7) +
    geom_point(size = 2.8) +
    scale_colour_manual(
      values = c("Top 5%" = "#2166ac", "Top 1%" = "#d6604d"),
      name = "Outcome") +
    scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                       name = "p < 0.05", labels = c("No", "Yes")) +
    scale_x_continuous(labels = label_percent(accuracy = 0.1)) +
    facet_wrap(~ outcome, scales = "free_x") +
    labs(
      title    = "Average Marginal Effects - Logistic Models (Set B)",
      subtitle = "Probability of being in the high-cost tail (95% CI)",
      x        = "Average Marginal Effect (percentage points)",
      y        = NULL,
      caption  = paste0("Source: MEPS 2022. Survey-weighted quasibinomial GLM.\n",
                        "Reference: Private insurance, High income, Age 18-34, ",
                        "0 chronic conditions, No SPD.")
    ) +
    theme_meps() +
    theme(legend.position = "right")

  ggsave("output/figures/fig2_forest_logistic_ames.png",
         fig2, width = 12, height = 9, dpi = 300)
  cat("  Saved: fig2_forest_logistic_ames.png\n")

} else {
  cat("  Skipped: AME files not found.\n")
}

# -----------------------------------------------------------------------------
# Figure 3: GLM forest (G-B, WEIGHTED AMEs)
# -----------------------------------------------------------------------------
cat("Fig 3: Forest plot GLM AMEs (dollars, weighted)...\n")

if (!is.null(ame_GB_df)) {

  ame_GB_clean <- ame_GB_df %>%
    mutate(
      label = if ("contrast" %in% names(.))
                paste0(term, ": ", contrast) else term,
      sig   = p.value < 0.05
    ) %>%
    filter(!grepl("Not observed|Intercept|(Intercept)",
                  label, ignore.case = TRUE)) %>%
    arrange(estimate)

  ame_GB_clean$label <- factor(ame_GB_clean$label,
                               levels = ame_GB_clean$label)

  fig3 <- ggplot(ame_GB_clean,
                 aes(x = estimate, y = label, colour = sig)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high),
                   height = 0.3, alpha = 0.7) +
    geom_point(size = 2.8) +
    scale_colour_manual(values = c("TRUE" = "#2166ac", "FALSE" = "grey60"),
                        name = "p < 0.05", labels = c("No", "Yes")) +
    scale_x_continuous(labels = label_dollar()) +
    labs(
      title    = "Average Marginal Effects on Total Expenditure - Gamma GLM (G-B)",
      subtitle = paste0("Conditional on positive expenditure. ",
                        "Survey-weighted AMEs in USD (95% CI)."),
      x        = "Average Marginal Effect (USD)",
      y        = NULL,
      caption  = paste0("Source: MEPS 2022. Survey-weighted Gamma GLM, log link.\n",
                        "Reference: Private insurance, High income, Age 18-34, ",
                        "0 chronic conditions, No SPD.")
    ) +
    theme_meps() +
    theme(legend.position = "right")

  ggsave("output/figures/fig3_forest_glm_ames.png",
         fig3, width = 11, height = 9, dpi = 300)
  cat("  Saved: fig3_forest_glm_ames.png\n")

} else {
  cat("  Skipped: GLM AME file not found.\n")
}

# -----------------------------------------------------------------------------
# Figure 4: Interaction (morbidity x insurance)
#   Prefers weighted GLM G-C; falls back to logistic C1 if unavailable.
# -----------------------------------------------------------------------------
cat("Fig 4: Interaction plot...\n")

# Prefer weighted GLM interaction AMEs; fall back to logistic if missing.
use_glm   <- !is.null(ame_GC_df) && "insurance_status" %in% names(ame_GC_df)
int_df    <- if (use_glm) ame_GC_df else ame_C1_df
int_ylab  <- if (use_glm) "AME on Expenditure (USD)" else "AME on Pr(Top 5%) (percentage points)"
int_title <- if (use_glm) {
  "Interaction: Morbidity x Insurance Type - Gamma GLM (G-C, weighted)"
} else {
  "Interaction: Morbidity x Insurance Type - Logistic Model C1 (Top 5%)"
}

if (!is.null(int_df) && "insurance_status" %in% names(int_df)) {

  int_df <- int_df %>%
    filter(!grepl("Not observed", contrast, ignore.case = TRUE)) %>%
    mutate(
      insurance_status = factor(insurance_status, levels = names(pal_ins)),
      contrast         = gsub("mean\\(|\\)", "", contrast)
    )

  y_fmt <- if (use_glm) label_dollar() else label_percent(accuracy = 0.1)

  fig4 <- ggplot(int_df,
                 aes(x = contrast, y = estimate,
                     colour = insurance_status, group = insurance_status)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_line(linewidth = 0.9, alpha = 0.85) +
    geom_point(size = 3.2) +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                  width = 0.2, alpha = 0.65) +
    scale_colour_manual(values = pal_ins, name = "Insurance type") +
    scale_y_continuous(labels = y_fmt) +
    labs(
      title    = int_title,
      subtitle = paste0("AME of each morbidity level vs. 0 chronic conditions, ",
                        "by insurance type (95% CI)"),
      x        = "Number of chronic conditions",
      y        = int_ylab,
      caption  = paste0("Source: MEPS 2022. Survey-weighted model. ",
                        "Reference: 0 chronic conditions, Private insurance.")
    ) +
    theme_meps()

  ggsave("output/figures/fig4_interaction_morbidity_insurance.png",
         fig4, width = 10, height = 6, dpi = 300)
  cat("  Saved: fig4_interaction_morbidity_insurance.png\n")

} else {
  cat("  Skipped: Interaction AME file not found.\n")
}

# -----------------------------------------------------------------------------
# Figure 5: Incremental DevR2
#   Now includes the participation model (Pr(Y > 0) | X) pseudo-R2 at
#   the 'Demo + Morbidity' specification level.
# -----------------------------------------------------------------------------
cat("Fig 5: Incremental DevR2...\n")

if (!is.null(r2_log) && !is.null(r2_glm)) {

  r2_log_clean <- r2_log %>%
    pivot_longer(cols = c(Top5_DevR2, Top1_DevR2),
                 names_to = "outcome", values_to = "DevR2") %>%
    mutate(
      outcome    = recode(outcome,
                          "Top5_DevR2" = "Logistic: Top 5%",
                          "Top1_DevR2" = "Logistic: Top 1%"),
      model_spec = sub("^[^:]+: ", "", Model)
    )

  r2_glm_clean <- r2_glm %>%
    mutate(
      outcome    = "Gamma GLM: Expenditure",
      model_spec = sub("^[^:]+: ", "", Model)
    ) %>%
    select(model_spec, outcome, DevR2)

  # Read participation pseudo-R2 from Step 06 output if available.
  # It is only available at the 'Demo + Morbidity' specification (single model).
  r2_participation <- tryCatch({
    ame_part <- read.csv("output/tables/ame_participation.csv")
    # No explicit R2 CSV; read from the saved participation model
    mod_P <- readRDS("output/models/participation_model.rds")
    dev_r2 <- 1 - mod_P$deviance / mod_P$null.deviance
    data.frame(
      model_spec = "Demo + Morbidity",
      outcome    = "Logistic: Pr(Y > 0)",
      DevR2      = round(dev_r2, 4)
    )
  }, error = function(e) NULL)

  r2_all <- bind_rows(r2_log_clean, r2_glm_clean, r2_participation) %>%
    mutate(model_spec = factor(
      model_spec,
      levels = c("Demographics only",
                 "Demo + Morbidity",
                 "Demo + Morbidity x Insurance")
    ))

  fig5 <- ggplot(r2_all,
                 aes(x = model_spec, y = DevR2,
                     colour = outcome, group = outcome)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 3.5) +
    geom_text(aes(label = sprintf("%.3f", DevR2)),
              vjust = -0.9, size = 3, show.legend = FALSE) +
    scale_colour_manual(
      values = c("Logistic: Top 5%"       = "#2166ac",
                 "Logistic: Top 1%"       = "#d6604d",
                 "Gamma GLM: Expenditure" = "#4dac26",
                 "Logistic: Pr(Y > 0)"    = "#8c510a"),
      name = "Model / Outcome"
    ) +
    scale_y_continuous(limits = c(0, NA),
                       labels = label_number(accuracy = 0.001)) +
    labs(
      title    = "Incremental Model Fit - Deviance-Based Pseudo-R2",
      subtitle = paste0("Pseudo-R2 across model specifications. ",
                        "Participation model estimated only under the ",
                        "Demo + Morbidity specification."),
      x        = "Model specification",
      y        = "Deviance-based pseudo-R2",
      caption  = paste0("Note: Pseudo-R2 = 1 - deviance/null.deviance; ",
                        "interpret incrementally, not as absolute fit. ",
                        "The participation series reports a single point ",
                        "because only one specification was fitted for that outcome.")
    ) +
    theme_meps() +
    theme(axis.text.x = element_text(angle = 15, hjust = 1))

  ggsave("output/figures/fig5_incremental_r2.png",
         fig5, width = 9, height = 5.5, dpi = 300)
  cat("  Saved: fig5_incremental_r2.png\n")

} else {
  cat("  Skipped: R2 files not found.\n")
}

# -----------------------------------------------------------------------------
# Figure 6: Two-part expected spending by profile
#   6 profiles: 2 morbidity levels x 3 insurance types.
#   Reference profile: Female, 55-64, White NH, College+, Middle Income,
#   No SPD, No cancer.
# -----------------------------------------------------------------------------
cat("Fig 6: Two-part expected spending by profile...\n")

if (!is.null(twopart_df) && all(c("total_E_Y", "morbidity_level") %in% names(twopart_df))) {

  # Post-P6 schema (12 profiles: 4 insurance x 3 z-score levels)
  twopart_plot <- twopart_df %>%
    mutate(
      insurance_status = factor(
        insurance_status,
        levels = c("Private", "Medicare/Other Public",
                   "Medicaid/Public Only", "Uninsured")
      ),
      morbidity_level  = factor(
        morbidity_level,
        levels = c("Low (z = -1)", "Medium (z = 0)", "High (z = +1)")
      )
    )

  fig6 <- ggplot(twopart_plot,
                 aes(x = insurance_status, y = total_E_Y,
                     fill = morbidity_level)) +
    geom_col(position = position_dodge(width = 0.8),
             width = 0.7, colour = "white") +
    geom_errorbar(aes(ymin = pmax(0, total_E_Y - total_E_Y_SE),
                      ymax = total_E_Y + total_E_Y_SE),
                  position = position_dodge(width = 0.8),
                  width = 0.2, alpha = 0.6) +
    geom_text(aes(label = sprintf("$%s",
                                  formatC(round(total_E_Y),
                                          format = "d", big.mark = ","))),
              position = position_dodge(width = 0.8),
              vjust = -1.0, size = 2.8, fontface = "bold") +
    scale_fill_manual(values = pal_morb,
                      name = "Morbidity (z-score)") +
    scale_y_continuous(
      labels = label_dollar(scale = 1e-3, suffix = "k"),
      expand = expansion(mult = c(0, 0.18))
    ) +
    labs(
      title    = "Expected Total Healthcare Spending by Insurance and Morbidity",
      subtitle = "Two-part model: E[Y] = Pr(Y > 0) x E[Y | Y > 0]. Error bars: +/- 1 SE.",
      x        = NULL,
      y        = "Expected total spending (USD)",
      caption  = paste0(
        "Reference profile: Female, age 45, South region, Middle Income (200-400% FPL), ",
        "Good self-rated health, no functional limitation, no 2022 inpatient stay.\n",
        "Source: MEPS 2022 (HC-243), survey-weighted two-part model. ",
        "Extensive margin omits any_inpatient_2022 to avoid endogeneity."
      )
    ) +
    theme_meps() +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 15, hjust = 1))

  ggsave("output/figures/fig6_twopart_expected_spending.png",
         fig6, width = 10, height = 6, dpi = 300)
  cat("  Saved: fig6_twopart_expected_spending.png\n")

} else {
  cat("  Skipped: two-part table missing required columns ",
      "(total_E_Y, morbidity_level). Run 06_twopart_profiles.R first.\n")
}

# -----------------------------------------------------------------------------
# Sanity checks for upstream figures from Productos 3-5
# These are produced by 02b_morbidity_score.R / 04_main_models.R / 05_sensitivity.R
# and not regenerated here. Just verify presence.
# -----------------------------------------------------------------------------
cat("\nUpstream figures check:\n")
upstream <- c(
  "P3 morbidity score distribution" = "output/figures/morbidity_score_distribution.png",
  "P4 insurance x morbidity AME"    = "output/figures/mfx_insurance_interaction.png",
  "P5 sensitivity summary (4 panels)"= "output/figures/sensitivity_summary.png"
)
for (lbl in names(upstream)) {
  path <- upstream[[lbl]]
  status <- if (file.exists(path)) "OK" else "MISSING (re-run upstream script)"
  cat(sprintf("  %-40s %s  [%s]\n", lbl, basename(path), status))
}

cat("\nStep 7 complete. All figures saved to output/figures/\n")
