# =============================================================================
# app.R — MEPS High-Cost Users Shiny Portfolio App
# =============================================================================
# Before first deploy, from R:
#   source("deploy.R")
# Or manually:
#   1. Make sure output/figures/*.png exist (run 07_figures.R)
#   2. This app copies figures to www/ at startup if needed
#   3. rsconnect::deployApp(appFiles = c("app.R", "www/", "output/tables/", ...))
# =============================================================================

# ── Packages ──────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(shiny)
  library(shinydashboard)
  library(ggplot2)
  library(dplyr)
  library(DT)
  library(scales)
  library(tidyr)
  library(readr)
})

# ── Startup: ensure www/ has all figures ──────────────────────────────────────
# Shiny only serves static images from www/, so we copy figures there if
# they aren't already. Safe to run on every startup (overwrite = TRUE).
local({
  if (!dir.exists("www")) dir.create("www", showWarnings = FALSE)
  fig_src <- list.files("output/figures", pattern = "\\.png$", full.names = TRUE)
  if (length(fig_src)) {
    file.copy(fig_src, "www/", overwrite = TRUE)
  }
})

# ── Helper: safe CSV loader ───────────────────────────────────────────────────
load_csv <- function(path) {
  if (file.exists(path)) {
    tryCatch(read_csv(path, show_col_types = FALSE),
             error = function(e) NULL)
  } else NULL
}

# ── Load pre-computed outputs ─────────────────────────────────────────────────
df_analytical <- if (file.exists("data/processed/meps_analytical.rds"))
                   readRDS("data/processed/meps_analytical.rds") else NULL

df_twopart    <- load_csv("output/tables/table_twopart_profiles.csv")
df_ame_b1     <- load_csv("output/tables/ame_B1_top5.csv")
df_ame_b2     <- load_csv("output/tables/ame_B2_top1.csv")

# Prefer weighted GLM AMEs; fall back to legacy unweighted
df_ame_glm    <- load_csv("output/tables/ame_GB_weighted.csv")
if (is.null(df_ame_glm))
  df_ame_glm  <- load_csv("output/tables/ame_glm_GB_main.csv")

df_ame_part   <- load_csv("output/tables/ame_participation.csv")
df_r2         <- load_csv("output/tables/deviance_r2_glm.csv")
df_thresholds <- load_csv("output/tables/thresholds.csv")
df_zero       <- load_csv("output/tables/zero_expenditure_stats.csv")

# Post-refactor outputs (Productos 3-5)
df_ame_p4_main    <- load_csv("output/tables/AME_morbidity_by_insurance_main.csv")
df_ame_p4_sens    <- load_csv("output/tables/AME_morbidity_by_insurance_sensitivity.csv")
df_wald_p4        <- load_csv("output/tables/wald_test_modulation.csv")
df_sens_outcomes  <- load_csv("output/tables/sensitivity_outcomes.csv")
df_sens_scores    <- load_csv("output/tables/sensitivity_scores.csv")
df_sens_subgroups <- load_csv("output/tables/sensitivity_subgroups.csv")
df_sens_extremes  <- load_csv("output/tables/sensitivity_extremes.csv")
df_cost_weights   <- load_csv("output/tables/morbidity_cost_weights.csv")

# Figure filenames — these must match the ones emitted by 07_figures.R.
# Images are served from www/, so paths here are relative to www/.
fig_paths <- list(
  hist        = "fig1_expenditure_distribution.png",
  forest_log  = "fig2_forest_logistic_ames.png",
  forest_glm  = "fig3_forest_glm_ames.png",
  interaction = "fig4_interaction_morbidity_insurance.png",
  r2          = "fig5_incremental_r2.png",
  twopart     = "fig6_twopart_expected_spending.png",
  # Post-refactor (Productos 3/4/5)
  mfx_p4      = "mfx_insurance_interaction.png",
  morb_dist   = "morbidity_score_distribution.png",
  sens_summ   = "sensitivity_summary.png"
)

# ── UI ────────────────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "blue",

  dashboardHeader(title = "MEPS 2022 — High-Cost Users"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Overview",          tabName = "overview",     icon = icon("info-circle")),
      menuItem("Expenditure",       tabName = "expenditure",  icon = icon("chart-bar")),
      menuItem("Logistic Models",   tabName = "logistic",     icon = icon("project-diagram")),
      menuItem("GLM Models",        tabName = "glm",          icon = icon("dollar-sign")),
      menuItem("Two-Part Profiles", tabName = "twopart",      icon = icon("users")),
      menuItem("Interaction",       tabName = "interaction",  icon = icon("exchange-alt")),
      menuItem("Insurance × Morbidity", tabName = "p4_interaction", icon = icon("link")),
      menuItem("Sensitivity",       tabName = "p5_sensitivity", icon = icon("flask")),
      menuItem("Morbidity Score",   tabName = "p3_morbidity",   icon = icon("calculator")),
      menuItem("Model Fit",         tabName = "fit",          icon = icon("chart-line")),
      menuItem("Sample Preview",    tabName = "sample",       icon = icon("table"))
    )
  ),

  dashboardBody(
    tabItems(

      # ── Overview ─────────────────────────────────────────────────────────────
      tabItem(tabName = "overview",
        fluidRow(
          box(width = 12,
              title = "Who pays for whom: the economics of concentrated healthcare spending",
              status = "primary", solidHeader = TRUE,
              p("U.S. healthcare spending is extremely concentrated: a small share of individuals accounts for a disproportionate share of annual expenditure. This project uses the 2022 MEPS Full-Year Consolidated File (HC-243) to ask two questions:"),
              tags$ol(
                tags$li("Who ends up in the high-cost tail (top 5% and top 1% of annual expenditure)?"),
                tags$li("How much does the system actually spend on them, conditional on their characteristics?")
              ),
              p("The analysis applies a two-part microeconometric framework: a survey-weighted quasibinomial logistic model for tail membership, a second logistic for participation (Pr(Y > 0)), and a Gamma GLM with log link for conditional expenditure. Average marginal effects are computed with survey weights explicitly passed through."),
              hr(),
              fluidRow(
                valueBoxOutput("box_n",    width = 3),
                valueBoxOutput("box_top5", width = 3),
                valueBoxOutput("box_top1", width = 3),
                valueBoxOutput("box_zero", width = 3)
              ),
              hr(),
              h4("Key findings (post-refactor, Productos 1-6)"),
              tags$ul(
                tags$li(strong("Insurance modulates the morbidity-cost gradient. "),
                        "Wald χ² = 11.04, df = 3, p = 0.012. The slope is ~1.8× steeper in Medicaid ($4,931/SD) than in Private ($2,723/SD). The Uninsured slope is statistically null (p = 0.114) — morbidity does not translate into observed spending without coverage."),
                tags$li(strong("Slope vs level distinction. "),
                        "Medicaid > Private in slope (AME); but Private > Medicaid in level (E[Y]) at moderate morbidity. The Medicaid finding is \"responds more,\" not \"pays more\" — they trace different but complementary patterns."),
                tags$li(strong("Uninsured access barrier quantified. "),
                        "Pr(Y > 0) drops to 50% at low morbidity for Uninsured vs ~87% for Private. The barrier sits at the extensive margin (system contact), not the intensive margin (per-encounter cost)."),
                tags$li(strong("Robust to 11 of 12 sensitivity tests. "),
                        "Outcome alternatives (top10, >$15k, >$50k), score variants, subpopulations, and outlier exclusion preserve the slope ordering Medicaid > Private > Uninsured.")
              )
          )
        )
      ),

      # ── Expenditure distribution ─────────────────────────────────────────────
      tabItem(tabName = "expenditure",
        fluidRow(
          box(width = 12, title = "Figure 1 — Distribution of Total Annual Healthcare Expenditure",
              status = "warning", solidHeader = TRUE,
              p("Log scale. Vertical lines mark the survey-weighted top 5% and top 1% thresholds."),
              uiOutput("fig_hist")
          )
        ),
        fluidRow(
          box(width = 6, title = "Expenditure Thresholds", status = "warning",
              DTOutput("table_thresholds")),
          box(width = 6, title = "Zero Expenditure", status = "warning",
              DTOutput("table_zero"))
        )
      ),

      # ── Logistic ─────────────────────────────────────────────────────────────
      tabItem(tabName = "logistic",
        fluidRow(
          box(width = 12, title = "Figure 2 — Average Marginal Effects: Logistic Models (Set B)",
              status = "success", solidHeader = TRUE,
              p("AMEs on the probability scale (percentage points). Filled points = p < 0.05. Reference: Private insurance, 0 chronic conditions, 18–34, Female, White Non-Hispanic, High Income, No SPD, No cancer."),
              uiOutput("fig_forest_log")
          )
        ),
        fluidRow(
          box(width = 12, title = "AME Table — Top 5% Model (B1)", status = "success",
              DTOutput("table_ame_b1"))
        ),
        fluidRow(
          box(width = 12, title = "AME Table — Top 1% Model (B2)", status = "success",
              DTOutput("table_ame_b2"))
        )
      ),

      # ── GLM ──────────────────────────────────────────────────────────────────
      tabItem(tabName = "glm",
        fluidRow(
          box(width = 12, title = "Figure 3 — Average Marginal Effects: Gamma GLM (Model G-B)",
              status = "danger", solidHeader = TRUE,
              p("Survey-weighted AMEs in dollars on conditional expenditure E[Y | Y > 0]. Restricted to observations with positive total expenditure."),
              uiOutput("fig_forest_glm")
          )
        ),
        fluidRow(
          box(width = 12, title = "AME Table — GLM Main Effects (G-B, weighted)",
              status = "danger",
              DTOutput("table_ame_glm"))
        )
      ),

      # ── Two-Part Profiles ────────────────────────────────────────────────────
      tabItem(tabName = "twopart",
        fluidRow(
          box(width = 12, title = "Two-Part Expected Spending by Profile",
              status = "primary", solidHeader = TRUE,
              p("Reference profile: Female, age 45, South region, Middle Income (200-400% FPL), Good self-rated health, no functional limitation, no 2022 inpatient stay."),
              p(strong("E[Y] = Pr(Y > 0) × E[Y | Y > 0]"),
                ". Morbidity varies in z-score units (low = -1, medium = 0, high = +1) of the AME-weighted score; insurance varies across all 4 categories."),
              p(em("Standard errors in parentheses (delta method, assumes independence between extensive and intensive stages).")),
              uiOutput("fig_twopart"),
              hr(),
              DTOutput("table_twopart")
          )
        )
      ),

      # ── Interaction (legacy: G-C with n_chronic_conditions_cat) ─────────────
      tabItem(tabName = "interaction",
        fluidRow(
          box(width = 12, title = "Figure 4 — Morbidity × Insurance Interaction (legacy)",
              status = "warning", solidHeader = TRUE,
              p("AME of each morbidity level vs. 0 chronic conditions, stratified by insurance type (survey-weighted Gamma GLM, legacy G-C with ", code("n_chronic_conditions_cat"), ")."),
              p(em("Refactor view: see \"Insurance × Morbidity\" tab for AME-weighted score interaction (Producto 4).")),
              uiOutput("fig_interaction")
          )
        )
      ),

      # ── P4: Insurance x Morbidity Interaction (refactor) ────────────────────
      tabItem(tabName = "p4_interaction",
        fluidRow(
          box(width = 12, title = "Insurance × Morbidity Interaction (Producto 4)",
              status = "primary", solidHeader = TRUE,
              p("Survey-weighted Gamma GLM (log link) on intensive margin (", code("total_exp > 0"),
                "). Outcome: conditional expenditure E[Y | Y > 0]. Predictor of interest: ",
                code("morbidity_score_std"), " (AME-weighted, standardized)."),
              p(strong("Wald test on interaction terms: "),
                htmlOutput("wald_inline", inline = TRUE)),
              imageOutput("fig_p4_mfx", height = "auto"),
              hr(),
              h4("AMEs by insurance type — PRIMARY (AME-weighted score, population-weighted)"),
              DTOutput("table_p4_main"),
              hr(),
              h4("AMEs by insurance type — SENSITIVITY (cost-mean log1p score)"),
              DTOutput("table_p4_sens"),
              hr(),
              p(strong("Slope vs level distinction. "),
                "AME (slope) and total expected spending (level) tell complementary stories. Medicaid has the steepest slope ($4,931/SD vs $2,723/SD for Private — 1.8×), but Private maintains higher overall E[Y] at moderate clinical complexity. See \"Two-Part Profiles\" tab for level estimates.")
          )
        )
      ),

      # ── P5: Sensitivity tests ──────────────────────────────────────────────
      tabItem(tabName = "p5_sensitivity",
        fluidRow(
          box(width = 12, title = "Sensitivity tests (Producto 5)",
              status = "info", solidHeader = TRUE,
              p("Four orthogonal robustness dimensions. ",
                strong("Slope ordering Medicaid > Private > Uninsured holds in 11 of 12 fits."),
                " The single failure is a 0.11 percentage-point inversion on top-5% binary outcome (clinically negligible)."),
              imageOutput("fig_p5_sens", height = "auto"),
              hr(),
              tabsetPanel(
                tabPanel("A) Outcome",        DTOutput("table_sens_outcomes")),
                tabPanel("B) Score",          DTOutput("table_sens_scores")),
                tabPanel("C) Subpopulation",  DTOutput("table_sens_subgroups")),
                tabPanel("D) Extremes",       DTOutput("table_sens_extremes"))
              ),
              hr(),
              p(em("Test D (excluding p99+ spenders) reduces insured AMEs by 22-25%; the Uninsured AME paradoxically increases by 42%. Slope ordering preserved; magnitudes should be read with this caveat."))
          )
        )
      ),

      # ── P3: Morbidity Score ────────────────────────────────────────────────
      tabItem(tabName = "p3_morbidity",
        fluidRow(
          box(width = 12, title = "Morbidity score construction (Producto 3)",
              status = "success", solidHeader = TRUE,
              p("Conditions are weighted by their average marginal effect (AME) on conditional expenditure, estimated via survey-weighted Gamma regression on the positive-expenditure subset, controlling for age (cubic spline df=4), sex, poverty, region, self-rated health, functional limitation, and any 2022 inpatient stay."),
              p("Weights are scaled to [0.5, 3.0]. Three conditions with non-positive AMEs after controls (cancer_active, hypertension, prevalent stroke) are clipped to the minimum weight of 0.5."),
              p(strong("Caveat: "),
                "After adding ", code("any_inpatient_2022"), " as a control, the AME of ", code("cancer_active"),
                " collapsed from $21,535 to −$123 (p = 0.90). The original effect was almost entirely the inpatient stay used to construct the active-cancer flag in the first place — i.e., not a finding about cancer per se, but about hospitalization."),
              imageOutput("fig_p3_dist", height = "auto"),
              hr(),
              h4("Cost weights (8 conditions)"),
              DTOutput("table_cost_weights")
          )
        )
      ),

      # ── Model Fit ────────────────────────────────────────────────────────────
      tabItem(tabName = "fit",
        fluidRow(
          box(width = 12, title = "Figure 5 — Deviance-Based Pseudo-R² Across Specifications",
              status = "info", solidHeader = TRUE,
              p("Pseudo-R² = 1 − deviance/null.deviance. Interpret incrementally, not as absolute fit. The participation series reports a single point because only one specification was fitted for that outcome."),
              uiOutput("fig_r2")
          )
        ),
        fluidRow(
          box(width = 12, title = "Model Fit Table (GLM)", status = "info",
              DTOutput("table_r2"))
        )
      ),

      # ── Sample preview ───────────────────────────────────────────────────────
      tabItem(tabName = "sample",
        fluidRow(
          box(width = 12, title = "Analytical Sample — First 500 Rows",
              status = "info", solidHeader = TRUE,
              p("Adults 18+ from MEPS 2022 (HC-243), analytical dataset. All modeling in this app uses survey-weighted estimation."),
              DTOutput("table_sample")
          )
        )
      )

    ) # end tabItems
  )   # end dashboardBody
)

# ── Server ────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  # ── Value boxes ────────────────────────────────────────────────────────────
  output$box_n <- renderValueBox({
    n <- if (!is.null(df_analytical)) nrow(df_analytical) else NA
    valueBox(
      if (is.na(n)) "—" else format(n, big.mark = ","),
      "Adults 18+ (n)", icon = icon("users"), color = "blue"
    )
  })
  output$box_top5 <- renderValueBox({
    val <- if (!is.null(df_thresholds) && "top5_threshold" %in% names(df_thresholds))
             dollar(df_thresholds$top5_threshold[1]) else "—"
    valueBox(val, "Top 5% threshold", icon = icon("arrow-up"), color = "yellow")
  })
  output$box_top1 <- renderValueBox({
    val <- if (!is.null(df_thresholds) && "top1_threshold" %in% names(df_thresholds))
             dollar(df_thresholds$top1_threshold[1]) else "—"
    valueBox(val, "Top 1% threshold", icon = icon("exclamation"), color = "red")
  })
  output$box_zero <- renderValueBox({
    val <- if (!is.null(df_zero)) {
      # Try several possible column names for robustness
      col <- intersect(c("p_zero", "pct_zero", "prop_zero"), names(df_zero))
      if (length(col)) sprintf("%.1f%%", df_zero[[col[1]]][1] * 100) else "14.5%"
    } else "—"
    valueBox(val, "Zero expenditure", icon = icon("times-circle"), color = "green")
  })

  # ── Figure renderer (from www/) ────────────────────────────────────────────
  render_fig <- function(filename) {
    renderUI({
      path_in_www <- file.path("www", filename)
      if (file.exists(path_in_www)) {
        tags$img(src = filename,
                 style = "width:100%; max-width:1100px; display:block; margin:auto;")
      } else {
        div(style = "color:gray; padding:20px; text-align:center;",
            em(paste("Figure not found:", filename)),
            br(),
            em("Run 07_figures.R to regenerate."))
      }
    })
  }

  output$fig_hist        <- render_fig(fig_paths$hist)
  output$fig_forest_log  <- render_fig(fig_paths$forest_log)
  output$fig_forest_glm  <- render_fig(fig_paths$forest_glm)
  output$fig_interaction <- render_fig(fig_paths$interaction)
  output$fig_r2          <- render_fig(fig_paths$r2)
  output$fig_twopart     <- render_fig(fig_paths$twopart)

  # ── Tables ─────────────────────────────────────────────────────────────────
  empty_msg <- function(msg) {
    datatable(data.frame(Message = msg),
              options = list(dom = "t"), rownames = FALSE)
  }

  fmt_ame <- function(df) {
    if (is.null(df)) return(empty_msg("Run pipeline to generate AMEs."))
    df %>%
      mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
      datatable(options = list(pageLength = 15, scrollX = TRUE),
                rownames = FALSE)
  }

  output$table_ame_b1  <- renderDT(fmt_ame(df_ame_b1))
  output$table_ame_b2  <- renderDT(fmt_ame(df_ame_b2))
  output$table_ame_glm <- renderDT(fmt_ame(df_ame_glm))

  output$table_thresholds <- renderDT({
    if (is.null(df_thresholds)) return(empty_msg("Run 02_clean.R"))
    datatable(df_thresholds, options = list(dom = "t"), rownames = FALSE)
  })

  output$table_zero <- renderDT({
    if (is.null(df_zero)) return(empty_msg("Run 03_descriptive.R"))
    datatable(df_zero, options = list(dom = "t"), rownames = FALSE)
  })

  output$table_r2 <- renderDT({
    if (is.null(df_r2)) return(empty_msg("Run 05_models_glm.R"))
    df_r2 %>%
      mutate(across(where(is.numeric), ~ round(.x, 4))) %>%
      datatable(options = list(dom = "t"), rownames = FALSE)
  })

  # ── Two-part table (post-P6 schema with legacy fallback) ───────────────────
  output$table_twopart <- renderDT({
    if (is.null(df_twopart)) return(empty_msg("Run 06_twopart_profiles.R"))
    if ("total_E_Y" %in% names(df_twopart)) {
      # Post-P6 schema: 12 profiles (4 insurance x 3 z-score levels)
      df_twopart %>%
        mutate(
          `Insurance`    = insurance_status,
          `Morbidity`    = morbidity_level,
          `Pr(Y > 0)`    = sprintf("%.1f%% (%.1f)",
                                   Pr_Y_pos * 100, Pr_Y_pos_SE * 100),
          `E[Y | Y > 0]` = sprintf("$%s ($%s)",
                                   formatC(round(E_Y_given_pos),
                                           format = "d", big.mark = ","),
                                   formatC(round(E_Y_given_pos_SE),
                                           format = "d", big.mark = ",")),
          `E[Y]`         = sprintf("$%s ($%s)",
                                   formatC(round(total_E_Y),
                                           format = "d", big.mark = ","),
                                   formatC(round(total_E_Y_SE),
                                           format = "d", big.mark = ","))
        ) %>%
        select(Insurance, Morbidity, `Pr(Y > 0)`, `E[Y | Y > 0]`, `E[Y]`) %>%
        datatable(options = list(dom = "t", pageLength = 12),
                  rownames = FALSE,
                  caption = "Standard errors in parentheses (delta method).")
    } else {
      # Legacy schema fallback
      df_twopart %>%
        mutate(
          `Insurance`    = insurance_status,
          `Morbidity`    = n_chronic_conditions_cat,
          `Pr(Y > 0)`    = sprintf("%.1f%%", pr_pos  * 100),
          `Pr(Top 5%)`   = sprintf("%.1f%%", pr_top5 * 100),
          `Pr(Top 1%)`   = sprintf("%.1f%%", pr_top1 * 100),
          `E[Y | Y > 0]` = dollar(e_y_pos),
          `E[Y]`         = dollar(e_y)
        ) %>%
        select(Insurance, Morbidity, `Pr(Y > 0)`,
               `Pr(Top 5%)`, `Pr(Top 1%)`, `E[Y | Y > 0]`, `E[Y]`) %>%
        datatable(options = list(dom = "t", pageLength = 10),
                  rownames = FALSE)
    }
  })

  # ── Sample preview ─────────────────────────────────────────────────────────
  output$table_sample <- renderDT({
    if (is.null(df_analytical)) return(empty_msg("Run 02_clean.R"))
    wanted <- c("insurance_status", "n_chronic_conditions_cat", "age_cat",
                "sex", "race_ethnicity", "education", "income_pct_fpl",
                "cancer_ever", "mental_health_cat", "top5", "top1", "total_exp")
    vars <- intersect(wanted, names(df_analytical))
    df_analytical[1:min(500, nrow(df_analytical)), vars, drop = FALSE] |>
      datatable(options = list(pageLength = 15, scrollX = TRUE),
                caption = "First 500 rows of the analytical dataset",
                rownames = FALSE)
  })

  # ── P4: Interaction tables, Wald test, image ────────────────────────────────
  fmt_p4 <- function(df) {
    if (is.null(df)) return(empty_msg("Run R/04_main_models.R"))
    df %>%
      mutate(
        `Insurance` = insurance_status,
        `AME (USD)` = dollar(estimate),
        `SE`        = dollar(std.error),
        `95% CI`    = sprintf("[%s; %s]",
                              dollar(conf.low), dollar(conf.high)),
        `p`         = formatC(p.value, format = "e", digits = 2)
      ) %>%
      select(Insurance, `AME (USD)`, SE, `95% CI`, p) %>%
      datatable(options = list(dom = "t"), rownames = FALSE)
  }
  output$table_p4_main <- renderDT(fmt_p4(df_ame_p4_main))
  output$table_p4_sens <- renderDT(fmt_p4(df_ame_p4_sens))

  output$wald_inline <- renderUI({
    if (is.null(df_wald_p4)) return(em("Wald test not available."))
    HTML(sprintf("χ² = %.2f, df = %.0f, <b>p = %.4f</b>",
                 df_wald_p4$chisq_statistic[1],
                 df_wald_p4$df[1],
                 df_wald_p4$p_value[1]))
  })

  output$fig_p4_mfx <- renderImage({
    src <- file.path("www", fig_paths$mfx_p4)
    if (!file.exists(src)) {
      return(list(src = "", contentType = "image/png", width = 0,
                  alt = "mfx_insurance_interaction.png missing"))
    }
    list(src = src, contentType = "image/png",
         width = 900, alt = "AME of morbidity by insurance type")
  }, deleteFile = FALSE)

  # ── P5: Sensitivity tables and image ────────────────────────────────────────
  fmt_sens <- function(df) {
    if (is.null(df)) return(empty_msg("Run R/05_sensitivity.R"))
    is_dollar <- any(abs(df$estimate) > 100, na.rm = TRUE)
    df %>%
      mutate(
        `Specification` = spec,
        `Insurance`     = insurance_status,
        `AME`           = if (is_dollar) dollar(estimate)
                          else sprintf("%.4f", estimate),
        `95% CI`        = if (is_dollar)
                            sprintf("[%s; %s]",
                                    dollar(conf.low), dollar(conf.high))
                          else
                            sprintf("[%.4f; %.4f]",
                                    conf.low, conf.high),
        `p`             = formatC(p.value, format = "e", digits = 2)
      ) %>%
      select(Specification, Insurance, AME, `95% CI`, p) %>%
      datatable(options = list(pageLength = 16), rownames = FALSE)
  }
  output$table_sens_outcomes  <- renderDT(fmt_sens(df_sens_outcomes))
  output$table_sens_scores    <- renderDT(fmt_sens(df_sens_scores))
  output$table_sens_subgroups <- renderDT(fmt_sens(df_sens_subgroups))
  output$table_sens_extremes  <- renderDT(fmt_sens(df_sens_extremes))

  output$fig_p5_sens <- renderImage({
    src <- file.path("www", fig_paths$sens_summ)
    if (!file.exists(src)) {
      return(list(src = "", contentType = "image/png", width = 0,
                  alt = "sensitivity_summary.png missing"))
    }
    list(src = src, contentType = "image/png",
         width = 1100, alt = "Sensitivity summary across 4 dimensions")
  }, deleteFile = FALSE)

  # ── P3: Morbidity score weights table and image ────────────────────────────
  output$table_cost_weights <- renderDT({
    if (is.null(df_cost_weights)) return(empty_msg("Run R/02b_morbidity_score.R"))
    df_cost_weights %>%
      mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
      datatable(options = list(dom = "t", pageLength = 8),
                rownames = FALSE,
                caption = "AMEs from Gamma GLM with full controls; weights scaled to [0.5, 3.0].")
  })

  output$fig_p3_dist <- renderImage({
    src <- file.path("www", fig_paths$morb_dist)
    if (!file.exists(src)) {
      return(list(src = "", contentType = "image/png", width = 0,
                  alt = "morbidity_score_distribution.png missing"))
    }
    list(src = src, contentType = "image/png",
         width = 800, alt = "Distribution of morbidity score (z-scored)")
  }, deleteFile = FALSE)
}

# ── Run ───────────────────────────────────────────────────────────────────────
shinyApp(ui, server)
