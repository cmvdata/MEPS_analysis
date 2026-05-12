# =============================================================================
# deploy.R — Deploy the MEPS Shiny app to shinyapps.io
#
# Before running this the FIRST time:
#   1. Install rsconnect:    install.packages("rsconnect")
#   2. Set account info once (from your shinyapps.io Tokens page):
#      rsconnect::setAccountInfo(name='cmvdata', token='...', secret='...')
#
# Then just:   source("deploy.R")
# =============================================================================

# ── Sanity: make sure we're in the project root ──────────────────────────────
if (!file.exists("app.R")) {
  stop("app.R not found. Set working directory to the project root first:\n",
       '   setwd("path/to/MEPS_analysis")')
}

# ── Prep www/: copy all figures so Shiny can serve them ──────────────────────
if (!dir.exists("www")) dir.create("www", showWarnings = FALSE)
fig_src <- list.files("output/figures", pattern = "\\.png$", full.names = TRUE)
if (!length(fig_src)) {
  stop("No figures found in output/figures/. Run R/07_figures.R first.")
}
n_copied <- sum(file.copy(fig_src, "www/", overwrite = TRUE))
cat(sprintf("Copied %d figure(s) to www/\n", n_copied))

# ── Build the explicit file list for deployment ──────────────────────────────
required_csvs <- c(
  "output/tables/table_twopart_profiles.csv",
  "output/tables/ame_B1_top5.csv",
  "output/tables/ame_B2_top1.csv",
  "output/tables/ame_GB_weighted.csv",      # preferred (weighted)
  "output/tables/ame_glm_GB_main.csv",      # fallback (unweighted)
  "output/tables/ame_participation.csv",
  "output/tables/deviance_r2_glm.csv",
  "output/tables/thresholds.csv",
  "output/tables/zero_expenditure_stats.csv"
)

required_data <- c("data/processed/meps_analytical.rds")

# Only include files that actually exist (missing ones simply show "—" in app)
app_files <- c(
  "app.R",
  list.files("www", pattern = "\\.png$", full.names = TRUE),
  Filter(file.exists, required_csvs),
  Filter(file.exists, required_data)
)

cat("\nFiles to deploy (", length(app_files), "):\n", sep = "")
cat(paste(" -", app_files), sep = "\n")

# ── Deploy ───────────────────────────────────────────────────────────────────
cat("\nDeploying to shinyapps.io...\n")
rsconnect::deployApp(
  appDir      = getwd(),
  appName     = "meps-high-cost-users",
  appFiles    = app_files,
  forceUpdate = TRUE,
  launch.browser = TRUE
)
