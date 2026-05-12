# R/_run_all.R
# =============================================================================
# MEPS High-Cost Healthcare Users — Master Pipeline
# =============================================================================
#
# Runs the full analytical pipeline in order. Run from project root:
#
#     Rscript R/_run_all.R
#
# Each step is wrapped in tryCatch so a failure in one script does not abort
# the rest. A summary of OK / FAIL status is printed at the end.
#
# Pipeline order:
#   01 download  -> 02 clean -> 02b morbidity score
#   03 descriptive
#   04 logistic LEGACY (B1/B2) + 04 main models (P4 interaction)
#   05 glm LEGACY (G-B/G-C) + 05 sensitivity (P5 robustness)
#   06 two-part profiles (post-P6 schema)
#   07 figures (legacy figs 1-6 + sanity check P3/P4/P5)
#   08 report (Rmd narrative; HTML render needs pandoc)
#
# Outputs land in:
#   data/processed/*.rds
#   output/tables/*.csv  output/tables/*.html
#   output/figures/*.png
#   output/models/*.rds  (gitignored)
#   docs/analytical_summary.{Rmd, md, html}
# =============================================================================

scripts <- c(
  "R/01_download.R",
  "R/02_clean.R",
  "R/02b_morbidity_score.R",
  "R/03_descriptive.R",
  "R/04_models_logistic_legacy.R",
  "R/04_main_models.R",
  "R/05_models_glm_legacy.R",
  "R/05_sensitivity.R",
  "R/06_twopart_profiles.R",
  "R/07_figures.R",
  "R/08_report.R"
)

results <- character(length(scripts))
names(results) <- scripts

for (i in seq_along(scripts)) {
  s <- scripts[i]
  cat("\n", strrep("-", 70), "\n", sep = "")
  cat(sprintf("[%d/%d] Running %s\n", i, length(scripts), s))
  cat(strrep("-", 70), "\n", sep = "")

  res <- tryCatch({
    if (!file.exists(s)) {
      stop(sprintf("File not found: %s", s))
    }
    source(s, echo = FALSE)
    "OK"
  }, error = function(e) {
    paste("FAILED:", conditionMessage(e))
  })

  results[s] <- res
  cat(sprintf("\n  -> %s\n", if (res == "OK") "[OK]" else paste("[FAIL]", res)))
}

# ---- Summary ---------------------------------------------------------------
cat("\n", strrep("=", 70), "\n", sep = "")
cat("PIPELINE SUMMARY\n")
cat(strrep("=", 70), "\n", sep = "")

n_ok <- sum(results == "OK")
for (i in seq_along(scripts)) {
  status_tag <- if (results[i] == "OK") "[OK]  " else "[FAIL]"
  msg        <- if (results[i] == "OK") ""
                else paste(" -", substr(results[i], 1, 80))
  cat(sprintf("  %s %s%s\n", status_tag, scripts[i], msg))
}
cat(sprintf("\n%d / %d scripts completed successfully.\n",
            n_ok, length(scripts)))

if (n_ok < length(scripts)) {
  warning(sprintf("%d script(s) failed. See summary above.",
                  length(scripts) - n_ok))
}
