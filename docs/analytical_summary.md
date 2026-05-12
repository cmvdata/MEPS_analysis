---
title: 'Cross-sectional characterization of high-cost healthcare users in the US'
subtitle: 'How the morbidity-cost gradient varies by insurance type in MEPS 2022'
author: 'Health Services Analyst Portfolio Project'
date: '`r Sys.Date()`'
output:
  html_document:
    theme: cosmo
    toc: true
    toc_float: true
---

## 1. Research question and design

U.S. healthcare spending is famously concentrated, but the structural reasons behind that concentration depend on who pays. This project uses the 2022 MEPS Full-Year Consolidated File (HC-243) to ask **how the morbidity-cost gradient varies by insurance type** - that is, the rate at which clinical complexity translates into observed healthcare spending.

The analysis is **cross-sectional by design**, not longitudinal. The MEPS panel rotates: linking 2022 respondents to 2021 (HC-233) loses ~40-50% of the sample and shifts the question from *"who are the high-cost users in 2022?"* to *"who became high-cost users between 2021 and 2022?"*. We chose the cross-sectional framing because it preserves statistical power (n = 17,909 adults), aligns with the original Maynou et al. (2023) design, and matches the policy question of interest: characterizing the standing relationship between morbidity, insurance, and spending in a given fiscal year.

## 2. Data and methods

### Sample and outcome
- **Dataset**: MEPS HC-243 2022 Full-Year Consolidated File
- **Sample**: 17,909 adults aged 18+
- **Primary outcome**: total annual healthcare expenditure (`TOTEXP22`), modeled both as a binary tail-membership indicator (Top 5%, Top 1%) and as a continuous variable conditional on positive expenditure (`Y > 0`)

### Covariates introduced beyond the legacy specification
- **Self-rated health** (`health_status_3`): RTHLTH31 collapsed to Excellent/Very good / Good / Fair-Poor
- **Functional limitation** (`any_functional_limitation`): need for help with ADLs or IADLs (round-1 measurement; HC-243 lacks a unified round-4 binary)
- **Concurrent inpatient utilization** (`any_inpatient_2022`): used as a covariate in the intensive margin only - see endogeneity note below
- **Cancer split**: `cancer_active` proxied as `cancer_ever AND any_inpatient_2022` (HC-242B refinement deferred); `cancer_remission` is the complement among `cancer_ever`

### Cost-weighted morbidity score (revised)
Conditions are weighted by their **average marginal effect** on conditional expenditure, estimated via survey-weighted Gamma regression (log link) on the positive-expenditure subset, controlling for age (cubic spline df=4), sex, poverty, region, self-rated health, functional limitation, and any inpatient stay in 2022. Weights are scaled to [0.5, 3.0] for numerical stability. The top three weights are assigned to **diabetes (3.00), arthritis (3.00), and asthma (2.60)** - chronic ambulatory conditions whose costs flow through pharmacy and outpatient management rather than acute hospitalization. Three conditions with non-positive AMEs after controls (`cancer_active`, hypertension, prevalent stroke) are clipped to the minimum weight of 0.5.

A key methodological finding emerges from the inpatient control: the AME of `cancer_active` collapses from **\$21,535** (uncontrolled) to **-\$123** (p = 0.90) once `any_inpatient_2022` enters the regression. The original effect was almost entirely attributable to the inpatient stay used to construct the active-cancer flag in the first place. Cross-sectional MEPS cannot identify "active cancer" net of "hospitalization" without HC-242B encounter-level diagnosis codes - a documented limitation rather than a finding.

### Two-part model: extensive and intensive margins
- **Extensive margin** (`Pr(Y > 0)`): survey-weighted quasibinomial logit. **`any_inpatient_2022` is omitted** from this stage to avoid circular endogeneity - an inpatient discharge implies positive expenditure by construction, so including it as a predictor of participation would estimate a tautology.
- **Intensive margin** (`E[Y | Y > 0]`): survey-weighted Gamma GLM with log link, restricted to `total_exp > 0`. `any_inpatient_2022` is included as a covariate to separate the cancer-active effect from the hospitalization effect.

This separation lets us decompose `E[Y | X] = Pr(Y > 0 | X) x E[Y | Y > 0, X]` and read each margin in its own clinical/policy register: extensive = access; intensive = intensity.

### Inference
All effects are reported as Average Marginal Effects (AMEs) on the response scale, computed via `marginaleffects::avg_slopes` with population weights (`PERWT22F`) explicitly passed via `wts = mod$prior.weights`. Standard errors use the delta method. The Wald test on interaction terms uses `car::linearHypothesis` on the svyglm object (returns Chi-squared, not F, for survey designs).

## 3. Findings

### 3a. Insurance modulates the morbidity-cost gradient (Wald p = 0.012)

The morbidity x insurance interaction in the intensive-margin Gamma GLM is statistically significant: **Wald Chi-squared = 11.04, df = 3, p = 0.012**. Per-insurance AMEs of a 1-SD increase in `morbidity_score_std` on conditional expenditure (population-weighted, USD):

| Insurance | AME (USD) | 95% CI | p |
|---|---|---|---|
| Private | 2,723 | [2,065; 3,382] | <0.001 |
| Medicare/Other Public | 2,123 | [1,383; 2,863] | <0.001 |
| **Medicaid/Public Only** | **4,931** | [3,310; 6,552] | <0.001 |
| Uninsured | 677 | [-163; 1,517] | 0.114 (n.s.) |

The slope of the morbidity-cost relationship is **roughly 1.8x steeper in Medicaid than in Private** (\$4,931 vs \$2,723 per SD). The Uninsured slope is statistically indistinguishable from zero in the primary specification (p = 0.114) and only marginal under the alternative cost-mean log1p score (p = 0.060).

**Plausible mechanisms for the Medicaid pattern**:

1. *Selection on residual severity.* Medicaid eligibility (low income, disability) selects individuals who, conditional on the same morbidity score, may carry residual severity not captured by the score. Each marginal SD of clinical complexity therefore translates into more spending.
2. *Reduced patient cost-sharing.* Medicaid covers copays and deductibles that Private plans typically leave to enrollees, removing demand-side friction that would otherwise dampen spending growth at high morbidity.
3. *Care management for high-need beneficiaries.* Medicaid-funded case management programs can mobilize utilization (specialist referrals, home-care, durable medical equipment) for clinically complex enrollees, expanding the spending response to morbidity.

These mechanisms are consistent with the data but cannot be separated in cross-section - only that the steeper slope is reproducible across specifications (see Robustness).

### 3b. The uninsured access barrier, quantified

The two-part decomposition exposes a feature that single-equation models obscure. For a representative profile (age 45, Female, middle income 200-400% FPL, South region, Good self-rated health, no functional limitation, no 2022 inpatient stay), predicted **`Pr(Y > 0)`** by insurance and morbidity:

| Morbidity (z) | Private | Medicare | Medicaid | **Uninsured** |
|---|---|---|---|---|
| Low (-1) | 0.869 | 0.797 | 0.774 | **0.499** |
| Medium (0) | 0.928 | 0.887 | 0.887 | **0.665** |
| High (+1) | 0.961 | 0.940 | 0.948 | **0.799** |

A **uninsured adult with low-to-medium morbidity has a 33-50% probability of incurring no observable healthcare spending in 2022**, against ~7-23% for the insured at the same morbidity. This is the empirical signature of an access barrier - not lower clinical need, but suppressed contact with the system at the extensive margin.

The slope of `Pr(Y > 0)` across morbidity levels is also markedly steeper for the uninsured (Delta = +0.30 from low to high) than for the insured (Delta ~ +0.09 Private, +0.14 Medicare, +0.17 Medicaid), suggesting that morbidity is a stronger driver of system contact when insurance does not provide a stable baseline of access. Among the insured, contact is near-universal; among the uninsured, contact is driven by clinical urgency.

### 3c. Slope vs level: two complementary readings

A common interpretive trap is to report "Medicaid > Private" without specifying which margin. The AME (slope) and total expected spending (level) tell different stories.

The intensive-margin AMEs from Section 3a (population-averaged across the sample) are:

- Private: **\$2,723 per 1-SD** increase in morbidity
- Medicaid: **\$4,931 per 1-SD** increase (~1.8x steeper)

At the reference profile of Section 3b, conditional spending E[Y | Y > 0] differs across morbidity extremes (z = -1 vs z = +1, a 2-SD range) by:

- Private: \$9,661 - \$5,182 = **\$4,479** (profile-specific local slope ~ \$2,240/SD)
- Medicaid: \$9,041 - \$4,252 = **\$4,789** (profile-specific local slope ~ \$2,395/SD)

The profile-specific slopes are *flatter* than the population AMEs because the reference profile (age 45, healthy, no functional limitation, no 2022 inpatient) sits well below the population mean E[Y]. For a Gamma log-link model the local marginal effect is `E[Y] x beta`, so the marginal dollar response shrinks at low-spending profiles. Both numbers are correct: the population AME describes the average patient; the profile-specific slope describes a young, healthy, modal-demographic patient.

For total expected spending **E[Y] = Pr(Y > 0) x E[Y | Y > 0]**:

| Morbidity | E[Y] Private | E[Y] Medicaid | Ranking | Gap |
|---|---|---|---|---|
| Low (z=-1) | \$4,506 | \$3,293 | Private > Medicaid | \$1,213 |
| Medium (z=0) | \$6,564 | \$5,502 | Private > Medicaid | \$1,062 |
| High (z=+1) | \$9,285 | \$8,568 | Private > Medicaid | \$717 |

At this profile, total expected spending is higher in Private than in Medicaid at every morbidity level, but the gap narrows steadily (from \$1,213 at low morbidity to \$717 at high). Both findings - Medicaid's steeper slope in AMEs and Private's higher level in profile predictions - are simultaneously true and complementary:

- **Private has a higher intercept**: better baseline access (`Pr(Y > 0)` 0.87 vs 0.77 at z = -1) and higher conditional spending at low morbidity.
- **Medicaid has a steeper slope**: each marginal SD of morbidity adds more dollars to expected spending in Medicaid than in Private. Total E[Y] in Medicaid grows by \$5,275 across the 2-SD morbidity range (\$8,568 - \$3,293), versus \$4,779 in Private (\$9,285 - \$4,506).

Reported jointly: *"Among the insured, the morbidity-cost gradient is steepest in Medicaid, but Private maintains higher overall expected spending at moderate clinical complexity."* The crossover between Private and Medicaid is not realized within the [-1, +1] z-score range at this profile, but the trajectory points toward convergence at very high morbidity.

The Medicaid-Private differential is *not* a "Medicaid pays more" story - it is a "Medicaid responds more" story, with policy implications that differ from a pure level interpretation.

### 3d. Robustness

Twelve sensitivity fits across four orthogonal dimensions. **The slope ordering Medicaid > Private > Uninsured holds in 11 / 12 fits.** The single failure is a 0.11 percentage-point inversion between Medicaid and Private on the top-5% binary outcome - clinically negligible.

| Dimension | Specifications tested | Pattern preserved |
|---|---|---|
| A) Outcome (binary) | top5, top10, >\$15k, >\$50k | 3 / 4 |
| B) Score | AME-weighted, cost-mean log1p, no-cancer | 3 / 3 |
| C) Subpopulation | <65, no functional limit, SAQ subsample | 3 / 3 |
| D) Outliers | full vs excluding p99+ | 2 / 2 |

Test D (outlier influence) is the most informative diagnostic. The AMEs of the insured groups are inflated by approximately **22-25% by the top 1% of spenders** (the threshold is \$88,317): excluding observations above p99 reduces Private's AME from \$2,723 to \$2,117, Medicare's from \$2,123 to \$1,735, and Medicaid's from \$4,931 to \$3,689. The Uninsured AME, paradoxically, *increases* from \$677 to \$962 - only a small base of extreme uninsured spenders, so their exclusion gives more weight to the mid-range gradient. The slope ordering is preserved, but **point estimates of magnitude should be read with the outlier caveat in mind**.

## 4. Implications

**For health policy.** Coverage by itself is not sufficient to flatten cost concentration. The modulation of the morbidity-cost gradient - not its existence - is what differs across insurance environments. Risk adjustment models calibrated against Private populations may underpay for the same morbidity score in Medicaid (whose response slope is roughly 1.8x steeper) and over-attribute "missing" spending to the uninsured (whose suppressed extensive margin reflects barriers, not low need).

**For care management strategy.** The steepest slope being in Medicaid is informative for plan design: the marginal high-cost member in Medicaid is not simply a Private-equivalent member with worse coverage - they are clinically and economically different in ways the morbidity score only partially captures. Capitation rate-setting that ignores the slope differential will systematically underfund high-morbidity Medicaid enrollees.

**For research.** Cross-sectional MEPS allows characterization but not causal identification. A longitudinal panel (HC-243 <-> HC-233) with appropriate IV or DiD design would be needed to identify whether the modulation reflects causal effects of coverage on spending response, or compositional differences across enrollee populations that coexist with coverage. Either way, the descriptive fact - that the morbidity-cost slope is steepest in Medicaid and flattest among the uninsured - is robust within 2022 and warrants attention.

## 5. Limitations

- **Cross-sectional design.** Findings characterize 2022 contemporaneous associations, not causal effects. Reverse causation (high morbidity -> enrollment in Medicaid via disability eligibility) is not identified.
- **Cancer active vs. remission proxy.** We use `cancer_ever AND any_inpatient_2022` as a proxy for active cancer. Approximately 60-70% of cancer treatment occurs in outpatient settings (NCI 2022), suggesting the misclassification may affect 30-50% of true active-cancer cases. The downward bias in `cancer_active` prevalence likely attenuates its weight in the morbidity score; HC-242B condition-level data would refine this, at the cost of substantial pipeline complexity.
- **Mental health complete-case attrition.** The K6-derived `mental_health` variable is observed only for the SAQ subsample (~60% of adults); the AME-weighted score uses NA -> 0 for individuals without K6, biasing the SPD weight downward. SPD-related coefficients should be interpreted as conservative.
- **MEPS underreporting.** OTC purchases, unreimbursed copays, and informal care are systematically underreported in MEPS. Total expenditure reflects administrative + patient-paid components, not the full economic burden of illness.
- **Outlier influence.** The intensive margin AMEs are inflated by approximately 22-25% by the top 1% of spenders (Test D); robust to exclusion in direction, sensitive in magnitude. Reported point estimates should be read as upper bounds for typical patterns.
- **No URBAN/MSA covariate.** `URBAN22` is not preserved in the analytical pipeline; rurality is partially absorbed via `REGION22` but a dedicated rural/urban control would strengthen identification.

## 6. Technical summary

- **Sample**: 17,909 adults (18+) from MEPS 2022 HC-243; 15,315 with positive expenditure; effective N for intensive Gamma fit = 15,155 (99% retention)
- **Survey design**: stratified multi-stage cluster (`VARSTR`, `VARPSU`, `PERWT22F`); `survey::svydesign` with `nest = TRUE`, `survey.lonely.psu = "adjust"`
- **Models**:
  - *Extensive*: `any_exp ~ insurance x morbidity_score_std + ns(AGE22X, df=4) + sex + POVCAT22 + REGION22 + health_status_3 + any_functional_limitation`, family = quasibinomial(logit). `any_inpatient_2022` excluded.
  - *Intensive*: as above plus `+ any_inpatient_2022`, family = Gamma(log), restricted to `total_exp > 0`.
- **Inference**: `marginaleffects::avg_slopes` with population weights (`wts = mod$prior.weights`), 95% CI via delta method. Wald test on interaction terms via `car::linearHypothesis` (returns Chi-squared, not F, for survey designs).
- **Score construction**: 8 conditions, weights from AMEs in the Gamma fit, scaled to [0.5, 3.0]; non-positive AMEs clipped to 0.5; standardized via `scale()` for analysis.
- **Two-part predictions**: manual via `model.matrix(delete.response(terms(model)), combined_data)` x `coef(model)` to bypass `predict.svyglm` issues with `splines::ns()` predvars in newdata. SE via `vcov` + delta method, assuming independence between extensive and intensive stages.
- **Stack**: R 4.5, `survey`, `marginaleffects`, `splines`, `gtsummary`, `gt`, `ggplot2`, `dplyr`.


### Supplementary table - full two-part profile predictions

*Reference profile: age 45, Female, Middle Income (200-400% FPL), South, Good self-rated health, no functional limitation, no 2022 inpatient stay. Standard errors in parentheses (delta method).*

| Insurance | Morbidity | Pr(Y>0) | E[Y\|Y>0] | E[Y] |
|---|---|---|---|---|
| Private | Low (z = -1) | 0.869 (0.013) | $5,182 ($375) | $4,506 ($332) |
| Private | Medium (z = 0) | 0.928 (0.007) | $7,076 ($467) | $6,564 ($436) |
| Private | High (z = +1) | 0.961 (0.006) | $9,661 ($733) | $9,285 ($707) |
| Medicare/Other Public | Low (z = -1) | 0.797 (0.043) | $5,017 ($501) | $3,998 ($455) |
| Medicare/Other Public | Medium (z = 0) | 0.887 (0.022) | $6,102 ($562) | $5,411 ($516) |
| Medicare/Other Public | High (z = +1) | 0.940 (0.016) | $7,422 ($701) | $6,976 ($669) |
| Medicaid/Public Only | Low (z = -1) | 0.774 (0.026) | $4,252 ($516) | $3,293 ($415) |
| Medicaid/Public Only | Medium (z = 0) | 0.887 (0.014) | $6,200 ($672) | $5,502 ($602) |
| Medicaid/Public Only | High (z = +1) | 0.948 (0.010) | $9,041 ($1,059) | $8,568 ($1,008) |
| Uninsured | Low (z = -1) | 0.499 (0.034) | $2,758 ($600) | $1,375 ($313) |
| Uninsured | Medium (z = 0) | 0.665 (0.032) | $3,312 ($463) | $2,203 ($326) |
| Uninsured | High (z = +1) | 0.799 (0.045) | $3,977 ($576) | $3,177 ($495) |



