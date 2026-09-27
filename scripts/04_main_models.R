# =============================================================================
# Dental status and MASLD / liver fibrosis — NHANES 2017–March 2020
# Script 04: main survey-weighted models (Aims 1 and 2) + dose-response splines
# -----------------------------------------------------------------------------
# Input : data/processed/nhanes_merged.rds   (created by Script 02)
# Output: output/tables/main_models.csv         (every OR from every model)
#         output/tables/main_models_primary.csv (Model 2 only — the headline table)
#         output/figures/spline_masld.png
#         output/figures/spline_fibrosis.png
#
# Outcomes : MASLD; significant fibrosis (LSM >= 8 kPa); cirrhosis range (>= 13.6)
# Exposures: tooth-loss category, missing teeth (per 5), fewer than 20 teeth;
#            among DENTATE adults only: any untreated caries, decayed teeth, DMFT
# Models   : 1 = age (spline) + sex + race
#            2 = Model 1 + education + income + smoking + alcohol + insurance  (PRIMARY)
#            3 = Model 2 + BMI + diabetes                     (possible mediators)
#            4 = Model 2 + log hs-CRP                         (possible mediator)
# =============================================================================

library(tidyverse)
library(survey)
library(splines)
library(here)

options(survey.lonely.psu = "adjust")

d <- readRDS(here("data/processed/nhanes_merged.rds"))

# ---- 1. Analysis variables --------------------------------------------------
KNOTS <- c(4, 12)          # pre-specified inner knots for missing teeth (0-28)

tl_basis <- ns(d$missing_teeth, knots = KNOTS, Boundary.knots = c(0, 28))

d <- d %>%
  mutate(
    pir_cat_m   = fct_na_value_to_level(pir_cat, level = "Missing"),
    teeth5      = missing_teeth / 5,           # OR per 5 additional missing teeth
    dmft5       = dmft / 5,                    # OR per 5 DMFT units
    tl_score    = as.integer(tooth_loss_cat) - 1,  # 0,1,2,3 for P-trend
    lt20_f      = factor(lt20_teeth, 0:1, c("20+ teeth", "<20 teeth")),
    caries_f    = factor(any_caries, 0:1, c("No", "Yes")),
    log_crp     = log(LBXHSCRP),
    dentate     = as.integer(edentulous == 0),
    ts1 = tl_basis[, 1], ts2 = tl_basis[, 2], ts3 = tl_basis[, 3]
  )

des <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMECPRP,
                 nest = TRUE, data = d %>% filter(!is.na(WTMECPRP), WTMECPRP > 0))

pops <- list(
  "All adults"          = subset(des, analytic == 1),
  "Adults with MASLD"   = subset(des, analytic == 1 & masld == 1),
  "Dentate adults"      = subset(des, analytic == 1 & dentate == 1)
)

# ---- 2. Covariate sets ------------------------------------------------------
m1 <- "ns(RIDAGEYR, df = 3) + sex + race"
m2 <- paste(m1, "+ education + pir_cat_m + smoking + alcohol_g_wk + insured")
m3 <- paste(m2, "+ BMXBMI + diabetes")
m4 <- paste(m2, "+ log_crp")
models <- c(`Model 1` = m1, `Model 2` = m2, `Model 3` = m3, `Model 4` = m4)

# ---- 3. Which analyses to run ----------------------------------------------
# exposure, the population it is analysed in, and how to label it
exposures <- tribble(
  ~exposure,         ~label,                         ~pop_all, ~pop_dentate,
  "tooth_loss_cat",  "Tooth loss (ref: 0 missing)",  TRUE,     FALSE,
  "tl_score",        "P for trend (tooth loss)",     TRUE,     FALSE,
  "teeth5",          "Per 5 missing teeth",          TRUE,     FALSE,
  "lt20_f",          "Fewer than 20 teeth",          TRUE,     FALSE,
  "caries_f",        "Any untreated caries",         FALSE,    TRUE,
  "decayed_teeth",   "Per decayed tooth",            FALSE,    TRUE,
  "dmft5",           "Per 5 DMFT",                   FALSE,    TRUE
)

analyses <- bind_rows(
  # Aim 1: MASLD in all adults (tooth loss) and dentate adults (caries)
  exposures %>% filter(pop_all)     %>% mutate(outcome = "masld",    population = "All adults"),
  exposures %>% filter(pop_dentate) %>% mutate(outcome = "masld",    population = "Dentate adults"),
  # Aim 2: fibrosis in all adults, in adults with MASLD, and caries in dentate adults
  exposures %>% filter(pop_all)     %>% mutate(outcome = "fib_sig",  population = "All adults"),
  exposures %>% filter(pop_all)     %>% mutate(outcome = "fib_sig",  population = "Adults with MASLD"),
  exposures %>% filter(pop_dentate) %>% mutate(outcome = "fib_sig",  population = "Dentate adults"),
  exposures %>% filter(pop_all)     %>% mutate(outcome = "fib_cirr", population = "All adults")
)

outcome_labels <- c(masld = "MASLD", fib_sig = "Significant fibrosis (LSM >= 8 kPa)",
                    fib_cirr = "Cirrhosis range (LSM >= 13.6 kPa)")

# ---- 4. Model fitting -------------------------------------------------------
fit_one <- function(outcome, exposure, model_name, pop) {
  f <- as.formula(paste(outcome, "~", exposure, "+", models[[model_name]]))
  design <- pops[[pop]]
  fit <- tryCatch(svyglm(f, design = design, family = quasibinomial()),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)

  est <- coef(fit)
  ci  <- confint(fit)
  p   <- summary(fit)$coefficients[, "Pr(>|t|)"]
  keep <- grepl(paste0("^", exposure), names(est))

  tibble(outcome = outcome, population = pop, exposure = exposure,
         term = sub(paste0("^", exposure), "", names(est)[keep]),
         model = model_name,
         n = length(fit$residuals),
         OR = exp(est[keep]), lower = exp(ci[keep, 1]), upper = exp(ci[keep, 2]),
         p = p[keep])
}

results <- pmap_dfr(analyses, function(exposure, outcome, population, ...) {
  map_dfr(names(models), ~fit_one(outcome, exposure, .x, population))
})

results <- results %>%
  left_join(exposures %>% select(exposure, label), by = "exposure") %>%
  mutate(
    outcome_label = outcome_labels[outcome],
    term = if_else(term == "", label, term),
    OR_CI = sprintf("%.2f (%.2f-%.2f)", OR, lower, upper),
    p_fmt = if_else(p < 0.001, "<0.001", sprintf("%.3f", p))
  )

# P for trend rows report only the P value
results <- results %>%
  mutate(OR_CI = if_else(exposure == "tl_score", "", OR_CI))

write_csv(results, here("output/tables/main_models.csv"))

# ---- 5. Headline table: Model 2 --------------------------------------------
primary <- results %>%
  filter(model == "Model 2") %>%
  transmute(Outcome = outcome_label, Population = population,
            Exposure = label, Level = term, n, `OR (95% CI)` = OR_CI, P = p_fmt)
write_csv(primary, here("output/tables/main_models_primary.csv"))

cat("\n===== PRIMARY RESULTS (Model 2) =====\n")
print(primary, n = Inf, width = Inf)

# Compact view across all four models for the key exposures
cat("\n===== ROBUSTNESS ACROSS MODELS (all adults) =====\n")
results %>%
  filter(population == "All adults", exposure %in% c("tooth_loss_cat", "teeth5")) %>%
  select(outcome, term, model, OR_CI) %>%
  pivot_wider(names_from = model, values_from = OR_CI) %>%
  print(n = Inf, width = Inf)

# ---- 6. Dose-response splines (Model 2 covariates) --------------------------
spline_curve <- function(outcome, pop = "All adults") {
  f   <- as.formula(paste(outcome, "~ ts1 + ts2 + ts3 +", m2))
  fit <- svyglm(f, design = pops[[pop]], family = quasibinomial())
  b   <- coef(fit)[c("ts1", "ts2", "ts3")]
  V   <- vcov(fit)[c("ts1", "ts2", "ts3"), c("ts1", "ts2", "ts3")]

  grid <- 0:28
  X    <- predict(tl_basis, grid) - matrix(predict(tl_basis, 0), length(grid), 3, byrow = TRUE)
  lp   <- as.vector(X %*% b)
  se   <- sqrt(rowSums((X %*% V) * X))

  # Overall and nonlinearity tests
  p_overall <- regTermTest(fit, ~ ts1 + ts2 + ts3)$p
  fit_lin   <- svyglm(as.formula(paste(outcome, "~ missing_teeth +", m2)),
                      design = pops[[pop]], family = quasibinomial())
  tibble(missing_teeth = grid, OR = exp(lp),
         lower = exp(lp - 1.96 * se), upper = exp(lp + 1.96 * se),
         p_overall = as.numeric(p_overall),
         p_nonlinear = as.numeric(regTermTest(fit, ~ ts2 + ts3)$p))
}

plot_curve <- function(cv, title, file) {
  lab <- sprintf("P overall %s; P nonlinearity %s",
                 ifelse(cv$p_overall[1] < 0.001, "<0.001", sprintf("= %.3f", cv$p_overall[1])),
                 ifelse(cv$p_nonlinear[1] < 0.001, "<0.001", sprintf("= %.3f", cv$p_nonlinear[1])))
  g <- ggplot(cv, aes(missing_teeth, OR)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#2B6CB0", alpha = 0.18) +
    geom_line(colour = "#2B6CB0", linewidth = 1) +
    scale_x_continuous(breaks = seq(0, 28, 4)) +
    scale_y_log10() +
    labs(x = "Number of missing teeth (of 28)",
         y = "Adjusted odds ratio (log scale, ref: 0 missing)",
         title = title, subtitle = lab) +
    theme_minimal(base_size = 12)
  ggsave(here("output/figures", file), g, width = 7, height = 5, dpi = 300)
  invisible(g)
}

cv_masld <- spline_curve("masld")
cv_fib   <- spline_curve("fib_sig")
plot_curve(cv_masld, "Missing teeth and MASLD", "spline_masld.png")
plot_curve(cv_fib, "Missing teeth and significant fibrosis", "spline_fibrosis.png")

cat("\n===== DOSE-RESPONSE (spline, Model 2) =====\n")
bind_rows(MASLD = cv_masld, Fibrosis = cv_fib, .id = "outcome") %>%
  filter(missing_teeth %in% c(0, 4, 8, 12, 16, 20, 24, 28)) %>%
  mutate(OR_CI = sprintf("%.2f (%.2f-%.2f)", OR, lower, upper)) %>%
  select(outcome, missing_teeth, OR_CI) %>%
  pivot_wider(names_from = outcome, values_from = OR_CI) %>%
  print(n = Inf)
cat("MASLD:    P overall =", signif(cv_masld$p_overall[1], 3),
    "; P nonlinearity =", signif(cv_masld$p_nonlinear[1], 3), "\n")
cat("Fibrosis: P overall =", signif(cv_fib$p_overall[1], 3),
    "; P nonlinearity =", signif(cv_fib$p_nonlinear[1], 3), "\n")

message("\nScript 04 finished. Send back the PRIMARY RESULTS, ROBUSTNESS, and DOSE-RESPONSE blocks,",
        " and look at the two figures in output/figures.")
