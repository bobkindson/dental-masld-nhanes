# =============================================================================
# Dental status and MASLD / liver fibrosis — NHANES 2017–March 2020
# Script 05: sensitivity analyses, Aim 3 (sex interaction + discrimination),
#            and the main figure (fibrosis spline among adults with MASLD)
# -----------------------------------------------------------------------------
# Input : data/processed/nhanes_merged.rds   (created by Script 02)
# Output: output/tables/sensitivity.csv
#         output/tables/aim3_interaction.csv
#         output/tables/aim3_auc.csv
#         output/figures/spline_fibrosis_masld.png   (main figure)
# Runtime: about 5-10 minutes (multiple imputation is the slow part)
# =============================================================================

library(tidyverse)
library(survey)
library(splines)
library(mice)
library(mitools)
library(pROC)
library(here)

options(survey.lonely.psu = "adjust")
set.seed(20260927)

d <- readRDS(here("data/processed/nhanes_merged.rds"))

# ---- 1. Analysis variables (same as Script 04, plus sensitivity versions) ---
KNOTS    <- c(4, 12)
tl_basis <- ns(d$missing_teeth, knots = KNOTS, Boundary.knots = c(0, 28))

teeth    <- sprintf("%02d", c(2:15, 18:31))
ctc      <- as.matrix(d[paste0("OHX", teeth, "CTC")]); ctc[is.na(ctc)] <- ""
dent_ok  <- d$OHDDESTS %in% 1

masld_at <- function(cut) case_when(
  is.na(d$LUXCAPM)               ~ NA_integer_,
  d$LUXCAPM < cut                ~ 0L,
  d$cmrf_any == 1                ~ 1L,
  d$cmrf_any == 0                ~ 0L,
  TRUE                           ~ NA_integer_)

d <- d %>%
  mutate(
    pir_cat_m  = fct_na_value_to_level(pir_cat, level = "Missing"),
    teeth5     = missing_teeth / 5,
    lt20_f     = factor(lt20_teeth, 0:1, c("20+ teeth", "<20 teeth")),
    caries_f   = factor(any_caries, 0:1, c("No", "Yes")),
    caries_z_f = factor(if_else(dent_ok, as.integer(rowSums(ctc == "Z") > 0), NA_integer_),
                        0:1, c("No", "Yes")),
    dentate    = as.integer(edentulous == 0),
    fib_adv    = as.integer(LUXSMED >= 9.7),
    masld_263  = masld_at(263),
    masld_274  = masld_at(274),
    ts1 = tl_basis[, 1], ts2 = tl_basis[, 2], ts3 = tl_basis[, 3],
    age60      = factor(if_else(RIDAGEYR >= 60, "60+", "20-59"), levels = c("20-59", "60+")),
    alc_g0      = replace_na(alcohol_g_wk, 0),          # for S2: unknown intake set to 0 ...
    alc_unknown = as.integer(is.na(alcohol_g_wk)),      # ... and flagged with an indicator

    # Every eligibility step except the final "MASLD classifiable" step
    base_ok = RIDAGEYR >= 20 & RIDSTATR == 2 & lux_complete %in% 1 & lux_reliable %in% 1 &
              dent_ok & pregnant == 0 & viral_hep == 0 & self_hep == 0,
    analytic_alc     = as.integer(base_ok & !(heavy_alcohol %in% 1) & !is.na(masld)),
    analytic_263     = as.integer(base_ok & heavy_alcohol %in% 0 & !is.na(masld_263)),
    analytic_274     = as.integer(base_ok & heavy_alcohol %in% 0 & !is.na(masld_274)),
    analytic_dentate = as.integer(analytic == 1 & dentate == 1),
    analytic_nocirr  = as.integer(analytic == 1 & LUXSMED < 13.6)
  )

des   <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMECPRP,
                   nest = TRUE, data = d %>% filter(!is.na(WTMECPRP), WTMECPRP > 0))
des_f <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTSAFPRP,
                   nest = TRUE, data = d %>% filter(!is.na(WTSAFPRP), WTSAFPRP > 0))

# ---- 2. Covariate sets ------------------------------------------------------
age_term <- "ns(RIDAGEYR, df = 3)"
m2_parts <- c(age_term, "sex", "race", "education", "pir_cat_m", "smoking",
              "alcohol_g_wk", "insured")
make_cov <- function(parts) paste(parts, collapse = " + ")
m2       <- make_cov(m2_parts)
m2_cot   <- make_cov(replace(m2_parts, m2_parts == "smoking", "log_cotinine"))
m2_alc   <- make_cov(replace(m2_parts, m2_parts == "alcohol_g_wk", "alc_g0 + alc_unknown"))
m2_nosex <- make_cov(setdiff(m2_parts, "sex"))

# ---- 3. Helpers -------------------------------------------------------------
fit_or <- function(outcome, exposure, covars, design) {
  f   <- as.formula(paste(outcome, "~", exposure, "+", covars))
  fit <- tryCatch(svyglm(f, design = design, family = quasibinomial()),
                  error = function(e) NULL)
  if (is.null(fit)) return(tibble(n = NA, OR = NA, lower = NA, upper = NA, p = NA))
  est  <- coef(fit); ci <- confint(fit)
  pv   <- summary(fit)$coefficients[, "Pr(>|t|)"]
  k    <- grep(paste0("^", exposure), names(est))[1]   # binary/continuous: one term
  tibble(n = length(fit$residuals), OR = exp(est[k]),
         lower = exp(ci[k, 1]), upper = exp(ci[k, 2]), p = pv[k])
}

fmt_or <- function(OR, lower, upper) sprintf("%.2f (%.2f-%.2f)", OR, lower, upper)
fmt_p  <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))

# The six "core" associations re-estimated in every sensitivity analysis
run_core <- function(analysis, design, flag = "analytic", covars = m2,
                     masld_var = "masld", fib_var = "fib_sig", caries_var = "caries_f") {
  v     <- design$variables
  all_d <- subset(design, v[[flag]] %in% 1)
  mas_d <- subset(design, v[[flag]] %in% 1 & v[[masld_var]] %in% 1)
  den_d <- subset(design, v[[flag]] %in% 1 & v$dentate %in% 1)

  specs <- tribble(
    ~outcome,  ~exposure,   ~population,         ~dsn,
    masld_var, "teeth5",    "All adults",        "all",
    masld_var, caries_var,  "Dentate adults",    "den",
    fib_var,   "teeth5",    "All adults",        "all",
    fib_var,   "teeth5",    "Adults with MASLD", "mas",
    fib_var,   "lt20_f",    "Adults with MASLD", "mas",
    fib_var,   caries_var,  "Dentate adults",    "den"
  )
  dsns <- list(all = all_d, mas = mas_d, den = den_d)

  pmap_dfr(specs, function(outcome, exposure, population, dsn) {
    r <- fit_or(outcome, exposure, covars, dsns[[dsn]])
    tibble(analysis, outcome = if_else(grepl("^masld", outcome), "MASLD", "Fibrosis"),
           exposure = case_when(exposure == "teeth5" ~ "Per 5 missing teeth",
                                exposure == "lt20_f" ~ "Fewer than 20 teeth",
                                TRUE ~ "Any untreated caries"),
           population) %>% bind_cols(r)
  })
}

# ---- 4. Sensitivity analyses ------------------------------------------------
cat("Running sensitivity analyses...\n")
sens <- bind_rows(
  run_core("Primary (Model 2)", des),
  run_core("S1: Caries = coronal lesions only (code Z)", des, caries_var = "caries_z_f"),
  run_core("S2: Unknown alcohol intake kept in", des, flag = "analytic_alc", covars = m2_alc),
  run_core("S3: Steatosis CAP >= 263 dB/m", des, flag = "analytic_263", masld_var = "masld_263"),
  run_core("S4: Steatosis CAP >= 274 dB/m", des, flag = "analytic_274", masld_var = "masld_274"),
  run_core("S5: Serum cotinine instead of self-reported smoking", des, covars = m2_cot),
  run_core("S6: Fasting subsample, fasting weights", des_f),
  run_core("S7: Advanced fibrosis (LSM >= 9.7 kPa)", des, fib_var = "fib_adv"),
  run_core("S8: Dentate adults only (edentulous excluded)", des, flag = "analytic_dentate"),
  run_core("S10: Cirrhosis-range stiffness excluded (LSM >= 13.6 kPa)", des, flag = "analytic_nocirr")
)

# S9: multiple imputation for income, education, insurance and smoking --------
cat("Running multiple imputation (20 datasets)...\n")
imp_vars <- c("pir_cat", "education", "insured", "smoking")
pred_vars <- c("RIDAGEYR", "sex", "race", "alcohol_g_wk", "masld", "fib_sig",
               "missing_teeth", "any_caries", "BMXBMI", "diabetes")
a_rows <- which(d$analytic == 1)
imp <- mice(d[a_rows, c(imp_vars, pred_vars)], m = 20, maxit = 10, printFlag = FALSE)

full <- d %>% filter(!is.na(WTMECPRP), WTMECPRP > 0)
pos  <- match(d$SEQN[a_rows], full$SEQN)
imp_list <- lapply(seq_len(imp$m), function(i) {
  x <- full; ci <- complete(imp, i)
  for (v in imp_vars) x[[v]][pos] <- ci[[v]]
  x
})
des_mi <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMECPRP,
                    nest = TRUE, data = imputationList(imp_list))
m2_mi <- sub("pir_cat_m", "pir_cat", m2)

mi_fit <- function(outcome, exposure, pop_expr) {
  sub_mi <- subset(des_mi, eval(pop_expr))
  fits <- with(sub_mi, svyglm(as.formula(paste(outcome, "~", exposure, "+", m2_mi)),
                              family = quasibinomial()))
  comb <- MIcombine(fits)
  k  <- grep(paste0("^", exposure), names(coef(comb)))[1]
  se <- sqrt(diag(vcov(comb)))[k]; b <- coef(comb)[k]
  tibble(n = sum(eval(pop_expr, full) %in% TRUE), OR = exp(b), lower = exp(b - 1.96 * se),
         upper = exp(b + 1.96 * se), p = 2 * pnorm(-abs(b / se)))
}
mi_rows <- bind_rows(
  tibble(outcome = "MASLD", exposure = "Per 5 missing teeth", population = "All adults") %>%
    bind_cols(mi_fit("masld", "teeth5", quote(analytic == 1))),
  tibble(outcome = "MASLD", exposure = "Any untreated caries", population = "Dentate adults") %>%
    bind_cols(mi_fit("masld", "caries_f", quote(analytic == 1 & edentulous == 0))),
  tibble(outcome = "Fibrosis", exposure = "Per 5 missing teeth", population = "All adults") %>%
    bind_cols(mi_fit("fib_sig", "teeth5", quote(analytic == 1))),
  tibble(outcome = "Fibrosis", exposure = "Per 5 missing teeth", population = "Adults with MASLD") %>%
    bind_cols(mi_fit("fib_sig", "teeth5", quote(analytic == 1 & masld == 1))),
  tibble(outcome = "Fibrosis", exposure = "Fewer than 20 teeth", population = "Adults with MASLD") %>%
    bind_cols(mi_fit("fib_sig", "lt20_f", quote(analytic == 1 & masld == 1))),
  tibble(outcome = "Fibrosis", exposure = "Any untreated caries", population = "Dentate adults") %>%
    bind_cols(mi_fit("fib_sig", "caries_f", quote(analytic == 1 & edentulous == 0)))
) %>% mutate(analysis = "S9: Multiple imputation of missing covariates (m = 20)")

sens <- bind_rows(sens, mi_rows) %>%
  mutate(OR_CI = fmt_or(OR, lower, upper), P = fmt_p(p))
write_csv(sens, here("output/tables/sensitivity.csv"))

cat("\n===== SENSITIVITY ANALYSES: OR (95% CI) =====\n")
sens %>%
  mutate(key = paste(outcome, "|", exposure, "|", population)) %>%
  select(analysis, key, OR_CI) %>%
  pivot_wider(names_from = key, values_from = OR_CI) %>%
  pivot_longer(-analysis, names_to = "association", values_to = "OR_CI") %>%
  pivot_wider(names_from = analysis, values_from = OR_CI) %>%
  print(n = Inf, width = Inf)

# ---- 5. Aim 3a: effect modification by sex (and age, exploratory) -----------
cat("\nTesting effect modification...\n")
all_d <- subset(des, analytic == 1)
mas_d <- subset(des, analytic == 1 & masld == 1)
den_d <- subset(des, analytic == 1 & dentate == 1)

int_specs <- tribble(
  ~outcome,  ~exposure,  ~population,          ~dsn,
  "fib_sig", "teeth5",   "Adults with MASLD",  "mas",
  "fib_sig", "lt20_f",   "Adults with MASLD",  "mas",
  "masld",   "caries_f", "Dentate adults",     "den",
  "masld",   "teeth5",   "All adults",         "all"
)
dsns <- list(all = all_d, mas = mas_d, den = den_d)

interaction_tbl <- pmap_dfr(int_specs, function(outcome, exposure, population, dsn) {
  dsn_obj <- dsns[[dsn]]
  map_dfr(c("sex", "age60"), function(mod) {
    covars <- if (mod == "sex") m2_nosex else m2
    f_int <- as.formula(paste(outcome, "~", exposure, "*", mod, "+", covars))
    fit   <- svyglm(f_int, design = dsn_obj, family = quasibinomial())
    p_int <- regTermTest(fit, as.formula(paste0("~", exposure, ":", mod)))$p
    lv    <- levels(dsn_obj$variables[[mod]])
    strata <- map_dfr(lv, function(l) {
      s <- subset(dsn_obj, dsn_obj$variables[[mod]] == l)
      r <- fit_or(outcome, exposure, if (mod == "sex") m2_nosex else m2, s)
      tibble(stratum = l, OR_CI = fmt_or(r$OR, r$lower, r$upper), n = r$n)
    })
    strata %>% mutate(outcome, exposure, population, modifier = mod,
                      P_interaction = fmt_p(as.numeric(p_int)))
  })
})
write_csv(interaction_tbl, here("output/tables/aim3_interaction.csv"))
cat("\n===== AIM 3a: EFFECT MODIFICATION (Model 2) =====\n")
print(interaction_tbl %>% select(outcome, exposure, population, modifier, stratum,
                                 n, OR_CI, P_interaction), n = Inf, width = Inf)

# ---- 6. Aim 3b: does tooth count add to fibrosis identification in MASLD? ---
cat("\nEstimating discrimination (AUC)...\n")
auc_vars <- c("fib_sig", "RIDAGEYR", "sex", "race", "BMXBMI", "diabetes",
              "missing_teeth", "WTMECPRP")
mas_cc_flag <- with(d, analytic == 1 & masld == 1 &
                      complete.cases(d[, auc_vars]))
des_auc <- subset(des, mas_cc_flag[match(des$variables$SEQN, d$SEQN)])
dat_auc <- d[mas_cc_flag, ]

base_f <- "fib_sig ~ ns(RIDAGEYR, df = 3) + sex + race + BMXBMI + diabetes"
fit_base  <- svyglm(as.formula(base_f), design = des_auc, family = quasibinomial())
fit_teeth <- svyglm(as.formula(paste(base_f, "+ missing_teeth")),
                    design = des_auc, family = quasibinomial())
p_base  <- as.numeric(predict(fit_base,  newdata = dat_auc, type = "response"))
p_teeth <- as.numeric(predict(fit_teeth, newdata = dat_auc, type = "response"))

weighted_auc <- function(pred, y, w) {
  g  <- factor(pred, levels = sort(unique(pred)))
  wc <- tapply(w * y, g, sum); wn <- tapply(w * (1 - y), g, sum)
  wc[is.na(wc)] <- 0; wn[is.na(wn)] <- 0
  below <- cumsum(wn) - wn
  sum(wc * (below + 0.5 * wn)) / (sum(wc) * sum(wn))
}
y <- dat_auc$fib_sig; w <- dat_auc$WTMECPRP

roc_b <- roc(y, p_base, quiet = TRUE); roc_t <- roc(y, p_teeth, quiet = TRUE)
dl    <- roc.test(roc_b, roc_t, method = "delong")

auc_tbl <- tibble(
  model = c("Clinical: age, sex, race, BMI, diabetes", "Clinical + number of missing teeth"),
  n = nrow(dat_auc), events = sum(y),
  weighted_AUC   = round(c(weighted_auc(p_base, y, w), weighted_auc(p_teeth, y, w)), 3),
  unweighted_AUC = round(c(auc(roc_b), auc(roc_t)), 3),
  AUC_95CI_unweighted = c(paste(round(ci.auc(roc_b)[c(1, 3)], 3), collapse = "-"),
                          paste(round(ci.auc(roc_t)[c(1, 3)], 3), collapse = "-")),
  P_DeLong = c("", fmt_p(dl$p.value))
)
write_csv(auc_tbl, here("output/tables/aim3_auc.csv"))
cat("\n===== AIM 3b: DISCRIMINATION FOR SIGNIFICANT FIBROSIS IN MASLD =====\n")
print(auc_tbl, width = Inf)

# ---- 7. Main figure: fibrosis spline among adults with MASLD ----------------
fit_sp <- svyglm(as.formula(paste("fib_sig ~ ts1 + ts2 + ts3 +", m2)),
                 design = mas_d, family = quasibinomial())
b <- coef(fit_sp)[c("ts1", "ts2", "ts3")]
V <- vcov(fit_sp)[c("ts1", "ts2", "ts3"), c("ts1", "ts2", "ts3")]
grid <- 0:28
X  <- predict(tl_basis, grid) - matrix(predict(tl_basis, 0), length(grid), 3, byrow = TRUE)
lp <- as.vector(X %*% b); se <- sqrt(rowSums((X %*% V) * X))
p_all <- as.numeric(regTermTest(fit_sp, ~ ts1 + ts2 + ts3)$p)
p_nl  <- as.numeric(regTermTest(fit_sp, ~ ts2 + ts3)$p)
cv <- tibble(missing_teeth = grid, OR = exp(lp),
             lower = exp(lp - 1.96 * se), upper = exp(lp + 1.96 * se))

# Distribution rug: where the MASLD participants actually are
rug_dat <- d %>% filter(analytic == 1, masld == 1) %>% count(missing_teeth)

g <- ggplot(cv, aes(missing_teeth, OR)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#2B6CB0", alpha = 0.18) +
  geom_line(colour = "#2B6CB0", linewidth = 1.1) +
  geom_rug(data = rug_dat, aes(x = missing_teeth), inherit.aes = FALSE,
           sides = "b", alpha = 0.4, colour = "grey30") +
  scale_x_continuous(breaks = seq(0, 28, 4)) +
  scale_y_log10() +
  labs(x = "Number of missing teeth (of 28, third molars excluded)",
       y = "Adjusted odds ratio of significant fibrosis\n(log scale; reference: 0 missing teeth)",
       title = "Missing teeth and significant liver fibrosis among US adults with MASLD",
       subtitle = sprintf("NHANES 2017-March 2020; P overall %s; P nonlinearity %s",
                          ifelse(p_all < 0.001, "<0.001", sprintf("= %.3f", p_all)),
                          ifelse(p_nl < 0.001, "<0.001", sprintf("= %.3f", p_nl))),
       caption = paste("Adjusted for age (spline), sex, race/ethnicity, education, income,",
                       "smoking, alcohol and insurance. Shaded band: 95% CI.")) +
  theme_minimal(base_size = 12) +
  theme(plot.caption = element_text(hjust = 0, colour = "grey30"))
ggsave(here("output/figures/spline_fibrosis_masld.png"), g, width = 8, height = 5.5, dpi = 300)

cat("\n===== MAIN FIGURE: fibrosis spline in MASLD =====\n")
cv %>% filter(missing_teeth %in% c(0, 4, 8, 12, 16, 20, 24, 28)) %>%
  mutate(OR_CI = fmt_or(OR, lower, upper)) %>% select(missing_teeth, OR_CI) %>% print(n = Inf)
cat("P overall =", signif(p_all, 3), "; P nonlinearity =", signif(p_nl, 3), "\n")

message("\nScript 05 finished. Attach sensitivity.csv, aim3_interaction.csv, aim3_auc.csv",
        " and spline_fibrosis_masld.png.")
