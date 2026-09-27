# =============================================================================
# Dental status and MASLD / liver fibrosis — NHANES 2017–March 2020
# Script 03: Table 1 (weighted characteristics) + covariate missingness
# -----------------------------------------------------------------------------
# Input : data/processed/nhanes_merged.rds   (created by Script 02)
# Output: output/tables/table1_by_toothloss.csv   (main Table 1)
#         output/tables/table1_by_masld.csv       (supplementary table)
#         output/tables/missingness_analytic.csv
#
# Format: continuous = weighted mean (SE); categorical = unweighted n (weighted %).
# P values: design-based F test (continuous), Rao-Scott chi-square (categorical).
# =============================================================================

library(tidyverse)
library(survey)
library(here)

options(survey.lonely.psu = "adjust")

d <- readRDS(here("data/processed/nhanes_merged.rds"))

# ---- 1. Prepare display variables -------------------------------------------
yn <- function(x) factor(if_else(x == 1, "Yes", "No"), levels = c("No", "Yes"))

d <- d %>%
  mutate(
    age_cat      = cut(RIDAGEYR, c(19, 39, 59, Inf), labels = c("20-39", "40-59", "60+")),
    diabetes_f   = yn(diabetes),
    masld_f      = factor(if_else(masld == 1, "MASLD", "No MASLD"),
                          levels = c("No MASLD", "MASLD")),
    fib_sig_f    = yn(fib_sig),
    fib_cirr_f   = yn(fib_cirr),
    any_caries_f = yn(any_caries),
    lt20_f       = yn(lt20_teeth),
    edent_f      = yn(edentulous)
  )

des <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMECPRP,
                 nest = TRUE, data = d %>% filter(!is.na(WTMECPRP), WTMECPRP > 0))
sub <- subset(des, analytic == 1)

# ---- 2. Rows of the table ---------------------------------------------------
# type: "cont" = mean (SE); "cat" = all levels; "bin" = show the "Yes" level only
rows <- tribble(
  ~var,            ~label,                               ~type,
  "RIDAGEYR",      "Age, years",                         "cont",
  "age_cat",       "Age group",                          "cat",
  "sex",           "Sex",                                "cat",
  "race",          "Race/ethnicity",                     "cat",
  "education",     "Education",                          "cat",
  "pir_cat",       "Family income-to-poverty ratio",     "cat",
  "insured",       "Health insurance",                   "cat",
  "smoking",       "Smoking status",                     "cat",
  "alcohol_g_wk",  "Alcohol, g/week",                    "cont",
  "BMXBMI",        "BMI, kg/m2",                         "cont",
  "BMXWAIST",      "Waist circumference, cm",            "cont",
  "sbp",           "Systolic blood pressure, mmHg",      "cont",
  "LBXGH",         "HbA1c, %",                           "cont",
  "LBDHDD",        "HDL cholesterol, mg/dL",             "cont",
  "LBXHSCRP",      "hs-CRP, mg/L",                       "cont",
  "diabetes_f",    "Diabetes",                           "bin",
  "missing_teeth", "Missing teeth (of 28)",              "cont",
  "lt20_f",        "Fewer than 20 teeth",                "bin",
  "edent_f",       "Edentulous",                         "bin",
  "decayed_teeth", "Decayed teeth",                      "cont",
  "any_caries_f",  "Any untreated caries",               "bin",
  "dmft",          "DMFT",                               "cont",
  "LUXCAPM",       "CAP, dB/m",                          "cont",
  "LUXSMED",       "Liver stiffness, kPa",               "cont",
  "masld_f",       "MASLD",                              "cat",
  "fib_sig_f",     "Significant fibrosis (LSM >= 8 kPa)", "bin",
  "fib_cirr_f",    "Cirrhosis range (LSM >= 13.6 kPa)",  "bin"
)

# ---- 3. Table builder -------------------------------------------------------
fmt_mean <- function(est, se) sprintf("%.1f (%.1f)", est, se)
fmt_pct  <- function(n, p)    sprintf("%s (%.1f)", format(n, big.mark = ","), 100 * p)
fmt_p    <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))

build_table1 <- function(design, group) {
  groups  <- levels(design$variables[[group]])
  designs <- c(list(Overall = design),
               setNames(lapply(groups, function(g)
                 subset(design, design$variables[[group]] == g)), groups))
  out <- list()

  # First row: unweighted sample size per column
  out[[1]] <- tibble(Characteristic = "Participants, n",
                     !!!setNames(lapply(designs, function(x) format(nrow(x), big.mark = ",")),
                                 names(designs)),
                     P = "")

  for (i in seq_len(nrow(rows))) {
    v <- rows$var[i]; lab <- rows$label[i]; type <- rows$type[i]
    if (v == group) next
    f <- as.formula(paste0("~", v))

    if (type == "cont") {
      cells <- sapply(designs, function(x) {
        m <- svymean(f, x, na.rm = TRUE)
        fmt_mean(coef(m), SE(m))
      })
      fit <- svyglm(as.formula(paste(v, "~", group)), design = design)
      p   <- tryCatch(regTermTest(fit, as.formula(paste0("~", group)))$p,
                      error = function(e) NA)
      out[[length(out) + 1]] <- tibble(Characteristic = paste0(lab, ", mean (SE)"),
                                       !!!as.list(cells), P = fmt_p(as.numeric(p)))
    } else {
      p <- tryCatch(svychisq(as.formula(paste0("~", v, "+", group)), design)$p.value,
                    error = function(e) NA)
      lv <- levels(design$variables[[v]])
      show_lv <- if (type == "bin") "Yes" else lv
      if (type == "cat")
        out[[length(out) + 1]] <- tibble(Characteristic = lab, P = fmt_p(p))
      for (l in show_lv) {
        cells <- sapply(designs, function(x) {
          vals <- x$variables[[v]]
          n_l  <- sum(vals == l, na.rm = TRUE)
          pr   <- svymean(f, x, na.rm = TRUE)
          fmt_pct(n_l, coef(pr)[paste0(v, l)])
        })
        name <- if (type == "bin") paste0(lab, ", n (%)") else paste0("   ", l)
        out[[length(out) + 1]] <- tibble(Characteristic = name, !!!as.list(cells),
                                         P = if (type == "bin") fmt_p(p) else "")
      }
    }
  }
  bind_rows(out) %>% mutate(across(everything(), ~replace_na(.x, "")))
}

# ---- 4. Build and save ------------------------------------------------------
t1 <- build_table1(sub, "tooth_loss_cat")
t1_masld <- build_table1(sub, "masld_f")

write_csv(t1, here("output/tables/table1_by_toothloss.csv"))
write_csv(t1_masld, here("output/tables/table1_by_masld.csv"))

cat("\n===== TABLE 1: by number of missing teeth =====\n")
print(t1, n = Inf, width = Inf)

# ---- 5. Missingness of candidate covariates (unweighted, analytic sample) ---
covars <- c("RIDAGEYR", "sex", "race", "education", "pir_cat", "insured", "smoking",
            "log_cotinine", "BMXBMI", "BMXWAIST", "diabetes", "LBXGH", "LBDHDD",
            "LBXTR", "LBXHSCRP", "alcohol_g_wk")
miss <- d %>%
  filter(analytic == 1) %>%
  summarise(across(all_of(covars), ~sum(is.na(.x)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing") %>%
  mutate(pct_missing = round(100 * n_missing / sum(d$analytic == 1), 1)) %>%
  arrange(desc(n_missing))

write_csv(miss, here("output/tables/missingness_analytic.csv"))
cat("\n===== MISSING DATA IN THE ANALYTIC SAMPLE =====\n")
print(miss, n = Inf)

message("\nScript 03 finished. Send back TABLE 1 and the MISSING DATA table.")
