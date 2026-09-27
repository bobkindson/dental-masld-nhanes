# =============================================================================
# Dental status and MASLD / liver fibrosis — NHANES 2017–March 2020
# Script 02: merge files, derive variables, apply exclusions, flowchart counts
# -----------------------------------------------------------------------------
# Input : data/raw/*.rds          (created by Script 01)
# Output: data/processed/nhanes_merged.rds   (ALL 15,560 participants + flags)
#         output/tables/flowchart_counts.csv
#         output/tables/dental_code_check.csv
#
# IMPORTANT: we keep every participant in the saved file and mark the analytic
# sample with a flag (analytic == 1). NHANES variance estimation must use the
# full design and subset() it later — never drop rows before building the design.
# =============================================================================

library(tidyverse)
library(haven)
library(survey)
library(here)

options(survey.lonely.psu = "adjust")

# ---- 0. Pre-specified analytic decisions (review before running) ------------
CAP_CUT      <- 285    # dB/m, steatosis threshold (matches Ogbu et al. 2026)
LSM_SIG      <- 8.0    # kPa, clinically significant fibrosis
LSM_CIRR     <- 13.6   # kPa, cirrhosis range (secondary outcome; NOT excluded)
MIN_VALID    <- 10     # minimum valid elastography measurements
MAX_IQR_MED  <- 30     # maximum IQR/median (%) for a reliable exam

# Clinical coding choices — for the dentist on the team to confirm:
COUNT_ROOT_FRAGMENTS_AS_PRESENT <- FALSE  # tooth count code 5 counted as missing
COUNT_ROOT_TIPS_AS_DECAYED      <- TRUE   # caries code J counted as decayed (WHO)
EXCLUDE_MISSING_ALCOHOL         <- TRUE   # unknown alcohol intake -> excluded

# ---- 1. Load files ----------------------------------------------------------
rd <- function(f) readRDS(here("data/raw", paste0(f, ".rds"))) %>% zap_labels()

pick <- function(df, vars, file) {
  missing_vars <- setdiff(vars, names(df))
  if (length(missing_vars) > 0)
    warning(file, " is missing expected variables: ",
            paste(missing_vars, collapse = ", "), call. = FALSE)
  df %>% select(SEQN, any_of(vars))
}

teeth   <- sprintf("%02d", c(2:15, 18:31))       # 28 teeth, third molars excluded
tc_vars <- paste0("OHX", teeth, "TC")            # tooth count
ctc_vars <- paste0("OHX", teeth, "CTC")          # coronal caries tooth code

demo <- pick(rd("P_DEMO"), c("RIDSTATR", "RIAGENDR", "RIDAGEYR", "RIDRETH3",
                             "DMDEDUC2", "INDFMPIR", "RIDEXPRG",
                             "WTMECPRP", "SDMVPSU", "SDMVSTRA"), "P_DEMO")
lux    <- pick(rd("P_LUX"), c("LUAXSTAT", "LUANMVGP", "LUXSMED",
                              "LUXSIQRM", "LUXCAPM"), "P_LUX")
ohx    <- pick(rd("P_OHXDEN"), c("OHDDESTS", tc_vars, ctc_vars), "P_OHXDEN")
bmx    <- pick(rd("P_BMX"), c("BMXBMI", "BMXWAIST"), "P_BMX")
bpxo   <- pick(rd("P_BPXO"), c(paste0("BPXOSY", 1:3), paste0("BPXODI", 1:3)), "P_BPXO")
bpq    <- pick(rd("P_BPQ"), c("BPQ050A", "BPQ100D"), "P_BPQ")
glu    <- pick(rd("P_GLU"), c("LBXGLU", "WTSAFPRP"), "P_GLU")
ghb    <- pick(rd("P_GHB"), c("LBXGH"), "P_GHB")
diq    <- pick(rd("P_DIQ"), c("DIQ010", "DIQ050", "DIQ070"), "P_DIQ")
trig   <- pick(rd("P_TRIGLY"), c("LBXTR"), "P_TRIGLY")
hdl    <- pick(rd("P_HDL"), c("LBDHDD"), "P_HDL")
hepb   <- pick(rd("P_HEPBD"), c("LBXHBC", "LBDHBG"), "P_HEPBD")
hepc   <- pick(rd("P_HEPC"), c("LBXHCR", "LBDHCI"), "P_HEPC")
alq    <- pick(rd("P_ALQ"), c("ALQ111", "ALQ121", "ALQ130"), "P_ALQ")
mcq    <- pick(rd("P_MCQ"), c("MCQ160L", "MCQ510D", "MCQ510E"), "P_MCQ")
smq    <- pick(rd("P_SMQ"), c("SMQ020", "SMQ040"), "P_SMQ")
cot    <- pick(rd("P_COT"), c("LBXCOT"), "P_COT")
hiq    <- pick(rd("P_HIQ"), c("HIQ011"), "P_HIQ")
crp    <- pick(rd("P_HSCRP"), c("LBXHSCRP"), "P_HSCRP")

d <- list(demo, lux, ohx, bmx, bpxo, bpq, glu, ghb, diq, trig, hdl,
          hepb, hepc, alq, mcq, smq, cot, hiq, crp) %>%
  reduce(left_join, by = "SEQN")

stopifnot(nrow(d) == 15560, !anyDuplicated(d$SEQN))

# Helper: criterion met if ANY component is TRUE; 0 if at least one component
# is observed and none is TRUE; NA if every component is missing.
any_true <- function(...) {
  m <- cbind(...)
  ifelse(rowSums(m == 1, na.rm = TRUE) > 0, 1L,
         ifelse(rowSums(!is.na(m)) > 0, 0L, NA_integer_))
}

# ---- 2. Dental exposures ----------------------------------------------------
tc  <- as.matrix(d[tc_vars])
ctc <- as.matrix(d[ctc_vars])
ctc[is.na(ctc)] <- ""

present_codes <- if (COUNT_ROOT_FRAGMENTS_AS_PRESENT) c(1, 2, 5) else c(1, 2)
decay_codes   <- if (COUNT_ROOT_TIPS_AS_DECAYED) c("Z", "J") else "Z"

dent_ok <- d$OHDDESTS == 1 & !is.na(d$OHDDESTS)

d <- d %>%
  mutate(
    teeth_present  = if_else(dent_ok, rowSums(matrix(tc %in% present_codes, nrow = nrow(tc))), NA_real_),
    missing_teeth  = 28 - teeth_present,
    edentulous     = as.integer(teeth_present == 0),
    lt20_teeth     = as.integer(teeth_present < 20),
    tooth_loss_cat = cut(missing_teeth, breaks = c(-Inf, 0, 8, 27, 28),
                         labels = c("0", "1-8", "9-27", "Edentulous")),
    decayed_teeth  = if_else(dent_ok, rowSums(matrix(ctc %in% decay_codes, nrow = nrow(ctc))), NA_real_),
    missing_dis    = if_else(dent_ok, rowSums(matrix(ctc %in% c("E", "P", "R"), nrow = nrow(ctc))), NA_real_),
    filled_teeth   = if_else(dent_ok, rowSums(matrix(ctc == "F", nrow = nrow(ctc))), NA_real_),
    dmft           = decayed_teeth + missing_dis + filled_teeth,
    any_caries     = as.integer(decayed_teeth > 0)
  )

# ---- 3. Liver outcomes ------------------------------------------------------
d <- d %>%
  mutate(
    lux_complete = as.integer(LUAXSTAT == 1 & !is.na(LUXCAPM) & !is.na(LUXSMED)),
    lux_reliable = as.integer(lux_complete == 1 & LUANMVGP >= MIN_VALID &
                                LUXSIQRM < MAX_IQR_MED),
    steatosis    = as.integer(LUXCAPM >= CAP_CUT),
    fib_sig      = as.integer(LUXSMED >= LSM_SIG),
    fib_cirr     = as.integer(LUXSMED >= LSM_CIRR)
  )

# ---- 4. Cardiometabolic risk factors (2023 MASLD consensus) -----------------
d <- d %>%
  mutate(
    male   = as.integer(RIAGENDR == 1),
    asian  = as.integer(RIDRETH3 == 6),
    sbp    = rowMeans(across(starts_with("BPXOSY")), na.rm = TRUE),
    dbp    = rowMeans(across(starts_with("BPXODI")), na.rm = TRUE),
    sbp    = if_else(is.nan(sbp), NA_real_, sbp),
    dbp    = if_else(is.nan(dbp), NA_real_, dbp),
    lipid_med = as.integer(BPQ100D == 1),

    cmrf_adipo = any_true(BMXBMI >= if_else(asian == 1, 23, 25),
                          BMXWAIST > if_else(male == 1, 94, 80)),
    cmrf_glyc  = any_true(LBXGLU >= 100, LBXGH >= 5.7,
                          DIQ010 == 1, DIQ050 == 1, DIQ070 == 1),
    cmrf_bp    = any_true(sbp >= 130, dbp >= 85, BPQ050A == 1),
    cmrf_tg    = any_true(LBXTR >= 150, lipid_med == 1),
    cmrf_hdl   = any_true(LBDHDD < if_else(male == 1, 40, 50), lipid_med == 1),
    cmrf_any   = any_true(cmrf_adipo, cmrf_glyc, cmrf_bp, cmrf_tg, cmrf_hdl),

    masld = case_when(
      steatosis == 0                  ~ 0L,
      steatosis == 1 & cmrf_any == 1  ~ 1L,
      steatosis == 1 & cmrf_any == 0  ~ 0L,   # steatosis without CMRF (cryptogenic)
      TRUE                            ~ NA_integer_
    ),
    cryptogenic_sld = as.integer(steatosis == 1 & cmrf_any == 0)
  )

# ---- 5. Exclusion variables -------------------------------------------------
freq_per_year <- c(`0` = 0, `1` = 365, `2` = 286, `3` = 182, `4` = 104, `5` = 52,
                   `6` = 30, `7` = 12, `8` = 9, `9` = 4.5, `10` = 1.5)

d <- d %>%
  mutate(
    pregnant   = as.integer(RIDEXPRG == 1 & !is.na(RIDEXPRG)),
    viral_hep  = as.integer((LBDHBG == 1 & !is.na(LBDHBG)) |
                              (LBXHCR == 1 & !is.na(LBXHCR))),
    self_hep   = as.integer(!is.na(MCQ510D) | !is.na(MCQ510E)),

    alq_freq   = unname(freq_per_year[as.character(ALQ121)]),
    alq_drinks = if_else(ALQ130 %in% c(777, 999), NA_real_, ALQ130),
    drinks_wk  = case_when(
      ALQ111 == 2               ~ 0,    # never had a drink
      ALQ121 == 0               ~ 0,    # none in past 12 months
      !is.na(alq_freq) & !is.na(alq_drinks) ~ alq_freq * alq_drinks / 52,
      TRUE                      ~ NA_real_
    ),
    alcohol_g_wk = drinks_wk * 14,
    heavy_alcohol = case_when(
      is.na(alcohol_g_wk)             ~ NA_integer_,
      male == 1 & alcohol_g_wk > 210  ~ 1L,
      male == 0 & alcohol_g_wk > 140  ~ 1L,
      TRUE                            ~ 0L
    )
  )

# ---- 6. Covariates ----------------------------------------------------------
d <- d %>%
  mutate(
    sex = factor(RIAGENDR, 1:2, c("Male", "Female")),
    race = factor(case_when(
      RIDRETH3 %in% c(1, 2) ~ "Hispanic",
      RIDRETH3 == 3 ~ "Non-Hispanic White",
      RIDRETH3 == 4 ~ "Non-Hispanic Black",
      RIDRETH3 == 6 ~ "Non-Hispanic Asian",
      RIDRETH3 == 7 ~ "Other/Multiracial"),
      levels = c("Non-Hispanic White", "Non-Hispanic Black", "Hispanic",
                 "Non-Hispanic Asian", "Other/Multiracial")),
    education = factor(case_when(
      DMDEDUC2 %in% c(1, 2) ~ "Less than high school",
      DMDEDUC2 == 3 ~ "High school/GED",
      DMDEDUC2 == 4 ~ "Some college",
      DMDEDUC2 == 5 ~ "College graduate"),
      levels = c("College graduate", "Some college", "High school/GED",
                 "Less than high school")),
    pir_cat = factor(case_when(
      INDFMPIR < 1.3 ~ "<1.3",
      INDFMPIR < 3.5 ~ "1.3-3.49",
      INDFMPIR >= 3.5 ~ ">=3.5"),
      levels = c(">=3.5", "1.3-3.49", "<1.3")),
    smoking = factor(case_when(
      SMQ020 == 2 ~ "Never",
      SMQ040 %in% c(1, 2) ~ "Current",
      SMQ040 == 3 ~ "Former"),
      levels = c("Never", "Former", "Current")),
    log_cotinine = log(LBXCOT + 0.011),
    insured = factor(case_when(HIQ011 == 1 ~ "Yes", HIQ011 == 2 ~ "No"),
                     levels = c("Yes", "No")),
    diabetes = any_true(DIQ010 == 1, LBXGH >= 6.5, LBXGLU >= 126)
  )

# ---- 7. Sequential exclusions and flowchart ---------------------------------
steps <- list(
  "All NHANES 2017-March 2020 participants" = rep(TRUE, nrow(d)),
  "Aged 20 years or older"                  = d$RIDAGEYR >= 20,
  "Examined at the MEC"                     = d$RIDSTATR == 2,
  "Completed liver elastography"            = d$lux_complete %in% 1,
  "Reliable elastography (>=10 valid, IQR/median <30%)" = d$lux_reliable %in% 1,
  "Complete dentition exam"                 = dent_ok,
  "Not pregnant"                            = d$pregnant == 0,
  "No hepatitis B (HBsAg+) or C (RNA+)"     = d$viral_hep == 0,
  "No self-reported viral/autoimmune hepatitis" = d$self_hep == 0,
  "No heavy alcohol use"                    = if (EXCLUDE_MISSING_ALCOHOL)
                                                d$heavy_alcohol %in% 0 else
                                                !(d$heavy_alcohol %in% 1),
  "MASLD status classifiable"               = !is.na(d$masld)
)

keep <- rep(TRUE, nrow(d))
flow <- tibble(step = character(), remaining = integer(), excluded = integer())
for (s in names(steps)) {
  before <- sum(keep)
  keep   <- keep & (steps[[s]] %in% TRUE)
  flow   <- add_row(flow, step = s, remaining = sum(keep), excluded = before - sum(keep))
}
d$analytic <- as.integer(keep)

cat("\n===== PARTICIPANT FLOW =====\n")
print(flow, n = Inf, width = Inf)

# Breakdown of the alcohol step (heavy vs unknown)
cat("\nAlcohol step detail (among those reaching it):\n")
at_alc <- Reduce(`&`, lapply(steps[1:9], function(x) x %in% TRUE))
print(table(heavy_alcohol = d$heavy_alcohol[at_alc], useNA = "ifany"))

# ---- 8. Sanity checks (weighted) --------------------------------------------
des <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMECPRP,
                 nest = TRUE, data = d %>% filter(!is.na(WTMECPRP), WTMECPRP > 0))
sub <- subset(des, analytic == 1)

cat("\n===== WEIGHTED CHECKS IN THE ANALYTIC SAMPLE =====\n")
cat("MASLD prevalence (%):\n");              print(round(100 * svymean(~masld, sub), 1))
cat("Significant fibrosis, LSM >= 8 (%):\n"); print(round(100 * svymean(~fib_sig, sub), 1))
cat("Cirrhosis range, LSM >= 13.6 (%):\n");   print(round(100 * svymean(~fib_cirr, sub), 1))
cat("Mean missing teeth (of 28):\n");         print(round(svymean(~missing_teeth, sub), 2))
cat("Edentulous (%):\n");                     print(round(100 * svymean(~edentulous, sub), 1))
cat("Any untreated caries (%):\n");           print(round(100 * svymean(~any_caries, sub), 1))
cat("Tooth-loss categories (%):\n");          print(round(100 * svymean(~tooth_loss_cat, sub), 1))
cat("Steatosis without any CMRF (n, unweighted):",
    sum(d$cryptogenic_sld[d$analytic == 1], na.rm = TRUE), "\n")

# Raw caries codes in the analytic sample (for the dentist to review)
code_check <- as_tibble(table(code = ctc[d$analytic == 1, ])) %>%
  arrange(desc(n))
cat("\nCoronal caries codes across all 28 teeth, analytic sample:\n")
print(code_check, n = Inf)

# ---- 9. Save ----------------------------------------------------------------
saveRDS(d, here("data/processed/nhanes_merged.rds"))
write_csv(flow, here("output/tables/flowchart_counts.csv"))
write_csv(code_check, here("output/tables/dental_code_check.csv"))

message("\nScript 02 finished. Send back the PARTICIPANT FLOW, the alcohol detail, ",
        "and the WEIGHTED CHECKS for review.")
