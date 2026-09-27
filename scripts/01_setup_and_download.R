# =============================================================================
# Dental status and MASLD / liver fibrosis — NHANES 2017–March 2020
# Script 01: project setup, package installation, and raw data download
# -----------------------------------------------------------------------------
# Run this ONCE from inside the RStudio Project folder (dental-masld.Rproj).
# It creates the folder structure, installs packages, downloads every NHANES
# file we need from the CDC, and saves each one as .rds for fast loading.
# =============================================================================

# ---- 1. Packages ------------------------------------------------------------
pkgs <- c(
  "tidyverse",  # data wrangling and plotting
  "haven",      # reads NHANES .xpt (SAS transport) files
  "survey",     # survey-weighted analysis (core of the whole project)
  "srvyr",      # tidy wrapper around survey
  "gtsummary",  # publication-ready Table 1 and regression tables
  "rms",        # restricted cubic splines for dose-response (Aim 1)
  "pROC",       # AUC / discrimination (Aim 3)
  "here",       # file paths relative to the project folder
  "janitor"     # quick tabulations and clean names
)
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

# ---- 2. Folder structure ----------------------------------------------------
dirs <- c("data/raw", "data/processed", "scripts",
          "output/tables", "output/figures", "docs")
for (d in dirs) dir.create(here(d), recursive = TRUE, showWarnings = FALSE)

# ---- 3. Files to download (2017–March 2020 pre-pandemic, prefix "P_") -------
files <- c(
  # Core
  P_DEMO   = "Demographics, survey weights, strata, PSU, pregnancy",
  P_OHXDEN = "Oral health exam: dentition (tooth count, caries) — EXPOSURE",
  P_LUX    = "Liver elastography: CAP and liver stiffness — OUTCOME",
  # MASLD cardiometabolic criteria
  P_BMX    = "Body measures: BMI, waist circumference",
  P_BPXO   = "Blood pressure (oscillometric)",
  P_BPQ    = "Blood pressure / cholesterol medication questionnaire",
  P_GLU    = "Fasting plasma glucose (fasting subsample)",
  P_GHB    = "Glycohemoglobin (HbA1c)",
  P_DIQ    = "Diabetes questionnaire (diagnosis, insulin, pills)",
  P_TRIGLY = "Triglycerides (fasting subsample)",
  P_HDL    = "HDL cholesterol",
  # Exclusions
  P_HEPBD  = "Hepatitis B core antibody and surface antigen",
  P_HEPC   = "Hepatitis C antibody and RNA",
  P_ALQ    = "Alcohol use questionnaire",
  P_MCQ    = "Medical conditions (liver conditions, autoimmune hepatitis)",
  # Covariates
  P_SMQ    = "Smoking questionnaire",
  P_COT    = "Serum cotinine (objective smoking exposure)",
  P_INQ    = "Income",
  P_HIQ    = "Health insurance",
  P_PAQ    = "Physical activity",
  P_HSCRP  = "High-sensitivity C-reactive protein",
  P_DR1TOT = "Dietary recall day 1: total nutrients (sugar, energy)"
)

# ---- 4. Download ------------------------------------------------------------
base_url <- "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/"
options(timeout = 600)  # some files are large; allow 10 minutes each

download_log <- tibble(file = names(files), status = NA_character_)

for (f in names(files)) {
  dest <- here("data/raw", paste0(f, ".xpt"))
  if (file.exists(dest)) {
    download_log$status[download_log$file == f] <- "already downloaded"
    next
  }
  url <- paste0(base_url, f, ".xpt")
  result <- tryCatch({
    download.file(url, dest, mode = "wb", quiet = TRUE)
    "downloaded"
  }, error = function(e) paste("FAILED:", conditionMessage(e)))
  download_log$status[download_log$file == f] <- result
  message(f, " -> ", result)
}

print(download_log, n = Inf)

# ---- 5. Read each file, check it, and save as .rds --------------------------
check_log <- tibble(file = character(), rows = integer(), cols = integer(),
                    seqn_unique = logical())

for (f in names(files)) {
  src <- here("data/raw", paste0(f, ".xpt"))
  if (!file.exists(src)) next
  dat <- read_xpt(src)
  saveRDS(dat, here("data/raw", paste0(f, ".rds")))
  check_log <- add_row(check_log,
                       file = f, rows = nrow(dat), cols = ncol(dat),
                       seqn_unique = !anyDuplicated(dat$SEQN))
}

print(check_log, n = Inf)

# ---- 6. Quick look at the two key files -------------------------------------
lux    <- readRDS(here("data/raw/P_LUX.rds"))
ohxden <- readRDS(here("data/raw/P_OHXDEN.rds"))
demo   <- readRDS(here("data/raw/P_DEMO.rds"))

cat("\nP_LUX variables:\n");    print(names(lux))
cat("\nP_OHXDEN variables (first 40):\n"); print(head(names(ohxden), 40))

cat("\nParticipants aged 20+ in DEMO:", sum(demo$RIDAGEYR >= 20, na.rm = TRUE), "\n")
cat("Of those, with an elastography record:",
    sum(demo$SEQN[demo$RIDAGEYR >= 20] %in% lux$SEQN), "\n")
cat("Of those, with a dentition exam record:",
    sum(demo$SEQN[demo$RIDAGEYR >= 20] %in% ohxden$SEQN), "\n")

# ---- 7. Record the R environment (needed for reproducibility) ---------------
writeLines(capture.output(sessionInfo()), here("docs/sessionInfo_script01.txt"))
write_csv(check_log, here("docs/download_check_script01.csv"))

message("\nScript 01 finished. Send the printed tables from steps 4, 5 and 6 back for review.")
