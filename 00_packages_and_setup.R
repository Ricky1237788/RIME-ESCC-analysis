# Run this script before all other analyses.

required_packages <- c("tidyverse", "survival", "survminer", "broom", "pROC", "timeROC", "patchwork")
missing_packages <- required_packages[!required_packages %in% rownames(installed.packages())]
if (length(missing_packages) > 0) install.packages(missing_packages)
invisible(lapply(required_packages, library, character.only = TRUE))

data_file <- "data/analysis_data.csv" # Do not upload identifiable patient-level data.
if (!file.exists(data_file)) stop("Place the de-identified analysis CSV at: ", data_file)

df <- readr::read_csv(data_file, show_col_types = FALSE)

# Input `rime` is the final baseline variable: ln(poly_RNA_ave) in the locked source workbook.
# Manuscript-defined exploratory cutoff, displayed to two decimals. The underlying
# cutoff is derived with maximally selected rank statistics in script 01.
rime_cutoff <- 2.23

df <- df %>%
  mutate(
    arm = factor(arm, levels = c("Control", "Experimental")),
    rime_group = factor(if_else(rime <= rime_cutoff, "RIME-low", "RIME-high"),
                        levels = c("RIME-low", "RIME-high")),
    pdl1_category = case_when(
      pdl1 < 1 ~ "<1%",
      pdl1 < 5 ~ "1-<5%",
      pdl1 < 10 ~ "5-<10%",
      TRUE ~ ">=10%"
    ),
    pdl1_category = factor(pdl1_category, levels = c("<1%", "1-<5%", "5-<10%", ">=10%")),
    pdl1_1 = if_else(pdl1 >= 1, "PD-L1-high", "PD-L1-low"),
    pdl1_5 = if_else(pdl1 >= 5, "PD-L1-high", "PD-L1-low"),
    pdl1_10 = if_else(pdl1 >= 10, "PD-L1-high", "PD-L1-low")
  )

exp_df <- filter(df, arm == "Experimental")
ctl_df <- filter(df, arm == "Control")

dir.create("output", showWarnings = FALSE)

save(df, exp_df, ctl_df, rime_cutoff, file = "output/analysis_objects.RData")
