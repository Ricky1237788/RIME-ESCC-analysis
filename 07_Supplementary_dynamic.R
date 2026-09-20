source("00_packages_and_setup.R")

# Expected file: one record per patient containing baseline and ONGOING raw
# poly(A)-normalized RIME values, plus outcome variables. The manuscript
# defines dynamic change as (ONGOING - Baseline) / Baseline.
dynamic_file <- "data/longitudinal_rime.csv"
if (!file.exists(dynamic_file)) stop("Provide a de-identified longitudinal data file at: ", dynamic_file)
long_df <- read_csv(dynamic_file, show_col_types = FALSE)

dynamic_cutoff <- -0.162
wide_df <- long_df %>%
  mutate(
    rime_change = (rime_post - rime_baseline) / rime_baseline,
    dynamic_group = factor(if_else(rime_change <= dynamic_cutoff,
                                   "RIME-decline", "RIME-increase"),
                           levels = c("RIME-decline", "RIME-increase"))
  )

km_pfs <- survfit(Surv(PFS_time, PFS_event) ~ dynamic_group, data = wide_df)
km_os <- survfit(Surv(OS_time, OS_event) ~ dynamic_group, data = wide_df)

pdf("output/FigureS2_dynamic_RIME.pdf", width = 10, height = 4.5)
print(ggsurvplot(km_pfs, data = wide_df, pval = TRUE, risk.table = FALSE,
                 title = "PFS | ONGOING-Baseline / Baseline"))
print(ggsurvplot(km_os, data = wide_df, pval = TRUE, risk.table = FALSE,
                 title = "OS | ONGOING-Baseline / Baseline"))
dev.off()

write_csv(wide_df, "output/dynamic_RIME_groups.csv")
