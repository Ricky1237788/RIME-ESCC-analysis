source("00_packages_and_setup.R")

rime_hr <- function(d, endpoint, label) {
  coxph(as.formula(paste0("Surv(", endpoint, "_time, ", endpoint, "_event) ~ rime_group")), data = d) %>%
    tidy(exponentiate = TRUE, conf.int = TRUE) %>% mutate(endpoint = endpoint, stratum = label)
}

fig5_forest <- bind_rows(
  map_dfr(c(1, 5, 10), ~ rime_hr(filter(exp_df, pdl1 >= .x), "OS", paste0("PD-L1 >=", .x, "%"))),
  map_dfr(c(1, 5, 10), ~ rime_hr(filter(exp_df, pdl1 >= .x), "PFS", paste0("PD-L1 >=", .x, "%")))
)
write.csv(fig5_forest, "output/Figure5_forest_data.csv", row.names = FALSE)

# Figure 5B defines PD-L1-high as >=10% and PD-L1-low as <10%.
joint_cutoff <- 10
joint_df <- exp_df %>% mutate(
  pdl1_binary = if_else(pdl1 >= joint_cutoff, "PD-L1-high", "PD-L1-low"),
  joint_group = factor(paste(rime_group, pdl1_binary, sep = " / "))
)
fit_joint <- survfit(Surv(OS_time, OS_event) ~ joint_group, data = joint_df)
print(ggsurvplot(fit_joint, data = joint_df, risk.table = TRUE, pval = TRUE, legend.title = NULL))
