source("00_packages_and_setup.R")

# Compare residual categorical PD-L1 information alone with PD-L1 + RIME in eligible positive strata.
fit_nested_survival <- function(d, endpoint) {
  fit_base <- coxph(as.formula(paste0("Surv(", endpoint, "_time, ", endpoint, "_event) ~ pdl1_category")), data = d)
  fit_full <- update(fit_base, . ~ . + rime_group)
  anova(fit_base, fit_full, test = "LRT")
}

fit_nested_orr <- function(d) {
  fit_base <- glm(ORR ~ pdl1_category, family = binomial, data = d)
  fit_full <- update(fit_base, . ~ . + rime_group)
  anova(fit_base, fit_full, test = "LRT")
}

d_ge1 <- filter(exp_df, pdl1 >= 1)
d_ge5 <- filter(exp_df, pdl1 >= 5)

results_s4 <- list(
  OS_ge1 = fit_nested_survival(d_ge1, "OS"), PFS_ge1 = fit_nested_survival(d_ge1, "PFS"), ORR_ge1 = fit_nested_orr(d_ge1),
  OS_ge5 = fit_nested_survival(d_ge5, "OS"), PFS_ge5 = fit_nested_survival(d_ge5, "PFS"), ORR_ge5 = fit_nested_orr(d_ge5)
)
saveRDS(results_s4, "output/Supplementary_FigureS4_nested_models.rds")
