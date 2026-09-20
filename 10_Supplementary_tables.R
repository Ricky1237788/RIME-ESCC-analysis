source("00_packages_and_setup.R")

# Table S1: cutoff sensitivity. Replace the endpoint/model if the final SAP specifies otherwise.
cutoff_models <- list(
  Analysis_derived = if_else(df$rime <= rime_cutoff, "RIME-low", "RIME-high"),
  Median = if_else(df$rime <= median(df$rime, na.rm = TRUE), "RIME-low", "RIME-high"),
  Tertile = cut(df$rime, breaks = quantile(df$rime, c(0, 1/3, 2/3, 1), na.rm = TRUE), include.lowest = TRUE)
)

# Table S2: BH-FDR correction across the 12 RIME-high versus RIME-low comparisons
# (ORR, PFS, and OS in PD-L1 <1%, >=1%, >=5%, and >=10% strata).
nominal_rime_test <- function(endpoint, d, label) {
  if (endpoint == "ORR") {
    p_value <- fisher.test(table(d$ORR, d$rime_group))$p.value
  } else {
    fit <- coxph(as.formula(paste0("Surv(", endpoint, "_time, ", endpoint, "_event) ~ rime_group")), data = d)
    p_value <- tidy(fit)$p.value[1]
  }
  data.frame(endpoint = endpoint, PDL1_stratum = label, p_value = p_value)
}

pdl1_strata <- list(
  "PD-L1 <1%" = filter(exp_df, pdl1 < 1),
  "PD-L1 >=1%" = filter(exp_df, pdl1 >= 1),
  "PD-L1 >=5%" = filter(exp_df, pdl1 >= 5),
  "PD-L1 >=10%" = filter(exp_df, pdl1 >= 10)
)
tableS2 <- imap_dfr(pdl1_strata, ~ map_dfr(c("ORR", "PFS", "OS"), nominal_rime_test, d = .x, label = .y)) %>%
  mutate(q_value_BH = p.adjust(p_value, method = "BH"))
write.csv(tableS2, "output/TableS2_BH_FDR.csv", row.names = FALSE)

# Table S3: multivariable models. Confirm covariate names and coding before running.
covariates <- "rime_group + age + sex + ECOG + stage + LN_status + metastatic_burden + prior_radiotherapy + factor(pdl1)"
fit_multi_os <- coxph(as.formula(paste("Surv(OS_time, OS_event) ~", covariates)), data = exp_df)
fit_multi_pfs <- coxph(as.formula(paste("Surv(PFS_time, PFS_event) ~", covariates)), data = exp_df)
fit_multi_orr <- glm(as.formula(paste("ORR ~", covariates)), family = binomial, data = exp_df)
tableS3 <- bind_rows(
  tidy(fit_multi_os, exponentiate = TRUE, conf.int = TRUE) %>% mutate(endpoint = "OS"),
  tidy(fit_multi_pfs, exponentiate = TRUE, conf.int = TRUE) %>% mutate(endpoint = "PFS"),
  tidy(fit_multi_orr, exponentiate = TRUE, conf.int = TRUE) %>% mutate(endpoint = "ORR")
)
write.csv(tableS3, "output/TableS3_multivariable_models.csv", row.names = FALSE)

# Table S4: safety outcomes. Expected variables: AE_grade, SAE, irAE, AE_discontinuation.
tableS4 <- exp_df %>% group_by(rime_group) %>% summarise(
  n = n(), grade3plus = sum(AE_grade >= 3, na.rm = TRUE),
  SAE = sum(SAE == 1, na.rm = TRUE), irAE = sum(irAE == 1, na.rm = TRUE),
  AE_discontinuation = sum(AE_discontinuation == 1, na.rm = TRUE), .groups = "drop"
)
write.csv(tableS4, "output/TableS4_safety.csv", row.names = FALSE)
