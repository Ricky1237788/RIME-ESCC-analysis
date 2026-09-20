source("00_packages_and_setup.R")

# Treatment-by-RIME interaction models. Control and RIME-low are reference levels.
fit_interaction <- function(d, endpoint) {
  if (endpoint %in% c("OS", "PFS")) {
    coxph(as.formula(paste0("Surv(", endpoint, "_time, ", endpoint, "_event) ~ arm * rime_group")), data = d)
  } else {
    glm(ORR ~ arm * rime_group, family = binomial, data = d)
  }
}

interaction_results <- bind_rows(
  tidy(fit_interaction(df, "OS"), exponentiate = TRUE, conf.int = TRUE) %>% mutate(endpoint = "OS"),
  tidy(fit_interaction(df, "PFS"), exponentiate = TRUE, conf.int = TRUE) %>% mutate(endpoint = "PFS"),
  tidy(fit_interaction(df, "ORR"), exponentiate = TRUE, conf.int = TRUE) %>% mutate(endpoint = "ORR")
)

# Repeat in protocol-defined PD-L1 strata if required for Figure S5.
interaction_by_stratum <- function(d, label) {
  bind_rows(
    tidy(fit_interaction(d, "OS"), exponentiate = TRUE, conf.int = TRUE) %>% mutate(endpoint = "OS"),
    tidy(fit_interaction(d, "PFS"), exponentiate = TRUE, conf.int = TRUE) %>% mutate(endpoint = "PFS"),
    tidy(fit_interaction(d, "ORR"), exponentiate = TRUE, conf.int = TRUE) %>% mutate(endpoint = "ORR")
  ) %>% mutate(PDL1_stratum = label)
}

figS5_data <- bind_rows(
  interaction_by_stratum(filter(df, pdl1 < 1), "PD-L1 <1%"),
  interaction_by_stratum(filter(df, pdl1 >= 1), "PD-L1 >=1%"),
  interaction_by_stratum(filter(df, pdl1 >= 5), "PD-L1 >=5%"),
  interaction_by_stratum(filter(df, pdl1 >= 10), "PD-L1 >=10%")
)
write.csv(figS5_data, "output/FigureS5_interaction_data.csv", row.names = FALSE)

