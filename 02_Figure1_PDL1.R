source("00_packages_and_setup.R")

# PD-L1 populations are nested and are analysed separately, not as mutually exclusive groups.
threshold_cox <- function(endpoint, threshold) {
  d <- exp_df %>% filter(pdl1 < 1 | pdl1 >= threshold) %>%
    mutate(pdl1_group = factor(if_else(pdl1 < 1, "PD-L1 <1%", paste0("PD-L1 >=", threshold, "%"))))
  coxph(as.formula(paste0("Surv(", endpoint, "_time, ", endpoint, "_event) ~ pdl1_group")), data = d) %>%
    tidy(exponentiate = TRUE, conf.int = TRUE) %>% mutate(endpoint = endpoint, threshold = threshold)
}

fig1_forest <- bind_rows(
  map_dfr(c(1, 5, 10), ~ threshold_cox("OS", .x)),
  map_dfr(c(1, 5, 10), ~ threshold_cox("PFS", .x))
)
write.csv(fig1_forest, "output/Figure1_forest_data.csv", row.names = FALSE)

fig1_orr <- map_dfr(c(1, 5, 10), function(x) {
  exp_df %>% filter(pdl1 < 1 | pdl1 >= x) %>%
    mutate(group = if_else(pdl1 < 1, "PD-L1 <1%", paste0("PD-L1 >=", x, "%"))) %>%
    group_by(group) %>% summarise(responders = sum(ORR == 1), total = n(), ORR = responders / total, .groups = "drop") %>%
    mutate(threshold = x)
})
write.csv(fig1_orr, "output/Figure1_ORR_data.csv", row.names = FALSE)

