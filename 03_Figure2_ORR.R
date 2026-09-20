source("00_packages_and_setup.R")

make_orr <- function(d, label) {
  d %>% group_by(rime_group) %>%
    summarise(responders = sum(ORR == 1, na.rm = TRUE), total = n(), ORR = responders / total, .groups = "drop") %>%
    mutate(stratum = label)
}

figure2_data <- bind_rows(
  make_orr(ctl_df, "Control overall"), make_orr(exp_df, "Experimental overall"),
  make_orr(filter(exp_df, pdl1 < 1), "PD-L1 <1%"),
  make_orr(filter(exp_df, pdl1 >= 1), "PD-L1 >=1%"),
  make_orr(filter(exp_df, pdl1 >= 5), "PD-L1 >=5%"),
  make_orr(filter(exp_df, pdl1 >= 10), "PD-L1 >=10%")
)
write.csv(figure2_data, "output/Figure2_ORR_data.csv", row.names = FALSE)

ggplot(figure2_data, aes(rime_group, ORR, fill = rime_group)) +
  geom_col(width = .65) + geom_text(aes(label = paste0(responders, "/", total)), vjust = -0.3) +
  facet_wrap(~stratum) + scale_y_continuous(labels = scales::percent, limits = c(0, 1)) + theme_classic()
ggsave("output/Figure2_ORR.pdf", width = 10, height = 6)
