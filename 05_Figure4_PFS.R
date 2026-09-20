source("00_packages_and_setup.R")

plot_km <- function(d, label) {
  fit <- survfit(Surv(PFS_time, PFS_event) ~ rime_group, data = d)
  g <- ggsurvplot(fit, data = d, risk.table = TRUE, pval = TRUE, conf.int = FALSE,
                  palette = c("#4DBBD5", "#E64B35"), title = label, xlab = "Time (months)")
  print(g)
}

pdf("output/Figure4_PFS.pdf", width = 11, height = 8)
plot_km(ctl_df, "Control arm")
plot_km(exp_df, "Experimental arm")
plot_km(filter(exp_df, pdl1 < 1), "PD-L1 <1%")
plot_km(filter(exp_df, pdl1 >= 1), "PD-L1 >=1%")
plot_km(filter(exp_df, pdl1 >= 5), "PD-L1 >=5%")
plot_km(filter(exp_df, pdl1 >= 10), "PD-L1 >=10%")
dev.off()

