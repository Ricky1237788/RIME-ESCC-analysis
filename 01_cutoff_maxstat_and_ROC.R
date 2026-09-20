source("00_packages_and_setup.R")

# The analysis grouping cutoff is derived using maximally selected rank
# statistics based on PFS time and event status. The reported/displayed value
# is rounded to two decimal places (2.23 in the locked analysis).
cutpoint_fit <- survminer::surv_cutpoint(
  df %>% dplyr::select(PFS_time, PFS_event, rime),
  time = "PFS_time", event = "PFS_event", variables = "rime", minprop = 0.10
)
cutoff_raw <- cutpoint_fit$cutpoint[["cutpoint"]][1]
cutoff_display <- round(cutoff_raw, 2)
write.csv(data.frame(cutoff_raw = cutoff_raw, cutoff_display = cutoff_display,
                     analysis_cutoff = rime_cutoff),
          "output/maxstat_cutoff.csv", row.names = FALSE)

if (!isTRUE(all.equal(cutoff_display, rime_cutoff))) {
  warning("The rounded data-derived cutoff differs from rime_cutoff in script 00.")
}

# Supplementary Figure S1: descriptive time-dependent ROC at the 6-month PFS milestone.
tdroc <- timeROC::timeROC(T = df$PFS_time, delta = df$PFS_event, marker = df$rime,
                           cause = 1, weighting = "marginal", times = 6, iid = TRUE)
t_index <- which(names(tdroc$AUC) == "t=6")
if (length(t_index) != 1) stop("Could not identify the 6-month ROC column.")

# Locate the timeROC coordinate corresponding to the displayed analysis cutoff.
markers <- sort(df$rime, decreasing = TRUE)
cutoff_index <- which.min(abs(markers - rime_cutoff))
sensitivity <- tdroc$TP[cutoff_index, t_index]
specificity <- 1 - tdroc$FP[cutoff_index, t_index]
roc_summary <- data.frame(AUC = tdroc$AUC[t_index], cutoff = rime_cutoff,
                          sensitivity = sensitivity, specificity = specificity)
write.csv(roc_summary, "output/FigureS1_ROC_summary.csv", row.names = FALSE)

pdf("output/FigureS1_time_dependent_ROC.pdf", width = 7.5, height = 7.5)
par(cex.axis = 0.85, cex.lab = 0.9)
plot(tdroc, time = 6, col = "#1F77B4", title = FALSE)
abline(a = 0, b = 1, lty = 2, col = "#FF7F0E")
points(1 - specificity, sensitivity, pch = 16, col = "#1F77B4", cex = 1.25)
title(main = "Time-dependent ROC anchored at 6-month PFS (PFS)", cex.main = 0.85)
legend("bottomright",
       legend = c(sprintf("AUC = %.3f", tdroc$AUC[t_index]),
                  sprintf("Cutoff = %.2f", rime_cutoff),
                  sprintf("Sens = %.2f, Spec = %.2f", sensitivity, specificity)),
       bty = "n")
dev.off()

# Bootstrap stability of the maximally selected rank-statistic cutoff.
set.seed(2026)
B <- 1000
boot_cutoffs <- replicate(B, {
  d <- df %>% dplyr::select(PFS_time, PFS_event, rime) %>% tidyr::drop_na()
  d <- d[sample.int(nrow(d), nrow(d), replace = TRUE), ]
  fit <- try(survminer::surv_cutpoint(d, time = "PFS_time", event = "PFS_event",
                                      variables = "rime", minprop = 0.10), silent = TRUE)
  if (inherits(fit, "try-error")) return(NA_real_)
  fit$cutpoint[["cutpoint"]][1]
})
boot_cutoffs <- boot_cutoffs[!is.na(boot_cutoffs)]
write.csv(data.frame(cutoff = boot_cutoffs), "output/bootstrap_cutoffs.csv", row.names = FALSE)

ggplot(data.frame(cutoff = boot_cutoffs), aes(cutoff)) +
  geom_histogram(bins = 30, fill = "#4DBBD5", colour = "white") +
  geom_vline(xintercept = cutoff_raw, linetype = 2, colour = "#E64B35") +
  theme_classic() +
  labs(x = "Maximally selected rank-statistic cutoff for RIME", y = "Bootstrap frequency")
ggsave("output/FigureS1_bootstrap_cutoff_stability.pdf", width = 6, height = 4)
