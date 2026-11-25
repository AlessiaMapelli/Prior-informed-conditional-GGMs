#################################################
## Compact ablation figure for the main text (Fig 2 panel): selective prior
## (prCGGM, baseline only) vs prior on all terms, AND, p = 100, q = 3.
## Reads results_ablation/summary/ablation_by_replication.csv (analyse_ablation.R).
## Usage (from Simulation_analysis_v2): Rscript --vanilla fig_ablation_main.R [results_ablation]
## Output: <ablation>/summary/Fig_ablation_main.{png,pdf}
#################################################
suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(trailingOnly = TRUE)
out <- file.path(if (length(args)) args[1] else "results_ablation", "summary")
d <- read.csv(file.path(out, "ablation_by_replication.csv"))
d <- d[d$symm == "AND", ]
arm_lab <- c(none = "No prior", noisy_baseline = "Noisy prior, baseline only (prCGGM)",
             noisy_all = "Noisy prior, all terms", perfect_baseline = "Perfect prior, baseline only (prCGGM)",
             perfect_all = "Perfect prior, all terms")
panels <- data.frame(component = c("Baseline", "x3_1", "x1", "x3_1"), metric = c("F1", "F1", "F1", "FPR_x1000"),
                     title = c("a. F1, baseline network", "b. F1, binary-effect network",
                               "c. F1, continuous-effect network", "d. FPR x 1000, binary-effect network"))
long <- do.call(rbind, lapply(seq_len(nrow(panels)), function(k) {
  x <- d[d$component == panels$component[k], ]
  data.frame(n = x$n, arm = x$arm, value = x[[panels$metric[k]]], panel = panels$title[k])
}))
long$panel <- factor(long$panel, panels$title)
long$Arm <- factor(arm_lab[long$arm], arm_lab)
long$N <- factor(paste0(long$n / 1000, "k"), paste0(sort(unique(long$n)) / 1000, "k"))
pal <- c(`No prior` = "#7f7f7f", `Noisy prior, baseline only (prCGGM)` = "#1f78b4", `Noisy prior, all terms` = "#a6cee3",
         `Perfect prior, baseline only (prCGGM)` = "#33a02c", `Perfect prior, all terms` = "#b2df8a")
g <- ggplot(long, aes(N, value, fill = Arm)) +
  geom_boxplot(outlier.size = 0.2, linewidth = 0.25, width = 0.8, position = position_dodge(0.85)) +
  facet_wrap(~ panel, scales = "free_y", nrow = 1) + scale_fill_manual(values = pal) +
  labs(x = "Sample size n", y = NULL, fill = NULL) +
  theme_bw(base_size = 8) + theme(legend.position = "bottom", panel.grid.minor = element_blank(),
                                  strip.text = element_text(hjust = 0)) +
  guides(fill = guide_legend(nrow = 2))
ggsave(file.path(out, "Fig_ablation_main.png"), g, width = 9, height = 3.2, dpi = 300, bg = "white")
ggsave(file.path(out, "Fig_ablation_main.pdf"), g, width = 9, height = 3.2)
cat("Written to", out, "\n")
