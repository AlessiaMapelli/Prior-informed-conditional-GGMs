#################################################
## A1. Global permutation result: T2D edge counts, observed vs null,
## 10-year endpoint (and the 5-year run as contrast).
## Run from the private project folder: Rscript --vanilla ../T2D_application_code/downstream/A1_global_null.R
#################################################
source(file.path(Sys.getenv("T2D_CODE_DIR", "../T2D_application_code"), "downstream", "00_common.R"))

runs <- c(`10-year (1,112 cases)` = file.path(PERM_DIR, "edge_counts_by_permutation.csv"),
          `5-year (about 270 cases)` = "results_bic_relaxed/permutation_null/edge_counts_by_permutation.csv")
d <- do.call(rbind, lapply(names(runs), function(r) {
  x <- read.csv(runs[[r]]); x <- x[x$component == "Diab_Incident", ]
  rbind(data.frame(endpoint = r, rule = "AND", perm = x$perm, count = x$AND),
        data.frame(endpoint = r, rule = "OR",  perm = x$perm, count = x$OR))
}))
d$endpoint <- factor(d$endpoint, levels = names(runs))

tab <- do.call(rbind, lapply(split(d, list(d$endpoint, d$rule)), function(x) {
  obs <- x$count[x$perm == 0]; nul <- x$count[x$perm > 0]
  data.frame(endpoint = x$endpoint[1], rule = x$rule[1], B = length(nul), observed = obs,
             null_mean = mean(nul), null_median = median(nul),
             null_q95 = unname(quantile(nul, 0.95)), null_max = max(nul),
             p = (1 + sum(nul >= obs)) / (length(nul) + 1),
             expected_false_share = min(1, mean(nul) / max(obs, 1)))
}))
rownames(tab) <- NULL
print(tab, digits = 3)
write.csv(tab, file.path(OUT_DIR, "A1_global_test.csv"), row.names = FALSE)

obs <- d[d$perm == 0, ]
lab <- merge(obs, tab[, c("endpoint", "rule", "p", "null_mean")])
lab$text <- sprintf("observed %d\nnull mean %.1f\np = %.3f", lab$count, lab$null_mean, lab$p)
g <- ggplot(d[d$perm > 0, ], aes(count)) +
  geom_histogram(binwidth = function(x) max(1, ceiling(diff(range(x)) / 30)),
                 fill = "grey70", colour = "white", boundary = 0) +
  geom_vline(data = obs, aes(xintercept = count), colour = "#d95f02", linewidth = 0.8) +
  geom_text(data = lab, aes(x = Inf, y = Inf, label = text), hjust = 1.05, vjust = 1.2, size = 2.6) +
  facet_grid(endpoint ~ rule, scales = "free") +
  labs(x = "Number of T2D-specific edges", y = "Permutations (B = 100)",
       title = "T2D component: observed edge count (orange) vs permutation null")
save_fig(g, "FigS_A1_permutation_null", 7, 5)
