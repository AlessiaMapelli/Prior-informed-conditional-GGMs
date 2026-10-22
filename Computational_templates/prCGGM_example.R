# Example of the prCGGM R package on the template example data
# (1000 subjects, 10 proteins, 1 binary covariate, perfect prior).
# It reproduces the README example of ggReg_main_functions.R with the
# package interface: prcggm() replaces GGReg_full_estimation() and
# predict() replaces predict_personalized_network().
#
# Install the package once with: install.packages("prCGGM")
# Run this script from the Computational_templates folder.

library(prCGGM)

load("Example_data.RData")  # expression_data, covariate_data, ppi_weight_matrix

#################################################
## Model estimation
#################################################

# set.seed() makes the cross-validation folds reproducible; with seed 2024
# the results match GGReg_full_estimation() with its default settings.
set.seed(2024)
fit <- prcggm(
  x          = expression_data,    # subjects x features
  covariates = covariate_data,     # subjects x covariates
  prior      = ppi_weight_matrix,  # features x features, values in [0, 1]
  verbose    = TRUE
)

fit
summary(fit)

# Main outputs (names in GGReg_full_estimation() in brackets)
fit$mean_effects           # covariate effects on the mean (Cov_effect)
names(fit$networks)        # adjacency matrices (Dic_adj_matrics)
fit$networks$Baseline
fit$networks$x1_1
fit$tuning                 # selected lambda, alpha, theta and BIC per node (optimal_params)

#################################################
## Personalized networks
#################################################

new_subjects <- data.frame(x1 = factor(c("0", "1")))
personal_net <- predict(fit, new_subjects)
names(personal_net)

#################################################
## Visualization and comparison
#################################################

plot(fit)                                   # baseline network
plot(fit, component = "x1_1")               # change associated with x1 = 1
plot_network(personal_net$Subject_1, main = "Subject 1")
plot_network_difference(personal_net$Subject_2, personal_net$Subject_1)
compare_networks(personal_net$Subject_1, personal_net$Subject_2,
                 method = "frobenius")      # also "difference", "correlation"

#################################################
## Larger problems
#################################################

# Parallel estimation over nodes on the local machine:
# fit <- prcggm(expression_data, covariate_data, ppi_weight_matrix, ncores = 4)
#
# SLURM array jobs (one job per node), see ?prcggm_slurm_prepare:
# job <- prcggm_slurm_prepare(expression_data, covariate_data,
#                             ppi_weight_matrix, dir = "results/prcggm")
# ... submit job$command from a shell, then once all jobs have completed:
# fit <- prcggm_slurm_collect("results/prcggm")
