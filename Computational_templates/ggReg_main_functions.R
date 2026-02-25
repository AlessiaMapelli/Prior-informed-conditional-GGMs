library(sparsegl)
library(glmnet) 
library(parallel)
library(doParallel)
library(foreach)

#################################################
## Full estimation function: mean + precision estimation
#################################################

#' GGReg_full_estimation
#'
#' Performs full estimation for the GGReg model, including parameter estimation and model fitting.
#'
#' @param x Expression data - data frame; nrow: subjects; ncol: features
#' @param known_ppi Previously available knowledge on interaction - matrix of values [0,1] that describe the confeidence of interaction existing; ncol: features, nrow: features. If NULL, no prior knowledge is used (default NULL)
#' @param covariates Covariates to regress on - data frame; nrow: subjects; ncol: covariates). If NULL, no covariates are included in the model (default NULL)
#' @param scr Logical indicating whether to perform screening based on correlation to speed up estimation - Logical; default TRUE
#' @param gamma Pearson correlation threshold for screening (if scr = TRUE) - Numeric; If NULL 20% quantile of absolute correlations (default NULL)
#' @param mean_estimation Logical indicating whether to perform mean estimation, can be set to false if the input data have a null mean  - Logical; (default TRUE)
#' @param lambda_mean Penalization term in Lasso for mean estimation (if mean_estimation = TRUE) - Numeric; If NULL, optimal lambda is selected via cross-validation (default NULL)
#' @param lambda_mean_type Type of lambda selection for mean estimation when using cross-validation for optimizing (if mean_estimation = TRUE) - Character; "1se" or "min" (default "1se")
#' @param lambda_prec Penalization term in Lasso for precision estimation - Numeric; If NULL, optimal lambda is selected via cross-validation (default NULL)
#' @param lambda_prec_type Type of lambda selection for precision estimation when using cross-validation for optimizing - Character; "1se" or "min" (default "1se")
#' @param tune_hyperparams Logical indicating whether to perform hyperparameter tuning (weight and aloha) for precision estimation - Logical; default TRUE
#' @param asparse_grid Grid of alpha values to consider for hyperparameter tuning (if tune_hyperparams = TRUE) - Numeric vector; default c(0.5, 0.75, 0.9, 0.95)
#' @param weight_grid Grid of weight values to consider for hyperparameter tuning (if tune_hyperparams = TRUE) - Numeric vector; default c(0.8, 1.0, 1.1, 1.3, 1.5)
#' @param random_hyper_search Logical indicating whether to perform random hyperparameter search instead of grid search for tuning (if tune_hyperparams = TRUE) - Logical; default FALSE
#' @param p.rand.hyper Proportion of random hyperparameter combinations to try if random_hyper_search = TRUE (if tune_hyperparams = TRUE) - Numeric; default NULL setting which means 50% of the combinations will be tried
#' @param K Number of folds for cross-validation when computing BIC for hyperparameter tuning (if tune_hyperparams = TRUE) - Numeric; default 5
#' @param use_slurm Logical indicating whether to use HPC (specifically SLURM) for parallelization of precision estimation - Logical; If FALSE computation is going to be done in parallel via foreach is the nodes number is grater then 150 otherwise sequentially (default TRUE)
#' @param slurm_script_path Path to the SLURM script to be used for parallelization (if use_slurm = TRUE) - Character; default "./slurm_ggReg_node.sbatch"
#' @param output_path Directory where SLURM job outputs will be saved (if use_slurm = TRUE) - Character; default "./results/"
#' @param name_output Base name for SLURM job output files (if use_slurm = TRUE) - Character; default "ggReg_result"
#' @param symm_method Method for symmetrization of the estimated precision matrix - Character; "OR" or "AND" (default "OR")
#' @param verbose Logical indicating whether to print detailed progress messages during estimation - Logical; default FALSE
#'
#' @return A list containing:
#'   \item{results}{A list with elements:}
#'     \itemize{
#'       \item{Cov_effect}{Estimated covariance effects from the mean regression.}
#'       \item{Dic_adj_matrics}{Dictionary of covariate-specific symmetrized weighted adjacency matrices estimation for each node reporting partial correlation between nodes.}
#'     }
#'   \item{additional_info}{A list with elements:}
#'     \itemize{
#'       \item{z}{Estimated latent variables from the mean regression.}
#'       \item{Prec_reg_matrix}{Aggregated coefficient regression matrix combining results from all nodes.}
#'       \item{No_sim_Delta_hat}{Aggregated non-symmetrized predicsion matrix combining results from all nodes.}
#'       \item{Sigma_hat}{Aggregated vector of estimated variances for each node.}
#'       \item{Dic_Delta_hat}{Dictionary of symmetrized precision matrices estimation for each node.}
#'       \item{Dic_adj_matrics}{Dictionary of symmetrized weighted adjacency matrices estimation for each node reporting partial correlation between nodes.}
#'       \item{optimal_params}{Data frame of optimal hyperparameters (asparse, weight) and BIC score for each node.}
#'       \item{computational_time_per_node}{Vector of computational time taken for estimation of each node.}
#'       \item{total_computational_time}{Total computational time taken for the entire estimation process.}
#'       \item{scaling_params}{List of scaling parameters (means and SDs) used for numerical covariates.}
#'       \item{dummy_params}{List of dummy coding parameters used for categorical covariates.}
#'     }
#' 
#'
#' @description
#' This function estimates the parameters of the GGReg model using the provided observed data and external covariates. 
#' It first performs mean estimation, and then estimates the precision matrix using either SLURM-based parallelization or a fallback method depending on the number of features. The function also includes options for hyperparameter tuning and screening to improve estimation efficiency.
#' It returns the estimated coefficients, fitted values, residuals, and convergence status.
#' 
#' @examples
#' set.seed(123)
#' n <- 1000  # number of samples
#' p <- 3   # number of features (nodes)
#' q <- 2    # number of covariates
#' # Generate synthetic expression data
#' X <- matrix(rnorm(n * p), n, p)
#' colnames(X) <- paste0("Gene_", 1:p)
#' # Generate synthetic covariates
#' covariates <- data.frame(
#'   age = rnorm(n, mean = 50, sd = 10),
#'   sex = rbinom(n, 1, 0.5)
#' )
#' covariates$sex <- as.factor(ifelse(covariates$sex,"Male", "Female"))
#' 
#' # Optional: Prior knowledge network (PPI)
#' known_ppi <- matrix(0, p, p)
#' known_ppi[1:2, 1:2] <- 0.3
#' diag(known_ppi) <- 0
#' 
#' source("ggReg_main_functions.R")
#' 
#' results <- GGReg_full_estimation(
#'   x = X,
#'   known_ppi = known_ppi,
#'   covariates = covariates)
#' 

GGReg_full_estimation <- function(
  x,
  known_ppi = NULL,
  covariates = NULL,
  scr = TRUE,
  gamma = NULL,
  mean_estimation = TRUE,
  lambda_mean = NULL,
  lambda_mean_type = "1se",
  lambda_prec = NULL,
  lambda_prec_type = "1se",
  tune_hyperparams = TRUE,
  asparse_grid = c(0.5, 0.75, 0.9, 0.95),
  weight_grid = c(0.8, 1.0, 1.1, 1.3, 1.5),
  random_hyper_search = FALSE,
  p.rand.hyper = NULL,
  K = 5,
  use_slurm = FALSE,
  slurm_script_path = "./slurm_ggReg_node.sbatch",
  output_path = "./results/",
  name_output = "ggReg_result",
  symm_method ="OR",
  verbose = FALSE) 
{
  if (verbose) {
    cat("=== Starting GGReg Full Estimation ===\n")
  }
  time.start <- Sys.time()
  if(mean_estimation){
    if (verbose) {cat("Estimating the mean\n")}
    # Step 1: Estimate mean
    res_mean_reg <- GGReg_mean_estimation(
      x = x,
      covariates = covariates,
      lambda_mean = lambda_mean,
      lambda_mean_type = lambda_mean_type,
      verbose = verbose)
    Z <- res_mean_reg$z
    if (verbose) {cat("Done estimating the mean\n")}
  }else{Z <- x}
  # Step 2: Estimate precision matrix with node-wise SLURM parallelization
  if (use_slurm) {
    if (verbose) {
      cat("Estimating the precision matrix with SLURM parallelization over nodes\n")
    }
    # Use SLURM array jobs - one job per node i
    res_cov_reg <- GGReg_cov_estimation_SLURM(
      Z0 = Z,
      known_ppi = known_ppi,
      covariates = covariates,
      scr = scr,
      gamma = gamma,
      lambda_prec = lambda_prec,
      lambda_prec_type = lambda_prec_type,
      tune_hyperparams = tune_hyperparams,
      asparse_grid = asparse_grid,
      weight_grid = weight_grid,
      random_hyper_search=random_hyper_search,
      p.rand.hyper = p.rand.hyper,
      K = K,
      slurm_script_path = slurm_script_path,
      output_path = output_path,
      name_output = name_output,
      symm_method = symm_method,
      verbose = verbose)
  } else { 
    if(ncol(Z)>150){
    if (verbose) {
      cat("Estimating the precision matrix with R parallelization via foreach over nodes\n")
    }
    # Fallback to the R parallelization via foreach
    res_cov_reg <- GGReg_cov_estimation_parallel(
      Z0 = Z,
      known_ppi = known_ppi,
      covariates = covariates,
      scr = scr,
      gamma = gamma,
      lambda_prec = lambda_prec,
      lambda_prec_type = lambda_prec_type,
      tune_hyperparams = tune_hyperparams,
      asparse_grid = asparse_grid,
      weight_grid = weight_grid,
      random_hyper_search=random_hyper_search,
      p.rand.hyper = p.rand.hyper,
      K = K,
      output_path = output_path,
      name_output = name_output,
      symm_method = symm_method,
      verbose = verbose) 
    } else {
      if (verbose) {
        cat("Estimating the precision matrix sequentially over nodes\n")
      }
      res_cov_reg <- GGReg_cov_estimation_sequential(
      Z0 = Z,
      known_ppi = known_ppi,
      covariates = covariates,
      scr = scr,
      gamma = gamma,
      lambda_prec = lambda_prec,
      lambda_prec_type = lambda_prec_type,
      tune_hyperparams = tune_hyperparams,
      asparse_grid = asparse_grid,
      weight_grid = weight_grid,
      random_hyper_search=random_hyper_search,
      p.rand.hyper = p.rand.hyper,
      K = K,
      output_path = output_path,
      name_output = name_output,
      symm_method = symm_method,
      verbose = verbose)
    }
  }
  if (verbose) {
    cat("Done estimating the precision matrix\n")
    cat("=== GGReg Full Estimation Complete ===\n")
  }
  if(mean_estimation){
    results <- list(
    Cov_effect = res_mean_reg$Cov_effect,
    Dic_adj_matrics = res_cov_reg$Dic_adj_matrics)
  }else{
    results <- list(
    Cov_effect = NULL,
    Dic_adj_matrics = res_cov_reg$Dic_adj_matrics)
  }
  additional_info <- list(
    z = Z,
    Prec_reg_matrix = res_cov_reg$Prec_reg_matrix,
    No_sim_Delta_hat = res_cov_reg$No_sim_Delta_hat,
    Sigma_hat = res_cov_reg$Sigma_hat,
    Dic_Delta_hat = res_cov_reg$Dic_Delta_hat,
    optimal_params = res_cov_reg$optimal_params,
    computational_time_per_node = res_cov_reg$computational_time_per_node,
    total_computational_time = difftime(Sys.time(), time.start, units = "mins"),
    scaling_params = res_cov_reg$scaling_params,
    dummy_params = res_cov_reg$dummy_params
  )
  return(list(results = results, additional_info = additional_info))
}

#################################################
## Predict Personalized Network for New Subject
#################################################

#' predict_personalized_network
#' 
#' This function predicts the personalized network structure for a new subject
#' based on their covariate values and the estimated covariate-specific effects
#' on the network structure.
#' 
#' @param Dic_Delta_hat Dictionary of covariate-specific adjacency matrices from GGReg estimation
#'   Should contain: $Baseline, and covariate-specific matrices
#' @param new_subject_covariates Data frame with covariates for the new subject (single row or multiple rows)
#' @param scaling_params List of scaling parameters (means and SDs) used for numerical covariates in training data
#' @param dummy_params List of dummy coding parameters used for categorical covariates in training data (including original column names and resulting dummy column names)
#' @param verbose Whether to print detailed information (default FALSE)
#' 
#' @return list of personalized adjacency matrices, one per subject
#'   
#' @examples
#' # Single subject prediction
#' new_subject <- data.frame(age = 45, sex = as.factor("Female"))
#' pred_net <- predict_personalized_network(
#'  Dic_Delta_hat = results$results$Dic_adj_matrics,
#'  new_subject_covariates = new_subject,
#'  scaling_params = results$additional_info$scaling_params,
#'  dummy_params = results$additional_info$dummy_params
#' )
#' 
#' # Multiple subjects
#' new_subjects <- data.frame(
#'   age = c(45, 30, 60),
#'   sex = as.factor(c("Female", "Male", "Female"))
#' )
#' pred_nets <- predict_personalized_network(
#'   Dic_Delta_hat = results$results$Dic_adj_matrics,
#'   new_subject_covariates = new_subjects,
#'   scaling_params = results$additional_info$scaling_params,
#'   dummy_params = results$additional_info$dummy_params
#' )
#' 
predict_personalized_network <- function(
  Dic_Delta_hat, 
  new_subject_covariates = NULL,
  scaling_params = NULL,
  dummy_params = NULL,
  verbose = FALSE) 
{ 
  if (!is.list(Dic_Delta_hat)) {
    stop("Dic_Delta_hat must be a list containing adjacency matrices")
  }

  baseline_key <- "Baseline"
  if (!(baseline_key %in% names(Dic_Delta_hat))) {
    stop("Dic_Delta_hat must contain the 'Baseline' network")
  }
  
  baseline_network <- Dic_Delta_hat[[baseline_key]]
  
  if (verbose) {
    cat("=== Predicting Personalized Network ===\n")
    cat("Baseline network:", baseline_key, "\n")
    cat("Available covariate effects:", setdiff(names(Dic_Delta_hat), baseline_key), "\n")
  }
  
  ## Case 1: No covariates - return baseline network
  
  if (is.null(new_subject_covariates) || is.null(scaling_params) || is.null(dummy_params)) {
    cat("No covariates provided - returning baseline network\n")
    cat("This represents the average network across all subjects\n")
    return(list(baseline_network))
  }else{
    n_subjects <- nrow(new_subject_covariates)
    if (verbose) {cat("Predicting for", n_subjects, "subject(s)\n")}
    missing_cols <- setdiff(dummy_params$original_colnames, colnames(new_subject_covariates))
    if(length(missing_cols) > 0) {
      warning("New subject covariates missing column: ", paste(missing_cols, collapse = ", "), ". Setting them to the mean/most numerous class.")
    }
    new_subject_covariates <- new_subject_covariates[, dummy_params$original_colnames, drop = FALSE]
    if (verbose && length(scaling_params) > 0) {
      cat("\nScaling parameters from training data:\n")
      for (col in names(scaling_params)) {
        cat("  ", col, ": mean =", round(scaling_params[[col]]$mean, 4), 
            ", sd =", round(scaling_params[[col]]$sd, 4), "\n")
      }
    }
    C_new <- data.frame(matrix(NA, n_subjects, length(dummy_params$dummy_colnames)))
    is_numeric <- rep(NA,length(dummy_params$dummy_colnames) )
    colnames(C_new) <- dummy_params$dummy_colnames
    names(is_numeric) <- dummy_params$dummy_colnames
    for (col in names(scaling_params)) {
      if (col %in% colnames(new_subject_covariates)) {
        C_new[[col]] <- (new_subject_covariates[[col]] - scaling_params[[col]]$mean) / 
                                      scaling_params[[col]]$sd
        is_numeric[col] <- TRUE
      }
    }
    if (verbose && length(scaling_params) > 0) {
      cat("\nScaled new subject covariates:\n")
      print(head(C_new))
    }
    col_to_dummy <- setdiff(dummy_params$original_colnames, dummy_params$dummy_colnames)
    for(col in col_to_dummy) {
      if (col %in% colnames(new_subject_covariates)) {
          dummies <- dummy_params$dummy_colnames[grepl(paste0("^", col), dummy_params$dummy_colnames)]
          for (dummy_col in dummies) {
            C_new[[dummy_col]] <- ifelse(new_subject_covariates[[col]] == sub(paste0("^", col),"",dummy_col), 1, 0)
            is_numeric[dummy_col] <- FALSE
          }
      }
    }
    if (verbose && length(col_to_dummy) > 0) {
      cat("\nDummy coded new subject covariates:\n")
      print(head(C_new))
    }

    # Predict pesonalized network for single or multiple subjects
    personalized_net <- list()
    for(j in 1:nrow(C_new)){
      if (verbose) {
        cat(" Estimating the personalized network for subject", j, "\n")
        }
      covariate_row <- C_new[j, ]
      signle_personalized_net <- baseline_network
      for (i in 1:length(covariate_row)) {
        covariate_name <- colnames(C_new)[i]
        covariate_value <- as.numeric(covariate_row[i])
        if (covariate_name %in% names(Dic_Delta_hat)) {
          if(is_numeric[covariate_name]){
            if (verbose) {
              cat("  Adding effect of numeric covariate", covariate_name, "with value", covariate_value, "\n")
            }
            signle_personalized_net <- signle_personalized_net + covariate_value * Dic_Delta_hat[[covariate_name]]
          } else {
            if (verbose && covariate_value == 1) {
              cat("  Adding effect of categorical covariate", covariate_name, "\n")
            }
            signle_personalized_net <- signle_personalized_net + covariate_value * Dic_Delta_hat[[covariate_name]]
          }
        } else {
          if (verbose) {
            warning("Covariate '", covariate_name, "' not found in Dic_Delta_hat. Skipping.")
          }
        }
        personalized_net[[j]] <- signle_personalized_net
      }
    }
    names(personalized_net) <- paste0("Subject_", 1:n_subjects)
    return(personalized_net)
  }
}

#################################################
## Pipeline Functions for the estimation
#################################################

#################################################
## STEP 1: Mean estimation
#################################################

#' GGReg_mean_estimation
#'
#' Estimates the mean structure in the GGReg model using the provided predictors and response.
#'
#' @param x Expression data - data frame; nrow: subjects; ncol: features
#' @param covariates Covariates to regress on - data frame; nrow: subjects; ncol: covariates). If NULL, no covariates are included in the model (default NULL)
#' @param lambda_mean Penalization term in Lasso for mean estimation (if mean_estimation = TRUE) - Numeric; If NULL, optimal lambda is selected via cross-validation (default NULL)
#' @param lambda_mean_type Type of lambda selection for mean estimation when using cross-validation for optimizing (if mean_estimation = TRUE) - Character; "1se" or "min" (default "1se")
#' @param verbose Logical indicating whether to print detailed progress messages during estimation - Logical; default FALSE
#'
#' @return A list containing:
#'   \item{z}{Estimated residuals.}
#'   \item{Cov_effect}{Estimated covariance effects from the mean regression.}
#'
#' @description
#' This function performs mean estimation for the GGReg model, returning the estimated covariance effects on the mean and the residuals used as input in the next step.

GGReg_mean_estimation <- function(
  x = NULL,
  covariates = NULL,
  lambda_mean = NULL,
  lambda_mean_type = "1se",
  verbose = FALSE) 
{
  if(is.null(x)){
    stop("Provide expression data")
  }
  n <- nrow(x)
  p <- ncol(x)
  z <- data.frame(matrix(NA, n, p))
  colnames(z) <-  colnames(x)
  if (n < 10) {
    warning("The sample size is too small! Network estimate may be unreliable!")
  }
  if (is.null(covariates)) {
    warning("No covariates included in the model")
    Cov_effect = matrix(0, p, 1)
    colnames(Cov_effect) <- "(Intercept)"
    for(i in 1:p){
      Cov_effect[i, ]= mean(x[, i])
      z[, i] = x[, i] - mean(x[, i])
    }
  } else {
    numeric_columns <- sapply(covariates, is.numeric)
    covariates[, numeric_columns] <- scale(covariates[, numeric_columns])
    C <- model.matrix( ~ ., covariates)
    X=x
    k <- ncol(C)
    Cov_effect = matrix(0, p, k)
    colnames(Cov_effect) <- colnames(C)
    C <- data.frame(C[,-1])
    colnames(C) <- colnames(Cov_effect)[-1]
    C <- as.matrix(C)
    for (i in 1:p) {
      Y <- matrix(X[, i], ncol = 1)
      mean_reg_coeff <- matrix(0, k, 1)
      if (k == 2){
        data <- data.frame(x = C, y = Y)
        colnames(data) <- c("x","y")
        mod <- lm(y ~ x, data)
        res.glm <- residuals(mod)
        Cov_effect[i, ] <- as.vector(coef(mod))
      } else {
        if (is.null(lambda_mean)) {
          mod <- cv.glmnet(x = C, y = Y)
          if (lambda_mean_type == "min"){
            mean_reg_coeff <- coef(mod, s = "lambda.min")[, 1]
          } else {
            mean_reg_coeff <- coef(mod, s = "lambda.1se")[, 1] 
          } 
          res.glm <- Y - cbind(1, C) %*% mean_reg_coeff
          Cov_effect[i, ] <- as.vector(mean_reg_coeff)
        } else {
          mod <- glmnet(x = C, y = Y, lambda = lambda_mean, relax = TRUE)
          mean_reg_coeff <- mod$beta[, 1]
          res.glm <- Y - cbind(1, C) %*% mean_reg_coeff
          Cov_effect[i, ] <- as.vector(mean_reg_coeff)
        }
      }
      z[, i] <- as.vector(res.glm)
    }
  }

  return(list(
    z = z,
    Cov_effect = Cov_effect
  ))
}

#################################################
## STEP 2: Precision Matrix Estimation -  SLURM-based parallelization over nodes
#################################################

#' GGReg_cov_estimation_SLURM
#'
#' Estimates the covariance structure in the GGReg model for each node, designed for distributed computation (e.g., SLURM clusters).
#'
#' @param Z0 A data frame of residual values from the mean estimation - data frame; nrow: subjects; ncol: features
#' @param known_ppi Previously available knowledge on interaction - matrix of values [0,1] that describe the confeidence of interaction existing; ncol: features, nrow: features. If NULL, no prior knowledge is used (default NULL)
#' @param covariates Covariates to regress on - data frame; nrow: subjects; ncol: covariates). If NULL, no covariates are included in the model (default NULL)
#' @param scr Logical indicating whether to perform screening based on correlation to speed up estimation - Logical; default TRUE
#' @param gamma Pearson correlation threshold for screening (if scr = TRUE) - Numeric; If NULL 20% quantile of absolute correlations (default NULL)
#' @param lambda_prec Penalization term in Lasso for precision estimation - Numeric; If NULL, optimal lambda is selected via cross-validation (default NULL)
#' @param lambda_prec_type Type of lambda selection for precision estimation when using cross-validation for optimizing - Character; "1se" or "min" (default "1se")
#' @param tune_hyperparams Logical indicating whether to perform hyperparameter tuning (weight and aloha) for precision estimation - Logical; default TRUE
#' @param asparse_grid Grid of alpha values to consider for hyperparameter tuning (if tune_hyperparams = TRUE) - Numeric vector; default c(0.5, 0.75, 0.9, 0.95)
#' @param weight_grid Grid of weight values to consider for hyperparameter tuning (if tune_hyperparams = TRUE) - Numeric vector; default c(0.8, 1.0, 1.1, 1.3, 1.5)
#' @param random_hyper_search Logical indicating whether to perform random hyperparameter search instead of grid search for tuning (if tune_hyperparams = TRUE) - Logical; default FALSE
#' @param p.rand.hyper Proportion of random hyperparameter combinations to try if random_hyper_search = TRUE (if tune_hyperparams = TRUE) - Numeric; default NULL setting which means 50% of the combinations will be tried
#' @param K Number of folds for cross-validation when computing BIC for hyperparameter tuning (if tune_hyperparams = TRUE) - Numeric; default 5
#' @param slurm_script_path Path to the SLURM script to be used for parallelization (if use_slurm = TRUE) - Character; default "./slurm_ggReg_node.sbatch"
#' @param output_path Directory where SLURM job outputs will be saved (if use_slurm = TRUE) - Character; default "./results/"
#' @param name_output Base name for SLURM job output files (if use_slurm = TRUE) - Character; default "ggReg_result"
#' @param symm_method Method for symmetrization of the estimated precision matrix - Character; "OR" or "AND" (default "OR")
#' @param verbose Logical indicating whether to print detailed progress messages during estimation - Logical; default FALSE
#' 
#' @return A list containing:
#'   \item{deltas}{Estimated delta parameters for each node.}
#'   \item{Sigma_hat}{Estimated covariance matrix.}
#'   \item{Dic_Delta_hat}{Dictionary of estimated Delta matrices for each node.}
#'   \item{optimal_params}{Optimal hyperparameters per node.}
#'
#' @description
#' This function performs node-wise covariance estimation for the GGReg model, supporting distributed execution (e.g., on SLURM clusters), and returns estimated parameters and matrices for each node.

GGReg_cov_estimation_SLURM <- function(
  Z0,
  known_ppi = NULL,
  covariates = NULL,
  scr = FALSE,
  gamma = NULL,
  lambda_prec = NULL,
  lambda_prec_type = "min",
  tune_hyperparams = TRUE,
  asparse_grid = c(0.5, 0.75, 0.9, 0.95),
  weight_grid = c(0.8, 1.0, 1.1, 1.3, 1.5),
  random_hyper_search = FALSE,
  p.rand.hyper = NULL,
  K = 5,
  slurm_script_path = "./slurm_ggReg_node.sbatch",
  output_path = "./results/",
  name_output = "ggReg_result",
  symm_method ="OR",
  verbose = FALSE) 
{
  n <- nrow(Z0)
  p <- ncol(Z0)
  if (verbose) {
    cat("Setting up SLURM array job for", p, "nodes\n")
    if (tune_hyperparams) {
      cat("Each node will tune hyperparameters: asparse =", length(asparse_grid), 
          "values, weight =", length(weight_grid), "values\n")
    }
  }
  if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
  }
  # Save all input data for each SLURM jobs
  input_data_path <- paste0(output_path, "input_data_nodes.rda")
  save(Z0, known_ppi, covariates, scr, gamma, lambda_prec, lambda_prec_type, 
       tune_hyperparams, asparse_grid, weight_grid, random_hyper_search, p.rand.hyper, K,
       file = input_data_path)
  if (verbose) {
    cat("Saved input data to:", input_data_path, "\n")
  }
  # Submit SLURM array job (one task per node i=1:p)
  job_cmd <- paste("sbatch --array=1-", p, " ", slurm_script_path, " ",
                  input_data_path, " ", output_path, " ", name_output, sep = "")
  cat("\n", rep("=", 60), "\n", sep = "")
  cat("PAUSED FOR MANUAL SUBMISSION\n")
  cat(rep("=", 60), "\n")
  cat("Set the right directory where your functions are stored: cd /your_folder/ \n")
  cat("Please run the following command manually:\n\n")
  cat("  ", job_cmd, "\n\n")
  cat("After running the command, you will get output like:\n")
  cat("  'Submitted batch job 12345'\n\n")
  cat("Monitor the job array via 'squeue -u $USER'\n")
  cat("Once the job is complete and enter the job ID below and press ENTER to continue:\n")
  cat(rep("=", 60), "\n")
  # Wait for user input
  job_id <- readline(prompt = "Job ID: ")
  # Validate input
  if (is.null(job_id) || job_id == "" || is.na(job_id)) {
    stop("No job ID provided. Please restart the function with a valid job ID.")
  }
  # Remove any whitespace
  job_id <- trimws(job_id)
  if (verbose) {
    cat("Received job ID:", job_id, "\n")
    cat("Collecting the results...\n")
  }
  # Collect results from all nodes
  result <- collect_node_results(p, output_path, name_output, symm_method,  verbose = verbose)
  return(result)
}


#################################################
## FALLBACK 1: Sequential GGReg_cov_estimation (no parallelization)
#################################################

#' GGReg_cov_estimation_sequential
#'
#' Estimates the covariance structure in the GGReg model for each node, designed for sequential computation (WARNING: it could be quite slow and may not be feasible for large p).
#'
#' @param Z0 A data frame of residual values from the mean estimation - data frame; nrow: subjects; ncol: features
#' @param known_ppi Previously available knowledge on interaction - matrix of values [0,1] that describe the confeidence of interaction existing; ncol: features, nrow: features. If NULL, no prior knowledge is used (default NULL)
#' @param covariates Covariates to regress on - data frame; nrow: subjects; ncol: covariates). If NULL, no covariates are included in the model (default NULL)
#' @param scr Logical indicating whether to perform screening based on correlation to speed up estimation - Logical; default TRUE
#' @param gamma Pearson correlation threshold for screening (if scr = TRUE) - Numeric; If NULL 20% quantile of absolute correlations (default NULL)
#' @param lambda_prec Penalization term in Lasso for precision estimation - Numeric; If NULL, optimal lambda is selected via cross-validation (default NULL)
#' @param lambda_prec_type Type of lambda selection for precision estimation when using cross-validation for optimizing - Character; "1se" or "min" (default "1se")
#' @param tune_hyperparams Logical indicating whether to perform hyperparameter tuning (weight and aloha) for precision estimation - Logical; default TRUE
#' @param asparse_grid Grid of alpha values to consider for hyperparameter tuning (if tune_hyperparams = TRUE) - Numeric vector; default c(0.5, 0.75, 0.9, 0.95)
#' @param weight_grid Grid of weight values to consider for hyperparameter tuning (if tune_hyperparams = TRUE) - Numeric vector; default c(0.8, 1.0, 1.1, 1.3, 1.5)
#' @param random_hyper_search Logical indicating whether to perform random hyperparameter search instead of grid search for tuning (if tune_hyperparams = TRUE) - Logical; default FALSE
#' @param p.rand.hyper Proportion of random hyperparameter combinations to try if random_hyper_search = TRUE (if tune_hyperparams = TRUE) - Numeric; default NULL setting which means 50% of the combinations will be tried
#' @param K Number of folds for cross-validation when computing BIC for hyperparameter tuning (if tune_hyperparams = TRUE) - Numeric; default 5
#' @param output_path Directory where SLURM job outputs will be saved (if use_slurm = TRUE) - Character; default "./results/"
#' @param name_output Base name for SLURM job output files (if use_slurm = TRUE) - Character; default "ggReg_result"
#' @param symm_method Method for symmetrization of the estimated precision matrix - Character; "OR" or "AND" (default "OR")
#' @param verbose Logical indicating whether to print detailed progress messages during estimation - Logical; default FALSE
#' 
#' @return A list containing:
#'   \item{deltas}{Estimated delta parameters for each node.}
#'   \item{Sigma_hat}{Estimated covariance matrix.}
#'   \item{Dic_Delta_hat}{Dictionary of estimated Delta matrices for each node.}
#'   \item{optimal_params}{Optimal hyperparameters per node.}
#'
#' @description
#' This function performs node-wise covariance estimation for the GGReg model sequentially and returns estimated parameters and matrices for each node.


GGReg_cov_estimation_sequential <- function(
  Z0,
  known_ppi = NULL,
  covariates = NULL,
  scr = TRUE,
  gamma = NULL,
  lambda_prec = NULL,
  lambda_prec_type = "min",
  tune_hyperparams = TRUE,
  asparse_grid = c(0.5, 0.75, 0.9, 0.95),
  weight_grid = c(0.8, 1.0, 1.1, 1.3, 1.5),
  random_hyper_search = FALSE,
  p.rand.hyper = NULL,
  K = 5,
  output_path = "./results/",
  name_output = "ggReg_result",
  symm_method ="OR",
  verbose = FALSE) 
{
  n <- nrow(Z0)
  p <- ncol(Z0)

  if (verbose && tune_hyperparams) {
    cat("Each node will tune hyperparameters: asparse =", length(asparse_grid), 
        "values, weight =", length(weight_grid), "values\n")
  }
  
  # Create output directory
  if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
  }
  
  # Save all input data for individual node runs
  input_data_path <- paste0(output_path, "input_data_nodes.rda")
  save(Z0, known_ppi, covariates, scr, gamma, lambda_prec, lambda_prec_type, 
       tune_hyperparams, asparse_grid, weight_grid, random_hyper_search, p.rand.hyper, K, 
       file = input_data_path)
  
  if (verbose) {
    cat("Saved input data to:", input_data_path, "\n")
    cat("Starting sequential processing of", p, "nodes...\n")
  }

  #################################################
  ##  Sequential processing of nodes
  #################################################
  
  for (i in 1:p) {
    if (verbose && (i %% 10 == 0 || i == 1)) {
      cat("Processing node", i, "of", p, "\n")
    }
    GGReg_cov_single_node_processing(
      input_data_path = input_data_path,
      node_index = i,
      output_path = output_path,
      name_output = name_output,
      verbose = verbose
    )
  }
  if (verbose) {
    cat("Done processing all nodes sequentially\n")
    cat("Collecting the results...\n")
  }
  result <- collect_node_results(p, output_path, name_output, symm_method,  verbose = verbose)
  return(result)
}

#################################################
## FALLBACK 1: Parallel GGReg_cov_estimation (based on forach in R, not distributed, but can speed up for moderate p)
#################################################

#' GGReg_cov_estimation_parallel
#'
#' Estimates the covariance structure in the GGReg model for each node, designed for parallel computation executed in R via foreach.
#'
#' @param Z0 A data frame of residual values from the mean estimation - data frame; nrow: subjects; ncol: features
#' @param known_ppi Previously available knowledge on interaction - matrix of values [0,1] that describe the confeidence of interaction existing; ncol: features, nrow: features. If NULL, no prior knowledge is used (default NULL)
#' @param covariates Covariates to regress on - data frame; nrow: subjects; ncol: covariates). If NULL, no covariates are included in the model (default NULL)
#' @param scr Logical indicating whether to perform screening based on correlation to speed up estimation - Logical; default TRUE
#' @param gamma Pearson correlation threshold for screening (if scr = TRUE) - Numeric; If NULL 20% quantile of absolute correlations (default NULL)
#' @param lambda_prec Penalization term in Lasso for precision estimation - Numeric; If NULL, optimal lambda is selected via cross-validation (default NULL)
#' @param lambda_prec_type Type of lambda selection for precision estimation when using cross-validation for optimizing - Character; "1se" or "min" (default "1se")
#' @param tune_hyperparams Logical indicating whether to perform hyperparameter tuning (weight and aloha) for precision estimation - Logical; default TRUE
#' @param asparse_grid Grid of alpha values to consider for hyperparameter tuning (if tune_hyperparams = TRUE) - Numeric vector; default c(0.5, 0.75, 0.9, 0.95)
#' @param weight_grid Grid of weight values to consider for hyperparameter tuning (if tune_hyperparams = TRUE) - Numeric vector; default c(0.8, 1.0, 1.1, 1.3, 1.5)
#' @param random_hyper_search Logical indicating whether to perform random hyperparameter search instead of grid search for tuning (if tune_hyperparams = TRUE) - Logical; default FALSE
#' @param p.rand.hyper Proportion of random hyperparameter combinations to try if random_hyper_search = TRUE (if tune_hyperparams = TRUE) - Numeric; default NULL setting which means 50% of the combinations will be tried
#' @param K Number of folds for cross-validation when computing BIC for hyperparameter tuning (if tune_hyperparams = TRUE) - Numeric; default 5
#' @param output_path Directory where SLURM job outputs will be saved (if use_slurm = TRUE) - Character; default "./results/"
#' @param name_output Base name for SLURM job output files (if use_slurm = TRUE) - Character; default "ggReg_result"
#' @param symm_method Method for symmetrization of the estimated precision matrix - Character; "OR" or "AND" (default "OR")
#' @param verbose Logical indicating whether to print detailed progress messages during estimation - Logical; default FALSE
#' 
#' @return A list containing:
#'   \item{deltas}{Estimated delta parameters for each node.}
#'   \item{Sigma_hat}{Estimated covariance matrix.}
#'   \item{Dic_Delta_hat}{Dictionary of estimated Delta matrices for each node.}
#'   \item{optimal_params}{Optimal hyperparameters per node.}
#'
#' @description
#' This function performs node-wise covariance estimation for the GGReg model in parallel exploiting R parallelization with foreach and returns estimated parameters and matrices for each node. Note the function automatically detects the core and uses all cores minus one for parallelization. For large p, this function may still be computationally intensive.

GGReg_cov_estimation_parallel <- function(
  Z0,
  known_ppi = NULL,
  covariates = NULL,
  scr = TRUE,
  gamma = NULL,
  lambda_prec = NULL,
  lambda_prec_type = "min",
  tune_hyperparams = TRUE,
  asparse_grid = c(0.5, 0.75, 0.9, 0.95),
  weight_grid = c(0.8, 1.0, 1.1, 1.3, 1.5),
  random_hyper_search = FALSE,
  p.rand.hyper = NULL,
  K = 5,
  output_path = "./results/",
  name_output = "ggReg_result",
  symm_method ="OR",
  verbose = FALSE) 
{
  n <- nrow(Z0)
  p <- ncol(Z0)

  if (verbose && tune_hyperparams) {
    cat("Each node will tune hyperparameters: asparse =", length(asparse_grid), 
        "values, weight =", length(weight_grid), "values\n")
  }
  
  # Create output directory
  if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
  }

  output_path <<- output_path
  name_output <<- name_output
  
  # Save all input data for individual node runs
  input_data_path <<- paste0(output_path, "input_data_nodes.rda")
  save(Z0, known_ppi, covariates, scr, gamma, lambda_prec, lambda_prec_type, 
       tune_hyperparams, asparse_grid, weight_grid, random_hyper_search, p.rand.hyper, K,
       file = input_data_path)
  
  if (verbose) {
    cat("Saved input data to:", input_data_path, "\n")
    cat("Starting parallel processing of", p, "nodes...\n")
  }

  #################################################
  ## Parallel processing of nodes with foreach
  #################################################
  
  num_cores <- detectCores() - 1
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  clusterExport(cl, c("GGReg_cov_single_node_processing", "input_data_path", "output_path", "name_output"))
  
  # Load libraries on workers
  clusterEvalQ(cl, {
    library(sparsegl)
    library(glmnet)
  })
  
  foreach(i = 1:p, .packages = c('glmnet', 'sparsegl')) %dopar% {
    GGReg_cov_single_node_processing(
      input_data_path = input_data_path,
      node_index = i,
      output_path = output_path,
      name_output = name_output,
      verbose = verbose
    )
  }
  stopCluster(cl)
  
  if (verbose) {
    cat("Done processing all nodes in parallel\n")
    cat("Collecting the results...\n")
  }
  
  result <- collect_node_results(p, output_path, name_output, symm_method, verbose = verbose)
  
  return(result)
}


#################################################
## Helper Functions 
#################################################

#' GGReg_cov_single_node_processing
#'
#' Individual node estimation function for the neighboorhood selection in the GGReg model with optional hyperparameter tuning, designed to be called for each node in parallel (e.g., via SLURM or foreach) or sequentially.
#'
#' @param input_data_path Path to the input data file containing Z0, known_ppi, covariates, scr, gamma, lambda_prec, lambda_prec_type, tune_hyperparams, asparse_grid, weight_grid, random_hyper_search, p.rand.hyper, K (saved as an RDA file) - Character
#' @param node_index Index of the node to be processed - Numeric (1-based index from 1 to p)
#' @param output_path Directory where the output for the node will be saved - Character
#' @param name_output Base name for the output file for the node - Character
#' @param verbose Logical indicating whether to print detailed progress messages during estimation - Logical; default FALSE
#'
#' @return
#'
#' @description
#' This function performs the neighboorhood estimation for a single node in the GGReg model, including setting up the design matrix, applying screening if specified, and running the regression to estimate the neighborhood of the node. It saves the results for the node to a specified output path.

GGReg_cov_single_node_processing <- function(
  input_data_path,
  node_index,
  output_path,
  name_output,
  verbose = FALSE)
{
  #rm(list = ls(all = TRUE))
  time.start <- Sys.time()
  if(verbose){
    cat("=== GGReg Node Estimation Job ===\n")
    cat("Input data:", input_data_path, "\n")
    cat("Processing node:", node_index, "\n")
    cat("Output path:", output_path, "\n")
    cat("Output name:", name_output, "\n")
  }

  #################################################
  ## Load Input Data
  #################################################
  if(verbose){cat("Loading input data...\n")}
  load(input_data_path)
  # This loads: Z0, known_ppi, covariates, scr, gamma, lambda_prec, lambda_prec_type, 
  #             tune_hyperparams, asparse_grid, weight_grid, K
  n <- nrow(Z0)
  p <- ncol(Z0)
  i <- node_index  # Current node index
  if(verbose){
    cat("Data dimensions: n =", n, ", p =", p, "\n")
    cat("Processing node", i, "of", p, "\n")
  }
  if (i > p) {
    stop("Node index ", i, " exceeds number of nodes ", p)
  }
  #################################################
  ## Setup Covariates and Design Matrix
  #################################################
  scaling_params <- list() 
  dummy_params <- NULL

  if (is.null(covariates)) {
    C <- data.frame(Zeros = rep(0, n))
    iU <- data.frame(rep(1, n))
    colnames(iU) <- "(Intercept)"
    q <- 0
  } else {
    numeric_columns <- sapply(covariates, is.numeric)
    covariates[, numeric_columns] <- scale(covariates[, numeric_columns])
    for (col in colnames(covariates)[numeric_columns]) {
        scaling_params[[col]] <- list(
          mean = mean(covariates[[col]], na.rm = TRUE),
          sd = sd(covariates[[col]], na.rm = TRUE)
        )
    }
    C <- model.matrix(~ ., covariates)
    dummy_params <- list(
      original_colnames = colnames(covariates),
      dummy_colnames = colnames(C)[-1]  # Exclude intercept
    )
    q <- ncol(C) - 1
    iU <- C
  }

  #################################################
  ## Setup Prior Knowledge and Screening
  #################################################

  if (!is.null(known_ppi)) {
    one <- matrix(FALSE, p, p)
    one[known_ppi > 0] <- TRUE
    rownames(one) <- colnames(Z0)
    penaltyfactor <- 1 - known_ppi
  } else {
    one <- matrix(FALSE, p, p)
    rownames(one) <- colnames(Z0)
    penaltyfactor <- matrix(1, p, p)
  }

  if (scr) {
    if(verbose){cat("Computing correlation matrix for screening...\n")}
    protdata_matrix <- as.matrix(Z0)
    Scor <- cor(protdata_matrix)
    scr_index <- matrix(1, p, p)
    if (is.null(gamma)) {
      gamma <- quantile(abs(Scor), probs = 0.2)
    }
    scr_index[abs(Scor) <= gamma] <- 0
    diag(scr_index) <- 0
    colnames(scr_index) <- colnames(Z0)
    rownames(scr_index) <- colnames(Z0)
    if(verbose){cat("Screening threshold (gamma):", gamma, "\n")}
  }

  #################################################
  ## Build Design Matrix
  #################################################

  if(verbose){cat("Building design matrix...\n")}
  interM <- data.frame(Intercept = rep(1, n))
  groups <- c(0, rep(1, ncol(Z0)))
  w_group <- c(0)

  for (j in 1:(q + 1)) {
    product <- Z0 * iU[, j]
    if (j != 1) {
      original_colnames <- colnames(product)
      iU_name <- colnames(iU)[j]
      new_colnames <- paste(iU_name, ":", original_colnames, sep = "")
      colnames(product) <- new_colnames
      groups <- c(groups, rep(j, ncol(product)))
      w_group <- c(w_group, 1)
    }
    interM <- cbind(interM, product)
  }
  # Save names of covariates for later use in output
  if (!is.null(covariates)) {
    U_names <- colnames(C)[-1]
  } else {
    U_names <- c()
  }
  U_names <- c("Intercept", "Baseline", U_names)

  if(verbose){cat("Design matrix built. Dimensions:", dim(interM), "\n")}

  #################################################
  ## Prepare Node-Specific Data
  #################################################

  Y <- matrix(Z0[, i], ncol = 1)
  protname <- colnames(Z0)[i]

  # Select features (exclude self-loops and screened features)
  pattern <- paste0(":", protname, "\\b")
  indexes_to_select <- !grepl(pattern, colnames(interM))
  indexes_to_select[colnames(interM) == protname] <- FALSE
  indexes_to_select_penalty <- rep(TRUE, ncol(Z0))
  indexes_to_select_penalty[colnames(Z0) == protname] <- FALSE

  if (scr) {
    zero_col_indices <- which(scr_index[i, ] == 0)
    zero_col_names <- colnames(scr_index)[zero_col_indices]
    for (j in zero_col_names) {
      pattern <- paste0(":", j, "\\b")
      indexes_to_select[grepl(pattern, colnames(interM))] <- FALSE
      indexes_to_select[colnames(interM) == j] <- FALSE
      indexes_to_select_penalty[colnames(Z0) == j] <- FALSE
    }
  }

  indexes_to_select_reg <- indexes_to_select
  indexes_to_select_reg[1] <- FALSE  # Remove intercept
  Xmat <- interM[, indexes_to_select_reg]
  Xgroup <- groups[indexes_to_select_reg]

  if(verbose){cat("Selected", ncol(Xmat), "features for node", i, "\n")}

  if(is.null(covariates)){
    if (is.null(known_ppi)) {
      set.seed(2024)
      if (is.null(lambda_prec)) {
        mod <- cv.glmnet(x = as.matrix(Xmat), y = Y, alpha = 1, family = c("gaussian"))
      } else{
        lambda_seq <- c(lambda_prec-0.01*lambda_prec, lambda_prec, lambda_prec+ 0.01*lambda_prec)
        mod <- cv.glmnet(x = as.matrix(Xmat), y = Y,lambda = lambda_seq, alpha = 1, family = c("gaussian"))
      }
    }else{
      w_sparsity <- rep(1, ncol(Xmat))
      temp_w_sparsity_g1 <- penaltyfactor[i, indexes_to_select_penalty]
      w_sparsity[Xgroup == 1] <- temp_w_sparsity_g1
      set.seed(2024)
      if (is.null(lambda_prec)) {
        mod <- cv.sparsegl(x = as.matrix(Xmat), y = Y, pf_sparse = w_sparsity, 
                          asparse = 1, intercept = FALSE, nfolds = K)
      } else {
        lambda_seq <- c(lambda_prec-0.01*lambda_prec, lambda_prec, lambda_prec+ 0.01*lambda_prec)
        mod <- cv.sparsegl(x = as.matrix(Xmat), y = Y, pf_sparse = w_sparsity, 
                          asparse = 1, lambda = lambda_seq, 
                          intercept = FALSE, nfolds = K)
      }
    }

    if (lambda_prec_type == "min") {
      prec_reg_coeff <- coef(mod, s = "lambda.min")[, 1]
      best_lambda_prec <- mod$lambda.min
    } else {
      prec_reg_coeff <- coef(mod, s = "lambda.1se")[, 1]
      best_lambda_prec <- mod$lambda.1se
    }

    # Compute BIC for this model
    res.model <- Y - cbind(1, as.matrix(Xmat)) %*% prec_reg_coeff
    df_j <- sum(prec_reg_coeff != 0)
    sigma <- sum(res.model^2) / (n - df_j)

    best_params <- list(best_lambda_prec = best_lambda_prec, bic_score = bic_score)

    # BIC = n * log(RSS/n) + log(n) * df
    bic_score <- n * log(sum(res.model^2) / n) + log(n) * df_j

    Prec_reg_matrix_row <- rep(0, ncol(interM))
    Prec_reg_matrix_row[indexes_to_select] <- as.vector(prec_reg_coeff)

    No_sim_Delta_hat_row <- rep(0, ncol(interM))
    No_sim_Delta_hat_row[indexes_to_select] <- -sigma * as.vector(prec_reg_coeff)

    if(verbose){
      cat("  Non-zero coefficients:", sum(abs(prec_reg_coeff) > 0), "\n")
      cat("  Residual variance:", sigma, "\n")}

    node_result <- list(
      node_index = i,
      Prec_reg_matrix_row = Prec_reg_matrix_row,
      No_sim_Delta_hat_row = No_sim_Delta_hat_row,
      sigma = sigma,
      optimal_params = best_params,
      computation_time = Sys.time() - time.start,
      # Save metadata for result collection
      col_names = colnames(interM),
      row_names = colnames(Z0),
      groups = groups,
      U_names = U_names,
      scaling_params = scaling_params,
      dummy_params = dummy_params
    )

    if(verbose){cat("Saving results for node", i, "...\n")}

    output_file <- paste0(output_path, name_output, "_node_", i, ".rda")
    save(node_result, file = output_file)

  }else{
    if (tune_hyperparams && !(length(asparse_grid) == 1 && length(weight_grid) == 1)) {
      if(verbose){
      cat("Starting hyperparameter tuning...\n")
      cat("Grid: asparse =", length(asparse_grid), "values, weight =", length(weight_grid), "values\n")}
      
      best_bic <- Inf
      best_params <- NULL
      best_result <- NULL

      param_combinations <- expand.grid(asparse = asparse_grid, weight = weight_grid)

      if(random_hyper_search){
        if(is.null(p.rand.hyper)){
          p.rand.hyper <- 0.5
        }
        set.seed(234+10*i)
        param_combinations <- param_combinations[sample(1:nrow(param_combinations), size = floor(p.rand.hyper * nrow(param_combinations))), ]
      }
      for (idx in 1:nrow(param_combinations)) {
        asparse_val <- param_combinations$asparse[idx]
        weight_val <- param_combinations$weight[idx]
        
        if (idx %% 5 == 1 && verbose) {
          cat("Testing asparse =", asparse_val, ", weight =", weight_val, 
              "(", idx, "of", nrow(param_combinations), ")\n")
        }
      
        # Setup penalty weights
        w_sparsity <- rep(weight_val, ncol(Xmat))
        temp_w_sparsity_g1 <- penaltyfactor[i, indexes_to_select_penalty]
        #temp_w_sparsity_g1[!one[i,indexes_to_select_penalty]] <- weight_val
        w_sparsity[Xgroup == 1] <- temp_w_sparsity_g1
        #w_group[which(w_group == 1)] <- sqrt(sum(Xgroup == 1))
        set.seed(2024)
        if (is.null(lambda_prec)) {
          mod <- cv.sparsegl(x = as.matrix(Xmat), y = Y, group = Xgroup, 
                            pf_group = w_group, pf_sparse = w_sparsity, 
                            asparse = asparse_val, intercept = FALSE, nfolds = K)
        } else {
          lambda_seq <- c(lambda_prec-0.01*lambda_prec, lambda_prec, lambda_prec+ 0.01*lambda_prec)
          mod <- cv.sparsegl(x = as.matrix(Xmat), y = Y, group = Xgroup, 
                            pf_group = w_group, pf_sparse = w_sparsity, 
                            asparse = asparse_val, lambda = lambda_seq, 
                            intercept = FALSE, nfolds = K)
        }
        if (lambda_prec_type == "min") {
          prec_reg_coeff <- coef(mod, s = "lambda.min")[, 1]
          best_lambda_prec <- mod$lambda.min
        } else {
          prec_reg_coeff <- coef(mod, s = "lambda.1se")[, 1]
          best_lambda_prec <- mod$lambda.1se
        }
        # Compute BIC for this parameter combination
        res.model <- Y - cbind(1, as.matrix(Xmat)) %*% prec_reg_coeff
        df_j <- sum(prec_reg_coeff != 0)
        sigma <- sum(res.model^2) / (n - df_j)
      
        # BIC = n * log(RSS/n) + log(n) * df
        bic_score <- n * log(sum(res.model^2) / n) + log(n) * df_j
      
        if (bic_score < best_bic) {
          best_bic <- bic_score
          best_params <- list(best_lambda_prec = best_lambda_prec, asparse = asparse_val, weight = weight_val, bic_score = bic_score)
          best_result <- list(prec_reg_coeff = prec_reg_coeff, sigma = sigma, mod = mod)

          if(verbose){

          cat("Current optimal parameters for node", i, ":\n")
          cat("  asparse =", best_params$asparse, "\n")
          cat("  weight =", best_params$weight, "\n")
          cat("  BIC score =", best_params$bic_score, "\n")

          cat("Preparing current best output for node", i, "...\n")}

          # Build full coefficient vector for this node
          Prec_reg_matrix_row <- rep(0, ncol(interM))
          Prec_reg_matrix_row[indexes_to_select] <- as.vector(prec_reg_coeff)

          No_sim_Delta_hat_row <- rep(0, ncol(interM))
          No_sim_Delta_hat_row[indexes_to_select] <- -sigma * as.vector(prec_reg_coeff)

          if(verbose){
          cat("  Non-zero coefficients:", sum(abs(prec_reg_coeff) > 0), "\n")
          cat("  Residual variance:", sigma, "\n")}

          node_result <- list(
            node_index = i,
            Prec_reg_matrix_row = Prec_reg_matrix_row,
            No_sim_Delta_hat_row = No_sim_Delta_hat_row,
            sigma = sigma,
            optimal_params = best_params,
            computation_time = Sys.time() - time.start,
            # Save metadata for result collection
            col_names = colnames(interM),
            row_names = colnames(Z0),
            groups = groups,
            U_names = U_names,
            scaling_params = scaling_params,
            dummy_params = dummy_params
          )

          cat("Saving best results so far for node", i, "...\n")

          output_file <- paste0(output_path, name_output, "_node_", i, ".rda")
          save(node_result, file = output_file)

        }
      
      }
      if (is.null(best_result)) {
        stop("No valid hyperparameter combination found for node ", i)
      }    
    } else {
      # No hyperparameter tuning - use first values from grids
      asparse_val <- asparse_grid[1]
      weight_val <- weight_grid[1]
      if(verbose){
      cat("Using fixed hyperparameters: asparse =", asparse_val, ", weight =", weight_val, "\n")}
    
      # Setup penalty weights
      w_sparsity <- rep(weight_val, ncol(Xmat))
      w_sparsity[Xgroup == 1] <- penaltyfactor[i, indexes_to_select_penalty]
      w_group[which(w_group == 1)] <- sqrt(sum(Xgroup == 1))
      set.seed(2024)
      if (is.null(lambda_prec)) {
        mod <- cv.sparsegl(x = as.matrix(Xmat), y = Y, group = Xgroup, 
                          pf_group = w_group, pf_sparse = w_sparsity, 
                          asparse = asparse_val, intercept = FALSE)
      } else {
        lambda_seq <- c(lambda_prec-0.01*lambda_prec, lambda_prec, lambda_prec+ 0.01*lambda_prec)
        mod <- cv.sparsegl(x = as.matrix(Xmat), y = Y, group = Xgroup, 
                          pf_group = w_group, pf_sparse = w_sparsity, 
                          asparse = asparse_val, lambda = lambda_seq, 
                          intercept = FALSE)
      }
      if (lambda_prec_type == "min") {
        prec_reg_coeff <- coef(mod, s = "lambda.min")[, 1]
        best_lambda_prec <- mod$lambda.min
      } else {
        prec_reg_coeff <- coef(mod, s = "lambda.1se")[, 1]
        best_lambda_prec <- mod$lambda.1se
      }

      # Compute residuals and variance
      res.model <- Y - cbind(1, as.matrix(Xmat)) %*% prec_reg_coeff
      df_j <- sum(prec_reg_coeff != 0)
      sigma <- sum(res.model^2) / (n - df_j)
    
      # BIC for record keeping
      bic_score <- n * log(sum(res.model^2) / n) + log(n) * df_j
    
      optimal_params <- list(best_lambda_prec = best_lambda_prec, asparse = asparse_val, weight = weight_val, bic_score = bic_score)

      if(verbose){cat("Preparing output for node", i, "...\n")}

      # Build full coefficient vector for this node
      Prec_reg_matrix_row <- rep(0, ncol(interM))
      Prec_reg_matrix_row[indexes_to_select] <- as.vector(prec_reg_coeff)

      No_sim_Delta_hat_row <- rep(0, ncol(interM))
      No_sim_Delta_hat_row[indexes_to_select] <- -sigma * as.vector(prec_reg_coeff)

      if(verbose){
      cat("Node", i, "estimation complete:\n")
      cat("  Non-zero coefficients:", sum(abs(prec_reg_coeff) > 0), "\n")
      cat("  Residual variance:", sigma, "\n")
      }

      if(verbose){cat("Saving results for node", i, "...\n")}

      node_result <- list(
        node_index = i,
        Prec_reg_matrix_row = Prec_reg_matrix_row,
        No_sim_Delta_hat_row = No_sim_Delta_hat_row,
        sigma = sigma,
        optimal_params = optimal_params,
        computation_time = difftime(Sys.time(), time.start, units = "mins"),
        # Save metadata for result collection
        col_names = colnames(interM),
        row_names = colnames(Z0),
        groups = groups,
        U_names = U_names,
        scaling_params = scaling_params,
        dummy_params = dummy_params
      )

      output_file <- paste0(output_path, name_output, "_node_", i, ".rda")
      save(node_result, file = output_file)
    }
  }
  #################################################
  ## Prepare Output for Node i
  #################################################

  cat("Final results saved to:", output_file, "\n")
  cat("Computation time for node", i, ":", difftime(Sys.time(), time.start, units = "mins"), "\n")
  cat("=== Node", i, "Complete ===\n")
}

#' collect_node_results
#'
#' Collects and aggregates the results from node-wise neighboorhod selection in the GGReg model.
#'
#' @param p Number of nodes (features) for which results are expected - Numeric
#' @param output_path Directory where SLURM job outputs are saved - Character
#' @param name_output Base name for SLURM job output files - Character
#' @param symm_method Method for symmetrization of the estimated precision matrix - Character; "OR" or "AND"
#' @param verbose Logical indicating whether to print detailed progress messages during result collection - Logical; default FALSE
#'
#' @return A list containing:
#'   \item{Prec_reg_matrix}{Aggregated coefficient regression matrix combining results from all nodes.}
#'   \item{No_sim_Delta_hat}{Aggregated non-symmetrized predicsion matrix combining results from all nodes.}
#'   \item{Sigma_hat}{Aggregated vector of estimated variances for each node.}
#'   \item{Dic_Delta_hat}{Dictionary of symmetrized precision matrices estimation for each node.}
#'   \item{Dic_adj_matrics}{Dictionary of symmetrized weighted adjacency matrices estimation for each node reporting partial correlation between nodes.}
#'   \item{optimal_params}{Data frame of optimal hyperparameters (asparse, weight) and BIC score for each node.}
#'   \item{computational_time_per_node}{Vector of computational times for each node.}
#'   \item{scaling_params}{List of scaling parameters (mean and sd) for numeric covariates used in the model.}
#'   \item{dummy_params}{Vector of names for dummy variables created from categorical covariates.}
#'
#' @description
#' This function aggregates the outputs from node-wise regression estimation procedures, combining them into unified results for further analysis in the GGReg model.
collect_node_results <- function(
  p,
  output_path,
  name_output,
  symm_method,
  verbose = FALSE)
{
  if (verbose) {
    cat("Collecting results from", p, "nodes...\n")
  }
  # Initialize result matrices
  first_result_file <- paste0(output_path, name_output, "_node_1.rda")
  if (!file.exists(first_result_file)) {
    stop("Result file for node 1 not found: ", first_result_file)
  }
  load(first_result_file)  # Should load 'node_result' object
  n_cols <- length(node_result$Prec_reg_matrix_row)
  Prec_reg_matrix <- matrix(0, p, n_cols)
  Sigma_hat <- rep(1, p)
  No_sim_Delta_hat <- matrix(0, p, n_cols)
  optimal_params <- data.frame(node = 1:p, asparse = NA, weight = NA, bic_score = NA)
  computational_time <- rep(0, p)
  scaling_params <- node_result$scaling_params
  dummy_params <- node_result$dummy_params
  successful_nodes <- 0
  for (i in 1:p) {
    result_file <- paste0(output_path, name_output, "_node_", i, ".rda")
    if (file.exists(result_file)) {
      load(result_file)   # Should load 'node_result' object
      if (exists("node_result") && !is.null(node_result)) {
        Prec_reg_matrix[i, ] <- node_result$Prec_reg_matrix_row
        Sigma_hat[i] <- node_result$sigma
        No_sim_Delta_hat[i, ] <- node_result$No_sim_Delta_hat_row
        
        if (!is.null(node_result$optimal_params)) {
          optimal_params[i, c("asparse", "weight", "bic_score")] <- 
            node_result$optimal_params[c("asparse", "weight", "bic_score")]
        }

        if (!is.null(node_result$computation_time)) {
          computational_time[i] <- node_result$computation_time
        }
        
        successful_nodes <- successful_nodes + 1
      } else {
        if (verbose) {
          cat("Warning: Invalid result for node", i, "\n")
        }
      }
    } else {
      if (verbose) {
        cat("Warning: Result file not found for node", i, ":", result_file, "\n")
      }
    }
  }
  if (successful_nodes == 0) {
    stop("No valid results found from any node!")
  }
  if (verbose) {
    cat("Successfully collected results from", successful_nodes, "of", p, "nodes\n")
  }
  
  # Set column and row names
  if (exists("node_result") && !is.null(node_result$col_names)) {
    colnames(Prec_reg_matrix) <- node_result$col_names
    colnames(No_sim_Delta_hat) <- node_result$col_names
  }
  if (exists("node_result") && !is.null(node_result$row_names)) {
    rownames(Prec_reg_matrix) <- node_result$row_names
    rownames(No_sim_Delta_hat) <- node_result$row_names
    names(Sigma_hat) <- node_result$row_names
  }
  
  # Build dictionary of Delta matrices (same as original)
  Dic_Delta_hat <- build_delta_dictionary(
    No_sim_Delta_hat,
    node_result$groups,
    node_result$U_names,
    symm_method,
    verbose = verbose)

  Dic_adj_matrics <- build_delta_dictionary(
    Prec_reg_matrix,
    node_result$groups,
    node_result$U_names,
    symm_method,
    verbose = verbose)

  return(list(
    Prec_reg_matrix = Prec_reg_matrix,
    No_sim_Delta_hat = No_sim_Delta_hat,
    Sigma_hat = Sigma_hat,
    Dic_Delta_hat = Dic_Delta_hat,
    Dic_adj_matrics = Dic_adj_matrics,
    optimal_params = optimal_params,
    computational_time_per_node = computational_time,
    scaling_params = scaling_params,
    dummy_params = dummy_params
  ))
}

#' build_delta_dictionary
#'
#' Constructs a dictionary of Delta matrices for the GGReg model based on provided parameters.
#'
#' @param No_sim_Delta_hat A non symmetrized matrix of coefficients derived from the node-wise regression estimation - matrix; nrow: features, ncol: features*(1+number of covariates)
#' @param groups A vector indicating the group membership of each column in No_sim_Delta_hat (e.g., which columns correspond to which covariates) - numeric vector; length: ncol(No_sim_Delta_hat)
#' @param U_names A vector of names corresponding to the groups in 'groups' (e.g., covariate names) - character vector; length: number of unique groups in 'groups'
#' @param symm_method Method for symmetrization of the estimated precision matrix - Character; "OR" or "AND"
#' @param verbose Logical indicating whether to print detailed progress messages during dictionary construction - Logical; default FALSE
#'
#' @return A list containing:
#'   \item{Dic_Delta_hat}{A List of matrices, where each matrix corresponds to a group in 'groups' and contains the symmetrized coefficients for that group. The names of the list elements correspond to 'U_names'.}
#'
#' @description
#' This function builds and returns a dictionary of Delta matrices for the GGReg model, applying the specified symmetrization method to the non-symmetrized coefficient matrix obtained from node-wise regression estimation. Each entry in the dictionary corresponds to a specific group of coefficients (e.g., related to a particular covariate) and contains the symmetrized values for that group.

build_delta_dictionary <- function(
  No_sim_Delta_hat,
  groups,
  U_names,
  symm_method,
  verbose = FALSE) 
{
  if (verbose) {
    cat(paste0("Building delta dictionary using ",symm_method," as symmetrization method \n"))
  }
  Dic_Delta_hat <- list()
  fact_groups <- as.factor(groups)
  
  for (k in 2:length(levels(fact_groups))) {
    id <- levels(fact_groups)[k]
    key <- U_names[k]
    value <- No_sim_Delta_hat[, fact_groups == id]
    if(symm_method == "OR"){
      choose_original <- abs(value) >= abs(t(value))
      sim_value <- ifelse(choose_original, value, t(value))
    } else {
      choose_original <- abs(value) < abs(t(value))
      sim_value <- ifelse(choose_original, value, t(value))
    }
    Dic_Delta_hat[[key]] <- sim_value
  }
  return(Dic_Delta_hat)
}

#################################################
## Network Visualization and Comparison Functions
#################################################

#' plot_personalized_network
#' 
#' @param network First network 
#' @param labels Optional node labels, if NULL, will use row/column names of the network or default to "Node_1", "Node_2", etc. - Character vector; length: nrow(network)
#' @param threshold Optional threshold to apply to the network for visualization (e.g., set values below threshold to zero) - Numeric; default 0
#' @param main Plot title
#' 
#' @description
#' This function visualizes a personalized network (e.g., the adj matrix for a specific subject) using a heatmap. If the 'corrplot' package is available, it uses 'corrplot' for a more polished visualization; otherwise, it falls back to a base R heatmap. The function allows for optional node labels and a custom title for the plot.
plot_personalized_network <- function(
  network,
  labels = NULL,
  threshold = 0,
  main = "Personalized Network") 
{
  # Set default labels
  if (is.null(labels)) {
    if (!is.null(rownames(network))) {
      labels <- rownames(network)
    } else if (!is.null(colnames(network))) {
      labels <- colnames(network)
    } else {
      labels <- paste0("Node_", 1:nrow(network))
    }
  }
  
  # Apply threshold if specified
  if (threshold > 0) {
    network_plot <- network
    network_plot[abs(network_plot) < threshold] <- 0
  } else {
    network_plot <- network
  }
  
  # Apply labels
  rownames(network_plot) <- labels
  colnames(network_plot) <- labels
  # Use heatmap or similar visualization
  if (requireNamespace("corrplot", quietly = TRUE)) {
    corrplot::corrplot(network, method = "color", 
                      is.corr = FALSE, 
                      tl.col = "black",
                      tl.srt = 45,
                      tl.cex = 0.8,
                      title = main,
                      mar = c(0, 0, 2, 0),
                      cl.pos = "r")
  } else {
    # Fallback to base R heatmap
    heatmap(network, main = main, 
            Rowv = NA, Colv = NA, 
            scale = "none",
            labRow = labels,
            labCol = labels,
            cexRow = 0.8,
            cexCol = 0.8,
            margins = c(10, 10),
            col = colorRampPalette(c("blue", "white", "red"))(100))
  }
}

#' compare_networks
#' 
#' @param network1 First adjacency matrix
#' @param network2 Second adjacency matrix
#' @param method Comparison method: "difference", "correlation", "frobenius"
#' 
#' @return Comparison metric
#' 
#' @description
#' This function compares two networks (adjacency matrices) using the specified method. The "difference" method computes the mean absolute difference between the two matrices, "correlation" computes the correlation of the upper triangular elements, and "frobenius" computes the Frobenius norm of the difference between the two matrices. The function checks that the input networks have the same dimensions before performing the comparison.
compare_networks <- function(
  network1,
  network2,
  method = "difference")
{
  
  if (!identical(dim(network1), dim(network2))) {
    stop("Networks must have same dimensions")
  }
  
  if (method == "difference") {
    # Mean absolute difference
    return(mean(abs(network1 - network2)))
    
  } else if (method == "correlation") {
    # Correlation of upper triangular elements
    upper_tri <- upper.tri(network1, diag = FALSE)
    return(cor(network1[upper_tri], network2[upper_tri]))
    
  } else if (method == "frobenius") {
    # Frobenius norm of difference
    return(norm(network1 - network2, "F"))
    
  } else {
    stop("Method must be 'difference', 'correlation', or 'frobenius'")
  }
}

#' plot_network_difference
#' 
#' @param network1 First network
#' @param network2 Second network  
#' @param labels Optional node labels, if NULL, will use row/column names of the network or default to "Node_1", "Node_2", etc. - Character vector; length: nrow(network1)
#' @param main Plot title
#' 
#' @description
#' This function visualizes the differences between two networks (e.g., personalized networks for two subjects) using a heatmap. The difference is computed as the element-wise difference between the two adjacency matrices, and the resulting matrix is visualized to highlight where the networks differ. If the 'corrplot' package is available, it uses 'corrplot' for a more polished visualization; otherwise, it falls back to a base R heatmap.
plot_network_difference <- function(
  network1,
  network2,
  labels = NULL,
  main = "Network Difference") 
{
  
  diff_network <- network1 - network2

  if (is.null(labels)) {
    if (!is.null(rownames(network1))) {
      labels <- rownames(network1)
    } else if (!is.null(colnames(network1))) {
      labels <- colnames(network1)
    } else {
      labels <- paste0("Node_", 1:nrow(network1))
    }
  }
  
  # Apply labels to matrix
  rownames(diff_network) <- labels
  colnames(diff_network) <- labels
  
  # Use heatmap or similar visualization
  if (requireNamespace("corrplot", quietly = TRUE)) {
    corrplot::corrplot(diff_network, 
                      method = "color", 
                      is.corr = FALSE,
                      tl.col = "black",      
                      tl.srt = 45,           
                      tl.cex = 0.8,       
                      title = main,
                      mar = c(0, 0, 2, 0),
                      addgrid.col = "gray",  
                      cl.pos = "r")         
  } else {
    # Fallback to base R heatmap
heatmap(diff_network, 
            main = main, 
            Rowv = NA,              
            Colv = NA,              
            scale = "none",
            labRow = labels,        
            labCol = labels,        
            cexRow = 0.8,           
            cexCol = 0.8,           
            margins = c(10, 10),    
            col = colorRampPalette(c("blue", "white", "red"))(100))  }
}
