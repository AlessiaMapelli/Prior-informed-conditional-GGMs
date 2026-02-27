library(igraph)
library(MASS)
library(matlib)
library(dplyr)
library(pROC)
library(ggplot2)
library(glmnet)
library(stabs)
library(matrixcalc)
library(fake)
library(reshape2)

#################################################
## NETWORK STRUCTURE GENERATION
#################################################

#' Generate adjacency matrix with specified network topology
#' @param n_nodes Number of nodes
#' @param method Network generation method ("SF" for scale-free, "ER" for Erdős-Rényi)
#' @param n_edges Target number of edges (for SF)
#' @param power Power law exponent (for SF, default 2.5)
#' @param edge_prob Edge probability (for ER)
#' @param interval Coefficient magnitude range
#' @param sign_options Sign options for coefficients
#' @return Symmetric adjacency matrix with coefficients
generate_network_structure <- function(
  n_nodes, 
  method = "SF", # or "ER" 
  n_edges = NULL, 
  power = 2.5, 
  edge_prob = NULL,
  interval = c(0.35, 0.5), 
  sign_options = c(-1, 1)
  ) 
{
  
  # Generate base graph structure
  if (method == "SF") {
    if (is.null(n_edges)) {
      stop("n_edges must be specified for scale-free networks")
    }
    # Use preferential attachment (Barabási-Albert model)
    g <- sample_pa(n = n_nodes, power = power, directed = FALSE)
    # If we have too many edges, remove some randomly
    current_edges <- gsize(g)
    if (current_edges > n_edges) {
      edges_to_remove <- sample(E(g), current_edges - n_edges)
      g <- delete_edges(g, edges_to_remove)
    }
  } else if (method == "ER") {
    if (is.null(edge_prob)) {
      stop("edge_prob must be specified for Erdős-Rényi networks")
    }
    g <- erdos.renyi.game(n_nodes, p = edge_prob, type = "gnp", directed = FALSE)
  }
  
  # Get adjacency matrix
  adj_matrix <- as.matrix(get.adjacency(g))
  
  # Generate coefficient magnitudes
  n_edges_actual <- sum(adj_matrix[upper.tri(adj_matrix)])
  if (n_edges_actual == 0) {
    warning("Generated network has no edges")
    return(matrix(0, n_nodes, n_nodes))
  }
  
  # Create coefficient matrix
  coef_matrix <- matrix(0, n_nodes, n_nodes)
  
  # Fill upper triangle with random coefficients
  upper_indices <- which(adj_matrix[upper.tri(adj_matrix)] == 1, arr.ind = FALSE)
  upper_coords <- which(upper.tri(adj_matrix), arr.ind = TRUE)
  edge_indices <- upper_coords[upper_indices, , drop = FALSE]
  
  for (i in 1:nrow(edge_indices)) {
    row_idx <- edge_indices[i, 1]
    col_idx <- edge_indices[i, 2]
    
    # Random magnitude and sign
    magnitude <- runif(1, interval[1], interval[2])
    sign <- sample(sign_options, 1)
    coef_matrix[row_idx, col_idx] <- sign * magnitude
  }
  
  # Make symmetric
  coef_matrix <- coef_matrix + t(coef_matrix)
  
  # Add small diagonal values to ensure positive definiteness
  diag(coef_matrix) <- 0.1
  
  # Check and enforce positive definiteness
  if (!is.positive.definite(coef_matrix)) {
    eigenvals <- eigen(coef_matrix, only.values = TRUE)$values
    min_eigenval <- min(eigenvals)
    if (min_eigenval <= 0) {
      diag(coef_matrix) <- diag(coef_matrix) + abs(min_eigenval) + 0.1
    }
  }
  
  return(coef_matrix)
}

#' Generate mean effect vectors for covariates
#' @param n_nodes Number of nodes
#' @param sparsity_prob Probability of non-zero effect
#' @param effect_mean Mean of non-zero effects
#' @param effect_sd Standard deviation of non-zero effects
#' @return Vector of mean effects
generate_mean_effects <- function(
  n_nodes, 
  sparsity_prob = 0.2, 
  effect_mean = 0, 
  effect_sd = 1) 
{
  
  # Generate sparsity pattern
  non_zero_mask <- rbinom(n_nodes, 1, sparsity_prob)
  
  # Generate effects
  effects <- rep(0, n_nodes)
  n_nonzero <- sum(non_zero_mask)
  
  if (n_nonzero > 0) {
    effects[non_zero_mask == 1] <- rnorm(n_nonzero, effect_mean, effect_sd)
  }
  
  return(effects)
}

#################################################
## PRIOR KNOWLEDGE GENERATION
#################################################

#' Generate prior knowledge matrix
#' @param true_adjacency True adjacency matrix (binary)
#' @param prior_type Type of prior ("perfect", "noisy", "none")
#' @param noise_params List with parameters for noisy prior
#' @return Prior weight matrix (values in [0,1])
generate_prior_knowledge <- function(
  true_adjacency, 
  prior_type = "perfect", # or "noisy", "none"
  noise_params = list(
  sensitivity = 0.8,  # Prob of detecting true edge
  specificity = 0.9,  # Prob of correctly identifying non-edge
  confidence_range = c(0.7, 1)  # Range for edge confidence
  )
  )
{
  
  p <- nrow(true_adjacency)
  
  if (prior_type == "none") {
    return(matrix(0, p, p))  # No prior knowledge
  }
  
  # Get true edge pattern (symmetric, no diagonal)
  true_edges <- (abs(true_adjacency) > 1e-8)
  diag(true_edges) <- FALSE
  
  if (prior_type == "perfect") {
    # Perfect prior: confidence 1 for true edges, 0 for non-edges
    prior_matrix <- matrix(0, p, p)
    prior_matrix[true_edges] <- 1
    return(prior_matrix)
  }
  
  if (prior_type == "noisy") {
    # Noisy prior: imperfect detection with confidence scores
    prior_matrix <- matrix(0, p, p)
    
    # For each potential edge
    for (i in 1:(p-1)) {
      for (j in (i+1):p) {
        if (true_edges[i, j]) {
          # True edge: detect with sensitivity
          if (runif(1) < noise_params$sensitivity) {
            confidence <- runif(1, noise_params$confidence_range[1], 
                              noise_params$confidence_range[2])
            prior_matrix[i, j] <- prior_matrix[j, i] <- confidence
          }
        } else {
          # Non-edge: false positive with (1 - specificity)
          if (runif(1) > noise_params$specificity) {
            # Add false positive with lower confidence
            confidence <- runif(1, 0.1, (noise_params$confidence_range[1]-0.1))
            prior_matrix[i, j] <- prior_matrix[j, i] <- confidence
          }
        }
      }
    }
    
    return(prior_matrix)
  }
}

#################################################
## DATA SIMULATION
#################################################

#' Simulate subject-level data
#' @param n_samples Number of samples
#' @param mean_effects List of mean effect vectors for each covariate
#' @param precision_matrices List of precision matrices for each covariate
#' @param covariate_values List of covariate values for each sample
#' @return List with simulated data and individual precision matrices
simulate_subject_data <- function(
  n_samples, 
  mean_effects, 
  precision_matrices, 
  covariate_values
  )
{
  
  p <- length(mean_effects[[1]])  # Number of variables
  q <- length(mean_effects) - 1   # Number of covariates (excluding intercept)
  
  # Initialize storage
  simulated_data <- matrix(0, n_samples, p)
  individual_precision_matrices <- array(0, dim = c(p, p, n_samples))
  
  for (i in 1:n_samples) {
    # Get covariate vector for this sample (including intercept)
    x_i <- c(1, unlist(covariate_values[i, ]))
    
    # Compute mean vector
    mean_i <- rep(0, p)
    if(!is.null(mean_effects)){
      for (h in 1:length(mean_effects)) {
        mean_i <- mean_i + mean_effects[[h]] * x_i[h]
      }
    }
    
    # Compute precision matrix
    precision_i <- matrix(0, p, p)
    for (h in 1:length(precision_matrices)) {
      precision_i <- precision_i + precision_matrices[[h]] * x_i[h]
    }
    
    # Ensure positive definiteness
    if (!is.positive.definite(precision_i)) {
      eigenvals <- eigen(precision_i, only.values = TRUE)$values
      min_eigenval <- min(eigenvals)
      if (min_eigenval <= 0) {
        diag(precision_i) <- diag(precision_i) + abs(min_eigenval) + 0.1
      }
    }
    
    # Compute covariance matrix
    covariance_i <- solve(precision_i)
    
    # Simulate data
    simulated_data[i, ] <- mvrnorm(1, mean_i, covariance_i)
    individual_precision_matrices[, , i] <- precision_i
  }
  
  return(list(
    data = simulated_data,
    precision_matrices = individual_precision_matrices
  ))
}

#################################################
## EVALUATION METRICS
#################################################

#' Calculate evaluation metrics for binary classification
#' @param true_matrix True binary adjacency matrix
#' @param estimated_matrix Estimated binary adjacency matrix
#' @return List of metrics (TPR, FPR, F1, etc.)
calculate_binary_metrics <- function(
  true_matrix,
  estimated_matrix
  )
{
  # Ensure binary and symmetric, no diagonal
  true_binary <- (abs(true_matrix) > 1e-8)
  est_binary <- (abs(estimated_matrix) > 1e-8)
  diag(true_binary) <- diag(est_binary) <- FALSE
  
  # Get upper triangle (since symmetric)
  true_edges <- true_binary[upper.tri(true_binary)]
  est_edges <- est_binary[upper.tri(est_binary)]
  
  # Confusion matrix elements
  TP <- sum(true_edges & est_edges)
  FP <- sum(!true_edges & est_edges)
  TN <- sum(!true_edges & !est_edges)
  FN <- sum(true_edges & !est_edges)
  
  # Calculate metrics
  TPR <- ifelse(TP + FN > 0, TP / (TP + FN), 0)  # Sensitivity/Recall
  FPR <- ifelse(FP + TN > 0, FP / (FP + TN), 0)  # 1 - Specificity
  Precision <- ifelse(TP + FP > 0, TP / (TP + FP), 0)
  F1 <- ifelse(Precision + TPR > 0, 2 * Precision * TPR / (Precision + TPR), 0)
  Accuracy <- (TP + TN) / (TP + TN + FP + FN)
  
  return(list(
    TP = TP, FP = FP, TN = TN, FN = FN,
    TPR = TPR, FPR = FPR, Precision = Precision,
    F1 = F1, Accuracy = Accuracy
  ))
}

#' Calculate evaluation metrics for binary classification (for vectors)
#' @param true_vector True binary vector
#' @param estimated_vector Estimated binary vector
#' @return List of metrics (TPR, FPR, F1, etc.)
calculate_binary_metrics_vector <- function(
  true_vector,
   estimated_vector
   )
{
  # Ensure binary
  true_binary <- (abs(true_vector) > 1e-8)
  est_binary <- (abs(estimated_vector) > 1e-8)
  
  # Confusion matrix elements
  TP <- sum(true_binary & est_binary)
  FP <- sum(!true_binary & est_binary)
  TN <- sum(!true_binary & !est_binary)
  FN <- sum(true_binary & !est_binary)
  
  # Calculate metrics
  TPR <- ifelse(TP + FN > 0, TP / (TP + FN), 0)  # Sensitivity/Recall
  FPR <- ifelse(FP + TN > 0, FP / (FP + TN), 0)  # 1 - Specificity
  Precision <- ifelse(TP + FP > 0, TP / (TP + FP), 0)
  F1 <- ifelse(Precision + TPR > 0, 2 * Precision * TPR / (Precision + TPR), 0)
  Accuracy <- (TP + TN) / (TP + TN + FP + FN)
  
  return(list(
    TP = TP, FP = FP, TN = TN, FN = FN,
    TPR = TPR, FPR = FPR, Precision = Precision,
    F1 = F1, Accuracy = Accuracy
  ))
}

#' Calculate Frobenius norm error
#' @param true_matrix True matrix
#' @param estimated_matrix Estimated matrix
#' @return Frobenius norm of difference
calculate_frobenius_error <- function(
  true_matrix,
  estimated_matrix
  ) 
{
  return(norm(true_matrix - estimated_matrix, type = "F"))
}

#' Calculate vector norm error (L2 norm)
#' @param true_vector True vector
#' @param estimated_vector Estimated vector
#' @return L2 norm of difference
calculate_vector_error <- function(
  true_vector,
  estimated_vector
  ) 
{
  return(sqrt(sum((true_vector - estimated_vector)^2)))
}

#' Evaluate mean effects estimation for each covariate
#' @param true_mean_effects List of true mean effect vectors for each covariate
#' @param estimated_mean_effects Matrix of estimated mean effects (p x q+1)
#' @return List of evaluation metrics for mean effects
evaluate_mean_effects <- function(
  true_mean_effects,
   estimated_mean_effects
   ) 
{
  
  mean_results <- list()
  
  # Get covariate names from true effects
  covariate_names <- names(true_mean_effects)
  
  # Convert estimated mean effects to list format if it's a matrix
  if (is.matrix(estimated_mean_effects)) {
    est_mean_list <- list()
    for (i in 1:ncol(estimated_mean_effects)) {
      cov_name <- colnames(estimated_mean_effects)[i]
      if (is.null(cov_name)) {
        cov_name <- covariate_names[i]  # Use true names if estimated names missing
      }
      est_mean_list[[cov_name]] <- estimated_mean_effects[, i]
    }
    estimated_mean_effects <- est_mean_list
  }
  
  # Evaluate each covariate's mean effects
  for (name in covariate_names) {
    if (name %in% names(estimated_mean_effects)) {
      true_vec <- true_mean_effects[[name]]
      est_vec <- estimated_mean_effects[[name]]
      
      # Ensure same length
      if (length(true_vec) != length(est_vec)) {
        warning(paste("Length mismatch for covariate", name, ": true =", 
                     length(true_vec), ", estimated =", length(est_vec)))
        next
      }
      
      # Binary metrics for sparsity pattern detection
      binary_metrics <- calculate_binary_metrics_vector(true_vec, est_vec)
      
      # Estimation error
      vector_error <- calculate_vector_error(true_vec, est_vec)
      
      # Relative error (normalized by true norm)
      true_norm <- sqrt(sum(true_vec^2))
      relative_error <- ifelse(true_norm > 1e-8, vector_error / true_norm, vector_error)
      
      mean_results[[paste0("Mean_", name)]] <- c(
        binary_metrics,
        list(
          Vector_Error = vector_error,
          Relative_Error = relative_error,
          True_Norm = true_norm
        )
      )
    }
  }
  
  # Overall mean effects metrics (average across covariates, excluding intercept for some metrics)
  if (length(mean_results) > 0) {
    all_results <- 
    metric_names <- c("TPR", "FPR", "F1", "Accuracy", "Vector_Error", "Relative_Error")
    overall_mean_metrics <- list()
        
    for (metric in metric_names) {
      # Overall metrics across all covariates
      all_values <- sapply(mean_results[-1], function(x) x[[metric]])
      overall_mean_metrics[[paste0("Overall_Mean_", metric)]] <- mean(all_values, na.rm = TRUE)
      
    }
    
    mean_results <- c(mean_results, overall_mean_metrics)
  }
  
  return(mean_results)
}

#' Evaluate Magnitude Preservation in Estimated Precision Matrix
#' @param true_matrix True binary adjacency matrix
#' @param estimated_matrix Estimated binary adjacency matrix
#' @return spearman correlation of the two
calculate_magnitude_preserved <- function(
  true_matrix, 
  estimated_matrix
  ) 
{

  # Get upper triangular elements (excluding diagonal)
  upper_tri_idx <- upper.tri(true_matrix, diag = FALSE)
  true_values <- true_matrix[upper_tri_idx]
  est_values <- estimated_matrix[upper_tri_idx]

  # Remove any NA or infinite values
  valid_idx <- is.finite(true_values) & is.finite(est_values)
  true_values <- true_values[valid_idx]
  est_values <- est_values[valid_idx]
  
  if(sum(true_values)==0 |sum(est_values)==0){
    magnitude_preserved =0
  }else{magnitude_preserved = cor(true_values, est_values, method = "spearman")}
  
  
  return(magnitude_preserved)
    
}

#' Evaluate estimation results
#' @param true_Delta_list List of true precision coefficient matrices
#' @param estimated_Delta_list List of estimated precision coefficient matrices
#' @param true_Gamma True mean coefficient matrix
#' @param estimated_Gamma Estimated mean coefficient matrix
#' @return List of evaluation metrics
evaluate_estimation <- function(
  true_Delta_list, 
  estimated_Delta_list, 
  true_mean_effects = NULL, 
  estimated_mean_effects = NULL
  ) 
{  
  results <- list()
  
  # Evaluate Delta matrices (precision coefficients)
  delta_names <- names(true_Delta_list)
  
  for (name in delta_names) {
    if (name %in% names(estimated_Delta_list)) {
      true_mat <- true_Delta_list[[name]]
      est_mat <- estimated_Delta_list[[name]]
      
      # Binary metrics
      binary_metrics <- calculate_binary_metrics(true_mat, est_mat)
      
      # Estimation error
      frobenius_error <- calculate_frobenius_error(true_mat, est_mat)
      magnitude_preserved <- calculate_magnitude_preserved(true_mat, est_mat)
      
      results[[paste0("Delta_", name)]] <- c(binary_metrics, 
                                           list(Frobenius_Error = frobenius_error, magnitude_preserved= magnitude_preserved ))
    }
  }
  
  # Overall Delta metrics (average across components)
  if (length(results) > 0) {
    metric_names <- c("TPR", "FPR", "F1", "Accuracy", "Frobenius_Error", "magnitude_preserved")
    overall_metrics <- list()
    
    for (metric in metric_names) {
      values <- sapply(results, function(x) x[[metric]])
      overall_metrics[[paste0("Overall_Delta_", metric)]] <- mean(values, na.rm = TRUE)
    }
    
    results <- c(results, overall_metrics)
  }
  
  if (!is.null(true_mean_effects) && !is.null(estimated_mean_effects)) {
    mean_evaluation <- evaluate_mean_effects(true_mean_effects, estimated_mean_effects)
    results <- c(results, mean_evaluation)
  }
  
  return(results)
}

#################################################
## SIMULATION UTILITIES
#################################################

#' Generate complete simulation dataset
#' @param n_samples Number of samples
#' @param n_nodes Number of nodes (proteins)
#' @param covariate_config List specifying covariate configuration
#' @param network_config List specifying network generation parameters
#' @param prior_config List specifying prior knowledge parameters
#' @param seed Random seed for reproducibility
#' @return Complete simulation dataset and metadata
generate_simulation_dataset <- function(
  n_samples, 
  n_nodes, 
  covariate_config,
  network_config,
  prior_config,
  seed = NULL, 
  output_path = "/group/diangelantonio/users/alessia_mapelli/Prot_graphs/UKB_data/APP_82779/Simulation_results",
  name_output = "ggReg_result"
  )
{
  
  if (!is.null(seed)) set.seed(seed)
  
  # Validate inputs
  if (n_samples < 10) stop("n_samples must be at least 10")
  if (n_nodes < 3) stop("n_nodes must be at least 3")
  
  # Extract configurations
  q <- covariate_config$n_continuous + covariate_config$n_binary
  
  # Generate covariate data
  covariate_data <- data.frame(row.names = 1:n_samples)
  
  # Add continuous covariates
  if (covariate_config$n_continuous > 0) {
    for (i in 1:covariate_config$n_continuous) {
      covariate_data[[paste0("x", i)]] <- runif(n_samples, 0, 1)
    }
  }
  
  # Add binary covariates
  if (covariate_config$n_binary > 0) {
    start_idx <- covariate_config$n_continuous + 1
    end_idx <- start_idx + covariate_config$n_binary - 1
    
    for (i in start_idx:end_idx) {
      covariate_data[[paste0("x", i)]] <- rbinom(n_samples, 1, 0.5)
    }
  }
  
  if(effective_mean){
      # Generate mean effects for each covariate component
     mean_effects <- list()
  
    # Intercept
    mean_effects[["(Intercept)"]] <- generate_mean_effects(
      n_nodes, 
      sparsity_prob = 0,  # All nodes have intercept
      effect_mean = 0, 
      effect_sd = 1
    )

    # Covariate effects
    for (i in 1:ncol(covariate_data)) {
      cov_name <- colnames(covariate_data)[i]
        mean_effects[[cov_name]] <- generate_mean_effects(
          n_nodes,
          sparsity_prob = covariate_config$mean_sparsity,
          effect_mean = 0.25,
          effect_sd = 0.001
        )
    }
  }else{
    mean_effects = NULL
  }
  
  # Generate precision coefficient matrices
  precision_matrices <- list()
  
  # Population-level precision matrix (Δ₀)
  precision_matrices[["Prot"]] <- generate_network_structure(
    n_nodes = n_nodes,
    method = network_config$population_method,
    n_edges = network_config$population_edges,
    power = network_config$population_power,
    edge_prob = network_config$population_prob,
    interval = network_config$coefficient_range,
    sign_options = network_config$sign_options
  )
  
  # Ensure positive definiteness with stronger diagonal
  min_diag <- max(0.5, max(abs(precision_matrices[["Prot"]][upper.tri(precision_matrices[["Prot"]])])) * 1.2)
  diag(precision_matrices[["Prot"]]) <- min_diag
  
  adj_df <- melt(abs(precision_matrices[["Prot"]]))
  colnames(adj_df) <- c("Row", "Col", "Value")
  plot <- ggplot(adj_df, aes(x = Col, y = Row, fill = Value)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "black") +
    theme_minimal() +
    coord_fixed() +
    scale_y_reverse() +  # to match matrix view
    labs(title = "Node adjacency Matrix", x = "", y = "")
  ggsave(plot, filename=paste(output_path, "Adjacency_matrix_",name_output ,"_", "Prot", ".png", sep=""), width=8, height=8, dpi=300)
  
  index_effective <- sample(1:ncol(covariate_data), covariate_config$effective_cov)

  # Covariate-dependent precision matrices
  for (i in 1:ncol(covariate_data)) {
    cov_name <- colnames(covariate_data)[i]

    if(i %in% index_effective){
      # Generate smaller covariate effects to maintain stability
    reduced_range <- network_config$coefficient_range * 0.5  # Reduce effect size
    
    precision_matrices[[cov_name]] <- generate_network_structure(
      n_nodes = n_nodes,
      method = network_config$covariate_method,
      n_edges = network_config$covariate_edges,
      power = network_config$power,
      edge_prob = network_config$covariate_prob,  # Sparser covariate effects
      interval = reduced_range,
      sign_options = network_config$sign_options
    )
    }else{
      precision_matrices[[cov_name]] <- matrix(0,n_nodes, n_nodes)
    }
    
    
    adj_df <- melt(abs(precision_matrices[[cov_name]]))
    colnames(adj_df) <- c("Row", "Col", "Value")
    plot <- ggplot(adj_df, aes(x = Col, y = Row, fill = Value)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "black") +
      theme_minimal() +
      coord_fixed() +
      scale_y_reverse() +  # to match matrix view
      labs(title = "Node adjacency Matrix", x = "", y = "")
    ggsave(plot, filename=paste(output_path, "Adjacency_matrix_",name_output ,"_", cov_name, ".png", sep=""), width=8, height=8, dpi=300)
    
    
  }
  
  # Simulate molecular data with error handling
  sim_result <- simulate_subject_data_robust(
    n_samples = n_samples,
    mean_effects = mean_effects,
    precision_matrices = precision_matrices,
    covariate_values = covariate_data
  )
  
  # Create molecular data frame
  molecular_data <- as.data.frame(sim_result$data)
  colnames(molecular_data) <- paste0("Prot", 1:n_nodes)
  
  # Scale molecular data to reasonable range
  molecular_data <- as.data.frame(scale(molecular_data))
  
  # Combine all data
  complete_data <- cbind(
    data.frame(pat_id = 1:n_samples),
    covariate_data,
    molecular_data
  )
  
  # Generate prior knowledge
  prior_knowledge <- generate_prior_knowledge(
    true_adjacency = precision_matrices[["Prot"]],
    prior_type = prior_config$type,
    noise_params = prior_config$noise_params
  )
  
  adj_df <- melt(prior_knowledge)
  colnames(adj_df) <- c("Row", "Col", "Value")
  plot <- ggplot(adj_df, aes(x = Col, y = Row, fill = Value)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "black") +
    theme_minimal() +
    coord_fixed() +
    scale_y_reverse() +  # to match matrix view
    labs(title = "Node adjacency Matrix", x = "", y = "")
  ggsave(plot, filename=paste(output_path, "Prior_knowledge_matrix_",name_output, ".png", sep=""), width=8, height=8, dpi=300)
  
  if(sum(prior_knowledge > 0) == 0){
    prior_knowledge = NULL
  }
  
  return(list(
    data = complete_data,
    true_mean_effects = mean_effects,
    true_precision_matrices = precision_matrices,
    prior_knowledge = prior_knowledge,
    individual_precision_matrices = sim_result$precision_matrices,
    config = list(
      n_samples = n_samples,
      n_nodes = n_nodes,
      covariate_config = covariate_config,
      network_config = network_config,
      prior_config = prior_config
    )
  ))
}

#' Robust simulation of subject-level data with error handling for precision matrix positive definiteness and numerical stability
#' @param n_samples Number of samples
#' @param mean_effects List of mean effect vectors for each covariate
#' @param precision_matrices List of precision matrices for each covariate
#' @param covariate_values List of covariate values for each sample
#' @return List with simulated data and individual precision matrices
simulate_subject_data_robust <- function(
  n_samples, 
  mean_effects, 
  precision_matrices, 
  covariate_values) 
{
  
  p <- length(mean_effects[[1]])  # Number of variables
  q <- length(mean_effects) - 1   # Number of covariates (excluding intercept)
  
  # Initialize storage
  simulated_data <- matrix(0, n_samples, p)
  individual_precision_matrices <- array(0, dim = c(p, p, n_samples))
  
  for (i in 1:n_samples) {
    # Get covariate vector for this sample (including intercept)
    x_i <- c(1, unlist(covariate_values[i, ]))
    
    # Compute mean vector
    mean_i <- rep(0, p)
    if(!is.null(mean_effects)){
      for (h in 1:length(mean_effects)) {
        mean_i <- mean_i + mean_effects[[h]] * x_i[h]
      }
     }
    
    # Compute precision matrix
    precision_i <- matrix(0, p, p)
    for (h in 1:length(precision_matrices)) {
      precision_i <- precision_i + precision_matrices[[h]] * x_i[h]
    }
    
    # Robust positive definiteness check
    max_attempts <- 3
    attempt <- 1
    
    while (attempt <= max_attempts) {
      if (!is.positive.definite(precision_i)) {
        eigenvals <- eigen(precision_i, symmetric = TRUE)
        min_eigenval <- min(eigenvals$values)
        if (min_eigenval <= 0) {
          correction <- abs(min_eigenval) + 0.1 * attempt
          diag(precision_i) <- diag(precision_i) + correction
        }
        attempt <- attempt + 1
      } else {
        break
      }
    }
    
    # Final check - if still not positive definite, use identity
    if (!is.positive.definite(precision_i)) {
      warning(paste("Using regularized precision matrix for sample", i))
      precision_i <- diag(p) * max(diag(precision_i))
    }
    
    # Compute covariance matrix with error handling
    tryCatch({
      covariance_i <- solve(precision_i)
    }, error = function(e) {
      # Use pseudo-inverse if solve fails
      covariance_i <- MASS::ginv(precision_i)
      warning(paste("Used pseudo-inverse for sample", i))
    })
    
    # Simulate data with error handling
    tryCatch({
      simulated_data[i, ] <- mvrnorm(1, mean_i, covariance_i)
    }, error = function(e) {
      # Fallback to independent normal if mvrnorm fails
      simulated_data[i, ] <- rnorm(p, mean_i, sqrt(diag(covariance_i)))
      warning(paste("Used independent sampling for sample", i))
    })
    
    individual_precision_matrices[, , i] <- precision_i
  }
  
  return(list(
    data = simulated_data,
    precision_matrices = individual_precision_matrices
  ))
}


#' Validate prediction performance of proteins values on covariates using LASSO and stability selection
#' @param sim_data Simulated dataset
#' @param true_precision_coef True precision coefficient matrix for outcome
#' @param outcome_var Name of outcome variable
#' @param feature_type Type of features to use ("proteins", "interactions", "both")
#' @return List with prediction performance AUC metrics for LASSO and stability selection
validate_data_prediction_performance <- function(
  sim_data, 
  true_precision_coef, 
  outcome_var = "x1",
  feature_type = "proteins" #or "interactions", "both"
  ) 
{
  
  # Get protein column names
  protein_cols <- colnames(sim_data)[grepl("^Prot", colnames(sim_data))]
  
  # Create interaction terms if needed
  if (feature_type %in% c("interactions", "both")) {
    true_adj <- (abs(true_precision_coef) > 1e-8)
    diag(true_adj) <- FALSE
    
    if (sum(true_adj) > 0) {
      g <- graph_from_adjacency_matrix(true_adj, mode = "undirected")
      edge_list <- as_edgelist(g)
      
      interaction_terms <- apply(edge_list, 1, function(edge) {
        paste0("Prot", edge[1], "*Prot", edge[2])
      })
    } else {
      interaction_terms <- character(0)
    }
  }
  
  # Build formula
  if (feature_type == "proteins") {
    predictor_terms <- protein_cols
  } else if (feature_type == "interactions") {
    predictor_terms <- interaction_terms
  } else {  # both
    predictor_terms <- c(protein_cols, interaction_terms)
  }
  
  if (length(predictor_terms) == 0) {
    return(list(auc_lasso = NA, auc_stable = NA))
  }
  
  formula_string <- paste(outcome_var, "~", paste(predictor_terms, collapse = " + "))
  
  # Create design matrix
  X <- model.matrix(as.formula(formula_string), sim_data)[, -1]  # Remove intercept
  y <- sim_data[[outcome_var]]
  
  # Determine family
  family_type <- if (length(unique(y)) == 2) "binomial" else "gaussian"
  
  # LASSO regression
  set.seed(123)
  cv_lasso <- cv.glmnet(X, y, alpha = 1, nfolds = 3, family = family_type)
  
  y_pred_lasso <- predict(cv_lasso, newx = X, s = "lambda.min", type = "response")
  
  if (family_type == "binomial") {
    roc_lasso <- roc(response = as.factor(y), predictor = as.numeric(y_pred_lasso),
                     levels = c(0, 1), smooth = FALSE, quiet = TRUE)
    auc_lasso <- auc(roc_lasso)
  } else {
    auc_lasso <- cor(y, y_pred_lasso)^2  # R-squared for continuous outcomes
  }
  
  # Stability selection (if package available and binary outcome)
  auc_stable <- NA
  if (family_type == "binomial" && requireNamespace("stabs", quietly = TRUE)) {
    tryCatch({
      stab_result <- stabs::stabsel(x = X, y = y, fitfun = stabs::glmnet.lasso, 
                                   cutoff = 0.9, PFER = 4)
      
      if (length(stab_result$selected) > 0) {
        selected_vars <- names(stab_result$selected)
        formula_stable <- paste("y ~", paste(selected_vars, collapse = " + "))
        
        temp_data <- as.data.frame(cbind(y = y, X[, selected_vars, drop = FALSE]))
        mod_glm <- glm(as.formula(formula_stable), family = binomial(link = "logit"), 
                      data = temp_data)
        
        y_pred_stable <- predict(mod_glm, type = "response")
        roc_stable <- roc(response = as.factor(y), predictor = y_pred_stable,
                         levels = c(0, 1), smooth = FALSE, quiet = TRUE)
        auc_stable <- auc(roc_stable)
      }
    }, error = function(e) {
      auc_stable <- NA
    })
  }
  
  return(list(auc_lasso = auc_lasso, auc_stable = auc_stable))
}



#################################################
## SIMULATION EXECUTION FUNCTIONS
#################################################

#' Run a single simulation replicate
#' @param config Single row from simulation grid
#' @param rep_id Replication ID
#' @param seed Random seed for this replicate
#' @return List of results for this replicate
run_single_simulation <- function(
  config,
  rep_id,
  seed,
  output_path = "/group/diangelantonio/users/alessia_mapelli/Prot_graphs/UKB_data/APP_82779/Simulation_results",
  name_output = "ggReg_result"
  ) 
{
  
  cat(sprintf("Running: n=%d, p=%d, q=%d, prior=%s, symm=%s, rep=%d\n",
              config$n_samples, config$n_nodes, config$n_covariates, 
              config$prior_type, config$symm_method, rep_id))
  
  start_time <- Sys.time()
  
  # Set up configurations
  network_config <- NETWORK_CONFIG
  network_config$population_edges <- round(config$n_nodes * (config$n_nodes - 1) / 2 * 
                                              network_config$population_density)
  
  covariate_config <- base_covariate_config(config$n_covariates)
  
  prior_config <- list(
    type = config$prior_type,
    noise_params = PRIOR_NOISE_PARAMS
  )
  
  method_params <- get_method_params(config$prior_type)
  
  # Generate simulation dataset
  sim_dataset <- generate_simulation_dataset(
    n_samples = config$n_samples,
    n_nodes = config$n_nodes,
    covariate_config = covariate_config,
    network_config = network_config,
    prior_config = prior_config,
    seed = seed,
    output_path = output_path,
    name_output = name_output

  )
  
  # Extract data components
  sim_data <- sim_dataset$data
  protein_cols <- colnames(sim_data)[grepl("^Prot", colnames(sim_data))]
  covariate_cols <- colnames(sim_data)[grepl("^x", colnames(sim_data))]
  
  # Prepare covariates dataframe
  if (length(covariate_cols) > 0) {
    covariates_df <- sim_data[, covariate_cols, drop = FALSE]
  } else {
    covariates_df <- NULL
  }
  
  if(is.null(sim_data$true_mean_effects)){
    mean_estimation=FALSE
  }else{mean_estimation=TRUE}
  
  # Run GGReg estimation
  ggReg_results <- GGReg_full_estimation(
    x = sim_data[, protein_cols],
    known_ppi = sim_dataset$prior_knowledge,
    covariates = covariates_df,
    scr = method_params$screening_procedure,
    mean_estimation = mean_estimation,
    lambda_prec_type = method_params$lambda_prec_type,
    tune_hyperparams = method_params$tune_hyperparams,
    asparse_grid = method_params$asparse_grid,
    weight_grid = method_params$weight_grid,
    random_hyper_search = method_params$random_hyper_search,
    p.rand.hyper = method_params$p.rand.hyper,
    K = method_params$K,
    use_slurm = FALSE,
    slurm_script_path = method_params$slurm_script_path,
    output_path = output_path,
    name_output = name_output,
    symm_method =config$symm_method,
    verbose = method_params$verbose
  )

  list_estimated_prec <- ggReg_results$additional_info$Dic_Delta_hat
  estimation_matrix_names <- names(list_estimated_prec)
  for(l in 1:length(estimation_matrix_names)){
    name_cov <- estimation_matrix_names[l]
    adj_df <- melt(abs(list_estimated_prec[[name_cov]]))
    colnames(adj_df) <- c("Row", "Col", "Value")
    adj_df$Value <- as.numeric(adj_df$Value)
    adj_df$Row <- as.numeric(sub("Prot", "", adj_df$Row))
    if(l==1){
      adj_df$Col <- as.numeric(sub("Prot", "", adj_df$Col))
    } else {
      adj_df$Col <- as.numeric(sub(paste(name_cov,":","Prot", sep=""), "", adj_df$Col))
    }
    # Plot
    plot <- ggplot(adj_df, aes(x = Col, y = Row, fill = Value)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "black") +
      theme_minimal() +
      coord_fixed() +
      scale_y_reverse() +  # to match matrix view
      labs(title = "Estimated adjacency Matrix", x = "", y = "")
    ggsave(plot, filename=paste(output_path, "Estimated_matrix_",name_output ,"_", name_cov, ".png", sep=""), width=8, height=8, dpi=300)

  }
  
  # Evaluate results
  if(is.null(sim_data$true_mean_effects)){
    evaluation_results <- evaluate_estimation(
      true_Delta_list = sim_dataset$true_precision_matrices,
      estimated_Delta_list = list_estimated_prec,
      true_mean_effects = NULL,
      estimated_mean_effects = NULL
    )
  }else{
    evaluation_results <- evaluate_estimation(
      true_Delta_list = sim_dataset$true_precision_matrices,
      estimated_Delta_list = list_estimated_prec,
      true_mean_effects = sim_dataset$true_mean_effects,
      estimated_mean_effects = ggReg_results$results$Cov_effect
    )
  }
  
  
  computation_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  # Compile results
  result <- list(
    config = config,
    rep_id = rep_id,
    evaluation = evaluation_results,
    computation_time = computation_time
  )
  
  return(result)
}

#' Run complete simulation study
#' @param output_file File to save results
#' @return List of all simulation results
run_complete_simulation <- function(
  SIMULATION_GRID = TEST_SIMULATION_GRID,
  N_REPLICATIONS = N_REPLICATIONS,
  output_folder= "/group/diangelantonio/users/alessia_mapelli/Prot_graphs/UKB_data/APP_82779/Simulation_results",
  output_file = "simulation_results.RData")
{
  
  cat("Starting comprehensive simulation study...\n")
  cat(sprintf("Total configurations: %d\n", nrow(SIMULATION_GRID)))
  cat(sprintf("Replications per configuration: %d\n", N_REPLICATIONS))
  cat(sprintf("Total simulations: %d\n", nrow(SIMULATION_GRID) * N_REPLICATIONS))
  
  all_results <- list()
  total_sims <- nrow(SIMULATION_GRID) * N_REPLICATIONS
  current_sim <- 0
  
  # Run simulations
  for (i in 1:nrow(SIMULATION_GRID)) {
    config <- SIMULATION_GRID[i, ]

    ggReg_output_folder <- sprintf("n%d_p%d_q%d_%s_%s", 
                                    config$n_samples, config$n_nodes, config$n_covariates,
                                    config$prior_type, config$symm_method)
    ggReg_output_path <- paste(output_folder,"/", ggReg_output_folder,"/", sep ="")
    dir.create(ggReg_output_path)
    
    for (rep in 1:N_REPLICATIONS) {
      current_sim <- current_sim + 1
      
      # Progress indicator
      if (current_sim %% 10 == 0) {
        cat(sprintf("Progress: %d/%d (%.1f%%)\n", 
                    current_sim, total_sims, 100 * current_sim / total_sims))
      }
      
      # Set unique seed for reproducibility
      seed <- 12345 + i * 1000 + rep
      
      ggReg_output_name <- sprintf("rep%d", rep)
      ggReg_output_path_rep <- paste(ggReg_output_path, ggReg_output_name, "/", sep="")
      dir.create(ggReg_output_path_rep)
      
      # Run simulation
      result <- run_single_simulation(config, rep, seed, ggReg_output_path_rep, ggReg_output_name)
      
      # Store result
      result_id <- sprintf("n%d_p%d_q%d_%s_%s_rep%d", 
                           config$n_samples, config$n_nodes, config$n_covariates,
                           config$prior_type, config$symm_method, rep)
      all_results[[result_id]] <- result
      
      # Save intermediate results every 50 simulations
      if (current_sim %% N_REPLICATIONS == 0) {
        cat("Saving intermediate results...\n")
        save(all_results, SIMULATION_GRID, file = paste0(output_folder, "/",ggReg_output_folder,"/", "setting_", i, "_", output_file))
      }
    }
  }
  
  # Save final results
  cat("Saving final results...\n")
  save(all_results, SIMULATION_GRID, file = paste(output_folder,ggReg_output_folder,output_file, sep="/"))
  
  return(all_results)
}

#################################################
## RESULTS ANALYSIS AND VISUALIZATION
#################################################

#' Process simulation results into analysis-ready format
#' @param results List of simulation results
#' @return Data frame with processed results
process_simulation_results <- function(
  results
  ) 
{
  
  cat("Processing simulation results...\n")
  
  processed_data <- data.frame()
  
  for (result_id in names(results)) {
    result <- results[[result_id]]
    
    if (is.null(result$evaluation)) {
      next  # Skip failed simulations
    }
    
    config <- result$config
    evaluation <- result$evaluation
    
    # Extract metrics for each Delta component
    delta_components <- grep("^Delta_", names(evaluation), value = TRUE)
    
    for (component in delta_components) {
      comp_name <- gsub("^Delta_", "", component)
      comp_metrics <- evaluation[[component]]
      
      if (is.list(comp_metrics)) {
        row_data <- data.frame(
          n_samples = config$n_samples,
          n_nodes = config$n_nodes,
          n_covariates = config$n_covariates,
          prior_type = config$prior_type,
          symm_method = config$symm_method,
          rep_id = result$rep_id,
          component = comp_name,
          component_type = "delta_individual",
          TP = comp_metrics$TP %||% NA,
          FP = comp_metrics$FP %||% NA,
          TN = comp_metrics$TN %||% NA,
          FN = comp_metrics$FN %||% NA,
          TPR = comp_metrics$TPR %||% NA,
          FPR = comp_metrics$FPR %||% NA,
          F1 = comp_metrics$F1 %||% NA,
          Accuracy = comp_metrics$Accuracy %||% NA,
          Frobenius_Error = comp_metrics$Frobenius_Error %||% NA,
          Magnitude_preserved = comp_metrics$magnitude_preserved %||% NA,
          Vector_Error = NA,  # Not applicable for matrices
          Relative_Error = NA,  # Not applicable for matrices
          computation_time = result$computation_time,
          stringsAsFactors = FALSE
        )
        
        processed_data <- rbind(processed_data, row_data)
      }
    }
    
    # Extract metrics for each Mean component
    mean_components <- grep("^Mean_", names(evaluation), value = TRUE)
    
    for (component in mean_components) {
      comp_name <- gsub("^Mean_", "", component)
      comp_metrics <- evaluation[[component]]
      
      if (is.list(comp_metrics)) {
        row_data <- data.frame(
          n_samples = config$n_samples,
          n_nodes = config$n_nodes,
          n_covariates = config$n_covariates,
          prior_type = config$prior_type,
          symm_method = config$symm_method,
          rep_id = result$rep_id,
          component = comp_name,
          component_type = "mean_individual",
          TP = comp_metrics$TP %||% NA,
          FP = comp_metrics$FP %||% NA,
          TN = comp_metrics$TN %||% NA,
          FN = comp_metrics$FN %||% NA,
          TPR = comp_metrics$TPR %||% NA,
          FPR = comp_metrics$FPR %||% NA,
          F1 = comp_metrics$F1 %||% NA,
          Accuracy = comp_metrics$Accuracy %||% NA,
          Frobenius_Error = NA,  # Not applicable for vectors
          Vector_Error = comp_metrics$Vector_Error %||% NA,
          Relative_Error = comp_metrics$Relative_Error %||% NA,
          computation_time = result$computation_time,
          stringsAsFactors = FALSE
        )
        
        processed_data <- rbind(processed_data, row_data)
      }
    }
    
    # Extract overall Delta metrics
    overall_delta_metrics <- evaluation[grep("^Overall_Delta_", names(evaluation))]
    
    if (length(overall_delta_metrics) > 0) {
      row_data <- data.frame(
        n_samples = config$n_samples,
        n_nodes = config$n_nodes,
        n_covariates = config$n_covariates,
        prior_type = config$prior_type,
        symm_method = config$symm_method,
        rep_id = result$rep_id,
        component = "Overall_Delta",
        component_type = "delta_overall",
        TP = NA,
        FP =  NA,
        TN =  NA,
        FN =  NA,
        TPR = overall_delta_metrics$Overall_Delta_TPR %||% NA,
        FPR = overall_delta_metrics$Overall_Delta_FPR %||% NA,
        F1 = overall_delta_metrics$Overall_Delta_F1 %||% NA,
        Accuracy = overall_delta_metrics$Overall_Delta_Accuracy %||% NA,
        Frobenius_Error = overall_delta_metrics$Overall_Delta_Frobenius_Error %||% NA,
        Magnitude_preserved = overall_delta_metrics$magnitude_preserved %||% NA,
        Vector_Error = NA,
        Relative_Error = NA,
        computation_time = result$computation_time,
        stringsAsFactors = FALSE
      )
      
      processed_data <- rbind(processed_data, row_data)
    }
    
    # Extract overall Mean metrics
    overall_mean_metrics <- evaluation[grep("^Overall_Mean_", names(evaluation))]
    
    if (length(overall_mean_metrics) > 0) {
      row_data <- data.frame(
        n_samples = config$n_samples,
        n_nodes = config$n_nodes,
        n_covariates = config$n_covariates,
        prior_type = config$prior_type,
        symm_method = config$symm_method,
        rep_id = result$rep_id,
        component = "Overall_Mean",
        component_type = "mean_overall",
        TP = NA,
        FP =  NA,
        TN =  NA,
        FN =  NA,
        TPR = overall_mean_metrics$Overall_Mean_TPR %||% NA,
        FPR = overall_mean_metrics$Overall_Mean_FPR %||% NA,
        F1 = overall_mean_metrics$Overall_Mean_F1 %||% NA,
        Accuracy = overall_mean_metrics$Overall_Mean_Accuracy %||% NA,
        Frobenius_Error = NA,
        Vector_Error = overall_mean_metrics$Overall_Mean_Vector_Error %||% NA,
        Relative_Error = overall_mean_metrics$Overall_Mean_Relative_Error %||% NA,
        computation_time = result$computation_time,
        stringsAsFactors = FALSE
      )
      
      processed_data <- rbind(processed_data, row_data)
    }
  }
  
  return(processed_data)
}

# Helper function for null coalescing
`%||%` <- function(x, y) if (is.null(x) || is.na(x)) y else x

#' Generate comprehensive analysis plots
#' @param processed_data Processed simulation results
#' @return List of ggplot objects
generate_analysis_plots <- function(
  processed_data
  )
{
  plots <- list()
  
  # 1. Overall Delta performance by prior type
  overall_delta_data <- processed_data[processed_data$component_type == "delta_overall", ]
  
  if (nrow(overall_delta_data) > 0) {
  # Aggregate by configuration
  delta_summary_data <- overall_delta_data %>%
    group_by(n_samples, n_nodes, n_covariates, prior_type, symm_method) %>%
    summarise(
      mean_TPR = mean(TPR, na.rm = TRUE),
      mean_FPR = mean(FPR, na.rm = TRUE),
      mean_F1 = mean(F1, na.rm = TRUE),
      mean_Frobenius = mean(Frobenius_Error, na.rm = TRUE),
      sd_TPR = sd(TPR, na.rm = TRUE),
      sd_FPR = sd(FPR, na.rm = TRUE),
      sd_F1 = sd(F1, na.rm = TRUE),
      sd_Frobenius = sd(Frobenius_Error, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # TPR comparison for Delta
  plots$delta_tpr_comparison <- ggplot(delta_summary_data, aes(x = factor(n_samples), y = mean_TPR, 
                                                               fill = prior_type, color = prior_type)) +
    geom_boxplot(position = position_dodge(0.8)) +
    geom_errorbar(aes(ymin = pmax(0, mean_TPR - sd_TPR), 
                      ymax = pmin(1, mean_TPR + sd_TPR), 
                      colour = prior_type),
                  position = position_dodge(0.8), width = 0.2) +
    facet_grid(n_nodes ~ n_covariates, labeller = label_both) +
    labs(title = "Delta Networks: True Positive Rate by Prior Type",
         x = "Sample Size", y = "True Positive Rate",
         fill = "Prior Type") +
    guides(colour = "none")+
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # F1 Score comparison for Delta
  plots$delta_f1_comparison <- ggplot(delta_summary_data, aes(x = factor(n_samples), y = mean_F1, 
                                                              fill = prior_type, color = prior_type)) +
    geom_boxplot(position = position_dodge(0.8)) +
    geom_errorbar(aes(ymin = pmax(0, mean_F1 - sd_F1), 
                      ymax = pmin(1, mean_F1 + sd_F1),
                      color = prior_type),
                  position = position_dodge(0.8), width = 0.2) +
    facet_grid(n_nodes ~ n_covariates, labeller = label_both) +
    labs(title = "Delta Networks: F1 Score by Prior Type",
         x = "Sample Size", y = "F1 Score",
         fill = "Prior Type") +
    guides(colour = "none")+
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Estimation Error comparison for Delta
  plots$delta_error_comparison <- ggplot(delta_summary_data, aes(x = factor(n_samples), y = mean_Frobenius, 
                                                                 fill = prior_type, color= prior_type)) +
    geom_boxplot(position = position_dodge(0.8)) +
    geom_errorbar(aes(ymin = pmax(0, mean_Frobenius - sd_Frobenius), 
                      ymax = mean_Frobenius + sd_Frobenius,
                      color= prior_type),
                  position = position_dodge(0.8), width = 0.2) +
    facet_grid(n_nodes ~ n_covariates, labeller = label_both, scales = "free_y") +
    labs(title = "Delta Networks: Estimation Error (Frobenius Norm) by Prior Type",
         x = "Sample Size", y = "Frobenius Error",
         fill = "Prior Type") +
    guides(colour = "none")+
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
  
  # 2. Overall Mean performance by prior type
  overall_mean_data <- processed_data[processed_data$component_type == "mean_overall", ]
  
  if (nrow(overall_mean_data) > 0) {
    # Aggregate by configuration
    mean_summary_data <- overall_mean_data %>%
      group_by(n_samples, n_nodes, n_covariates, prior_type, symm_method) %>%
      summarise(
        mean_TPR = mean(TPR, na.rm = TRUE),
        mean_FPR = mean(FPR, na.rm = TRUE),
        mean_F1 = mean(F1, na.rm = TRUE),
        mean_Vector_Error = mean(Vector_Error, na.rm = TRUE),
        mean_Relative_Error = mean(Relative_Error, na.rm = TRUE),
        sd_TPR = sd(TPR, na.rm = TRUE),
        sd_FPR = sd(FPR, na.rm = TRUE),
        sd_F1 = sd(F1, na.rm = TRUE),
        sd_Vector_Error = sd(Vector_Error, na.rm = TRUE),
        sd_Relative_Error = sd(Relative_Error, na.rm = TRUE),
        .groups = 'drop'
      )
    
    # TPR comparison for Mean effects
    plots$mean_tpr_comparison <- ggplot(mean_summary_data, aes(x = factor(n_samples), y = mean_TPR, 
                                                      fill = prior_type, color = prior_type)) +
      geom_boxplot(position = position_dodge(0.8)) +
      geom_errorbar(aes(ymin = pmax(0, mean_TPR - sd_TPR), 
                       ymax = pmin(1, mean_TPR + sd_TPR,  color = prior_type)),
                   position = position_dodge(0.8), width = 0.2) +
      facet_grid(n_nodes ~ n_covariates, labeller = label_both) +
      labs(title = "Mean Effects: True Positive Rate by Prior Type",
           x = "Sample Size", y = "True Positive Rate",
           fill = "Prior Type") +
      guides(colour = "none")+
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    # F1 Score comparison for Mean effects
    plots$mean_f1_comparison <- ggplot(mean_summary_data, aes(x = factor(n_samples), y = mean_F1, 
                                                     fill = prior_type,  color = prior_type)) +
      geom_boxplot(position = position_dodge(0.8)) +
      geom_errorbar(aes(ymin = pmax(0, mean_F1 - sd_F1), 
                       ymax = pmin(1, mean_F1 + sd_F1,  color = prior_type)),
                   position = position_dodge(0.8), width = 0.2) +
      facet_grid(n_nodes ~ n_covariates, labeller = label_both) +
      labs(title = "Mean Effects: F1 Score by Prior Type",
           x = "Sample Size", y = "F1 Score",
           fill = "Prior Type") +
      guides(colour = "none")+
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    # Relative Error comparison for Mean effects
    plots$mean_error_comparison <- ggplot(mean_summary_data, aes(x = factor(n_samples), y = mean_Relative_Error, 
                                                        fill = prior_type,  color = prior_type)) +
      geom_boxplot(position = position_dodge(0.8)) +
      geom_errorbar(aes(ymin = pmax(0, mean_Relative_Error - sd_Relative_Error), 
                       ymax = mean_Relative_Error + sd_Relative_Error,  color = prior_type),
                   position = position_dodge(0.8), width = 0.2) +
      facet_grid(n_nodes ~ n_covariates, labeller = label_both, scales = "free_y") +
      labs(title = "Mean Effects: Relative Estimation Error by Prior Type",
           x = "Sample Size", y = "Relative Error",
           fill = "Prior Type") +
      guides(colour = "none")+
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  
  # 3. Component-wise analysis for Delta
  delta_component_data <- processed_data[processed_data$component_type == "delta_individual", ]
  
  if (nrow(delta_component_data) > 0) {
    delta_comp_summary <- delta_component_data %>%
      group_by(component, prior_type, n_nodes) %>%
      summarise(
        mean_TPR = mean(TPR, na.rm = TRUE),
        mean_F1 = mean(F1, na.rm = TRUE),
        mean_Frobenius = mean(Frobenius_Error, na.rm = TRUE),
        mean_accuracy = mean(Accuracy, na.rm = TRUE),
        sd_TPR = sd(TPR, na.rm = TRUE),
        sd_F1 = sd(F1, na.rm = TRUE),
        sd_Frobenius = sd(Frobenius_Error, na.rm = TRUE),
        sd_accuracy = sd(Accuracy, na.rm = TRUE),
        .groups = 'drop'
      )
    
    plots$delta_component_tpr <- ggplot(delta_comp_summary, aes(x = component, y = mean_TPR, 
                                                       fill = prior_type)) +
      geom_bar(stat = "identity", position = position_dodge()) +
      geom_errorbar(aes(ymin = pmax(0, mean_TPR - sd_TPR), 
                       ymax = mean_TPR + sd_TPR),
                   position = position_dodge(0.8), width = 0.2) +
      facet_wrap(~ n_nodes, labeller = label_both) +
      labs(title = "Delta Components: TPR by Component and Prior Type",
           x = "Component", y = "Mean TPR",
           fill = "Prior Type") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    plots$delta_component_F1 <- ggplot(delta_comp_summary, aes(x = component, y = mean_F1, 
                                                       fill = prior_type)) +
      geom_bar(stat = "identity", position = position_dodge()) +
      geom_errorbar(aes(ymin = pmax(0, mean_F1 - sd_F1), 
                       ymax = mean_F1 + sd_F1),
                   position = position_dodge(0.8), width = 0.2) +
      facet_wrap(~ n_nodes, labeller = label_both) +
      labs(title = "Delta Components: F1 by Component and Prior Type",
           x = "Component", y = "Mean F1",
           fill = "Prior Type") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

      plots$delta_component_Acc <- ggplot(delta_comp_summary, aes(x = component, y = mean_accuracy, 
                                                       fill = prior_type)) +
      geom_bar(stat = "identity", position = position_dodge()) +
      geom_errorbar(aes(ymin = pmax(0, mean_accuracy - sd_accuracy), 
                       ymax = mean_accuracy + sd_accuracy),
                   position = position_dodge(0.8), width = 0.2) +
      facet_wrap(~ n_nodes, labeller = label_both) +
      labs(title = "Delta Components: Accuracy by Component and Prior Type",
           x = "Component", y = "Mean accuracy",
           fill = "Prior Type") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

      plots$delta_component_error <- ggplot(delta_comp_summary, aes(x = component, y = mean_Frobenius, 
                                                       fill = prior_type)) +
      geom_bar(stat = "identity", position = position_dodge()) +
      geom_errorbar(aes(ymin = pmax(0, mean_Frobenius - sd_Frobenius), 
                       ymax = mean_Frobenius + sd_Frobenius),
                   position = position_dodge(0.8), width = 0.2) +
      facet_wrap(~ n_nodes, labeller = label_both) +
      labs(title = "Delta Components: Frobenius by Component and Prior Type",
           x = "Component", y = "Mean Frobenius",
           fill = "Prior Type") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  
  # 4. Component-wise analysis for Mean effects
  mean_component_data <- processed_data[processed_data$component_type == "mean_individual", ]
  
  if (nrow(mean_component_data) > 0) {
    mean_comp_summary <- mean_component_data %>%
      group_by(component, prior_type, n_nodes) %>%
      summarise(
        mean_TPR = mean(TPR, na.rm = TRUE),
        mean_F1 = mean(F1, na.rm = TRUE),
        mean_Relative_Error = mean(Relative_Error, na.rm = TRUE),
        sd_TPR = sd(TPR, na.rm = TRUE),
        sd_F1 = sd(F1, na.rm = TRUE),
        sd_Relative_Error = sd(Relative_Error, na.rm = TRUE),

        .groups = 'drop'
      )
    
    plots$mean_component_tpr <- ggplot(mean_comp_summary, aes(x = component, y = mean_TPR, 
                                                     fill = prior_type)) +
      geom_bar(stat = "identity", position = position_dodge()) +
      geom_errorbar(aes(ymin = pmax(0, mean_TPR - sd_TPR), 
                       ymax = mean_TPR + sd_TPR),
                   position = position_dodge(0.8), width = 0.2) +
      facet_wrap(~ n_nodes, labeller = label_both) +
      labs(title = "Mean Effects Components: TPR by Component and Prior Type",
           x = "Component", y = "Mean TPR",
           fill = "Prior Type") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    plots$mean_component_error <- ggplot(mean_comp_summary, aes(x = component, y = mean_Relative_Error, 
                                                       fill = prior_type)) +
      geom_bar(stat = "identity", position = position_dodge()) +
      geom_errorbar(aes(ymin = pmax(0, mean_Relative_Error - sd_Relative_Error), 
                       ymax = mean_Relative_Error + sd_Relative_Error),
                   position = position_dodge(0.8), width = 0.2) +
      facet_wrap(~ n_nodes, labeller = label_both, scales = "free_y") +
      labs(title = "Mean Effects Components: Relative Error by Component and Prior Type",
           x = "Component", y = "Mean Relative Error",
           fill = "Prior Type") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  
  # 5. Symmetrization method comparison (Delta and Mean)
  if (length(unique(processed_data$symm_method)) > 1) {
    # Delta symmetrization comparison
    delta_symm_data <- overall_delta_data %>%
      group_by(symm_method, prior_type, n_samples, n_nodes) %>%
      summarise(
        mean_F1 = mean(F1, na.rm = TRUE),
        sd_F1 = sd(F1, na.rm = TRUE),
        .groups = 'drop'
      )
    
    plots$delta_symmetrization_comparison <- ggplot(delta_symm_data, aes(x = symm_method, y = mean_F1, 
                                                            fill = prior_type, color= prior_type)) +
      geom_boxplot(position = position_dodge(0.8)) +
      geom_errorbar(aes(ymin = pmax(0, mean_F1 - sd_F1), 
                       ymax = mean_F1 + sd_F1, color= prior_type),
                   position = position_dodge(0.8), width = 0.2) +
      facet_grid(n_samples ~ n_nodes, labeller = label_both) +
      labs(title = "Delta Networks: F1 Score by Symmetrization Method",
           x = "Symmetrization Method", y = "Mean F1 Score",
           fill = "Prior Type") +
      guides(colour = "none")+
      theme_minimal()
  }
  
  # 6. Combined Delta vs Mean performance comparison
  combined_data <- processed_data[processed_data$component_type %in% c("delta_overall", "mean_overall"), ]
  
  if (nrow(combined_data) > 0) {
    combined_summary <- combined_data %>%
      mutate(metric_type = ifelse(component_type == "delta_overall", "Delta Networks", "Mean Effects")) %>%
      group_by(metric_type, prior_type, n_samples, n_nodes) %>%
      summarise(
        mean_F1 = mean(F1, na.rm = TRUE),
        mean_TPR = mean(TPR, na.rm = TRUE),
        sd_F1 = sd(F1, na.rm = TRUE),
        sd_TPR = sd(TPR, na.rm = TRUE),
        .groups = 'drop'
      )
    
    plots$combined_f1_comparison <- ggplot(combined_summary, aes(x = factor(n_samples), y = mean_F1, 
                                                        fill = prior_type, color=prior_type)) +
      geom_boxplot(position = position_dodge(0.8)) +
      geom_errorbar(aes(ymin = pmax(0, mean_F1 - sd_F1), 
                       ymax = mean_F1 + sd_F1, color= prior_type),
                   position = position_dodge(0.8), width = 0.2) +
      facet_grid(metric_type ~ n_nodes, labeller = label_both, scales = "free_y") +
      labs(title = "F1 Score Comparison: Delta Networks vs Mean Effects",
           x = "Sample Size", y = "Mean F1 Score",
           fill = "Prior Type") +
      guides(colour = "none")+
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  
  return(plots)
}

#' Generate summary tables
#' @param processed_data Processed simulation results
#' @return List of summary tables
generate_summary_tables <- function(processed_data) {
  
  tables <- list()
  
  # 1. Overall Delta performance table
  overall_delta_data <- processed_data[processed_data$component_type == "delta_overall", ]
  
  if (nrow(overall_delta_data) > 0) {
    tables$overall_delta_summary <- overall_delta_data %>%
      group_by(n_samples, n_nodes, n_covariates, prior_type) %>%
      summarise(
        Mean_TPR = sprintf("%.3f (%.3f)", mean(TPR, na.rm = TRUE), sd(TPR, na.rm = TRUE)),
        Mean_FPR = sprintf("%.3f (%.3f)", mean(FPR, na.rm = TRUE), sd(FPR, na.rm = TRUE)),
        Mean_F1 = sprintf("%.3f (%.3f)", mean(F1, na.rm = TRUE), sd(F1, na.rm = TRUE)),
        Mean_Error = sprintf("%.3f (%.3f)", mean(Frobenius_Error, na.rm = TRUE), 
                            sd(Frobenius_Error, na.rm = TRUE)),
        Mean_Time = sprintf("%.1f", mean(computation_time, na.rm = TRUE)),
        .groups = 'drop'
      ) %>%
      arrange(n_nodes, n_samples, prior_type)
  }
  
  # 2. Overall Mean effects performance table
  overall_mean_data <- processed_data[processed_data$component_type == "mean_overall", ]
  
  if (nrow(overall_mean_data) > 0) {
    tables$overall_mean_summary <- overall_mean_data %>%
      group_by(n_samples, n_nodes, n_covariates, prior_type) %>%
      summarise(
        Mean_TPR = sprintf("%.3f (%.3f)", mean(TPR, na.rm = TRUE), sd(TPR, na.rm = TRUE)),
        Mean_FPR = sprintf("%.3f (%.3f)", mean(FPR, na.rm = TRUE), sd(FPR, na.rm = TRUE)),
        Mean_F1 = sprintf("%.3f (%.3f)", mean(F1, na.rm = TRUE), sd(F1, na.rm = TRUE)),
        Mean_Vector_Error = sprintf("%.3f (%.3f)", mean(Vector_Error, na.rm = TRUE), 
                                   sd(Vector_Error, na.rm = TRUE)),
        Mean_Relative_Error = sprintf("%.3f (%.3f)", mean(Relative_Error, na.rm = TRUE), 
                                     sd(Relative_Error, na.rm = TRUE)),
        Mean_Time = sprintf("%.1f", mean(computation_time, na.rm = TRUE)),
        .groups = 'drop'
      ) %>%
      arrange(n_nodes, n_samples, prior_type)
  }
  
  # 3. Prior type comparison table (Delta networks)
  if (nrow(overall_delta_data) > 0) {
    tables$delta_prior_comparison <- overall_delta_data %>%
      group_by(prior_type) %>%
      summarise(
        Overall_TPR = sprintf("%.3f ± %.3f", mean(TPR, na.rm = TRUE), sd(TPR, na.rm = TRUE)),
        Overall_FPR = sprintf("%.3f ± %.3f", mean(FPR, na.rm = TRUE), sd(FPR, na.rm = TRUE)),
        Overall_F1 = sprintf("%.3f ± %.3f", mean(F1, na.rm = TRUE), sd(F1, na.rm = TRUE)),
        Overall_Error = sprintf("%.3f ± %.3f", mean(Frobenius_Error, na.rm = TRUE), 
                               sd(Frobenius_Error, na.rm = TRUE)),
        Avg_Time = sprintf("%.1f s", mean(computation_time, na.rm = TRUE)),
        .groups = 'drop'
      )
  }
  
  # 4. Prior type comparison table (Mean effects)
  if (nrow(overall_mean_data) > 0) {
    tables$mean_prior_comparison <- overall_mean_data %>%
      group_by(prior_type) %>%
      summarise(
        Overall_TPR = sprintf("%.3f ± %.3f", mean(TPR, na.rm = TRUE), sd(TPR, na.rm = TRUE)),
        Overall_FPR = sprintf("%.3f ± %.3f", mean(FPR, na.rm = TRUE), sd(FPR, na.rm = TRUE)),
        Overall_F1 = sprintf("%.3f ± %.3f", mean(F1, na.rm = TRUE), sd(F1, na.rm = TRUE)),
        Overall_Vector_Error = sprintf("%.3f ± %.3f", mean(Vector_Error, na.rm = TRUE), 
                                      sd(Vector_Error, na.rm = TRUE)),
        Overall_Relative_Error = sprintf("%.3f ± %.3f", mean(Relative_Error, na.rm = TRUE), 
                                        sd(Relative_Error, na.rm = TRUE)),
        Avg_Time = sprintf("%.1f s", mean(computation_time, na.rm = TRUE)),
        .groups = 'drop'
      )
  }
  
  # 5. Component-wise performance for Delta
  delta_component_data <- processed_data[processed_data$component_type == "delta_individual", ]
  
  if (nrow(delta_component_data) > 0) {
    tables$delta_component_summary <- delta_component_data %>%
      group_by(component, prior_type) %>%
      summarise(
        TPR = sprintf("%.3f ± %.3f", mean(TPR, na.rm = TRUE), sd(TPR, na.rm = TRUE)),
        FPR = sprintf("%.3f ± %.3f", mean(FPR, na.rm = TRUE), sd(FPR, na.rm = TRUE)),
        F1 = sprintf("%.3f ± %.3f", mean(F1, na.rm = TRUE), sd(F1, na.rm = TRUE)),
        Accuracy = sprintf("%.3f ± %.3f", mean(Accuracy, na.rm = TRUE), sd(Accuracy, na.rm = TRUE)),
        Error = sprintf("%.3f ± %.3f", mean(Frobenius_Error, na.rm = TRUE), 
                       sd(Frobenius_Error, na.rm = TRUE)),
        .groups = 'drop'
      ) %>%
      pivot_wider(names_from = prior_type, values_from = c(TPR, F1, Error),
                 names_sep = "_")

      tables$delta_component_complete_summary <- delta_component_data %>%
      group_by(n_samples, n_nodes, n_covariates, prior_type, symm_method, component) %>%
      summarise(
        TPR = sprintf("%.3f ± %.3f", mean(TPR, na.rm = TRUE), sd(TPR, na.rm = TRUE)),
        FPR = sprintf("%.3f ± %.3f", mean(FPR, na.rm = TRUE), sd(FPR, na.rm = TRUE)),
        F1 = sprintf("%.3f ± %.3f", mean(F1, na.rm = TRUE), sd(F1, na.rm = TRUE)),
        Accuracy = sprintf("%.3f ± %.3f", mean(Accuracy, na.rm = TRUE), sd(Accuracy, na.rm = TRUE)),
        Error = sprintf("%.3f ± %.3f", mean(Frobenius_Error, na.rm = TRUE), 
                        sd(Frobenius_Error, na.rm = TRUE)),
        .groups = 'drop'
      )
  }
  
  # 6. Component-wise performance for Mean effects
  mean_component_data <- processed_data[processed_data$component_type == "mean_individual", ]
  
  if (nrow(mean_component_data) > 0) {
    tables$mean_component_summary <- mean_component_data %>%
      group_by(component, prior_type) %>%
      summarise(
        TPR = sprintf("%.3f ± %.3f", mean(TPR, na.rm = TRUE), sd(TPR, na.rm = TRUE)),
        F1 = sprintf("%.3f ± %.3f", mean(F1, na.rm = TRUE), sd(F1, na.rm = TRUE)),
        Vector_Error = sprintf("%.3f ± %.3f", mean(Vector_Error, na.rm = TRUE), 
                              sd(Vector_Error, na.rm = TRUE)),
        Relative_Error = sprintf("%.3f ± %.3f", mean(Relative_Error, na.rm = TRUE), 
                                sd(Relative_Error, na.rm = TRUE)),
        .groups = 'drop'
      ) %>%
      pivot_wider(names_from = prior_type, values_from = c(TPR, F1, Vector_Error, Relative_Error),
                 names_sep = "_")
  }
  
  # 7. Combined comparison table (Delta vs Mean overall performance)
  combined_data <- processed_data[processed_data$component_type %in% c("delta_overall", "mean_overall"), ]
  
  if (nrow(combined_data) > 0) {
    tables$combined_comparison <- combined_data %>%
      mutate(metric_type = ifelse(component_type == "delta_overall", "Delta", "Mean")) %>%
      group_by(metric_type, prior_type) %>%
      summarise(
        TPR = sprintf("%.3f ± %.3f", mean(TPR, na.rm = TRUE), sd(TPR, na.rm = TRUE)),
        F1 = sprintf("%.3f ± %.3f", mean(F1, na.rm = TRUE), sd(F1, na.rm = TRUE)),
        Primary_Error = sprintf("%.3f ± %.3f", 
                               mean(ifelse(is.na(Frobenius_Error), Vector_Error, Frobenius_Error), na.rm = TRUE),
                               sd(ifelse(is.na(Frobenius_Error), Vector_Error, Frobenius_Error), na.rm = TRUE)),
        .groups = 'drop'
      ) %>%
      pivot_wider(names_from = prior_type, values_from = c(TPR, F1, Primary_Error),
                 names_sep = "_")
  }
  
  return(tables)
}

#################################################
## SIMULATION EXECUTION FUNCTIONS - Parrelized on HPC
#################################################

#' Generate all simulation datasets across all configurations, replications
#' @param SIMULATION_GRID Dataframe with simulation configurations
#' @param N_REPLICATIONS Number of replications per configuration
#' @param output_folder Base folder to save simulation input data
#' @param output_file Name of the output file to save the list of input data for computations
#' @return Dataframe with paths to input data for computations
generate_input_datasets_simulation <- function(
  SIMULATION_GRID,
  N_REPLICATIONS,
  output_folder,
  output_file
  )
{
  cat("Starting comprehensive simulation study...\n")
  cat(sprintf("Total configurations: %d\n", nrow(SIMULATION_GRID)))
  cat(sprintf("Replications per configuration: %d\n", N_REPLICATIONS))
  cat(sprintf("Total simulations: %d\n", nrow(SIMULATION_GRID) * N_REPLICATIONS))
  
  total_sims <- nrow(SIMULATION_GRID) * N_REPLICATIONS
  current_sim <- 0
  
  input_computation <- data.frame(p=NA,slurm_script_path=NA, input_data_path=NA, output_path=NA, name_output=NA)
  
  # Create input data for each simulation configuration and replication
  for (i in 1:nrow(SIMULATION_GRID)) {
    config <- SIMULATION_GRID[i, ]
    ggReg_output_folder <- sprintf("n%d_p%d_q%d_%s_%s", 
                                   config$n_samples, config$n_nodes, config$n_covariates,
                                   config$prior_type, config$symm_method)
    ggReg_output_path <- paste(output_folder,"/", ggReg_output_folder,"/", sep ="")
    dir.create(ggReg_output_path)

    for (rep in 1:N_REPLICATIONS) {
      current_sim <- current_sim + 1

      if (current_sim %% 10 == 0) {
        cat(sprintf("Progress: %d/%d (%.1f%%)\n", 
                    current_sim, total_sims, 100 * current_sim / total_sims))
      }
      
      seed <- 12345 + i * 1000 + rep

      ggReg_output_name <- sprintf("rep%d", rep)
      ggReg_output_path_rep <- paste(ggReg_output_path, ggReg_output_name, "/", sep="")
      dir.create(ggReg_output_path_rep)
    
      cat(sprintf("Running: n=%d, p=%d, q=%d, prior=%s, symm=%s, rep=%d\n",
                    config$n_samples, config$n_nodes, config$n_covariates, 
                    config$prior_type, config$symm_method, rep))
      
        
      network_config <- NETWORK_CONFIG
      network_config$population_edges <- round(config$n_nodes * (config$n_nodes - 1) / 2 * 
                                                     network_config$population_density)
          
      covariate_config <- base_covariate_config(config$n_covariates)
          
      prior_config <- list(
        type = config$prior_type,
        noise_params = PRIOR_NOISE_PARAMS
      )
          
      method_params <- get_method_params(config$prior_type)
          
      # Generate simulation dataset
      sim_dataset <- generate_simulation_dataset(
            n_samples = config$n_samples,
            n_nodes = config$n_nodes,
            covariate_config = covariate_config,
            network_config = network_config,
            prior_config = prior_config,
            seed = seed,
            output_path = ggReg_output_path_rep,
            name_output = ggReg_output_name
            
      )
      rep_id=rep
      save(sim_dataset, config, rep_id, file=paste0(ggReg_output_path_rep, "sim_dataset_full.RData"))
      
      sim_data <- sim_dataset$data
      protein_cols <- colnames(sim_data)[grepl("^Prot", colnames(sim_data))]
      covariate_cols <- colnames(sim_data)[grepl("^x", colnames(sim_data))]
          
      # Prepare covariates dataframe
      if (length(covariate_cols) > 0) {
          covariates_df <- sim_data[, covariate_cols, drop = FALSE]
      } else {
          covariates_df <- NULL
      }

      if(is.null(sim_data$true_mean_effects)){
        mean_estimation=FALSE
      }else{mean_estimation=TRUE}

      x = sim_data[, protein_cols]
      known_ppi = sim_dataset$prior_knowledge
      covariates = covariates_df
      scr = method_params$screening_procedure
      gamma = NULL
      mean_estimation = mean_estimation
      lambda_mean = NULL
      lambda_mean_type = "1se"
      lambda_prec = NULL
      lambda_prec_type = method_params$lambda_prec_type
      tune_hyperparams = method_params$tune_hyperparams
      asparse_grid = method_params$asparse_grid
      weight_grid = method_params$weight_grid
      random_hyper_search = method_params$random_hyper_search
      p.rand.hyper = method_params$p.rand.hyper
      K = method_params$K
      use_slurm = FALSE
      slurm_script_path = method_params$slurm_script_path
      output_path = output_path
      name_output = name_output
      symm_method =config$symm_method
      verbose = method_params$verbose

      if(mean_estimation){
        res_mean_reg <- GGReg_mean_estimation(
          x = x,
          covariates = covariates,
          lambda_mean = lambda_mean,
          lambda_mean_type = lambda_mean_type,
          verbose = verbose)
        Z <- res_mean_reg$z
      }else {
         Z <- x
      }
      n <- nrow(Z)
      p <- ncol(Z)

      input_data_path <- paste0(ggReg_output_path_rep, "input_data_nodes.rda")
      save(Z, known_ppi, covariates, scr, gamma, lambda_prec, lambda_prec_type, 
          tune_hyperparams, asparse_grid, weight_grid, random_hyper_search, p.rand.hyper, K,
          file = input_data_path)

      if (verbose) {
        cat("Saved input data to:", input_data_path, "\n")
      }
        
      temp_input_computation <- data.frame(p=p,slurm_script_path=slurm_script_path, input_data_path=input_data_path, output_path=ggReg_output_path_rep, name_output=ggReg_output_name)
        
      input_computation <- rbind(input_computation, temp_input_computation)
    }
  }

  write.csv(input_computation[-1,], file = paste0(output_folder, "/","input_computation_file.csv"), row.names = F)
  cat("The input dataset have been generated and saved successfully in: ", paste0(output_folder, "/","input_computation_file.csv"), "\n")
  return(input_computation[-1,])
}

#' Collect results from GGReg parallel computations and evaluate performance
#' @param p Number of nodes
#' @param output_path Path to the output folder where GGReg results are saved
#' @param name_output Name of the output file to identify the results
#' @param symm_method Symmetrization method used in GGReg estimation
#' @return List with evaluation results and computational time
#' 
collect_and_evaluate_resuts <- function(
  p,
  output_path,
  name_output,
  symm_method = "OR"
  )
{
  ggReg_results <- collect_node_results(p, output_path, name_output, symm_method)
  
  list_estimated_prec <- ggReg_results$additional_info$Dic_Delta_hat
  estimation_matrix_names <- names(list_estimated_prec)
  for(l in 1:length(estimation_matrix_names)){
    name_cov <- estimation_matrix_names[l]
    adj_df <- melt(abs(list_estimated_prec[[name_cov]]))
    colnames(adj_df) <- c("Row", "Col", "Value")
    adj_df$Value <- as.numeric(adj_df$Value)
    adj_df$Row <- as.numeric(sub("Prot", "", adj_df$Row))
    if(l==1){
      adj_df$Col <- as.numeric(sub("Prot", "", adj_df$Col))
    } else {
      adj_df$Col <- as.numeric(sub(paste(name_cov,":","Prot", sep=""), "", adj_df$Col))
    }
    # Plot
    plot <- ggplot(adj_df, aes(x = Col, y = Row, fill = Value)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "black") +
      theme_minimal() +
      coord_fixed() +
      scale_y_reverse() +  # to match matrix view
      labs(title = "Estimated adjacency Matrix", x = "", y = "")
    ggsave(plot, filename=paste(output_path, "Estimated_matrix_",name_output ,"_", name_cov, ".png", sep=""), width=8, height=8, dpi=300)
  }


  # Evaluate results
  evaluation_results <- evaluate_estimation(
    true_Delta_list = sim_dataset$true_precision_matrices,
    estimated_Delta_list = list_estimated_prec,
    true_mean_effects = NULL,
    estimated_mean_effects = NULL
  )

  computation_time <- as.numeric(mean(ggReg_results$additional_info$computational_time_per_node))
  load(paste0(output_path,"sim_dataset_full.RData"))

  # Compile results
  result <- list(
    config = config,
    rep_id = rep_id,
    evaluation = evaluation_results,
    computation_time = computation_time
  )
  return(result)
}

