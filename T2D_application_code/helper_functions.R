# Helper functions of Network_estimation_case_study_T2D_UKB.R: confounder
# analysis (variance of the protein PCs explained by candidate covariates,
# association with incident T2D) and the cohort split.
library(dplyr)
library(pheatmap)
library(tidyr)
library(caret)
library(knitr)
library(kableExtra)

# Function to calculate R-squared for PC ~ confounding factor
calculate_pc_rsquared <- function(
    pc_name, 
    conf_factor, 
    data) 
  {
  tryCatch({
    # Handle different types of confounding factors
    if (is.character(data[[conf_factor]]) || is.factor(data[[conf_factor]])) {
      # For categorical variables
      formula_str <- paste(pc_name, "~", conf_factor)
    } else {
      # For continuous variables
      formula_str <- paste(pc_name, "~", conf_factor)
    }
    
    # Fit linear model
    model <- lm(as.formula(formula_str), data = data)
    model_summary <- summary(model)
    
    # Return adjusted R-squared (multiply by 100 for percentage)
    return(model_summary$r.squared * 100)
    
  }, error = function(e) {
    cat("Error calculating R² for", pc_name, "~", conf_factor, ":", e$message, "\n")
    return(NA)
  })
}

# COMPUTING PRINCIPAL COMPONENTS 
perform_confounder_analysis <- function(
    mol_data, 
    pheno_data, 
    col_id,  
    conf_factor) 
  {
moldata <- mol_data[, !names(mol_data) %in% col_id]
pca_result <- prcomp(moldata, center = TRUE, scale. = TRUE)
pc_scores <- as.data.frame(pca_result$x)

# Calculate variance explained by each PC
variance_explained <- (pca_result$sdev^2 / sum(pca_result$sdev^2))*100
names(variance_explained) <- colnames(pc_scores)
significant_pcs <- names(variance_explained)[variance_explained > 1]
cat("\nPCs explaining >1% variance:", length(significant_pcs), "\n")

# Add participant IDs
pc_scores <- pc_scores[,significant_pcs]
pc_scores[[col_id]] <- mol_data[[col_id]]
dtjoin<-dtjoin <- pc_scores %>%
  left_join(pheno_data %>%
              dplyr::select(dplyr::all_of(c(col_id, conf_factor))),
            by = col_id
  )

# Initialize results matrix
results_matrix <- matrix(NA, 
                         nrow = length(significant_pcs), 
                         ncol = length(conf_factor),
                         dimnames = list(significant_pcs, conf_factor))

cat("Calculating R² values for each PC-confounding factor combination...\n")

# Calculate R-squared for each combination
for (pc in significant_pcs) {
  for (factor in conf_factor) {
    cat("Processing:", pc, "~", factor, "\n")
    rsq_value <- calculate_pc_rsquared(pc, factor, dtjoin)
    results_matrix[pc, factor] <- rsq_value
  }
}

# Convert to data frame for easier handling
results_df <- as.data.frame(results_matrix)

# Prepare data for heatmap (transpose to match your example format)
heatmap_data <- t(results_matrix)

# Create a better color palette
# Use RdBu (Red-Blue) palette similar to your example
colors <- colorRampPalette(c( "white", "orange", "darkorange", "red"))(100)

# Create the main heatmap using pheatmap
pheatmap(
  heatmap_data,
  color = colors,
  breaks = seq(0, 
               max(abs(heatmap_data), na.rm = TRUE), 
               length.out = 101),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.2f",
  fontsize_number = 8,
  main = "PC Confounding Factor Analysis\nR² Values (%)",
  fontsize = 10,
  angle_col = 45,
  cellwidth = 25,
  cellheight = 20,
  filename = "Conf_analysis_heatmap_plot.png",
  na_col = "grey90"
)

return(results_df)
}

split_dataset <- function(
    data,
    sex_col,
    age_col,
    endpoint_col,
    incident_name,
    non_incident_name,
    p_sel = 0.8,
    p_case_train = 0.5,
    clinical_var = NULL,
    controls_per_case = 4,
    seed = 123,
    output_dir = "."
) {
  set.seed(seed)
  
  # ------------------------------------------------------------------
  # 1. Split into selection and remaining data
  #    Stratify jointly by endpoint and sex
  # ------------------------------------------------------------------
  
  stratification_group <- interaction(
    data[[endpoint_col]],
    data[[sex_col]],
    drop = TRUE
  )
  
  index_1 <- caret::createDataPartition(
    y = stratification_group,
    p = p_sel,
    list = FALSE,
    times = 1
  )
  
  selection_data <- data[index_1, , drop = FALSE]
  remaining_data <- data[-index_1, , drop = FALSE]
  
  # ------------------------------------------------------------------
  # 2. Separate incident and non-incident participants
  # ------------------------------------------------------------------
  
  incident_participants <- remaining_data %>%
    dplyr::filter(.data[[endpoint_col]] == incident_name)
  
  non_incident_participants <- remaining_data %>%
    dplyr::filter(.data[[endpoint_col]] == non_incident_name)
  
  if (nrow(incident_participants) == 0) {
    stop("No incident participants were found in the remaining data.")
  }
  
  if (nrow(non_incident_participants) == 0) {
    stop("No non-incident participants were found in the remaining data.")
  }
  
  # ------------------------------------------------------------------
  # 3. Split incident participants into training and testing
  #    Stratified by sex
  # ------------------------------------------------------------------
  
  set.seed(seed + 1)
  
  index_2 <- caret::createDataPartition(
    y = incident_participants[[sex_col]],
    p = p_case_train,
    list = FALSE,
    times = 1
  )
  
  incident_training_data <- incident_participants[
    index_2,
    ,
    drop = FALSE
  ]
  
  incident_testing_data <- incident_participants[
    -index_2,
    ,
    drop = FALSE
  ]
  
  # ------------------------------------------------------------------
  # 4. Match controls to incident training participants
  #    Exact matching on sex and age
  # ------------------------------------------------------------------
  
  matched_controls <- vector(
    mode = "list",
    length = nrow(incident_training_data)
  )
  
  available_controls <- non_incident_participants
  
  set.seed(seed + 2)
  
  for (i in seq_len(nrow(incident_training_data))) {
    
    incident_row <- incident_training_data[i, , drop = FALSE]
    
    eligible_controls <- available_controls %>%
      dplyr::filter(
        .data[[age_col]] == incident_row[[age_col]][1],
        .data[[sex_col]] == incident_row[[sex_col]][1]
      )
    
    if (nrow(eligible_controls) == 0) {
      warning(
        "No matching controls found for incident participant ID ",
        incident_row$id[1],
        "."
      )
      
      matched_controls[[i]] <- NULL
      next
    }
    
    n_to_sample <- min(
      controls_per_case,
      nrow(eligible_controls)
    )
    
    selected_controls <- eligible_controls %>%
      dplyr::slice_sample(
        n = n_to_sample,
        replace = FALSE
      )
    
    matched_controls[[i]] <- selected_controls
    
    available_controls <- available_controls %>%
      dplyr::filter(!id %in% selected_controls$id)
  }
  
  matched_controls <- dplyr::bind_rows(matched_controls)
  
  train_data <- dplyr::bind_rows(
    incident_training_data,
    matched_controls
  ) %>%
    dplyr::distinct(id, .keep_all = TRUE)
  
  # All remaining participants not used for training enter the test set
  test_data <- remaining_data %>%
    dplyr::filter(!id %in% train_data$id)
  
  # ------------------------------------------------------------------
  # 5. Ensure categorical levels in test are represented in training
  # ------------------------------------------------------------------
  
  if (!is.null(clinical_var) && length(clinical_var) > 0) {
    
    is_categorical <- function(x) {
      is.factor(x) || is.character(x) || is.logical(x)
    }
    
    train_data_id <- train_data$id
    test_data_id <- test_data$id
    
    ids_to_move <- character(0)
    
    for (v in clinical_var) {
      
      x_train <- train_data[[v]]
      x_test <- test_data[[v]]
      
      if (!is_categorical(x_train) && !is_categorical(x_test)) {
        next
      }
      
      train_values <- unique(
        as.character(x_train[!is.na(x_train)])
      )
      
      test_values <- as.character(x_test)
      
      unseen_mask <- !is.na(test_values) &
        nzchar(test_values) &
        !(test_values %in% train_values)
      
      if (any(unseen_mask)) {
        
        unseen_ids <- test_data$id[unseen_mask]
        
        message(
          "Moving ",
          length(unique(unseen_ids)),
          " participant(s) because of unseen levels in `",
          v,
          "`."
        )
        
        ids_to_move <- c(ids_to_move, unseen_ids)
      }
    }
    
    ids_to_move <- unique(ids_to_move)
    
    if (length(ids_to_move) == 0) {
      
      message(
        "All categorical levels in TEST are present in TRAIN. ",
        "No IDs moved."
      )
      
    } else {
      
      train_data_id <- unique(c(train_data_id, ids_to_move))
      test_data_id <- setdiff(test_data_id, ids_to_move)
      
      train_data <- remaining_data %>%
        dplyr::filter(id %in% train_data_id)
      
      test_data <- remaining_data %>%
        dplyr::filter(id %in% test_data_id)
    }
  }
  
  # ------------------------------------------------------------------
  # 6. Summary
  # ------------------------------------------------------------------
  
  message("Selection set: ", nrow(selection_data), " participants")
  message("Training set: ", nrow(train_data), " participants")
  message("Testing set: ", nrow(test_data), " participants")
  
  # ------------------------------------------------------------------
  # 7. Save IDs
  # ------------------------------------------------------------------
  
  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  write.table(
    selection_data$id,
    file.path(output_dir, "Selection_data_id.txt"),
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
  
  write.table(
    train_data$id,
    file.path(output_dir, "Train_data_id.txt"),
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
  
  write.table(
    test_data$id,
    file.path(output_dir, "Test_data_id.txt"),
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
  
  return(
    list(
      selection_data = selection_data,
      train_data = train_data,
      test_data = test_data
    )
  )
}


summary_subset <- function(df,
                           endpoint_col,
                           incident_name,
                           sex_col,
                           age_col) {
  
  n <- nrow(df)
  
  incidence <- mean(df[[endpoint_col]] == incident_name, na.rm = TRUE) * 100
  
  age_mean <- mean(df[[age_col]], na.rm = TRUE)
  age_sd <- sd(df[[age_col]], na.rm = TRUE)
  
  sex_tab <- table(df[[sex_col]])
  
  n_female <- if ("Female" %in% names(sex_tab))
    sex_tab["Female"]
  else if ("F" %in% names(sex_tab))
    sex_tab["F"]
  else
    NA
  
  n_male <- if ("Male" %in% names(sex_tab))
    sex_tab["Male"]
  else if ("M" %in% names(sex_tab))
    sex_tab["M"]
  else
    NA
  
  tibble(
    n = n,
    `T2D Inc. (%)` = round(incidence, 1),
    `Age (mean ± sd)` = sprintf(
      "%.2f ± %.2f",
      age_mean,
      age_sd
    ),
    `Sex (F/M)` = paste0(n_female, "/", n_male)
  )
}

summary_cohort <- function(list_df, endpoint_col, incident_name, sex_col, age_col){
summary_table <- bind_rows(
  summary_subset(
    list_df$selection_data,
    endpoint_col = endpoint_col,
    incident_name = incident_name,
    sex_col = sex_col,
    age_col = age_col
  ),
  summary_subset(
    list_df$train_data,
    endpoint_col = endpoint_col,
    incident_name = incident_name,
    sex_col = sex_col,
    age_col = age_col
  ),
  summary_subset(
    list_df$test_data,
    endpoint_col = endpoint_col,
    incident_name = incident_name,
    sex_col = sex_col,
    age_col = age_col
  ),
  summary_subset(
    bind_rows(list_df$selection_data, list_df$train_data, list_df$test_data),
    endpoint_col = endpoint_col,
    incident_name = incident_name,
    sex_col = sex_col,
    age_col = age_col
  )
) %>%
  mutate(
    Subset = c(
      "Estimation",
      "Training (balanced)",
      "Test",
      "Total"
    )
  ) %>%
  relocate(Subset)
return(summary_table)
}

# Association between candidate covariates and 5-year incident T2D
perform_t2d_confounder_analysis <- function(
    model_data,
    pheno_data = NULL,
    col_id = "id",
    conf_factor,
    outcome = "surv_status_5y",
    adjust_vars = c("sex", "Age_at_recruitment"),
    categorical_vars = c("sex", "Smoking_status",
                         "Cholesterol_lowering_medication",
                         "Month_blood_sample_collected"),
    increments = c(Age_at_recruitment = 5, BMI = 5),
    output_dir = "T2D_confounder_analysis",
    make_plots = TRUE) {
  
  # 1. Match missing candidate variables by participant ID.
  data <- as.data.frame(model_data)
  needed <- unique(c(conf_factor, adjust_vars))
  if (!all(c(col_id, outcome) %in% names(data))) {
    stop("model_data must contain the ID and outcome columns.")
  }
  if (anyDuplicated(data[[col_id]])) stop("Duplicate participant IDs in model_data.")
  missing_vars <- setdiff(needed, names(data))
  if (length(missing_vars)) {
    if (is.null(pheno_data)) {
      stop("pheno_data is needed for: ", paste(missing_vars, collapse = ", "))
    }
    ph <- as.data.frame(pheno_data)
    if (!col_id %in% names(ph) || anyDuplicated(ph[[col_id]])) {
      stop("pheno_data must have unique participant IDs.")
    }
    unavailable <- setdiff(missing_vars, names(ph))
    if (length(unavailable)) {
      stop("Variables not found: ", paste(unavailable, collapse = ", "))
    }
    index <- match(as.character(data[[col_id]]), as.character(ph[[col_id]]))
    for (v in missing_vars) data[[v]] <- ph[[v]][index]
  }
  
  # 2. Require explicit binary incident-T2D coding.
  y <- toupper(trimws(as.character(data[[outcome]])))
  data$.outcome <- ifelse(y %in% c("1", "TRUE"), 1L,
                          ifelse(y %in% c("0", "FALSE"), 0L, NA_integer_))
  unexpected <- !is.na(data[[outcome]]) & is.na(data$.outcome)
  if (any(unexpected)) {
    stop("Unexpected outcome coding. Inspect unique(model_data[[outcome]]).")
  }
  
  for (v in intersect(categorical_vars, names(data))) {
    data[[v]] <- as.factor(data[[v]])
  }
  # Preserve existing factor/character variables as categorical.
  for (v in needed) {
    if (is.character(data[[v]]) || is.logical(data[[v]])) {
      data[[v]] <- as.factor(data[[v]])
    }
  }
  
  # 3. Fit nested logistic regressions on the same complete-case subset.
  #    LRT p-values are omnibus for multi-level categorical predictors.
  summary_rows <- list()
  coefficient_rows <- list()
  failures <- character(0)
  for (v in conf_factor) {
    controls <- setdiff(adjust_vars, v)
    columns <- unique(c(".outcome", v, controls))
    dd <- data[stats::complete.cases(data[, columns, drop = FALSE]),
               columns, drop = FALSE]
    dd <- droplevels(dd)
    if (nrow(dd) < 30L || length(unique(dd$.outcome)) < 2L ||
        length(unique(dd[[v]])) < 2L) {
      failures <- c(failures, paste(v, "insufficient data or variation"))
      next
    }
    # A constant adjustment factor should not cause contrasts errors.
    controls <- controls[vapply(controls, function(z) {
      length(unique(dd[[z]])) > 1L
    }, logical(1))]
    
    dd$.candidate <- dd[[v]]
    scale_label <- "category vs reference"
    reference <- NA_character_
    if (is.numeric(dd$.candidate)) {
      if (v %in% names(increments)) {
        increment <- unname(increments[[v]])
        scale_label <- paste0("per ", increment, " units")
      } else {
        increment <- stats::sd(dd$.candidate)
        scale_label <- "per 1 SD"
      }
      if (!is.finite(increment) || increment <= 0) {
        failures <- c(failures, paste(v, "invalid numeric scale"))
        next
      }
      dd$.candidate <- dd$.candidate / increment
    } else {
      dd$.candidate <- droplevels(as.factor(dd$.candidate))
      reference <- levels(dd$.candidate)[1]
    }
    
    base_formula <- stats::reformulate("1", response = ".outcome")
    single_formula <- stats::reformulate(".candidate", response = ".outcome")
    adjusted_base <- stats::reformulate(
      if (length(controls)) controls else "1", response = ".outcome")
    adjusted_full <- stats::reformulate(
      c(controls, ".candidate"), response = ".outcome")
    
    result <- tryCatch({
      m0 <- stats::glm(base_formula, data = dd, family = stats::binomial())
      m1 <- stats::glm(single_formula, data = dd, family = stats::binomial())
      a0 <- stats::glm(adjusted_base, data = dd, family = stats::binomial())
      a1 <- stats::glm(adjusted_full, data = dd, family = stats::binomial())
      p_lrt <- function(small, large) {
        tab <- stats::anova(small, large, test = "Chisq")
        as.numeric(tab[["Pr(>Chi)"]][2])
      }
      list(unadjusted = m1, adjusted = a1,
           p_unadjusted = p_lrt(m0, m1),
           p_adjusted = p_lrt(a0, a1),
           converged = m1$converged && a1$converged)
    }, error = function(e) {
      failures <<- c(failures, paste(v, conditionMessage(e)))
      NULL
    })
    if (is.null(result)) next
    
    summary_rows[[v]] <- data.frame(
      factor = v, scale = scale_label, reference = reference,
      n = nrow(dd), incident_cases = sum(dd$.outcome),
      incidence_pct = 100 * mean(dd$.outcome),
      p_unadjusted = result$p_unadjusted,
      p_adjusted = result$p_adjusted,
      converged = result$converged, stringsAsFactors = FALSE)
    
    for (model_name in c("unadjusted", "adjusted")) {
      tab <- summary(result[[model_name]])$coefficients
      rows <- grep("^\\.candidate", rownames(tab))
      for (j in rows) {
        beta <- tab[j, "Estimate"]
        se <- tab[j, "Std. Error"]
        label <- if (is.na(reference)) scale_label else
          paste0(sub("^\\.candidate", "", rownames(tab)[j]), " vs ", reference)
        coefficient_rows[[length(coefficient_rows) + 1L]] <- data.frame(
          factor = v, model = model_name, contrast = label,
          OR = exp(beta), CI_low = exp(beta - 1.96 * se),
          CI_high = exp(beta + 1.96 * se),
          p_coefficient = tab[j, "Pr(>|z|)"],
          scale = scale_label, stringsAsFactors = FALSE)
      }
    }
  }
  
  if (!length(summary_rows)) stop("No candidate factors could be analysed.")
  summary_df <- do.call(rbind, summary_rows)
  rownames(summary_df) <- NULL
  summary_df$FDR_unadjusted <- stats::p.adjust(summary_df$p_unadjusted, "BH")
  summary_df$FDR_adjusted <- stats::p.adjust(summary_df$p_adjusted, "BH")
  summary_df <- summary_df[order(summary_df$p_adjusted), , drop = FALSE]
  rownames(summary_df) <- NULL
  coef_df <- if (length(coefficient_rows)) do.call(rbind, coefficient_rows) else data.frame()
  
  # Compact table: an OR is shown only when the factor has one coefficient.
  summary_df$OR_adjusted_95CI <- NA_character_
  if (nrow(coef_df)) {
    adj <- coef_df[coef_df$model == "adjusted", , drop = FALSE]
    for (v in summary_df$factor) {
      one <- adj[adj$factor == v, , drop = FALSE]
      if (nrow(one) == 1L && all(is.finite(unlist(one[1, c("OR", "CI_low", "CI_high")])))) {
        summary_df$OR_adjusted_95CI[summary_df$factor == v] <- sprintf(
          "%.2f (%.2f-%.2f)", one$OR, one$CI_low, one$CI_high)
      }
    }
  }
  
  # 4. Export tables and heatmap, following perform_confounder_analysis().
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(summary_df, file.path(output_dir, "T2D_factor_association_summary.csv"),
                   row.names = FALSE)
  utils::write.csv(coef_df, file.path(output_dir, "T2D_factor_association_OR.csv"),
                   row.names = FALSE)
  
  heatmap_data <- cbind(
    "Unadjusted" = -log10(pmax(summary_df$p_unadjusted, .Machine$double.xmin)),
    "Age/sex adjusted" = -log10(pmax(summary_df$p_adjusted, .Machine$double.xmin)))
  rownames(heatmap_data) <- summary_df$factor
  heatmap_data[!is.finite(heatmap_data)] <- NA_real_
  
  if (make_plots) {
    if (!requireNamespace("pheatmap", quietly = TRUE)) {
      stop("Install the pheatmap package to create the heatmap.")
    }
    colors <- grDevices::colorRampPalette(c("white", "orange", "darkorange", "red"))(100)
    pheatmap::pheatmap(
      heatmap_data, color = colors,
      breaks = seq(0, max(1, max(heatmap_data, na.rm = TRUE)), length.out = 101),
      cluster_rows = FALSE, cluster_cols = FALSE,
      display_numbers = matrix(sprintf("%.1f", heatmap_data),
                               nrow = nrow(heatmap_data), dimnames = dimnames(heatmap_data)),
      fontsize_number = 9, fontsize = 10,
      main = "Association with 5-year incident T2D\n-log10(p-value)",
      cellwidth = 110, cellheight = 23,
      filename = file.path(output_dir, "T2D_factor_association_heatmap.png"),
      na_col = "grey90")
    
    # Optional forest plot for all estimable adjusted OR contrasts.
    if (requireNamespace("ggplot2", quietly = TRUE) && nrow(coef_df)) {
      plot_df <- coef_df[coef_df$model == "adjusted" &
                           is.finite(coef_df$OR) & is.finite(coef_df$CI_low) &
                           is.finite(coef_df$CI_high) & coef_df$CI_low > 0, , drop = FALSE]
      if (nrow(plot_df)) {
        plot_df$label <- ifelse(plot_df$scale == "category vs reference",
                                paste0(plot_df$factor, ": ", plot_df$contrast),
                                paste0(plot_df$factor, " (", plot_df$scale, ")"))
        plot_df$label <- factor(plot_df$label, levels = rev(unique(plot_df$label)))
        fig <- ggplot2::ggplot(plot_df, ggplot2::aes(x = OR, y = label)) +
          ggplot2::geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55") +
          ggplot2::geom_errorbarh(ggplot2::aes(xmin = CI_low, xmax = CI_high),
                                  height = 0.2) +
          ggplot2::geom_point(size = 2) +
          ggplot2::scale_x_log10() +
          ggplot2::labs(x = "Odds ratio (95% CI), age/sex adjusted",
                        y = NULL, title = "Candidate factors and incident T2D") +
          ggplot2::theme_minimal(base_size = 11)
        ggplot2::ggsave(file.path(output_dir, "T2D_factor_association_forest.png"),
                        fig, width = 9, height = max(5, 0.34 * nrow(plot_df) + 1.8), dpi = 250)
      }
    }
  }
  
  message("T2D analysis complete: ", nrow(summary_df), " factors; ",
          sum(data$.outcome == 1, na.rm = TRUE), " incident cases in full input.")
  if (length(failures)) warning("Skipped/failed: ", paste(failures, collapse = "; "))
  message("Interpretation: outcome=0 must denote adequate event-free 5-year follow-up.")
  return(list(summary = summary_df, coefficients = coef_df,
              heatmap = heatmap_data, warnings = failures))
}




