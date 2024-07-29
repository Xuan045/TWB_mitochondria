library(stats)

perform_phewas_analysis <- function(data, phenotypes, genotypes, covariates, min_cases = 20, age_sex_cov = FALSE) {
  results <- list()
  
  for (pheno in phenotypes) {
    subset_data <- na.omit(data[c(pheno, covariates)])
    
    # Check each covariate in the subset data
    # Some phenotypes may require excluding specific covariates due to lack of variability
    valid_covariates <- covariates
    for (cov in covariates) {
      if (length(unique(subset_data[[cov]])) < 2) {
        valid_covariates <- valid_covariates[valid_covariates != cov]
      }
    }
    
    for (geno in genotypes) {
      # Create the formula for the model
      formula_str <- paste(pheno, "~", geno)
      # Add valid covariates
      if (length(valid_covariates) > 0) {
        formula_str <- paste(formula_str, "+", paste(valid_covariates, collapse=" + "))
      }
      # Check if 'age' and 'sex' are in valid_covariates and append interaction terms if so
      if ('AGE' %in% valid_covariates && 'SEX' %in% valid_covariates && age_sex_cov == TRUE) {
        formula_str <- paste(formula_str, "+ AGE * SEX")
      }
      # Convert string to formula
      formula <- as.formula(formula_str)
      
      # Initialize note for recording any modifications
      note <- ""
      dropped_covs <- covariates[!covariates %in% valid_covariates]
      if (length(dropped_covs) > 0) {
        note <- paste0("[Note: Column(s) dropped due to lack of variability in subset: ", paste(dropped_covs, collapse=", "), "]")
      }
      
      # Model fitting
      if (class(data[[pheno]]) %in% c("factor", "logical")) {
        model_type <- "logistic"
        model <- glm(formula, data = data, family = binomial())
        n_cases <- sum(data[[pheno]] == 1, na.rm = TRUE)
        n_controls <- sum(data[[pheno]] == 0, na.rm = TRUE)
        
        # Check if the number of cases is below the threshold
        if (n_cases < min_cases) {
          note <- paste(note, "Not enough cases (n < ", min_cases, ")")
        }
      } else {
        model_type <- "linear"
        model <- glm(formula, data = data)
        n_cases <- NA
        n_controls <- NA
      }
      
      # Check if the model converged
      if (!model$converged) {
        next
      }
      
      # Results processing
      n_total <- nrow(na.omit(data[c(pheno, geno, valid_covariates)]))
      
      if (!is.null(model) && (geno %in% rownames(summary(model)$coefficients))) {
        summary_model <- summary(model)
        coef_info <- summary_model$coefficients[geno, ]
        if (model_type == "logistic") {
          # Calculate OR for logistic regression
          or <- exp(coef_info["Estimate"])
          or_ci_lower <- exp(coef_info["Estimate"] - 1.96 * coef_info["Std. Error"])
          or_ci_upper <- exp(coef_info["Estimate"] + 1.96 * coef_info["Std. Error"])
        } else {
          # Set OR to NA for linear regression
          or <- NA
          or_ci_lower <- NA
          or_ci_upper <- NA
        }
        
        # Create result data frame
        results[[paste(pheno, geno, sep="_")]] <- data.frame(
          phenotype = pheno,
          snp = geno,
          Beta = coef_info["Estimate"],
          SE = coef_info["Std. Error"],
          Statistic = if (model_type == "linear") {
            coef_info["t value"]
          } else {
            coef_info["z value"]
          },
          p = if (model_type == "linear") {
            coef_info["Pr(>|t|)"]
          } else {
            coef_info["Pr(>|z|)"]
          },
          OR = or,
          OR_CI_lower = or_ci_lower,
          OR_CI_upper = or_ci_upper,
          type = model_type,
          n_total = n_total,
          n_cases = n_cases,
          n_controls = n_controls,
          formula = formula_str,
          note = note,
          stringsAsFactors = FALSE
        )
      } else {
        next  # Skip to the next iteration if the coefficient for geno is not found
      }
    }
  }
  
  results_df <- do.call(rbind, results)
  results_df$FDR <- p.adjust(results_df$p, method = "BH")  # Benjamini-Hochberg method
  return(results_df)
}
