#helper functions for data_cleaning.R
#Written by Jacob S. Berv, jacob.berv@gmail.com

# find_string_matches: Searches for occurrences of a specified substring within a vector of strings.  
# Returns a list containing:  
# - A logical vector indicating which elements contain the substring.  
# - A vector of the matching elements.  
# - A vector of the indices of the matching elements.  
find_string_matches <- function(string_vector, search_string) {
  # Check each element for the presence of search_string
  matches <- grepl(search_string, string_vector)
  
  # Check if any element contains the specified string
  any_match <- any(matches)
  
  # Prepare the result
  result <- list(
    matches = matches,
    matching_elements = character(0),  # Empty character vector if no matches
    matching_indices = integer(0)      # Empty integer vector if no matches
  )
  
  # If there are matching elements, add them to the result
  if (any_match) {
    result$matching_elements <- string_vector[matches]
    result$matching_indices <- which(matches)
  }
  
  return(result)
}

# calculate_na_percentage: Computes the percentage of NA (missing) entries in a data frame.  
# Prints the total number of non-NA values, NA values, and overall values.  
# Returns the percentage of NA values in the data frame. 
calculate_na_percentage <- function(data_frame) {
  # Check for NA values in the data frame
  na_logical <- is.na(data_frame)
  
  # Calculate the total number of NA entries
  total_na_count <- sum(na_logical)
  
  # Calculate the total number of elements in the data frame
  total_elements <- length(unlist(data_frame))
  
  # Calculate the total number of non-NA entries
  total_non_na_count <- total_elements - total_na_count
  
  # Calculate the percentage of NA entries
  percentage_na <- (total_na_count / total_elements) * 100
  
  # Print the total number of non-NA and NA values
  cat("Total non-NA values:", total_non_na_count, "\n")
  cat("Total NA values:", total_na_count, "\n")
  cat("Total values:", total_na_count+total_non_na_count, "\n")
  
  
  return(percentage_na)
}

# calculate_specimen_statistics: Estimates statistical measures of specimen counts per species.  
# Groups data by species and calculates variance, standard deviation, range, mean, and median  
# of the number of specimens per species. Returns these statistics as a list.  
calculate_specimen_statistics <- function(data) {
  # Ensure dplyr is loaded
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("dplyr package is not installed. Please install it using install.packages('dplyr').")
  }
  
  # Load dplyr for data manipulation
  library(dplyr)
  
  # Group by species and count the number of specimens per species
  specimens_per_species <- data %>%
    group_by(phylo) %>%
    summarise(specimen_count = n())
  
  # Calculate the variance, standard deviation, range, mean, and median of the specimen counts
  variance_specimens <- var(specimens_per_species$specimen_count)
  std_dev_specimens <- sd(specimens_per_species$specimen_count)
  range_specimens <- range(specimens_per_species$specimen_count)
  mean_specimens <- mean(specimens_per_species$specimen_count)
  median_specimens <- median(specimens_per_species$specimen_count)
  
  # Create a list to store and return the results
  results <- list(
    variance = variance_specimens,
    standard_deviation = std_dev_specimens,
    range_min = range_specimens[1],
    range_max = range_specimens[2],
    mean = mean_specimens,
    median = median_specimens
  )
  
  return(results)
}

# calculate_missing_percentage: Estimates the percentage of missing numerical data in a dataset,  
# excluding mass values, and associates it with specimen mass.  
# Returns a data frame with ID, species, log-transformed mass, and missing data percentage.  
# Optionally generates a violin plot of mass distribution by missing data percentage.  
calculate_missing_percentage <- function(data, plot = FALSE) {
  # Filter rows where mass is not NA
  data_with_mass <- data[!is.na(data$vertnet_mass), ]
  
  # Select numerical columns (excluding mass)
  num_cols <- sapply(data_with_mass, is.numeric)
  num_cols["vertnet_mass"] <- FALSE
  numerical_data <- data_with_mass[, num_cols]
  
  # Calculate % of missing data in numerical columns
  missing_percentage <- apply(numerical_data, 1, function(row) {
    mean(is.na(row)) * 100
  })
  
  # Create the output data frame
  data_with_mass$missing_percentage <- missing_percentage
  result <- data.frame(
    ID = data_with_mass$ID,
    phylo = data_with_mass$phylo,
    mass = log(data_with_mass$vertnet_mass),
    missing_percentage = missing_percentage
  )
  
  # Plotting, if requested
  if (plot) {
    library(vioplot)
    library(viridis)  # Assuming you're using the viridis package for color scales
    
    # Split the mass data based on unique values of missing percentage
    split_mass <- split(result$mass, result$missing_percentage)
    
    # Create the violin plot
    op <- par(mar = c(5, 4, 4, 2) + 0.1)
    
    # Round the missing percentages in the data frame
    result$rounded_missing_percentage <- round(result$missing_percentage, 0)
    
    # Split the mass data based on rounded values of missing percentage
    split_mass_rounded <- split(result$mass, result$rounded_missing_percentage)
    
    # Calculate the sample size for each group and normalize
    group_sizes <- sapply(split_mass_rounded, length)
    wex_values <- group_sizes / max(group_sizes)
    
    vioplot(split_mass_rounded, names = names(split_mass_rounded), col = viridis(length(split_mass_rounded)), wex = wex_values + 0.15, pchMed = 21, colMed = 'black', colMed2 = 'white', cex = 2)
    title(main = "Violin Plot of Mass by % of Missing Data", xlab = "% of Missing Data", ylab = "Mass")
    
    # Add text for sample sizes
    for (i in seq_along(split_mass_rounded)) {
      text(i, max(split_mass_rounded[[i]]) + 0.025, labels = group_sizes[i], cex = 0.8, pos = 3)
    }
    
    par(op)
  }
  
  return(result)
}

# calculate_r_squared_matching_columns: Computes the R-squared values for matching numeric columns  
# between two data frames by fitting linear models.  
# Returns a named vector where each element represents the R-squared value for a matching column.  
calculate_r_squared_matching_columns <- function(df1, df2) {
  # Find matching column names in both data frames
  matching_names <- intersect(names(df1), names(df2))
  
  # Initialize a vector to store the R^2 values
  r_squared_values <- numeric(length(matching_names))
  
  # Calculate R^2 for each pair of matching columns
  for (i in seq_along(matching_names)) {
    col_name <- matching_names[i]
    fit <- lm(df2[[col_name]] ~ df1[[col_name]])
    r_squared_values[i] <- summary(fit)$r.squared
  }
  
  # Create a named vector to make the output more informative
  names(r_squared_values) <- paste("R^2 for", matching_names)
  
  return(r_squared_values)
}

# calculate_correlations_matching_columns: Computes Pearson correlation coefficients for matching  
# numeric columns between two data frames.  
# Returns a named vector where each element represents the correlation coefficient for a matching column. 
calculate_correlations_matching_columns <- function(df1, df2) {
  # Find matching column names in both data frames
  matching_names <- intersect(names(df1), names(df2))
  
  # Initialize a vector to store the correlation coefficients
  correlations <- numeric(length(matching_names))
  
  # Calculate correlation for each pair of matching columns
  for (i in seq_along(matching_names)) {
    col_name <- matching_names[i]
    correlations[i] <- cor(df1[[col_name]], df2[[col_name]], use = "complete.obs", method = "pearson")
  }
  
  # Create a named vector to make the output more informative
  names(correlations) <- paste("Correlation between", matching_names)
  
  return(correlations)
}

# create_datasets: Generates three datasets by applying a mask to reconstructed data  
# based on missing values in the original data.  
# Returns a list containing:  
# - filtered_dataset: Original data with ID column.  
# - imputed_dataset: Reconstructed data with missing values retained and species column adjusted.  
# - completed_dataset: Fully reconstructed data with ID and species columns.  
create_datasets <- function(ID, original_data, reconstructed_data) {
  # Ensure ID is a data frame with a single column
  if (!is.data.frame(ID) || ncol(ID) != 1) {
    stop("ID should be a data frame with a single column.")
  }
  
  # Creating the mask
  mask <- is.na(original_data)
  
  # Apply mask to reconstructed data to create tmp
  tmp <- reconstructed_data
  tmp[mask] <- NA
  
  # Check for equivalence with tolerance (optional)
  equivalence_result <- check_full_equivalence_with_tolerance(original_data, tmp)
  print(equivalence_result)
  
  names <- ID
  
  # Creating filtered_dataset
  filtered_dataset <- cbind(ID = names[,1], original_data)
  
  # Inverse mask to mask the measured data
  tmp <- reconstructed_data
  tmp[!mask] <- NA
  imputed_dataset <- cbind(ID = names[,1], tmp)
  colnames(imputed_dataset)[2]<-'species'
  imputed_dataset[,2] <- filtered_dataset$species
  
  # Creating completed_dataset
  completed_dataset <- cbind(ID = names[,1], reconstructed_data)
  colnames(completed_dataset)[2]<-'species'
  
  # Return a list of the datasets
  list(
    filtered_dataset = filtered_dataset,
    imputed_dataset = imputed_dataset,
    completed_dataset = completed_dataset
    #equivalence_result = equivalence_result
  )
}

# check_full_equivalence_with_tolerance: Compares two data frames for full equivalence,  
# allowing for small numerical differences within a specified tolerance.  
# Returns TRUE if the data frames are equivalent, otherwise prints the mismatch and returns FALSE.  
check_full_equivalence_with_tolerance <- function(df1, df2, tolerance = 1e-8) {
  if (!all(dim(df1) == dim(df2))) {
    print("Dimensions of data frames do not match.")
    return(FALSE)
  }
  
  for (row in 1:nrow(df1)) {
    for (col in 1:ncol(df1)) {
      val1 <- df1[row, col]
      val2 <- df2[row, col]
      
      # Check for NA and NaN equivalence
      if (!((is.na(val1) && is.na(val2)) || (is.nan(val1) && is.nan(val2)))) {
        # If not both NA or NaN, check for value mismatch with tolerance for numeric values
        if (is.numeric(val1) && is.numeric(val2)) {
          if (abs(val1 - val2) > tolerance) {
            print(paste("Mismatch found at row", row, "column", col, 
                        ": values", val1, "and", val2, 
                        "differ by more than tolerance", tolerance))
            return(FALSE)
          }
        } else if (val1 != val2) {
          print(paste("Mismatch found at row", row, "column", col, 
                      ": values", val1, "and", val2, "are not equal"))
          return(FALSE)
        }
      }
    }
  }
  
  print("Data frames are equivalent within specified tolerance.")
  return(TRUE)
}

# construct_CI_df: Constructs confidence intervals, standard errors, variances,  
# and reconstructed values data frames from given reconstructed values and variances.  
# Returns a list containing:  
# - Recon: Data frame of reconstructed values.  
# - CI: Data frame with lower and upper 95% confidence intervals.  
# - SE: Data frame of standard errors.  
# - Var: Data frame of variances.  
construct_CI_df <- function(recon_values, variances) {
  # Calculate the standard errors (square root of variances)
  se = sqrt(variances)
  
  # Extract species names from the row names of recon_values
  species <- rownames(recon_values)
  
  # Initialize the data frames with the same row names and number of rows as recon_values
  n_rows <- nrow(recon_values)
  n_cols <- ncol(recon_values)
  
  ci_df <- data.frame(Species = species, matrix(nrow = n_rows, ncol = 2 * n_cols))
  se_df <- data.frame(Species = species, matrix(nrow = n_rows, ncol = n_cols))
  var_df <- data.frame(Species = species, matrix(nrow = n_rows, ncol = n_cols))
  recon_df <- data.frame(Species = species, matrix(nrow = n_rows, ncol = n_cols))
  
  # Assign column names and populate the data frames
  for (i in 1:n_cols) {
    # Calculate CI values
    lower_ci <- recon_values[, i] - (se[, i] * 1.96)
    upper_ci <- recon_values[, i] + (se[, i] * 1.96)
    
    # Assign CI values to ci_df
    ci_df[,(2 * i)] <- lower_ci
    ci_df[,(2 * i + 1)] <- upper_ci
    
    # Assign column names for CI
    colnames(ci_df)[(2 * i)] <- paste0(colnames(recon_values)[i], "_LowerCI")
    colnames(ci_df)[(2 * i + 1)] <- paste0(colnames(recon_values)[i], "_UpperCI")
    
    # Populate other data frames
    se_df[, i + 1] <- se[, i]
    var_df[, i + 1] <- variances[, i]
    recon_df[, i + 1] <- recon_values[, i]
    
    # Assign column names for other data frames
    colnames(se_df)[i + 1] <- paste0(colnames(recon_values)[i], "_SE")
    colnames(var_df)[i + 1] <- paste0(colnames(recon_values)[i], "_Var")
    colnames(recon_df)[i + 1] <- colnames(recon_values)[i]
  }
  
  # Return a list containing the CI, SE, Variance, and Reconstructed Values data frames
  return(list(Recon = recon_df, CI = ci_df, SE = se_df, Var = var_df))
}

# mask_and_track_percentage: Applies random masking to a specified percentage of non-missing values  
# in each column of a dataset, tracking masked positions and sample sizes.  
# Runs for a specified number of repetitions, excluding specified columns from masking.  
# Returns a list where each element contains:  
# - masked_data: Data frame with randomly masked values.  
# - masked_positions: Boolean matrix indicating masked positions.  
# - sample_size: Vector tracking the number of masked values per column.  
mask_and_track_percentage <- function(data, mask_percentage, repetitions, exclude_cols = NULL) {
  # Initialize a list to store results for each repetition
  results_list <- vector("list", repetitions)
  
  # Loop for each repetition
  for (rep in 1:repetitions) {
    cat("Repetition", rep, ":\n")
    
    # Initialize dataframes for masking, tracking, and storing sample sizes
    masked_data <- data
    masked_positions <- data.frame(matrix(FALSE, nrow = nrow(data), ncol = ncol(data)))
    sample_sizes <- numeric(ncol(data))
    names(masked_positions) <- names(data)
    names(sample_sizes) <- names(data)
    
    # Loop through each column, excluding specified columns
    for (col in setdiff(names(data), exclude_cols)) {
      non_missing_indices <- which(!is.na(data[[col]]))
      num_non_missing = length(non_missing_indices)
      
      # Calculate the number of values to mask based on the percentage
      num_to_mask = ceiling(num_non_missing * mask_percentage / 100)
      
      if (num_non_missing >= num_to_mask) {
        indices_to_mask <- sample(non_missing_indices, num_to_mask, replace = FALSE)
        masked_data[[col]][indices_to_mask] <- NA
        masked_positions[[col]][indices_to_mask] <- TRUE
        sample_sizes[col] <- length(indices_to_mask)
        
        # Print details of masking for this column
        actual_mask_percentage <- (length(indices_to_mask) / num_non_missing) * 100
        cat("  Column", col, ": Intended to mask", mask_percentage, "%, actually masked", actual_mask_percentage, "% of non-missing data (", length(indices_to_mask), "out of", num_non_missing, "cells)\n")
      } else {
        cat("  Column", col, ": Not enough non-missing data to mask\n")
      }
    }
    
    # Store the results of this repetition
    results_list[[rep]] <- list(masked_data = masked_data, masked_positions = masked_positions, sample_size = sample_sizes)
    cat("\n")
  }
  
  return(results_list)
}

# compute_metrics: Evaluates imputation performance by computing RMSE and percentage bias  
# for each numeric column in the dataset, using masked positions for comparison.  
# Returns a list containing:  
# - rmse: Average RMSE per column across imputation repetitions.  
# - se_rmse: Standard error of RMSE.  
# - ci_rmse: 95% confidence interval for RMSE.  
# - p_bias: Average percentage bias per column.  
# - se_p_bias: Standard error of p-bias.  
# - ci_p_bias: 95% confidence interval for p-bias. 
compute_metrics <- function(original_data, imputed_data_list, masked_positions_list) {
  # Exclude the first column (species names) and check if remaining columns are numeric
  original_data_numeric <- original_data[-1]
  if (!all(sapply(original_data_numeric, is.numeric))) {
    stop("Not all columns (excluding the first one) are numeric in the original data.")
  }
  
  # Initialize variables for cumulative RMSE and p-bias, and their squared differences
  cum_rmse_values <- setNames(numeric(ncol(original_data_numeric)), names(original_data_numeric))
  cum_p_bias_values <- setNames(numeric(ncol(original_data_numeric)), names(original_data_numeric))
  rmse_squared_diffs <- vector("list", ncol(original_data_numeric))
  p_bias_squared_diffs <- vector("list", ncol(original_data_numeric))
  names(rmse_squared_diffs) <- names(original_data_numeric)
  names(p_bias_squared_diffs) <- names(original_data_numeric)
  
  repetitions <- length(imputed_data_list)
  
  # Loop through each repetition of imputed data
  for (rep in 1:repetitions) {
    imputed_data <- imputed_data_list[[rep]][-1]
    masked_positions <- masked_positions_list[[rep]][-1]
    
    # Check if column names match across original, imputed, and masked datasets
    if (!identical(names(original_data_numeric), names(imputed_data)) || !identical(names(original_data_numeric), names(masked_positions))) {
      stop("Column names do not match across the datasets.")
    }
    
    # Loop through each column to compute metrics
    for (col in names(original_data_numeric)) {
      original_col <- original_data_numeric[[col]]
      imputed_col <- imputed_data[[col]]
      masked_col <- masked_positions[[col]]
      
      # Skip column if no positions are masked or if any NA values are found
      if (sum(masked_col) == 0 || any(is.na(original_col[masked_col])) || any(is.na(imputed_col[masked_col]))) {
        next
      }
      
      # Compute RMSE and p-bias for the current column
      rmse_val <- sqrt(mean((original_col[masked_col] - imputed_col[masked_col])^2))
      p_bias_val <- sum((imputed_col[masked_col] - original_col[masked_col]) / original_col[masked_col]) / sum(masked_col) * 100
      
      # Accumulate the sum of RMSE and p-bias values
      cum_rmse_values[col] <- cum_rmse_values[col] + rmse_val
      cum_p_bias_values[col] <- cum_p_bias_values[col] + p_bias_val
      
      # Accumulate squared differences for RMSE and p-bias for standard error calculation
      rmse_squared_diffs[[col]] <- c(rmse_squared_diffs[[col]], (rmse_val - (cum_rmse_values[col] / rep))^2)
      p_bias_squared_diffs[[col]] <- c(p_bias_squared_diffs[[col]], (p_bias_val - (cum_p_bias_values[col] / rep))^2)
    }
  }
  
  # Calculate average RMSE and p-bias, their standard errors, and 95% confidence intervals
  avg_rmse_values <- cum_rmse_values / repetitions
  avg_p_bias_values <- cum_p_bias_values / repetitions
  se_rmse_values <- sapply(rmse_squared_diffs, function(x) sqrt(sum(x) / (length(x) - 1)) / sqrt(length(x)))
  se_p_bias_values <- sapply(p_bias_squared_diffs, function(x) sqrt(sum(x) / (length(x) - 1)) / sqrt(length(x)))
  ci_rmse_values <- 1.96 * se_rmse_values
  ci_p_bias_values <- 1.96 * se_p_bias_values
  
  return(list(rmse = avg_rmse_values, se_rmse = se_rmse_values, ci_rmse = ci_rmse_values,
              p_bias = avg_p_bias_values, se_p_bias = se_p_bias_values, ci_p_bias = ci_p_bias_values))
}

# plot_metrics: Visualizes RMSE and percentage bias (p-bias) with error bars,  
# using either standard error or 95% confidence intervals.  
# Generates two plots:  
# - RMSE plot with error bars.  
# - p-bias plot with error bars. 
plot_metrics <- function(metrics, error_type = "se") {
  num_traits <- length(metrics$rmse)
  par(mfrow = c(2, 1), mar = c(7, 4, 4, 2))
  
  # Select error type (standard error or confidence interval)
  error_rmse <- if (error_type == "se") metrics$se_rmse else metrics$ci_rmse
  error_p_bias <- if (error_type == "se") metrics$se_p_bias else metrics$ci_p_bias
  
  # Titles based on error type
  rmse_title <- if (error_type == "se") "RMSE with Standard Error" else "RMSE with 95% Confidence Interval"
  p_bias_title <- if (error_type == "se") "p-bias with Standard Error" else "p-bias with 95% Confidence Interval"
  
  # Determine y-axis limits dynamically for RMSE plot
  rmse_upper_limit <- max(metrics$rmse + error_rmse, na.rm = TRUE)
  rmse_lower_limit <- min(metrics$rmse - error_rmse, na.rm = TRUE)
  
  # RMSE plot
  #lot(1:num_traits, metrics$rmse, xaxt = 'n', pch = 19, ylim = c(rmse_lower_limit, rmse_upper_limit), xlab = "", ylab = "Mean RMSE log(mm)", main = rmse_title)
  plot(1:num_traits, metrics$rmse, xaxt = 'n', pch = 19, ylim = c(0, 0.2), xlab = "", ylab = "Mean RMSE log(mm)", main = rmse_title)
  
  axis(1, at = 1:num_traits, labels = names(metrics$rmse), las = 2, cex.axis=0.75)
  for (i in 1:num_traits) {
    lines(c(i, i), c(metrics$rmse[i] - error_rmse[i], metrics$rmse[i] + error_rmse[i]))
    lines(c(i - 0.1, i + 0.1), c(metrics$rmse[i] + error_rmse[i], metrics$rmse[i] + error_rmse[i]))
    lines(c(i - 0.1, i + 0.1), c(metrics$rmse[i] - error_rmse[i], metrics$rmse[i] - error_rmse[i]))
  }
  
  # Determine y-axis limits dynamically for p-bias plot
  p_bias_upper_limit <- max(metrics$p_bias + error_p_bias, na.rm = TRUE)
  p_bias_lower_limit <- min(metrics$p_bias - error_p_bias, na.rm = TRUE)
  
  # p-bias plot
  #plot(1:num_traits, metrics$p_bias, xaxt = 'n', pch = 19, ylim = c(p_bias_lower_limit, p_bias_upper_limit), xlab = "", ylab = "Mean p-bias (%)", main = p_bias_title)
  plot(1:num_traits, metrics$p_bias, xaxt = 'n', pch = 19, ylim = c(-5, 5), xlab = "", ylab = "Mean p-bias (%)", main = p_bias_title)
  axis(1, at = 1:num_traits, labels = names(metrics$p_bias), las = 2, cex.axis=0.75)
  for (i in 1:num_traits) {
    lines(c(i, i), c(metrics$p_bias[i] - error_p_bias[i], metrics$p_bias[i] + error_p_bias[i]))
    lines(c(i - 0.1, i + 0.1), c(metrics$p_bias[i] + error_p_bias[i], metrics$p_bias[i] + error_p_bias[i]))
    lines(c(i - 0.1, i + 0.1), c(metrics$p_bias[i] - error_p_bias[i], metrics$p_bias[i] - error_p_bias[i]))
  }
}

# plot_rmse_and_pbias: Generates two plots to visualize RMSE and percentage bias (p-bias)  
# across different masking percentages for multiple traits.  
# Includes legends to differentiate traits and allows custom x-axis limits.  
plot_rmse_and_pbias <- function(rmse_dat, pbias_dat, x_values, xlim_rmse = NULL, xlim_pbias = NULL) {
  library(viridis) # Load the viridis package
  
  # Check that the number of rows in the matrices matches the length of x_values
  if (nrow(rmse_dat) != length(x_values)) {
    stop("Number of rows in rmse_dat must be equal to the length of x_values")
  }
  if (nrow(pbias_dat) != length(x_values)) {
    stop("Number of rows in pbias_dat must be equal to the length of x_values")
  }
  
  # Generate a set of colors and point styles
  num_traits <- ncol(rmse_dat)
  set.seed(123)  # Set a seed for reproducibility
  colors <- sample(viridis(num_traits), num_traits)
  pch_styles <- sample(1:25, num_traits, replace = TRUE)
  
  # Define layout matrix for the plots and legends
  layout_matrix <- rbind(c(1, 2), c(3, 4))
  
  # Define the relative widths of the plot and legend areas
  layout_widths <- c(5, 1) # The plot area is 5 times the width of the legend area
  
  # Define common margins for both plots
  plot_margins <- c(5, 4, 4, 2) # Bottom, left, top, right margins
  
  # Setup the layout for the plots and legends
  layout(layout_matrix, widths = layout_widths)
  
  # Set the x-axis limits if not specified
  if (is.null(xlim_rmse)) {
    xlim_rmse <- range(x_values)
  }
  if (is.null(xlim_pbias)) {
    xlim_pbias <- range(x_values)
  }
  
  # RMSE plot in the first cell
  par(mar = plot_margins)
  plot(x_values, rmse_dat[, 1], type = 'b', pch = pch_styles[1], col = colors[1],
       ylim = c(0, max(rmse_dat)),
       xlim = xlim_rmse,
       xaxt = 'n', yaxt = 'n',
       xlab = "% Masked", ylab = "RMSE log(mm)", 
       main = "Mean RMSE (10 replicates per trait)", bty = 'n')
  for (i in 2:num_traits) {
    lines(x_values, rmse_dat[, i], type = 'b', pch = pch_styles[i], col = colors[i])
  }
  axis(1, at = x_values, labels = x_values)
  axis(2, las = 2, at = seq(0, max(rmse_dat), length.out = 9), labels = round(seq(0, max(rmse_dat), length.out = 9), 2))
  
  # Add RMSE Legend in the second cell
  par(mar = c(plot_margins[1], 0, plot_margins[3], plot_margins[4])) # Adjust left margin to 0 for legend
  plot.new()
  legend("topright", inset = c(0.30, 0), legend = colnames(rmse_dat), col = colors, 
         pch = pch_styles, lty = 1, cex = 0.8, xpd = NA, bty = 'n')
  
  # p-bias plot in the third cell
  par(mar = plot_margins)
  plot(x_values, pbias_dat[, 1], type = 'b', pch = pch_styles[1], col = colors[1],
       ylim = c(-1, 1),
       xlim = xlim_pbias,
       xaxt = 'n', yaxt = 'n',
       xlab = "% Masked", ylab = "p-bias (%)", 
       main = "Mean p-bias (10 replicates per trait)", bty = 'n')
  for (i in 2:num_traits) {
    lines(x_values, pbias_dat[, i], type = 'b', pch = pch_styles[i], col = colors[i])
  }
  axis(1, at = x_values, labels = x_values)
  axis(2, las = 2, at = seq(-1, 1, length.out = 9), labels = seq(-100, 100, length.out = 9))
  
  # Add p-bias Legend in the fourth cell
  par(mar = c(plot_margins[1], 0, plot_margins[3], plot_margins[4]))
  plot.new()
  legend("topright", inset = c(0.30, 0), legend = colnames(pbias_dat), col = colors, 
         pch = pch_styles, lty = 1, cex = 0.8, xpd = NA, bty = 'n')
  
  # Reset to default par settings
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
}