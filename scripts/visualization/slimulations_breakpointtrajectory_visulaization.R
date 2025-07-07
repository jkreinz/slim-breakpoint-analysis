#!/usr/bin/env Rscript
library(dplyr)
library(ggplot2)
library(viridis)
library(tidyr)
library(patchwork)

# Read summary results for all three simulation types
summary_additive <- read.csv("/Users/juliakreiner/simulation_outputs/additive_n200_summary_updated.csv")
summary_additive_n70 <- read.csv("/Users/juliakreiner/simulation_outputs/additive_n70_summary_updated.csv")

# Handle recessive summary file (check if it exists and has the right name)
recessive_file <- "/Users/juliakreiner/simulation_outputs_recessive/updated_breakpoint_summary_recessive.csv"

# Read recessive data and fix column names if needed
if (file.exists(recessive_file)) {
  summary_recessive <- read.csv(recessive_file, header = TRUE)
  
  # If no header, read again and assign column names
  if (ncol(summary_recessive) != ncol(summary_additive) ||
      !all(names(summary_recessive) == names(summary_additive))) {
    summary_recessive <- read.csv(recessive_file, header = FALSE)
    names(summary_recessive) <- names(summary_additive)
  }
} else {
  cat("Warning: Could not find recessive summary file. Proceeding with only additive data.\n")
  summary_recessive <- data.frame() # Empty dataframe
}

# Define consistent viridis colors for each simulation type
sim_colors <- c("Additive (n=200)" = "#440154FF",    # Deep purple (viridis start)
                "Additive (n=70)" = "#21908CFF",     # Teal (viridis middle)
                "Recessive (n=200)" = "#FDC100Ff")   # Green (better contrast than yellow)

# True breakpoints and selection coefficients (based on your SLiM design)
true_breakpoints <- c(50, 100)
true_s <- c(0.05, -0.025, 0.1)  # Periods 1, 2, 3

# Process first 100 iterations of each simulation type
max_iter <- 100
cat("Processing summary data for first", max_iter, "iterations...\n")

# Filter to first max_iter iterations and add simulation type labels
summary_additive_filtered <- summary_additive %>% 
  dplyr::filter(iteration <= max_iter) %>%
  dplyr::mutate(simulation_type = "Additive (n=200)")

summary_additive_n70_filtered <- summary_additive_n70 %>% 
  dplyr::filter(iteration <= max_iter) %>%
  dplyr::mutate(simulation_type = "Additive (n=70)")

# Only add recessive data if it exists
if (nrow(summary_recessive) > 0) {
  summary_recessive_filtered <- summary_recessive %>% 
    dplyr::filter(iteration <= max_iter) %>%
    dplyr::mutate(simulation_type = "Recessive (n=200)")
  
  # Combine all summary data
  all_summary <- rbind(summary_additive_filtered, summary_additive_n70_filtered, summary_recessive_filtered)
} else {
  # Only use additive data
  all_summary <- rbind(summary_additive_filtered, summary_additive_n70_filtered)
  cat("Note: Only using additive simulation data.\n")
}

# ANALYSIS 1: SELECTION COEFFICIENT PRECISION
# Updated column names: s_period_1 → s_1, etc.
s_diff_data <- all_summary %>%
  dplyr::filter(best_bic_bp == 2) %>%  # Only include iterations with 2 breakpoints
  dplyr::mutate(
    period_1_diff = s_period_1 - true_s[1],      # Changed from s_period_1 to s_1
    period_2_diff = s_period_2 - true_s[2],      # Changed from s_period_2 to s_2  
    period_3_diff = s_period_3 - true_s[3]       # Changed from s_period_3 to s_3
  ) %>%
  dplyr::select(simulation_type, iteration, period_1_diff, period_2_diff, period_3_diff) %>%
  tidyr::pivot_longer(cols = c(period_1_diff, period_2_diff, period_3_diff),
                      names_to = "period", values_to = "s_difference") %>%
  dplyr::mutate(period = dplyr::case_when(
    period == "period_1_diff" ~ "Period 1\ns=+0.05",
    period == "period_2_diff" ~ "Period 2\ns=-0.025", 
    period == "period_3_diff" ~ "Period 3\ns=+0.10"
  ))

# Filter out simulation types with insufficient data (< 5 iterations with 2 breakpoints)
valid_sim_types <- s_diff_data %>%
  dplyr::group_by(simulation_type) %>%
  dplyr::summarise(n_valid = length(unique(iteration)), .groups = 'drop') %>%
  dplyr::filter(n_valid >= 5) %>%
  dplyr::pull(simulation_type)

s_diff_data <- s_diff_data %>%
  dplyr::filter(simulation_type %in% valid_sim_types)

# ANALYSIS 2: BREAKPOINT ACCURACY  
# Updated column names: breakpoint_1 → bp_1, etc.
bp_data <- all_summary %>%
  dplyr::filter(best_bic_bp == 2) %>%  # Only iterations with 2 breakpoints detected
  dplyr::select(simulation_type, iteration, dplyr::contains("breakpoint_")) %>%
  # Handle the new column naming convention
  {
    if ("bp_1" %in% names(.)) {
      tidyr::pivot_longer(., cols = c("bp_1", "bp_2"),
                          names_to = "bp_number", values_to = "breakpoint")
    } else if ("breakpoint_1" %in% names(.)) {
      tidyr::pivot_longer(., cols = c("breakpoint_1", "breakpoint_2"),
                          names_to = "bp_number", values_to = "breakpoint")
    } else {
      # If no recognizable breakpoint columns, create empty dataframe
      data.frame(simulation_type = character(0), iteration = numeric(0),
                 bp_number = character(0), breakpoint = numeric(0))
    }
  } %>%
  dplyr::mutate(
    bp_number = as.numeric(gsub("breakpoint_|bp_", "", bp_number)),
    true_bp = ifelse(bp_number == 1, true_breakpoints[1], true_breakpoints[2]),
    bp_difference = breakpoint - true_bp,
    bp_label = paste("Breakpoint", bp_number)
  ) %>%
  dplyr::filter(!is.na(breakpoint))

# ANALYSIS 3: DETECTION POWER
# This shows ALL iterations (not just those with 2 breakpoints)
success_rates <- all_summary %>%
  dplyr::group_by(simulation_type) %>%
  dplyr::summarise(
    total_iterations = dplyr::n(),
    correct_detections = sum(best_bic_bp == 2, na.rm = TRUE),
    success_rate = correct_detections / total_iterations,
    success_rate_se = sqrt((success_rate * (1 - success_rate)) / total_iterations),
    .groups = 'drop'
  )

# CREATE COMBINED PLOT
# Panel A: Detection power (point plot with error bars)
p4 <- ggplot(success_rates, aes(x = simulation_type, y = success_rate, color = simulation_type)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = pmax(0, success_rate - 1.96*success_rate_se),
                    ymax = pmin(1, success_rate + 1.96*success_rate_se)),
                width = 0.2, size = 1) +
  geom_text(aes(label = paste0(correct_detections, "/", total_iterations)),
            vjust = 4, size = 3, color = "black") +
  scale_color_manual(values = sim_colors) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
  labs(#title = "A) Detection Power",
       x = "Simulation Type",
       y = "Success Rate\n(2 BP detected)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "none")

# Panel B: Selection coefficient precision (violin plots with means)
p2 <- ggplot(s_diff_data, aes(x = simulation_type, y = s_difference, fill = simulation_type)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.2, alpha = 0.9, outlier.shape = NA) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.7) +
  facet_grid(. ~ period, scales = "free_y") +
  scale_fill_manual(values = sim_colors) +
  coord_cartesian(ylim = c(-0.075, .05)) +
  labs(#title = "B) Selection Coefficient Precision",
       x = "",
       y = "Bias\n(Estimated s - True s)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "none",
        strip.text = element_text(size = 10))

# Panel C: Breakpoint accuracy (violin plots)
if(nrow(bp_data) > 0) {
  p3 <- ggplot(bp_data, aes(x = simulation_type, y = bp_difference, fill = simulation_type)) +
    geom_violin(alpha = 0.7, trim = FALSE) +
    geom_boxplot(width = 0.2, alpha = 0.9, outlier.shape = NA) +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.7) +
    facet_grid(. ~ bp_label, scales = "free_y") +
    scale_fill_manual(values = sim_colors) +
    labs(#title = "C) Breakpoint Accuracy",
         x = "",
         y = "Bias\n(Estimated BP - True BP)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          legend.position = "none",
          strip.text = element_text(size = 10)) +
    coord_cartesian(ylim = c(-15, 20))
} else {
  p3 <- ggplot() +
    labs() +
    theme_minimal() +
    coord_cartesian(ylim = c(-15, 20))
}

# Combine summary-based panels horizontally - Power first, then Precision and Accuracy
combined_summary_analysis <- p4 | p2 | p3
combined_summary_analysis

ggsave("/Users/juliakreiner/updated_summary_method_performance.png", combined_summary_analysis, 
       width = 18, height = 6, dpi = 300)

# Display summary statistics
cat("\nSummary Statistics:\n")
print(success_rates)

print(combined_summary_analysis)

#########################

# Function to extract sample counts per year from individual files
extract_sample_counts <- function(summary_data, output_dir, file_prefix, sim_type, max_iter = 50) {
  all_sample_counts <- data.frame()
  
  # Limit to first max_iter iterations
  iterations_to_process <- summary_data$iteration[summary_data$iteration <= max_iter]
  
  for(iter in iterations_to_process) {
    # Only process iterations where 2 breakpoints were detected
    best_bp <- summary_data$best_bic_bp[summary_data$iteration == iter]
    if(best_bp != 2) {
      cat("Skipping", sim_type, "iteration", iter, "- detected", best_bp, "breakpoints instead of 2\n")
      next
    }
    
    cat("Processing", sim_type, "iteration", iter, "for sample counts\n")
    
    # Read individual names
    inds <- read.table(paste0(output_dir, "/iteration_", iter, "/", file_prefix, iter, ".012.indv"))
    
    # Extract year information (same logic as in your script)
    year <- gsub("i","",gsub("t","",(gsub("_sample_i0","",inds$V1))))
    year2 <- as.numeric(unlist(lapply(strsplit(year,"_"), `[[`, 1)))
    
    # Count samples per year for this iteration
    year_counts <- table(year2)
    
    # Convert to data frame
    iter_counts <- data.frame(
      year = as.numeric(names(year_counts)),
      count = as.numeric(year_counts),
      iteration = iter,
      simulation_type = sim_type
    )
    
    all_sample_counts <- rbind(all_sample_counts, iter_counts)
  }
  
  return(all_sample_counts)
}

# Extract sample counts for each simulation type
cat("Extracting sample counts for additive (n=200) simulations...\n")
additive_sample_counts <- extract_sample_counts(
  summary_additive, 
  "/Users/juliakreiner/simulation_outputs", 
  "m2_allsamples_nomulti_iter", 
  "Additive (n=200)",
  max_iter = 100
)

cat("Extracting sample counts for additive (n=70) simulations...\n")
additive_n70_sample_counts <- extract_sample_counts(
  summary_additive_n70, 
  "/Users/juliakreiner/simulation_outputs", 
  "m2_allsamples_n70_nomulti_iter", 
  "Additive (n=70)",
  max_iter = 100
)

if (nrow(summary_recessive) > 0) {
  cat("Extracting sample counts for recessive (n=200) simulations...\n")
  recessive_sample_counts <- extract_sample_counts(
    summary_recessive, 
    "/Users/juliakreiner/simulation_outputs_recessive", 
    "m2_allsamples_nomulti_recessive_iter", 
    "Recessive (n=200)",
    max_iter = 100
  )
  
  # Combine all sample counts
  all_sample_counts <- rbind(
    additive_sample_counts,
    additive_n70_sample_counts,
    recessive_sample_counts
  )
} else {
  # Only additive data
  all_sample_counts <- rbind(
    additive_sample_counts,
    additive_n70_sample_counts
  )
}

# Calculate the final sampling year from our sample count data
final_sampling_year_samples <- max(all_sample_counts$year, na.rm = TRUE)
cat("Final sampling year detected from sample data:", final_sampling_year_samples, "\n")

# Calculate mean sample counts per year across iterations
sample_counts_summary <- all_sample_counts %>%
  group_by(simulation_type, year) %>%
  summarise(
    mean_count = mean(count, na.rm = TRUE),
    sd_count = sd(count, na.rm = TRUE),
    n_iterations = n(),
    .groups = 'drop'
  ) %>%
  mutate(
    se_count = sd_count / sqrt(n_iterations),
    years_before_present = final_sampling_year_samples - year  # Convert to years before present
  )


# Create stacked bar plot
p_sample_counts <- ggplot(sample_counts_summary, 
                          aes(x = years_before_present, y = mean_count, fill = simulation_type)) +
  geom_col(position = "stack", alpha = 0.8) +
  scale_fill_manual(name = "Simulation Type", values = sim_colors) +
  scale_x_reverse() +  # Reverse x-axis to match your trajectory plot
  labs(
    title = "Sample Distribution Across Time",
    x = "Years Before Present",
    y = "Mean Number of Samples"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p_sample_counts

# ==========================================
# DETAILED INDIVIDUAL FILE ANALYSIS (SLOW)
# ==========================================

cat("\n=== DETAILED INDIVIDUAL FILE ANALYSIS ===\n")

# Function to process iterations for detailed trajectory analysis
process_iterations <- function(summary_data, output_dir, file_prefix, sim_type, max_iter = 50) {
  all_data <- list()
  all_predictions <- list()
  all_breakpoints <- data.frame()
  
  # Limit to first max_iter iterations
  iterations_to_process <- summary_data$iteration[summary_data$iteration <= max_iter]
  
  for(iter in iterations_to_process) {
    # Only process iterations where 2 breakpoints were detected
    best_bp <- summary_data$best_bic_bp[summary_data$iteration == iter]
    if(best_bp != 2) {
      cat("Skipping", sim_type, "iteration", iter, "- detected", best_bp, "breakpoints instead of 2\n")
      next
    }
    
    cat("Processing", sim_type, "iteration", iter, "- correctly detected 2 breakpoints\n")
    
    # Read genotype data
    genos <- read.table(paste0(output_dir, "/iteration_", iter, "/", file_prefix, iter, ".012"))
    inds <- read.table(paste0(output_dir, "/iteration_", iter, "/", file_prefix, iter, ".012.indv"))
    
    # Extract year information
    year <- gsub("i","",gsub("t","",(gsub("_sample_i0","",inds$V1))))
    year2 <- as.numeric(unlist(lapply(strsplit(year,"_"), `[[`, 1)))
    
    # Convert to long format
    genos_long <- data.table::melt(genos[,-1], id.vars = NULL)
    genos_long$year <- rep(year2, ncol(genos)-1)
    genos_long$locus <- rep(1:(ncol(genos)-1), each = nrow(genos))
    genos_long$individual <- rep(1:nrow(genos), ncol(genos)-1)
    genos_long$value <- genos_long$value/2
    genos_long <- genos_long[!is.na(genos_long$value),]
    genos_long$iteration <- iter
    genos_long$simulation_type <- sim_type
    
    # Store raw data
    all_data[[iter]] <- genos_long
    
    # Fit base linear model
    lm_base <- glm(value ~ year, family="binomial", data=genos_long)
    
    # Fit segmented model with 2 breakpoints
    tryCatch({
      seg_model <- segmented(lm_base, seg.Z = ~year, npsi = 2)
    }, error = function(e) {
      cat("Error fitting segmented model for", sim_type, "iteration", iter, "\n")
      seg_model <<- lm_base
    })
    
    # Create prediction data
    year_seq <- seq(min(genos_long$year), max(genos_long$year), length.out = 50)
    newdata <- data.frame(year = year_seq)
    pred_data <- data.frame(
      year = year_seq,
      predicted = predict(seg_model, newdata, type = "response"),
      iteration = iter,
      simulation_type = sim_type
    )
    all_predictions[[iter]] <- pred_data
    
    # Store breakpoints if segmented model worked
    if(class(seg_model)[1] == "segmented") {
      breakpoints_est <- seg_model$psi[,"Est."]
      bp_data <- data.frame(
        iteration = iter,
        breakpoint = breakpoints_est,
        bp_number = 1:length(breakpoints_est),
        simulation_type = sim_type
      )
      all_breakpoints <- rbind(all_breakpoints, bp_data)
    }
  }
  
  return(list(
    data = all_data,
    predictions = all_predictions,
    breakpoints = all_breakpoints
  ))
}

# Process all simulation types for detailed trajectory analysis
cat("Processing additive (n=200) simulations for trajectory plot...\n")
additive_results <- process_iterations(
  summary_additive, 
  "/Users/juliakreiner/simulation_outputs", 
  "m2_allsamples_nomulti_iter", 
  "Additive (n=200)",
  max_iter = 100
)

cat("Processing additive (n=70) simulations for trajectory plot...\n")
additive_n70_results <- process_iterations(
  summary_additive_n70, 
  "/Users/juliakreiner/simulation_outputs", 
  "m2_allsamples_n70_nomulti_iter", 
  "Additive (n=70)",
  max_iter = 100
)

if (nrow(summary_recessive) > 0) {
  cat("Processing recessive (n=200) simulations for trajectory plot...\n")
  recessive_results <- process_iterations(
    summary_recessive, 
    "/Users/juliakreiner/simulation_outputs_recessive", 
    "m2_allsamples_nomulti_recessive_iter", 
    "Recessive (n=200)",
    max_iter = 100
  )
  
  # Combine all data
  combined_data_detailed <- rbind(
    do.call(rbind, additive_results$data),
    do.call(rbind, additive_n70_results$data),
    do.call(rbind, recessive_results$data)
  )
  
  combined_predictions_detailed <- rbind(
    do.call(rbind, additive_results$predictions),
    do.call(rbind, additive_n70_results$predictions),
    do.call(rbind, recessive_results$predictions)
  )
  
  combined_breakpoints_detailed <- rbind(
    additive_results$breakpoints,
    additive_n70_results$breakpoints,
    recessive_results$breakpoints
  )
} else {
  # Only additive data
  combined_data_detailed <- rbind(
    do.call(rbind, additive_results$data),
    do.call(rbind, additive_n70_results$data)
  )
  
  combined_predictions_detailed <- rbind(
    do.call(rbind, additive_results$predictions),
    do.call(rbind, additive_n70_results$predictions)
  )
  
  combined_breakpoints_detailed <- rbind(
    additive_results$breakpoints,
    additive_n70_results$breakpoints
  )
}

# Calculate mean and CI of estimated selection coefficients for detailed plot annotations
calculate_s_summary_detailed <- function(summary_data, sim_type, valid_iterations) {
  summary_data %>%
    dplyr::filter(iteration %in% valid_iterations) %>%
    dplyr::summarise(
      simulation_type = sim_type,
      period_1_mean = mean(s_period_1, na.rm = TRUE),
      period_1_ci_lower = mean(s_period_1, na.rm = TRUE) - 1.96 * mean(s_period_1_se, na.rm = TRUE),
      period_1_ci_upper = mean(s_period_1, na.rm = TRUE) + 1.96 * mean(s_period_1_se, na.rm = TRUE),
      period_2_mean = mean(s_period_2, na.rm = TRUE),
      period_2_ci_lower = mean(s_period_2, na.rm = TRUE) - 1.96 * mean(s_period_2_se, na.rm = TRUE),
      period_2_ci_upper = mean(s_period_2, na.rm = TRUE) + 1.96 * mean(s_period_2_se, na.rm = TRUE),
      period_3_mean = mean(s_period_3, na.rm = TRUE),
      period_3_ci_lower = mean(s_period_3, na.rm = TRUE) - 1.96 * mean(s_period_3_se, na.rm = TRUE),
      period_3_ci_upper = mean(s_period_3, na.rm = TRUE) + 1.96 * mean(s_period_3_se, na.rm = TRUE)
    )
}

s_summary_additive_detailed <- calculate_s_summary_detailed(
  summary_additive, 
  "Additive (n=200)", 
  unique(additive_results$breakpoints$iteration)
)

s_summary_additive_n70_detailed <- calculate_s_summary_detailed(
  summary_additive_n70, 
  "Additive (n=70)", 
  unique(additive_n70_results$breakpoints$iteration)
)

if (nrow(summary_recessive) > 0) {
  s_summary_recessive_detailed <- calculate_s_summary_detailed(
    summary_recessive, 
    "Recessive (n=200)", 
    unique(recessive_results$breakpoints$iteration)
  )
}

# Find the maximum year in your data to determine the final sampling year
final_sampling_year <- max(c(
  if(exists("combined_data_detailed")) max(combined_data_detailed$year, na.rm = TRUE) else 0,
  if(exists("combined_predictions_detailed")) max(combined_predictions_detailed$year, na.rm = TRUE) else 0
), na.rm = TRUE)

cat("Final sampling year detected as:", final_sampling_year, "\n")

# Convert data to years before present
if(exists("combined_data_detailed")) {
  combined_data_detailed <- combined_data_detailed %>%
    mutate(years_before_present = final_sampling_year - year)
}

if(exists("combined_predictions_detailed")) {
  combined_predictions_detailed <- combined_predictions_detailed %>%
    mutate(years_before_present = final_sampling_year - year)
}

# Convert true breakpoints to years before present
true_breakpoints_ybp <- final_sampling_year - true_breakpoints

# Convert breakpoint ranges to years before present
if(exists("breakpoint_ranges")) {
  breakpoint_ranges_ybp <- breakpoint_ranges %>%
    mutate(
      min_bp_ybp = final_sampling_year - max_bp,  # Note: min/max swap because we're reversing
      max_bp_ybp = final_sampling_year - min_bp
    )
}

# Create the updated plot with Years Before Present
p1_ybp <- ggplot() +
  # Fitted lines for each iteration  
  geom_line(data = combined_predictions_detailed,
            aes(x = years_before_present, y = predicted, 
                group = paste(iteration, simulation_type),
                color = simulation_type),
            alpha = 0.6, size = 0.4) +
  
  # True breakpoints (solid black lines) - converted to years before present
  geom_vline(xintercept = true_breakpoints_ybp,
             color = "black", linetype = "solid", size = 1.2) +
  
  # Breakpoint ranges as translucent rectangles - converted to years before present
  geom_rect(data = breakpoint_ranges_ybp,
            aes(xmin = min_bp_ybp, xmax = max_bp_ybp,
                ymin = -Inf, ymax = Inf,
                fill = simulation_type),
            alpha = 0.25) +
  
  # Manual viridis color scale for both color and fill
  scale_color_manual(name = "Simulation Type", values = sim_colors) +
  scale_fill_manual(name = "Simulation Type", values = sim_colors) +
  
  # Reverse x-axis so it goes from past (high values) to present (0)
  scale_x_reverse() +
  
  # Coordinate limits and theme
  coord_cartesian(ylim = c(0, 0.22)) +
  labs(
  #  title = "Allele Frequency Trajectories Over Time",
  #  subtitle = "Herbarium sampling perspective: from historical specimens to present",
    x = "Years Before Present",
    y = "Allele Frequency"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

# Display the updated plot
p1_ybp

# Display final summary statistics
cat("\nDetailed Analysis Summary Statistics (First", max_iter, "Iterations):\n")
cat("Additive (n=200) simulations: ", length(unique(additive_results$breakpoints$iteration)), " iterations with 2 breakpoints\n")
cat("Additive (n=70) simulations: ", length(unique(additive_n70_results$breakpoints$iteration)), " iterations with 2 breakpoints\n")
if (nrow(summary_recessive) > 0) {
  cat("Recessive (n=200) simulations: ", length(unique(recessive_results$breakpoints$iteration)), " iterations with 2 breakpoints\n")
}

# Save detailed trajectory plot
ggsave("/Users/juliakreiner/detailed_trajectory_plot_first50.png", p1, width = 16, height = 12, dpi = 300)

# Display both plots
print(combined_summary_analysis)
print(p1)
# Save detailed trajectory plot
ggsave("/Users/juliakreiner/detailed_trajectory_plot_first50.png", p1, width = 16, height = 12, dpi = 300)

# Display both plots
print(combined_summary_analysis)
print(p1)
cat("\nSelection Coefficient Bias Summary:\n")
s_bias_summary <- s_diff_data %>%
  group_by(simulation_type, period) %>%
  summarise(
    mean_bias = mean(s_difference, na.rm = TRUE),
    rmse = sqrt(mean(s_difference^2, na.rm = TRUE)),
    .groups = 'drop'
  )
print(s_bias_summary)

if(nrow(bp_data) > 0) {
  cat("\nBreakpoint Accuracy Summary:\n")
  bp_accuracy_summary <- bp_data %>%
    group_by(simulation_type, bp_label) %>%
    summarise(
      mean_bias = mean(bp_difference, na.rm = TRUE),
      rmse = sqrt(mean(bp_difference^2, na.rm = TRUE)),
      .groups = 'drop'
    )
  print(bp_accuracy_summary)
}

###################################

