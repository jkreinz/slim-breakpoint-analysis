#!/usr/bin/env Rscript
# Streamlined VCF Analysis Script - Updated for new directory structure
# Produces only trajectory plots and bias plots

library(vcfR)
library(dplyr)
library(ggplot2)
library(purrr)
library(tidyr)

# Updated color scheme for three simulation types
sim_colors <- c("Additive (n=15)" = "#21908CFF",     # Teal (viridis middle)
                "Additive (n=40)" = "#440154FF",     # Deep purple (viridis start)
                "Recessive (n=40)" = "#FDC100FF")    # Yellow-orange

# Define true selection parameters
TRUE_SELECTION <- data.frame(
  shift = c("shift1", "shift2", "shift3"),
  s_before = c(0.0, 0.05, -0.025),
  s_after = c(0.05, -0.025, 0.1),
  generations_between = c(4, 4, 4)
)

# Function to read and process a single VCF file
process_vcf <- function(vcf_path) {
  vcf <- read.vcfR(vcf_path, verbose = FALSE)
  gt_matrix <- extract.gt(vcf, element = "GT")
  fix_data <- as.data.frame(getFIX(vcf))
  
  allele_freqs <- apply(gt_matrix, 1, function(gt_row) {
    allele_counts <- sapply(strsplit(gt_row, "[|/]"), function(x) sum(as.numeric(x)))
    total_alleles <- length(allele_counts) * 2
    sum(allele_counts) / total_alleles
  })
  
  results <- data.frame(
    mutation_id = paste0(fix_data$CHROM, "_", fix_data$POS, "_", fix_data$REF, "_", fix_data$ALT),
    allele_freq = allele_freqs,
    sample_size = ncol(gt_matrix),
    file = basename(vcf_path)
  )
  
  return(results %>% distinct(mutation_id, .keep_all = TRUE))
}

# Updated function to process all VCF files for a single iteration
process_iteration <- function(iteration_dir, dominance_type = "additive_n15") {
  iteration_num <- basename(iteration_dir)
  # Only look for files that start with m2_ and end with .vcf
  all_m2_files <- list.files(iteration_dir, pattern = "^m2_.*\\.vcf$", full.names = TRUE)
  
  vcf_files <- character(0)
  
  if (dominance_type == "additive_n15") {
    # Look for m2_ files ending with _sample_n15_h5.vcf (including m2_last_sample_n15_h5.vcf)
    vcf_files <- all_m2_files[grepl("_sample_n15_h5\\.vcf$", all_m2_files)]
  } else if (dominance_type == "additive_n40") {
    # Look for m2_ files ending with _sample_n40_h5.vcf (including m2_last_sample_n40_h5.vcf)
    vcf_files <- all_m2_files[grepl("_sample_n40_h5\\.vcf$", all_m2_files)]
  } else if (dominance_type == "recessive_n40") {
    # Look for m2_ files ending with _sample_n40_h015.vcf (including m2_last_sample_n40_h015.vcf)
    vcf_files <- all_m2_files[grepl("_sample_n40_h015\\.vcf$", all_m2_files)]
  }
  
  if (length(vcf_files) == 0) {
    cat("    No", dominance_type, "files found in", iteration_dir, "\n")
    return(NULL)
  }
  
  cat("    Found", length(vcf_files), dominance_type, "files\n")
  
  all_data <- map_dfr(vcf_files, process_vcf)
  
  # Extract timepoint information from filename
  if (dominance_type == "additive_n15") {
    # Handle both regular and last files
    all_data$timepoint <- ifelse(
      grepl("m2_last_sample_n15_h5\\.vcf$", all_data$file),
      "last",
      # Extract everything between m2_ and _sample_n15_h5.vcf
      gsub("^m2_(.+)_sample_n15_h5\\.vcf$", "\\1", all_data$file)
    )
  } else if (dominance_type == "additive_n40") {
    all_data$timepoint <- ifelse(
      grepl("m2_last_sample_n40_h5\\.vcf$", all_data$file),
      "last", 
      # Extract everything between m2_ and _sample_n40_h5.vcf
      gsub("^m2_(.+)_sample_n40_h5\\.vcf$", "\\1", all_data$file)
    )
  } else if (dominance_type == "recessive_n40") {
    all_data$timepoint <- ifelse(
      grepl("m2_last_sample_n40_h015\\.vcf$", all_data$file),
      "last",
      # Extract everything between m2_ and _sample_n40_h015.vcf
      gsub("^m2_(.+)_sample_n40_h015\\.vcf$", "\\1", all_data$file)
    )
  }
  
  all_data$iteration <- iteration_num
  all_data$dominance_type <- dominance_type
  
  return(all_data)
}

# Function to calculate instantaneous selection coefficients
calculate_instantaneous_s <- function(freq_data) {
  before_after_pairs <- list(
    shift1 = c("shift1_before", "shift1_after"),
    shift2 = c("shift2_before", "shift2_after"),
    shift3 = c("shift3_before", "shift3_after")
  )
  
  selection_results <- map_dfr(names(before_after_pairs), function(shift_name) {
    timepoints <- before_after_pairs[[shift_name]]
    
    before_data <- freq_data %>% 
      filter(timepoint == timepoints[1]) %>%
      distinct(mutation_id, iteration, dominance_type, .keep_all = TRUE) %>%
      dplyr::select(mutation_id, iteration, dominance_type, sample_size, p_before = allele_freq)
    
    after_data <- freq_data %>% 
      filter(timepoint == timepoints[2]) %>%
      distinct(mutation_id, iteration, dominance_type, .keep_all = TRUE) %>%
      dplyr::select(mutation_id, iteration, dominance_type, p_after = allele_freq)
    
    merged <- inner_join(before_data, after_data, by = c("mutation_id", "iteration", "dominance_type"))
    
    merged %>%
      mutate(
        shift = shift_name,
        delta_p = p_after - p_before,
        q_before = 1 - p_before,
        pq_before = p_before * q_before,
        s_instantaneous = ifelse(pq_before > 0.001, delta_p / pq_before, NA),
        valid_mutation = p_before > 0.05 & p_before < 0.95 & !is.na(s_instantaneous)
      ) %>%
      filter(valid_mutation)
  })
  
  return(selection_results)
}

# Function to add true selection values
add_true_selection <- function(selection_data) {
  selection_data %>%
    left_join(TRUE_SELECTION, by = "shift") %>%
    mutate(
      s_true_average = (s_before + s_after) / 2
    )
}

# Function to create summary statistics
summarize_selection <- function(selection_data) {
  summary_stats <- selection_data %>%
    group_by(shift, iteration, dominance_type) %>%
    summarise(
      mean_s_instantaneous = mean(s_instantaneous, na.rm = TRUE),
      s_true_average = first(s_true_average),
      sample_size = first(sample_size),
      .groups = "drop"
    )
  
  return(summary_stats)
}

# Updated main analysis function - now takes base directory as parameter
run_analysis <- function(base_dirs = NULL, max_iterations = NULL) {
  
  # If no base directories provided, try to find them automatically
  if (is.null(base_dirs)) {
    # Look for directories with "instants" in the name
    possible_dirs <- list.dirs(".", recursive = FALSE)
    base_dirs <- possible_dirs[grepl("instants", basename(possible_dirs))]
    
    if (length(base_dirs) == 0) {
      # Fallback to old naming
      base_dirs <- "instants_simulation_outputs"
    }
    
    cat("Auto-detected base directories:", paste(base_dirs, collapse = ", "), "\n")
  }
  
  # Process all base directories
  all_freq_data <- data.frame()
  
  for (base_dir in base_dirs) {
    if (!dir.exists(base_dir)) {
      cat("Warning: Directory", base_dir, "does not exist. Skipping.\n")
      next
    }
    
    cat("Processing directory:", base_dir, "\n")
    
    iteration_dirs <- list.dirs(base_dir, recursive = FALSE)
    if (!is.null(max_iterations)) {
      iteration_dirs <- iteration_dirs[1:min(max_iterations, length(iteration_dirs))]
    }
    
    cat("Processing", length(iteration_dirs), "iterations for all three simulation types\n")
    
    # Process all three types with progress reporting
    cat("Processing additive (n=15) files...\n")
    additive_n15_data <- map_dfr(iteration_dirs, function(dir) {
      iter_num <- basename(dir)
      cat("  Iteration", iter_num, "(additive n=15)\n")
      # Debug: show what m2_ files exist
      all_m2_files <- list.files(dir, pattern = "^m2_.*\\.vcf$", full.names = FALSE)
      if (length(all_m2_files) > 0) {
        cat("    All m2_ files:", paste(head(all_m2_files, 3), collapse = ", "), 
            ifelse(length(all_m2_files) > 3, "...", ""), "\n")
      }
      process_iteration(dir, "additive_n15")
    })
    
    cat("Processing additive (n=40) files...\n")
    additive_n40_data <- map_dfr(iteration_dirs, function(dir) {
      iter_num <- basename(dir)
      cat("  Iteration", iter_num, "(additive n=40)\n")
      process_iteration(dir, "additive_n40")
    })
    
    cat("Processing recessive (n=40) files...\n")
    recessive_n40_data <- map_dfr(iteration_dirs, function(dir) {
      iter_num <- basename(dir)
      cat("  Iteration", iter_num, "(recessive n=40)\n")
      process_iteration(dir, "recessive_n40")
    })
    
    # Combine data from this base directory
    dir_data <- rbind(additive_n15_data, additive_n40_data, recessive_n40_data)
    all_freq_data <- rbind(all_freq_data, dir_data)
    
    cat("Data summary for", base_dir, ":\n")
    cat("- Additive (n=15):", nrow(additive_n15_data), "rows\n")
    cat("- Additive (n=40):", nrow(additive_n40_data), "rows\n") 
    cat("- Recessive (n=40):", nrow(recessive_n40_data), "rows\n")
  }
  
  cat("Combined data summary:\n")
  cat("- Total:", nrow(all_freq_data), "rows\n")
  
  cat("Calculating selection coefficients...\n")
  selection_data <- calculate_instantaneous_s(all_freq_data)
  cat("Selection data rows:", nrow(selection_data), "\n")
  
  selection_data <- add_true_selection(selection_data)
  cat("After adding true selection:", nrow(selection_data), "\n")
  
  summaries <- summarize_selection(selection_data)
  cat("Summary data rows:", nrow(summaries), "\n")
  
  # Debug: check what timepoints we actually have
  cat("Unique timepoints in data:", paste(unique(all_freq_data$timepoint), collapse = ", "), "\n")
  cat("Timepoint counts:\n")
  print(table(all_freq_data$timepoint))
  
  cat("Analysis complete!\n")
  
  return(list(
    raw_frequencies = all_freq_data,
    selection_coefficients = selection_data,
    summaries = summaries
  ))
}

# Updated function to create bias plot
create_bias_plot <- function(results) {
  bias_data <- results$summaries %>%
    mutate(
      simulation_type = case_when(
        dominance_type == "additive_n15" ~ "Additive (n=15)",
        dominance_type == "additive_n40" ~ "Additive (n=40)",
        dominance_type == "recessive_n40" ~ "Recessive (n=40)",
        TRUE ~ "Unknown"
      ),
      # Set the factor levels to control ordering
      simulation_type = factor(simulation_type, levels = c("Additive (n=40)", "Additive (n=15)", "Recessive (n=40)")),
      bias = mean_s_instantaneous - s_true_average,
      shift_label = case_when(
        shift == "shift1" ~ "Around BP 1: S 0 → 0.05\n(true avg S: 0.025)",
        shift == "shift2" ~ "Around BP 2: S +0.05 → -0.025\n(true avg S: 0.0125)",
        shift == "shift3" ~ "Around BP 3: S -0.025 → 0.1\n(true avg S: 0.0375)"
      )
    ) %>%
    filter(simulation_type %in% names(sim_colors), !is.na(bias))
  
  # Debug: check if we have data
  cat("Bias plot data summary:\n")
  cat("- Total rows:", nrow(bias_data), "\n")
  cat("- Simulation types:", paste(unique(bias_data$simulation_type), collapse = ", "), "\n")
  cat("- Shifts:", paste(unique(bias_data$shift_label), collapse = ", "), "\n")
  
  if (nrow(bias_data) == 0) {
    cat("No data for bias plot!\n")
    return(NULL)
  }
  
  ggplot(bias_data, aes(x = simulation_type, y = bias, fill = simulation_type)) +
    geom_violin(alpha = 0.7, trim = FALSE) +
    geom_boxplot(width = 0.1, alpha = 0.3, outlier.shape = NA) +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black") +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 1) +
    facet_wrap(~shift_label, scales = "free_x") +
    scale_fill_manual(values = sim_colors) +
    labs(
      y = "Bias\n(Estimated - True)",
      fill = "Simulation Type"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      legend.position = "bottom",
      strip.text = element_text(size = 10)
    )
}

# Function to create trajectory plot
create_trajectory_plot <- function(results, max_iterations = 50) {
  trajectory_data <- results$selection_coefficients %>%
    mutate(
      simulation_type = case_when(
        dominance_type == "additive_n15" ~ "Additive (n=15)",
        dominance_type == "additive_n40" ~ "Additive (n=40)",
        dominance_type == "recessive_n40" ~ "Recessive (n=40)",
        TRUE ~ "Unknown"
      )
    ) %>%
    filter(simulation_type %in% names(sim_colors), !is.na(delta_p)) %>%
    group_by(iteration, shift, simulation_type) %>%
    summarise(
      mean_p_before = mean(p_before, na.rm = TRUE),
      mean_p_after = mean(p_after, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      time_before = case_when(
        shift == "shift1" ~ 143,
        shift == "shift2" ~ 93,
        shift == "shift3" ~ 43
      ),
      time_after = case_when(
        shift == "shift1" ~ 139,
        shift == "shift2" ~ 89,
        shift == "shift3" ~ 39
      )
    )
  
  # Sample iterations if too many
  if (length(unique(trajectory_data$iteration)) > max_iterations) {
    sampled_iterations <- sample(unique(trajectory_data$iteration), max_iterations)
    trajectory_data <- trajectory_data %>% filter(iteration %in% sampled_iterations)
  }
  
  # Convert to long format for plotting
  trajectory_data_long <- trajectory_data %>%
    tidyr::pivot_longer(
      cols = c(mean_p_before, mean_p_after),
      names_to = "timepoint",
      values_to = "mean_frequency"
    ) %>%
    mutate(
      years_before_present = case_when(
        timepoint == "mean_p_before" ~ time_before,
        timepoint == "mean_p_after" ~ time_after
      ),
      point_type = case_when(
        timepoint == "mean_p_before" ~ "Before",
        timepoint == "mean_p_after" ~ "After"
      )
    )
  
  # Create data for inter-breakpoint connections (grey lines)
  inter_data <- trajectory_data %>%
    arrange(iteration, simulation_type, shift) %>%
    group_by(iteration, simulation_type) %>%
    filter(n() == 3) %>%
    do({
      conn1 <- data.frame(
        iteration = .$iteration[1],
        simulation_type = .$simulation_type[1],
        x = c(.$time_after[.$shift == "shift1"], .$time_before[.$shift == "shift2"]),
        y = c(.$mean_p_after[.$shift == "shift1"], .$mean_p_before[.$shift == "shift2"]),
        connection = "1to2"
      )
      conn2 <- data.frame(
        iteration = .$iteration[1],
        simulation_type = .$simulation_type[1],
        x = c(.$time_after[.$shift == "shift2"], .$time_before[.$shift == "shift3"]),
        y = c(.$mean_p_after[.$shift == "shift2"], .$mean_p_before[.$shift == "shift3"]),
        connection = "2to3"
      )
      rbind(conn1, conn2)
    }) %>%
    ungroup()
  
  # Add final timepoint data for grey line connection only
  final_data <- results$raw_frequencies %>%
    filter(timepoint == "last") %>%
    mutate(
      simulation_type = case_when(
        dominance_type == "additive_n15" ~ "Additive (n=15)",
        dominance_type == "additive_n40" ~ "Additive (n=40)",
        dominance_type == "recessive_n40" ~ "Recessive (n=40)",
        TRUE ~ "Unknown"
      )
    ) %>%
    filter(simulation_type %in% names(sim_colors)) %>%
    group_by(iteration, simulation_type) %>%
    summarise(
      mean_p_final = mean(allele_freq, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Add connection from shift3_after to final timepoint (at time 0) - same simple approach
  final_connections <- trajectory_data %>%
    filter(shift == "shift3") %>%
    inner_join(final_data, by = c("iteration", "simulation_type")) %>%
    group_by(iteration, simulation_type) %>%
    do({
      data.frame(
        iteration = .$iteration[1],
        simulation_type = .$simulation_type[1],
        x = c(.$time_after[1], 0),  # From shift3_after time to present (0)
        y = c(.$mean_p_after[1], .$mean_p_final[1]),  # From shift3_after freq to final freq
        connection = "3tofinal"
      )
    }) %>%
    ungroup()
  
  # Combine all connection data
  all_inter_data <- rbind(inter_data, final_connections)
  
  # Create the plot
  ggplot(trajectory_data_long, aes(x = years_before_present, y = mean_frequency)) +
    geom_vline(xintercept = c(141, 91, 41), color = "black", linetype = "solid", size = 1.2) +
    geom_line(data = all_inter_data,
              aes(x = x, y = y, group = paste(iteration, simulation_type, connection)),
              color = "grey70", alpha = 0.4, size = 0.3) +
    geom_line(aes(group = paste(iteration, shift, simulation_type), color = simulation_type),
              alpha = 0.6, size = 0.4) +
    geom_point(aes(color = simulation_type, shape = point_type),
               alpha = 0.7, size = 1.5) +
    scale_x_reverse(limits = c(145, 0)) +
    scale_color_manual(name = "Simulation Type", values = sim_colors) +
    scale_shape_manual(name = "Timepoint", values = c("Before" = 16, "After" = 17)) +
    coord_cartesian(ylim = c(0, 0.35)) +
    labs(
      x = "Years Before Present",
      y = "Mean Allele Frequency"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
}

# USAGE - Updated to be more flexible
cat("=== UPDATED STREAMLINED VCF ANALYSIS SCRIPT ===\n")
cat("Updated for new directory structure.\n")
cat("Will auto-detect directories with 'instants' in the name.\n")
cat("File patterns:\n")
cat("- *_sample_n15_h5.vcf (Additive n=15, teal)\n")
cat("- *_sample_n40_h5.vcf (Additive n=40, purple)\n")
cat("- *_sample_n40_h015.vcf (Recessive n=40, yellow)\n\n")

# You can specify directories manually or let it auto-detect
# Example manual specification:
# base_directories <- c("breakpoint_instants_additiven15_outputs", 
#                       "breakpoint_instants_additiven40_outputs", 
#                       "breakpoint_instants_recessiven40_outputs")
# results <- run_analysis(base_dirs = base_directories, max_iterations = 100)

cat("Running analysis with auto-detection...\n")
results <- run_analysis(max_iterations = 100)

cat("Creating bias plot...\n")
bias_plot <- create_bias_plot(results)
if (!is.null(bias_plot)) {
  print(bias_plot)
}

cat("Creating trajectory plot...\n")
trajectory_plot <- create_trajectory_plot(results)
print(trajectory_plot)

cat("\nAnalysis complete!\n")
cat("Updated color scheme:\n")
cat("- Additive (n=15): Teal\n")
cat("- Additive (n=40): Deep purple\n")
cat("- Recessive (n=40): Yellow\n")

# Calculate sample counts per timepoint across iterations
if (nrow(results$raw_frequencies) > 0) {
  sample_counts_summary <- results$raw_frequencies %>%
    filter(timepoint != "last") %>%  # Exclude the final timepoint for now
    # Only keep additive types (exclude recessive)
    filter(dominance_type %in% c("additive_n15", "additive_n40")) %>%
    group_by(simulation_type = case_when(
      dominance_type == "additive_n15" ~ "Additive (n=15)",
      dominance_type == "additive_n40" ~ "Additive (n=40)",
      TRUE ~ "Unknown"
    ), timepoint) %>%
    summarise(
      mean_sample_size = mean(sample_size, na.rm = TRUE),
      n_iterations = length(unique(iteration)),
      .groups = 'drop'
    ) %>%
    # Map timepoints to years before present (same as trajectory plot)
    mutate(
      years_before_present = case_when(
        timepoint == "shift1_before" ~ 143,
        timepoint == "shift1_after" ~ 139,
        timepoint == "shift2_before" ~ 93,
        timepoint == "shift2_after" ~ 89,
        timepoint == "shift3_before" ~ 43,
        timepoint == "shift3_after" ~ 39
      )
    )
  
  # Create the sample distribution plot
  if (require(ggpattern, quietly = TRUE)) {
    p_sample_counts <- ggplot(sample_counts_summary, 
                              aes(x = years_before_present, y = mean_sample_size, 
                                  fill = simulation_type, pattern = simulation_type,
                                  pattern_fill = simulation_type)) +
      geom_col_pattern(position = "stack", alpha = 0.8, width = 1.5,
                       pattern_density = 0.6, pattern_spacing = 0.1,
                       colour = NA, pattern_colour = NA) +
      scale_fill_manual(name = "Sample Type", 
                        values = c("Additive (n=15)" = "#21908CFF", "Additive (n=40)" = "#440154FF")) +
      scale_pattern_fill_manual(name = "Sample Type",
                                values = c("Additive (n=15)" = "#21908CFF", "Additive (n=40)" = "#FDC100FF")) +
      scale_pattern_manual(name = "Sample Type", 
                           values = c("Additive (n=15)" = "none", "Additive (n=40)" = "stripe"), 
                           guide = "none") +
      scale_x_reverse(limits = c(145, 0)) +
      labs(
        title = "Sample Distribution Across Time",
        x = "Years Before Present",
        y = "Mean Sample Size"
      ) +
      theme_minimal() +
      theme(
        legend.position = "bottom"
      )
    
    print(p_sample_counts)
  } else {
    cat("ggpattern package not available, skipping sample distribution plot\n")
  }
}
