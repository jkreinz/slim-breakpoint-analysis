#!/usr/bin/env Rscript
library(data.table)
library(dplyr)
library(segmented)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 2) {
  cat("Usage: Rscript summarise_results.R <iter_num> <file_prefix>\n")
  cat("Example: Rscript summarise_results.R 1 m2_allsamples_nomulti_iter\n")
  quit(status = 1)
}

iter_num <- args[1]
file_prefix <- args[2]

# Extract simulation type from file prefix for output naming
if(grepl("n70", file_prefix)) {
  sim_type <- "additive_n70"
} else if(grepl("recessive", file_prefix)) {
  sim_type <- "recessive_n200"
} else {
  sim_type <- "additive_n200"
}

# Read data files
genos_file <- paste0(file_prefix, iter_num, ".012")
inds_file <- paste0(file_prefix, iter_num, ".012.indv")
snps_file <- paste0(file_prefix, iter_num, ".012.pos")
selcoef_file <- paste0("selection_coefficients_iter", iter_num, ".txt")

# Check if required files exist
if(!file.exists(genos_file) || !file.exists(inds_file)) {
  cat("Error: Required data files not found for iteration", iter_num, "\n")
  quit(status = 1)
}

# Read data
genos <- read.table(genos_file)
inds <- read.table(inds_file)
if(file.exists(snps_file)) {
  snps <- read.table(snps_file)
}
if(file.exists(selcoef_file)) {
  selcoef <- read.table(selcoef_file)
}

# Extract year information - handles new VCF naming patterns
# Updated for patterns like: t128_sample_h50_n200_i0, t130_sample_h50_n200_i0, etc.
year <- gsub("i","",gsub("t","",(gsub("_sample_i0","",inds$V1))))
year2 <- as.numeric(unlist(lapply(strsplit(year,"_"), `[[`, 1)))

# Convert genotype data to long format
genos_long <- data.table::melt(genos[,-1], id.vars = NULL)
genos_long$year <- rep(year2, ncol(genos)-1)
genos_long$locus <- rep(1:(ncol(genos)-1), each = nrow(genos))
genos_long$value <- genos_long$value/2
genos_long <- genos_long[!is.na(genos_long$value),]

# Fit base linear model
lm_base <- glm(value ~ year, family="binomial", data=genos_long)

# Fit all models (0-4 breakpoints)
results <- data.frame(n_breakpoints = 0:4, aic = NA, bic = NA, converged = FALSE)
results$aic[1] <- AIC(lm_base)
results$bic[1] <- BIC(lm_base)
results$converged[1] <- TRUE

# Store all models
all_models <- list(lm_base)

# Fit segmented models with 1-4 breakpoints
for(i in 1:4) {
  tryCatch({
    seg_model <- segmented(lm_base, seg.Z = ~year, npsi = i)
    results$aic[i+1] <- AIC(seg_model)
    results$bic[i+1] <- BIC(seg_model)
    results$converged[i+1] <- TRUE
    all_models[[i+1]] <- seg_model
  }, error = function(e) {
    results$converged[i+1] <<- FALSE
    all_models[[i+1]] <<- NULL
  })
}

# Select best models
best_aic <- which.min(results$aic[results$converged])
best_bic <- which.min(results$bic[results$converged])
best_aic_bp <- results$n_breakpoints[best_aic]
best_bic_bp <- results$n_breakpoints[best_bic]

# USE BIC MODEL for extracting breakpoints and slopes
best_model <- all_models[[best_bic]]

# Extract breakpoints and slopes from BIC-selected model
if(best_bic_bp > 0) {
  breakpoints_result <- as.data.frame(best_model$psi)
  breakpoints <- breakpoints_result$Est.
  breakpoint_ses <- breakpoints_result$St.Err
  slopes_result <- as.data.frame(slope(best_model)$year)
  slopes <- slopes_result$Est. * 2
  slope_ses <- slopes_result$St.Err. * 2
} else {
  breakpoints <- NA
  breakpoint_ses <- NA
  slopes <- coef(best_model)[2] * 2
  slope_ses <- summary(best_model)$coefficients[2,2] * 2
}

# Create standardized output with up to 4 breakpoints and 5 slopes
# Initialize all as NA
bp_1 <- bp_2 <- bp_3 <- bp_4 <- NA
bp_1_se <- bp_2_se <- bp_3_se <- bp_4_se <- NA
s_1 <- s_2 <- s_3 <- s_4 <- s_5 <- NA
s_1_se <- s_2_se <- s_3_se <- s_4_se <- s_5_se <- NA

# Fill in breakpoint values based on what we have
if(length(breakpoints) >= 1 && !is.na(breakpoints[1])) {
  bp_1 <- round(breakpoints[1], 2)
  bp_1_se <- round(breakpoint_ses[1], 3)
}
if(length(breakpoints) >= 2 && !is.na(breakpoints[2])) {
  bp_2 <- round(breakpoints[2], 2)
  bp_2_se <- round(breakpoint_ses[2], 3)
}
if(length(breakpoints) >= 3 && !is.na(breakpoints[3])) {
  bp_3 <- round(breakpoints[3], 2)
  bp_3_se <- round(breakpoint_ses[3], 3)
}
if(length(breakpoints) >= 4 && !is.na(breakpoints[4])) {
  bp_4 <- round(breakpoints[4], 2)
  bp_4_se <- round(breakpoint_ses[4], 3)
}

# Fill in slope values
if(length(slopes) >= 1 && !is.na(slopes[1])) {
  s_1 <- round(slopes[1], 4)
  s_1_se <- round(slope_ses[1], 4)
}
if(length(slopes) >= 2 && !is.na(slopes[2])) {
  s_2 <- round(slopes[2], 4)
  s_2_se <- round(slope_ses[2], 4)
}
if(length(slopes) >= 3 && !is.na(slopes[3])) {
  s_3 <- round(slopes[3], 4)
  s_3_se <- round(slope_ses[3], 4)
}
if(length(slopes) >= 4 && !is.na(slopes[4])) {
  s_4 <- round(slopes[4], 4)
  s_4_se <- round(slope_ses[4], 4)
}
if(length(slopes) >= 5 && !is.na(slopes[5])) {
  s_5 <- round(slopes[5], 4)
  s_5_se <- round(slope_ses[5], 4)
}

# Create summary line with standardized columns
summary_line <- paste(iter_num, best_aic_bp, best_bic_bp,
                     bp_1, bp_1_se, bp_2, bp_2_se, bp_3, bp_3_se, bp_4, bp_4_se,
                     s_1, s_1_se, s_2, s_2_se, s_3, s_3_se, s_4, s_4_se, s_5, s_5_se,
                     sep = ",")

# Write summary file with simulation type in filename
output_file <- paste0("summary_", sim_type, "_iter", iter_num, ".csv")
write.table(summary_line, output_file, row.names = FALSE, col.names = FALSE, quote = FALSE)
