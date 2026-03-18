#!/usr/bin/env Rscript
# scSeuratBuilder 10X File Check Examples
# This script demonstrates how to use the 10X file validation functions

library(scSeuratBuilder)

# ==============================================================================
# Example 1: Check a single 10X directory
# ==============================================================================

cat("========================================\n")
cat("Example 1: Single Directory Check\n")
cat("========================================\n\n")

# Replace with your actual data directory
data_dir <- "path/to/your/10x/data/"

# Run the check
result <- check_10x_triplets(
  data_dir = data_dir,
  feature_col = 2,        # Use gene symbol (column 2)
  make_unique_names = TRUE,
  verbose = TRUE
)

# Check if valid
if (result$valid) {
  cat("\n==> Data is valid! Ready for Seurat object creation.\n")

  # Access the validated matrix
  mat <- result$matrix
  cat("Matrix dimensions:", dim(mat), "\n")
} else {
  cat("\n==> Errors found:\n")
  print(result$errors)
}

# ==============================================================================
# Example 2: Check Seurat object consistency
# ==============================================================================

cat("\n========================================\n")
cat("Example 2: Seurat Object Consistency\n")
cat("========================================\n\n")

# Assuming you have a Seurat object
# obj <- create_seurat_obj(data_dir)

# Verify the object matches original files
# is_consistent <- check_seurat_consistency(obj, data_dir, verbose = TRUE)

# ==============================================================================
# Example 3: Batch check multiple directories
# ==============================================================================

cat("\n========================================\n")
cat("Example 3: Batch Check\n")
cat("========================================\n\n")

# Example: Check all sample directories
# sample_dirs <- c(
#   "data/sample1/",
#   "data/sample2/",
#   "data/sample3/"
# )

# batch_result <- batch_check_10x(sample_dirs, verbose = TRUE)

# View results
# print(batch_result)

# Get only failed directories
# failed_dirs <- batch_result$dir[!batch_result$valid]
# if (length(failed_dirs) > 0) {
#   cat("\nFailed directories:\n")
#   print(failed_dirs)
# }

# ==============================================================================
# Example 4: Check with Ensembl IDs instead of symbols
# ==============================================================================

cat("\n========================================\n")
cat("Example 4: Use Ensembl IDs\n")
cat("========================================\n\n")

# result_ensembl <- check_10x_triplets(
#   data_dir = data_dir,
#   feature_col = 1,        # Use Ensembl ID (column 1)
#   verbose = TRUE
# )

# ==============================================================================
# Example 5: Access detailed summary
# ==============================================================================

cat("\n========================================\n")
cat("Example 5: Detailed Summary\n")
cat("========================================\n\n")

# Print custom summary
# print_summary <- function(check_result) {
#   cat("Check Summary for:", check_result$dir, "\n")
#   cat("  Valid:", check_result$valid, "\n")
#   cat("  Features:", check_result$summary$mtx$header_rows, "\n")
#   cat("  Barcodes:", check_result$summary$mtx$header_cols, "\n")
#   cat("  Non-zero entries:", check_result$summary$mtx$actual_nnz, "\n")
#   cat("  Duplicated features:", check_result$summary$name_quality$dup_features, "\n")
#   cat("  Duplicated barcodes:", check_result$summary$name_quality$dup_barcodes, "\n")
#   cat("  Errors:", length(check_result$errors), "\n")
#   cat("  Warnings:", length(check_result$warnings), "\n")
# }

# print_summary(result)

cat("\n========================================\n")
cat("Examples completed!\n")
cat("========================================\n")
