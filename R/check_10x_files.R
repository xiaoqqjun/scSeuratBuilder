#' Check 10X Data File Alignment
#'
#' Comprehensive validation of 10X Genomics triplet files (barcodes.tsv.gz,
#' features.tsv.gz, matrix.mtx.gz) to ensure structural consistency before
#' Seurat object creation.
#'
#' @param data_dir Character. Path to directory containing 10X files
#' @param feature_col Integer. Which column in features.tsv.gz to use as gene
#'   names (default: 2, typically gene symbol). Use 1 for Ensembl ID.
#' @param make_unique_names Logical. Whether to make duplicate names unique
#'   using \code{make.unique()} (default: TRUE)
#' @param verbose Logical. Print detailed progress (default: TRUE)
#'
#' @return A list with components:
#' \describe{
#'   \item{valid}{Logical. TRUE if all checks passed}
#'   \item{matrix}{The sparse matrix with row/col names assigned}
#'   \item{features}{Data frame of feature metadata}
#'   \item{barcodes}{Data frame of barcode metadata}
#'   \item{summary}{List with validation statistics}
#'   \item{errors}{Character vector of error messages}
#'   \item{warnings}{Character vector of warning messages}
#' }
#'
#' @details
#' This function performs the following checks:
#' \enumerate{
#'   \item File existence verification
#'   \item Dimension alignment (features == matrix rows, barcodes == matrix cols)
#'   \item MTX header consistency
#'   \item Matrix index bounds checking (no out-of-range indices)
#'   \item Name quality checks (NA, empty, duplicated)
#' }
#'
#' @examples
#' \dontrun{
#' # Check a single sample directory
#' result <- check_10x_triplets("path/to/sample/")
#'
#' if (result$valid) {
#'   cat("All checks passed!\n")
#' } else {
#'   cat("Errors found:\n")
#'   print(result$errors)
#' }
#'
#' # Access the validated matrix
#' mat <- result$matrix
#' }
#'
#' @export
check_10x_triplets <- function(data_dir,
                               feature_col = 2,
                               make_unique_names = TRUE,
                               verbose = TRUE) {

  # ============================================================================
  # File Path Setup
  # ============================================================================

  barcodes_file <- file.path(data_dir, "barcodes.tsv.gz")
  features_file <- file.path(data_dir, "features.tsv.gz")
  matrix_file <- file.path(data_dir, "matrix.mtx.gz")

  # Initialize result
  result <- list(
    valid = TRUE,
    matrix = NULL,
    features = NULL,
    barcodes = NULL,
    summary = list(),
    errors = character(),
    warnings = character()
  )

  # ============================================================================
  # 1. File Existence Check
  # ============================================================================

  if (verbose) cat("=== 10X Triplet File Check ===\n")
  if (verbose) cat("Directory:", data_dir, "\n\n")

  missing_files <- character()
  if (!file.exists(barcodes_file)) missing_files <- c(missing_files, "barcodes.tsv.gz")
  if (!file.exists(features_file)) missing_files <- c(missing_files, "features.tsv.gz")
  if (!file.exists(matrix_file)) missing_files <- c(missing_files, "matrix.mtx.gz")

  if (length(missing_files) > 0) {
    result$valid <- FALSE
    result$errors <- c(result$errors,
                      paste("Missing files:", paste(missing_files, collapse = ", ")))
    stop("Missing required 10X files: ", paste(missing_files, collapse = ", "))
  }

  if (verbose) cat("[1/6] File existence check... OK\n")

  # ============================================================================
  # 2. Read Files
  # ============================================================================

  if (verbose) cat("[2/6] Reading files...\n")

  barcodes <- tryCatch({
    read.table(
      gzfile(barcodes_file),
      sep = "\t",
      header = FALSE,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    result$errors <<- c(result$errors, paste("Failed to read barcodes:", e$message))
    NULL
  })

  features <- tryCatch({
    read.table(
      gzfile(features_file),
      sep = "\t",
      header = FALSE,
      stringsAsFactors = FALSE,
      quote = ""
    )
  }, error = function(e) {
    result$errors <<- c(result$errors, paste("Failed to read features:", e$message))
    NULL
  })

  mat <- tryCatch({
    Matrix::readMM(gzfile(matrix_file))
  }, error = function(e) {
    result$errors <<- c(result$errors, paste("Failed to read matrix:", e$message))
    NULL
  })

  if (is.null(barcodes) || is.null(features) || is.null(mat)) {
    result$valid <- FALSE
    stop("File reading failed. Check error messages.")
  }

  if (verbose) {
    cat("  Matrix dim:", nrow(mat), "x", ncol(mat), "\n")
    cat("  Features rows:", nrow(features), "\n")
    cat("  Barcodes rows:", nrow(barcodes), "\n")
  }

  # ============================================================================
  # 3. Dimension Alignment Check
  # ============================================================================

  if (verbose) cat("[3/6] Dimension alignment check...\n")

  dim_errors <- character()

  if (nrow(features) != nrow(mat)) {
    msg <- sprintf("Features rows (%d) != Matrix rows (%d)",
                   nrow(features), nrow(mat))
    dim_errors <- c(dim_errors, msg)
    result$errors <- c(result$errors, msg)
  }

  if (nrow(barcodes) != ncol(mat)) {
    msg <- sprintf("Barcodes rows (%d) != Matrix cols (%d)",
                   nrow(barcodes), ncol(mat))
    dim_errors <- c(dim_errors, msg)
    result$errors <- c(result$errors, msg)
  }

  if (length(dim_errors) > 0) {
    result$valid <- FALSE
    if (verbose) {
      cat("  ERROR: Dimension mismatch!\n")
      for (e in dim_errors) cat("   -", e, "\n")
    }
  } else {
    if (verbose) cat("  OK: Dimensions aligned\n")
  }

  # ============================================================================
  # 4. Name Quality Check
  # ============================================================================

  if (verbose) cat("[4/6] Name quality check...\n")

  if (feature_col > ncol(features)) {
    result$valid <- FALSE
    result$errors <- c(result$errors,
                      paste("feature_col =", feature_col,
                            "exceeds features columns =", ncol(features)))
    stop("feature_col exceeds features table columns")
  }

  feature_names <- as.character(features[, feature_col])
  barcode_names <- as.character(barcodes[, 1])

  # Count issues
  na_features <- sum(is.na(feature_names))
  empty_features <- sum(feature_names == "" | is.na(feature_names))
  dup_features <- sum(duplicated(feature_names))

  na_barcodes <- sum(is.na(barcode_names))
  empty_barcodes <- sum(barcode_names == "" | is.na(barcode_names))
  dup_barcodes <- sum(duplicated(barcode_names))

  if (verbose) {
    cat("  Feature names:\n")
    cat("    NA:", na_features, "\n")
    cat("    Empty:", empty_features, "\n")
    cat("    Duplicated:", dup_features, "\n")
    cat("  Barcode names:\n")
    cat("    NA:", na_barcodes, "\n")
    cat("    Empty:", empty_barcodes, "\n")
    cat("    Duplicated:", dup_barcodes, "\n")
  }

  # Store in summary
  result$summary$name_quality <- list(
    na_features = na_features,
    empty_features = empty_features,
    dup_features = dup_features,
    na_barcodes = na_barcodes,
    empty_barcodes = empty_barcodes,
    dup_barcodes = dup_barcodes
  )

  # Warnings for quality issues
  if (na_features > 0 || empty_features > 0) {
    msg <- sprintf("Feature names have %d NA/empty values", na_features + empty_features)
    result$warnings <- c(result$warnings, msg)
  }
  if (dup_features > 0) {
    msg <- sprintf("Feature names have %d duplicates", dup_features)
    result$warnings <- c(result$warnings, msg)
  }
  if (na_barcodes > 0 || empty_barcodes > 0) {
    msg <- sprintf("Barcode names have %d NA/empty values", na_barcodes + empty_barcodes)
    result$warnings <- c(result$warnings, msg)
  }

  # ============================================================================
  # 5. MTX Index Bounds Check
  # ============================================================================

  if (verbose) cat("[5/6] MTX index bounds check...\n")

  mtx_check <- .check_mtx_indices(matrix_file, nrow(features), nrow(barcodes))

  if (!mtx_check$valid) {
    result$valid <- FALSE
    result$errors <- c(result$errors, mtx_check$errors)
    if (verbose) cat("  ERROR: Matrix index out of bounds!\n")
  } else {
    if (verbose) cat("  OK: Matrix indices within bounds\n")
  }

  # Store MTX info
  result$summary$mtx <- list(
    header_rows = mtx_check$header_rows,
    header_cols = mtx_check$header_cols,
    header_nnz = mtx_check$header_nnz,
    actual_nnz = mtx_check$actual_nnz
  )

  # ============================================================================
  # 6. Assign Names to Matrix
  # ============================================================================

  if (verbose) cat("[6/6] Assigning names to matrix...\n")

  if (make_unique_names) {
    feature_names <- make.unique(feature_names)
    barcode_names <- make.unique(barcode_names)
    if (verbose) cat("  Names made unique\n")
  }

  rownames(mat) <- feature_names
  colnames(mat) <- barcode_names

  # ============================================================================
  # Finalize Result
  # ============================================================================

  result$matrix <- mat
  result$features <- features
  result$barcodes <- barcodes

  if (verbose) {
    cat("\n=== Summary ===\n")
    cat("Valid:", result$valid, "\n")
    cat("Errors:", length(result$errors), "\n")
    cat("Warnings:", length(result$warnings), "\n")
    if (result$valid) {
      cat("\nMatrix ready for Seurat object creation.\n")
    } else {
      cat("\nPlease fix errors before creating Seurat object.\n")
    }
  }

  # Add class
  class(result) <- c("check_10x_result", "list")

  return(result)
}


#' Check MTX Matrix Indices
#'
#' Internal function to verify that indices in the MTX file are within
#' valid bounds.
#'
#' @param matrix_file Path to matrix.mtx.gz file
#' @param n_features Expected number of features (rows)
#' @param n_barcodes Expected number of barcodes (columns)
#'
#' @return List with validation results
#'
#' @keywords internal
.check_mtx_indices <- function(matrix_file, n_features, n_barcodes) {

  result <- list(
    valid = TRUE,
    errors = character(),
    header_rows = NA,
    header_cols = NA,
    header_nnz = NA,
    actual_nnz = NA
  )

  con <- gzfile(matrix_file, open = "rt")
  lines <- readLines(con)
  close(con)

  # Remove comment lines
  lines2 <- lines[!grepl("^%", lines)]

  if (length(lines2) < 1) {
    result$valid <- FALSE
    result$errors <- "Empty MTX file"
    return(result)
  }

  # Parse header
  header <- strsplit(lines2[1], "\\s+")[[1]]
  header <- header[header != ""]

  if (length(header) < 3) {
    result$valid <- FALSE
    result$errors <- "Invalid MTX header format"
    return(result)
  }

  result$header_rows <- as.integer(header[1])
  result$header_cols <- as.integer(header[2])
  result$header_nnz <- as.integer(header[3])

  # Check header consistency
  if (result$header_rows != n_features) {
    result$errors <- c(result$errors,
                       sprintf("MTX header rows (%d) != features count (%d)",
                               result$header_rows, n_features))
  }
  if (result$header_cols != n_barcodes) {
    result$errors <- c(result$errors,
                       sprintf("MTX header cols (%d) != barcodes count (%d)",
                               result$header_cols, n_barcodes))
  }

  # Parse non-zero entries
  if (length(lines2) > 1) {
    dat <- data.table::fread(
      text = paste(lines2[-1], collapse = "\n"),
      header = FALSE,
      showProgress = FALSE
    )
    data.table::setnames(dat, c("i", "j", "x"))

    result$actual_nnz <- nrow(dat)

    # Check bounds
    out_of_bounds_rows <- any(dat$i < 1 | dat$i > n_features)
    out_of_bounds_cols <- any(dat$j < 1 | dat$j > n_barcodes)

    if (out_of_bounds_rows) {
      result$valid <- FALSE
      result$errors <- c(result$errors,
                         "Matrix contains row indices out of feature bounds")
    }
    if (out_of_bounds_cols) {
      result$valid <- FALSE
      result$errors <- c(result$errors,
                         "Matrix contains column indices out of barcode bounds")
    }
  }

  return(result)
}


#' Check Seurat Object Consistency with 10X Files
#'
#' Verify that a Seurat object matches the original 10X triplet files.
#' Useful for validating that object creation was successful.
#'
#' @param seurat_obj A Seurat object
#' @param data_dir Character. Path to original 10X data directory
#' @param assay Character. Assay name to check (default: "RNA")
#' @param slot Character. Slot to check (default: "counts")
#' @param verbose Logical. Print detailed comparison (default: TRUE)
#'
#' @return Logical. TRUE if object matches original files
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(scSeuratBuilder)
#'
#' # Create object
#' obj <- create_seurat_obj("data/")
#'
#' # Verify consistency
#' check_seurat_consistency(obj, "data/")
#' }
#'
#' @export
check_seurat_consistency <- function(seurat_obj,
                                     data_dir,
                                     assay = "RNA",
                                     slot = "counts",
                                     verbose = TRUE) {

  if (verbose) cat("=== Seurat Object Consistency Check ===\n\n")

  # ============================================================================
  # Get Object Matrix
  # ============================================================================

  obj_mat <- Seurat::GetAssayData(seurat_obj, assay = assay, slot = slot)

  if (verbose) {
    cat("Object matrix dim:", dim(obj_mat), "\n")
    cat("Object assay:", assay, "\n")
    cat("Object slot:", slot, "\n\n")
  }

  # ============================================================================
  # Read Original Files
  # ============================================================================

  check_result <- check_10x_triplets(data_dir, verbose = FALSE, make_unique_names = FALSE)

  if (!check_result$valid) {
    warning("Original 10X files have errors. Cannot compare.")
    return(FALSE)
  }

  orig_mat <- check_result$matrix

  # ============================================================================
  # Compare Dimensions
  # ============================================================================

  if (verbose) cat("[1/3] Dimension comparison...\n")

  dim_match <- all(dim(obj_mat) == dim(orig_mat))

  if (verbose) {
    cat("  Original dim:", dim(orig_mat), "\n")
    cat("  Object dim:   ", dim(obj_mat), "\n")
    cat("  Match:", dim_match, "\n\n")
  }

  # ============================================================================
  # Compare Row Names
  # ============================================================================

  if (verbose) cat("[2/3] Row name comparison...\n")

  # Handle make.unique differences
  obj_rownames <- rownames(obj_mat)
  orig_rownames <- rownames(orig_mat)

  # If original has duplicates, object should have make.unique applied
  if (any(duplicated(orig_rownames))) {
    expected_rownames <- make.unique(orig_rownames)
  } else {
    expected_rownames <- orig_rownames
  }

  rownames_match <- identical(obj_rownames, expected_rownames)

  if (verbose) {
    cat("  Original rows:", length(orig_rownames), "\n")
    cat("  Object rows:   ", length(obj_rownames), "\n")
    cat("  Match:", rownames_match, "\n")
    if (!rownames_match && verbose) {
      cat("  First mismatch at position:",
          which(obj_rownames != expected_rownames)[1], "\n")
    }
    cat("\n")
  }

  # ============================================================================
  # Compare Column Names
  # ============================================================================

  if (verbose) cat("[3/3] Column name comparison...\n")

  obj_colnames <- colnames(obj_mat)
  orig_colnames <- colnames(orig_mat)

  if (any(duplicated(orig_colnames))) {
    expected_colnames <- make.unique(orig_colnames)
  } else {
    expected_colnames <- orig_colnames
  }

  colnames_match <- identical(obj_colnames, expected_colnames)

  if (verbose) {
    cat("  Original cols:", length(orig_colnames), "\n")
    cat("  Object cols:   ", length(obj_colnames), "\n")
    cat("  Match:", colnames_match, "\n")
    if (!colnames_match && verbose) {
      cat("  First mismatch at position:",
          which(obj_colnames != expected_colnames)[1], "\n")
    }
    cat("\n")
  }

  # ============================================================================
  # Final Result
  # ============================================================================

  all_match <- dim_match && rownames_match && colnames_match

  if (verbose) {
    cat("=== Result ===\n")
    if (all_match) {
      cat("OK: Seurat object matches original 10X files.\n")
    } else {
      cat("WARNING: Seurat object does NOT match original files.\n")
    }
  }

  return(all_match)
}


#' Batch Check 10X Directories
#'
#' Check multiple 10X data directories in batch.
#'
#' @param data_dirs Character vector of directory paths
#' @param feature_col Integer. Which column in features.tsv.gz to use (default: 2)
#' @param verbose Logical. Print detailed progress (default: TRUE)
#'
#' @return Data frame with check results for each directory
#'
#' @examples
#' \dontrun{
#' dirs <- c("sample1/", "sample2/", "sample3/")
#' results <- batch_check_10x(dirs)
#' print(results)
#'
#' # Get only failed directories
#' failed <- results$dir[!results$valid]
#' }
#'
#' @export
batch_check_10x <- function(data_dirs,
                            feature_col = 2,
                            verbose = TRUE) {

  if (verbose) {
    cat("=== Batch 10X Check ===\n")
    cat("Directories:", length(data_dirs), "\n\n")
  }

  results_list <- lapply(data_dirs, function(dir) {
    if (verbose) cat("Checking:", dir, "\n")

    result <- tryCatch({
      check <- check_10x_triplets(dir, feature_col = feature_col, verbose = FALSE)
      list(
        dir = dir,
        valid = check$valid,
        n_features = nrow(check$features),
        n_barcodes = nrow(check$barcodes),
        nnz = check$summary$mtx$actual_nnz,
        dup_features = check$summary$name_quality$dup_features,
        dup_barcodes = check$summary$name_quality$dup_barcodes,
        n_errors = length(check$errors),
        n_warnings = length(check$warnings)
      )
    }, error = function(e) {
      list(
        dir = dir,
        valid = FALSE,
        n_features = NA,
        n_barcodes = NA,
        nnz = NA,
        dup_features = NA,
        dup_barcodes = NA,
        n_errors = 1,
        n_warnings = 0,
        error_msg = conditionMessage(e)
      )
    })

    if (verbose) {
      if (result$valid) {
        cat("  OK\n")
      } else {
        cat("  FAILED (", result$n_errors, " errors)\n")
      }
    }

    result
  })

  # Convert to data frame
  results_df <- do.call(rbind, lapply(results_list, as.data.frame))
  rownames(results_df) <- NULL

  if (verbose) {
    cat("\n=== Summary ===\n")
    cat("Total:", nrow(results_df), "\n")
    cat("Passed:", sum(results_df$valid), "\n")
    cat("Failed:", sum(!results_df$valid), "\n")
  }

  return(results_df)
}


#' Print Check Result Summary
#'
#' @param x A check_10x_result object
#' @param ... Additional arguments (unused)
#'
#' @export
print.check_10x_result <- function(x, ...) {
  cat("=== 10X Triplet Check Result ===\n\n")

  cat("Status:", ifelse(x$valid, "VALID", "INVALID"), "\n\n")

  cat("Dimensions:\n")
  cat("  Features:", x$summary$mtx$header_rows, "\n")
  cat("  Barcodes:", x$summary$mtx$header_cols, "\n")
  cat("  Non-zero entries:", x$summary$mtx$actual_nnz, "\n\n")

  cat("Name Quality:\n")
  cat("  Duplicated features:", x$summary$name_quality$dup_features, "\n")
  cat("  Duplicated barcodes:", x$summary$name_quality$dup_barcodes, "\n\n")

  if (length(x$errors) > 0) {
    cat("Errors:\n")
    for (e in x$errors) cat("  -", e, "\n")
    cat("\n")
  }

  if (length(x$warnings) > 0) {
    cat("Warnings:\n")
    for (w in x$warnings) cat("  -", w, "\n")
    cat("\n")
  }

  invisible(x)
}
