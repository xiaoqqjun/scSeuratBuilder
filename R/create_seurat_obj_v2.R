#' Create Seurat Object with Enhanced Validation and Per-Sample Processing
#'
#' Enhanced version of create_seurat_obj with:
#' - Per-sample processing loop
#' - 10X triplet file validation
#' - Barcode format consistency check
#' - Feature column detection (2 vs 3 columns)
#' - Intermediate file saving for resumability
#' - Per-sample gene synchronization
#' - No scale.data generation (saves storage)
#'
#' @param data_dir Character. Path to data directory
#' @param metadata_file Character. Path to metadata file (xlsx/xls). If NULL, will look for metadata.xlsx/xls in data_dir
#' @param species Character. Species for gene synchronization. Options: "homo"/"human"/"hsa" (human),
#'   "mus"/"mouse"/"mmu" (mouse), "rat"/"rattus"/"rno" (rat). Required when enable_gene_sync = TRUE.
#' @param enable_gene_sync Logical. Enable gene name synchronization (default: FALSE)
#' @param min_cells Integer. Minimum cells per feature (default: 3)
#' @param min_features Integer. Minimum features per cell (default: 100)
#' @param verbose Logical. Print detailed progress (default: TRUE)
#' @param output_dir Character. Path to output directory for final object. If NULL, object is not saved (default: NULL)
#' @param save_intermediate Logical. Save per-sample RDS files (default: FALSE)
#' @param intermediate_dir Character. Directory for intermediate files. If NULL, uses data_dir/intermediate (default: NULL)
#' @param resume Logical. Skip already processed samples (default: FALSE)
#' @param skip_samples Character vector. Samples to skip (default: NULL)
#' @param on_sample_error Character. How to handle sample errors: "stop", "skip", or "warn" (default: "stop")
#' @param check_10x Logical. Validate 10X triplet files before processing (default: TRUE)
#' @param skip_scale_data Logical. Skip scale.data generation to save storage (default: TRUE)
#'
#' @return A Seurat object with integrated metadata
#'
#' @details
#' ## Workflow
#' 1. Validate metadata file (mandatory)
#' 2. Detect data structure and format
#' 3. Validate 10X files (if check_10x = TRUE)
#' 4. Match samples with metadata
#' 5. Process each sample:
#'    - Read data
#'    - Create Seurat object
#'    - Gene synchronization (per-sample)
#'    - Save intermediate (if save_intermediate = TRUE)
#' 6. Merge all objects
#' 7. Integrate metadata
#' 8. Return final object
#'
#' ## Per-Sample Processing Benefits
#' - Progress tracking for each sample
#' - Resume from last processed sample
#' - Skip problematic samples
#' - Save intermediate results
#' - Better error handling
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' obj <- create_seurat_obj_v2(
#'   data_dir = "path/to/Data",
#'   metadata_file = "path/to/metadata.xlsx"
#' )
#'
#' # With intermediate saving (for resumability)
#' obj <- create_seurat_obj_v2(
#'   data_dir = "path/to/Data",
#'   save_intermediate = TRUE,
#'   intermediate_dir = "path/to/intermediate"
#' )
#'
#' # Resume from previous run
#' obj <- create_seurat_obj_v2(
#'   data_dir = "path/to/Data",
#'   save_intermediate = TRUE,
#'   resume = TRUE
#' )
#'
#' # With gene synchronization
#' obj <- create_seurat_obj_v2(
#'   data_dir = "path/to/Data",
#'   species = "human",
#'   enable_gene_sync = TRUE
#' )
#' }
#'
#' @export
#'
create_seurat_obj_v2 <- function(data_dir,
                                 metadata_file = NULL,
                                 species = NULL,
                                 enable_gene_sync = FALSE,
                                 min_cells = 3,
                                 min_features = 100,
                                 verbose = TRUE,
                                 output_dir = NULL,
                                 save_intermediate = FALSE,
                                 intermediate_dir = NULL,
                                 resume = FALSE,
                                 skip_samples = NULL,
                                 on_sample_error = c("stop", "skip", "warn"),
                                 check_10x = TRUE,
                                 skip_scale_data = TRUE) {

  # ============================================================================
  # Step 0: Initial Validation
  # ============================================================================

  start_time <- Sys.time()

  if (missing(data_dir)) {
    stop("data_dir is required!\n",
         "Example: create_seurat_obj_v2(data_dir = \"path/to/Data\")")
  }

  on_sample_error <- match.arg(on_sample_error)

  # Remove trailing slash
  data_dir <- gsub("[/\\\\]$", "", data_dir)

  if (!dir.exists(data_dir)) {
    stop("Data directory does not exist: ", data_dir)
  }

  # Set up intermediate directory
  if (is.null(intermediate_dir)) {
    intermediate_dir <- file.path(data_dir, ".scSeuratBuilder_intermediate")
  }

  if (verbose) {
    cat("\n")
    cat(paste(rep("=", 70), collapse = ""), "\n")
    cat("scSeuratBuilder v2: Enhanced Seurat Object Builder\n")
    cat("Mode: Per-sample processing with validation\n")
    cat(paste(rep("=", 70), collapse = ""), "\n")
    cat("\n")
    cat("Data directory:", data_dir, "\n")
    cat("Intermediate directory:", intermediate_dir, "\n")
    cat("Save intermediate:", save_intermediate, "\n")
    cat("Resume mode:", resume, "\n")
    cat("Check 10X files:", check_10x, "\n")
    cat("\n")
  }

  # ============================================================================
  # Step 1: Validate Metadata File (MANDATORY)
  # ============================================================================

  cat("[Step 1/9] Validating metadata file...\n")

  # Look for metadata file if not specified
  if (is.null(metadata_file)) {
    all_files <- list.files(data_dir, ignore.case = TRUE)
    metadata_idx <- grep("^metadata\\\\.xls", all_files)

    if (length(metadata_idx) > 0) {
      metadata_file <- file.path(data_dir, all_files[metadata_idx[1]])
    } else {
      xlsx_idx <- grep("\\.xlsx?$", all_files)
      if (length(xlsx_idx) == 1) {
        metadata_file <- file.path(data_dir, all_files[xlsx_idx[1]])
      }
    }
  }

  if (is.null(metadata_file) || !file.exists(metadata_file)) {
    cat("\n")
    cat("========================================\n")
    cat("ERROR: Metadata file not found!\n")
    cat("========================================\n")
    cat("\n")
    stop("Metadata file is required")
  }

  metadata <- readxl::read_excel(metadata_file)

  if (!("samples" %in% tolower(names(metadata)))) {
    stop("Invalid metadata structure: missing 'samples' column")
  }

  names(metadata) <- tolower(names(metadata))

  if (!("groups" %in% names(metadata))) {
    stop("Invalid metadata structure: missing 'groups' column")
  }

  metadata_samples <- as.character(metadata$samples)
  n_metadata_samples <- length(metadata_samples)

  cat("  ✓ Metadata file found:", basename(metadata_file), "\n")
  cat("  ✓ Samples in metadata:", n_metadata_samples, "\n")
  cat("\n")

  # ============================================================================
  # Step 2: Detect Data Structure and Format
  # ============================================================================

  cat("[Step 2/9] Detecting data structure...\n")

  structure_info <- detect_data_structure_enhanced(data_dir, verbose = verbose)

  data_type <- structure_info$type
  data_format <- structure_info$format

  cat("  Detected structure:", data_type, "\n")
  cat("  Detected format:", data_format, "\n")

  if (data_format == "h5") {
    cat("\n")
    cat("========================================\n")
    cat("H5 Format Detected\n")
    cat("========================================\n")
    cat("\n")
    cat("H5 files are detected but not fully supported in v2.\n")
    cat("Please use create_seurat_obj() for H5 files, or\n")
    cat("convert H5 to 10X triplet format first.\n")
    cat("\n")
    stop("H5 format not supported in create_seurat_obj_v2()")
  }

  cat("\n")

  # ============================================================================
  # Step 3: Validate 10X Triplet Files (if applicable)
  # ============================================================================

  if (check_10x && data_type %in% c("subfolders", "direct_files")) {
    cat("[Step 3/9] Checking 10X triplet files...\n")

    validation_result <- validate_all_10x_files(
      data_dir = data_dir,
      data_type = data_type,
      metadata_samples = metadata_samples,
      verbose = verbose
    )

    # Show warnings if any (now non-blocking)
    if (length(validation_result$warnings) > 0) {
      cat("\n")
      cat("  ⚠ Warnings:\n")
      for (w in validation_result$warnings) {
        cat("    -", w, "\n")
      }
      cat("  Continuing with processing...\n")
    }

    cat("  ✓ Feature columns:", validation_result$n_feature_cols, "\n")
    cat("  ✓ Barcode format:", validation_result$barcode_format, "\n")
    cat("\n")
  } else {
    cat("[Step 3/9] Skipping 10X validation...\n\n")
  }

  # ============================================================================
  # Step 4: Match Samples with Metadata
  # ============================================================================

  cat("[Step 4/9] Matching samples with metadata...\n")

  match_result <- match_samples_with_metadata(
    data_dir = data_dir,
    data_type = data_type,
    metadata_samples = metadata_samples,
    verbose = verbose
  )

  if (!match_result$success) {
    cat("\n")
    cat("========================================\n")
    cat("ERROR: Sample name mismatch!\n")
    cat("========================================\n")
    cat("\n")
    stop("Sample name mismatch detected")
  }

  samples_to_process <- match_result$matched_samples

  # Apply skip_samples filter
  if (!is.null(skip_samples)) {
    samples_skipped <- samples_to_process %in% skip_samples
    if (any(samples_skipped)) {
      cat("  Skipping samples:", paste(samples_to_process[samples_skipped], collapse = ", "), "\n")
      samples_to_process <- samples_to_process[!samples_skipped]
    }
  }

  n_samples <- length(samples_to_process)

  cat("  ✓ Samples to process:", n_samples, "\n")
  cat("\n")

  # ============================================================================
  # Step 5: Create Intermediate Directory
  # ============================================================================

  if (save_intermediate) {
    cat("[Step 5/9] Setting up intermediate directory...\n")

    if (!dir.exists(intermediate_dir)) {
      dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)
    }

    cat("  ✓ Intermediate directory:", intermediate_dir, "\n")
    cat("\n")
  } else {
    cat("[Step 5/9] Skipping intermediate directory...\n\n")
  }

  # ============================================================================
  # Step 6: Process Each Sample (Per-Sample Loop)
  # ============================================================================

  cat("[Step 6/9] Processing samples (per-sample loop)...\n\n")

  seurat_list <- list()
  processed_samples <- character()
  failed_samples <- character()
  skipped_samples <- character()

  # Track gene sync info
  gene_sync_info <- list()

  for (i in seq_along(samples_to_process)) {
    sample_name <- samples_to_process[i]

    cat(sprintf("  [%d/%d] %s", i, n_samples, sample_name), "... ")

    # Check if already processed (resume mode)
    intermediate_file <- file.path(intermediate_dir, paste0(sample_name, ".rds"))

    if (resume && save_intermediate && file.exists(intermediate_file)) {
      cat("SKIPPED (already processed)")
      if (verbose) cat(" - loading from", intermediate_file)
      cat("\n")

      seurat_list[[sample_name]] <- readRDS(intermediate_file)
      skipped_samples <- c(skipped_samples, sample_name)
      next
    }

    # Process sample
    sample_result <- tryCatch({
      process_single_sample(
        sample_name = sample_name,
        data_dir = data_dir,
        data_type = data_type,
        file_mapping = match_result$file_mapping,
        min_cells = min_cells,
        min_features = min_features,
        species = species,
        enable_gene_sync = enable_gene_sync,
        skip_scale_data = skip_scale_data,
        verbose = verbose
      )
    }, error = function(e) {
      list(
        success = FALSE,
        error = conditionMessage(e),
        obj = NULL
      )
    })

    # Handle result
    if (!sample_result$success) {
      if (on_sample_error == "stop") {
        cat("ERROR!\n")
        cat("    Error:", sample_result$error, "\n")
        stop("Sample processing failed. Use on_sample_error = 'skip' to continue.")
      } else if (on_sample_error == "warn") {
        cat("ERROR (continuing)...\n")
        cat("    Error:", sample_result$error, "\n")
        failed_samples <- c(failed_samples, sample_name)
      } else {
        # skip
        cat("SKIPPED (error)\n")
        failed_samples <- c(failed_samples, sample_name)
      }
      next
    }

    obj <- sample_result$obj

    # Save intermediate
    if (save_intermediate) {
      saveRDS(obj, intermediate_file)
      if (verbose) cat(" + saved")
    }

    # Store gene sync info
    if (!is.null(sample_result$gene_sync_summary)) {
      gene_sync_info[[sample_name]] <- sample_result$gene_sync_summary
    }

    seurat_list[[sample_name]] <- obj
    processed_samples <- c(processed_samples, sample_name)

    # Show summary
    n_genes <- nrow(obj)
    n_cells <- ncol(obj)
    cat(sprintf("OK (%d genes x %d cells)", n_genes, n_cells))

    if (verbose && !is.null(sample_result$gene_sync_summary)) {
      cat(sprintf(" [synced: %d -> %d genes]",
                  sample_result$gene_sync_summary$before,
                  sample_result$gene_sync_summary$after))
    }

    cat("\n")
  }

  cat("\n")
  cat("  Processed:", length(processed_samples), "\n")
  if (length(skipped_samples) > 0) {
    cat("  Skipped (resume):", length(skipped_samples), "\n")
  }
  if (length(failed_samples) > 0) {
    cat("  Failed:", length(failed_samples), "\n")
  }
  cat("\n")

  if (length(processed_samples) == 0) {
    stop("No samples were successfully processed.")
  }

  # ============================================================================
  # Step 7: Merge Objects
  # ============================================================================

  cat("[Step 7/9] Merging Seurat objects...\n")

  # Filter to only successfully processed samples
  seurat_list <- seurat_list[processed_samples]

  merged_obj <- merge_seurat_objects_safe(seurat_list, verbose = verbose)

  cat("  ✓ Objects merged\n")
  cat("\n")

  # ============================================================================
  # Step 8: Integrate Metadata
  # ============================================================================

  cat("[Step 8/9] Integrating metadata...\n")

  merged_obj <- integrate_metadata_safe(merged_obj, metadata, verbose = verbose)

  cat("  ✓ Metadata integrated\n")
  cat("\n")

  # ============================================================================
  # Step 9: Final Summary and Return
  # ============================================================================

  cat("[Step 9/9] Finalizing...\n")

  end_time <- Sys.time()
  build_time <- difftime(end_time, start_time, units = "secs")

  # Save if output_dir is specified
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    output_file <- file.path(output_dir, "seurat_object_v2.rds")
    saveRDS(merged_obj, file = output_file, compress = FALSE)
    cat("  ✓ Object saved:", output_file, "\n")
  } else {
    cat("  ✓ Object created (not saved)\n")
  }

  cat("\n")

  if (verbose) {
    cat("\n")
    cat(paste(rep("=", 70), collapse = ""), "\n")
    cat("✓ Pipeline completed successfully!\n")
    cat(paste(rep("=", 70), collapse = ""), "\n")
    cat("\n")
    cat("Seurat object: merged_obj\n")
    cat("  Dimensions:", nrow(merged_obj), "genes x", ncol(merged_obj), "cells\n")
    cat("  Samples:", length(processed_samples), "\n")
    cat("  Build time:", round(as.numeric(build_time), 2), "seconds\n")
    cat("\n")

    # Show groups
    if ("groups" %in% names(merged_obj@meta.data)) {
      cat("Groups:\n")
      groups_table <- table(merged_obj@meta.data$groups)
      for (g in names(groups_table)) {
        cat(sprintf("  - %s: %d cells\n", g, groups_table[g]))
      }
    }

    # Show failed samples
    if (length(failed_samples) > 0) {
      cat("\nFailed samples:\n")
      for (s in failed_samples) {
        cat("  -", s, "\n")
      }
    }
    cat("\n")
  }

  return(merged_obj)
}


#' Process Single Sample
#'
#' Internal function to process a single sample through all steps
#'
#' @param sample_name Sample name
#' @param data_dir Data directory
#' @param data_type Data structure type
#' @param file_mapping File mapping (for direct_files type)
#' @param min_cells Minimum cells per feature
#' @param min_features Minimum features per cell
#' @param species Species for gene sync
#' @param enable_gene_sync Enable gene synchronization
#' @param skip_scale_data Skip scale.data generation
#' @param verbose Print progress
#'
#' @return List with success, obj, and gene_sync_summary
#'
#' @keywords internal
#'
process_single_sample <- function(sample_name,
                                   data_dir,
                                   data_type,
                                   file_mapping = NULL,
                                   min_cells = 3,
                                   min_features = 100,
                                   species = NULL,
                                   enable_gene_sync = FALSE,
                                   skip_scale_data = TRUE,
                                   verbose = FALSE) {

  # Wrap entire processing in suppressWarnings to handle non-unique features
  obj <- suppressWarnings({

    # Determine sample path
    if (data_type == "subfolders") {
      sample_path <- file.path(data_dir, sample_name)
      format_info <- detect_folder_format(sample_path)

      obj <- switch(format_info$format,
        "threebrothers" = read_threebrothers(sample_path, sample_name,
                                               min_cells, min_features),
        "rds" = read_seurat_file(file.path(sample_path,
                                            list.files(sample_path)[grepl("\\.rds$", list.files(sample_path), ignore.case = TRUE)][1])),
        "rdata" = read_seurat_file(file.path(sample_path,
                                              list.files(sample_path)[grepl("\\.(rda|rdata)$", list.files(sample_path), ignore.case = TRUE)][1]))
      )

    } else if (data_type == "direct_files") {
      file_name <- file_mapping[[sample_name]]
      file_path <- file.path(data_dir, file_name)

      ext <- tolower(tools::file_ext(file_name))

      obj <- switch(ext,
        "rds" = read_seurat_file(file_path),
        "rda" = read_seurat_file(file_path),
        "rdata" = read_seurat_file(file_path),
        stop("Unsupported file format: ", ext)
      )
    }
  })

  if (is.null(obj) || !inherits(obj, "Seurat")) {
    return(list(
      success = FALSE,
      error = "Failed to create Seurat object",
      obj = NULL
    ))
  }

  # Gene synchronization (per-sample)
  gene_sync_summary <- NULL
  if (enable_gene_sync && !is.null(species)) {
    n_before <- nrow(obj)

    obj <- tryCatch({
      geneSync::gene_sync_add_to_obj(
        obj,
        species = species,
        update_features = TRUE,
        verbose = FALSE
      )
    }, error = function(e) {
      warning("Gene sync failed for ", sample_name, ": ", conditionMessage(e))
      obj
    })

    n_after <- nrow(obj)
    gene_sync_summary <- list(
      before = n_before,
      after = n_after
    )
  }

  # Skip scale.data to save storage (V4/V5 compatible)
  if (skip_scale_data) {
    tryCatch({
      default_assay <- Seurat::DefaultAssay(obj)
      assay_obj <- obj[[default_assay]]

      # Detect Seurat version: V5 has @layers, V4 has @slots
      has_layers <- tryCatch({
        exists("layers", inherits = "Assay", mode = function(env) {
          is(env, "environment") || is(env, "list")
        }) && !is.null(assay_obj@layers)
      }, error = function(e) FALSE)

      if (has_layers && length(assay_obj@layers) > 0) {
        # Seurat V5: Remove scale layer
        if ("scale" %in% names(assay_obj@layers)) {
          # Use slot assignment for V5
          if (is(assay_obj@layers[["scale"]], "dgCMatrix") ||
              is(assay_obj@layers[["scale"]], "ANY") ||
              is.null(assay_obj@layers[["scale"]])) {
            assay_obj@layers[["scale"]] <- NULL
          }
          obj[[default_assay]] <- assay_obj
        }
      } else {
        # Seurat V4: Handle scale.data slot
        # First check if scale.data slot exists and can be modified
        if (.hasSlot(assay_obj, "scale.data")) {
          scale_data <- assay_obj@scale.data
          # Only try to remove if it's a matrix (not S4 object)
          if (is.matrix(scale_data) || is.null(scale_data)) {
            assay_obj@scale.data <- NULL
            obj[[default_assay]] <- assay_obj
          }
        }
      }

      # Remove ScaleData command to prevent regeneration
      if ("ScaleData" %in% names(obj@commands)) {
        obj@commands$ScaleData <- NULL
      }

    }, error = function(e) {
      # Silently ignore - if we can't remove scale.data, it's okay
    })
  }

  return(list(
    success = TRUE,
    obj = obj,
    gene_sync_summary = gene_sync_summary
  ))
}


#' Detect Data Structure and Format (Enhanced)
#'
#' Enhanced detection that also identifies file format
#'
#' @param data_dir Data directory
#' @param verbose Print progress
#'
#' @return List with type and format
#'
#' @keywords internal
#'
detect_data_structure_enhanced <- function(data_dir, verbose = TRUE) {

  all_items <- list.files(data_dir)
  is_dir <- dir.exists(file.path(data_dir, all_items))
  subfolders <- all_items[is_dir]
  files <- all_items[!is_dir]

  n_subfolders <- length(subfolders)
  n_files <- length(files)

  # Check for H5 files
  has_h5 <- any(grepl("\\.h5$", files, ignore.case = TRUE))

  # Check for RDS/RDATA files
  has_rds <- any(grepl("\\.(rds|rda|rdata)$", files, ignore.case = TRUE))

  # Determine type
  if (n_subfolders > 0) {
    type <- "subfolders"
  } else if (n_files > 0) {
    type <- "direct_files"
  } else {
    type <- "unknown"
  }

  # Determine format
  if (has_h5) {
    format <- "h5"
  } else if (has_rds && n_files == 1) {
    format <- "seurat_object"
  } else {
    format <- "10x_triplet"
  }

  return(list(
    type = type,
    format = format,
    n_subfolders = n_subfolders,
    n_files = n_files
  ))
}


#' Validate All 10X Files
#'
#' Check barcode format consistency and feature columns across all samples
#'
#' @param data_dir Data directory
#' @param data_type Data structure type
#' @param metadata_samples Sample names from metadata
#' @param verbose Print progress
#'
#' @return List with validation results
#'
#' @keywords internal
#'
validate_all_10x_files <- function(data_dir,
                                    data_type,
                                    metadata_samples,
                                    verbose = TRUE) {

  result <- list(
    valid = TRUE,
    warnings = character(),
    n_feature_cols = NULL,
    barcode_format = NULL
  )

  # Collect barcode formats and feature column counts
  barcode_formats <- character()
  feature_col_counts <- integer()

  for (sample in metadata_samples) {
    if (data_type == "subfolders") {
      sample_path <- file.path(data_dir, sample)
    } else {
      # For direct_files, find the corresponding file
      files <- list.files(data_dir)
      sample_file <- files[grepl(paste0("^", sample, "[._]"), files, ignore.case = TRUE)][1]
      if (is.na(sample_file)) {
        sample_file <- files[grep(sample, files, ignore.case = TRUE)][1]
      }
      sample_path <- file.path(data_dir, sample_file)
    }

    # Check barcode format
    barcode_file <- find_file_in_path(sample_path, "barcode")
    if (!is.null(barcode_file)) {
      barcodes <- tryCatch({
        read.table(barcode_file, header = FALSE, stringsAsFactors = FALSE)[, 1]
      }, error = function(e) NULL)

      if (!is.null(barcodes) && length(barcodes) > 0) {
        # Check barcode format: standard (no _) or multipart (has _)
        has_underscore <- any(grepl("_", barcodes))
        barcode_formats <- c(barcode_formats, ifelse(has_underscore, "multipart", "standard"))
      }
    }

    # Check feature columns
    feature_file <- find_file_in_path(sample_path, "feature")
    if (!is.null(feature_file)) {
      features <- tryCatch({
        read.table(feature_file, header = FALSE, stringsAsFactors = FALSE, quote = "")
      }, error = function(e) NULL)

      if (!is.null(features)) {
        feature_col_counts <- c(feature_col_counts, ncol(features))
      }
    }
  }

  # Check barcode format consistency - now warning only
  if (length(barcode_formats) > 0) {
    unique_formats <- unique(barcode_formats)
    result$barcode_format <- paste(unique_formats, collapse = "/")

    if (length(unique_formats) > 1) {
      result$warnings <- c(result$warnings,
                          "Inconsistent barcode formats across samples")
    }
  }

  # Check feature column consistency - now warning only
  if (length(feature_col_counts) > 0) {
    unique_cols <- unique(feature_col_counts)
    result$n_feature_cols <- unique_cols

    if (length(unique_cols) > 1) {
      result$warnings <- c(result$warnings,
                          sprintf("Feature columns vary: %s",
                                  paste(unique_cols, collapse = ", ")))
    }
  }

  # Always return valid (warnings don't stop processing)
  result$valid <- TRUE

  return(result)
}


#' Find File in Path (Internal Helper)
#'
#' @param path Directory to search
#' @param pattern File name pattern
#'
#' @return File path or NULL
#'
#' @keywords internal
#'
find_file_in_path <- function(path, pattern) {
  if (!dir.exists(path)) return(NULL)

  files <- list.files(path, ignore.case = TRUE)

  # Direct match
  matches <- files[grepl(pattern, files, ignore.case = TRUE)]

  if (length(matches) == 0) return(NULL)

  # Check for standard 10X subdirectories
  for (subdir in c("filtered_feature_bc_matrix", "raw_feature_bc_matrix")) {
    subdir_path <- file.path(path, subdir)
    if (dir.exists(subdir_path)) {
      subdir_files <- list.files(subdir_path, ignore.case = TRUE)
      subdir_matches <- subdir_files[grepl(pattern, subdir_files, ignore.case = TRUE)]
      if (length(subdir_matches) > 0) {
        return(file.path(subdir_path, subdir_matches[1]))
      }
    }
  }

  # Return first match
  return(file.path(path, matches[1]))
}


#' Match Samples with Metadata (Internal)
#'
#' @param data_dir Data directory
#' @param data_type Data structure type
#' @param metadata_samples Sample names from metadata
#' @param verbose Print progress
#'
#' @return List with match results
#'
#' @keywords internal
#'
match_samples_with_metadata <- function(data_dir,
                                         data_type,
                                         metadata_samples,
                                         verbose = TRUE) {

  if (data_type == "subfolders") {
    return(match_samples_subfolders(data_dir, metadata_samples, verbose))
  } else if (data_type == "direct_files") {
    return(match_samples_direct_files(data_dir, metadata_samples, verbose))
  } else {
    # single_file - always succeed
    return(list(
      success = TRUE,
      matched_samples = metadata_samples,
      file_mapping = NULL
    ))
  }
}
