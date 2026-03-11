#' Create Seurat Object with Strict Metadata Validation
#'
#' Main function to build Seurat objects with mandatory metadata validation.
#' Supports three data structures:
#' - Type A: Subfolder per sample
#' - Type B: Direct files in data_dir
#' - Type C: Single merged object
#'
#' @param data_dir Character. Path to data directory
#' @param metadata_file Character. Path to metadata file (xlsx/xls). If NULL, will look for metadata.xlsx/xls in data_dir
#' @param species Character. Species for gene synchronization. Options: "homo"/"human"/"hsa" (human),
#'   "mus"/"mouse"/"mmu" (mouse), "rat"/"rattus"/"rno" (rat). Required when enable_gene_sync = TRUE.
#' @param enable_gene_sync Logical. Enable gene name synchronization (default: FALSE)
#' @param min_cells Integer. Minimum cells per feature (default: 3)
#' @param min_features Integer. Minimum features per cell (default: 100)
#' @param verbose Logical. Print detailed progress (default: TRUE)
#'
#' @return A Seurat object with integrated metadata
#'
#' @details
#' ## Workflow
#' 1. Validate metadata file (mandatory)
#' 2. Detect data structure (Type A/B/C)
#' 3. Match samples with metadata
#' 4. Read and create Seurat objects
#' 5. (Optional) Gene synchronization
#' 6. Merge and integrate metadata
#' 7. Return final object
#'
#' ## Metadata Requirements
#' - Column 1: samples (must match sample names)
#' - Column 2: groups (sample grouping)
#' - Additional columns: optional
#'
#' @export
#'
create_seurat_obj <- function(data_dir,
                               metadata_file = NULL,
                               species = NULL,
                               enable_gene_sync = FALSE,
                               min_cells = 3,
                               min_features = 100,
                               verbose = TRUE) {

  # ============================================================================
  # Step 0: Initial Validation
  # ============================================================================

  if (missing(data_dir)) {
    stop("data_dir is required!\n",
         "Example: create_seurat_obj(data_dir = \"path/to/Data\")")
  }

  # Remove trailing slash
  data_dir <- gsub("[/\\\\]$", "", data_dir)

  if (!dir.exists(data_dir)) {
    stop("Data directory does not exist: ", data_dir)
  }

  if (verbose) {
    cat("\n")
    cat(paste(rep("=", 70), collapse = ""), "\n")
    cat("scSeuratBuilder: Intelligent Seurat Object Builder\n")
    cat("Version: 1.0.0\n")
    cat(paste(rep("=", 70), collapse = ""), "\n")
    cat("\n")
    cat("Data directory:", data_dir, "\n")
    cat("\n")
  }

  # ============================================================================
  # Step 1: Validate Metadata File (MANDATORY)
  # ============================================================================

  cat("[Step 1/7] Validating metadata file...\n")

  # Look for metadata file if not specified
  if (is.null(metadata_file)) {
    # List all files and look for metadata.xlsx/xls
    all_files <- list.files(data_dir, ignore.case = TRUE)

    # Find metadata file
    metadata_idx <- grep("^metadata\\\\.xls", all_files)

    if (length(metadata_idx) > 0) {
      metadata_file <- file.path(data_dir, all_files[metadata_idx[1]])
    } else {
      # Try to find any Excel file
      xlsx_idx <- grep("\\.xlsx?$", all_files)
      if (length(xlsx_idx) == 1) {
        # Only one Excel file, assume it's metadata
        metadata_file <- file.path(data_dir, all_files[xlsx_idx[1]])
      }
    }
  }

  # Debug info
  cat("  Looking for metadata file in:", data_dir, "\n")
  cat("  Files found:", paste(list.files(data_dir), collapse = ", "), "\n")
  if (!is.null(metadata_file)) {
    cat("  Found:", metadata_file, "\n")
  }

  # Validate metadata file exists
  if (is.null(metadata_file) || !file.exists(metadata_file)) {
    cat("\n")
    cat("========================================\n")
    cat("ERROR: Metadata file not found!\n")
    cat("========================================\n")
    cat("\n")
    cat("Please provide a metadata file in Excel format (.xlsx or .xls)\n")
    cat("The metadata file MUST contain:\n")
    cat("  - Column 1: samples (sample names)\n")
    cat("  - Column 2: groups (sample grouping)\n")
    cat("\n")
    cat("Place the file in the data directory with one of these names:\n")
    cat("  - metadata.xlsx\n")
    cat("  - metadata.xls\n")
    cat("\n")
    cat("Or specify the file path:\n")
    cat("  create_seurat_obj(data_dir, metadata_file = \"path/to/metadata.xlsx\")\n")
    cat("\n")
    stop("Metadata file is required")
  }

  # Read and validate metadata structure
  metadata <- readxl::read_excel(metadata_file)

  if (!("samples" %in% tolower(names(metadata)))) {
    cat("\n")
    cat("========================================\n")
    cat("ERROR: Invalid metadata structure!\n")
    cat("========================================\n")
    cat("\n")
    cat("The metadata file must have:\n")
    cat("  - First column named 'samples'\n")
    cat("  - Second column named 'groups'\n")
    cat("\n")
    cat("Current columns:", paste(names(metadata), collapse = ", "), "\n")
    cat("\n")
    stop("Invalid metadata structure: missing 'samples' column")
  }

  # Normalize column names
  names(metadata) <- tolower(names(metadata))

  if (!("groups" %in% names(metadata))) {
    cat("\n")
    cat("========================================\n")
    cat("ERROR: Invalid metadata structure!\n")
    cat("========================================\n")
    cat("\n")
    cat("The metadata file must have a 'groups' column (second column)\n")
    cat("\n")
    cat("Current columns:", paste(names(metadata), collapse = ", "), "\n")
    cat("\n")
    stop("Invalid metadata structure: missing 'groups' column")
  }

  # Extract sample names from metadata
  metadata_samples <- as.character(metadata$samples)
  n_metadata_samples <- length(metadata_samples)

  cat("  ✓ Metadata file found:", basename(metadata_file), "\n")
  cat("  ✓ Samples in metadata:", n_metadata_samples, "\n")
  cat("\n")

  # ============================================================================
  # Step 2: Detect Data Structure
  # ============================================================================

  cat("[Step 2/7] Detecting data structure...\n")

  structure_info <- detect_data_structure(data_dir, verbose = verbose)

  data_type <- structure_info$type

  cat("  Detected structure:", data_type, "\n")
  cat("\n")

  # ============================================================================
  # Step 3: Extract and Match Samples
  # ============================================================================

  cat("[Step 3/7] Extracting and matching samples...\n")

  match_result <- switch(data_type,
    "subfolders" = match_samples_subfolders(data_dir, metadata_samples, verbose),
    "direct_files" = match_samples_direct_files(data_dir, metadata_samples, verbose),
    "single_file" = match_samples_single_file(data_dir, metadata_samples, verbose)
  )

  # Check if matching succeeded
  # - Type A/B: strict matching required, fail if mismatch
  # - Type C: merged object, show warning but continue
  if (!match_result$success) {
    cat("\n")
    cat("========================================\n")
    cat("ERROR: Sample name mismatch!\n")
    cat("========================================\n")
    cat("\n")

    if (length(match_result$missing_in_data) > 0) {
      cat("Samples in metadata but NOT in data:\n")
      for (s in match_result$missing_in_data) {
        cat("  -", s, "\n")
      }
      cat("\n")
    }

    if (length(match_result$missing_in_metadata) > 0) {
      cat("Samples in data but NOT in metadata:\n")
      for (s in match_result$missing_in_metadata) {
        cat("  -", s, "\n")
      }
      cat("\n")
    }

    cat("This mismatch will affect downstream analysis and grouping.\n")
    cat("Please check:\n")
    cat("  1. Sample names in metadata (first column)\n")
    cat("  2. Folder names / file names in data directory\n")
    cat("\n")
    cat("Sample names must match exactly!\n")
    cat("\n")
    stop("Sample name mismatch detected")
  }

  # For Type C (merged object), show warning if metadata doesn't fully match
  if (data_type == "single_file" && !match_result$metadata_matched) {
    cat("\n")
    cat("⚠ Warning: Sample name mismatch detected!\n")
    cat("\n")

    if (length(match_result$missing_in_data) > 0) {
      cat("Samples in metadata but NOT in object:\n")
      for (s in match_result$missing_in_data) {
        cat("  -", s, "\n")
      }
      cat("\n")
    }

    if (length(match_result$missing_in_metadata) > 0) {
      cat("Samples in object but NOT in metadata:\n")
      for (s in match_result$missing_in_metadata) {
        cat("  -", s, "\n")
      }
      cat("\n")
    }

    cat("Proceeding with available samples...\n")
    cat("  Matched samples:", length(match_result$matched_samples), "\n")
    cat("\n")
  } else if (data_type != "single_file") {
    # Type A/B: show success message
    cat("  ✓ All samples matched successfully!\n")
    cat("  Matched samples:", length(match_result$matched_samples), "\n")
    cat("\n")
  } else {
    # Type C with full match
    cat("  ✓ All samples matched successfully!\n")
    cat("  Matched samples:", length(match_result$matched_samples), "\n")
    cat("\n")
  }

  # ============================================================================
  # Step 4: Read Data and Create Seurat Objects
  # ============================================================================

  cat("[Step 4/7] Reading data and creating Seurat objects...\n")

  seurat_list <- switch(data_type,
    "subfolders" = read_data_from_subfolders(
      data_dir, match_result$matched_samples,
      min_cells, min_features, verbose
    ),
    "direct_files" = read_data_from_direct_files(
      data_dir, match_result$matched_samples, match_result$file_mapping,
      min_cells, min_features, verbose
    ),
    "single_file" = read_merged_object(
      data_dir, verbose
    )
  )

  if (length(seurat_list) == 0) {
    stop("No Seurat objects were created")
  }

  cat("  ✓ Created", length(seurat_list), "Seurat object(s)\n")
  cat("\n")

  # ============================================================================
  # Step 5: Gene Synchronization (Optional)
  # ============================================================================

  if (enable_gene_sync && !is.null(species)) {
    cat("[Step 5/7] Gene name synchronization...\n")

    # IMPORTANT: Always sync when enable_gene_sync = TRUE!
    # Purpose: Eliminate Cell Ranger version differences (e.g., Gene-old vs Gene-new)
    # Even if genes are symbols, sync will standardize them using authority symbols

    # Show sample gene names for verification
    first_obj <- seurat_list[[1]]
    sample_genes <- head(rownames(first_obj), 10)
    cat("  Sample gene names:", paste(sample_genes, collapse = ", "), "\n")

    cat("  Purpose: Unify gene names across samples processed with different Cell Ranger versions\n")
    cat("  Syncing all samples -> standardized gene names ...\n")

    seurat_list <- sync_genes_for_list(
      seurat_list, species, verbose
    )

    cat("  ✓ Gene synchronization completed\n")
    cat("\n")
  } else {
    cat("[Step 5/7] Gene synchronization skipped\n")
    cat("\n")
  }

  # ============================================================================
  # Step 6: Merge Objects and Integrate Metadata
  # ============================================================================

  cat("[Step 6/7] Merging objects and integrating metadata...\n")

  if (data_type == "single_file") {
    # For single file, object may already be merged
    merged_obj <- seurat_list[[1]]

    # Special warning for Type C
    cat("\n")
    cat("╔══════════════════════════════════════════════════════════════╗\n")
    cat("║  IMPORTANT: You loaded a pre-merged object!                      ║\n")
    cat("╠══════════════════════════════════════════════════════════════╣\n")
    cat("║  Please verify:                                                   ║\n")
    cat("║    1. Sample names (orig.ident) match your expectations          ║\n")
    cat("║    2. Group information is correct                               ║\n")
    cat("║    3. Data comes from a reliable source                          ║\n")
    cat("╚══════════════════════════════════════════════════════════════╝\n")
    cat("\n")

    # Check if metadata was successfully matched
    if (match_result$metadata_matched) {
      cat("  ✓ Metadata matched successfully\n")
    } else {
      cat("  ⚠ Metadata matching had issues - please review!\n")
    }

  } else {
    # Merge multiple objects
    merged_obj <- merge_seurat_objects_safe(seurat_list, verbose)

    cat("  ✓ Objects merged\n")

    # Integrate metadata
    merged_obj <- integrate_metadata_safe(merged_obj, metadata, verbose)

    cat("  ✓ Metadata integrated\n")
  }

  cat("\n")

  # ============================================================================
  # Step 7: Save and Return
  # ============================================================================

  cat("[Step 7/7] Saving and returning...\n")

  # Create output directory
  rundate <- gsub("-", "_", Sys.Date())
  output_dir <- paste0("scSeuratBuilder_output_", rundate)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }

  # Save object with progress indication
  output_file <- file.path(output_dir, "seurat_object.rds")

  cat("  Saving Seurat object...\n")
  cat("    Saving...\r")

  saveRDS(merged_obj, file = output_file, compress = FALSE)

  cat("    Done!\n")

  file_size <- file.size(output_file)
  size_mb <- round(file_size / 1024 / 1024, 2)
  cat(sprintf("    File size: %.2f MB\n", size_mb))

  # Save summary
  summary_file <- file.path(output_dir, "summary.txt")
  write_lines(summary_file, capture_output_text({
    cat("scSeuratBuilder Summary\n")
    cat("=====================\n")
    cat("Date:", Sys.Date(), "\n")
    cat("Data directory:", data_dir, "\n")
    cat("Data type:", data_type, "\n")
    cat("Total samples:", n_metadata_samples, "\n")
    cat("Final object:", nrow(merged_obj), "genes x", ncol(merged_obj), "cells\n")
    cat("Gene sync:", ifelse(enable_gene_sync, "Yes", "No"), "\n")
  }))

  cat("  ✓ Object saved:", output_file, "\n")
  cat("\n")

  # ============================================================================
  # Final Summary
  # ============================================================================

  if (verbose) {
    cat("\n")
    cat(paste(rep("=", 70), collapse = ""), "\n")
    cat("✓ Pipeline completed successfully!\n")
    cat(paste(rep("=", 70), collapse = ""), "\n")
    cat("\n")
    cat("Output directory:", output_dir, "\n")
    cat("Seurat object: merged_obj\n")
    cat("  Dimensions:", nrow(merged_obj), "genes x", ncol(merged_obj), "cells\n")
    cat("\n")

    # Show groups
    if ("groups" %in% names(merged_obj@meta.data)) {
      cat("Groups:\n")
      groups_table <- table(merged_obj@meta.data$groups)
      for (g in names(groups_table)) {
        cat(sprintf("  - %s: %d cells\n", g, groups_table[g]))
      }
      cat("\n")
    }
  }

  return(merged_obj)
}


#' Detect Data Structure Type
#'
#' Detects whether data_dir contains:
#' - Type A: Subfolders (each = one sample)
#' - Type B: Direct files
#' - Type C: Single file (merged object)
#'
#' @param data_dir Character. Path to data directory
#' @param verbose Logical. Print progress
#'
#' @return List with type and details
#'
#' @keywords internal
#'
detect_data_structure <- function(data_dir, verbose = TRUE) {

  # Get all items in data_dir
  all_items <- list.files(data_dir)

  # Separate files and directories
  is_dir <- dir.exists(file.path(data_dir, all_items))
  subfolders <- all_items[is_dir]
  files <- all_items[!is_dir]

  # Count by type
  n_subfolders <- length(subfolders)
  n_files <- length(files)

  if (verbose) {
    cat("  Scanning directory...\n")
    cat("    Subfolders:", n_subfolders, "\n")
    cat("    Files:", n_files, "\n")
  }

  # Type C: Single Seurat object file (pre-merged object)
  # Exclude metadata files first
  data_files <- files[!grepl("^metadata\\\\.", files, ignore.case = TRUE)]
  data_files <- data_files[!tolower(tools::file_ext(data_files)) %in% c("xlsx", "xls", "xlsb")]
  n_data_files <- length(data_files)

  # If only 1 data file and it's RDS/RDATA, treat as Type C (merged object)
  if (n_data_files == 1 && n_subfolders == 0) {
    file_ext <- tolower(tools::file_ext(data_files[1]))
    if (file_ext %in% c("rds", "rda", "rdata")) {
      return(list(
        type = "single_file",
        file_path = file.path(data_dir, data_files[1]),
        description = "Single merged object"
      ))
    }
  }

  # Type A: Subfolders (each = one sample)
  if (n_subfolders > 0) {
    # Check if subfolders contain data files
    has_data_in_subfolders <- FALSE
    for (sub in subfolders) {
      sub_path <- file.path(data_dir, sub)
      sub_files <- list.files(sub_path)
      # Check for data files (not metadata)
      data_files <- sub_files[!grepl("^metadata\\\\.", sub_files, ignore.case = TRUE)]
      if (length(data_files) > 0) {
        has_data_in_subfolders <- TRUE
        break
      }
    }

    if (has_data_in_subfolders) {
      return(list(
        type = "subfolders",
        subfolders = subfolders,
        n_samples = n_subfolders,
        description = "Subfolder-per-sample structure"
      ))
    }
  }

  # Type B: Direct files in data_dir
  # (data_files already computed above for Type C check)

  if (n_data_files > 0) {
    return(list(
      type = "direct_files",
      files = data_files,
      n_files = n_data_files,
      description = "Direct files in data_dir"
    ))
  }

  # Fallback
  stop("Unable to detect data structure. Please check your data directory.")
}


#' Match Samples with Metadata (Type A: Subfolders)
#'
#' @keywords internal
#'
match_samples_subfolders <- function(data_dir, metadata_samples, verbose) {

  subfolders <- list.files(data_dir)

  # Exclude potential metadata folders
  subfolders <- subfolders[!grepl("^metadata", subfolders, ignore.case = TRUE)]

  data_samples <- subfolders

  # Match
  only_in_metadata <- setdiff(metadata_samples, data_samples)
  only_in_data <- setdiff(data_samples, metadata_samples)
  common <- intersect(metadata_samples, data_samples)

  if (verbose) {
    cat("  Matching subfolder names with metadata...\n")
    cat("    Metadata samples:", length(metadata_samples), "\n")
    cat("    Subfolders:", length(subfolders), "\n")
    cat("    Matched:", length(common), "\n")
  }

  if (length(only_in_metadata) > 0 || length(only_in_data) > 0) {
    return(list(
      success = FALSE,
      matched_samples = character(),
      missing_in_data = only_in_metadata,
      missing_in_metadata = only_in_data,
      common = common
    ))
  }

  return(list(
    success = TRUE,
    matched_samples = common,
    missing_in_data = character(),
    missing_in_metadata = character(),
    common = common
  ))
}


#' Match Samples with Metadata (Type B: Direct Files)
#'
#' @keywords internal
#'
match_samples_direct_files <- function(data_dir, metadata_samples, verbose) {

  files <- list.files(data_dir)

  # Exclude metadata files (better pattern)
  data_files <- files[!grepl("^metadata\\.xls", files, ignore.case = TRUE)]

  # Also exclude by extension check
  data_files <- data_files[!tolower(tools::file_ext(data_files)) %in% c("xlsx", "xls", "xlsb")]

  # Extract sample names from file names
  file_samples <- extract_sample_names_from_filenames(data_files)

  # Debug
  if (verbose) {
    cat("    All files:", paste(files, collapse = ", "), "\n")
    cat("    Data files:", paste(data_files, collapse = ", "), "\n")
    cat("    Extracted samples:", paste(file_samples, collapse = ", "), "\n")
  }

  # Fuzzy match - pass original file names too
  match_result <- fuzzy_match_samples(metadata_samples, file_samples, data_files)

  if (verbose) {
    cat("  Matching file names with metadata...\n")
    cat("    Metadata samples:", length(metadata_samples), "\n")
    cat("    Data files:", length(data_files), "\n")
    cat("    Matched:", length(match_result$matched), "\n")
  }

  if (!match_result$complete_match) {
    return(list(
      success = FALSE,
      matched_samples = match_result$matched,
      missing_in_data = match_result$unmatched_in_metadata,
      missing_in_metadata = match_result$unmatched_in_files,
      file_mapping = match_result$mapping
    ))
  }

  return(list(
    success = TRUE,
    matched_samples = match_result$matched,
    missing_in_data = character(),
    missing_in_metadata = character(),
    file_mapping = match_result$mapping
  ))
}


#' Match Samples with Metadata (Type C: Single File)
#'
#' @keywords internal
#'
match_samples_single_file <- function(data_dir, metadata_samples, verbose) {

  files <- list.files(data_dir)
  data_file <- files[!grepl("^metadata\\\\.", files, ignore.case = TRUE)][1]

  # Load object and extract orig.ident
  ext <- tolower(tools::file_ext(data_file))

  obj <- if (ext == "rds") {
    readRDS(file.path(data_dir, data_file))
  } else {
    env <- new.env()
    load(file.path(data_dir, data_file), envir = env)
    get(ls(env)[1], envir = env)
  }

  if (!inherits(obj, "Seurat")) {
    stop("File is not a valid Seurat object")
  }

  # Extract sample names from orig.ident
  if ("orig.ident" %in% colnames(obj@meta.data)) {
    data_samples <- as.character(unique(obj@meta.data$orig.ident))
  } else {
    # Try to extract from cell names
    cell_names <- colnames(obj)
    # Extract prefix before first _ or -
    data_samples <- unique(gsub("[-_].*$", "", cell_names))
  }

  # Match
  only_in_metadata <- setdiff(metadata_samples, data_samples)
  only_in_data <- setdiff(data_samples, metadata_samples)
  common <- intersect(metadata_samples, data_samples)

  if (verbose) {
    cat("  Extracting samples from object...\n")
    cat("    Metadata samples:", length(metadata_samples), "\n")
    cat("    Object samples:", length(data_samples), "\n")
    cat("    Matched:", length(common), "\n")
  }

  # For Type C, we allow partial match but warn user
  metadata_matched <- length(only_in_metadata) == 0 && length(only_in_data) == 0

  return(list(
    success = TRUE,  # Always succeed for Type C, user will verify
    matched_samples = common,
    missing_in_data = only_in_metadata,
    missing_in_metadata = only_in_data,
    common = common,
    metadata_matched = metadata_matched,
    seurat_object = obj
  ))
}


#' Extract Sample Names from File Names
#'
#' @keywords internal
#'
extract_sample_names_from_filenames <- function(file_names) {

  # Remove extensions
  base_names <- tools::file_path_sans_ext(file_names)

  # Extract sample names (before first _ or .)
  samples <- gsub("[-_.].*$", "", base_names)

  return(samples)
}


#' Fuzzy Match Samples
#'
#' @keywords internal
#'
fuzzy_match_samples <- function(metadata_samples, file_samples, original_files = NULL) {

  matched <- character()
  unmatched_in_metadata <- character()
  unmatched_in_files <- character()
  mapping <- list()

  # If original file names not provided, use file_samples
  if (is.null(original_files)) {
    original_files <- file_samples
  }

  for (meta_samp in metadata_samples) {
    # Try exact match
    if (meta_samp %in% file_samples) {
      matched <- c(matched, meta_samp)
      # Find the original file name
      file_idx <- which(file_samples == meta_samp)[1]
      mapping[[meta_samp]] <- original_files[file_idx]
      next
    }

    # Try fuzzy match (case insensitive, underscore/space variations)
    lower_meta <- tolower(meta_samp)
    lower_files <- tolower(file_samples)

    # Remove underscores and spaces for comparison
    meta_clean <- gsub("[_ ]", "", lower_meta)
    files_clean <- gsub("[_ ]", "", lower_files)

    match_idx <- which(files_clean == meta_clean)

    if (length(match_idx) > 0) {
      matched <- c(matched, meta_samp)
      # Store the ORIGINAL file name, not the extracted sample name
      mapping[[meta_samp]] <- original_files[match_idx[1]]
    } else {
      unmatched_in_metadata <- c(unmatched_in_metadata, meta_samp)
    }
  }

  # Check for files not matched to any metadata
  all_matched_files <- unlist(mapping)
  unmatched_in_files <- setdiff(original_files, all_matched_files)

  complete_match <- length(unmatched_in_metadata) == 0 &&
                       length(unmatched_in_files) == 0

  return(list(
    matched = matched,
    unmatched_in_metadata = unmatched_in_metadata,
    unmatched_in_files = unmatched_in_files,
    complete_match = complete_match,
    mapping = mapping
  ))
}


#' Validate Metadata File
#'
#' @param file_path Path to metadata file
#' @return List with validation results
#'
#' @export
#'
validate_metadata_file <- function(file_path) {

  if (!file.exists(file_path)) {
    return(list(
      valid = FALSE,
      message = "File not found"
    ))
  }

  # Try to read
  tryCatch({
    metadata <- readxl::read_excel(file_path)
    names(metadata) <- tolower(names(metadata))

    has_samples <- "samples" %in% names(metadata)
    has_groups <- "groups" %in% names(metadata)

    if (!has_samples || !has_groups) {
      return(list(
        valid = FALSE,
        message = "Missing required columns (samples, groups)",
        columns = names(metadata)
      ))
    }

    return(list(
      valid = TRUE,
      message = "Valid metadata file",
      n_samples = nrow(metadata),
      columns = names(metadata)
    ))

  }, error = function(e) {
    return(list(
      valid = FALSE,
      message = paste("Failed to read:", e$message)
    ))
  })
}


#' Capture Output Text
#'
#' @keywords internal
#'
capture_output_text <- function(expr) {
  tmp <- tempfile()
  sink(tmp)
  on.exit(sink())
  expr
  sink()
  tmp_text <- readLines(tmp)
  paste(tmp_text, collapse = "\n")
}


#' Write Lines to File
#'
#' @keywords internal
#'
write_lines <- function(file, ...) {
  text <- capture.output(...)

  # Remove empty lines
  text <- text[nchar(trimws(text)) > 0]

  writeLines(text, file)
}
