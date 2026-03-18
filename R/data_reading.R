#' Read Data from Subfolders (Type A)
#'
#' @keywords internal
#'
read_data_from_subfolders <- function(data_dir, sample_names,
                                      min_cells, min_features, verbose) {

  seurat_list <- list()

  for (i in seq_along(sample_names)) {
    sample_name <- sample_names[i]
    sample_path <- file.path(data_dir, sample_name)

    cat(sprintf("  [%d/%d] Reading: %s\n", i, length(sample_names), sample_name))

    tryCatch({
      # Detect format in subfolder
      format_info <- detect_folder_format(sample_path)

      obj <- switch(format_info$format,
        "threebrothers" = read_threebrothers(sample_path, sample_name,
                                               min_cells, min_features),
        "h5" = read_h5_file(sample_path, sample_name,
                             min_cells, min_features),
        "rds" = read_seurat_file(sample_path),
        "rdata" = read_seurat_file(sample_path)
      )

      if (!is.null(obj)) {
        seurat_list[[sample_name]] <- obj
        cat("    ✓ Loaded:", nrow(obj), "genes x", ncol(obj), "cells\n")
      }

    }, error = function(e) {
      cat("    ✗ Error:", e$message, "\n")
    })
  }

  return(seurat_list)
}


#' Read Data from Direct Files (Type B)
#'
#' @keywords internal
#'
read_data_from_direct_files <- function(data_dir, sample_names, file_mapping,
                                         min_cells, min_features, verbose) {

  seurat_list <- list()

  for (i in seq_along(sample_names)) {
    sample_name <- sample_names[i]
    file_name <- file_mapping[[sample_name]]
    file_path <- file.path(data_dir, file_name)

    cat(sprintf("  [%d/%d] Reading: %s\n", i, length(sample_names), file_name))

    tryCatch({
      # Detect format from file extension
      ext <- tolower(tools::file_ext(file_name))

      cat(sprintf("    File type: %s\n", ext))
      cat(sprintf("    File path: %s\n", file_path))

      obj <- switch(ext,
        "h5" = read_h5_file_direct(file_path, sample_name, min_cells, min_features),
        "rds" = read_seurat_file(file_path),
        "rda" = read_seurat_file(file_path),
        "rdata" = read_seurat_file(file_path),
        "mtx" = read_threebrothers_from_file(file_path, sample_name,
                                              file_name, min_cells, min_features)
      )

      if (!is.null(obj) && inherits(obj, "Seurat")) {
        seurat_list[[sample_name]] <- obj
        cat(sprintf("    ✓ Loaded: %d genes x %d cells\n", nrow(obj), ncol(obj)))
      } else {
        cat("    ✗ Failed to create Seurat object\n")
      }

    }, error = function(e) {
      cat(sprintf("    ✗ Error: %s\n", e$message))
      cat(sprintf("    File: %s\n", file_path))
    })
  }

  return(seurat_list)
}


#' Read Merged Object (Type C)
#'
#' @keywords internal
#'
read_merged_object <- function(data_dir, verbose) {

  files <- list.files(data_dir)
  data_file <- files[!grepl("^metadata\\\\.", files, ignore.case = TRUE)][1]
  file_path <- file.path(data_dir, data_file)

  cat("  Loading merged object:", data_file, "\n")

  ext <- tolower(tools::file_ext(data_file))

  obj <- if (ext == "rds") {
    readRDS(file_path)
  } else if (ext %in% c("rda", "rdata")) {
    env <- new.env()
    load(file_path, envir = env)
    get(ls(env)[1], envir = env)
  } else {
    stop("Unsupported file format: ", ext)
  }

  if (!inherits(obj, "Seurat")) {
    stop("File is not a valid Seurat object")
  }

  cat("    ✓ Loaded:", nrow(obj), "genes x", ncol(obj), "cells\n")

  return(list(merged = obj))
}


#' Detect Folder Format
#'
#' @keywords internal
#'
detect_folder_format <- function(folder_path) {

  files <- list.files(folder_path)

  # Check for H5
  if (any(grepl("\\.h5$", files, ignore.case = TRUE))) {
    return(list(format = "h5"))
  }

  # Check for RDS/RDATA
  if (any(grepl("\\.(rds|rda|rdata)$", files, ignore.case = TRUE))) {
    ext <- tolower(tools::file_ext(files[grepl("\\.(rds|rda|rdata)$", files, ignore.case = TRUE)][1]))
    return(list(format = if (ext == "rds") "rds" else "rdata"))
  }

  # Check for three brothers (MTX)
  has_barcodes <- any(grepl("barcode", files, ignore.case = TRUE))
  has_features <- any(grepl("feature|gene", files, ignore.case = TRUE))
  has_matrix <- any(grepl("matrix", files, ignore.case = TRUE))

  if (has_barcodes && has_features && has_matrix) {
    return(list(format = "threebrothers"))
  }

  # Check for standard 10X directories
  subdirs <- list.files(folder_path)
  if (any(subdirs %in% c("filtered_feature_bc_matrix", "raw_feature_bc_matrix"))) {
    return(list(format = "threebrothers"))
  }

  stop("Unable to detect data format in folder: ", folder_path)
}


#' Read Three Brothers Format (MTX)
#'
#' @keywords internal
#'
read_threebrothers <- function(folder_path, sample_name, min_cells, min_features) {

  # Check for standard 10X subdirectories
  possible_dirs <- c("filtered_feature_bc_matrix", "raw_feature_bc_matrix")

  for (subdir in possible_dirs) {
    full_path <- file.path(folder_path, subdir)
    if (dir.exists(full_path)) {
      counts <- Seurat::Read10X(data.dir = full_path)

      # Suppress warning about non-unique features
      obj <- suppressWarnings({
        Seurat::CreateSeuratObject(
          counts = counts,
          project = sample_name,
          min.cells = min_cells,
          min.features = min_features
        )
      })

      return(validate_and_fix_orig_ident(obj, sample_name))
    }
  }

  # Try reading with flexible detection
  counts <- read_with_flexible_detection(folder_path)

  # Suppress warning about non-unique features
  obj <- suppressWarnings({
    Seurat::CreateSeuratObject(
      counts = counts,
      project = sample_name,
      min.cells = min_cells,
      min.features = min_features
    )
  })

  # Validate and fix orig.ident if needed
  return(validate_and_fix_orig_ident(obj, sample_name))
}


#' Read Three Brothers from File
#'
#' @keywords internal
#'
read_threebrothers_from_file <- function(file_path, sample_name, original_filename,
                                          min_cells, min_features) {

  # Get directory and detect other files
  dir_path <- dirname(file_path)

  counts <- read_with_flexible_detection(dir_path)

  # Suppress warning about non-unique features
  obj <- suppressWarnings({
    Seurat::CreateSeuratObject(
      counts = counts,
      project = sample_name,
      min.cells = min_cells,
      min.features = min_features
    )
  })

  # Validate and fix orig.ident if needed
  return(validate_and_fix_orig_ident(obj, sample_name))
}


#' Read with Flexible Detection
#'
#' @keywords internal
#'
read_with_flexible_detection <- function(data_path) {

  files <- list.files(data_path)

  # Find files by pattern
  barcode_file <- grep("barcode", files, value = TRUE, ignore.case = TRUE)[1]
  feature_file <- grep("feature|gene", files, value = TRUE, ignore.case = TRUE)[1]
  matrix_file <- grep("matrix", files, value = TRUE, ignore.case = TRUE)[1]

  if (is.na(barcode_file) || is.na(feature_file) || is.na(matrix_file)) {
    stop("Could not identify all required files (barcodes, features, matrix)")
  }

  barcode_path <- file.path(data_path, barcode_file)
  feature_path <- file.path(data_path, feature_file)
  matrix_path <- file.path(data_path, matrix_file)

  # Read files
  barcodes <- tryCatch({
    data.table::fread(barcode_path, header = FALSE)$V1
  }, error = function(e) {
    read.delim(barcode_path, header = FALSE)[,1]
  })

  features <- tryCatch({
    data.table::fread(feature_path, header = FALSE)
  }, error = function(e) {
    read.delim(feature_path, header = FALSE)
  })

  mat <- Matrix::readMM(matrix_path)

  # Set row and column names
  if (ncol(features) >= 2) {
    rownames(mat) <- features[, 2]
  } else {
    rownames(mat) <- features[, 1]
  }

  colnames(mat) <- barcodes

  # Detect non-standard barcode formats
  if (any(grepl("_", barcodes))) {
    n_with_underscore <- sum(grepl("_", barcodes))
    n_total <- length(barcodes)

    message(sprintf("  ⚠ Warning: Non-standard barcode format detected"))
    message(sprintf("    - %d out of %d barcodes contain underscores", n_with_underscore, n_total))
    message(sprintf("    - Sample may use custom barcode format (e.g., multi-part barcodes)"))
    message(sprintf("    - This is common in datasets from certain platforms (e.g., some GEO datasets)"))
    message(sprintf("    - orig.ident will be set to the sample name to ensure correct grouping"))

    # Analyze barcode structure for debugging
    sample_barcodes <- head(barcodes, min(5, n_total))
    parts_list <- strsplit(sample_barcodes, "_")
    n_parts <- unique(sapply(parts_list, length))

    if (length(n_parts) == 1 && n_parts > 1) {
      message(sprintf("    - Barcode structure: %d segments (separated by '_')", n_parts))
    }
  }

  return(mat)
}


#' Validate and Fix orig.ident
#'
#' Ensures that orig.ident is correctly set to the sample name.
#' This is crucial for samples with non-standard barcode formats (e.g., multi-part barcodes).
#'
#' @param obj Seurat object
#' @param expected_sample_name Character. Expected sample name for orig.ident
#'
#' @return Seurat object with corrected orig.ident
#'
#' @keywords internal
#'
validate_and_fix_orig_ident <- function(obj, expected_sample_name) {

  # Check if orig.ident matches expected sample name
  current_orig_idents <- unique(obj@meta.data$orig.ident)

  if (length(current_orig_idents) == 1 && current_orig_idents[1] == expected_sample_name) {
    # Everything is correct
    return(obj)
  }

  # orig.ident is not correct - fix it
  if (length(current_orig_idents) > 1) {
    message(sprintf("  ⚠ Warning: orig.ident contains %d different values instead of 1", length(current_orig_idents)))
    message(sprintf("    Expected: '%s'", expected_sample_name))
    message(sprintf("    Found: %s", paste(head(current_orig_idents, 3), collapse = ", ")))
    message(sprintf("    Fixing: Setting all cells to orig.ident = '%s'", expected_sample_name))

  } else if (current_orig_idents[1] != expected_sample_name) {
    message(sprintf("  ⚠ Warning: orig.ident is '%s' but expected '%s'", current_orig_idents[1], expected_sample_name))
    message(sprintf("    Fixing: Setting orig.ident to '%s'", expected_sample_name))
  }

  # Fix orig.ident
  obj@meta.data$orig.ident <- expected_sample_name

  # Verify the fix
  if (length(unique(obj@meta.data$orig.ident)) == 1 && unique(obj@meta.data$orig.ident) == expected_sample_name) {
    message(sprintf("    ✓ Fixed: All %d cells now have orig.ident = '%s'", ncol(obj), expected_sample_name))
  }

  return(obj)
}


#' Read H5 File
#'
#' @keywords internal
#'
read_h5_file <- function(folder_path, sample_name, min_cells, min_features) {

  files <- list.files(folder_path)
  h5_files <- files[grepl("\\.h5$", files, ignore.case = TRUE)]

  if (length(h5_files) == 0) {
    stop("No H5 file found in folder")
  }

  h5_path <- file.path(folder_path, h5_files[1])

  counts <- Seurat::Read10X_h5(filename = h5_path)

  return(Seurat::CreateSeuratObject(
    counts = counts,
    project = sample_name,
    min.cells = min_cells,
    min.features = min_features
  ))
}


#' Read H5 File Direct
#'
#' @keywords internal
#'
read_h5_file_direct <- function(file_path, sample_name, min_cells, min_features) {

  cat(sprintf("    Reading H5 file: %s\n", basename(file_path)))

  # Check if file exists
  if (!file.exists(file_path)) {
    stop("File does not exist: ", file_path)
  }

  # Read H5 file
  counts <- Seurat::Read10X_h5(filename = file_path)

  # Handle multiple modalities (Gene Expression + others)
  if (is.list(counts)) {
    cat("    Multiple modalities detected, using 'Gene Expression' data\n")
    if ("Gene Expression" %in% names(counts)) {
      counts <- counts[["Gene Expression"]]
    } else if ("RNA" %in% names(counts)) {
      counts <- counts[["RNA"]]
    } else {
      # Use first available modality
      modalities <- names(counts)
      cat(sprintf("    Available modalities: %s\n", paste(modalities, collapse = ", ")))
      cat(sprintf("    Using first modality: %s\n", modalities[1]))
      counts <- counts[[1]]
    }
  }

  # Check if counts is valid
  if (!is.matrix(counts) && !inherits(counts, "dgCMatrix")) {
    stop("Unable to extract valid count matrix from H5 file")
  }

  cat(sprintf("    Counts: %d genes x %d cells\n", nrow(counts), ncol(counts)))

  obj <- Seurat::CreateSeuratObject(
    counts = counts,
    project = sample_name,
    min.cells = min_cells,
    min.features = min_features
  )

  return(obj)
}


#' Read Seurat File (RDS/RDATA)
#'
#' @keywords internal
#'
read_seurat_file <- function(file_path) {

  ext <- tolower(tools::file_ext(file_path))

  if (ext == "rds") {
    return(readRDS(file_path))
  } else {
    env <- new.env()
    load(file_path, envir = env)
    return(get(ls(env)[1], envir = env))
  }
}
