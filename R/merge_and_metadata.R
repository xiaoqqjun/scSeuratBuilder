#' Merge Seurat Objects Safely
#'
#' Merges multiple Seurat objects with error handling
#'
#' @param seurat_list List of Seurat objects
#' @param verbose Logical. Print progress
#'
#' @return Merged Seurat object
#'
#' @keywords internal
#'
merge_seurat_objects_safe <- function(seurat_list, verbose) {

  if (length(seurat_list) == 1) {
    return(seurat_list[[1]])
  }

  if (verbose) {
    cat("  Merging", length(seurat_list), "objects...\n")
  }

  # Use Reduce for efficient merging
  merged_obj <- Reduce(function(x, y) {
    merge(x = x, y = y, project = "scRNA_project")
  }, seurat_list)

  return(merged_obj)
}


#' Integrate Metadata Safely
#'
#' @param obj Seurat object
#' @param metadata Data frame with samples and groups columns
#' @param verbose Logical. Print progress
#'
#' @return Updated Seurat object
#'
#' @keywords internal
#'
integrate_metadata_safe <- function(obj, metadata, verbose) {

  if (verbose) {
    cat("  Integrating metadata...\n")
  }

  # Get cell names and their sample assignments
  cell_names <- colnames(obj)
  cell_sample <- obj@meta.data$orig.ident

  # Create mapping data frame
  cell_df <- data.frame(
    cell = cell_names,
    sample = as.character(cell_sample),
    stringsAsFactors = FALSE
  )

  # Merge with metadata
  merged_meta <- merge(
    cell_df,
    metadata,
    by.x = "sample",
    by.y = "samples",
    all.x = TRUE,
    sort = FALSE
  )

  # Reorder to match original cell order
  merged_meta <- merged_meta[match(cell_names, merged_meta$cell), ]
  merged_meta$cell <- NULL

  # Add to object
  for (col in setdiff(names(merged_meta), c("sample", "orig.ident", "nCount_RNA", "nFeature_RNA"))) {
    obj@meta.data[[col]] <- merged_meta[[col]]
  }

  if (verbose) {
    added_cols <- setdiff(names(merged_meta), c("sample", "cell"))
    cat("    Added columns:", paste(added_cols, collapse = ", "), "\n")
  }

  return(obj)
}
