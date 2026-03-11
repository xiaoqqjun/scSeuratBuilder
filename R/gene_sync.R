#' Auto Synchronize Genes for a List of Seurat Objects
#'
#' Synchronizes gene names using geneSync to eliminate Cell Ranger version differences.
#' IMPORTANT: Always sync when enable_gene_sync = TRUE, even if genes appear to be symbols!
#'
#' @param seurat_list List of Seurat objects
#' @param species Character. Species (e.g., "homo", "mus")
#' @param verbose Logical. Print progress
#'
#' @return Updated list of Seurat objects
#'
#' @keywords internal
#'
sync_genes_for_list <- function(seurat_list, species, verbose) {

  if (!requireNamespace("geneSync", quietly = TRUE)) {
    stop("Package 'geneSync' is required. Install from local path or use:\n",
         "  devtools::install_local('path/to/geneSync_V3.2.0')")
  }

  cat("  Gene synchronization purpose:\n")
  cat("    - Eliminate Cell Ranger version differences\n")
  cat("    - Unify gene names across samples processed with different versions\n")
  cat("    - e.g., 'Gene-old' vs 'Gene-new' → both become 'Gene'\n")
  cat("    - Standardize to authority symbols when available\n")
  cat("\n")

  # Map user input species to geneSync species names
  species_mapping <- c(
    "homo" = "homo",
    "human" = "homo",
    "hsa" = "homo",
    "mus" = "mmu",
    "mouse" = "mmu",
    "mmu" = "mmu",
    "rat" = "rat",
    "rattus" = "rat",
    "rno" = "rat"
  )

  species_lower <- tolower(species)
  if (!species_lower %in% names(species_mapping)) {
    stop("Unknown species: ", species, "\n",
         "Supported: homo/human/hsa, mus/mouse/mmu, rat/rattus/rno")
  }

  geneSync_species <- species_mapping[species_lower]

  # Detect gene type from first object to determine 'from' parameter
  first_obj <- seurat_list[[1]]
  gene_type <- detect_gene_id_type(first_obj, verbose = verbose)

  # Determine 'from' parameter for geneSync
  from_type <- switch(gene_type,
    "ensembl_id_version" = "ensembl_id",
    "ensembl_id" = "ensembl_id",
    "entrez_id" = "gene_id",
    "symbol" = "symbol"
  )

  cat("  Species:", species, "-> geneSync:", geneSync_species, "\n")
  cat("  From:", from_type, "-> To: symbol (final_symbol)\n")
  cat("\n")

  # Sync each object
  for (name in names(seurat_list)) {
    cat("  [", which(names(seurat_list) == name), "/", length(seurat_list), "] ",
        name, "... ", sep = "")

    tryCatch({
      # Get gene count before sync
      n_before <- nrow(seurat_list[[name]])

      # Call geneSync with correct parameters
      seurat_list[[name]] <- geneSync::gene_sync_add_to_obj(
        seurat_list[[name]],
        species = geneSync_species,
        from = from_type,
        to = "symbol",
        symbol_type = "final_symbol",  # Use authority symbols when available
        update_features = TRUE,
        verbose = FALSE
      )

      # Get gene count after sync
      n_after <- nrow(seurat_list[[name]])

      # Show change in gene count
      if (n_after != n_before) {
        cat(sprintf("OK (%d -> %d genes)\n", n_before, n_after))
      } else {
        cat(sprintf("OK (%d genes)\n", n_after))
      }

    }, error = function(e) {
      cat("ERROR:", conditionMessage(e), "\n")
    })
  }

  cat("\n")
  return(seurat_list)
}


#' Detect Gene ID Type from Seurat Object
#'
#' @param obj Seurat object
#' @param verbose Logical. Print details
#'
#' @return Character: "symbol", "ensembl_id", "ensembl_id_version", "entrez_id", or "unknown"
#'
#' @keywords internal
#'
detect_gene_id_type <- function(obj, verbose = TRUE) {

  gene_names <- rownames(obj)
  n_genes <- length(gene_names)

  # Sample for faster detection
  sample_size <- min(1000, n_genes)
  sample_genes <- sample(gene_names, sample_size)

  # Patterns
  ensembl_version_pattern <- "^ENS[GTP]\\d{11}\\.\\d+$"
  ensembl_pattern <- "^ENS[GTP]\\d{11}$"
  entrez_pattern <- "^\\d+$"

  n_total <- length(sample_genes)

  n_ensembl_version <- sum(grepl(ensembl_version_pattern, sample_genes))
  n_ensembl <- sum(grepl(ensembl_pattern, sample_genes))
  n_entrez <- sum(grepl(entrez_pattern, sample_genes))

  prop_ensembl_version <- n_ensembl_version / n_total
  prop_ensembl <- n_ensembl / n_total
  prop_entrez <- n_entrez / n_total

  threshold <- 0.8

  if (verbose) {
    cat("    Gene ID type detection:\n")
    cat("      Total genes:", n_genes, "\n")
    cat("      Sample size:", sample_size, "\n")
    cat("      Ensembl ID (with version):", n_ensembl_version, sprintf("(%.1f%%)", prop_ensembl_version * 100), "\n")
    cat("      Ensembl ID (no version):", n_ensembl, sprintf("(%.1f%%)", prop_ensembl * 100), "\n")
    cat("      Entrez ID (numbers only):", n_entrez, sprintf("(%.1f%%)", prop_entrez * 100), "\n")
  }

  if (prop_ensembl_version >= threshold) {
    if (verbose) cat("    → Detected: ensembl_id_version\n")
    return("ensembl_id_version")
  } else if (prop_ensembl >= threshold) {
    if (verbose) cat("    → Detected: ensembl_id\n")
    return("ensembl_id")
  } else if (prop_entrez >= threshold) {
    if (verbose) cat("    → Detected: entrez_id\n")
    return("entrez_id")
  } else {
    if (verbose) cat("    → Detected: symbol (default)\n")
    return("symbol")
  }
}
