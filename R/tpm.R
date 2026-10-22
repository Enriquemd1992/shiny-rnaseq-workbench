prepare_gene_lengths <- function(data, gene_id_column, length_column) {
  if (!all(c(gene_id_column, length_column) %in% names(data))) stop("Select valid gene ID and length columns.", call. = FALSE)
  ids <- trimws(as.character(data[[gene_id_column]]))
  lengths <- suppressWarnings(as.numeric(data[[length_column]]))
  if (any(!nzchar(ids) | is.na(ids))) stop("Gene length identifiers contain missing values.", call. = FALSE)
  if (anyDuplicated(ids)) stop("Gene length identifiers must be unique.", call. = FALSE)
  if (any(!is.finite(lengths) | lengths <= 0)) stop("Gene lengths must be positive finite values in base pairs.", call. = FALSE)
  stats::setNames(lengths, ids)
}

compute_tpm_export <- function(counts, annotations = NULL, gene_lengths,
                               min_count = 10, min_samples = 2,
                               gene_types = character()) {
  if (!is.numeric(min_count) || min_count < 0 || min_samples < 1 || min_samples > ncol(counts)) stop("Invalid expression filter.", call. = FALSE)
  matched <- match(rownames(counts), names(gene_lengths))
  if (anyNA(matched)) stop(sprintf("Gene lengths are missing for %d count-matrix gene(s). Use lengths from the matching annotation release.", sum(is.na(matched))), call. = FALSE)
  lengths <- as.numeric(gene_lengths[matched])
  keep <- rowSums(counts >= min_count) >= min_samples
  if (length(gene_types)) {
    if (is.null(annotations) || !"gene_type" %in% names(annotations)) stop("gene_type annotation is required for this filter.", call. = FALSE)
    keep <- keep & as.character(annotations[rownames(counts), "gene_type"]) %in% gene_types
  }
  if (!any(keep)) stop("Filters retain no genes.", call. = FALSE)
  filtered_counts <- counts[keep, , drop = FALSE]
  filtered_lengths <- lengths[keep]
  rpk <- filtered_counts / (filtered_lengths / 1000)
  scaling <- colSums(rpk) / 1e6
  if (any(!is.finite(scaling) | scaling <= 0)) stop("TPM cannot be calculated for a sample with zero total RPK.", call. = FALSE)
  tpm <- sweep(rpk, 2, scaling, "/")
  ann <- if (is.null(annotations)) data.frame(row.names = rownames(counts)) else annotations[rownames(counts), , drop = FALSE]
  list(counts = filtered_counts, tpm = tpm,
    annotations = ann[keep, , drop = FALSE], lengths_bp = stats::setNames(filtered_lengths, rownames(filtered_counts)),
    parameters = list(min_count = min_count, min_samples = min_samples,
      gene_types = if (length(gene_types)) paste(gene_types, collapse = ", ") else "All",
      genes_input = nrow(counts), genes_retained = sum(keep), samples = paste(colnames(counts), collapse = ", ")))
}

tpm_export_table <- function(result, matrix_name = c("tpm", "counts")) {
  matrix_name <- match.arg(matrix_name); matrix <- result[[matrix_name]]
  data.frame(gene_id = rownames(matrix), result$annotations, length_bp = unname(result$lengths_bp[rownames(matrix)]),
    matrix, row.names = NULL, check.names = FALSE)
}
