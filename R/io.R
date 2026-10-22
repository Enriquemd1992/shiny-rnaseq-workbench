allowed_extensions <- c("csv", "tsv", "txt", "xlsx")

safe_extension <- function(path, original_name = basename(path)) {
  ext <- tolower(tools::file_ext(original_name))
  if (!ext %in% allowed_extensions) {
    stop("Unsupported file type. Use CSV, TSV/TXT, or XLSX.", call. = FALSE)
  }
  ext
}

read_tabular_file <- function(path, original_name = basename(path), max_bytes = 100 * 1024^2) {
  info <- file.info(path)
  if (is.na(info$size) || info$size > max_bytes) {
    stop("The file is missing or exceeds the 100 MB limit.", call. = FALSE)
  }
  ext <- safe_extension(path, original_name)
  out <- switch(ext,
    csv = utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE),
    tsv = utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE),
    txt = utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE),
    xlsx = as.data.frame(readxl::read_excel(path, .name_repair = "minimal"))
  )
  if (!nrow(out) || ncol(out) < 2L) stop("The file must contain rows and at least two columns.", call. = FALSE)
  out
}

prepare_counts <- function(data, gene_column) {
  prepare_count_data(data, gene_column, annotation_columns = character())$counts
}

is_numeric_column <- function(x) {
  converted <- suppressWarnings(as.numeric(x))
  all(is.na(x) | trimws(as.character(x)) == "" | !is.na(converted))
}

suggest_annotation_columns <- function(data, gene_column) {
  candidates <- setdiff(names(data), gene_column)
  candidates[!vapply(data[candidates], is_numeric_column, logical(1))]
}

prepare_count_data <- function(data, gene_column, annotation_columns = character()) {
  if (!gene_column %in% names(data)) stop("Select a valid gene identifier column.", call. = FALSE)
  if (any(!annotation_columns %in% setdiff(names(data), gene_column))) stop("Select valid annotation columns.", call. = FALSE)
  genes <- trimws(as.character(data[[gene_column]]))
  if (any(!nzchar(genes) | is.na(genes))) stop(sprintf("The selected gene identifier column '%s' contains empty values.", gene_column), call. = FALSE)
  if (anyDuplicated(genes)) {
    duplicated_ids <- sum(duplicated(genes))
    stop(sprintf("The selected gene identifier column '%s' is not unique (%s duplicated rows). For annotated matrices, use a unique key such as ensembl_id and keep gene_symbol as an annotation.", gene_column, format(duplicated_ids, big.mark=",")), call. = FALSE)
  }
  values <- data[, setdiff(names(data), c(gene_column, annotation_columns)), drop = FALSE]
  if (!ncol(values)) stop("No sample count columns remain after excluding annotations.", call. = FALSE)
  numeric_values <- lapply(values, function(x) suppressWarnings(as.numeric(x)))
  matrix_values <- as.matrix(as.data.frame(numeric_values, check.names = FALSE))
  rownames(matrix_values) <- genes
  storage.mode(matrix_values) <- "numeric"
  annotations <- data[, annotation_columns, drop = FALSE]
  rownames(annotations) <- genes
  list(counts = matrix_values, annotations = annotations, gene_id_column = gene_column, annotation_columns = annotation_columns)
}

prepare_metadata <- function(data, sample_column) {
  if (!sample_column %in% names(data)) stop("Select a valid sample identifier column.", call. = FALSE)
  samples <- trimws(as.character(data[[sample_column]]))
  out <- data[, setdiff(names(data), sample_column), drop = FALSE]
  rownames(out) <- samples
  out
}

select_sample_ids <- function(metadata, mode = "all", variable = NULL, values = character(), manual_ids = character()) {
  ids <- rownames(metadata)
  selected <- switch(mode,
    all = ids,
    metadata = {
      if (is.null(variable) || !variable %in% names(metadata)) stop("Select a valid metadata variable for sample filtering.",call.=FALSE)
      if (!length(values)) stop("Select at least one metadata value.",call.=FALSE)
      ids[as.character(metadata[[variable]]) %in% as.character(values)]
    },
    manual = {
      if (!length(manual_ids)) stop("Select at least one sample.",call.=FALSE)
      intersect(ids,manual_ids)
    },
    stop("Unknown sample-selection mode.",call.=FALSE)
  )
  selected <- ids[ids %in% selected]
  if (length(selected) < 2L) stop("Sample selection must retain at least two samples.",call.=FALSE)
  selected
}

subset_samples <- function(counts, metadata, sample_ids) {
  sample_ids <- colnames(counts)[colnames(counts) %in% sample_ids]
  if (length(sample_ids) < 2L) stop("At least two matched samples are required.",call.=FALSE)
  if (!all(sample_ids %in% rownames(metadata))) stop("Selected samples are missing from metadata.",call.=FALSE)
  list(counts=counts[,sample_ids,drop=FALSE],metadata=metadata[sample_ids,,drop=FALSE])
}
