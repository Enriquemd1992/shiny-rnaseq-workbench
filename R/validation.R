issue <- function(severity, code, message) data.frame(severity, code, message, stringsAsFactors = FALSE)

validate_counts <- function(counts) {
  out <- list()
  add <- function(x) out[[length(out) + 1L]] <<- x
  if (anyNA(counts)) add(issue("error", "missing_counts", "Counts contain missing or non-numeric values."))
  if (any(counts < 0, na.rm = TRUE)) add(issue("error", "negative_counts", "Counts must be non-negative."))
  if (any(counts != floor(counts), na.rm = TRUE)) add(issue("error", "non_integer_counts", "Counts must be integer-valued raw counts."))
  if (any(!nzchar(rownames(counts)))) add(issue("error", "empty_gene_id", "Some gene identifiers are empty."))
  if (anyDuplicated(rownames(counts))) add(issue("error", "duplicate_genes", "Gene identifiers must be unique."))
  if (anyDuplicated(colnames(counts))) add(issue("error", "duplicate_samples_counts", "Sample names in counts must be unique."))
  if (!length(out)) return(data.frame(severity=character(), code=character(), message=character()))
  do.call(rbind, out)
}

validate_metadata <- function(metadata) {
  out <- list()
  add <- function(x) out[[length(out) + 1L]] <<- x
  if (any(!nzchar(rownames(metadata)))) add(issue("error", "empty_sample_id", "Some sample identifiers are empty."))
  if (anyDuplicated(rownames(metadata))) add(issue("error", "duplicate_samples_metadata", "Sample identifiers in metadata must be unique."))
  if (any(!nzchar(names(metadata)))) add(issue("error", "empty_column_name", "Metadata column names cannot be empty."))
  if (anyDuplicated(names(metadata))) add(issue("error", "duplicate_column_names", "Metadata column names must be unique."))
  for (nm in names(metadata)) {
    x <- metadata[[nm]]
    observed <- unique(x[!is.na(x) & trimws(as.character(x)) != ""])
    if (anyNA(x) || any(trimws(as.character(x)) == "", na.rm = TRUE)) add(issue("warning", "missing_metadata", sprintf("Metadata column '%s' contains missing or empty values.", nm)))
    if (length(observed) <= 1L) add(issue("warning", "constant_column", sprintf("Metadata column '%s' is constant or empty.", nm)))
    if (is.factor(x) && any(table(x) == 0)) add(issue("warning", "empty_levels", sprintf("Factor '%s' contains unused levels.", nm)))
    if ((is.factor(x) || is.character(x)) && any(table(x, useNA = "no") < 2L)) add(issue("warning", "low_replication", sprintf("Column '%s' has a group with fewer than 2 samples.", nm)))
  }
  if (!length(out)) return(data.frame(severity=character(), code=character(), message=character()))
  do.call(rbind, out)
}

validate_sample_match <- function(counts, metadata) {
  count_ids <- colnames(counts); meta_ids <- rownames(metadata)
  if (identical(count_ids, meta_ids)) return(issue("ok", "samples_aligned", "Samples match and are in the same order."))
  if (setequal(count_ids, meta_ids)) return(issue("warning", "sample_order", "Sample sets match, but metadata order differs. Use explicit reordering."))
  missing_meta <- setdiff(count_ids, meta_ids); extra_meta <- setdiff(meta_ids, count_ids)
  issue("error", "sample_mismatch", paste0("Sample sets differ. Missing in metadata: ", paste(missing_meta, collapse=", "), "; absent from counts: ", paste(extra_meta, collapse=", "), "."))
}

validate_annotations <- function(annotations) {
  if (is.null(annotations) || !ncol(annotations)) return(data.frame(severity=character(), code=character(), message=character()))
  out <- list()
  for (nm in names(annotations)) {
    x <- as.character(annotations[[nm]])
    present <- !is.na(x) & nzchar(trimws(x))
    duplicated_values <- unique(x[present & (duplicated(x) | duplicated(x, fromLast = TRUE))])
    identifier_like <- sum(present) > 0 && length(unique(x[present])) / sum(present) >= 0.5
    if (length(duplicated_values) && identifier_like) out[[length(out)+1L]] <- issue("warning", "duplicate_annotation", sprintf("Annotation '%s' has %s duplicated value(s). They will be preserved and never collapsed automatically.", nm, format(length(duplicated_values), big.mark=",")))
  }
  if (!length(out)) return(data.frame(severity=character(), code=character(), message=character()))
  do.call(rbind,out)
}

validate_all <- function(counts, metadata, annotations = NULL) rbind(validate_counts(counts), validate_metadata(metadata), validate_sample_match(counts, metadata), validate_annotations(annotations))
has_errors <- function(issues) any(issues$severity == "error")

align_metadata <- function(counts, metadata) {
  if (!setequal(colnames(counts), rownames(metadata))) stop("Samples cannot be reordered because the sets differ.", call. = FALSE)
  metadata[colnames(counts), , drop = FALSE]
}

validate_design <- function(metadata, variables) {
  if (!length(variables)) return(issue("error", "empty_design", "Select at least one design variable."))
  if (any(!variables %in% names(metadata))) return(issue("error", "unknown_variable", "The design contains an unknown metadata column."))
  formula <- stats::reformulate(variables)
  mm <- try(stats::model.matrix(formula, metadata), silent = TRUE)
  if (inherits(mm, "try-error")) return(issue("error", "invalid_design", "The design matrix could not be constructed; check missing values and variable types."))
  if (qr(mm)$rank < ncol(mm)) return(issue("error", "rank_deficient", "The design is not full rank. Variables may be constant or completely confounded."))
  issue("ok", "valid_design", sprintf("Design %s is full rank (%d coefficients).", paste(deparse(formula), collapse=""), ncol(mm)))
}
