`%||%` <- function(x, y) if (is.null(x)) y else x

assert_metadata_column <- function(metadata, column) {
  if (length(column) != 1L || !nzchar(column) || !column %in% names(metadata)) stop("Select a valid metadata column.", call. = FALSE)
  invisible(TRUE)
}

validate_new_column_name <- function(metadata, name, current = NULL) {
  name <- trimws(name)
  if (!nzchar(name)) stop("Column names cannot be empty.", call. = FALSE)
  if (name %in% setdiff(names(metadata), current %||% character())) stop(sprintf("A column named '%s' already exists.", name), call. = FALSE)
  name
}

convert_metadata_column <- function(metadata, column, type) {
  assert_metadata_column(metadata, column); old <- metadata[[column]]
  converted <- switch(type,
    character = as.character(old),
    factor = factor(as.character(old), levels = unique(as.character(old)[!is.na(old) & nzchar(trimws(as.character(old)))])),
    numeric = suppressWarnings(as.numeric(as.character(old))),
    stop("Unsupported metadata type.", call. = FALSE)
  )
  if (identical(type, "numeric")) {
    newly_missing <- is.na(converted) & !is.na(old) & nzchar(trimws(as.character(old)))
    if (any(newly_missing)) stop(sprintf("Numeric conversion would lose these values: %s.", paste(utils::head(unique(as.character(old)[newly_missing]), 5), collapse = ", ")), call. = FALSE)
  }
  metadata[[column]] <- converted; metadata
}

rename_metadata_column <- function(metadata, column, new_name) {
  assert_metadata_column(metadata, column); new_name <- validate_new_column_name(metadata, new_name, current = column)
  names(metadata)[match(column, names(metadata))] <- new_name; metadata
}

add_metadata_column <- function(metadata, name, value = "") {
  name <- validate_new_column_name(metadata, name); metadata[[name]] <- rep(value, nrow(metadata)); metadata
}

delete_metadata_column <- function(metadata, column) {
  assert_metadata_column(metadata, column)
  if (ncol(metadata) <= 1L) stop("Keep at least one metadata variable.", call. = FALSE)
  metadata[[column]] <- NULL; metadata
}

edit_metadata_cell <- function(metadata, row, column, value) {
  assert_metadata_column(metadata, column)
  if (length(row) != 1L || row < 1L || row > nrow(metadata)) stop("Invalid metadata row.", call. = FALSE)
  x <- metadata[[column]]
  if (is.factor(x)) {
    value <- as.character(value); if (!value %in% levels(x)) levels(x) <- c(levels(x), value); x[row] <- value
  } else if (is.numeric(x)) {
    parsed <- suppressWarnings(as.numeric(value))
    if (is.na(parsed) && nzchar(trimws(as.character(value)))) stop("Enter a valid numeric value.", call. = FALSE)
    x[row] <- parsed
  } else x[row] <- as.character(value)
  metadata[[column]] <- x; metadata
}

assign_metadata_value <- function(metadata, rows, column, value) {
  assert_metadata_column(metadata, column); rows <- unique(as.integer(rows))
  if (!length(rows) || anyNA(rows) || any(rows < 1L | rows > nrow(metadata))) stop("Select one or more valid samples in the table.", call. = FALSE)
  for (row in rows) metadata <- edit_metadata_cell(metadata, row, column, value)
  metadata
}

set_metadata_reference <- function(metadata, column, reference) {
  assert_metadata_column(metadata, column)
  if (is.numeric(metadata[[column]])) stop("Convert this numeric column to a factor before setting a reference.", call. = FALSE)
  values <- as.character(metadata[[column]])
  observed <- unique(values[!is.na(values) & nzchar(trimws(values))])
  if (!reference %in% observed) stop("The reference must be an observed value.", call. = FALSE)
  metadata[[column]] <- stats::relevel(factor(values, levels = observed), ref = reference); metadata
}

metadata_type_label <- function(x) if (is.factor(x)) "factor" else if (is.numeric(x)) "numeric" else "character"
