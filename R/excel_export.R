safe_excel_table_name <- function(name, index) {
  clean <- gsub("[^A-Za-z0-9_]", "_", name)
  clean <- sub("^[^A-Za-z_]", "T_", clean)
  paste0(substr(clean, 1, 220), "_", index)
}

write_filterable_workbook <- function(sheets, file) {
  if (!length(sheets) || is.null(names(sheets)) || any(!nzchar(names(sheets)))) {
    stop("Workbook sheets must be a named list.", call. = FALSE)
  }

  workbook <- openxlsx::createWorkbook()
  header_style <- openxlsx::createStyle(
    fontColour = "#FFFFFF", fgFill = "#2C3E50", textDecoration = "bold",
    halign = "center", valign = "center"
  )
  numeric_style <- openxlsx::createStyle(numFmt = "0.0000")
  probability_style <- openxlsx::createStyle(numFmt = "0.00E+00")
  integer_style <- openxlsx::createStyle(numFmt = "0")
  wrapped_style <- openxlsx::createStyle(wrapText = TRUE, valign = "top")

  for (index in seq_along(sheets)) {
    sheet_name <- substr(names(sheets)[index], 1, 31)
    data <- as.data.frame(sheets[[index]], check.names = FALSE, stringsAsFactors = FALSE)
    openxlsx::addWorksheet(workbook, sheet_name)
    openxlsx::writeDataTable(
      workbook, sheet_name, data, startRow = 1, startCol = 1,
      tableName = safe_excel_table_name(sheet_name, index),
      tableStyle = "TableStyleMedium2", withFilter = TRUE,
      headerStyle = header_style
    )
    openxlsx::freezePane(workbook, sheet_name, firstRow = TRUE)

    if (ncol(data)) {
      widths <- rep(14, ncol(data))
      widths[grepl("gene_id|ensembl", names(data), ignore.case = TRUE)] <- 24
      widths[grepl("gene_symbol", names(data), ignore.case = TRUE)] <- 18
      widths[grepl("gene_type", names(data), ignore.case = TRUE)] <- 30
      widths[names(data) == "parameter"] <- 30
      widths[names(data) == "value"] <- 70
      openxlsx::setColWidths(workbook, sheet_name, cols = seq_len(ncol(data)), widths = widths)

      numeric_cols <- which(vapply(data, is.numeric, logical(1)))
      probability_cols <- which(names(data) %in% c("pvalue", "padj"))
      integer_cols <- which(names(data) %in% c("length_bp"))
      ordinary_numeric <- setdiff(numeric_cols, c(probability_cols, integer_cols))
      if (nrow(data) && length(ordinary_numeric)) openxlsx::addStyle(workbook, sheet_name, numeric_style, rows = 2:(nrow(data) + 1), cols = ordinary_numeric, gridExpand = TRUE)
      if (nrow(data) && length(probability_cols)) openxlsx::addStyle(workbook, sheet_name, probability_style, rows = 2:(nrow(data) + 1), cols = probability_cols, gridExpand = TRUE)
      if (nrow(data) && length(integer_cols)) openxlsx::addStyle(workbook, sheet_name, integer_style, rows = 2:(nrow(data) + 1), cols = integer_cols, gridExpand = TRUE)
      if (nrow(data) && "value" %in% names(data)) openxlsx::addStyle(workbook, sheet_name, wrapped_style, rows = 2:(nrow(data) + 1), cols = match("value", names(data)), gridExpand = TRUE)
    }
  }
  openxlsx::saveWorkbook(workbook, file, overwrite = TRUE)
  invisible(file)
}

bundled_gene_lengths_path <- function() {
  working_directory <- normalizePath(getwd(), mustWork = FALSE)
  search_roots <- unique(c(
    working_directory,
    dirname(working_directory),
    dirname(dirname(working_directory))
  ))
  candidates <- c(
    getOption("labrnaseq.gene_lengths", ""),
    Sys.getenv("LABRNASEQ_GENE_LENGTHS", unset = ""),
    file.path(search_roots, "data", "gencode_v22_exonic_gene_lengths.tsv")
  )
  candidates <- candidates[nzchar(candidates)]
  existing <- candidates[file.exists(candidates)]
  if (length(existing)) normalizePath(existing[1]) else ""
}
