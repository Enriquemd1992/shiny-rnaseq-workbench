validate_deseq_request <- function(counts, metadata, variables, contrast_variable, numerator, denominator) {
  count_issues <- validate_counts(counts)
  if (has_errors(count_issues)) stop(count_issues$message[count_issues$severity == "error"][1], call. = FALSE)
  if (length(variables) && any(!stats::complete.cases(metadata[, variables, drop = FALSE]))) stop("Design variables contain missing values. Complete the metadata before running DESeq2.", call. = FALSE)
  if (length(variables) && any(vapply(metadata[, variables, drop = FALSE], function(x) any(trimws(as.character(x)) == ""), logical(1)))) stop("Design variables contain empty values. Complete the metadata before running DESeq2.", call. = FALSE)
  issues <- validate_design(metadata, variables)
  if (has_errors(issues)) stop(issues$message[1], call. = FALSE)
  if (!contrast_variable %in% variables) stop("The contrast variable must be included in the design.", call. = FALSE)
  if (is.numeric(metadata[[contrast_variable]])) stop("The contrast variable must be categorical.", call. = FALSE)
  if (identical(numerator, denominator)) stop("Numerator and denominator must be different.", call. = FALSE)
  observed <- as.character(metadata[[contrast_variable]])
  if (!all(c(numerator, denominator) %in% observed)) stop("Both contrast levels must be present in metadata.", call. = FALSE)
  group_sizes <- table(observed[observed %in% c(numerator, denominator)])
  if (any(group_sizes < 2L)) stop("Each contrast group must contain at least two samples.", call. = FALSE)
  if (!identical(colnames(counts), rownames(metadata))) stop("Counts and metadata must be aligned before running DESeq2.", call. = FALSE)
  invisible(TRUE)
}

gene_filter_mask <- function(counts, annotations = NULL, gene_types = character(), min_count = 10, min_samples = 2) {
  if (!is.numeric(min_count) || min_count < 0 || !is.numeric(min_samples) || min_samples < 1 || min_samples > ncol(counts)) stop("Invalid low-expression filter settings.", call. = FALSE)
  keep <- rowSums(counts >= min_count) >= min_samples
  if (length(gene_types)) {
    if (is.null(annotations) || !"gene_type" %in% names(annotations)) stop("A gene_type annotation is required for this filter.", call. = FALSE)
    keep <- keep & as.character(annotations[rownames(counts),"gene_type"]) %in% gene_types
  }
  keep
}

run_deseq_analysis <- function(counts, metadata, variables, contrast_variable, numerator, denominator,
                               alpha = 0.05, shrinkage = "none", annotations = NULL,
                               gene_types = character(), min_count = 10, min_samples = 2, sample_selection = NULL) {
  validate_deseq_request(counts, metadata, variables, contrast_variable, numerator, denominator)
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) stop("Alpha must be between 0 and 1.", call. = FALSE)
  if (is.null(annotations)) annotations <- data.frame(row.names=rownames(counts))
  annotations <- annotations[rownames(counts),,drop=FALSE]
  keep <- gene_filter_mask(counts,annotations,gene_types,min_count,min_samples)
  if (sum(keep) < 2L) stop("Gene filters retain fewer than two genes. Relax the filters.", call. = FALSE)
  genes_input <- nrow(counts); counts <- counts[keep,,drop=FALSE]; annotations <- annotations[keep,,drop=FALSE]
  metadata <- as.data.frame(metadata)
  for (variable in variables) if (is.character(metadata[[variable]])) metadata[[variable]] <- factor(metadata[[variable]])
  metadata[[contrast_variable]] <- stats::relevel(factor(metadata[[contrast_variable]]), ref = denominator)
  design_formula <- stats::reformulate(variables)
  count_matrix <- round(as.matrix(counts))
  storage.mode(count_matrix) <- "integer"
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = count_matrix, colData = metadata, design = design_formula)
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  contrast <- c(contrast_variable, numerator, denominator)
  raw <- DESeq2::results(dds, contrast = contrast, alpha = alpha)
  annotation_table <- data.frame(gene_id=rownames(annotations),annotations,row.names=NULL,check.names=FALSE)
  table <- cbind(annotation_table,as.data.frame(raw)[annotation_table$gene_id,,drop=FALSE],row.names=NULL)
  table$shrunkenLog2FoldChange <- NA_real_
  shrinkage_used <- "none"
  if (identical(shrinkage, "normal")) {
    shrunk <- DESeq2::lfcShrink(dds, contrast = contrast, res = raw, type = "normal")
    table$shrunkenLog2FoldChange <- shrunk$log2FoldChange
    shrinkage_used <- "normal"
  }
  transformed <- DESeq2::varianceStabilizingTransformation(dds, blind = FALSE)
  list(
    dds = dds,
    results = table,
    normalized_counts = DESeq2::counts(dds, normalized = TRUE),
    transformed = SummarizedExperiment::assay(transformed),
    metadata = metadata,
    annotations = annotations,
    parameters = list(
      design = paste(deparse(design_formula), collapse = ""), variables = variables,
      contrast = contrast, contrast_variable = contrast_variable, numerator = numerator, denominator = denominator,
      alpha = alpha, shrinkage = shrinkage_used,
      gene_type_filter = if(length(gene_types)) paste(gene_types,collapse=", ") else "All gene types",
      min_count = min_count, min_samples = min_samples, genes_input = genes_input, genes_retained = sum(keep),
      sample_selection = if(is.null(sample_selection)) "All samples" else if(identical(sample_selection$mode,"all")) "All samples" else if(identical(sample_selection$mode,"metadata")) paste0(sample_selection$variable," in {",paste(sample_selection$values,collapse=", "),"}") else "Manual selection",
      included_samples = paste(colnames(counts),collapse=", "),
      excluded_samples = if(is.null(sample_selection)||!length(sample_selection$excluded)) "None" else paste(sample_selection$excluded,collapse=", "),
      executed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
    )
  )
}

classify_deseq_results <- function(table, alpha = 0.05, lfc_threshold = 1, use_shrunken = TRUE) {
  out <- table
  lfc <- if (use_shrunken && "shrunkenLog2FoldChange" %in% names(out) && any(!is.na(out$shrunkenLog2FoldChange))) out$shrunkenLog2FoldChange else out$log2FoldChange
  out$displayLog2FoldChange <- lfc
  out$direction <- "Not significant"
  valid <- !is.na(out$padj) & !is.na(lfc)
  out$direction[valid & out$padj <= alpha & lfc >= lfc_threshold] <- "Up"
  out$direction[valid & out$padj <= alpha & lfc <= -lfc_threshold] <- "Down"
  out
}

rank_deseq_results <- function(table, method = "padj") {
  score <- switch(method,
    padj = ifelse(is.na(table$padj), Inf, table$padj),
    lfc = -abs(table$displayLog2FoldChange),
    combined = -(abs(table$displayLog2FoldChange) * -log10(pmax(table$padj, .Machine$double.xmin))),
    stop("Unknown ranking method.", call. = FALSE)
  )
  table[order(score, na.last = TRUE), , drop = FALSE]
}

deseq_summary <- function(table) c(up = sum(table$direction == "Up"), down = sum(table$direction == "Down"), not_significant = sum(table$direction == "Not significant"))

volcano_plot <- function(table, alpha, lfc_threshold) {
  table$minus_log10_padj <- -log10(pmax(table$padj, .Machine$double.xmin))
  label_source <- if("gene_symbol" %in% names(table)) ifelse(is.na(table$gene_symbol)|!nzchar(table$gene_symbol),table$gene_id,table$gene_symbol) else table$gene_id
  table$.label <- ""; candidates <- which(table$direction!="Not significant" & !is.na(table$padj)); candidates <- candidates[order(table$padj[candidates])]
  table$.label[utils::head(candidates,10)] <- label_source[utils::head(candidates,10)]
  ggplot2::ggplot(table, ggplot2::aes(displayLog2FoldChange, minus_log10_padj, color = direction)) +
    ggplot2::geom_point(alpha = .55, size = 1.4, na.rm = TRUE) +
    ggplot2::geom_text(ggplot2::aes(label=.label),check_overlap=TRUE,vjust=-.6,size=3,show.legend=FALSE,na.rm=TRUE) +
    ggplot2::geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), linetype = 2, color = "grey45") +
    ggplot2::geom_hline(yintercept = -log10(alpha), linetype = 2, color = "grey45") +
    ggplot2::scale_color_manual(values = c(Down = "#2c7fb8", `Not significant` = "grey70", Up = "#d7301f")) +
    ggplot2::labs(x = "log2 fold change", y = "-log10 adjusted p-value", color = NULL) + ggplot2::theme_minimal()
}

ma_plot <- function(table, alpha, lfc_threshold) {
  ggplot2::ggplot(table, ggplot2::aes(log10(baseMean + 1), displayLog2FoldChange, color = direction)) +
    ggplot2::geom_point(alpha = .55, size = 1.4, na.rm = TRUE) +
    ggplot2::geom_hline(yintercept = c(-lfc_threshold, 0, lfc_threshold), linetype = c(2, 1, 2), color = "grey45") +
    ggplot2::scale_color_manual(values = c(Down = "#2c7fb8", `Not significant` = "grey70", Up = "#d7301f")) +
    ggplot2::labs(x = "log10 mean normalized count + 1", y = "log2 fold change", color = NULL) + ggplot2::theme_minimal()
}
