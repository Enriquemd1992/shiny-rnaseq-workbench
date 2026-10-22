run_response_comparison <- function(counts, metadata, annotations = NULL,
                                    group_variable, group_a, group_b,
                                    treatment_variable, treated, control,
                                    alpha = 0.05, lfc_threshold = 1,
                                    min_count = 10, min_samples = 2) {
  required <- c(group_variable, treatment_variable)
  if (any(!required %in% names(metadata))) stop("Select valid group and treatment variables.", call. = FALSE)
  if (identical(group_variable, treatment_variable)) stop("Group and treatment variables must be different.", call. = FALSE)
  if (identical(group_a, group_b) || identical(treated, control)) stop("Select two different levels for each variable.", call. = FALSE)
  if (!all(c(group_a, group_b) %in% as.character(metadata[[group_variable]]))) stop("Both group levels must be present.", call. = FALSE)
  if (!all(c(treated, control) %in% as.character(metadata[[treatment_variable]]))) stop("Both treatment levels must be present.", call. = FALSE)

  keep_samples <- as.character(metadata[[group_variable]]) %in% c(group_a, group_b) &
    as.character(metadata[[treatment_variable]]) %in% c(treated, control)
  md <- metadata[keep_samples, , drop = FALSE]
  count_matrix <- counts[, rownames(md), drop = FALSE]
  md[[group_variable]] <- stats::relevel(factor(md[[group_variable]]), group_a)
  md[[treatment_variable]] <- stats::relevel(factor(md[[treatment_variable]]), control)
  cell_sizes <- table(md[[group_variable]], md[[treatment_variable]])
  if (any(cell_sizes < 2L)) stop("Each group-by-treatment cell must contain at least two samples.", call. = FALSE)
  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1) stop("Alpha must be between 0 and 1.", call. = FALSE)
  if (!is.numeric(lfc_threshold) || lfc_threshold < 0) stop("The LFC threshold must be non-negative.", call. = FALSE)

  gene_keep <- rowSums(count_matrix >= min_count) >= min_samples
  if (sum(gene_keep) < 2L) stop("Gene filters retain fewer than two genes.", call. = FALSE)
  count_matrix <- round(count_matrix[gene_keep, , drop = FALSE]); storage.mode(count_matrix) <- "integer"
  if (is.null(annotations)) annotations <- data.frame(row.names = rownames(counts))
  ann <- annotations[rownames(count_matrix), , drop = FALSE]

  design <- stats::as.formula(sprintf("~ `%s` * `%s`", group_variable, treatment_variable))
  dds <- DESeq2::DESeqDataSetFromMatrix(count_matrix, md, design)
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  coefficient_names <- DESeq2::resultsNames(dds)
  treatment_coefficient <- grep(paste0("^", treatment_variable, ".*_vs_"), coefficient_names, value = TRUE)
  interaction_coefficient <- grep(paste0(group_variable, ".*\\.", treatment_variable), coefficient_names, value = TRUE)
  if (length(treatment_coefficient) != 1L || length(interaction_coefficient) != 1L)
    stop("Could not identify the treatment and interaction coefficients.", call. = FALSE)

  result_a <- DESeq2::results(dds, name = treatment_coefficient, alpha = alpha)
  result_b <- DESeq2::results(dds, contrast = list(c(treatment_coefficient, interaction_coefficient)), alpha = alpha)
  interaction <- DESeq2::results(dds, name = interaction_coefficient, alpha = alpha)
  genes <- rownames(dds)
  table <- data.frame(gene_id = genes, ann, baseMean = result_a$baseMean,
    lfc_group_a = result_a$log2FoldChange, pvalue_group_a = result_a$pvalue, padj_group_a = result_a$padj,
    lfc_group_b = result_b$log2FoldChange, pvalue_group_b = result_b$pvalue, padj_group_b = result_b$padj,
    interaction_lfc = interaction$log2FoldChange, interaction_pvalue = interaction$pvalue,
    interaction_padj = interaction$padj, check.names = FALSE, row.names = NULL)

  significant_a <- !is.na(table$padj_group_a) & table$padj_group_a <= alpha & abs(table$lfc_group_a) >= lfc_threshold
  significant_b <- !is.na(table$padj_group_b) & table$padj_group_b <= alpha & abs(table$lfc_group_b) >= lfc_threshold
  direction_a <- sign(table$lfc_group_a); direction_b <- sign(table$lfc_group_b)
  category <- rep("Not detected", nrow(table))
  category[significant_a & significant_b & direction_a > 0 & direction_b > 0] <- "Shared up"
  category[significant_a & significant_b & direction_a < 0 & direction_b < 0] <- "Shared down"
  category[significant_a & significant_b & direction_a > 0 & direction_b < 0] <- paste0(group_a, " up / ", group_b, " down")
  category[significant_a & significant_b & direction_a < 0 & direction_b > 0] <- paste0(group_a, " down / ", group_b, " up")
  category[significant_a & !significant_b & direction_a > 0] <- paste0("Detected only in ", group_a, " — up")
  category[significant_a & !significant_b & direction_a < 0] <- paste0("Detected only in ", group_a, " — down")
  category[!significant_a & significant_b & direction_b > 0] <- paste0("Detected only in ", group_b, " — up")
  category[!significant_a & significant_b & direction_b < 0] <- paste0("Detected only in ", group_b, " — down")
  table$response_category <- category
  table$response_difference_significant <- !is.na(table$interaction_padj) & table$interaction_padj <= alpha
  table <- table[order(table$interaction_padj, pmin(table$padj_group_a, table$padj_group_b, na.rm = TRUE), na.last = TRUE), , drop = FALSE]

  transformed <- SummarizedExperiment::assay(DESeq2::varianceStabilizingTransformation(dds, blind = FALSE))
  list(table = table, dds = dds, transformed = transformed, metadata = md,
    parameters = list(group_variable = group_variable, group_a = group_a, group_b = group_b,
      treatment_variable = treatment_variable, treated = treated, control = control,
      alpha = alpha, lfc_threshold = lfc_threshold, min_count = min_count,
      min_samples = min_samples, genes_retained = nrow(dds), cell_sizes = cell_sizes,
      treatment_coefficient = treatment_coefficient, interaction_coefficient = interaction_coefficient))
}

response_category_summary <- function(result) {
  as.data.frame(table(result$table$response_category), stringsAsFactors = FALSE) |>
    stats::setNames(c("category", "genes"))
}

response_scatter_plot <- function(result) {
  x <- result$table; p <- result$parameters
  shown <- x$response_category != "Not detected"
  x$response_direction <- ifelse(x$lfc_group_a > 0 & x$lfc_group_b > 0, "Up in both",
    ifelse(x$lfc_group_a < 0 & x$lfc_group_b < 0, "Down in both",
      ifelse(x$lfc_group_a > 0 & x$lfc_group_b < 0, paste0(p$group_a, " up / ", p$group_b, " down"),
        paste0(p$group_a, " down / ", p$group_b, " up"))))
  x$detection_status <- ifelse(grepl("^Shared| / ", x$response_category), "DEG in both",
    ifelse(grepl(paste0("only in ", p$group_a), x$response_category, fixed = TRUE), paste0("DEG only in ", p$group_a),
      ifelse(grepl(paste0("only in ", p$group_b), x$response_category, fixed = TRUE), paste0("DEG only in ", p$group_b), "Not detected")))
  labels <- rep("", nrow(x)); candidates <- which(x$response_difference_significant)
  candidates <- head(candidates[order(x$interaction_padj[candidates])], 12)
  symbol <- if ("gene_symbol" %in% names(x)) as.character(x$gene_symbol) else x$gene_id
  symbol[is.na(symbol) | !nzchar(symbol)] <- x$gene_id[is.na(symbol) | !nzchar(symbol)]
  labels[candidates] <- symbol[candidates]
  ggplot2::ggplot(x, ggplot2::aes(lfc_group_a, lfc_group_b)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey80") + ggplot2::geom_vline(xintercept = 0, colour = "grey80") +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey45") +
    ggplot2::geom_point(colour = "grey82", alpha = .22, size = 1) +
    ggplot2::geom_point(data = x[shown, , drop = FALSE],
      ggplot2::aes(colour = response_direction, shape = detection_status), alpha = .8, size = 2) +
    ggplot2::geom_point(data = x[x$response_difference_significant, , drop = FALSE],
      shape = 21, fill = NA, colour = "black", stroke = 1.2, size = 3.6) +
    ggplot2::geom_text(ggplot2::aes(label = labels), check_overlap = TRUE, size = 3, vjust = -.6) +
    ggplot2::scale_colour_manual(values = c("Up in both" = "#d7301f", "Down in both" = "#2c7fb8",
      stats::setNames("#7a0177", paste0(p$group_a, " up / ", p$group_b, " down")),
      stats::setNames("#1b9e77", paste0(p$group_a, " down / ", p$group_b, " up")))) +
    ggplot2::scale_shape_manual(values = stats::setNames(c(16, 17, 15),
      c("DEG in both", paste0("DEG only in ", p$group_a), paste0("DEG only in ", p$group_b)))) +
    ggplot2::coord_equal() + ggplot2::labs(
      x = paste0("log2FC ", p$group_a, ": ", p$treated, " vs ", p$control),
      y = paste0("log2FC ", p$group_b, ": ", p$treated, " vs ", p$control),
      colour = "Estimated direction", shape = "DEG threshold") +
    ggplot2::theme_minimal() + ggplot2::theme(legend.position = "right")
}

interaction_volcano_plot <- function(result) {
  x <- result$table; p <- result$parameters
  x$minus_log10_padj <- -log10(pmax(x$interaction_padj, .Machine$double.xmin))
  x$status <- ifelse(x$response_difference_significant, "Different response", "Not significant")
  ggplot2::ggplot(x, ggplot2::aes(interaction_lfc, minus_log10_padj, colour = status)) +
    ggplot2::geom_point(alpha = .55, size = 1.3, na.rm = TRUE) +
    ggplot2::geom_hline(yintercept = -log10(p$alpha), linetype = 2, colour = "grey45") +
    ggplot2::scale_colour_manual(values = c(`Different response` = "#d7301f", `Not significant` = "grey72")) +
    ggplot2::labs(x = paste0("Interaction log2FC (", p$group_b, " response minus ", p$group_a, " response)"),
      y = "-log10 interaction adjusted p-value", colour = NULL) + ggplot2::theme_minimal()
}
