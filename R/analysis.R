library_sizes <- function(counts) colSums(counts)
detected_genes <- function(counts) colSums(counts > 0)

log_cpm <- function(counts, prior_count = 0.5) {
  sizes <- colSums(counts)
  log2(t(t(counts + prior_count) / (sizes + prior_count * nrow(counts))) * 1e6)
}

compute_pca <- function(counts, metadata) {
  transformed <- log_cpm(counts)
  keep <- apply(transformed, 1, stats::var) > 0
  if (sum(keep) < 2L) stop("At least two variable genes are needed for PCA.", call. = FALSE)
  fit <- stats::prcomp(t(transformed[keep, , drop = FALSE]), scale. = FALSE)
  scores <- cbind(sample_id = rownames(fit$x), as.data.frame(fit$x[, 1:2, drop = FALSE]), metadata[rownames(fit$x), , drop = FALSE])
  attr(scores, "variance") <- summary(fit)$importance[2, 1:2] * 100
  scores
}

prepare_qc_transformation <- function(counts, method = "vst", min_count = 10, min_samples = 2, top_n = 500) {
  if (!method %in% c("vst","rlog","logcpm")) stop("Unknown QC transformation.",call.=FALSE)
  if (min_count < 0 || min_samples < 1 || min_samples > ncol(counts)) stop("Invalid QC expression filter.",call.=FALSE)
  keep <- rowSums(counts >= min_count) >= min_samples
  if (sum(keep) < 2L) stop("QC filters retain fewer than two genes.",call.=FALSE)
  filtered <- counts[keep,,drop=FALSE]
  if (method == "logcpm") {
    transformed <- log_cpm(filtered)
  } else {
    integer_counts <- round(filtered); storage.mode(integer_counts) <- "integer"
    dds <- DESeq2::DESeqDataSetFromMatrix(integer_counts,data.frame(row.names=colnames(integer_counts)),~1)
    object <- if(method=="vst") DESeq2::varianceStabilizingTransformation(dds,blind=TRUE) else DESeq2::rlog(dds,blind=TRUE)
    transformed <- SummarizedExperiment::assay(object)
  }
  variances <- apply(transformed,1,stats::var)
  variable <- names(sort(variances,decreasing=TRUE)); variable <- variable[is.finite(variances[variable]) & variances[variable]>0]
  if (!length(variable)) stop("No variable genes remain for PCA.",call.=FALSE)
  n_use <- if(is.infinite(top_n)) length(variable) else min(as.integer(top_n),length(variable))
  list(matrix=transformed,pca_matrix=transformed[variable[seq_len(n_use)],,drop=FALSE],genes_retained=nrow(transformed),genes_used=n_use,method=method,min_count=min_count,min_samples=min_samples)
}

compute_matrix_pca <- function(matrix, metadata) {
  if (!identical(colnames(matrix),rownames(metadata))) stop("Transformed matrix and metadata are not aligned.",call.=FALSE)
  fit <- stats::prcomp(t(matrix),center=TRUE,scale.=FALSE)
  scores <- cbind(sample_id=rownames(fit$x),as.data.frame(fit$x),metadata[rownames(fit$x),,drop=FALSE])
  attr(scores,"variance") <- summary(fit)$importance[2,]*100
  scores
}
