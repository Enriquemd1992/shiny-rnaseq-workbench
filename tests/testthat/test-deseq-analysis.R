make_deseq_fixture <- function() {
  set.seed(42); n <- 120; samples <- paste0("s",1:8)
  metadata <- data.frame(batch=factor(rep(c("A","B"),4)), condition=factor(rep(c("control","treated"),each=4)), row.names=samples)
  mu <- matrix(40,n,8); mu[1:15,5:8] <- 120
  counts <- matrix(stats::rnbinom(n*8,mu=mu,size=12),nrow=n,dimnames=list(paste0("gene",1:n),samples))
  list(counts=counts,metadata=metadata)
}

test_that("DESeq2 wrapper matches an explicit manual contrast", {
  x <- make_deseq_fixture()
  fit <- run_deseq_analysis(x$counts,x$metadata,c("batch","condition"),"condition","treated","control",alpha=.05,min_count=0,min_samples=1)
  md <- x$metadata; md$condition <- stats::relevel(md$condition,"control")
  dds <- DESeq2::DESeqDataSetFromMatrix(x$counts,md,~batch+condition); dds <- DESeq2::DESeq(dds,quiet=TRUE)
  manual <- DESeq2::results(dds,contrast=c("condition","treated","control"),alpha=.05)
  expect_equal(fit$results$log2FoldChange,manual$log2FoldChange,tolerance=1e-10)
  expect_equal(fit$results$padj,manual$padj,tolerance=1e-10)
})

test_that("invalid contrasts are blocked", {
  x <- make_deseq_fixture()
  expect_error(validate_deseq_request(x$counts,x$metadata,"batch","condition","treated","control"),"included")
  expect_error(validate_deseq_request(x$counts,x$metadata,"condition","condition","control","control"),"different")
  x$metadata$condition[1] <- NA
  expect_error(validate_deseq_request(x$counts,x$metadata,"condition","condition","treated","control"),"missing values")
})

test_that("classification and ranking do not alter raw statistics", {
  tab <- data.frame(gene_id=c("a","b","c"),baseMean=1:3,log2FoldChange=c(2,-2,.2),lfcSE=1,stat=1,pvalue=c(.001,.002,.5),padj=c(.01,.02,.8),shrunkenLog2FoldChange=NA_real_)
  classified <- classify_deseq_results(tab,.05,1)
  expect_equal(classified$direction,c("Up","Down","Not significant")); expect_equal(classified$log2FoldChange,tab$log2FoldChange)
  expect_equal(rank_deseq_results(classified,"padj")$gene_id,c("a","b","c"))
})

test_that("analysis setup runs and invalidates stale results", {
  x <- make_deseq_fixture(); data_source <- shiny::reactiveVal(x)
  suppressWarnings(shiny::testServer(mod_analysis_setup_server,args=list(analysis_data=data_source),{
    analysis_result <- session$getReturned()
    session$setInputs(variables=c("batch","condition"),contrast_var="condition",numerator="treated",denominator="control",alpha=.05,shrinkage="none",gene_types=character(),min_count=0,min_samples=1,run=1)
    session$flushReact()
    expect_false(is.null(analysis_result()))
    expect_equal(analysis_result()$parameters$contrast,c("condition","treated","control"))
    session$setInputs(alpha=.1)
    session$flushReact()
    expect_null(analysis_result())
  }))
})
