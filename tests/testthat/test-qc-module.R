test_that("QC PCA plots render for full data and a subset", {
  set.seed(7)
  counts <- matrix(stats::rpois(600,lambda=30),nrow=100,dimnames=list(paste0("g",1:100),paste0("s",1:6)))
  metadata <- data.frame(treatment=c("NC","NC","NC","NC","A","A"),source_group=c(1,1,2,2,1,2),row.names=colnames(counts))
  suppressWarnings(shiny::testServer(mod_qc_server,args=list(analysis_data=shiny::reactiveVal(list(counts=counts,metadata=metadata,annotations=data.frame(row.names=rownames(counts))))),{
    session$setInputs(color="treatment",shape="source_group",label="",subset_var="treatment",subset_level="NC")
    session$flushReact()
    expect_false(is.null(output$pca))
    expect_false(is.null(output$subset_pca))
    expect_false(is.null(output$download_pca))
  }))
})
