test_that("QC transformation filters genes and selects top variable genes", {
  set.seed(2); counts<-matrix(stats::rpois(1200,30),nrow=200,dimnames=list(paste0("g",1:200),paste0("s",1:6)));counts[1:10,]<-0
  q<-prepare_qc_transformation(counts,"vst",10,2,50)
  expect_equal(q$genes_retained,190);expect_equal(q$genes_used,50);expect_equal(nrow(q$pca_matrix),50)
  md<-data.frame(group=rep(c("A","B"),each=3),row.names=colnames(counts));scores<-compute_matrix_pca(q$pca_matrix,md)
  expect_equal(nrow(scores),6);expect_true(length(attr(scores,"variance"))>=2)
})

test_that("log-CPM and all-variable-gene modes are supported", {
  counts<-matrix(c(20,30,40,50,60,70,5,6,7,8,9,10),nrow=2,byrow=TRUE,dimnames=list(c("g1","g2"),paste0("s",1:6)))
  q<-prepare_qc_transformation(counts,"logcpm",0,1,Inf)
  expect_equal(q$genes_used,2);expect_equal(dim(q$matrix),c(2,6))
})
