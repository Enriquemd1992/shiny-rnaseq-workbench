test_that("samples can be selected by metadata while preserving count order", {
  counts<-matrix(1:24,nrow=4,dimnames=list(paste0("g",1:4),paste0("s",1:6)))
  metadata<-data.frame(experiment=c("A","A","B","B","B","B"),replicate=1:6,row.names=colnames(counts))
  ids<-select_sample_ids(metadata,"metadata","experiment","B")
  selected<-subset_samples(counts,metadata,ids)
  expect_identical(colnames(selected$counts),paste0("s",3:6));expect_identical(rownames(selected$metadata),colnames(selected$counts))
})

test_that("manual selection is explicit and guarded", {
  metadata<-data.frame(group=c("A","A","B"),row.names=c("s1","s2","s3"))
  expect_equal(select_sample_ids(metadata,"manual",manual_ids=c("s3","s1")),c("s1","s3"))
  expect_error(select_sample_ids(metadata,"metadata","group","missing"),"at least two")
})
