test_that("annotation columns are excluded from counts and retained", {
  x <- data.frame(ensembl_id=c("e1","e2"),gene_symbol=c("A","A"),gene_type=c("protein_coding","lncRNA"),s1=c(1,2),s2=c(3,4),check.names=FALSE)
  prepared <- prepare_count_data(x,"ensembl_id",c("gene_symbol","gene_type"))
  expect_equal(suggest_annotation_columns(x,"ensembl_id"),c("gene_symbol","gene_type"))
  expect_equal(colnames(prepared$counts),c("s1","s2"))
  expect_equal(prepared$annotations$gene_symbol,c("A","A"))
  expect_true("duplicate_annotation" %in% validate_annotations(prepared$annotations)$code)
  expect_error(prepare_count_data(x,"gene_symbol",c("ensembl_id","gene_type")),"not unique")
})

test_that("gene filters combine expression and optional type", {
  counts <- matrix(c(20,20,20,0,5,5),nrow=3,byrow=TRUE,dimnames=list(c("g1","g2","g3"),c("s1","s2")))
  annotations <- data.frame(gene_type=c("protein_coding","protein_coding","lncRNA"),row.names=rownames(counts))
  expect_equal(unname(gene_filter_mask(counts,annotations,character(),10,1)),c(TRUE,TRUE,FALSE))
  expect_equal(unname(gene_filter_mask(counts,annotations,"protein_coding",10,2)),c(TRUE,FALSE,FALSE))
})
