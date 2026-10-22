test_that("TPM export uses lengths and preserves selected samples", {
  counts <- matrix(c(100,200,300,400,50,50), nrow = 2, byrow = TRUE,
    dimnames = list(c("g1","g2"), c("s1","s2","s3")))
  ann <- data.frame(gene_symbol=c("A","B"),gene_type="protein_coding",row.names=rownames(counts))
  result <- compute_tpm_export(counts,ann,c(g1=1000,g2=2000),min_count=0,min_samples=1)
  expect_equal(unname(colSums(result$tpm)),rep(1e6,3),tolerance=1e-6)
  expect_identical(colnames(result$tpm),colnames(counts))
  expect_true(all(c("gene_id","gene_symbol","length_bp") %in% names(tpm_export_table(result))))
})

test_that("TPM export rejects incomplete or invalid lengths", {
  counts <- matrix(10,2,2,dimnames=list(c("g1","g2"),c("s1","s2")))
  expect_error(compute_tpm_export(counts,gene_lengths=c(g1=1000),min_count=0,min_samples=1),"missing")
  expect_error(prepare_gene_lengths(data.frame(id="g1",length=0),"id","length"),"positive")
})
