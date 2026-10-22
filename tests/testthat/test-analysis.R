test_that("QC metrics use raw counts", {
  counts <- matrix(c(0,2,3,4),nrow=2,dimnames=list(c("g1","g2"),c("s1","s2")))
  expect_equal(unname(library_sizes(counts)),c(2,7))
  expect_equal(unname(detected_genes(counts)),c(1,2))
})
