test_that("valid integer counts pass", {
  counts <- matrix(c(1,2,3,4), nrow=2, dimnames=list(c("g1","g2"),c("s1","s2")))
  expect_equal(nrow(validate_counts(counts)), 0)
})

test_that("invalid count values are rejected", {
  counts <- matrix(c(1,-1,2.5,NA), nrow=2, dimnames=list(c("g1","g2"),c("s1","s2")))
  expect_true(all(c("missing_counts","negative_counts","non_integer_counts") %in% validate_counts(counts)$code))
})

test_that("sample order mismatch is explicit and reversible", {
  counts <- matrix(1:4,nrow=2,dimnames=list(c("g1","g2"),c("s1","s2")))
  metadata <- data.frame(group=c("B","A"),row.names=c("s2","s1"))
  expect_equal(validate_sample_match(counts,metadata)$code,"sample_order")
  aligned <- align_metadata(counts,metadata)
  expect_identical(rownames(aligned),colnames(counts))
})

test_that("confounded design is detected", {
  metadata <- data.frame(batch=factor(c("A","A","B","B")),condition=factor(c("X","X","Y","Y")))
  expect_equal(validate_design(metadata,c("batch","condition"))$code,"rank_deficient")
})
