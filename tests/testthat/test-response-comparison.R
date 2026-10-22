test_that("joint response comparison returns simple effects and interaction", {
  set.seed(91)
  md <- expand.grid(cell_model = c("A", "B"), treatment = c("UT", "IFNy"), replicate = 1:2,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  rownames(md) <- paste0("s", seq_len(nrow(md)))
  counts <- matrix(rnbinom(200 * nrow(md), mu = 60, size = 15), 200,
    dimnames = list(paste0("g", 1:200), rownames(md)))
  counts[1:15, md$treatment == "IFNy"] <- counts[1:15, md$treatment == "IFNy"] + 250
  counts[16:25, md$treatment == "IFNy" & md$cell_model == "B"] <- counts[16:25, md$treatment == "IFNy" & md$cell_model == "B"] + 400
  ann <- data.frame(gene_symbol = paste0("G", 1:200), row.names = rownames(counts))
  result <- run_response_comparison(counts, md, ann, "cell_model", "A", "B", "treatment", "IFNy", "UT",
    min_count = 0, min_samples = 1)
  expect_true(all(c("lfc_group_a", "lfc_group_b", "interaction_lfc", "response_category") %in% names(result$table)))
  expect_equal(sum(response_category_summary(result)$genes), nrow(result$table))
  expect_equal(as.integer(result$parameters$cell_sizes), rep(2L, 4))
})

test_that("response comparison requires replicated complete cells", {
  md <- data.frame(group = c("A", "A", "B"), treatment = c("UT", "X", "UT"), row.names = paste0("s", 1:3))
  counts <- matrix(10L, 20, 3, dimnames = list(paste0("g", 1:20), rownames(md)))
  expect_error(run_response_comparison(counts, md, group_variable = "group", group_a = "A", group_b = "B",
    treatment_variable = "treatment", treated = "X", control = "UT", min_count = 0, min_samples = 1), "at least two")
})
