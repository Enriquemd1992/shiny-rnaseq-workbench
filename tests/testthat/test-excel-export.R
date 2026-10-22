test_that("filterable workbook contains Excel tables and preserves logical columns", {
  output <- tempfile(fileext = ".xlsx")
  data <- data.frame(
    gene_id = c("ENSG1", "ENSG2"),
    log2FoldChange = c(1.25, -0.75),
    padj = c(0.01, 0.2),
    significant = c(TRUE, FALSE)
  )
  write_filterable_workbook(
    list(Full_results = data, Parameters = data.frame(parameter = "filter", value = "count >= 10")),
    output
  )
  expect_true(file.exists(output))
  expect_gt(file.info(output)$size, 0)
  workbook <- openxlsx::loadWorkbook(output)
  expect_equal(names(workbook), c("Full_results", "Parameters"))
  imported <- openxlsx::read.xlsx(output, sheet = "Full_results", colNames = TRUE)
  expect_equal(imported$significant, c(TRUE, FALSE))
  table_xml <- unzip(output, list = TRUE)$Name
  expect_true(any(grepl("xl/tables/table", table_xml)))
})

test_that("bundled GENCODE length path is discoverable", {
  path <- bundled_gene_lengths_path()
  expect_true(nzchar(path))
  expect_true(file.exists(path))
})
