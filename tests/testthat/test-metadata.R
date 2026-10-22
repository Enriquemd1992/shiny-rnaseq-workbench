metadata_fixture <- function() data.frame(condition=c("control","control","treated","treated"), batch=c("A","B","A","B"), score=c("1","2","3","4"), row.names=paste0("sample_",1:4), stringsAsFactors=FALSE)

test_that("metadata columns convert without losing values", {
  x <- convert_metadata_column(metadata_fixture(), "score", "numeric")
  expect_type(x$score, "double"); expect_equal(x$score, 1:4)
  bad <- metadata_fixture(); bad$score[2] <- "bad"
  expect_error(convert_metadata_column(bad, "score", "numeric"), "would lose")
})

test_that("factor references are explicit", {
  x <- set_metadata_reference(metadata_fixture(), "condition", "treated")
  expect_true(is.factor(x$condition)); expect_equal(levels(x$condition)[1], "treated")
  expect_error(set_metadata_reference(transform(metadata_fixture(), score=as.numeric(score)), "score", "1"), "numeric column")
})

test_that("factor editing safely adds a level", {
  x <- convert_metadata_column(metadata_fixture(), "condition", "factor")
  x <- edit_metadata_cell(x, 1, "condition", "vehicle")
  expect_equal(as.character(x$condition[1]), "vehicle"); expect_true("vehicle" %in% levels(x$condition))
})

test_that("structural edits and bulk assignment are guarded", {
  x <- rename_metadata_column(metadata_fixture(), "batch", "sequencing_batch")
  expect_error(rename_metadata_column(x, "score", "condition"), "already exists")
  x <- add_metadata_column(x, "sex", "unknown"); x <- assign_metadata_value(x, c(1,3), "sex", "female")
  expect_equal(x$sex, c("female","unknown","female","unknown"))
  x <- delete_metadata_column(x, "score"); expect_false("score" %in% names(x))
})

test_that("metadata validation reports empty cells", {
  x <- metadata_fixture(); x$batch[2] <- ""
  expect_true("missing_metadata" %in% validate_metadata(x)$code)
})

test_that("draft edits do not propagate until Apply", {
  imported_metadata <- metadata_fixture()
  suppressWarnings(shiny::testServer(mod_metadata_server, args = list(imported = list(accepted = shiny::reactiveVal(list(metadata = imported_metadata)))), {
    session$flushReact()
    applied <- session$getReturned()
    expect_type(applied()$score, "character")

    session$setInputs(column = "score", type = "numeric", convert = 1)
    session$flushReact()
    expect_type(applied()$score, "character")

    session$setInputs(undo = 1)
    session$flushReact()
    expect_type(draft()$score, "character")

    session$setInputs(convert = 2)
    session$flushReact()

    session$setInputs(apply = 1)
    session$flushReact()
    expect_type(applied()$score, "double")
  }))
})
