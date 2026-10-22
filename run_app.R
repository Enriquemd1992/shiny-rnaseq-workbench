source_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)

if (!is.null(source_file)) {
  app_dir <- dirname(normalizePath(source_file, mustWork = TRUE))
} else if (length(file_arg)) {
  app_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE))
} else {
  app_dir <- normalizePath(getwd(), mustWork = TRUE)
}

activate <- file.path(app_dir, "renv", "activate.R")
if (file.exists(activate)) source(activate, local = FALSE)

required <- c("shiny", "bslib", "DT", "readxl", "ggplot2", "scales")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) {
  stop(
    "Missing packages: ", paste(missing, collapse=", "),
    ". Run the GitHub installer again or execute install.R from the app folder.",
    call. = FALSE
  )
}
message("Starting Lab RNA-seq in your default browser...")
shiny::runApp(app_dir, launch.browser = TRUE)
