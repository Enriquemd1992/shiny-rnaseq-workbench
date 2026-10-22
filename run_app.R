required <- c("shiny", "bslib", "DT", "readxl", "ggplot2", "scales")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse=", "), ". Run: Rscript install.R", call.=FALSE)
message("Starting Lab RNA-seq in your default browser...")
shiny::runApp(".", launch.browser = TRUE)
