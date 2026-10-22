if (!requireNamespace("renv", quietly=TRUE)) install.packages("renv", repos="https://cloud.r-project.org")
if (file.exists("renv.lock")) {
  renv::restore(prompt=FALSE)
} else {
  message("renv.lock not found; installing declared dependencies.")
  renv::install()
  renv::snapshot(prompt=FALSE)
}
message("Installation complete. Start with: Rscript run_app.R")
