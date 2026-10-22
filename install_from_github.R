# One-command installer for shiny-rnaseq-workbench.
# Run with:
# source("https://raw.githubusercontent.com/Enriquemd1992/shiny-rnaseq-workbench/main/install_from_github.R")

repo_archive <- paste0(
  "https://github.com/Enriquemd1992/",
  "shiny-rnaseq-workbench/archive/refs/heads/main.zip"
)
configured_dir <- Sys.getenv("SHINY_RNASEQ_INSTALL_DIR", unset = "")
install_dir <- if (nzchar(configured_dir)) {
  normalizePath(configured_dir, mustWork = FALSE)
} else {
  path.expand("~/shiny-rnaseq-workbench")
}
required_r <- "4.3"

current_r <- paste(R.version$major, strsplit(R.version$minor, ".", fixed = TRUE)[[1]][1], sep = ".")
if (!identical(current_r, required_r)) {
  stop(
    "This release requires R 4.3.x to reproduce its Bioconductor environment. ",
    "You are running R ", getRversion(), ".\n",
    "Install R 4.3.2, select it in RStudio, restart RStudio, and run this command again.",
    call. = FALSE
  )
}

is_existing_install <-
  dir.exists(install_dir) &&
  file.exists(file.path(install_dir, "renv.lock")) &&
  file.exists(file.path(install_dir, "run_app.R"))

if (dir.exists(install_dir) && !is_existing_install) {
  stop(
    "The destination already exists but is not a valid app installation:\n  ",
    install_dir, "\nRename or remove that folder, then try again.",
    call. = FALSE
  )
}

if (!is_existing_install) {
  message("Downloading shiny-rnaseq-workbench from GitHub...")
  archive <- tempfile(fileext = ".zip")
  unpack_dir <- tempfile("shiny-rnaseq-workbench-")
  dir.create(unpack_dir)

  download.file(repo_archive, archive, mode = "wb", quiet = TRUE)
  utils::unzip(archive, exdir = unpack_dir)

  extracted <- list.dirs(unpack_dir, recursive = FALSE, full.names = TRUE)
  if (length(extracted) != 1L) {
    stop("The downloaded repository could not be identified.", call. = FALSE)
  }

  dir.create(install_dir, recursive = TRUE)
  contents <- list.files(
    extracted,
    all.files = TRUE,
    full.names = TRUE,
    no.. = TRUE
  )
  copied <- file.copy(contents, install_dir, recursive = TRUE, copy.mode = TRUE)
  if (!all(copied)) {
    stop("Some application files could not be copied to: ", install_dir, call. = FALSE)
  }
} else {
  message("An existing installation was found; checking its R environment...")
}

old_dir <- setwd(install_dir)
on.exit(setwd(old_dir), add = TRUE)

options(repos = c(CRAN = "https://cloud.r-project.org"))
source(file.path("renv", "activate.R"), local = FALSE)
renv::restore(prompt = FALSE)

message("\nInstallation complete and environment restored.")
message("App folder: ", normalizePath(install_dir, winslash = "/"))
message("\nTo start the app now or later, run:")
message('source("~/shiny-rnaseq-workbench/run_app.R")')
