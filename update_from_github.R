# One-command updater for an existing shiny-rnaseq-workbench installation.
# Close the running app, then run:
# source("https://raw.githubusercontent.com/Enriquemd1992/shiny-rnaseq-workbench/main/update_from_github.R")

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

required_existing <- c("app.R", "run_app.R", "renv.lock")
missing_existing <- required_existing[
  !file.exists(file.path(install_dir, required_existing))
]
if (!dir.exists(install_dir) || length(missing_existing)) {
  stop(
    "No valid installation was found at:\n  ", install_dir, "\n",
    "Run the installation command first.",
    call. = FALSE
  )
}

message("Downloading the latest application code from GitHub...")
archive <- tempfile(fileext = ".zip")
unpack_dir <- tempfile("shiny-rnaseq-workbench-update-")
backup_dir <- tempfile("shiny-rnaseq-workbench-backup-")
dir.create(unpack_dir)
dir.create(backup_dir)

download.file(repo_archive, archive, mode = "wb", quiet = TRUE)
utils::unzip(archive, exdir = unpack_dir)

extracted <- list.dirs(unpack_dir, recursive = FALSE, full.names = TRUE)
if (length(extracted) != 1L) {
  stop("The downloaded repository could not be identified.", call. = FALSE)
}

required_download <- c(
  "app.R", "run_app.R", "renv.lock", "install_from_github.R"
)
missing_download <- required_download[
  !file.exists(file.path(extracted, required_download))
]
if (length(missing_download)) {
  stop(
    "The downloaded update is incomplete: ",
    paste(missing_download, collapse = ", "),
    call. = FALSE
  )
}

source_files <- list.files(
  extracted,
  all.files = TRUE,
  recursive = TRUE,
  full.names = TRUE,
  include.dirs = FALSE,
  no.. = TRUE
)
relative_files <- substring(source_files, nchar(extracted) + 2L)
destination_files <- file.path(install_dir, relative_files)
existing <- file.exists(destination_files)

# Keep temporary copies of files that will be replaced. If copying fails, the
# previous application files are restored. Extra user-created files are never
# selected or removed by this updater.
if (any(existing)) {
  backup_files <- file.path(backup_dir, relative_files[existing])
  invisible(lapply(unique(dirname(backup_files)), dir.create, recursive = TRUE, showWarnings = FALSE))
  backed_up <- file.copy(
    destination_files[existing],
    backup_files,
    overwrite = TRUE,
    copy.mode = TRUE
  )
  if (!all(backed_up)) {
    stop("The existing installation could not be backed up safely.", call. = FALSE)
  }
}

invisible(lapply(unique(dirname(destination_files)), dir.create, recursive = TRUE, showWarnings = FALSE))
copied <- file.copy(
  source_files,
  destination_files,
  overwrite = TRUE,
  copy.mode = TRUE
)

if (!all(copied)) {
  if (any(existing)) {
    file.copy(
      file.path(backup_dir, relative_files[existing]),
      destination_files[existing],
      overwrite = TRUE,
      copy.mode = TRUE
    )
  }
  stop(
    "The update could not be copied completely. Existing files were restored. ",
    "Close the app and try again.",
    call. = FALSE
  )
}

message("Application files updated. Checking the reproducible R environment...")
source(file.path(install_dir, "install_from_github.R"), local = FALSE)

message("\nUpdate complete.")
message("Start the app with:")
message('source("~/shiny-rnaseq-workbench/run_app.R")')
