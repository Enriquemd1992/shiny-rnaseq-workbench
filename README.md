# Lab RNA-seq app

Local R Shiny application for safe, reproducible bulk RNA-seq workflows. The current version includes validated import, safe metadata editing, QC, explicit DESeq2 contrasts, interactive results, plots, and exports.

## Privacy

Files are read only by the local R process. The app contains no upload-to-cloud code. Never commit private matrices, metadata, FASTQ, BAM, credentials, or saved projects. The included `example_data/` files are synthetic.

## Install from Git

Prerequisites: R >= 4.3, Git, and build tools where a binary R package is unavailable (Rtools on Windows, Xcode Command Line Tools on macOS, or compilers and common development libraries on Linux).

```sh
git clone <URL>
cd <repository>
Rscript install.R
Rscript run_app.R
```

## Install from ZIP

Download the repository ZIP, extract it, open a terminal in the extracted folder, and run:

```sh
Rscript install.R
Rscript run_app.R
```

`renv` restores the package versions recorded in `renv.lock`. Bioconductor packages such as DESeq2 may require system libraries and a compatible R/Bioconductor release.

## Current workflow

1. Import counts and metadata; select identifier columns.
   Gene annotation columns such as `gene_symbol` and `gene_type` can be selected separately and are excluded from the numeric matrix.
2. Preview and validate the complete dataset, then optionally retain a cohort by metadata values or by manually selecting samples. Counts and metadata are subset together and kept in the same order.
3. Accept the selected data only after blocking errors are resolved. Changing the sample selection invalidates the previous acceptance and downstream analyses.
4. Edit metadata in a draft with undo, bulk assignment, safe type conversion, factor references, validation, and explicit Apply/Discard controls.
5. Explore library sizes, detected genes, and PCA. No automatic outlier removal occurs.
6. Configure an ordered design and an explicit numerator-versus-denominator contrast in **Analysis setup**.
7. Run DESeq2 and explore the complete results, selectable Top N, volcano, MA plot, heatmap, and per-gene normalized expression.
8. Export the complete table, displayed Top N, parameters, and plots. Display thresholds never overwrite the original DESeq2 statistics.

Real datasets, study-specific scripts, result workbooks, logs, and internal validation notes are intentionally excluded from the repository. Use the synthetic files in `example_data/` to test the complete import workflow.

Run tests with:

```sh
Rscript tests/testthat.R
```

## Scope and roadmap

The next milestones are correlation/distance QC heatmaps, a complete reproducibility bundle and equivalent R script, and local project persistence. Interactions, LRT, enrichment analysis, and automatic outlier removal remain out of scope for this MVP.

This repository is intended for private internal evaluation. Do not commit scientific datasets, generated results, credentials, or identifying metadata.
