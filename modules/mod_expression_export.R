mod_expression_export_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::fluidPage(
    shiny::h2("Expression export"),
    shiny::p("Generate filtered raw-count and TPM tables for the currently accepted samples. TPM is for expression reporting; DESeq2 must continue to use raw counts."),
    bslib::layout_sidebar(
      sidebar = bslib::sidebar(
        shiny::radioButtons(ns("length_source"), "Gene lengths", choices = c("Bundled GENCODE v22" = "bundled", "Upload another annotation" = "upload"), selected = "bundled"),
        shiny::conditionalPanel(
          sprintf("input['%s'] === 'upload'", ns("length_source")),
          shiny::fileInput(ns("length_file"), "Gene lengths (CSV/TSV/XLSX)", accept = c(".csv", ".tsv", ".txt", ".xlsx")),
          shiny::selectInput(ns("length_gene_id"), "Gene ID column", choices = NULL),
          shiny::selectInput(ns("length_value"), "Length column (bp)", choices = NULL)
        ),
        shiny::uiOutput(ns("length_status")),
        shiny::selectizeInput(ns("gene_types"), "Gene types (empty = all)", choices = NULL, multiple = TRUE),
        shiny::numericInput(ns("min_count"), "Minimum raw count", 10, min = 0, step = 1),
        shiny::numericInput(ns("min_samples"), "Minimum samples", 2, min = 1, step = 1),
        shiny::actionButton(ns("generate"), "Generate TPM", class = "btn-primary"),
        shiny::hr(),
        shiny::downloadButton(ns("download_tpm"), "Download TPM CSV"),
        shiny::downloadButton(ns("download_counts"), "Download filtered counts CSV"),
        shiny::downloadButton(ns("download_workbook"), "Download workbook")
      ),
      bslib::navset_card_tab(
        bslib::nav_panel("Summary", shiny::uiOutput(ns("summary"))),
        bslib::nav_panel("TPM preview", DT::DTOutput(ns("tpm_preview"))),
        bslib::nav_panel("Filtered counts", DT::DTOutput(ns("counts_preview")))
      )
    )
  )
}

mod_expression_export_server <- function(id, analysis_data) shiny::moduleServer(id, function(input, output, session) {
  raw_lengths <- shiny::reactiveVal(NULL); export_result <- shiny::reactiveVal(NULL)
  observeEvent(input$length_file, {
    x <- tryCatch(read_tabular_file(input$length_file$datapath, input$length_file$name),
      error = function(e) { shiny::showNotification(e$message, type = "error"); NULL })
    raw_lengths(x)
    if (!is.null(x)) {
      id_guess <- intersect(c("ensembl_id", "gene_id", "Geneid"), names(x)); id_guess <- if(length(id_guess)) id_guess[1] else names(x)[1]
      length_guess <- grep("length", names(x), ignore.case = TRUE, value = TRUE); length_guess <- if(length(length_guess)) length_guess[1] else names(x)[min(2, ncol(x))]
      shiny::updateSelectInput(session, "length_gene_id", choices = names(x), selected = id_guess)
      shiny::updateSelectInput(session, "length_value", choices = names(x), selected = length_guess)
    }
  })
  observeEvent(analysis_data(), {
    req(analysis_data()); types <- if("gene_type" %in% names(analysis_data()$annotations)) sort(unique(as.character(analysis_data()$annotations$gene_type))) else character()
    shiny::updateSelectizeInput(session, "gene_types", choices = types, selected = character(), server = TRUE)
    shiny::updateNumericInput(session, "min_samples", max = ncol(analysis_data()$counts), value = min(2, ncol(analysis_data()$counts)))
    export_result(NULL)
  })
  output$length_status <- shiny::renderUI({
    if (!identical(input$length_source, "bundled")) return(NULL)
    path <- bundled_gene_lengths_path()
    if (nzchar(path)) shiny::div(class = "alert alert-info py-2", "Using bundled GENCODE v22 non-overlapping exonic gene lengths.")
    else shiny::div(class = "alert alert-warning py-2", "Bundled GENCODE v22 lengths were not found. Select Upload another annotation.")
  })
  observeEvent(input$generate, {
    x <- tryCatch({
      req(analysis_data())
      if (identical(input$length_source, "bundled")) {
        path <- bundled_gene_lengths_path()
        if (!nzchar(path)) stop("Bundled GENCODE v22 lengths were not found. Upload a matching gene-length table instead.", call. = FALSE)
        length_data <- read_tabular_file(path)
        lengths <- prepare_gene_lengths(length_data, "ensembl_id", "length_bp")
      } else {
        req(raw_lengths(), input$length_gene_id, input$length_value)
        lengths <- prepare_gene_lengths(raw_lengths(), input$length_gene_id, input$length_value)
      }
      compute_tpm_export(analysis_data()$counts, analysis_data()$annotations, lengths,
        input$min_count, input$min_samples, input$gene_types %||% character())
    }, error = function(e) { shiny::showNotification(e$message, type = "error", duration = 10); NULL })
    if (!is.null(x)) { export_result(x); shiny::showNotification("TPM and filtered counts generated.", type = "message") }
  })
  output$summary <- shiny::renderUI({ req(export_result()); p <- export_result()$parameters
    shiny::div(class = "alert alert-success", sprintf("Retained %s of %s genes across %s selected samples.", p$genes_retained, p$genes_input, ncol(export_result()$tpm)), shiny::br(),
      "Each TPM sample sums to approximately 1,000,000. Raw counts remain unchanged for DESeq2.") })
  output$tpm_preview <- DT::renderDT({ req(export_result()); DT::datatable(utils::head(tpm_export_table(export_result(), "tpm"), 100), rownames = FALSE, options = list(scrollX = TRUE, pageLength = 15)) })
  output$counts_preview <- DT::renderDT({ req(export_result()); DT::datatable(utils::head(tpm_export_table(export_result(), "counts"), 100), rownames = FALSE, options = list(scrollX = TRUE, pageLength = 15)) })
  output$download_tpm <- shiny::downloadHandler(function() "filtered_TPM.csv", function(file) write.csv(tpm_export_table(export_result(), "tpm"), file, row.names = FALSE, na = ""))
  output$download_counts <- shiny::downloadHandler(function() "filtered_raw_counts.csv", function(file) write.csv(tpm_export_table(export_result(), "counts"), file, row.names = FALSE, na = ""))
  output$download_workbook <- shiny::downloadHandler(function() "expression_export.xlsx", function(file) {
    p <- export_result()$parameters
    p$length_source <- if (identical(input$length_source, "bundled")) "GENCODE v22 non-overlapping exonic gene length" else input$length_file$name
    params <- data.frame(parameter = names(p), value = vapply(p, function(z) paste(z, collapse = ", "), character(1)))
    write_filterable_workbook(list(TPM = tpm_export_table(export_result(), "tpm"), Filtered_counts = tpm_export_table(export_result(), "counts"), Parameters = params), file)
  })
  export_result
})
