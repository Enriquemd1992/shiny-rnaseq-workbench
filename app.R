options(shiny.maxRequestSize = 100 * 1024^2)

source("R/io.R")
source("R/validation.R")
source("R/metadata.R")
source("R/analysis.R")
source("R/deseq_analysis.R")
source("R/response_comparison.R")
source("R/tpm.R")
source("R/excel_export.R")
source("modules/mod_import.R")
source("modules/mod_metadata.R")
source("modules/mod_design.R")
source("modules/mod_qc.R")
source("modules/mod_results.R")
source("modules/mod_response_comparison.R")
source("modules/mod_expression_export.R")

ui <- bslib::page_navbar(
  title = "Lab RNA-seq",
  theme = bslib::bs_theme(version = 5, bootswatch = "flatly"),
  bslib::nav_panel("1. Import", mod_import_ui("import")),
  bslib::nav_panel("2. Metadata", mod_metadata_ui("metadata")),
  bslib::nav_panel("3. QC", mod_qc_ui("qc")),
  bslib::nav_panel("4. Analysis setup", mod_analysis_setup_ui("setup")),
  bslib::nav_panel("5. Results", mod_results_ui("results")),
  bslib::nav_panel("6. Response comparison", mod_response_comparison_ui("response")),
  bslib::nav_panel("7. Expression export", mod_expression_export_ui("expression")),
  bslib::nav_spacer(),
  bslib::nav_item(shiny::span("Local only", class = "badge text-bg-success")),
  footer = shiny::div(class = "container-fluid small text-muted py-2",
    "All files stay in this R session. No network upload is performed.")
)

server <- function(input, output, session) {
  imported <- mod_import_server("import")
  edited_metadata <- mod_metadata_server("metadata", imported)
  analysis_data <- shiny::reactive({
    req(imported$accepted())
    list(counts = imported$accepted()$counts, annotations = imported$accepted()$annotations, metadata = edited_metadata(), sample_selection = imported$accepted()$sample_selection)
  })
  mod_qc_server("qc", analysis_data)
  analysis_result <- mod_analysis_setup_server("setup", analysis_data)
  mod_results_server("results", analysis_result)
  mod_response_comparison_server("response", analysis_data)
  mod_expression_export_server("expression", analysis_data)
}

shiny::shinyApp(ui, server)
