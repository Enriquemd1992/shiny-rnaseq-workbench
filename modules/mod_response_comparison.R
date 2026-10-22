mod_response_comparison_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::fluidPage(
    shiny::h2("Response comparison"),
    shiny::p("Compare one treatment response across two biological groups using a joint DESeq2 interaction model."),
    bslib::layout_sidebar(
      sidebar = bslib::sidebar(
        shiny::selectInput(ns("group_variable"), "Group variable", choices = NULL),
        shiny::selectInput(ns("group_a"), "Group A", choices = NULL),
        shiny::selectInput(ns("group_b"), "Group B", choices = NULL),
        shiny::selectInput(ns("treatment_variable"), "Treatment variable", choices = NULL),
        shiny::selectInput(ns("treated"), "Treated level", choices = NULL),
        shiny::selectInput(ns("control"), "Control level", choices = NULL),
        shiny::numericInput(ns("alpha"), "FDR threshold", .05, min = .001, max = .25, step = .01),
        shiny::numericInput(ns("lfc_threshold"), "Absolute log2FC threshold", 1, min = 0, step = .25),
        shiny::numericInput(ns("min_count"), "Minimum count", 10, min = 0, step = 1),
        shiny::numericInput(ns("min_samples"), "Minimum samples", 2, min = 1, step = 1),
        shiny::actionButton(ns("run"), "Run response comparison", class = "btn-primary"),
        shiny::hr(), shiny::downloadButton(ns("download_table"), "Download classified table"),
        shiny::downloadButton(ns("download_scatter"), "Download scatter PNG"),
        shiny::downloadButton(ns("download_interaction"), "Download interaction PNG")
      ),
      bslib::navset_card_tab(
        bslib::nav_panel("Response scatter", shiny::plotOutput(ns("scatter"), height = "700px")),
        bslib::nav_panel("Categories",
          shiny::div(class="alert alert-warning",
            shiny::strong("Descriptive threshold overlap: "),
            "'Detected only' does not prove a different response. Use interaction_padj to test the between-group difference."),
          DT::DTOutput(ns("category_table"))),
        bslib::nav_panel("Interaction volcano", shiny::plotOutput(ns("interaction"), height = "650px")),
        bslib::nav_panel("Genes", DT::DTOutput(ns("genes")))
      )
    )
  )
}

mod_response_comparison_server <- function(id, analysis_data) shiny::moduleServer(id, function(input, output, session) {
  result <- shiny::reactiveVal(NULL)
  observeEvent(analysis_data(), {
    req(analysis_data()); fields <- names(analysis_data()$metadata)
    shiny::updateSelectInput(session, "group_variable", choices = fields, selected = if ("cell_model" %in% fields) "cell_model" else fields[1])
    shiny::updateSelectInput(session, "treatment_variable", choices = fields, selected = if ("treatment" %in% fields) "treatment" else fields[min(2, length(fields))])
    result(NULL)
  })
  observeEvent(input$group_variable, {
    req(analysis_data(), input$group_variable); levels <- unique(as.character(analysis_data()$metadata[[input$group_variable]])); levels <- levels[!is.na(levels) & nzchar(levels)]
    shiny::updateSelectInput(session, "group_a", choices = levels, selected = levels[1])
    shiny::updateSelectInput(session, "group_b", choices = levels, selected = levels[min(2, length(levels))])
  })
  observeEvent(input$treatment_variable, {
    req(analysis_data(), input$treatment_variable); levels <- unique(as.character(analysis_data()$metadata[[input$treatment_variable]])); levels <- levels[!is.na(levels) & nzchar(levels)]
    shiny::updateSelectInput(session, "treated", choices = levels, selected = if ("IFN-y" %in% levels) "IFN-y" else levels[min(2, length(levels))])
    shiny::updateSelectInput(session, "control", choices = levels, selected = if ("UT" %in% levels) "UT" else levels[1])
  })
  observeEvent(input$run, {
    req(analysis_data(), input$group_a, input$group_b, input$treated, input$control)
    x <- tryCatch(run_response_comparison(
      analysis_data()$counts, analysis_data()$metadata, analysis_data()$annotations,
      input$group_variable, input$group_a, input$group_b, input$treatment_variable, input$treated, input$control,
      input$alpha, input$lfc_threshold, input$min_count, input$min_samples),
      error = function(e) { shiny::showNotification(e$message, type = "error", duration = 10); NULL })
    if (!is.null(x)) { result(x); shiny::showNotification("Response comparison completed.", type = "message") }
  })
  output$scatter <- shiny::renderPlot({ req(result()); response_scatter_plot(result()) })
  output$interaction <- shiny::renderPlot({ req(result()); interaction_volcano_plot(result()) })
  output$category_table <- DT::renderDT({ req(result()); DT::datatable(response_category_summary(result()), rownames = FALSE, options = list(pageLength = 15)) })
  output$genes <- DT::renderDT({ req(result()); x <- result()$table; keep <- intersect(c("gene_id", "gene_symbol", "gene_type", "lfc_group_a", "padj_group_a", "lfc_group_b", "padj_group_b", "interaction_lfc", "interaction_padj", "response_category", "response_difference_significant"), names(x)); DT::datatable(x[, keep, drop = FALSE], rownames = FALSE, filter = "top", options = list(scrollX = TRUE, pageLength = 25)) })
  output$download_table <- shiny::downloadHandler(function() "response_comparison.csv", function(file) utils::write.csv(result()$table, file, row.names = FALSE, na = ""))
  png_download <- function(filename, plotter) shiny::downloadHandler(function() filename, function(file) { grDevices::png(file, 1800, 1500, res = 200); on.exit(grDevices::dev.off()); print(plotter()) })
  output$download_scatter <- png_download("response_scatter.png", function() response_scatter_plot(result()))
  output$download_interaction <- png_download("interaction_volcano.png", function() interaction_volcano_plot(result()))
  result
})
