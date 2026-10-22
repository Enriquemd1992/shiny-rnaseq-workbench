mod_analysis_setup_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::fluidPage(
    shiny::h2("Analysis setup"),
    shiny::p("Configure one explicit DESeq2 contrast. Adjustment variables should precede the biological variable."),
    bslib::layout_sidebar(
      sidebar = bslib::sidebar(
        shiny::selectizeInput(ns("variables"), "Design variables (ordered)", choices=NULL, multiple=TRUE),
        shiny::selectInput(ns("contrast_var"), "Contrast variable", choices=NULL),
        shiny::selectInput(ns("numerator"), "Numerator", choices=NULL),
        shiny::selectInput(ns("denominator"), "Denominator / reference", choices=NULL),
        shiny::numericInput(ns("alpha"), "DESeq2 alpha", value=.05, min=.001, max=.25, step=.01),
        shiny::selectInput(ns("shrinkage"), "LFC shrinkage", choices=c("None"="none", "Normal prior (built-in)"="normal")),
        shiny::selectizeInput(ns("gene_types"), "Gene types (empty = all)", choices=NULL, multiple=TRUE),
        shiny::numericInput(ns("min_count"), "Minimum raw count", value=10, min=0, step=1),
        shiny::numericInput(ns("min_samples"), "Required number of samples", value=2, min=1, step=1),
        shiny::actionButton(ns("run"), "Run DESeq2", class="btn-primary")
      ),
      bslib::card(
        bslib::card_header("Analysis specification"),
        shiny::verbatimTextOutput(ns("formula")),
        shiny::uiOutput(ns("validation")),
        shiny::uiOutput(ns("filter_summary")),
        shiny::uiOutput(ns("interpretation")),
        shiny::uiOutput(ns("run_state"))
      )
    )
  )
}

mod_analysis_setup_server <- function(id, analysis_data) shiny::moduleServer(id, function(input, output, session) {
  result <- shiny::reactiveVal(NULL)

  observeEvent(analysis_data(), {
    md <- analysis_data()$metadata
    eligible <- names(md)[vapply(md, function(x) length(unique(x[!is.na(x)])) > 1L, logical(1))]
    shiny::updateSelectizeInput(session,"variables",choices=eligible,selected=NULL,server=TRUE)
    shiny::updateSelectInput(session,"contrast_var",choices=character())
    gene_types <- if("gene_type" %in% names(analysis_data()$annotations)) sort(unique(as.character(analysis_data()$annotations$gene_type))) else character()
    shiny::updateSelectizeInput(session,"gene_types",choices=gene_types,selected=character(),server=TRUE)
    shiny::updateNumericInput(session,"min_samples",max=ncol(analysis_data()$counts),value=min(2,ncol(analysis_data()$counts)))
  })
  observeEvent(input$variables, {
    req(analysis_data()); md <- analysis_data()$metadata
    categorical <- input$variables[input$variables %in% names(md) & vapply(md[input$variables], function(x) !is.numeric(x), logical(1))]
    shiny::updateSelectInput(session,"contrast_var",choices=categorical,selected=utils::tail(categorical,1))
  })
  observeEvent(input$contrast_var, {
    req(analysis_data(),input$contrast_var)
    lev <- unique(as.character(analysis_data()$metadata[[input$contrast_var]])); lev <- lev[!is.na(lev) & nzchar(trimws(lev))]
    shiny::updateSelectInput(session,"numerator",choices=lev,selected=if(length(lev)>1) lev[2] else lev[1])
    shiny::updateSelectInput(session,"denominator",choices=lev,selected=lev[1])
  })

  signature <- shiny::reactive(list(data=analysis_data(),variables=input$variables,contrast=input$contrast_var,numerator=input$numerator,denominator=input$denominator,alpha=input$alpha,shrinkage=input$shrinkage,gene_types=input$gene_types,min_count=input$min_count,min_samples=input$min_samples))
  observeEvent(signature(), { if (!is.null(result())) result(NULL) }, ignoreInit=TRUE)

  request_error <- shiny::reactive({
    req(analysis_data())
    tryCatch({ validate_deseq_request(analysis_data()$counts,analysis_data()$metadata,input$variables,input$contrast_var,input$numerator,input$denominator); NULL }, error=function(e)e$message)
  })

  output$formula <- shiny::renderText({ if(!length(input$variables)) "Design: not configured" else paste("Design:",paste(deparse(stats::reformulate(input$variables)),collapse="")) })
  output$validation <- shiny::renderUI({
    req(analysis_data()); error <- request_error()
    if(is.null(error)) shiny::div(class="alert alert-success","Design and contrast are valid.") else shiny::div(class="alert alert-danger",error)
  })
  output$filter_summary <- shiny::renderUI({
    req(analysis_data(),input$min_count,input$min_samples)
    keep <- tryCatch(gene_filter_mask(analysis_data()$counts,analysis_data()$annotations,input$gene_types %||% character(),input$min_count,input$min_samples),error=function(e)e)
    if(inherits(keep,"error")) shiny::div(class="alert alert-danger",keep$message)
    else shiny::div(class="alert alert-secondary",sprintf("Pre-DESeq2 filter: raw count ≥ %s in at least %s sample(s); gene type: %s. Retains %s of %s genes.",
      input$min_count,input$min_samples,if(length(input$gene_types))paste(input$gene_types,collapse=", ") else "all",format(sum(keep),big.mark=","),format(length(keep),big.mark=",")))
  })
  output$interpretation <- shiny::renderUI({
    req(input$contrast_var,input$numerator,input$denominator)
    shiny::div(class="alert alert-info",shiny::strong(input$numerator)," versus ",shiny::strong(input$denominator)," for ",input$contrast_var,
      ". Positive log2 fold changes mean higher expression in the numerator.")
  })
  output$run_state <- shiny::renderUI({
    if(is.null(result())) shiny::p(class="text-muted","No current DESeq2 result. Changes to metadata or this configuration invalidate the previous run.")
    else shiny::div(class="alert alert-success",shiny::strong("Completed: "),result()$parameters$executed_at)
  })

  observeEvent(input$run, {
    req(analysis_data())
    error <- request_error(); if(!is.null(error)) return(shiny::showNotification(error,type="error",duration=8))
    fit <- tryCatch(shiny::withProgress(message="Running DESeq2",value=.2,{
      ans <- run_deseq_analysis(analysis_data()$counts,analysis_data()$metadata,input$variables,input$contrast_var,input$numerator,input$denominator,input$alpha,input$shrinkage,
        annotations=analysis_data()$annotations,gene_types=input$gene_types %||% character(),min_count=input$min_count,min_samples=input$min_samples,sample_selection=analysis_data()$sample_selection)
        shiny::incProgress(.8,detail="Preparing results"); ans
    }),error=function(e){shiny::showNotification(paste("DESeq2 failed:",e$message),type="error",duration=NULL);NULL})
    if(!is.null(fit)){ result(fit); shiny::showNotification("DESeq2 analysis completed. Open Results.",type="message",duration=8) }
  })
  shiny::reactive(result())
})
