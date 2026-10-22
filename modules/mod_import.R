mod_import_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::fluidPage(
    shiny::h2("Import and validate"),
    shiny::p("Load raw integer counts and sample metadata. Files remain on this computer."),
    bslib::layout_sidebar(
      sidebar = bslib::sidebar(
        shiny::fileInput(ns("counts_file"), "Raw counts", accept = c(".csv", ".tsv", ".txt", ".xlsx")),
        shiny::selectInput(ns("gene_col"), "Gene identifier column", choices = NULL),
        shiny::selectizeInput(ns("annotation_cols"), "Gene annotation columns (excluded from counts)", choices = NULL, multiple = TRUE),
        shiny::fileInput(ns("metadata_file"), "Metadata", accept = c(".csv", ".tsv", ".txt", ".xlsx")),
        shiny::selectInput(ns("sample_col"), "Sample identifier column", choices = NULL),
        shiny::actionButton(ns("validate"), "Validate", class = "btn-primary"),
        shiny::actionButton(ns("reorder"), "Reorder metadata", class = "btn-outline-secondary"),
        shiny::actionButton(ns("accept"), "Accept data", class = "btn-success")
      ),
      bslib::navset_card_tab(
        id = ns("import_tabs"),
        bslib::nav_panel("Counts preview", value = "counts", DT::DTOutput(ns("counts_preview"))),
        bslib::nav_panel("Metadata preview", value = "metadata", DT::DTOutput(ns("metadata_preview"))),
        bslib::nav_panel("Sample selection", value = "samples",
          shiny::p("Optional: define the analysis cohort after validation. Counts and metadata are always filtered together."),
          shiny::radioButtons(ns("selection_mode"),"Samples to use",choices=c("All samples"="all","Filter by metadata"="metadata","Select manually"="manual"),selected="all",inline=TRUE),
          bslib::layout_columns(
            shiny::selectInput(ns("filter_variable"),"Metadata variable",choices=NULL),
            shiny::div(
              shiny::div(class="d-flex justify-content-between align-items-center",
                shiny::strong("Values to include"),
                shiny::div(
                  shiny::actionLink(ns("select_all_values"), "Select all"),
                  " · ",
                  shiny::actionLink(ns("clear_filter_values"), "Clear")
                )
              ),
              shiny::div(
                style="max-height: 240px; overflow-y: auto; padding: 0.5rem 0.75rem; border: 1px solid #dee2e6; border-radius: 0.375rem;",
                shiny::checkboxGroupInput(ns("filter_values"), label=NULL, choices=NULL, selected=character())
              )
            ),
            col_widths=c(6,6)
          ),
          shiny::selectizeInput(ns("manual_samples"),"Manual sample selection",choices=NULL,multiple=TRUE),
          shiny::div(class="d-flex gap-2 mb-3",shiny::actionButton(ns("apply_selection"),"Apply sample selection",class="btn-primary"),shiny::actionButton(ns("reset_selection"),"Reset to all samples",class="btn-outline-secondary")),
          shiny::uiOutput(ns("selection_summary")),DT::DTOutput(ns("selection_preview"))
        ),
        bslib::nav_panel("Validation", value = "validation", shiny::uiOutput(ns("issues")))
      )
    )
  )
}

mod_import_server <- function(id) shiny::moduleServer(id, function(input, output, session) {
  raw_counts <- shiny::reactiveVal(); raw_metadata <- shiny::reactiveVal()
  prepared_full <- shiny::reactiveVal(); prepared <- shiny::reactiveVal(); accepted <- shiny::reactiveVal()
  observeEvent(input$counts_file, {
    x <- tryCatch(read_tabular_file(input$counts_file$datapath, input$counts_file$name), error = function(e) { shiny::showNotification(e$message, type="error"); NULL })
    raw_counts(x); if (!is.null(x)) {
      gene_guess <- if ("ensembl_id" %in% names(x)) "ensembl_id" else names(x)[1]
      shiny::updateSelectInput(session, "gene_col", choices=names(x), selected=gene_guess)
      shiny::updateSelectizeInput(session,"annotation_cols",choices=setdiff(names(x),gene_guess),selected=suggest_annotation_columns(x,gene_guess),server=TRUE)
    }
  })
  observeEvent(input$gene_col, {
    req(raw_counts(),input$gene_col)
    shiny::updateSelectizeInput(session,"annotation_cols",choices=setdiff(names(raw_counts()),input$gene_col),selected=suggest_annotation_columns(raw_counts(),input$gene_col),server=TRUE)
  },ignoreInit=TRUE)
  observeEvent(input$metadata_file, {
    x <- tryCatch(read_tabular_file(input$metadata_file$datapath, input$metadata_file$name), error = function(e) { shiny::showNotification(e$message, type="error"); NULL })
    raw_metadata(x); if (!is.null(x)) shiny::updateSelectInput(session, "sample_col", choices=names(x), selected=names(x)[1])
  })
  output$counts_preview <- DT::renderDT({ req(raw_counts()); DT::datatable(utils::head(raw_counts(), 20), options=list(scrollX=TRUE)) })
  output$metadata_preview <- DT::renderDT({ req(raw_metadata()); DT::datatable(raw_metadata(), options=list(scrollX=TRUE, pageLength=10)) })
  observeEvent(input$validate, {
    req(raw_counts(), raw_metadata(), input$gene_col, input$sample_col)
    result <- tryCatch({
      prepared_counts <- prepare_count_data(raw_counts(), input$gene_col, input$annotation_cols %||% character())
      counts <- prepared_counts$counts
      metadata <- prepare_metadata(raw_metadata(), input$sample_col)
      list(counts=counts, annotations=prepared_counts$annotations, metadata=metadata, gene_id_column=input$gene_col,
        annotation_columns=input$annotation_cols %||% character(),sample_selection=list(mode="all",variable="",values=character(),included=colnames(counts),excluded=character()),
        issues=validate_all(counts, metadata, prepared_counts$annotations))
    }, error=function(e) { shiny::showNotification(e$message, type="error"); NULL })
    prepared_full(result);prepared(result)
    if (!is.null(result)) {
      fields<-names(result$metadata);shiny::updateSelectInput(session,"filter_variable",choices=fields,selected=if("experiment"%in%fields)"experiment" else fields[1])
      shiny::updateSelectizeInput(session,"manual_samples",choices=rownames(result$metadata),selected=rownames(result$metadata),server=TRUE)
      bslib::nav_select("import_tabs", "validation", session = session)
      n_errors <- sum(result$issues$severity == "error")
      n_warnings <- sum(result$issues$severity == "warning")
      if (n_errors > 0) {
        shiny::showNotification(sprintf("Validation finished: %d error(s), %d warning(s).", n_errors, n_warnings), type = "error")
      } else if (n_warnings > 0) {
        shiny::showNotification(sprintf("Validation finished: no errors, %d warning(s).", n_warnings), type = "warning")
      } else {
        shiny::showNotification("Validation passed: no errors or warnings.", type = "message")
      }
    }
  })
  observeEvent(input$filter_variable,{
    req(prepared_full(),input$filter_variable);x<-prepared_full()$metadata[[input$filter_variable]];values<-unique(as.character(x));values<-values[!is.na(values)];counts<-table(as.character(x),useNA="no");labels<-paste0(values," (",as.integer(counts[values])," samples)")
    shiny::updateCheckboxGroupInput(session,"filter_values",choices=stats::setNames(values,labels),selected=character())
  })
  observeEvent(input$select_all_values, {
    req(prepared_full(), input$filter_variable)
    x <- prepared_full()$metadata[[input$filter_variable]]
    values <- unique(as.character(x)); values <- values[!is.na(values)]
    counts <- table(as.character(x), useNA="no")
    labels <- paste0(values, " (", as.integer(counts[values]), " samples)")
    shiny::updateCheckboxGroupInput(session,"filter_values",choices=stats::setNames(values,labels),selected=values)
  })
  observeEvent(input$clear_filter_values, {
    shiny::updateCheckboxGroupInput(session,"filter_values",selected=character())
  })
  observeEvent(input$apply_selection,{
    req(prepared_full());full<-prepared_full()
    selected_ids<-tryCatch(select_sample_ids(full$metadata,input$selection_mode,input$filter_variable,input$filter_values %||% character(),input$manual_samples %||% character()),error=function(e){shiny::showNotification(e$message,type="error");NULL})
    if(is.null(selected_ids))return()
    subset<-subset_samples(full$counts,full$metadata,selected_ids);x<-full;x$counts<-subset$counts;x$metadata<-subset$metadata
    x$sample_selection<-list(mode=input$selection_mode,variable=if(input$selection_mode=="metadata")input$filter_variable else "",values=if(input$selection_mode=="metadata")input$filter_values else character(),included=colnames(x$counts),excluded=setdiff(colnames(full$counts),colnames(x$counts)))
    x$issues<-validate_all(x$counts,x$metadata,x$annotations);prepared(x);accepted(NULL);bslib::nav_select("import_tabs","samples",session=session)
    shiny::showNotification(sprintf("Sample selection applied: %d included, %d excluded. Validate the summary, then Accept data.",ncol(x$counts),length(x$sample_selection$excluded)),type="message",duration=8)
  })
  observeEvent(input$reset_selection,{
    req(prepared_full());x<-prepared_full();x$sample_selection<-list(mode="all",variable="",values=character(),included=colnames(x$counts),excluded=character());prepared(x);accepted(NULL);shiny::updateRadioButtons(session,"selection_mode",selected="all");shiny::showNotification("All matched samples restored.",type="message")
  })
  observeEvent(input$reorder, {
    req(prepared())
    x <- tryCatch({ prepared()$metadata <- align_metadata(prepared()$counts, prepared()$metadata); prepared() }, error=function(e) { shiny::showNotification(e$message, type="error"); NULL })
    if (!is.null(x)) { x$issues <- validate_all(x$counts, x$metadata, x$annotations); prepared(x);if(identical(x$sample_selection$mode,"all"))prepared_full(x) }
  })
  observeEvent(input$accept, {
    req(prepared())
    if (has_errors(prepared()$issues)) shiny::showNotification("Resolve validation errors before accepting.", type="error")
    else { accepted(prepared()); shiny::showNotification("Data accepted. Continue to Metadata.", type="message") }
  })
  output$issues <- shiny::renderUI({
    req(prepared()); x <- prepared()$issues
    shiny::tagList(lapply(seq_len(nrow(x)), function(i) {
      cls <- switch(x$severity[i], error="alert-danger", warning="alert-warning", ok="alert-success")
      shiny::div(class=paste("alert", cls), shiny::strong(toupper(x$severity[i])), " — ", x$message[i])
    }))
  })
  output$selection_summary<-shiny::renderUI({
    req(prepared());s<-prepared()$sample_selection
    description<-if(s$mode=="all")"All matched samples" else if(s$mode=="metadata")paste0(s$variable," ∈ {",paste(s$values,collapse=", "),"}") else "Manual selection"
    shiny::div(class="alert alert-info",shiny::strong(description),shiny::br(),sprintf("Included: %d. Excluded: %d.",length(s$included),length(s$excluded)),if(length(s$excluded))shiny::tagList(shiny::br(),"Excluded: ",paste(s$excluded,collapse=", ")))
  })
  output$selection_preview<-DT::renderDT({req(prepared());DT::datatable(cbind(sample_id=rownames(prepared()$metadata),prepared()$metadata),rownames=FALSE,options=list(scrollX=TRUE,pageLength=10))})
  list(accepted=accepted)
})
