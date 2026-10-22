mod_metadata_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::fluidPage(
    shiny::h2("Metadata editor"),
    shiny::p("Edit a draft safely, review validation, then apply it to Design and QC."),
    bslib::layout_sidebar(
      sidebar = bslib::sidebar(
        shiny::selectInput(ns("column"), "Selected column", choices = NULL),
        shiny::div(class = "d-grid gap-2 mb-3",
          shiny::actionButton(ns("undo"), "Undo last change", class = "btn-outline-secondary"),
          shiny::actionButton(ns("discard"), "Discard pending changes", class = "btn-outline-secondary"),
          shiny::actionButton(ns("reset"), "Reset to imported metadata", class = "btn-outline-warning")
        ),
        shiny::hr(), shiny::h5("Column type"),
        shiny::selectInput(ns("type"), "Convert to", choices = c("Text" = "character", "Factor" = "factor", "Numeric" = "numeric")),
        shiny::actionButton(ns("convert"), "Convert type"),
        shiny::hr(), shiny::h5("Factor reference"),
        shiny::selectInput(ns("reference"), "Reference level", choices = NULL),
        shiny::actionButton(ns("set_reference"), "Set reference"),
        shiny::hr(), shiny::h5("Rename selected column"),
        shiny::textInput(ns("rename_to"), "New name", placeholder = "e.g. condition"),
        shiny::actionButton(ns("rename"), "Rename"),
        shiny::hr(), shiny::h5("Assign selected samples"),
        shiny::textInput(ns("bulk_value"), "Value"),
        shiny::actionButton(ns("bulk_assign"), "Assign to selected rows"),
        shiny::hr(), shiny::h5("Add or remove columns"),
        shiny::textInput(ns("new_name"), "New column name"),
        shiny::textInput(ns("new_value"), "Initial value"),
        shiny::actionButton(ns("add"), "Add column"),
        shiny::actionButton(ns("delete"), "Delete selected", class = "btn-outline-danger")
      ),
      bslib::navset_card_tab(
        bslib::nav_panel("Editor", shiny::uiOutput(ns("state_banner")), DT::DTOutput(ns("table")),
          shiny::div(class = "d-flex gap-2 mt-3", shiny::actionButton(ns("apply"), "Apply changes", class = "btn-success"), shiny::uiOutput(ns("summary")))),
        bslib::nav_panel("Validation", shiny::uiOutput(ns("validation"))),
        bslib::nav_panel("Column details", DT::DTOutput(ns("details")))
      )
    )
  )
}

mod_metadata_server <- function(id, imported) shiny::moduleServer(id, function(input, output, session) {
  original <- shiny::reactiveVal(); committed <- shiny::reactiveVal(); draft <- shiny::reactiveVal(); history <- shiny::reactiveVal(list())
  is_dirty <- shiny::reactive({ req(draft(), committed()); !identical(draft(), committed()) })

  update_column_controls <- function(selected = NULL) {
    x <- draft(); if (is.null(x)) return(invisible())
    choices <- names(x); selected <- selected %||% input$column %||% choices[1]
    if (!selected %in% choices) selected <- choices[1]
    shiny::updateSelectInput(session, "column", choices = choices, selected = selected)
  }
  push_draft <- function(value, selected = NULL) {
    history(c(history(), list(draft()))); draft(value); update_column_controls(selected)
  }
  notify_error <- function(expr) tryCatch(expr, error = function(e) { shiny::showNotification(e$message, type = "error", duration = 7); NULL })

  observeEvent(imported$accepted(), {
    x <- imported$accepted()$metadata
    original(x); committed(x); draft(x); history(list()); update_column_controls(names(x)[1])
  }, ignoreNULL = TRUE)

  observeEvent(list(input$column, draft()), {
    req(draft(), input$column); if (!input$column %in% names(draft())) return()
    x <- draft()[[input$column]]; values <- unique(as.character(x)[!is.na(x) & nzchar(trimws(as.character(x)))])
    reference <- if (is.factor(x) && length(levels(x))) levels(x)[1] else values[1]
    shiny::updateSelectInput(session, "reference", choices = values, selected = reference)
    shiny::updateSelectInput(session, "type", selected = metadata_type_label(x))
  }, ignoreInit = TRUE)

  output$table <- DT::renderDT({
    req(draft())
    DT::datatable(cbind(sample_id = rownames(draft()), draft()), rownames = FALSE,
      selection = list(mode = "multiple", target = "row"), editable = list(target = "cell", disable = list(columns = 0)),
      options = list(scrollX = TRUE, pageLength = 15))
  })

  observeEvent(input$table_cell_edit, {
    req(draft()); info <- input$table_cell_edit; column_index <- info$col
    if (column_index < 1L || column_index > ncol(draft())) return()
    column <- names(draft())[column_index]; updated <- notify_error(edit_metadata_cell(draft(), info$row, column, info$value))
    if (!is.null(updated)) push_draft(updated, column)
  })
  observeEvent(input$convert, {
    req(draft(), input$column, input$type); updated <- notify_error(convert_metadata_column(draft(), input$column, input$type))
    if (!is.null(updated)) push_draft(updated, input$column)
  })
  observeEvent(input$set_reference, {
    req(draft(), input$column, input$reference); updated <- notify_error(set_metadata_reference(draft(), input$column, input$reference))
    if (!is.null(updated)) { push_draft(updated, input$column); shiny::showNotification(sprintf("'%s' is now the reference for '%s'.", input$reference, input$column), type = "message") }
  })
  observeEvent(input$rename, {
    req(draft(), input$column); updated <- notify_error(rename_metadata_column(draft(), input$column, input$rename_to))
    if (!is.null(updated)) { new_name <- trimws(input$rename_to); push_draft(updated, new_name); shiny::updateTextInput(session, "rename_to", value = "") }
  })
  observeEvent(input$bulk_assign, {
    req(draft(), input$column); updated <- notify_error(assign_metadata_value(draft(), input$table_rows_selected, input$column, input$bulk_value))
    if (!is.null(updated)) push_draft(updated, input$column)
  })
  observeEvent(input$add, {
    req(draft()); updated <- notify_error(add_metadata_column(draft(), input$new_name, input$new_value))
    if (!is.null(updated)) { new_name <- trimws(input$new_name); push_draft(updated, new_name); shiny::updateTextInput(session, "new_name", value = ""); shiny::updateTextInput(session, "new_value", value = "") }
  })
  observeEvent(input$delete, {
    req(draft(), input$column); updated <- notify_error(delete_metadata_column(draft(), input$column))
    if (!is.null(updated)) push_draft(updated, names(updated)[1])
  })
  observeEvent(input$undo, {
    stack <- history(); if (!length(stack)) return(shiny::showNotification("There is no change to undo.", type = "warning"))
    draft(stack[[length(stack)]]); history(stack[-length(stack)]); update_column_controls()
  })
  observeEvent(input$discard, {
    req(committed()); draft(committed()); history(list()); update_column_controls(); shiny::showNotification("Pending changes discarded.", type = "message")
  })
  observeEvent(input$reset, {
    req(original()); push_draft(original(), names(original())[1]); shiny::showNotification("Imported metadata restored in the draft. Click Apply changes to use it.", type = "warning")
  })
  observeEvent(input$apply, {
    req(draft()); issues <- validate_metadata(draft())
    if (has_errors(issues)) return(shiny::showNotification("Resolve metadata errors before applying changes.", type = "error"))
    committed(draft()); history(list()); n_warnings <- sum(issues$severity == "warning")
    message <- if (n_warnings) sprintf("Metadata applied with %d warning(s).", n_warnings) else "Metadata applied successfully."
    shiny::showNotification(message, type = if (n_warnings) "warning" else "message", duration = 7)
  })

  output$state_banner <- shiny::renderUI({
    req(draft(), committed())
    if (is_dirty()) shiny::div(class = "alert alert-warning", shiny::strong("Pending changes."), " Design and QC still use the last applied metadata.")
    else shiny::div(class = "alert alert-success", "Metadata is applied and up to date.")
  })
  output$summary <- shiny::renderUI({ req(draft()); shiny::span(class = "align-self-center text-muted", sprintf("%d samples × %d variables", nrow(draft()), ncol(draft()))) })
  output$validation <- shiny::renderUI({
    req(draft()); issues <- validate_metadata(draft())
    if (!nrow(issues)) return(shiny::div(class = "alert alert-success", "No metadata warnings or errors."))
    shiny::tagList(lapply(seq_len(nrow(issues)), function(i) {
      cls <- switch(issues$severity[i], error = "alert-danger", warning = "alert-warning", "alert-info")
      shiny::div(class = paste("alert", cls), shiny::strong(toupper(issues$severity[i])), " — ", issues$message[i])
    }))
  })
  output$details <- DT::renderDT({
    req(draft())
    details <- do.call(rbind, lapply(names(draft()), function(nm) {
      x <- draft()[[nm]]; observed <- unique(as.character(x)[!is.na(x) & nzchar(trimws(as.character(x)))])
      data.frame(column=nm, type=metadata_type_label(x), reference=if (is.factor(x) && length(levels(x))) levels(x)[1] else "—",
        unique_values=length(observed), missing=sum(is.na(x) | trimws(as.character(x)) == ""), values=paste(utils::head(observed, 8), collapse=", "), stringsAsFactors=FALSE)
    }))
    DT::datatable(details, rownames = FALSE, options = list(scrollX = TRUE, pageLength = 15))
  })

  shiny::reactive({ req(committed()); committed() })
})
