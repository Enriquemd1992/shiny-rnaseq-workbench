mod_results_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::fluidPage(
    shiny::h2("DESeq2 results"),
    shiny::uiOutput(ns("empty")),
    shiny::conditionalPanel(sprintf("output['%s']",ns("has_results")),
      bslib::layout_sidebar(
        sidebar=bslib::sidebar(
          shiny::numericInput(ns("display_alpha"),"FDR threshold",.05,min=.001,max=.25,step=.01),
          shiny::numericInput(ns("lfc_threshold"),"Absolute log2FC threshold",1,min=0,step=.25),
          shiny::numericInput(ns("top_n"),"Top genes to display",50,min=1,max=1000,step=1),
          shiny::selectInput(ns("ranking"),"Rank by",choices=c("Adjusted p-value"="padj","Absolute log2FC"="lfc","Combined evidence"="combined")),
          shiny::checkboxInput(ns("significant_only"),"Only significant genes",FALSE),
          shiny::selectizeInput(
            ns("gene"), "Search gene symbol or Ensembl ID", choices=NULL,
            options=list(
              placeholder="Type a symbol (for example MYC) or an Ensembl ID...",
              openOnFocus=FALSE, maxOptions=50, closeAfterSelect=TRUE,
              plugins=list("clear_button")
            )
          ),
          shiny::helpText("The selected value is always the unique gene ID. Duplicate symbols are shown separately and are never collapsed."),
          shiny::hr(),
          shiny::downloadButton(ns("download_csv"),"Download full CSV"),
          shiny::downloadButton(ns("download_xlsx"),"Download workbook"),
          shiny::downloadButton(ns("download_record"),"Download analysis record")
        ),
        bslib::navset_card_tab(
          bslib::nav_panel("Summary",shiny::uiOutput(ns("summary")),shiny::verbatimTextOutput(ns("specification"))),
          bslib::nav_panel("Top genes",DT::DTOutput(ns("table"))),
          bslib::nav_panel("Volcano",shiny::plotOutput(ns("volcano"),height="600px"),shiny::downloadButton(ns("download_volcano"),"Download PNG")),
          bslib::nav_panel("MA plot",shiny::plotOutput(ns("ma"),height="600px"),shiny::downloadButton(ns("download_ma"),"Download PNG")),
          bslib::nav_panel("Heatmap",shiny::plotOutput(ns("heatmap"),height="700px"),shiny::downloadButton(ns("download_heatmap"),"Download PNG")),
          bslib::nav_panel("Gene expression",shiny::uiOutput(ns("gene_prompt")),shiny::plotOutput(ns("gene_plot"),height="550px"),shiny::downloadButton(ns("download_gene"),"Download PNG"))
        )
      )
    )
  )
}

mod_results_server <- function(id, analysis_result) shiny::moduleServer(id,function(input,output,session){
  output$has_results <- shiny::reactive(!is.null(analysis_result())); shiny::outputOptions(output,"has_results",suspendWhenHidden=FALSE)
  output$empty <- shiny::renderUI({ if(is.null(analysis_result())) shiny::div(class="alert alert-info","Run a valid contrast in Analysis setup to populate this page.") })
  observeEvent(analysis_result(),{
    req(analysis_result()); p <- analysis_result()$parameters
    shiny::updateNumericInput(session,"display_alpha",value=p$alpha)
    tab<-analysis_result()$results; genes<-tab$gene_id
    labels<-if("gene_symbol"%in%names(tab)) paste0(ifelse(is.na(tab$gene_symbol)|!nzchar(tab$gene_symbol),"(no symbol)",tab$gene_symbol)," — ",genes) else genes
    shiny::updateSelectizeInput(session,"gene",choices=stats::setNames(genes,labels),selected=character(),server=TRUE)
  })
  classified <- shiny::reactive({ req(analysis_result()); classify_deseq_results(analysis_result()$results,input$display_alpha,input$lfc_threshold,TRUE) })
  ranked <- shiny::reactive({
    x <- classified(); if(isTRUE(input$significant_only)) x <- x[x$direction!="Not significant",,drop=FALSE]
    utils::head(rank_deseq_results(x,input$ranking),input$top_n)
  })
  output$summary <- shiny::renderUI({
    s <- deseq_summary(classified())
    bslib::layout_columns(
      bslib::value_box("Upregulated",s[["up"]],theme="danger"),
      bslib::value_box("Downregulated",s[["down"]],theme="primary"),
      bslib::value_box("Not significant",s[["not_significant"]],theme="secondary"),col_widths=c(4,4,4))
  })
  output$specification <- shiny::renderText({ p<-analysis_result()$parameters; paste("Design:",p$design,"\nContrast:",p$numerator,"versus",p$denominator,"\nShrinkage:",p$shrinkage,"\nExecuted:",p$executed_at) })
  output$table <- DT::renderDT({
    x<-ranked(); annotation_cols<-setdiff(names(analysis_result()$annotations),"gene_id"); keep<-intersect(c("gene_id",annotation_cols,"baseMean","log2FoldChange","shrunkenLog2FoldChange","lfcSE","stat","pvalue","padj","direction"),names(x))
    widget <- DT::datatable(x[,keep,drop=FALSE],rownames=FALSE,filter="top",options=list(scrollX=TRUE,pageLength=25))
    DT::formatSignif(widget,columns=intersect(c("baseMean","log2FoldChange","shrunkenLog2FoldChange","lfcSE","stat","pvalue","padj"),keep),digits=4)
  })
  volcano <- shiny::reactive(volcano_plot(classified(),input$display_alpha,input$lfc_threshold))
  ma <- shiny::reactive(ma_plot(classified(),input$display_alpha,input$lfc_threshold))
  output$volcano <- shiny::renderPlot(volcano())
  output$ma <- shiny::renderPlot(ma())
  heatmap_draw <- function(){
    req(analysis_result()); genes<-intersect(ranked()$gene_id,rownames(analysis_result()$transformed)); genes<-utils::head(genes,100)
    shiny::validate(shiny::need(length(genes)>=2,"At least two genes are required for the heatmap."))
    mat<-analysis_result()$transformed[genes,,drop=FALSE]
    if("gene_symbol"%in%names(analysis_result()$annotations)){symbols<-as.character(analysis_result()$annotations[genes,"gene_symbol"]);symbols[is.na(symbols)|!nzchar(symbols)]<-genes[is.na(symbols)|!nzchar(symbols)];rownames(mat)<-make.unique(paste0(symbols," | ",genes))}
    ann<-analysis_result()$metadata[,analysis_result()$parameters$contrast_variable,drop=FALSE]
    pheatmap::pheatmap(mat,scale="row",annotation_col=ann,show_rownames=length(genes)<=60,fontsize_row=7)
  }
  output$heatmap <- shiny::renderPlot(heatmap_draw())
  output$gene_prompt <- shiny::renderUI({
    if (is.null(input$gene) || !nzchar(input$gene))
      shiny::div(class="alert alert-info", "Search for a gene by symbol or Ensembl ID in the left panel to draw its expression plot.")
  })
  gene_plot <- shiny::reactive({
    req(analysis_result(),input$gene); p<-analysis_result()$parameters
    shiny::validate(shiny::need(input$gene %in% rownames(analysis_result()$normalized_counts),"Select a valid gene."))
    df<-data.frame(sample=colnames(analysis_result()$normalized_counts),expression=as.numeric(analysis_result()$normalized_counts[input$gene,]),group=analysis_result()$metadata[[p$contrast_variable]])
    symbol<-if("gene_symbol"%in%names(analysis_result()$annotations))as.character(analysis_result()$annotations[input$gene,"gene_symbol"]) else ""
    title<-if(!is.na(symbol)&&nzchar(symbol))paste0(symbol," | ",input$gene) else input$gene
    ggplot2::ggplot(df,ggplot2::aes(group,expression,color=group))+ggplot2::geom_boxplot(outlier.shape=NA,alpha=.15)+ggplot2::geom_point(position=ggplot2::position_jitter(width=.12),size=2.5)+ggplot2::labs(title=title,x=p$contrast_variable,y="DESeq2 normalized count",color=NULL)+ggplot2::theme_minimal()
  })
  output$gene_plot <- shiny::renderPlot(gene_plot())

  output$download_csv <- shiny::downloadHandler(function()"deseq2_full_results.csv",function(file)utils::write.csv(classified(),file,row.names=FALSE,na=""))
  output$download_xlsx <- shiny::downloadHandler(function()"deseq2_results.xlsx",function(file){
    p<-analysis_result()$parameters; params<-data.frame(parameter=names(p)[!vapply(p,is.list,logical(1))],value=vapply(p[!vapply(p,is.list,logical(1))],function(x)paste(x,collapse=", "),character(1)))
    write_filterable_workbook(list(Full_results=classified(),Displayed_top=ranked(),Parameters=params),file)
  })
  output$download_record <- shiny::downloadHandler(function()"analysis_record.txt",function(file){
    p<-analysis_result()$parameters; lines<-c("Shiny RNA-seq workbench analysis record",paste0("Generated: ",format(Sys.time(),"%Y-%m-%d %H:%M:%S %Z")),"",vapply(names(p),function(nm)paste0(nm,": ",paste(p[[nm]],collapse=", ")),character(1)),"",capture.output(utils::sessionInfo()))
    writeLines(lines,file)
  })
  png_download <- function(filename,plotter,width=1800,height=1400) shiny::downloadHandler(function()filename,function(file){grDevices::png(file,width=width,height=height,res=180);on.exit(grDevices::dev.off());print(plotter())})
  output$download_volcano <- png_download("volcano.png",volcano)
  output$download_ma <- png_download("ma_plot.png",ma)
  output$download_gene <- png_download("gene_expression.png",gene_plot)
  output$download_heatmap <- shiny::downloadHandler(function()"heatmap.png",function(file){grDevices::png(file,width=1800,height=1800,res=180);on.exit(grDevices::dev.off());heatmap_draw()})
})
