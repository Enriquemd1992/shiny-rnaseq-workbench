mod_qc_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::fluidPage(
    shiny::h2("Exploration and QC"),
    shiny::p("QC settings affect visualizations only. No sample is removed and DESeq2 settings are not changed."),
    bslib::card(
      bslib::card_header("PCA and sample-similarity settings"),
      bslib::layout_columns(
        shiny::selectInput(ns("method"),"Transformation",choices=c("DESeq2 VST (recommended)"="vst","DESeq2 rlog (slower)"="rlog","log-CPM (exploratory)"="logcpm"),selected="vst"),
        shiny::numericInput(ns("qc_min_count"),"Minimum raw count",10,min=0,step=1),
        shiny::numericInput(ns("qc_min_samples"),"Required samples",2,min=1,step=1),
        shiny::selectInput(ns("top_genes"),"Most variable genes for PCA",choices=c("500"="500","1,000"="1000","2,000"="2000","All filtered genes"="Inf"),selected="500"),
        col_widths=c(3,3,3,3)
      ),
      bslib::layout_columns(
        shiny::selectInput(ns("x_pc"),"X axis",choices=paste0("PC",1:10),selected="PC1"),
        shiny::selectInput(ns("y_pc"),"Y axis",choices=paste0("PC",1:10),selected="PC2"),
        shiny::selectInput(ns("color"),"Color PCA by",choices=NULL),
        shiny::selectInput(ns("shape"),"Shape PCA by",choices=NULL),
        shiny::selectInput(ns("label"),"Label PCA by",choices=NULL),
        col_widths=c(2,2,3,3,2)
      ),
      shiny::actionButton(ns("apply_qc"),"Apply QC settings",class="btn-primary"),
      shiny::uiOutput(ns("qc_status"))
    ),
    bslib::layout_columns(
      bslib::card(bslib::card_header("Library size"),shiny::plotOutput(ns("library"))),
      bslib::card(bslib::card_header("Detected genes"),shiny::plotOutput(ns("detected"))),col_widths=c(6,6)
    ),
    bslib::navset_card_tab(
      bslib::nav_panel(
        "PCA",
        shiny::plotOutput(ns("pca"),height="620px"),
        shiny::downloadButton(ns("download_pca"),"Export PCA",class="btn-secondary mb-3")
      ),
      bslib::nav_panel("Variance explained",shiny::plotOutput(ns("scree"),height="500px")),
      bslib::nav_panel("Sample correlation",shiny::plotOutput(ns("correlation"),height="700px")),
      bslib::nav_panel("Sample distances",shiny::plotOutput(ns("distance"),height="700px")),
      bslib::nav_panel("Subset exploration",
        bslib::layout_columns(shiny::selectInput(ns("subset_var"),"Subset variable",choices=NULL),shiny::selectInput(ns("subset_level"),"Subset level",choices=NULL),col_widths=c(6,6)),
        bslib::navset_card_tab(
          bslib::nav_panel("Subset PCA",shiny::plotOutput(ns("subset_pca"),height="600px")),
          bslib::nav_panel("Subset correlation",shiny::plotOutput(ns("subset_correlation"),height="650px"))
        ))
    )
  )
}

mod_qc_server <- function(id,analysis_data) shiny::moduleServer(id,function(input,output,session){
  applied <- shiny::reactiveVal()
  default_settings <- function() list(method="vst",min_count=10,min_samples=min(2,ncol(analysis_data()$counts)),top_n=500)
  observeEvent(analysis_data(),{
    fields<-names(analysis_data()$metadata);applied(default_settings())
    shiny::updateNumericInput(session,"qc_min_samples",max=ncol(analysis_data()$counts),value=min(2,ncol(analysis_data()$counts)))
    shiny::updateSelectInput(session,"color",choices=c("None"="",fields),selected=if("treatment"%in%fields)"treatment" else "")
    shiny::updateSelectInput(session,"shape",choices=c("None"="",fields),selected="")
    shiny::updateSelectInput(session,"label",choices=c("No labels"="","Sample ID"=".sample_id",fields),selected=".sample_id")
    shiny::updateSelectInput(session,"subset_var",choices=fields,selected=if("treatment"%in%fields)"treatment" else fields[1])
    pcs<-paste0("PC",seq_len(min(10,nrow(analysis_data()$metadata)-1)));shiny::updateSelectInput(session,"x_pc",choices=pcs,selected="PC1");shiny::updateSelectInput(session,"y_pc",choices=pcs,selected="PC2")
  })
  observeEvent(input$apply_qc,{
    req(analysis_data(),input$method,input$qc_min_count,input$qc_min_samples,input$top_genes)
    top_n<-if(input$top_genes=="Inf") Inf else as.numeric(input$top_genes)
    candidate<-list(method=input$method,min_count=input$qc_min_count,min_samples=input$qc_min_samples,top_n=top_n)
    error<-tryCatch({gene_filter_mask(analysis_data()$counts,analysis_data()$annotations,character(),candidate$min_count,candidate$min_samples);NULL},error=function(e)e$message)
    if(!is.null(error)) shiny::showNotification(error,type="error") else applied(candidate)
  })
  observeEvent(input$subset_var,{
    req(analysis_data(),input$subset_var);values<-unique(as.character(analysis_data()$metadata[[input$subset_var]]));values<-values[!is.na(values)]
    shiny::updateSelectInput(session,"subset_level",choices=values,selected=if("NC"%in%values)"NC" else values[1])
  })

  qc <- shiny::reactive({
    req(analysis_data(),applied());s<-applied()
    shiny::withProgress(message="Preparing QC transformation",value=.3,{q<-prepare_qc_transformation(analysis_data()$counts,s$method,s$min_count,s$min_samples,s$top_n);shiny::incProgress(.7);q})
  })
  scores <- shiny::reactive(compute_matrix_pca(qc()$pca_matrix,analysis_data()$metadata))
  output$qc_status <- shiny::renderUI({
    q<-qc();method<-switch(q$method,vst="DESeq2 VST (blind)",rlog="DESeq2 rlog (blind)",logcpm="log-CPM")
    shiny::div(class="alert alert-info mt-3",sprintf("Applied to QC only: %s; count ≥ %s in at least %s sample(s); %s genes retained; %s most-variable genes used for PCA.",method,q$min_count,q$min_samples,format(q$genes_retained,big.mark=","),format(q$genes_used,big.mark=",")))
  })

  qc_frame<-shiny::reactive({req(analysis_data());data.frame(sample=colnames(analysis_data()$counts),library=library_sizes(analysis_data()$counts),detected=detected_genes(analysis_data()$counts))})
  output$library<-shiny::renderPlot({ggplot2::ggplot(qc_frame(),ggplot2::aes(sample,library))+ggplot2::geom_col(fill="#2c7fb8")+ggplot2::scale_y_continuous(labels=scales::label_comma())+ggplot2::labs(x=NULL,y="Raw count sum")+ggplot2::theme_minimal()+ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45,hjust=1))})
  output$detected<-shiny::renderPlot({ggplot2::ggplot(qc_frame(),ggplot2::aes(sample,detected))+ggplot2::geom_col(fill="#41ab5d")+ggplot2::labs(x=NULL,y="Genes with count > 0")+ggplot2::theme_minimal()+ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45,hjust=1))})

  pca_plot<-shiny::reactive({
    x<-scores();variance<-attr(x,"variance");xpc<-input$x_pc %||% "PC1";ypc<-input$y_pc %||% "PC2"
    shiny::validate(shiny::need(all(c(xpc,ypc)%in%names(x)),"Select valid principal components."))
    color_var<-input$color %||% "";shape_var<-input$shape %||% "";label_var<-input$label %||% ".sample_id"
    if(nzchar(shape_var)&&shape_var%in%names(x))x$.shape<-factor(as.character(x[[shape_var]]))
    x$.label<-if(identical(label_var,".sample_id"))x$sample_id else if(nzchar(label_var)&&label_var%in%names(x))as.character(x[[label_var]]) else ""
    has_color<-nzchar(color_var)&&color_var%in%names(x);has_shape<-".shape"%in%names(x)
    mapping<-if(has_color&&has_shape)ggplot2::aes(x=.data[[xpc]],y=.data[[ypc]],color=.data[[color_var]],shape=.shape) else if(has_color)ggplot2::aes(x=.data[[xpc]],y=.data[[ypc]],color=.data[[color_var]]) else if(has_shape)ggplot2::aes(x=.data[[xpc]],y=.data[[ypc]],shape=.shape) else ggplot2::aes(x=.data[[xpc]],y=.data[[ypc]])
    plot<-ggplot2::ggplot(x,mapping)+ggplot2::geom_point(size=3)
    if(any(nzchar(x$.label)))plot<-plot+ggplot2::geom_text(ggplot2::aes(label=.label),vjust=-1.15,check_overlap=TRUE,show.legend=FALSE)
    if(".shape"%in%names(x)){
      shape_values<-rep(c(16,17,15,18,3,8,7,9,10,12,13,14),length.out=nlevels(x$.shape))
      plot<-plot+ggplot2::scale_shape_manual(values=shape_values,name=shape_var)
    }
    plot+ggplot2::labs(x=sprintf("%s (%.1f%%)",xpc,variance[as.integer(sub("PC","",xpc))]),y=sprintf("%s (%.1f%%)",ypc,variance[as.integer(sub("PC","",ypc))]))+ggplot2::theme_minimal()+ggplot2::coord_cartesian(clip="off")+ggplot2::theme(plot.margin=ggplot2::margin(12,12,12,12))
  })
  output$pca<-shiny::renderPlot(pca_plot())
  output$download_pca<-shiny::downloadHandler(
    filename=function()"pca.png",
    contentType="image/png",
    content=function(file){
      grDevices::png(file,width=1800,height=1400,res=180)
      on.exit(grDevices::dev.off(),add=TRUE)
      print(pca_plot())
    }
  )
  output$scree<-shiny::renderPlot({v<-attr(scores(),"variance");n<-min(10,length(v));df<-data.frame(component=factor(paste0("PC",seq_len(n)),levels=paste0("PC",seq_len(n))),variance=v[seq_len(n)]);ggplot2::ggplot(df,ggplot2::aes(component,variance))+ggplot2::geom_col(fill="#756bb1")+ggplot2::geom_text(ggplot2::aes(label=sprintf("%.1f%%",variance)),vjust=-.3)+ggplot2::labs(x=NULL,y="Variance explained (%)")+ggplot2::theme_minimal()})
  output$correlation<-shiny::renderPlot({pheatmap::pheatmap(stats::cor(qc()$matrix),main=paste("Pearson correlation —",toupper(qc()$method)),border_color=NA)})
  output$distance<-shiny::renderPlot({pheatmap::pheatmap(as.matrix(stats::dist(t(qc()$matrix))),main=paste("Euclidean distance —",toupper(qc()$method)),border_color=NA)})

  subset_matrix <- shiny::reactive({
    req(input$subset_var,input$subset_level);keep<-as.character(analysis_data()$metadata[[input$subset_var]])==input$subset_level;shiny::validate(shiny::need(sum(keep)>=3,"At least three samples are required in the selected subset."))
    mat<-qc()$matrix[,keep,drop=FALSE];vars<-sort(apply(mat,1,stats::var),decreasing=TRUE);n<-min(qc()$genes_used,sum(vars>0));list(all=mat,pca=mat[names(vars)[seq_len(n)],,drop=FALSE],keep=keep)
  })
  output$subset_pca<-shiny::renderPlot({
    sm<-subset_matrix();x<-compute_matrix_pca(sm$pca,analysis_data()$metadata[sm$keep,,drop=FALSE]);variance<-attr(x,"variance");color_var<-if("source_group"%in%names(x))"source_group" else input$color;x$.label<-x$sample_id
    if(nzchar(color_var %||% "")&&color_var%in%names(x)){x$.subset_color<-factor(as.character(x[[color_var]]));plot<-ggplot2::ggplot(x,ggplot2::aes(PC1,PC2,color=.subset_color,label=.label))}else plot<-ggplot2::ggplot(x,ggplot2::aes(PC1,PC2,label=.label))
    plot+ggplot2::geom_point(size=4)+ggplot2::geom_text(nudge_y=.15)+ggplot2::labs(title=paste(input$subset_var,"=",input$subset_level),x=sprintf("PC1 (%.1f%%)",variance[1]),y=sprintf("PC2 (%.1f%%)",variance[2]),color=color_var)+ggplot2::theme_minimal()
  })
  output$subset_correlation<-shiny::renderPlot({sm<-subset_matrix();ann<-analysis_data()$metadata[sm$keep,intersect(c("source_group",input$subset_var),names(analysis_data()$metadata)),drop=FALSE];pheatmap::pheatmap(stats::cor(sm$all),annotation_col=ann,main=paste("Correlation:",input$subset_var,"=",input$subset_level),border_color=NA)})
})
