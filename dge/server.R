library(shiny)
library(limma)
library(ggplot2)
library(dplyr)
library(tibble)
library(stringr)
library(reshape2)
library(plotly)

# TODO: interactive plots and text
# TODO: normalise RNA-seq data better
# TODO: rearrange gene columns (symbol and entrez) order
# 
# Lots of repeated code... should refactor sometime


# For dev/debugging
#example_csv_filename <- "example_4.csv"

############################################################
## Load and arrange data
############################################################

load("data/data.RData")
original_rna_phenotype <- rna_phenotype %>% 
  select(sample_id, study:afs_score) #%>%
  # filter(! is.na(age))
original_array_phenotype <- array_phenotype

############################################################
## Functions
############################################################

# get_slider_samples <- function(phenotype, slider) {
#   phenotype %>% filter(between(bmi, slider[1], slider[2])) %>% .$sample_id
# }
# 
# get_checkbox_samples <- function(phenotype, checkbox) {
#   phenotype %>% filter(class %in% checkbox) %>% .$sample_id
# }

############################################################
## Server
############################################################


shinyServer(function(input, output) {
  
  # Set reactive values
  rv <- reactiveValues(
    exprs = NA,
    rna_phenotype = original_rna_phenotype,
    array_phenotype = original_array_phenotype,
    selection_method = NULL,
    rna_selected_gene = NA,
    array_selected_gene = NA,
    rna_results = NULL,
    array_results = NULL,
    rna_plot_list = NULL,
    array_plot_list = NULL
  )
  
  ############################################################
  ## PANEL 1
  ############################################################
  
  observeEvent({input$file1}, {
    message("CSV file uploaded")
    tryCatch(
      {
        df <- read.csv(input$file1$datapath,
                       header=FALSE,
                       sep=",",
                       stringsAsFactors=FALSE)
        colnames(df) <- c("sample_id", "group")
        print(class(df$group))
        if (is.numeric(df$group) & length(unique(df$group)) > 2) {
          message("Detected numerical covariate")
          colnames(df) <- c("sample_id", "variable")
        } else {
          message("Detected categorical covariate")
          df$group <- factor(df$group)
          stopifnot(length(levels(df$group)) == 2)
        }
      },
      error = function(e) {
        # return a safeError if a parsing error occurs
        stop(safeError(e))
      }
    )
    rv$rna_phenotype <- merge(df, original_rna_phenotype, by="sample_id")
    rv$array_phenotype <- merge(df, original_array_phenotype, by="sample_id")
  })
  
  output$rna_phenotype <- renderDataTable({
    datatable(rv$rna_phenotype,
              options=list(pageLength=10), 
              selection="single")
  })
  
  output$rna_sample_list_text <- renderText({
    if ("group" %in% colnames(rv$rna_phenotype)) {
      group_names <- levels(rv$rna_phenotype$group)
      txt <- sprintf("DGE comparison: %s (n=%d) vs %s (n=%d)",
                     group_names[1], sum(rv$rna_phenotype$group == group_names[1]),
                     group_names[2], sum(rv$rna_phenotype$group == group_names[2]))
    } else if ("variable" %in% colnames(rv$rna_phenotype)) {
      txt <- sprintf("Numerical factor (n=%d)", nrow(rv$rna_phenotype))
    } else {
      txt <- paste0("Upload a CSV file in the left panel to define groups for a ",
                    "differential expression analysis. Available RNA-seq samples ",
                    "are shown below.")
    }
    return(txt)
  })
  
  output$array_phenotype <- renderDataTable({
    datatable(rv$array_phenotype,
              options=list(pageLength=10), 
              selection="single")
  })
  
  output$array_sample_list_text <- renderText({
    if ("group" %in% colnames(rv$array_phenotype)) {
      group_names <- levels(rv$array_phenotype$group)
      txt <- sprintf("DGE comparison: %s (n=%d) vs %s (n=%d)",
                     group_names[1], sum(rv$array_phenotype$group == group_names[1]),
                     group_names[2], sum(rv$array_phenotype$group == group_names[2]))
    } else if ("variable" %in% colnames(rv$array_phenotype)) {
      txt <- sprintf("Numerical factor (n=%d)", nrow(rv$array_phenotype))
    } else {
      txt <- paste0("Upload a CSV file in the left panel to define groups for a ",
                    "differential expression analysis. Available microarray samples ",
                    "are shown below.")
    }
    return(txt)
  })
  
  
  ############################################################
  ## PANEL 2
  ############################################################
  
  # Perform DGE analysis and get top table
  observeEvent({rv$rna_phenotype}, {
    message("Performing RNA-seq DGE analysis...")
    
    # Check there's at least 3 samples in each group
    if (! "variable" %in% colnames(rv$rna_phenotype)) {
      group_names <- levels(rv$rna_phenotype$group)
      if (sum(rv$rna_phenotype$group == group_names[1]) < 3 |
          sum(rv$rna_phenotype$group == group_names[2]) < 3) {
        message("Too few samples in groups.")
        rv$rna_results <- list(top_table=data.frame())
        rv$rna_selected_gene <- NULL
        return()
      }
    }
    if (! "group" %in% colnames(rv$rna_phenotype)) {
      if (nrow(rv$rna_phenotype) < 3) {
        message("Too few samples in groups.")
        rv$rna_results <- list(top_table=data.frame())
        rv$rna_selected_gene <- NULL
        return()
      }
    }
    
    # Setup design
    # TODO: add age option
    if ("group" %in% colnames(rv$rna_phenotype)) {
      design <- model.matrix(~group, data=rv$rna_phenotype)
    } else {
      design <- model.matrix(~variable, data=rv$rna_phenotype)
    }
    
    #design <- model.matrix(~group+age, data=rv$rna_phenotype)
    
    # Fit model
    lm <- lmFit(rna_cc[,rv$rna_phenotype$sample_id], design)
    # lm <- lmFit(rna_bc[,rv$rna_phenotype$sample_id], design) # Sanity check
    fit <- eBayes(lm)
    rna_top_table <- topTable(fit, coef=2, n=Inf) %>%
      tibble::rownames_to_column(var="ensembl_id") %>%
      select(-t, -B) %>%
      merge(rna_gene_info, by="ensembl_id") %>%
      arrange(P.Value)
    
    rv$rna_results <- list(top_table=rna_top_table)
    rv$rna_selected_gene <- rna_top_table[1,"ensembl_id"]
    
  })

  # Output top table
  output$rna_top_table <- renderDataTable({
    datatable(rv$rna_results$top_table %>% head(input$n_top_table) %>% 
                format(digits=3),
              options=list(pageLength=5), selection="single")
  })

  # Update selected gene when row is clicked
  observeEvent({input$rna_top_table_cell_clicked}, {
    if (!is.null(input$rna_top_table_cell_clicked$row)) {
      message(paste0("Row clicked: ", input$rna_top_table_cell_clicked$row))
      rv$rna_selected_gene <- rv$rna_results$top_table[input$rna_top_table_cell_clicked$row,"ensembl_id"]
    }
  })

  # Get data for plotting expression
  observeEvent({rv$rna_results; rv$rna_selected_gene}, {
    message("Updating plot data: ", rv$rna_selected_gene)

    # Get gene name
    gene_name <- rna_gene_info %>%
      filter(ensembl_id == rv$rna_selected_gene) %>%
      .$symbol

    # Get expression and phenotype data
    dat <- rv$rna_phenotype %>%
      mutate(exprs=rna_cc[rv$rna_selected_gene,sample_id])
    
    # Get expression for non-cycle-corrected data
    cycle_dat <- original_rna_phenotype %>%
      mutate(exprs=rna_bc[rv$rna_selected_gene,sample_id]) %>%
      merge(rv$rna_phenotype %>% select(1,2), all=TRUE)

    rv$rna_plot_list <- list(dat=dat, cycle_dat=cycle_dat, 
                         gene=rv$rna_selected_gene, gene_name=gene_name)
  })

  # Plot expression with ggplot
  output$ggplot_rna_exprs <- renderPlot({
    if (! "group" %in% colnames(rv$rna_plot_list$dat) &
        ! "variable" %in% colnames(rv$rna_plot_list$dat)) {
      return()
    }
    g <- ggplot(rv$rna_plot_list$dat, 
                aes_string(x=colnames(rv$rna_phenotype)[2], y="exprs", 
                           color=colnames(rv$rna_phenotype)[2])) +
      labs(title=paste(rv$rna_plot_list$gene, " - ", rv$rna_plot_list$gene_name)) + 
      theme_bw()
    if ("group" %in% colnames(rv$rna_plot_list$dat)) {
      # Add jitter & boxplot if categorical groups
      g <- g + 
        geom_jitter(height=0, width=0.3, size=2) +
        geom_boxplot(alpha=0.1, outlier.alpha=0)
    } else {
      # Scatterplot & linear fit if numerical
      g <- g + 
        geom_point(size=2) +
        geom_smooth(method="lm", alpha=0.2)
    }
    return(g)
  })
  
  # Plot cycle expression with ggplot
  output$ggplot_rna_cycle <- renderPlot({
    if (! "group" %in% colnames(rv$rna_plot_list$dat) &
        ! "variable" %in% colnames(rv$rna_plot_list$dat)) {
      return()
    }
    ggplot(rv$rna_plot_list$cycle_dat, 
             aes_string(x="cycle_stage", y="exprs", 
                        color=colnames(rv$rna_phenotype)[2])) +
      geom_jitter(height=0, width=0.3, size=2) +
      labs(title=paste(rv$rna_plot_list$gene, " - ", rv$rna_plot_list$gene_name)) + 
      theme_bw()
  })
  
  # P-value histogram
  output$ggplot_rna_pval <- renderPlot({
    if (nrow(rv$rna_results$top_table) == 0) {
      return()
    }
    print(head(rv$rna_results$top_table))
    ggplot(rv$rna_results$top_table, aes(x=P.Value)) +
      geom_histogram(binwidth=0.05, fill="white", color="black") +
      labs(title="P-value histogram") + 
      theme_bw()
  })
  
  ############################################################
  ## PANEL 3
  ############################################################
  
  # Perform DGE analysis and get top table
  observeEvent({rv$array_phenotype}, {
    message("Performing microarray DGE analysis...")
    
    # Check there's at least 3 samples in each group
    if (! "variable" %in% colnames(rv$array_phenotype)) {
      group_names <- levels(rv$array_phenotype$group)
      if (sum(rv$array_phenotype$group == group_names[1]) < 3 |
          sum(rv$array_phenotype$group == group_names[2]) < 3) {
        message("Too few samples in groups.")
        rv$array_results <- list(top_table=data.frame())
        rv$array_selected_gene <- NULL
        return()
      }
    }
    if (! "group" %in% colnames(rv$array_phenotype)) {
      if (nrow(rv$array_phenotype) < 3) {
        message("Too few samples in groups.")
        rv$array_results <- list(top_table=data.frame())
        rv$array_selected_gene <- NULL
        return()
      }
    }
    
    # Setup design
    # TODO: add age option
    if ("group" %in% colnames(rv$array_phenotype)) {
      design <- model.matrix(~group, data=rv$array_phenotype)
    } else {
      design <- model.matrix(~variable, data=rv$array_phenotype)
    }
    
    #design <- model.matrix(~group+age, data=rv$array_phenotype)
    
    # Fit model
    lm <- lmFit(array_cc[,rv$array_phenotype$sample_id], design)
    # lm <- lmFit(array_bc[,rv$array_phenotype$sample_id], design) # Sanity check
    fit <- eBayes(lm)
    array_top_table <- topTable(fit, coef=2, n=Inf) %>%
      tibble::rownames_to_column(var="illumina_id") %>%
      select(-t, -B) %>%
      merge(array_probe_info, by="illumina_id") %>%
      arrange(P.Value)
    
    rv$array_results <- list(top_table=array_top_table)
    rv$array_selected_gene <- array_top_table[1,"illumina_id"]
    
  })
  
  # Output top table
  output$array_top_table <- renderDataTable({
    datatable(rv$array_results$top_table %>% head(input$n_top_table) %>% 
                format(digits=3),
              options=list(pageLength=5), selection="single")
  })
  
  # Update selected gene when row is clicked
  observeEvent({input$array_top_table_cell_clicked}, {
    if (!is.null(input$array_top_table_cell_clicked$row)) {
      message(paste0("Row clicked: ", input$array_top_table_cell_clicked$row))
      rv$array_selected_gene <- rv$array_results$top_table[input$array_top_table_cell_clicked$row,"illumina_id"]
    }
  })
  
  # Get data for plotting expression
  observeEvent({rv$array_results; rv$array_selected_gene}, {
    message("Updating plot data: ", rv$array_selected_gene)
    
    # Get gene name
    gene_name <- array_probe_info %>%
      filter(illumina_id == rv$array_selected_gene) %>%
      .$symbol
    
    # Get expression and phenotype data
    dat <- rv$array_phenotype %>%
      mutate(exprs=array_cc[rv$array_selected_gene,sample_id])
    
    # Get expression for non-cycle-corrected data
    cycle_dat <- original_array_phenotype %>%
      mutate(exprs=array_bc[rv$array_selected_gene,sample_id]) %>%
      merge(rv$array_phenotype %>% select(1,2), all=TRUE)
    
    rv$array_plot_list <- list(dat=dat, cycle_dat=cycle_dat, 
                             gene=rv$array_selected_gene, gene_name=gene_name)
  })
  
  # Plot expression with ggplot
  output$ggplot_array_exprs <- renderPlot({
    if (! "group" %in% colnames(rv$array_plot_list$dat) &
        ! "variable" %in% colnames(rv$array_plot_list$dat)) {
      return()
    }
    g <- ggplot(rv$array_plot_list$dat, 
                aes_string(x=colnames(rv$array_phenotype)[2], y="exprs", 
                           color=colnames(rv$array_phenotype)[2])) +
      labs(title=paste(rv$array_plot_list$gene, " - ", rv$array_plot_list$gene_name)) + 
      theme_bw()
    if ("group" %in% colnames(rv$array_plot_list$dat)) {
      # Add jitter & boxplot if categorical groups
      g <- g + 
        geom_jitter(height=0, width=0.3, size=2) +
        geom_boxplot(alpha=0.1, outlier.alpha=0)
    } else {
      # Scatterplot & linear fit if numerical
      g <- g + 
        geom_point(size=2) +
        geom_smooth(method="lm", alpha=0.2)
    }
    return(g)
  })
  
  # Plot cycle expression with ggplot
  output$ggplot_array_cycle <- renderPlot({
    if (! "group" %in% colnames(rv$array_plot_list$dat) &
        ! "variable" %in% colnames(rv$array_plot_list$dat)) {
      return()
    }
    ggplot(rv$array_plot_list$cycle_dat, 
           aes_string(x="model_day", y="exprs", 
                      color=colnames(rv$array_phenotype)[2])) +
      geom_jitter(height=0, width=0.3, size=2) +
      labs(title=paste(rv$array_plot_list$gene, " - ", rv$array_plot_list$gene_name)) + 
      theme_bw()
  })
  
  # P-value histogram
  output$ggplot_array_pval <- renderPlot({
    if (nrow(rv$array_results$top_table) == 0) {
      return()
    }
    print(head(rv$array_results$top_table))
    ggplot(rv$array_results$top_table, aes(x=P.Value)) +
      geom_histogram(binwidth=0.05, fill="white", color="black") +
      labs(title="P-value histogram") + 
      theme_bw()
  })
  
  
  
  
  
})