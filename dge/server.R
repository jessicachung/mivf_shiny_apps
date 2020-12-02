library(shiny)
library(limma)
library(ggplot2)
library(dplyr)
library(tibble)
library(stringr)
library(reshape2)
library(plotly)

############################################################
## Load and arrange data
############################################################

load("data/data.RData")

stopifnot(colnames(array_bc) == array_phenotype$sample_id)
stopifnot(colnames(array_cc) == array_phenotype$sample_id)

original_rna_phenotype <- rna_phenotype %>% 
  mutate(model_time=round(model_time, 2))
original_array_phenotype <- array_phenotype %>% 
  mutate(model_time=round(model_time, 2))

############################################################
## Functions
############################################################

# Common functions for both RNA and microarray analysis

ok_sample_numbers <- function(pheno) {
  # Check there's at least 3 samples in each group
  if (! "variable" %in% colnames(pheno)) {
    group_names <- levels(pheno$group)
    if (sum(pheno$group == group_names[1]) < 3 |
        sum(pheno$group == group_names[2]) < 3) {
      return(FALSE)
    }
  }
  if (! "group" %in% colnames(pheno)) {
    if (nrow(pheno) < 3) {
      return(FALSE)
    }
  }
  return(TRUE)
}

perform_dge <- function(exprs, pheno, gene_info, gene_column) {
  if (! ok_sample_numbers(pheno)) {
    message("Too few samples in groups.")
    return(NULL)
  }
  
  # Setup design
  if ("group" %in% colnames(pheno)) {
    design <- model.matrix(~group, data=pheno)
  } else {
    design <- model.matrix(~variable, data=pheno)
  }
  # TODO: add age option
  # design <- model.matrix(~group+age, data=pheno)
  
  # Fit model
  lm <- lmFit(exprs[,pheno$sample_id], design)
  fit <- eBayes(lm)
  top_table <- topTable(fit, coef=2, n=Inf) %>%
    tibble::rownames_to_column(var=gene_column) %>%
    select(-t, -B) %>%
    merge(gene_info, by=gene_column) %>%
    arrange(P.Value)
  return(top_table)
}

get_plot_data <- function(gene_id, gene_info, bc_exprs, cc_exprs, phenotype) {
  # Get gene name
  gene_name <- gene_info[gene_info[,1] == gene_id, "symbol"]
  
  # Get expression and phenotype data
  dat <- phenotype %>%
    mutate(exprs=cc_exprs[gene_id,sample_id])
  
  # Get expression for non-cycle-corrected data
  cycle_dat <- phenotype %>%
    mutate(exprs=bc_exprs[gene_id,sample_id],
           cc_exprs=cc_exprs[gene_id,sample_id]) %>%
    merge(phenotype %>% select(1,2), all=TRUE)
  
  return(list(dat=dat, cycle_dat=cycle_dat, 
              gene=gene_id, gene_name=gene_name))
}

expression_plot <- function(plot_list, phenotype, plot_type="ggplot") {
  if (plot_type == "ggplot") {
    size <- 2
  } else {
    size <- 1
  }
  if (! "group" %in% colnames(plot_list$dat) &
      ! "variable" %in% colnames(plot_list$dat)) {
    return()
  }
  g <- ggplot(plot_list$dat, 
              aes_string(x=colnames(phenotype)[2], y="exprs", 
                         color=colnames(phenotype)[2])) +
    labs(title=paste(plot_list$gene, " - ", plot_list$gene_name)) + 
    theme_bw()
  if ("group" %in% colnames(plot_list$dat)) {
    # Add jitter & boxplot if categorical groups
    g <- g + 
      geom_jitter(aes(text=sample_id), height=0, width=0.3, size=size) +
      geom_boxplot(alpha=0.1, outlier.alpha=0)
  } else {
    # Scatterplot & linear fit if numerical
    g <- g + 
      geom_point(aes(text=sample_id), size=size) +
      geom_smooth(method="lm", alpha=0.2)
  }
  
  if (plot_type == "plotly" & "group" %in% colnames(plot_list$dat)) {
    # Need to manually edit plotly points to make outlier transparent
    # (although points are still able to be hovered over)
    g <- plotly_build(g)
    for(i in 1:length(g$x$data)) {
      if (is.na(g$x$data[[i]]$marker$opacity)) {
        g$x$data[[i]]$marker$opacity <- 0
      }
    }
  }
  return(g)
}

cycle_plot <- function(plot_list, phenotype, x="model_time", y="exprs",
                       plot_type="ggplot") {
  if (plot_type == "ggplot") {
    size <- 2
  } else {
    size <- 1
  }
  if (! "group" %in% colnames(plot_list$dat) &
      ! "variable" %in% colnames(plot_list$dat)) {
    return()
  }
  ggplot(plot_list$cycle_dat, 
         aes_string(x=x, y=y, 
                    color=colnames(phenotype)[2],
                    text="sample_id")) +
    geom_point(size=size) +
    labs(title=paste(plot_list$gene, " - ", plot_list$gene_name)) + 
    theme_bw()
}

histogram_plot <- function(top_table) {
  if (nrow(top_table) == 0) {
    return()
  }
  ggplot(top_table, aes(x=P.Value)) +
    geom_histogram(binwidth=0.05, boundary=0, fill="white", color="black") +
    labs(title="P-value histogram") + 
    theme_bw()
}


############################################################
## Server
############################################################

shinyServer(function(input, output) {
  
  # Set reactive values
  rv <- reactiveValues(
    exprs = NA,
    rna_phenotype = original_rna_phenotype,
    array_phenotype = original_array_phenotype,
    filtered_rna_phenotype = original_rna_phenotype,
    filtered_array_phenotype = original_array_phenotype,
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
  
  observeEvent({input$file1; input$csv_select}, {
    tryCatch(
      {
        if (input$csv_select == "age") {
          df <- rbind(
            original_rna_phenotype %>% select(sample_id, age),
            original_array_phenotype %>% select(sample_id, age)
          ) %>%
            filter(! is.na(age)) %>%
            mutate(age=as.numeric(age)) %>%
            .[! duplicated(.),]
        } else if (input$csv_select == "endo") {
          df <- rbind(
            original_rna_phenotype %>% select(sample_id, endo),
            original_array_phenotype %>% select(sample_id, endo)
          ) %>%
            filter(! is.na(endo)) %>%
            .[! duplicated(.),]
        } else {
          if (! is.null(input$file1$datapath)) {
            df <- read.csv(input$file1$datapath,
                           header=FALSE,
                           sep=",",
                           stringsAsFactors=FALSE)
            message("CSV file uploaded")
          } else {
            return(0)
          }
        }
        colnames(df) <- c("sample_id", "group")
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
    rv$filtered_rna_phenotype <- rv$rna_phenotype
    rv$array_phenotype <- merge(df, original_array_phenotype, by="sample_id")
    rv$filtered_array_phenotype <- rv$array_phenotype
  })
  
  output$rna_phenotype <- renderDataTable({
    datatable(rv$filtered_rna_phenotype,
              options=list(pageLength=10), 
              selection="single")
  })
  
  output$rna_sample_list_text <- renderText({
    if ("group" %in% colnames(rv$filtered_rna_phenotype)) {
      group_names <- levels(rv$filtered_rna_phenotype$group)
      txt <- sprintf("DGE comparison: %s (n=%d) vs %s (n=%d)",
                     group_names[1], sum(rv$filtered_rna_phenotype$group == group_names[1]),
                     group_names[2], sum(rv$filtered_rna_phenotype$group == group_names[2]))
    } else if ("variable" %in% colnames(rv$filtered_rna_phenotype)) {
      txt <- sprintf("Numerical factor (n=%d)", nrow(rv$filtered_rna_phenotype))
    } else {
      txt <- paste0("Upload a CSV file in the left panel to define groups for a ",
                    "differential expression analysis. Available RNA-seq samples ",
                    "are shown below.")
    }
    return(txt)
  })
  
  output$array_phenotype <- renderDataTable({
    datatable(rv$filtered_array_phenotype,
              options=list(pageLength=10), 
              selection="single")
  })
  
  output$array_sample_list_text <- renderText({
    if ("group" %in% colnames(rv$filtered_array_phenotype)) {
      group_names <- levels(rv$filtered_array_phenotype$group)
      txt <- sprintf("DGE comparison: %s (n=%d) vs %s (n=%d)",
                     group_names[1], 
                     sum(rv$filtered_array_phenotype$group == group_names[1]),
                     group_names[2], 
                     sum(rv$filtered_array_phenotype$group == group_names[2]))
    } else if ("variable" %in% colnames(rv$filtered_array_phenotype)) {
      txt <- sprintf("Numerical factor (n=%d)", nrow(rv$filtered_array_phenotype))
    } else {
      txt <- paste0("Upload a CSV file in the left panel to define groups for a ",
                    "differential expression analysis. Available microarray samples ",
                    "are shown below.")
    }
    return(txt)
  })
  
  
  # Update samples when filtered on cycle stage or study
  observeEvent({rv$array_phenotype; rv$rna_phenotype; input$cycle_checkbox; 
    input$study_checkbox}, {
      message("Updating phenotype dataframe...")
      if (input$cycle_select == "all") {
        ok_cycles <- c("s1", "s2", "s3", "s4", "s5", "s6", "s7")
      } else {
        ok_cycles <- input$cycle_checkbox
      }
      if (input$study_select == "all") {
        ok_studies <- c("study_1", "study_2", "HMB")
      } else {
        ok_studies <- input$study_checkbox
      }
      
      rv$filtered_rna_phenotype <- rv$rna_phenotype %>%
        filter(study %in% ok_studies)
      rv$filtered_array_phenotype <- rv$array_phenotype %>%
        filter(study %in% ok_studies)
  })
  
  ############################################################
  ## PANEL 2
  ############################################################
  
  # Perform DGE analysis and get top table
  observeEvent({rv$filtered_rna_phenotype}, {
    message("Performing RNA-seq DGE analysis...")
    top_table <- perform_dge(exprs=rna_cc, pheno=rv$filtered_rna_phenotype,
                             gene_info=rna_gene_info, gene_column="ensembl_id")
    if (! is.null(top_table)) {
      rv$rna_results <- list(top_table=top_table)
      rv$rna_selected_gene <- top_table[1,"ensembl_id"]
    } else {
      rv$rna_results <- list(top_table=data.frame())
      rv$rna_selected_gene <- NULL
    }
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
      rv$rna_selected_gene <- rv$rna_results$top_table[
        input$rna_top_table_cell_clicked$row,"ensembl_id"]
    }
  })

  # Get data for plotting expression
  observeEvent({rv$rna_results; rv$rna_selected_gene}, {
    message("Updating plot data: ", rv$rna_selected_gene)
    rv$rna_plot_list <- get_plot_data(gene_id=rv$rna_selected_gene, 
                                      gene_info=rna_gene_info,
                                      bc_exprs=rna_bc, cc_exprs=rna_cc, 
                                      phenotype=rv$filtered_rna_phenotype)
    
  })

  # Plot expression data
  output$ggplot_rna_exprs <- renderPlot({
    expression_plot(plot_list=rv$rna_plot_list,
                    phenotype=rv$filtered_rna_phenotype,
                    plot_type="ggplot")
  })
  
  output$plotly_rna_exprs <- renderPlotly({
    expression_plot(plot_list=rv$rna_plot_list,
                    phenotype=rv$filtered_rna_phenotype,
                    plot_type="plotly")
  })
  
  # Plot cycle expression
  output$ggplot_rna_cycle <- renderPlot({
    cycle_plot(plot_list=rv$rna_plot_list, phenotype=rv$filtered_rna_phenotype,
               x="model_time", plot_type="ggplot")
  })
  
  output$plotly_rna_cycle <- renderPlotly({
    cycle_plot(plot_list=rv$rna_plot_list, phenotype=rv$filtered_rna_phenotype,
               x="model_time", plot_type="plotly")
  })
  
  # Plot cycle-corrected expression
  output$ggplot_rna_corrected <- renderPlot({
    cycle_plot(plot_list=rv$rna_plot_list, phenotype=rv$filtered_rna_phenotype,
               x="model_time", y="cc_exprs", plot_type="ggplot")
  })
  
  output$plotly_rna_corrected <- renderPlotly({
    cycle_plot(plot_list=rv$rna_plot_list, phenotype=rv$filtered_rna_phenotype,
               x="model_time", y="cc_exprs", plot_type="plotly")
  })
  
  # Plot P-value histogram
  output$ggplot_rna_pval <- renderPlot({
    histogram_plot(rv$rna_results$top_table)
  })
  
  ############################################################
  ## PANEL 3
  ############################################################
  
  # Perform DGE analysis and get top table
  observeEvent({rv$filtered_array_phenotype}, {
    message("Performing microarray DGE analysis...")
    top_table <- perform_dge(exprs=array_cc, pheno=rv$filtered_array_phenotype,
                             gene_info=array_probe_info, gene_column="illumina_id")
    if (! is.null(top_table)) {
      rv$array_results <- list(top_table=top_table)
      rv$array_selected_gene <- top_table[1,"illumina_id"]
    } else {
      rv$array_results <- list(top_table=data.frame())
      rv$array_selected_gene <- NULL
    }
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
      rv$array_selected_gene <- rv$array_results$top_table[
        input$array_top_table_cell_clicked$row,"illumina_id"]
    }
  })
  
  # Get data for plotting expression
  observeEvent({rv$array_results; rv$array_selected_gene}, {
    message("Updating plot data: ", rv$array_selected_gene)
    rv$array_plot_list <- get_plot_data(gene_id=rv$array_selected_gene, 
                                        gene_info=array_probe_info,
                                        bc_exprs=array_bc, cc_exprs=array_cc, 
                                        phenotype=rv$filtered_array_phenotype)
    
  })
  
  # Plot expression data
  output$ggplot_array_exprs <- renderPlot({
    expression_plot(plot_list=rv$array_plot_list,
                    phenotype=rv$filtered_array_phenotype,
                    plot_type="ggplot")
  })
  
  output$plotly_array_exprs <- renderPlotly({
    expression_plot(plot_list=rv$array_plot_list,
                    phenotype=rv$filtered_array_phenotype,
                    plot_type="plotly")
  })
  
  # Plot cycle expression
  output$ggplot_array_cycle <- renderPlot({
    cycle_plot(plot_list=rv$array_plot_list, phenotype=rv$filtered_array_phenotype,
               x="model_time", plot_type="ggplot")
  })
  
  output$plotly_array_cycle <- renderPlotly({
    cycle_plot(plot_list=rv$array_plot_list, phenotype=rv$filtered_array_phenotype,
               x="model_time", plot_type="plotly")
  })
  
  # Plot cycle-corrected expression
  output$ggplot_array_corrected <- renderPlot({
    cycle_plot(plot_list=rv$array_plot_list, phenotype=rv$filtered_array_phenotype,
               x="model_time", y="cc_exprs", plot_type="ggplot")
  })
  
  output$plotly_array_corrected <- renderPlotly({
    cycle_plot(plot_list=rv$array_plot_list, phenotype=rv$filtered_array_phenotype,
               x="model_time", y="cc_exprs", plot_type="plotly")
  })
  
  
  # Plot P-value histogram
  output$ggplot_array_pval <- renderPlot({
    histogram_plot(rv$array_results$top_table)
  })
  
})