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

load("data/bmi_rnaseq.RData")

# Dimensions
message(paste0("Expression matrix dimensions: ", 
               nrow(cc_exprs)), " x ", ncol(cc_exprs))

# # Subset expression matrix with prior genes
# subset_exprs <- cc_exprs[prior_genes, phenotype$sample_id]
# message(paste0("Expression matrix dimensions with prior genes: ", 
#                nrow(subset_exprs)), " x ", ncol(subset_exprs))

############################################################
## Functions
############################################################

get_slider_samples <- function(phenotype, slider) {
  phenotype %>% filter(between(bmi, slider[1], slider[2])) %>% .$sample_id
}

get_checkbox_samples <- function(phenotype, checkbox) {
  phenotype %>% filter(class %in% checkbox) %>% .$sample_id
}

############################################################
## Server
############################################################


shinyServer(function(input, output) {
  
  # Set reactive values
  rv <- reactiveValues(
    exprs = cc_exprs,
    phenotype = NULL,
    selection_method = NULL,
    group_1_samples = NULL,
    group_2_samples = NULL,
    selected_gene = NA,
    results = NULL,
    plot_list = NULL
  )
  
  # # Update exprs matrix when selected
  # observeEvent({input$exprs_selection}, {
  #   message("Updating exprs matrix...")
  #   if (input$exprs_selection == "good") {
  #     rv$exprs <- good_exprs
  #   } else if (input$exprs_selection == "prior") {
  #     rv$exprs <- subset_exprs
  #   } else if (input$exprs_selection == "high") {
  #     rv$exprs <- high_exprs
  #   }
  # })
  
  observeEvent({input$patient_selection}, {
    message("Updating phenotype dataframe...")
    if (input$patient_selection == "all") {
      rv$phenotype <- phenotype
    } else if (input$patient_selection == "endo") {
      rv$phenotype <- phenotype %>% filter(endo == 1)
    } else if (input$patient_selection == "control") {
      rv$phenotype <- phenotype %>% filter(endo == 0)
    } else if (input$patient_selection == "woi") {
      rv$phenotype <- phenotype %>% 
        mutate(cycle_stage=as.numeric(cycle_stage)) %>% 
        filter(between(cycle_stage, 4, 6))
    } else if(input$patient_selection == "afs1-5") {
      rv$phenotype <- phenotype %>% 
        mutate(afs_score=as.numeric(afs_score)) %>% 
        filter(between(afs_score, 1, 5))
    }
  })
  
  # Update rv groups when ranges/checkboxes change
  observeEvent({input$selection; input$slider_1; input$slider_2;
    input$checkbox_1; input$checkbox_2; rv$phenotype}, {
      message("Updating groups...")
      if (input$selection == "slider") {
        rv$group_1_samples <- get_slider_samples(rv$phenotype, input$slider_1)
        rv$group_2_samples <- get_slider_samples(rv$phenotype, input$slider_2)
      } else {
        rv$group_1_samples <- get_checkbox_samples(rv$phenotype, input$checkbox_1)
        rv$group_2_samples <- get_checkbox_samples(rv$phenotype, input$checkbox_2)
      }
    }
  )
  
  # Update text with number of patients in each group
  output$n_group_1 <- renderText({
    paste("Group 1 samples:", length(rv$group_1_samples))
  })
  output$n_group_2 <- renderText({
    paste("Group 2 samples:", length(rv$group_2_samples))
  })
  
  # Perform DGE analysis and get top table
  observeEvent({rv$exprs; rv$phenotype; rv$group_1_samples; rv$group_2_samples; 
    input$n_top_table}, {
      message("Performing DGE analysis...")
      
      # Subset data
      pheno <- rv$phenotype %>% 
        filter(sample_id %in% c(rv$group_1_samples, rv$group_2_samples)) %>%
        mutate(group=ifelse(sample_id %in% rv$group_1_samples, "Group 1", "Group 2"))
      
      # Check there's at least 3 samples in each group
      if (sum(pheno$group == "Group 1") < 3 | 
          sum(pheno$group == "Group 2") < 3) {
        message("Too few samples in groups.")
        rv$results <- NULL
        rv$selected_gene <- NULL
        
      } else {
        # Setup design
        design <- model.matrix(~group, data = pheno)
        
        # Fit model
        lm <- lmFit(rv$exprs[,pheno$sample_id], design)
        fit <- eBayes(lm)
        top_table <- topTable(fit, coef=2, n=input$n_top_table) %>% 
          tibble::rownames_to_column(var="ensembl_id") %>%
          select(-t, -B) %>%
          merge(qld_counts$gene_info, by="ensembl_id") %>% 
          arrange(P.Value)
        
        rv$results <- list(top_table=top_table, pheno=pheno)
        rv$selected_gene <- top_table[1,"ensembl_id"]
      }
    })
  
  # Output top table
  output$top_table <- renderDataTable({
    datatable(rv$results$top_table %>% format(digits=3), 
              options=list(pageLength=5), selection="single")
  })
  
  # Update selected gene when row is clicked
  observeEvent({input$top_table_cell_clicked}, {
    if (!is.null(input$top_table_cell_clicked$row)) {
      message(paste0("Row clicked: ", input$top_table_cell_clicked$row))
      rv$selected_gene <- rv$results$top_table[input$top_table_cell_clicked$row,"ensembl_id"]
    }
  })
  
  # Get data for plotting expression
  observeEvent({rv$results; rv$selected_gene}, {
    message("Updating plot data: ", rv$selected_gene)
    
    # Get gene name
    gene_name <- qld_counts$gene_info %>% 
      filter(ensembl_id == rv$selected_gene) %>%
      .$symbol

    # Get expression and phenotype data
    dat <- data.frame(sample_id=colnames(rv$exprs),
                      exprs=rv$exprs[rv$selected_gene,]) %>%
      merge(rv$phenotype %>% select(sample_id, bmi, text), by="sample_id", all.x=TRUE) %>%
      merge(rv$results$pheno %>% select(sample_id, group), by="sample_id", all.x=TRUE)

    rv$plot_list <- list(dat=dat, gene=rv$selected_gene, gene_name=gene_name)
  })

  # Plot BMI vs expression with ggplot
  output$ggplot <- renderPlot({
    ggplot(rv$plot_list$dat, aes(x=bmi, y=exprs, color=group)) +
      geom_point(size=2) +
      labs(title=paste(rv$plot_list$gene, " - ", rv$plot_list$gene_name))
  })

  # Plot BMI vs expression with plotly
  output$plotly <- renderPlotly({
    ggplot(rv$plot_list$dat, aes(x=bmi, y=exprs, color=group, text=text)) +
      geom_point(size=1) +
      labs(title=paste(rv$plot_list$gene, " - ", rv$plot_list$gene_name))
  })
  
})
