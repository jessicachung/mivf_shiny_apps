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

load("data/cc_combat_2batch_cycle_exprs.rda")
all_phenotype <- readRDS("data/phenotype.rds")
probe_info <- readRDS("data/probe_info.rds")

# Get samples with BMI values
phenotype <- all_phenotype %>% filter(include)

# Filter probes
good_probes <- probe_info %>% 
  filter(pass_threshold & pass_quality) %>% .$IlluminaID

# Subset expression matrix
exprs <- cc_combat_2batch_cycle_exprs[good_probes, phenotype$sample_id]

# Dimensions
message(paste0("Expression matrix dimensions: ", 
               nrow(exprs)), " x ", ncol(exprs))

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
    selection_method = NULL,
    group_1_samples = NULL,
    group_2_samples = NULL,
    selected_probe = NA,
    results = NULL,
    plot_list = NULL
  )
  
  # Update rv groups when ranges/checkboxes change
  observeEvent({input$selection; input$slider_1; input$slider_2;
    input$checkbox_1; input$checkbox_2}, {
      message("Updating groups...")
      if (input$selection == "slider") {
        rv$group_1_samples <- get_slider_samples(phenotype, input$slider_1)
        rv$group_2_samples <- get_slider_samples(phenotype, input$slider_2)
      } else {
        rv$group_1_samples <- get_checkbox_samples(phenotype, input$checkbox_1)
        rv$group_2_samples <- get_checkbox_samples(phenotype, input$checkbox_2)
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
  observeEvent({rv$group_1_samples; rv$group_2_samples; input$n_top_table}, {
    message("Performing DGE analysis...")
    
    # Subset data
    pheno <- phenotype %>% select(sample_id:subclass, batch, text) %>%
      filter(sample_id %in% c(rv$group_1_samples, rv$group_2_samples)) %>%
      mutate(group=ifelse(sample_id %in% rv$group_1_samples, "Group 1", "Group 2"))
    
    # Check there's at least 3 samples in each group
    if (sum(pheno$group == "Group 1") < 3 | 
        sum(pheno$group == "Group 2") < 3) {
      message("Too few samples in groups.")
      rv$results <- NULL
      rv$selected_probe <- NULL
      
    } else {
      # Setup design
      design <- model.matrix(~group, data = pheno)
      
      # Fit model
      lm <- lmFit(exprs[,pheno$sample_id], design)
      fit <- eBayes(lm)
      top_table <- topTable(fit, coef=2, n=input$n_top_table) %>% 
        tibble::rownames_to_column(var="IlluminaID") %>%
        select(-t, -B) %>%
        merge(probe_info %>% select(IlluminaID, SymbolReannotated), 
              by="IlluminaID") %>% 
        arrange(P.Value)
      
      rv$results <- list(top_table=top_table, pheno=pheno)
      rv$selected_probe <- top_table[1,"IlluminaID"]
    }
  })

  # Output top table
  output$top_table <- renderDataTable({
    datatable(rv$results$top_table %>% format(digits=3), 
              options=list(pageLength=5), selection="single")
  })
  
  # Update selected probe when row is clicked
  observeEvent({input$top_table_cell_clicked}, {
    if (!is.null(input$top_table_cell_clicked$row)) {
      message(paste0("Row clicked: ", input$top_table_cell_clicked$row))
      rv$selected_probe <- rv$results$top_table[input$top_table_cell_clicked$row,"IlluminaID"]
    }
  })
  
  # Get data for plotting expression
  observeEvent({rv$results; rv$selected_probe}, {
    message("Updating plot data: ", rv$selected_probe)
    
    # Get probe name
    probe_name <- probe_info %>% filter(IlluminaID == rv$selected_probe) %>%
      .$SymbolReannotated
    
    # Get expression and phenotype data
    dat <- data.frame(sample_id=colnames(exprs), 
                      exprs=exprs[rv$selected_probe,]) %>%
      merge(phenotype %>% select(sample_id, bmi, text), by="sample_id", all.x=TRUE) %>%
      merge(rv$results$pheno %>% select(sample_id, group), by="sample_id", all.x=TRUE)
    
    rv$plot_list <- list(dat=dat, probe=rv$selected_probe, probe_name=probe_name)
  })
  
  # Plot BMI vs expression with ggplot
  output$ggplot <- renderPlot({
    ggplot(rv$plot_list$dat, aes(x=bmi, y=exprs, color=group)) +
      geom_point(size=2) +
      labs(title=paste(rv$plot_list$probe, " - ", rv$plot_list$probe_name))
  })

  # Plot BMI vs expression with plotly
  output$plotly <- renderPlotly({
    ggplot(rv$plot_list$dat, aes(x=bmi, y=exprs, color=group, text=text)) +
      geom_point(size=1) +
      labs(title=paste(rv$plot_list$probe, " - ", rv$plot_list$probe_name))
  })
  
})
