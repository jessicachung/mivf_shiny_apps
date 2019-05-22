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

# Separate m/z and retention time
mz <- feature_data$mz_rt %>% str_split("_") %>% sapply(function(x) x[1]) %>%
  as.numeric
rt <- feature_data$mz_rt %>% str_split("_") %>% sapply(function(x) x[2]) %>%
  as.numeric
feature_data <- feature_data %>% 
  mutate(pass_threshold=(! is.na(shapiro_pval)),
         # shapiro_pval=ifelse(is.na(shapiro_pval), Inf, shapiro_pval),
         # batch_pval=ifelse(is.na(batch_pval), Inf, batch_pval),
         mz=mz, rt=rt) %>%
  select(mz_rt, mz, rt, pass_threshold, everything())
spike_in_controls <- 
  c("487.3456_2.96", "518.3144_2.67", "529.3917_2.87", "610.5316_10.75", 
    "675.5366_7.45", "675.6492_18.15", "712.6343_14.7", "753.5942_9.11", 
    "755.5653_10.61", "781.4483_13.88", "829.7874_17.56", "846.744_16.99"
  )

############################################################
## Functions
############################################################

############################################################
## Server
############################################################


shinyServer(function(input, output) {
  
  # Set reactive values
  rv <- reactiveValues(
    exprs = imputed_matrix,
    phenotype = sample_data,
    feature_data = feature_data,
    selected_lipid = "399.9738_1.25",
    plot_data = NULL
  )
  
  # Update exprs matrix when selected
  observeEvent({input$exprs_dataset}, {
    message("Updating exprs matrix...")
    if (input$exprs_dataset == "raw") {
      rv$exprs <- ms_matrix
    } else if (input$exprs_dataset == "qnorm") {
      rv$exprs <- imputed_matrix
    } else if (input$exprs_dataset == "combat") {
      rv$exprs <- combat_cc
    } else if (input$exprs_dataset == "metablo") {
      rv$exprs <- metablo_matrix
    }
  })
  
  # Update feature dataframe when selected
  observeEvent({input$lipid_subset}, {
    message("Updating feature table...")
    if (input$lipid_subset == "all") {
      rv$feature_data <- feature_data
    } else if (input$lipid_subset == "random_100") {
      random_indices <- sample(seq_len(nrow(feature_data)), 100) %>% sort
      rv$feature_data <- feature_data[random_indices,]
    } else if (input$lipid_subset == "random_1000") {
      random_indices <- sample(seq_len(nrow(feature_data)), 1000) %>% sort
      rv$feature_data <- feature_data[random_indices,]
    } else if (input$lipid_subset == "controls") {
      rv$feature_data <- feature_data %>% filter(mz_rt %in% spike_in_controls)
    }
  })
  
  # Output feature table
  observeEvent({rv$feature_data}, {
    output$feature_table <- renderDataTable({
      datatable(rv$feature_data %>% format(digits=3), 
                options=list(pageLength=10), selection="single")
    })
  })

  # Update selected lipid when row is clicked
  observeEvent({input$feature_table_cell_clicked}, {
    if (!is.null(input$feature_table_cell_clicked$row)) {
      message(paste0("Row clicked: ", input$feature_table_cell_clicked$row))
      rv$selected_lipid <- rv$feature_data[input$feature_table_cell_clicked$row,"mz_rt"]
    }
  })

  observeEvent({input$refresh}, {
    # Sometime the server hangs when plotting
    message("Refresh")
    rv$plot_data <- rv$plot_data
  })
  
  # Get data for plotting expression
  observeEvent({rv$selected_lipid; input$exprs_dataset}, {
    message("Updating plot data: ", rv$selected_lipid)
    # Get expression and phenotype data
    if (! rv$selected_lipid %in% rownames(rv$exprs)) {
      rv$plot_data <- NULL
      return(NULL)
    }
    dat <- data.frame(sample_id=colnames(rv$exprs),
                      exprs=rv$exprs[rv$selected_lipid,colnames(rv$exprs)]) %>%
      merge(rv$phenotype, by="sample_id", all.x=TRUE)

    rv$plot_data <- dat
  })

  # Plot 1 with ggplot
  output$plot_1_ggplot <- renderPlot({
    ggplot(rv$plot_data, aes(y=exprs)) +
      geom_point(aes_string(x=input$plot_1_x_axis, color=input$plot_1_colour), size=2) +
      labs(title=paste(rv$selected_lipid))
  })
  
  # Plot 2 with ggplot
  output$plot_2_ggplot <- renderPlot({
    ggplot(rv$plot_data, aes(x=exprs)) +
      geom_histogram(aes_string(fill=input$plot_1_colour), bins=50) +
      labs(title=paste(rv$selected_lipid))
  })

  # Plot 1 with plotly
  output$plot_1_plotly <- renderPlotly({
    ggplot(rv$plot_data, aes(y=exprs)) +
      geom_point(aes_string(x=input$plot_1_x_axis, color=input$plot_1_colour, text=input$sample_label)) +
      labs(title=paste(rv$selected_lipid))
  })

})
