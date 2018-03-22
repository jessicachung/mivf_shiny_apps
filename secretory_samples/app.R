############################################################
## Secretory samples analysis
############################################################

library(shiny)
library(DT)
library(ggplot2)
library(dplyr)
library(stringr)
library(reshape2)
library(plotly)
library(stats)
library(glmnet)
library(splines)
library(mgcv)

############################################################
## Load data

load("data/shiny.RData")

model_days <- phenotype$model_estimate %>% setNames(phenotype$sample_id)
path_days <- phenotype$avg_path %>% setNames(phenotype$sample_id)
batch <- phenotype$batch %>% factor
endo <- phenotype$endo %>% factor
endo_group <- phenotype$endo_group %>% factor
urs <- phenotype$urs %>% factor


############################################################
## Functions

fit_gam_model <- function(exprs, days, probe, 
                          formula="exprs ~ s(day_cycle, k=3)", 
                          gam_method="GCV.Cp", weights=NULL, 
                          newdata_covariates=list(), ...) {
  kwargs <- list(...)
  cycle_range <- seq(min(days), max(days), by=0.25)
  # Check if probe is valid and samples are valid
  stopifnot(probe %in% rownames(exprs))
  stopifnot(names(days) %in% colnames(exprs))
  
  # Get data into a tidy data frame
  dat <- data.frame(sample_id=names(days),
                    exprs=exprs[probe, names(days)],
                    day_cycle=days,
                    batch=batch)
  for (x in names(kwargs)) {
    dat[[x]] <- kwargs[[x]]
  }
  
  # Fit model
  tryCatch({
    fit <- gam(eval(parse(text=formula)), data=dat, family=gaussian(), 
               method=gam_method, weights=weights)
  }, error=function(e) {
    fit <- NULL
    warning(probe, " not fitted.\n")
  })
  if (is.null(fit)) {
    return(list(predict=NULL, residuals=NULL, coef=NULL, pseudo_r2=NULL))
  }  
  
  # Get coefficients
  fit_summary <- summary(fit)
  coef <- cbind(coef=fit_summary$p.coeff, p=fit_summary$p.pv)
  
  # Get predictions
  newdata <- data.frame(day_cycle=cycle_range)
  for (x in names(newdata_covariates)) {
    newdata[[x]] <- newdata_covariates[[x]]
  }
  predict <- predict(fit, newdata=newdata)
  names(predict) <- cycle_range
  residuals <- fit$residuals
  pseudo_r2 <- 1 - (fit$deviance/fit$null.deviance)
  return(list(predict=predict, residuals=residuals, coef=coef, pseudo_r2=pseudo_r2))
}

plot_model <- function(exprs, days, probe, predict, r2,
                       point_size=1, color="batch", ...) {
  kwargs <- list(...)
  cycle_range <- seq(min(days), max(days), by=0.25)
  dat <- data.frame(sample_id=names(days),
                    exprs=exprs[probe, names(days)],
                    day_cycle=days)
  for (x in names(kwargs)) {
    dat[[x]] <- kwargs[[x]]
  }
  pred <- data.frame(day_cycle=predict %>% names %>% as.numeric,
                     pred=predict)
  title <- paste0(probe, " (R^2 = ", round(r2, 3), ")")
  g1 <- ggplot(dat, aes(x=day_cycle, y=exprs)) +
    geom_point(aes_string(color=color, label="sample_id"), size=point_size) +
    geom_line(data=pred, aes(x=day_cycle, y=pred)) + 
    labs(title=title)
  g2 <- ggplot(dat, aes_string(x=color, y="exprs", color=color)) +
    geom_jitter(width=0.2, height=0) +
    labs(title=probe)
  return(list(g1=g1, g2=g2))
}

plot_points <- function(exprs, days, probe, point_size=1, color="batch", ...) {
  kwargs <- list(...)
  dat <- data.frame(sample_id=names(days),
                    exprs=exprs[probe, names(days)],
                    day_cycle=days)
  for (x in names(kwargs)) {
    dat[[x]] <- kwargs[[x]]
  }
  g1 <- ggplot(dat, aes(x=day_cycle, y=exprs)) +
    geom_point(aes_string(color=color), size=point_size) +
    labs(title=probe)
  g2 <- ggplot(dat, aes_string(x=color, y="exprs", color=color)) +
    geom_jitter(width=0.2, height=0) +
    labs(title=probe)
  return(list(g1=g1, g2=g2))
}


############################################################
## UI

ui <- fluidPage(
  titlePanel("Secretory samples"),
  
  fluidRow(
    column(3, wellPanel(
      # Static or interactive plot
      radioButtons("plot_type",
                   label = h4("Plot type:"),
                   choices = list("Static (ggplot2)" = "ggplot",
                                  "Interactive (plotly)" = "plotly"),
                   selected = "ggplot"),
      
      # Path day or model day
      radioButtons("day_type",
                   label = h4("Post ovulation day:"),
                   choices = list("Average pathology" = "path",
                                  "Molecular model estimate" = "model"),
                   selected = "path"),
      
      # Spline k
      sliderInput("spline_k",
                  label="Spline k (smooth term)",
                  min=2,
                  max=10,
                  value=4,
                  step=1),
      
      # Colour
      radioButtons("point_colour",
                   label = h4("Colour groups:"),
                   choices = list("URS" = "urs",
                                  "Endo" = "endo",
                                  "Endo group" = "endo_group",
                                  "Batch" = "batch"),
                   selected = "urs"),
      
      # Choose input type
      radioButtons("input_type",
                   label = h4("Input type:"),
                   choices = list("Probe ID from list" = 1,
                                  "Enter probe ID" = 2,
                                  "Random probe" = 3),
                   selected = 3),
      
      # Toggle input type for drop down menu
      conditionalPanel(
        condition="input.input_type == 1",
        selectInput("pa_group_select", 
                    label = h4("Illumina probe ID"), 
                    choices = list("Tissue Enriched (5)"="tissue_enriched",
                                   "Group Enriched (42)"="group_enriched",
                                   "Tissue Enhanced (136)"="tissue_enhanced"),
                    selected = "tissue_enriched"),
        conditionalPanel(
          condition="input.pa_group_select == 'tissue_enriched'",
          selectInput("probe_select_tissue_enriched", 
                      label = NULL, 
                      choices = pa_list[["tissue_enriched"]],
                      selected =  pa_list[["tissue_enriched"]][[1]])
        ),
        conditionalPanel(
          condition="input.pa_group_select == 'group_enriched'",
          selectInput("probe_select_group_enriched", 
                      label = NULL, 
                      choices =  pa_list[["group_enriched"]],
                      selected = pa_list[["group_enriched"]][[1]])
        ),
        conditionalPanel(
          condition="input.pa_group_select == 'tissue_enhanced'",
          selectInput("probe_select_tissue_enhanced", 
                      label = NULL, 
                      choices = pa_list[["tissue_enhanced"]],
                      selected = pa_list[["tissue_enhanced"]][[1]])
        )
      ),
      
      # Toggle input type for probe text input
      conditionalPanel(
        condition="input.input_type == 2",
        textInput("probe_text", 
                  label = h4("Illumina probe ID"), 
                  value = "ILMN_1774828")
      ),
      
      # Toggle for random
      conditionalPanel(
        condition="input.input_type == 3",
        sliderInput("r2_range",
                    label="Limit R^2 range",
                    min=0,
                    max=1,
                    value=c(0,1),
                    step=0.01),
        textInput("random_text", 
                  label = h4("Illumina probe ID"), 
                  value = "")
      ),
      
      # Dynamic submit button
      uiOutput("submit")
      
    ) # end wellPanel
    ), # end column
    
    column(
      9,
      column(8,
      # Plot
      conditionalPanel(
        condition="input.plot_type == 'ggplot'",
        plotOutput("ggplot_1a")
      ),
      
      conditionalPanel(
        condition="input.plot_type == 'plotly'",
        plotlyOutput("plotly_1a")
      ),
      
      hr(),
      
      # Cycle day corrected plots
      conditionalPanel(
        condition="input.plot_type == 'ggplot'",
        plotOutput("ggplot_2a")
      ),
      
      conditionalPanel(
        condition="input.plot_type == 'plotly'",
        plotlyOutput("plotly_2a")
      ),
      
      # Coefficients output
      conditionalPanel(
        condition="input.show_coefs == true",
        h4("Coefficients:"),
        tableOutput("coef")
      )
      ), # end column
      
      column(4,
             plotOutput("ggplot_1b"),
             hr(),
             plotOutput("ggplot_2b")
      )
        
    ) # end colum
  ), # end fluidRow
    
  hr(),
  
  # Probe info table
  h3("Selected probe:"),
  tableOutput("probe_info"),
  
  hr(),
  
  # Given gene name, search for probes
  textInput("gene_search", 
            label = h4("Enter a gene name to search for probes:"), 
            value = "VEZT"),
  
  # Probe search table
  tableOutput("gene_search_table")
  
)

############################################################
## Server

server <- function(input, output, session){
  
  # Set reactive values
  rv <- reactiveValues(
    probe_name = "ILMN_1748004",
    days=path_days
  )

  # Update action button text
  output$submit <- renderUI({
    label <- ifelse(input$input_type == "3", "Randomise", "Submit")
    actionButton("submit", label = label)
  })

  # Update probe name when submit button is pressed
  observeEvent({input$submit; input$input_type; input$pa_group_select; 
                input$probe_select_tissue_enriched; 
                input$probe_select_group_enriched; 
                input$probe_select_tissue_enhanced}, {
    if (input$input_type == "1") {
      # Probe from list
      rv$probe_name <- input[[paste0("probe_select_", input$pa_group_select)]]
    } else if (input$input_type == "2") {
      # Probe from text box
      rv$probe_name <- input$probe_text
    } else {
      # Random probe (overwrite seed from jitter)
      set.seed(Sys.time() %>% as.numeric)
      in_range <- names(r2)[between(r2, input$r2_range[1], input$r2_range[2])]
      rv$probe_name <- sample(in_range, 1)
      # Update text box
      updateTextInput(session, "random_text",
                      value = paste(rv$probe_name)
      )
    }
  })
  
  # Update POD selection
  observeEvent({input$day_type}, {
    if (input$day_type == "path") {
      rv$days <- path_days
    } else {
      rv$days <- model_days
    }
  })
  
  # Update plots when probe_name changes
  observeEvent({rv$probe_name; rv$days; input$point_colour; input$plot_type;
                input$spline_k}, {
    # Set point size smaller if plotting with Plotly
    if (input$plot_type == "ggplot") {
      ps <- 1.5
    } else {
      ps <- 0.8
    }
    message("plot: ", rv$probe_name)
    formula <- sprintf("exprs ~ s(day_cycle, k=%d)", input$spline_k)
    model <- fit_gam_model(exprs, days=rv$days, probe=rv$probe_name,
                  formula=formula, gam_method="REML")
    
    # Calculated day-corrected exprs (Always use model day to correct, not path day)
    exp <- model$predict[as.character(model_days)]
    probe_mean <- mean(exp)
    obs <- exprs[rv$probe_name,names(rv$days)]
    day_corrected <- obs - exp
    day_corrected <- matrix(day_corrected + probe_mean, nrow=1) %>%
      magrittr::set_colnames(names(rv$days)) %>%
      magrittr::set_rownames(rv$probe_name)
    
    # Plot expression across the cycle
    x <- plot_model(exprs=exprs, days=rv$days, probe=rv$probe_name, point_size=ps,
                    predict=model$predict, r2=model$pseudo_r2, 
                    batch=batch, endo=endo, endo_group=endo_group, 
                    urs=urs, color=input$point_colour)
    rv$plot_1a <- x$g1
    rv$plot_1b <- x$g2
    
    # Plot day-corrected expression
    y <- plot_points(exprs=day_corrected, days=rv$days, probe=rv$probe_name, 
                     point_size=ps, batch=batch, endo=endo, 
                     endo_group=endo_group, urs=urs, color=input$point_colour)
    rv$plot_2a <- y$g1
    rv$plot_2b <- y$g2
  })
  
  ## ggplot
  output$ggplot_1a <- renderPlot({
    rv$plot_1a
  })

  output$ggplot_1b <- renderPlot({
    rv$plot_1b
  })
  
  output$ggplot_2a <- renderPlot({
    rv$plot_2a
  })
  
  output$ggplot_2b <- renderPlot({
    rv$plot_2b
  })
  
  # plotly
  output$plotly_1a <- renderPlotly({
    rv$plot_1a
  })
  
  output$plotly_2a <- renderPlotly({
    rv$plot_2a
  })
  
  # Probe info table
  output$probe_info <- renderTable({
    probe_info %>% filter(IlluminaID == rv$probe_name)
  })

  # Search table
  output$gene_search_table <- renderTable({
    if (nchar(input$gene_search) > 1) {
      probe_info %>% filter(str_detect(Symbol, input$gene_search))
    } else {
      NULL
    }
  })
  
}

############################################################
## Run

shinyApp(ui=ui, server=server)
