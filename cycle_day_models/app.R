############################################################
## Model expression along 28 day cycle
############################################################

# Plot expression a given probe along cycle days and fit a line with a
# polynomial model or spline.

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

# setwd("~/Work/2016_mivf/shiny_apps/cycle_day_models")

# Load data
load("data/cycle_data.RData")
probes <- rownames(combat_2batch_exprs)
phenotype <- phenotype %>% mutate(pathology_day=day_cycle)

############################################################
## Functions

get_tidy_cycle_data <- function(exprs, pheno, probe) {
  # Melt data
  dat <- exprs[probe,,drop=FALSE] %>% melt %>% 
    S4Vectors::rename(Var2="sample_id")
  stopifnot(dat[,"sample_id"] == pheno[,"sample_id"])
  original_dat <- cbind(dat, pheno %>% dplyr::select(-sample_id)) %>% 
    filter(!is.na(day_cycle)) %>% mutate(original=TRUE)
  return(original_dat)
}

extend_days_data <- function(dat, extend_days) {
  # Extend days backwards and forwards
  if (extend_days > 0 && extend_days <= 28) {
    forwards <- dat %>% filter(day_cycle < extend_days) %>% 
      mutate(day_cycle=day_cycle + 28, original=FALSE)
    backwards <- dat %>% filter(28 - day_cycle < extend_days) %>% 
      mutate(day_cycle=day_cycle - 28, original=FALSE)
    extended_dat <- rbind(backwards, dat, forwards)
    return(extended_dat)
  } else {
    stop("extend_days must be between 0 and 28")
  }
}

get_outlier_samples <- function(residuals, pheno, sd, band_size) {
  # Get data frame of outlier samples
  stopifnot(pheno[,"sample_id"] == names(residuals))
  df <- pheno %>% mutate(residuals=round(residuals,4)) %>%
    filter(abs(residuals) - band_size * sd > 0) %>% 
    arrange(day_cycle) %>% 
    dplyr::select(-text, -afs_score_log, -sample_type, -sample_section, -day_cycle) %>%
    dplyr::select(sample_id, contains("day"), everything())
  return(df)
}

calculate_R2 <- function(residuals, dat) {
  # Calculate R^2
  SSR <- residuals^2 %>% sum
  SST <- (dat$value - mean(dat$value))^2 %>% sum
  R2 <- round(1 - SSR/SST, digits=3)
  return(R2)
}

expression_plot <- function(dat, predict, band_size, sd, jitter_scale, 
                            color_str, point_size, probe, R2, cycle_range) {
  # Plot expression and curve with ggplot
  if (jitter_scale > 0) {
    set.seed(0)
    # Add jitter to avoid overplotting
    jitter <- runif(nrow(dat),jitter_scale * -1,jitter_scale)
    dat <- mutate(dat, day_cycle=day_cycle+jitter)
  }
  predict <- data.frame(day_cycle=cycle_range, predict=predict, 
                        ymin=predict-band_size*sd, ymax=predict+band_size*sd)
  title <- paste0(probe, " (R^2 = ", R2, ")")
  g <- ggplot(dat, aes(x=day_cycle, y=value)) +
    geom_point(aes_string(text="text", color=color_str), size=point_size) + 
    geom_line(data=predict, aes(x=day_cycle, y=predict), alpha=0.3, size=1.5) +
    geom_ribbon(data=predict, aes(x=day_cycle, y=predict, ymin=ymin, ymax=ymax), alpha=0.1) +
    labs(title=title)
  return(g)
}

no_model <- function(exprs, pheno, probe, jitter_scale=1, color_str="endo",
                     point_size=1) {
  dat <- get_tidy_cycle_data(exprs, pheno, probe)
  if (jitter_scale > 0) {
    set.seed(0)
    # Add jitter to avoid overplotting
    jitter <- runif(nrow(dat),jitter_scale * -1,jitter_scale)
    dat <- mutate(dat, day_cycle=day_cycle+jitter)
  }
  g <- ggplot(dat, aes(x=day_cycle, y=value)) +
    geom_point(aes_string(text="text", color=color_str), size=point_size) + 
    labs(title=probe)
  return(list(plot=g, outliers=NULL, coef=NULL, R2=NULL))
}

fit_polynomial_model <- function(exprs, pheno, probe, extend_days, poly_degree,
    poly_raw, cycle_range=seq(0,28,by=0.5), jitter_scale=1, band_size=2, 
    color_str="endo", point_size=1) {
  # Check if probe is valid
  if (! probe %in% rownames(exprs)) return(NULL)
  
  # Get data into a tidy data frame
  original_dat <- get_tidy_cycle_data(exprs, pheno, probe)
  if (extend_days) {
    extended_dat <- extend_days_data(original_dat, extend_days)
  } else {
    extended_dat <- original_dat
  }
  
  # Fit model
  x <- poly(extended_dat$day_cycle, poly_degree, raw=poly_raw)
  poly_coefs <- attr(x,"coefs")
  fit <- lm(extended_dat$value ~ x)
  
  # Prediction at each day in cycle
  x <- poly(cycle_range, poly_degree, coefs=poly_coefs, raw=poly_raw)
  predict <- predict(fit, x)
  
  # Get coefficients
  coef <- fit$coefficients %>% t %>% as.data.frame
  
  # Get residuals
  residuals <- fit$residuals[which(extended_dat$original)]
  names(residuals) <- extended_dat %>% filter(original) %>% .$sample_id
  
  # Calculate R^2
  sd <- sd(original_dat$value)
  R2 <- calculate_R2(residuals=residuals, dat=original_dat)
  
  # Get outlier samples
  outliers <- get_outlier_samples(residuals=residuals, pheno=pheno, sd=sd,
                                  band_size=band_size)
  
  # Plot
  g <- expression_plot(dat=original_dat, predict=predict, band_size=band_size,
                       sd=sd, jitter_scale=jitter_scale, color_str=color_str, 
                       point_size=point_size, probe=probe, R2=R2, 
                       cycle_range=cycle_range)
  
  return(list(plot=g, outliers=outliers, coef=coef, R2=R2))
}

fit_en_polynomial_model <- function(exprs, pheno, probe, extend_days, poly_degree,
    poly_raw, use_weights = FALSE, elastic_alpha, cycle_range=seq(0,28,by=0.5), jitter_scale=1, 
    band_size=2, color_str="endo", point_size=1) {
  # Check if probe is valid
  if (! probe %in% rownames(exprs)) return(NULL)
  
  # Get data into a tidy data frame
  original_dat <- get_tidy_cycle_data(exprs, pheno, probe)
  if (extend_days) {
    extended_dat <- extend_days_data(original_dat, extend_days)
  } else {
    extended_dat <- original_dat
  }
  
  # Get weights
  if (use_weights) {
    day_cycle_weights <- extended_dat %>% group_by(day_cycle) %>% summarise(weight=1/n())
    extended_dat <- merge(extended_dat, day_cycle_weights, by="day_cycle")
  } else {
    extended_dat <- extended_dat %>% mutate(weight=1)
  }
  
  # Fit model
  x <- poly(extended_dat$day_cycle, poly_degree, raw=poly_raw)
  poly_coefs <- attr(x,"coefs")
  cv_fit <- cv.glmnet(x, extended_dat$value, weights=extended_dat$weight,
                      family="gaussian", alpha=elastic_alpha)
  
  # Prediction at each day in cycle
  x <- poly(cycle_range, poly_degree, coefs=poly_coefs, raw=poly_raw)
  predict <- predict(cv_fit, newx=x, s="lambda.min")[,1]
  
  # Get coefficients
  coef <- coef(cv_fit, s = "lambda.min") %>% as.matrix %>% t %>% as.data.frame
  coef[coef == 0] <- "."
  
  # Get residuals
  x <- poly(original_dat$day_cycle, poly_degree, coefs=poly_coefs, raw=poly_raw)
  residuals <- predict(cv_fit, newx=x, s="lambda.min")[,1] - original_dat$value
  names(residuals) <- original_dat %>% .$sample_id
  
  # Calculate R^2
  sd <- sd(original_dat$value)
  R2 <- calculate_R2(residuals=residuals, dat=original_dat)
  
  # Get outlier samples
  outliers <- get_outlier_samples(residuals=residuals, pheno=pheno, sd=sd,
                                  band_size=band_size)
  
  # Plot
  g <- expression_plot(dat=original_dat, predict=predict, band_size=band_size,
                       sd=sd, jitter_scale=jitter_scale, color_str=color_str, 
                       point_size=point_size, probe=probe, R2=R2,
                       cycle_range=cycle_range)
  
  return(list(plot=g, outliers=outliers, coef=coef, R2=R2))
}


fit_spline_model <- function(exprs, pheno, probe, extend_days, spline_df,
    cycle_range=seq(0,28,by=0.5), jitter_scale=1, band_size=2, color_str="endo", 
    point_size=1) {
  # Check if probe is valid
  if (! probe %in% rownames(exprs)) return(NULL)
  
  # Get data into a tidy data frame
  original_dat <- get_tidy_cycle_data(exprs, pheno, probe)
  if (extend_days) {
    extended_dat <- extend_days_data(original_dat, extend_days)
  } else {
    extended_dat <- original_dat
  }
  
  # Fit model
  fit <- lm(value ~ ns(day_cycle, df=spline_df), data=extended_dat)
  
  # Prediction at each day in cycle
  pred <- data.frame(day_cycle=cycle_range)
  predict <- predict(fit, pred)
  
  # Get coefficients
  coef <- fit$coefficients %>% t %>% as.data.frame
  colnames(coef) <- colnames(coef) %>% str_replace("ns([^0-9]+)", "x")
  
  # Get residuals
  residuals <- fit$residuals[which(extended_dat$original)]
  names(residuals) <- extended_dat %>% filter(original) %>% .$sample_id
  
  # Calculate R^2
  sd <- sd(original_dat$value)
  R2 <- calculate_R2(residuals=residuals, dat=original_dat)
  
  # Get outlier samples
  outliers <- get_outlier_samples(residuals=residuals, pheno=pheno, sd=sd,
                                  band_size=band_size)
  
  # Plot
  g <- expression_plot(dat=original_dat, predict=predict, band_size=band_size,
                       sd=sd, jitter_scale=jitter_scale, color_str=color_str, 
                       point_size=point_size, probe=probe, R2=R2,
                       cycle_range=cycle_range)
  
  return(list(plot=g, outliers=outliers, coef=coef, R2=R2))
}

############################################################
## UI

ui <- fluidPage(
  titlePanel("Model expression across 28 day cycle"),
  
  fluidRow(
    column(4, wellPanel(
      # Static or interactive plot
      radioButtons("plot_type",
                   label = h4("Plot type:"),
                   choices = list("Static (ggplot2)" = 1,
                                  "Interactive (plotly)" = 2),
                   selected = 1),
      
      hr(),
      
      # Type of model to fit
      radioButtons("model_type",
                   label = h4("Model:"),
                   choices = list("Polynomial" = 1,
                                  "Polynomial with elastic net" = 2,
                                  "B-spline" = 3,
                                  "None" = 4),
                   selected = 2),
      
      # Polynomial degrees
      conditionalPanel(
        condition="input.model_type == 1 || input.model_type == 2",
        sliderInput("poly_degree",
                    label = "Polynomial degree:",
                    min=1,
                    max=20,
                    step=1,
                    value=13)
      ),
      
      # Elastic net
      conditionalPanel(
        condition="input.model_type == 2",
        sliderInput("elastic_alpha",
                    label = "Elastic net alpha:",
                    min=0,
                    max=1,
                    step=0.05,
                    value=0.5)
      ),
      
      # Spline parameters
      conditionalPanel(
        condition="input.model_type == 3",
        sliderInput("spline_df",
                    label = "Spline df:",
                    min=1,
                    max=20,
                    step=1,
                    value=6)
      ),
      
      # Advanced options
      checkboxInput("advanced",
                    label = "Show advanced options",
                    value = FALSE),
      
      # Polynomial raw
      conditionalPanel(
        condition="input.advanced == true && 
                   (input.model_type == 1 || input.model_type == 2)",
        checkboxInput("poly_raw",
                      label = "Raw polynomials (not orthogonal)",
                      value=FALSE)
      ),
      
      # Use weights
      conditionalPanel(
        condition="input.advanced == true && input.model_type == 2",
        checkboxInput("use_weights",
                      label = "Weighted equally for each cycle day",
                      value = FALSE)
      ),
      
      # Show table of coefficients
      conditionalPanel(condition="input.advanced == true",
        checkboxInput("show_coefs",
                      label = "Display coefficients",
                      value = TRUE)
      ),
      
      # Expression values
      conditionalPanel(condition="input.advanced == true",
        selectInput("expression_values",
                    label = "Expression values:",
                    choices = list("Combat 2batch expression" = 1,
                                   "Combat 2batch cycle expression" = 2),
                    selected = 2)
      ),
      
      # 28-day cycle values
      conditionalPanel(condition="input.advanced == true",
                       selectInput("day_cycle_values",
                                   label = "28-day cycle values:",
                                   choices = list("Updated cycle values" = 1,
                                                  "Original pathology" = 2,
                                                  "SVM classification" = 3),
                                   selected = 1)
      ),
      
      # Colour
      conditionalPanel(
        condition="input.advanced == true",
        selectInput("point_color",
                    label = "Point colours:",
                    choices = list("Endo status" = "endo",
                                   "Study" = "study",
                                   "AFS score (log)" = "afs_score_log"),
                    selected = 1)
        ),
      
      # Curve band size
      conditionalPanel(condition="input.advanced == true",
        sliderInput("band_size",
                    label = "Band size (std dev):",
                    min=0,
                    max=3,
                    step=0.1,
                    value=2)
      ),
      
      # Jitter to avoid overplotting
      conditionalPanel(condition="input.advanced == true",
        sliderInput("jitter_scale",
                    label = "Cycle day jitter:",
                    min=0,
                    max=1,
                    step=0.05,
                    value=0.1)
      ),
      
      # Extend days
      conditionalPanel(condition="input.advanced == true",
        sliderInput("extend_days",
                    label="Extend days forward and backwards:",
                    min=0,
                    max=28,
                    step=0.5,
                    value=12)
      ),
      
      # Show outlier table
      checkboxInput("show_outliers",
                    label = "Display outlier samples in table",
                    value = TRUE)
      
    )),
    
    column(7, wellPanel(
      fluidRow(
        column(6,
          # Choose input type
          radioButtons("input_type",
                       label = h4("Input type:"),
                       choices = list("Probe ID from list" = 1,
                                      "Enter probe ID" = 2,
                                      "Random probe" = 3),
                       selected = 3)
        ),
        
        column(5,
          # Toggle input type for drop down menu
          conditionalPanel(
            condition="input.input_type == 1",
            selectInput("probe_select", 
                        label = h4("Illumina probe ID"), 
                        choices = probe_list,
                        selected = probe_list[[1]])
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
            textInput("random_text", 
                      label = h4("Illumina probe ID"), 
                      value = "")
          ),
          
          # Dynamic submit button
          uiOutput("submit")
          
        )
      )),
      
      hr(),
      
      # Plot
      conditionalPanel(
        condition="input.plot_type == 1",
        plotOutput("ggplot")
      ),
      
      conditionalPanel(
        condition="input.plot_type == 2",
        plotlyOutput("plotly")
      ),
      
      hr(),
      
      # Coefficients output
      conditionalPanel(
        condition="input.show_coefs == true",
        h4("Coefficients:"),
        tableOutput("coef")
      )
    )
    
  ),
  
  hr(),
  
  # Outlier samples
  conditionalPanel(
    condition="input.show_outliers == true",
    h3("Outlier samples:"),
    dataTableOutput("outlier_samples"),
    hr()
  ),
  
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
    probe_name = probe_list[[1]],
    exprs = combat_2batch_cycle_exprs,
    phenotype = phenotype,
    model = NA
  )
  
  # Update action button text
  output$submit <- renderUI({
    label <- ifelse(input$input_type == "3", "Randomise", "Submit")
    actionButton("submit", label = label)
  })
  
  # Update probe name when submit button is pressed
  observeEvent({input$submit; input$probe_select}, {
    if (input$input_type == "1") {
      # Probe from list
      rv$probe_name <- input$probe_select
    } else if (input$input_type == "2") {
      # Probe from text box
      rv$probe_name <- input$probe_text
    } else {
      # Random probe (overwrite seed from jitter)
      set.seed(Sys.time() %>% as.numeric)
      rv$probe_name <- sample(probes, 1)
      # Update text box
      updateTextInput(session, "random_text",
                      value = paste(rv$probe_name)
      )
    }
  })
  
  # Update expression values when changed
  observeEvent({input$expression_values}, {
    if (input$expression_values == "1") {
      rv$exprs <- combat_2batch_exprs
    } else {
      rv$exprs <- combat_2batch_cycle_exprs
    }
  })
  
  # Update day cycle when changed
  observeEvent({input$day_cycle_values}, {
    if (input$day_cycle_values == "1") {
      rv$phenotype$day_cycle <- rv$phenotype$day_cycle_v2
    } else if (input$day_cycle_values == "2") {
      rv$phenotype$day_cycle <- rv$phenotype$pathology_day
    } else {
      rv$phenotype$day_cycle <- rv$phenotype$svm_predicted_day
    }
  })
  
  # Update model when arguments update
  observeEvent({rv$probe_name; rv$exprs; rv$phenotype; input$plot_type; input$model_type; 
    input$poly_degree; input$poly_raw; input$elastic_alpha; input$spline_df; 
    input$jitter_scale; input$band_size; input$point_color; input$extend_days;
    input$use_weights}, {
      # Set point size smaller if plotting with Plotly
      if (input$plot_type == "1") {
        ps <- 1
      } else {
        ps <- 0.5
      }
      if (input$model_type == "1") {
        # Polynomial
        rv$model <- fit_polynomial_model(exprs=rv$exprs, pheno=rv$phenotype, probe=rv$probe_name,
                        extend_days=input$extend_days, poly_degree=input$poly_degree,
                        poly_raw=input$poly_raw, jitter_scale=input$jitter_scale, 
                        band_size=input$band_size, color_str=input$point_color, point_size=ps)
      } else if (input$model_type == "2") {
        # Polynomial with elastic net
        rv$model <- fit_en_polynomial_model(exprs=rv$exprs, pheno=rv$phenotype, probe=rv$probe_name,
                        extend_days=input$extend_days, poly_degree=input$poly_degree,
                        poly_raw=input$poly_raw, elastic_alpha=input$elastic_alpha,
                        jitter_scale=input$jitter_scale, band_size=input$band_size, 
                        color_str=input$point_color, point_size=ps, use_weights=input$use_weights)
      } else if (input$model_type == "3") {
        # Splines
        rv$model <- fit_spline_model(exprs=rv$exprs, pheno=rv$phenotype, probe=rv$probe_name,
                        extend_days=input$extend_days, spline_df=input$spline_df,
                        jitter_scale=input$jitter_scale, band_size=input$band_size, 
                        color_str=input$point_color, point_size=ps)
      } else {
        rv$model <- no_model(exprs=rv$exprs, pheno=rv$phenotype, probe=rv$probe_name,
                             color_str=input$point_color, point_size=ps)
      }
    })
  
  output$ggplot <- renderPlot({
    message("ggplot: ", rv$probe_name)
    rv$model$plot
  })
  
  output$plotly <- renderPlotly({
    message("plotly: ", rv$probe_name)
    rv$model$plot
  })
  
  # Coefficients output
  output$coef <- renderTable({
    rv$model$coef %>% format(digits=3)
  })
  
  # Sample outliers table
  output$outlier_samples <- renderDataTable({
    datatable(rv$model$outliers, options=list(pageLength=15))
  })
  
  # Probe info table
  output$probe_info <- renderTable({
    probe_df %>% filter(IlluminaID == rv$probe_name)
  })
  
  # Search table
  output$gene_search_table <- renderTable({
    probe_df %>% filter(Symbol == input$gene_search)
  })
  
}

############################################################
## Run

shinyApp(ui=ui, server=server)
