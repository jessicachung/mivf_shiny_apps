############################################################
## Expression along 28 day cycle
############################################################

# Plot expression a given probe along cycle days with LOESS smoothing

## TODO: 
## Add random button
## R^2 function

library(shiny)
library(ggplot2)
library(dplyr)
library(stringr)
library(reshape2)
library(plotly)
library(stats)

# setwd("~/Work/2016_mivf/shiny_apps/cycle_day_plots")

# Expression data
load("data/combat_exprs.rda")

# Phenotype info
phenotype <- read.table("data/cycle_phenotype_2016-08-11.tsv", header=TRUE, sep="\t", 
                        stringsAsFactors=FALSE) %>%
  mutate(endo=factor(endo), afs_score_log=log2(afs_score+0.01))

# Probe info
probe_df <- read.table("data/illumina_v4_annotation_with_detection.tsv", header=TRUE, sep="\t", 
                       stringsAsFactors=FALSE) %>% 
  S4Vectors::rename(X10.="10%", X50.="50%", X90.="90%") %>%
  dplyr::select(IlluminaID, SymbolReannotated, GenomicLocation, ProbeQuality, MeanDetectionPVal:`90%`) %>%
  S4Vectors::rename(SymbolReannotated="Symbol")

# Subset expression data
combat_exprs <- combat_exprs[,phenotype[,"sample_id"]]
cycle <- phenotype[,"day_cycle"]

# Probes of interest from Jane
probes_raw <- "ILMN_1764096
ILMN_3194087
ILMN_3238259
ILMN_3272768
ILMN_1734552
ILMN_1652431
ILMN_1696183
ILMN_1743290
ILMN_1786015
ILMN_1807529
ILMN_2091454
ILMN_2169736
ILMN_2367126
ILMN_1740706
ILMN_2060719
ILMN_1784217
ILMN_1729033
ILMN_1782743"

# Make list of probes of interest for drop down menu
probes_of_interest <- probes_raw %>% str_split("\n") %>% unlist
probe_list <- list()
for (p in probes_of_interest) {
  probe_list[[p]] <- p
}

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

get_outlier_samples <- function(fit, dat, pheno, sd, band_size) {
  # Get outliers
  residuals <- fit$residuals[between(fit$x, 0.1, 28)]
  names(residuals) <- dat[,"sample_id"]
  indices <- which(abs(residuals) - band_size * sd > 0)
  outliers <- pheno[indices,] %>% arrange(day_cycle) %>% 
    dplyr::select(-text, -afs_score_log)
  return(outliers)
}

calculate_R2 <- function() {
  #TODO
  return("NA")
}

fit_loess_model <- function(exprs, pheno, probe, extend_days, loess_span=0.6, 
                            jitter_scale=1, band_size=2, color_str="endo", 
                            point_size=1) {
  if (! probe %in% rownames(exprs)) return(NULL)
  
  # Get data into a tidy data frame
  original_dat <- get_tidy_cycle_data(exprs, pheno, probe)
  if (extend_days) {
    extended_dat <- extend_days_data(original_dat, extend_days)
    cycle_range <- seq(0,28,by=0.5)
  } else {
    extended_dat <- original_dat
    cycle_range <- seq(min(cycle),max(cycle),by=0.5)
  }
  
  # Fit loess model
  fit <- loess(value ~ day_cycle, extended_dat, span=loess_span, model=TRUE)
  
  # Prediction at each day in cycle
  predict <- predict(fit, cycle_range)
  
  # Calculate R^2
  sd <- sd(original_dat[,"value"])
  R2 <- calculate_R2()
  
  # Get outlier samples
  outliers <- get_outlier_samples(fit, dat=original_dat, pheno=pheno, sd=sd, 
                                  band_size=band_size)
  
  # Jitter
  if (jitter_scale > 0) {
    set.seed(0)
    jitter <- runif(nrow(original_dat),jitter_scale * -1,jitter_scale)
    original_dat <- mutate(original_dat, day_cycle=day_cycle+jitter)
  }
  
  predict <- data.frame(day_cycle=cycle_range, predict=predict, 
                        ymin=predict-band_size*sd, ymax=predict+band_size*sd)
  title <- paste0(probe, " (R^2 = ", R2, ")")
  g <- ggplot(original_dat, aes(x=day_cycle, y=value)) +
    geom_point(aes_string(text="text", color=color_str), size=point_size) + 
    geom_line(data=predict, aes(x=day_cycle, y=predict), alpha=0.3, size=1.5) +
    geom_ribbon(data=predict, aes(x=day_cycle, y=predict, ymin=ymin, ymax=ymax), alpha=0.1) +
    labs(title=title)
  
  return(list(plot=g, outliers=outliers, R2=R2))
}


############################################################
## UI

ui <- fluidPage(
  titlePanel("Expression across 28 day cycle"),
  
  sidebarLayout(
    sidebarPanel(
      
      # Static or interactive plot
      radioButtons("plot_type",
                   label = h4("Plot type:"),
                   choices = list("Static (ggplot2)" = 1,
                                  "Interactive (plotly)" = 2),
                   selected = 1),
      
      # Advanced options
      # h4("Advanced Options:")
      checkboxInput("advanced",
                    label = "Show advanced options",
                    value = FALSE),
      
      # Colour
      conditionalPanel(
        condition="input.advanced == true", {
        selectInput("point_color",
                    label = "Point colours:",
                    choices = list("Endo status" = "endo",
                                   "Study" = "study",
                                   "AFS score (log)" = "afs_score_log"),
                    selected = 1)
      }),
      
      # Smoothing
      conditionalPanel(
        condition="input.advanced == true",
        sliderInput("loess_span",
                    label = "LOESS alpha (smoothing):",
                    min=0.2,
                    max=2,
                    step=0.05,
                    value=0.4)
      ),
      
      # Curve band size
      conditionalPanel(
        condition="input.advanced == true",
        sliderInput("band_size",
                    label = "Band size (std dev):",
                    min=0,
                    max=3,
                    step=0.1,
                    value=2)
      ),
      
      # Jitter to avoid overplotting
      conditionalPanel(
        condition="input.advanced == true",
        sliderInput("jitter_scale",
                    label = "Cycle day jitter:",
                    min=0,
                    max=2,
                    step=0.1,
                    value=0.5)
      ),
      
      # Extend days
      conditionalPanel(
        condition="input.advanced == true",
        sliderInput("extend_days",
                    label="Extend days forward and backwards:",
                    min=0,
                    max=28,
                    step=0.5,
                    value=10)
      ),
      
      # Choose input type
      radioButtons("input_type",
                   label = h4("Input type:"),
                   choices = list("Probe ID from list" = 1,
                                  "Enter probe ID" = 2),
                   selected = 2),
      
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
      
      # Submit button
      actionButton("submit", 
                   label="Sumbit")
      
    ),
    mainPanel(
      # Plot
      conditionalPanel(
        condition="input.plot_type == 1",
        plotOutput("ggplot")
      ),
      
      conditionalPanel(
        condition="input.plot_type == 2",
        plotlyOutput("plotly")
      ),
      
      hr()
      
    )
    
  ),
  
  hr(),
  
  # Outlier samples
  h3("Outlier samples:"),
  tableOutput("outlier_samples"),
  
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

server <- function(input, output){
  
  # Set reactive values
  rv <- reactiveValues(
    probe_name = probe_list[[1]],
    model = NA
  )
  
  # Update probe name when submit button is pressed
  observeEvent({input$submit; input$probe_select}, {
    if (input$input_type == "1") {
      rv$probe_name <- input$probe_select
    } else {
      rv$probe_name <- input$probe_text
    }
  })
  
  # Update model when arguments update
  observeEvent({rv$probe_name; input$plot_type; input$loess_span; 
    input$jitter_scale; input$band_size; input$point_color; input$extend_days}, {
      if (input$plot_type == "1") {
        rv$model <- fit_loess_model(exprs=combat_exprs, pheno=phenotype, probe=rv$probe_name,
                                    loess_span=input$loess_span, jitter_scale=input$jitter_scale,
                                    band_size=input$band_size, color_str=input$point_color,
                                    point_size=1, extend_days=input$extend_days)
      } else {
        # Use smaller point_size with plotly output
        rv$model <- fit_loess_model(exprs=combat_exprs, pheno=phenotype, probe=rv$probe_name,
                                    loess_span=input$loess_span, jitter_scale=input$jitter_scale,
                                    band_size=input$band_size, color_str=input$point_color,
                                    point_size=0.5, extend_days=input$extend_days)
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
  
  # Sample outliers table
  output$outlier_samples <- renderTable({
    rv$model$outliers
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
