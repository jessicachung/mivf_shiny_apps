############################################################
## Expression along 28 day cycle
############################################################

# Plot expression a given probe along cycle days

## TODO: 
## Test other models
## Change plots to reactive

library(shiny)
library(ggplot2)
library(dplyr)
library(stringr)
library(reshape2)
library(plotly)
library(stats)
#library(splines) # splines
#library(MASS)    # lm.ridge
#library(glmnet)  # lasso (alpha=1), ridge (alpha=0)

# setwd("~/Work/2016_mivf/shiny_apps/cycle_day_plots")
load("data/combat_exprs.rda")
phenotype <- read.table("data/combat_phenotype_with_day_cycle.tsv", header=TRUE, sep="\t", 
                           stringsAsFactors=FALSE) %>%
  mutate(endo=factor(endo),
         afs_score_log=log2(afs_score+0.01))
probe_df <- read.table("data/illumina_v4_annotation_with_detection.tsv", header=TRUE, sep="\t", 
                       stringsAsFactors=FALSE) %>% 
  S4Vectors::rename(X10.="10%", X50.="50%", X90.="90%")

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

combat_exprs <- combat_exprs[,phenotype[,"sample_id"]]

cycle <- phenotype[,"day_cycle"]

expression_cycle <- function(exprs, pheno, probe, loess_span=0.75, jitter_scale=1, 
                             band_size=2, color_str="endo", point_size=1, title="",
                             cycle_data=TRUE) {
  # Melt data
  dat <- exprs[probe,,drop=FALSE] %>% melt %>% 
    S4Vectors::rename(Var2="sample_id")
  stopifnot(dat[,"sample_id"] == pheno[,"sample_id"])
  dat <- cbind(dat, pheno %>% select(-sample_id))
  sd <- sd(dat[,"value"])
  
  if (cycle_data) {
    # Repeat data for cyclic nature
    cyclic_dat <- rbind(mutate(dat, day_cycle=day_cycle - 28),
                        dat,
                        mutate(dat, day_cycle=day_cycle + 28))
    
    # Fit loess model
    fit <- loess(value ~ day_cycle, cyclic_dat, span=loess_span, model=TRUE)
    cycle_range <- seq(0,28,by=0.5)
  } else {
    # Fit loess model
    fit <- loess(value ~ day_cycle, dat, span=loess_span, model=TRUE)
    cycle_range <- seq(min(cycle),max(cycle),by=0.5)
  }
  
  # Get curve of loess model
  predict <- predict(fit, cycle_range)
  curve <- data.frame(day_cycle=cycle_range, loess=predict, 
                      ymin=predict-band_size*sd, ymax=predict+band_size*sd)
  
  # Get outliers
  residuals <- fit$residuals[between(fit$x, 0.1, 28)]
  names(residuals) <- dat[,"sample_id"]
  indices <- which(abs(residuals) - band_size * sd > 0)
  outliers <- pheno[indices,] %>% arrange(day_cycle) %>% dplyr::select(-text)
  
  # Jitter
  if (jitter_scale > 0) {
    set.seed(0)
    jitter <- runif(nrow(dat),jitter_scale * -1,jitter_scale)
    dat <- mutate(dat, day_cycle=day_cycle+jitter)
  }
  
  g <- ggplot(dat, aes(x=day_cycle, y=value)) +
    geom_point(aes_string(text="text", color=color_str), size=point_size) + 
    geom_line(data=curve, aes(x=day_cycle, y=loess), alpha=0.3, size=1.5) +
    geom_ribbon(data=curve, aes(x=day_cycle, y=loess, ymin=ymin, ymax=ymax), alpha=0.1) +
    labs(title=title)
  
  return(list(plot=g, outliers=outliers, cycle_means=NA))
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
      
      # Cycle data or not
      conditionalPanel(
        condition="input.advanced == true",
        checkboxInput("cycle_data",
                      label = "Loop 28 day data",
                      value=TRUE)
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
            value = "e.g. VEZT"),
  
  # Probe search table
  tableOutput("gene_search_table")
  
)

############################################################
## Server

server <- function(input, output){
  
  # Set reactive values
  rv <- reactiveValues(
    probe_name = probe_list[[1]]
  )
  
  # Update probe name when submit button is pressed
  observeEvent({input$submit; input$probe_select}, {
    if (input$input_type == "1") {
      rv$probe_name <- input$probe_select
    } else {
      rv$probe_name <- input$probe_text
    }
  })
  
  output$ggplot <- renderPlot({
    print(rv$probe_name)
    expression_cycle(exprs=combat_exprs, pheno=phenotype, probe=rv$probe_name,
                     loess_span=input$loess_span, jitter_scale=input$jitter_scale,
                     band_size=input$band_size, color_str=input$point_color,
                     cycle_data=input$cycle_data,
                     title=rv$probe_name)
  })
  
  output$plotly <- renderPlotly({
    expression_cycle(exprs=combat_exprs, pheno=phenotype, probe=rv$probe_name,
                     loess_span=input$loess_span, jitter_scale=input$jitter_scale,
                     band_size=input$band_size, color_str=input$point_color, 
                     point_size=0.5, cycle_data=input$cycle_data,
                     title=rv$probe_name)$plot
  })
  
  output$outlier_samples <- renderTable({
    ### WIP, change to reactive
    expression_cycle(exprs=combat_exprs, pheno=phenotype, probe=rv$probe_name,
                     loess_span=input$loess_span, jitter_scale=input$jitter_scale,
                     band_size=input$band_size, color_str=input$point_color, 
                     point_size=0.5, cycle_data=input$cycle_data,
                     title=rv$probe_name)$outliers %>% dplyr::select(-afs_score_log)
  })
  
  # Probe info table
  output$probe_info <- renderTable({
    probe_df %>% filter(IlluminaID == rv$probe_name) %>%
      dplyr::select(IlluminaID, SymbolReannotated, GenomicLocation, ProbeQuality, MeanDetectionPVal:`90%`) %>%
      S4Vectors::rename(SymbolReannotated="Symbol")
  })
  
  # Search table
  output$gene_search_table <- renderTable({
    probe_df %>% filter(SymbolReannotated == input$gene_search) %>%
      dplyr::select(IlluminaID, SymbolReannotated, GenomicLocation, ProbeQuality, MeanDetectionPVal:`90%`) %>%
      S4Vectors::rename(SymbolReannotated="Symbol")
  })
  
}

############################################################
## Run

shinyApp(ui=ui, server=server)
