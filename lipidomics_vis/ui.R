library(shiny)
library(plotly)
library(DT)

shinyUI(fluidPage(
  
  titlePanel("MIVF Lipidomics Visualisation"),
  
  fluidRow(
    column(4, wellPanel(
      
      h4(em("Visualisation options:")),
      
      selectInput("plot_type",
                  label = "Plot type:",
                  choices = list("Static" = "ggplot",
                                 "Interactive" = "plotly"),
                  selected = "ggplot"),
      
      
      selectInput("exprs_dataset",
                  label = "Dataset:",
                  choices = list("Raw data" = "raw",
                                 "Quantile normalised (+imputed)" = "qnorm",
                                 "Combat normalised (batch + cycle)" = "combat"),
                  selected = "qnorm"),
      
      selectInput("lipid_subset",
                  label = "Lipids to display in table:",
                  choices = list("All lipids (can be slow!)" = "all",
                                 "Random 100" = "random_100",
                                 "Random 1000" = "random_1000"),
                  selected = "random_100"),
      
      selectInput("sample_label",
                  label = "Sample ID label (with interactive plots):",
                  choices = list("QMIR (X210###)" = "sample_id",
                                 "RMH ID" = "index",
                                 "Mass spec ID (P##)" = "ms_sample_id"),
                  selected = "sample_id"),
      
      selectInput("plot_1_x_axis",
                  label = "X-axis:",
                  choices = list("Run number" = "run_number",
                                 "Endo group" = "endo_group",
                                 "Cycle stage" = "path_menst_cycle"),
                  selected = "run"),
      
      selectInput("plot_1_colour",
                  label = "Colour:",
                  choices = list("Endo group" = "endo_group",
                                 "Cycle stage" = "path_menst_cycle",
                                 "Batch" = "batch",
                                 "Sample type" = "sample_type"),
                  selected = "endo")
      
      
      )
      
    ),
    
    column(8, wellPanel(
      # Feature dataframe
      dataTableOutput("feature_table")
    )
    )),
  
  column(8, 
         conditionalPanel(
           condition="input.plot_type == 'ggplot'",
           plotOutput("plot_1_ggplot", height="500px")),
         
         conditionalPanel(
           condition="input.plot_type == 'plotly'",
           plotlyOutput("plot_1_plotly", height="500px"))
  ),
  column(4, 
         plotOutput("plot_2_ggplot", height="500px")
  )
  
))
