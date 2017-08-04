library(shiny)
library(plotly)
library(DT)

shinyUI(fluidPage(
  
  # Application title
  titlePanel("MIVF BMI Models"),
  
  fluidRow(
    column(4, wellPanel(
      
      h4(em("Analysis options:")),
      
      selectInput("plot_type",
                  label = "Plot type:",
                  choices = list("Static" = "ggplot",
                                 "Interactive" = "plotly"),
                  selected = "ggplot"),
      
      selectInput("exprs_selection",
                  label = "Gene subset:",
                  choices = list("All good genes" = "good",
                                 "Genes with priors" = "prior")),
      
      selectInput("patient_selection",
                  label = "Patient subset:",
                  choices = list("All patients" = "all",
                                 "Endo patients" = "endo",
                                 "Control patients" = "control",
                                 "WOI patients (16-24)" = "woi",
                                 "Cycle stage 4" = "cycle4")),
      
      selectInput("n_top_table",
                  label = "Number of genes to list:",
                  choices = list("100" = 100,
                                 "All" = Inf),
                  selected = 100),
      
      selectInput("selection",
                  label = "Selection method:",
                  choices = list("Slider" = "slider",
                                 "Checkbox" = "checkbox"),
                  selected = "slider"),
      
      # Group 1 & 2 slider input
      conditionalPanel(
        condition="input.selection == 'slider'",
        sliderInput("slider_1",
                    label = "Group 1 BMI:",
                    min = 15,
                    max = 52,
                    value = c(15, 25)),
        sliderInput("slider_2",
                    label = "Group 2 BMI:",
                    min = 15,
                    max = 52,
                    value = c(25, 52))
      ),
      
      # Group 1 & 2 checkbox input
      conditionalPanel(
        condition="input.selection == 'checkbox'",
        checkboxGroupInput("checkbox_1",
                           label = "Group 1 BMI classes:",
                           choices = list("Underweight" = "Underweight",
                                          "Normal" = "Normal",
                                          "Pre-obese" = "Pre-obese",
                                          "Obese" = "Obese"),
                           selected = c("Underweight", "Normal")
        ),
        checkboxGroupInput("checkbox_2",
                           label = "Group 2 BMI classes:",
                           choices = list("Underweight" = "Underweight",
                                          "Normal" = "Normal",
                                          "Pre-obese" = "Pre-obese",
                                          "Obese" = "Obese"),
                           selected = c("Pre-obese", "Obese")
        )
      ),
      
      br(),
      textOutput("n_group_1"),
      textOutput("n_group_2"),
      br(),
      
      actionButton("refresh", label = "Refresh")
      
    )),
    
    column(8, wellPanel(
      # Top table of results
      dataTableOutput("top_table")
    ),
    
    conditionalPanel(
      condition="input.plot_type == 'ggplot'",
      plotOutput("ggplot", height="500px")),
    
    conditionalPanel(
      condition="input.plot_type == 'plotly'",
      plotlyOutput("plotly", height="500px"))
    
    ))
))
