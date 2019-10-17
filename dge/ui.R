library(shiny)
library(plotly)
library(DT)

shinyUI(fluidPage(
  
  # Application title
  titlePanel("Endometrial RNA-Seq and Microarray DGE"),
  h4("Peter Rogers"),
  
  fluidRow(
    column(3, wellPanel(
      
      h4(em("Upload sample groups:")),
      
      fileInput("file1", "Choose CSV file:",
                multiple = FALSE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),

      h4(em("Analysis options:")),
      
      selectInput("cycle_select",
                  label = "Filter by cycle stages:",
                  choices = list("Use samples from all stages" = "all",
                                 "Select specific stages" = "checkbox"),
                  selected = "all"),
      
      conditionalPanel(
        condition="input.cycle_select == 'checkbox'",
        wellPanel(
          checkboxGroupInput("cycle_checkbox",
                             label = "Included cycle stages (WIP)",
                             choices = list("Stage 1" = "s1",
                                            "Stage 2" = "s2",
                                            "Stage 3" = "s3",
                                            "Stage 4" = "s4",
                                            "Stage 5" = "s5",
                                            "Stage 6" = "s6",
                                            "Stage 7" = "s7"),
                             selected = c("s1", "s2", "s3", "s4", "s5", "s6", "s7"))
        )
      ),
      
      
      selectInput("study_select",
                  label = "Filter by study:",
                  choices = list("Use samples from all studies" = "all",
                                 "Select specific study samples" = "checkbox"),
                  selected = "all"),
      
      conditionalPanel(
        condition="input.study_select == 'checkbox'",
        wellPanel(
          checkboxGroupInput("study_checkbox",
                             label = "Included study samples",
                             choices = list("Study 1" = "study_1",
                                            "Study 2" = "study_2",
                                            "HMB" = "HMB"),
                             selected = c("study_1", "study_2", "HMB"))
        )
      ),
      
      # selectInput("filter_select",
      #             label = "Filter counts (WIP):",
      #             choices = list("Normal" = "normal",
      #                            "Lenient" = "lenient"),
      #             selected = "normal"),
      
      br(),
      
      h4(em("Display options:")),
      selectInput("plot_type",
                  label = "Plot type:",
                  choices = list("Static" = "ggplot",
                                 "Interactive" = "plotly"),
                  selected = "ggplot"),

      selectInput("n_top_table",
                  label = "Number of genes to list:",
                  choices = list("100" = 100,
                                 "All" = Inf),
                  selected = Inf),
      br(),

      actionButton("refresh", label = "Refresh")

    )),

    column(9, 
    
      mainPanel(
        tabsetPanel(
          tabPanel("Sample list", 
                   br(),
                   tabsetPanel(
                     tabPanel("RNA-seq samples", br(),
                              wellPanel(h5(textOutput("rna_sample_list_text"))), 
                              dataTableOutput("rna_phenotype")),
                     tabPanel("Microarray samples", br(),
                              wellPanel(h5(textOutput("array_sample_list_text"))), 
                              dataTableOutput("array_phenotype"))
                   )
          ),
                   
          tabPanel("RNA-seq DGE", br(),
                   # Top table of results
                   dataTableOutput("rna_top_table"), br(),
                   
                   tabsetPanel(
                     tabPanel("Expression plot", br(),
                              conditionalPanel(
                                condition="input.plot_type == 'ggplot'",
                                plotOutput("ggplot_rna_exprs", height="500px")),
                              
                              conditionalPanel(
                                condition="input.plot_type == 'plotly'",
                                plotlyOutput("plotly_rna_exprs", height="500px"))
                              ),
                     tabPanel("Cycle stage plot", br(),
                              conditionalPanel(
                                condition="input.plot_type == 'ggplot'",
                                plotOutput("ggplot_rna_cycle", height="500px")),
                              
                              conditionalPanel(
                                condition="input.plot_type == 'plotly'",
                                plotlyOutput("plotly_rna_cycle", height="500px"))
                              ),
                     tabPanel("P-value histogram", br(),
                              plotOutput("ggplot_rna_pval", height="500px"))
                   )
          ),
          
          tabPanel("Microarray DGE", br(),
                   # Top table of results
                   dataTableOutput("array_top_table"), br(),
                   
                   tabsetPanel(
                     tabPanel("Expression plot", br(),
                              conditionalPanel(
                                condition="input.plot_type == 'ggplot'",
                                plotOutput("ggplot_array_exprs", height="500px")),
                              
                              conditionalPanel(
                                condition="input.plot_type == 'plotly'",
                                plotlyOutput("plotly_array_exprs", height="500px"))
                     ),
                     tabPanel("Cycle stage plot", br(),
                              conditionalPanel(
                                condition="input.plot_type == 'ggplot'",
                                plotOutput("ggplot_array_cycle", height="500px")),
                              
                              conditionalPanel(
                                condition="input.plot_type == 'plotly'",
                                plotlyOutput("plotly_array_cycle", height="500px"))
                     ),
                     tabPanel("P-value histogram", br(),
                              plotOutput("ggplot_array_pval", height="500px"))
                   )
          ),
          
          tabPanel("...", br(), h4("Notes:"),
                   p("This app is still under development. Some of the options on the left panel don't work yet."),
                   p("The uploaded CSV file should have two columns. The first column should be sample ID and the second column is group name. There should be two groups only for A vs B analysis, or a numerical column for linear analysis."),
                   p("DGE analysis will not run if fewer than three samples per gorup."),
                   p("Currently, the DGE model doesn't include age. #TODO: add option to include age in model."),
                   p("'Expression plot' plots the expression of cycle-corrected data for the highlighted genes for the defined groups."),
                   p("'Cycle stage plot' plots the batch-corrected data across cycle stage."),
                   p("Currently the RNA data is normalised by 7-stage cycle predicted from the molecular data and the microarray data is normalised by 28-day cycle predicted from the molecular data."),
                   p("Filtering is lenient at the moment. Currently not filtering genes based on which samples are included in analysis. The current method is CPM > 0.5 for at least 20% of all samples in RNA-seq, and detection p-value < 0.05 in at least 20% of all samples in microarray."),
                   p(""))
        )
      )
      
    
    ))
))
