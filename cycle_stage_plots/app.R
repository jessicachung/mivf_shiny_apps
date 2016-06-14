############################################################
## Expression along cycle stage
############################################################

# Plot expression a given probe along cycle stage

library(shiny)
library(ggplot2)
library(dplyr)
library(stringr)
library(reshape2)

load("data/combat_exprs.rda")
phenotype_df <- read.table("data/combat_phenotype.tsv", header=TRUE, sep="\t", 
                           stringsAsFactors=FALSE)
probe_df <- read.table("data/illumina_v4_annotation.tsv", header=TRUE, sep="\t", 
                       stringsAsFactors=FALSE)

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

# Remove samples with cycle_stage 8,9,10
phenotype_df <- phenotype_df %>% 
  # filter(cycle_stage %in% 1:7) %>%
  mutate(cycle_stage=factor(cycle_stage))
combat_exprs <- combat_exprs[,phenotype_df[,"sample_id"]]

# TODO: Replace with real p-values later
probe_df <- probe_df %>% mutate(mean_pval=NA)


expression_cycle <- function(exprs, pheno, probe, title="") {
  dat <- exprs[probe,,drop=FALSE] %>% melt %>% 
    S4Vectors::rename(Var2="sample_id")
  stopifnot(dat[,"sample_id"] == pheno[,"sample_id"])
  dat <- cbind(dat, cycle=pheno[,"cycle_stage"])
  dat_means <- dat %>% group_by(cycle) %>% summarise(mean=mean(value))
  g <- ggplot(dat, aes(x=cycle, y=value, color=cycle)) + 
    geom_point() + 
    geom_point(data=dat_means, aes(x=cycle, y=mean), col="red", size=6, alpha=0.5) +
    geom_text(data=dat_means, size=3,
              aes(x=cycle, y=mean, label=sprintf("              %0.2f", mean))) +
    labs(title=title)
  return(list(plot=g, cycle_means=dat_means))
}


############################################################
## UI

ui <- fluidPage(
  titlePanel("Expression across cycle stage"),
  
  sidebarLayout(
    sidebarPanel(
      
      # Choose input type (radioButtons or drop down menu?)
      radioButtons("input_type",
                   label = h4("Input type:"),
                   choices = list("Probe ID from list" = 1,
                                  "Enter probe ID" = 2),
                   selected = 1),
      
      # Toggle input type for drop down menu
      conditionalPanel(
        condition="input.input_type == 1",
        radioButtons("probe_select", 
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
      plotOutput("plot"),
      
      hr(),
      
      # Probe info table
      h3("Selected probe:"),
      tableOutput("probe_info")
    )
    
  ),
  
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
  
  output$plot <- renderPlot({
    print(rv$probe_name)
    expression_cycle(exprs=combat_exprs, pheno=phenotype_df, probe=rv$probe_name,
                     title=rv$probe_name)
  })
  
  # Probe info table
  output$probe_info <- renderTable({
    probe_df %>% filter(IlluminaID == rv$probe_name) %>%
      dplyr::select(IlluminaID, SymbolReannotated, GenomicLocation, ProbeQuality, mean_pval) %>%
      S4Vectors::rename(SymbolReannotated="Symbol", mean_pval="MeanDetectionPValue")
  })
  
  # Search table
  output$gene_search_table <- renderTable({
    probe_df %>% filter(SymbolReannotated == input$gene_search) %>%
      dplyr::select(IlluminaID, SymbolReannotated, GenomicLocation, ProbeQuality, mean_pval) %>%
      S4Vectors::rename(SymbolReannotated="Symbol", mean_pval="MeanDetectionPValue")
  })
  
}

############################################################
## Run

shinyApp(ui=ui, server=server)
