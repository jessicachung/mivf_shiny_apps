############################################################
## Boxplot of expression for SNP genotypes
############################################################

# Plot expression for different genotypes for a given SNP, for any given probe
# Notes:
#  - Removed samples not in cycle stage 1-7
#  - Removed probes weren't labelled good or perfect
#  - Expression values are actually residuals after fitting the expression to a
#    linear model with cycle stage as covariates
#  - Original study samples (analysis 1): n=121
#  - Extended set of endometriosis patients (analysis 2): n=155

library(shiny)
library(ggplot2)
library(dplyr)

load("data/data.rda")

plot_expressions <- function(y, snp, probe, geno_df) {
  samples <- colnames(y)
  m <- match(samples, geno_df %>% .$Sample_ID)
  genotype <- as.data.frame(geno_df)[m,snp]
  dat <- data.frame(genotype=genotype, exprs=y[probe,]) %>% filter(genotype != "NA")
  g <- ggplot(dat, aes(x=genotype, y=exprs)) +
    geom_boxplot() + geom_jitter(width=0.3) +
    labs(title=paste0(snp, "\n", probe))
  return(g)
}

# Make snp list for drop down box
snp_list <- list()
for (i in seq_len(nrow(snp_locs))) {
  snp_name <- sprintf("%s (%s:%s)", snp_locs[i, "snp"], snp_locs[i, "chr"], 
                      as.character(snp_locs[i, "loc"]))
  snp_list[[snp_name]] <- i
}

############################################################
## UI

ui <- fluidPage(
  titlePanel("Expression boxplots for SNP genotypes"),
  
  sidebarLayout(
    sidebarPanel(
      
      # Samples (original or all)
      radioButtons("sample_set", 
                   label = h4("Samples to include:"),
                   choices = list("Original samples (n=121)" = 1, 
                                  "All study 1 samples (n=155)" = 2), 
                   selected = 1),
      
      # Text 1
      selectInput("snp", 
                  label = h4("SNP:"),
                  choices = snp_list,
                  selected = 38),
      
      # Probe
      textInput("probe", 
                label = h4("Illumina probe ID"), 
                value = "ILMN_1774828"),
      
      # Submit button
      actionButton("submit", 
                   label="Sumbit")
    ),
    mainPanel(
      # Boxplot
      plotOutput("plot")
    )
    
  ),
  
  hr(),
  
  # Probe info table
  h3("Selected probe:"),
  tableOutput("probe_info"),
  
  # SNP info table
  h3("Selected SNP:"),
  tableOutput("snp_info"),
  
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
    probe_name = "ILMN_1774828"
  )
  
  
  # Update probe name when submit button is pressed
  observeEvent(input$submit, {
    rv$probe_name <- input$probe
  })
  
  output$plot <- renderPlot({
    # Convert input snp number to name
    snp_name <- snp_locs[as.numeric(input$snp),"snp"] %>% as.character()
    print(snp_name)
    print(rv$probe_name)
    
    if (input$sample_set == "1") {
      # If analysis 1
      plot_expressions(y=residuals_1, snp=snp_name, 
                            probe=rv$probe_name, geno_df=genotype_df)
    } else {
      # If analysis 2
      plot_expressions(y=residuals_2, snp=snp_name, 
                            probe=rv$probe_name, geno_df=genotype_df)
    }
  })
  
  # SNP info table
  output$snp_info <- renderTable({
    snp_locs[as.numeric(input$snp),]
  })
  
  # Probe info table
  output$probe_info <- renderTable({
    good_probes %>% filter(IlluminaID == rv$probe_name) %>%
      dplyr::select(IlluminaID, SymbolReannotated, GenomicLocation, ProbeQuality, mean_pval) %>%
      S4Vectors::rename(SymbolReannotated="Symbol", mean_pval="MeanDetectionPValue")
  })
  
  # Search table
  output$gene_search_table <- renderTable({
    good_probes %>% filter(SymbolReannotated == input$gene_search) %>%
      dplyr::select(IlluminaID, SymbolReannotated, GenomicLocation, ProbeQuality, mean_pval) %>%
      S4Vectors::rename(SymbolReannotated="Symbol", mean_pval="MeanDetectionPValue")
  })
  
}

############################################################
## Run

shinyApp(ui=ui, server=server)
