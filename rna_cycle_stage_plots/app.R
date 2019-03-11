############################################################
## Expression along cycle stage
############################################################

# Plot expression a given gene along cycle stage

library(shiny)
library(ggplot2)
library(dplyr)
library(DT)
library(stringr)
library(reshape2)
library(edgeR)

qld_counts <- readRDS("data/qld_counts_list.rds")
sample_info <- readRDS("data/sample_info_2019-02-08.rds")
all_annotations <- readRDS("data/all_annotations.rds")

# Remove samples with cycle_stage 8,9,10 (note that URS samples don't have cycle stage)
phenotype_df <- sample_info %>% 
  filter(cycle_stage %in% 1:7) %>%
  mutate(cycle_stage=factor(cycle_stage),
         endo=factor(ifelse(afs_score > 0, 1, 0)))

# Remove replicate samples (those with names ending in _2)
# replicated <- phenotype_df %>% filter(str_detect(sample_id, "_"))
# stopifnot(replicated$sample_base_id %in% phenotype_df$sample_id)
phenotype_df <- phenotype_df %>% 
  filter(! str_detect(sample_id, "_2$"))

# Get counts
counts <- qld_counts$counts[,phenotype_df[,"sample_id"]]

# Filter counts
# ok <- filterByExpr(counts)
ok <- rowMeans(cpm(counts) > 0.5) > 0.5
counts <- counts[ok,]
gene_info <- merge(qld_counts$gene_info[ok,], 
             all_annotations %>% select(-strand), 
             by.x="ensembl_id", by.y="ensembl_gene_id") %>%
  arrange(ensembl_id)

# Get normalised counts
y <- DGEList(counts)
y <- calcNormFactors(y)
nc <- cpm(y, normalized.lib.sizes=TRUE, log=TRUE)

expression_cycle <- function(exprs, pheno, gene_id, title="", subtitle="") {
  dat <- exprs[gene_id,,drop=FALSE] %>% melt %>% 
    S4Vectors::rename(Var2="sample_id")
  stopifnot(dat[,"sample_id"] == pheno[,"sample_id"])
  dat <- cbind(dat, cycle=pheno[,"cycle_stage"], endo=pheno[,"endo"])
  dat_means <- dat %>% group_by(cycle) %>% summarise(mean=mean(value))
  g <- ggplot(dat, aes(x=cycle, y=value)) + 
    geom_jitter(aes(color=endo), height=0, width=0.1) + 
    geom_point(data=dat_means, aes(x=cycle, y=mean), col="orange", size=8, alpha=0.5) +
    geom_text(data=dat_means, size=3,
              aes(x=cycle, y=mean, label=sprintf("                   %0.2f", mean))) +
    labs(title=title, subtitle=subtitle) +
    theme_bw()
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
                   choices = list("Select from table" = 1,
                                  "Ensembl gene ID" = 2,
                                  "Gene symbol" = 3),
                   selected = 1),
    
      # Toggle input type for drop down menu
      conditionalPanel(
        condition="input.input_type == 2",
        textInput("gene_id", 
                  label = h4("Ensembl gene ID"), 
                  value = "ENSG00000082175")
      ),
      
      # Toggle input type for probe text input
      conditionalPanel(
        condition="input.input_type == 3",
        textInput("gene_symbol", 
                  label = h4("Gene symbol"), 
                  value = "PGR")
      ),
    
      # Submit button
      conditionalPanel(
        condition="input.input_type != 1",
        actionButton("submit", 
                     label="Sumbit")
      )
      
    ),
    
    mainPanel(
      # Plot
      plotOutput("plot")
    )
    
  ),
   
  hr(),

  # Probe search table
  dataTableOutput("gene_search_table"),
  
  # Search gene symbol
  textInput("gene_search",
            label = h4("Search for a gene name with regular expressions:"),
            value = ""),
  
  hr()
)

############################################################
## Server

server <- function(input, output){
  
  # Set reactive values
  rv <- reactiveValues(
    gene_id = "ENSG00000082175"
  )
  
  # Update probe name when submit button is pressed
  observeEvent({input$submit; input$input_type}, {
    if (input$input_type == "1") {
      rv$gene_id <- input$gene_id
    } else {
      rv$gene_id <- qld_counts$gene_info %>% filter(symbol == input$gene_symbol) %>% .$ensembl_id
    }
  })
  
  # Expression plot across cycle stage
  output$plot <- renderPlot({
    print(rv$gene_id)
    gene_symbol <- qld_counts$gene_info %>% filter(ensembl_id == rv$gene_id) %>% .$symbol
    expression_cycle(exprs=nc, pheno=phenotype_df, gene_id=rv$gene_id,
                     title=gene_symbol, subtitle=rv$gene_id)
  })
  
  # Table of genes
  output$gene_search_table <- renderDataTable({
    datatable(gene_info %>% filter(str_detect(symbol, input$gene_search)),
              options=list(pageLength=5), selection="single")
  })
  
  # Update selected probe when row is clicked
  observeEvent({input$gene_search_table_cell_clicked}, {
    if (!is.null(input$gene_search_table_cell_clicked$row)) {
      message(paste0("Row clicked: ", input$gene_search_table_cell_clicked$row))
      rv$gene_id <- gene_info %>% filter(str_detect(symbol, input$gene_search)) %>%
        .[input$gene_search_table_cell_clicked$row,"ensembl_id"]
    }
  })
  
}

############################################################
## Run

shinyApp(ui=ui, server=server)
