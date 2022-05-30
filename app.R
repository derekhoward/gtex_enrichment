library(shiny)
library(ggplot2)
library(readr)
library(magrittr)
library(plotly)
library(tibble)
library(tidyr)
library(dplyr)
library(pheatmap)
library(purrr)
#library(shinyjs)
source("./AUCFunction.R")
source("./string_processing.R")

apply_MWU <- function(column, targetIndices) {
  wilcox.test(column[targetIndices], column[!targetIndices], conf.int = F)$p.value
}

ui <- fluidPage(
  shinyjs::useShinyjs(),
  #tags$head(includeHTML("google-analytics.html")),
  # App title ----
  titlePanel("Polygenic tester for GTEx data"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    # Sidebar panel for inputs ----
    sidebarPanel(
      textAreaInput(
        inputId = "genelist",
        label = "Input your gene list:",
        value = 'THRA\nRTN1\nTUBA1A\nSTMN2\nCRMP1\nTUBB3\nISLR2',
        rows = 10
      ),
      textAreaInput(
        inputId = "background_genelist",
        label = "Background gene list (optional):",
        value = '',
        rows = 2
      ),
      selectInput(
        inputId = 'species',
        label = 'Species of input genes:',
        choices = c('Human', 'Mouse', 'Rhesus Macaque')
      ),
      actionButton(inputId = "submit",
                   label = "Submit"),
      br(),
      br(),
      downloadButton(outputId = "download_data", label = "Download results as .csv"),
      downloadButton(outputId = "download_heatmap", label = "Download heatmap"),
      hr(),
      tags$b("Data was made available by the Genotype-Tissue Expression (GTEx) project and is available from: "),
      #br(),
      tags$a(href="https://gtexportal.org/home/datasets", "GTEx Portal data"),
      hr(),
      tags$a(href="https://github.com/derekhoward/gtex_enrichment", "Source code")
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      div(
        id = "main",
        # Output: Verbatim text for data summary ----
        verbatimTextOutput("summary"),
        br(),
        dataTableOutput("view"),
        br(),
        #plotlyOutput("dotplot"),
        verbatimTextOutput("info")
      )
    )
  )
)


# Define server logic process and output top cortical layers/zones ----
server <- function(input, output) {
  output$summary <- renderPrint({
    cat("\nResults will load here when complete")
    cat("\n")
    #print(gc())
    #print(Sys.info()['nodename'])
  })
  
  observeEvent(input$submit, {
    start <- Sys.time()
    
    cleaned_gene_list <- isolate(process_input_genes(input$genelist))
    background_cleaned_gene_list <- isolate(process_input_genes(input$background_genelist))
    
    
    gtex_expression_ranks <- read_csv(file = './data/processed/gtex_processed_ranks.csv')
    
    cleaned_gene_list <- convert2human(input_genes = cleaned_gene_list, in_species = input$species)
    background_cleaned_gene_list <- convert2human(input_genes = background_cleaned_gene_list, in_species = input$species)
    
    first_cleaned_gene_list <- cleaned_gene_list #for printing out the original size
    
    if (length(background_cleaned_gene_list) == 1 && background_cleaned_gene_list == "") {
      gene_universe <- read_table('./data/gene_universe.txt', col_names = F) %>% .$X1
    } else { #if given a background
      gene_universe <- background_cleaned_gene_list
    }
    gtex_expression_ranks %<>% filter(gene_symbol %in% gene_universe)
    cleaned_gene_list <- intersect(cleaned_gene_list, gene_universe)
    #re rank based on new background
    gtex_expression_ranks %<>% mutate_if(is.numeric, rank)
    
    print(paste0("Before time taken:", Sys.time() - start))
    
    #for indices - use dplyr for ease
    forIndices <- as_tibble(gtex_expression_ranks$gene_symbol)
    names(forIndices) <- 'gene_symbol'
    forIndices %<>% mutate(isTargetGene = gene_symbol %in% cleaned_gene_list)
    targetIndices <- forIndices$isTargetGene
    
    # only columns from cortical zones remain in df
    df <- gtex_expression_ranks %>%
      select(-gene_symbol)
    
    AUROC <- map_df(df, auroc_analytic, targetIndices)
    wilcox_tests <- map_df(df, apply_MWU, targetIndices)
    
    # group results together in a single table
    table <- bind_cols(gather(AUROC, key = Tissue, value = AUROC), 
                       gather(wilcox_tests, value = pValue)) %>%
      select(-key)
    
    print(paste0("Wilcox time taken:", Sys.time() - start))
    
    # these are the values for the results table
    table %<>% arrange(-AUROC)
    table %<>% mutate(pValue = signif(pValue, digits = 3), 
                      AUROC = signif(AUROC, digits = 3),
                      adjusted_P = signif(p.adjust(pValue, method = "bonferroni"), digits = 3))
    
    
    output$summary <- renderPrint({
      #count of intersection of submitted genes with total gene list
      cat(paste("Time taken:", round(Sys.time() - start), "seconds"))
      cat(paste(
        "\nGenes found in data:",
        length(cleaned_gene_list),
        "of",
        length(first_cleaned_gene_list)
      ))
      cat(paste(
        "\nBackground genes:", nrow(gtex_expression_ranks)
      ))
    })
    
    output$view <- renderDataTable(    table, escape = FALSE, options = list(
      paging =FALSE,
      pageLength =  54 
    ))
    
    #code duplication with HPA 
    output$download_heatmap <-
      downloadHandler(
        filename = "polygenic_GTEx_heatmap.pdf",
        content = function(file) {
          df <- as.data.frame(gtex_expression_ranks %>% filter(gene_symbol %in% cleaned_gene_list))
          rownames(df) <- df$gene_symbol
          df <- df[ , (names(df) != "gene_symbol"), drop = T]
          
          colnames(df) <- gsub("[|]", " ", colnames(df))
          
          plot_out <- pheatmap(df, main = paste0("GTEx heatmap for ", length(cleaned_gene_list), " genes\ncolor scale represents specific expression rank (higher is more specific)"))
          
          pdf(file, width=1+(23/170) * ncol(df), height=5+nrow(df) *.2)
          grid::grid.newpage()
          grid::grid.draw(plot_out$gtable)
          dev.off()
        }
      )
    output$download_data <-
      downloadHandler(
        filename = "polygenic_GTEx_results.csv",
        content = function(file) {
          write_csv(table, file)
        }
      )
    
  })
}

shinyApp(ui, server)