# The server.R file creates the server side of the Shiny App.
# This Shiny App was written by Kasun Buddika to distribute data from Buddika et al., 2021.

## Load required libraries
library(DT)
library(ggiraph)
library(tidyverse)

## RNA-seq data import and pre-processing
patr1_data <- read_csv("https://github.com/jkkbuddika/buddika_2021_analysis/raw/main/shinyApp/data/DESeq2_Results_Patr-1.csv")
esg_data <- read_csv("https://github.com/jkkbuddika/buddika_2021_analysis/raw/main/shinyApp/data/DESeq2_Results_esg.csv")
rad21_data <- read_csv("https://github.com/jkkbuddika/buddika_2021_analysis/raw/main/shinyApp/data/DESeq2_Results_rad21.csv")

dataSets <- list(patr1_data, esg_data, rad21_data)
names(dataSets) <- c("Patr-1 Data", "esg Data", "Rad21 Data")

## Screen data import and pre-processing
screen_data <- read_csv("https://github.com/jkkbuddika/buddika_2021_analysis/raw/main/shinyApp/data/Screen_data.csv")

function(input, output, session) {
  
  ## Enable switching between datasets
  updateSelectInput(session,
                    "select_dataset",
                    choices = names(dataSets)
  )
  
  updateSelectInput(session,
                    "select_dataset_plot",
                    choices = names(dataSets)
  )
  
  ## Enable downloading data tables
  output$download_datatables <- downloadHandler(
    filename = function() {
      paste0(Sys.Date(), "_dge_data.csv")
    },
    content = function(file) {
      dataSets[[input$select_dataset]] %>%
        write_csv(file)
    }
  )
  
  output$download_screenData <- downloadHandler(
    filename = function() {
      paste0(Sys.Date(), "_screen_data.csv")
    },
    content = function(file) {
      screen_data %>%
        write_csv(file)
    }
  )
  
  ## Render data tables with the DT package
  output$dgeData <- renderDT({
    dataSets[[input$select_dataset]] %>% 
      select(Symbol:padj) %>%
      arrange(!desc(padj)) %>%
      datatable(rownames = FALSE)
  })
  
  ## Render volcano plot visualization
  output$plot <- renderGirafe({
    
    ## To remove single quotation marks from symbols (if any)
    dataSets[[input$select_dataset_plot]]$Symbol <- str_replace_all(
      dataSets[[input$select_dataset_plot]]$Symbol, "'", "")
    
    ## ggplot object
    gg_point <- dataSets[[input$select_dataset_plot]] %>%
      mutate(Significance = ifelse(log2FoldChange > 1 & padj < 0.05, "Up",
                                   ifelse(log2FoldChange < -1 & padj < 0.05, "Down", "NS"))) %>%
      
      ggplot(aes(x = log2FoldChange, y = -log(padj), color = Significance)) +
      geom_point_interactive(aes(data_id = Symbol, tooltip = Symbol), 
                             size = 1.5, alpha = 0.6) +
      scale_color_manual(values = c("#008B8B", "#808080", "#DC143C")) +
      geom_hline(yintercept = -log(0.05), linetype = "dashed") +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
      theme_light() +
      theme(panel.grid = element_blank())
    
    ## Produce a graphic using girafe
    intPlot <- girafe(code = print(gg_point),
                options = list(
                  opts_hover(css = "fill:#800080;stroke:black;cursor:pointer;", reactive = TRUE),
                  opts_selection(
                    type = "multiple", css = "fill:#FF3333;stroke:black;")
                ))
    intPlot
  })
  
  output$screenData <- renderDT({
    screen_data %>% 
      arrange(!desc(Symbol)) %>%
      datatable(rownames = FALSE)
  })
  
}