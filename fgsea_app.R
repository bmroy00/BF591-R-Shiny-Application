library(shiny)
library(bslib)
library(ggplot2)
library(tidyverse)
library(colourpicker)
library(DT)
library(shinythemes)

# Define UI for application that draws a histogram
ui <- fluidPage(
  theme = shinytheme("simplex"),
  # Title page
  titlePanel("GSEA"),
  
  
  sidebarLayout(
    # All in the sidebar
    sidebarPanel(
      # Accept only csv files as input
      fileInput(inputId = "fgsea_file", "Input fgsea File:",
                accept=c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv")
      ),
      
      # Button to generate plots and table
      submitButton(text = 'Submit', icon = NULL)
    ),
    
    # Organize where our plots and table exist
    mainPanel(
      tabsetPanel(
        tabPanel("Top Results", 
                 sliderInput(inputId = 'num_paths', min = 1, max = 40,
                             label = 'Top Results by adjusted p-value',
                             value = 5,
                             step = 1),
                 plotOutput('nes_plot')),
        
        tabPanel("Table",
                 sliderInput(inputId = 'table_padj_slider', min = 1, max = 40,
                             label = '-log(P-value) to Filter Table',
                             value = 25,
                             step = 1),
                 radioButtons(inputId = 'nes_path_button', 'Select NES Pathway(s)',
                              choices = c("Positive",
                                          "Negative",
                                          "All"),
                              selected = 'All'),
                 downloadButton(outputId = "downloadTable", "Download Table"),
                 DTOutput('nes_table')),
        
        tabPanel("Plots", 
                 sliderInput(inputId = 'scatter_nes_slider', min = 1, max = 40,
                            label = '-log(P-value) for Plot Coloring',
                            value = 3,
                            step = 1),
                 plotOutput('nes_scatter_plot'))
      )
    )
    
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  options(shiny.maxRequestSize = 100*1024^2)
  #' load_Data
  load_data <- reactive({
    df <- read.csv(input$fgsea_file$datapath)
    
    return(df)
  })
  
  #' 
  
  plot_top_pathways <- function(fgsea_results, num_paths){
    # Sort by NES
    fgsea_results <- arrange(fgsea_results,padj)
    
    # Take top 10
    top_vals <- head(fgsea_results,num_paths)

    # Plot top ten and bottom ten NES
    bars <- ggplot(top_vals, aes(x = reorder(pathway,NES),y= NES, fill = NES > 0)) + 
      geom_bar(stat='identity') +
      labs(title = 'fgsea results for Gene Ontology gene set', 
           y = 'Normalized Enrichment Score (NES)') + 
      theme(axis.text = element_text(size = 5),
            axis.title = element_text(size = 5)) + coord_flip()
    
    return(bars)
  }
  
  #' Summary Table
  
  draw_fgsea_table <- function(dataf) {
    # Total number of samples
    table <- as.data.frame(dataf)
    table <- filter(table, -log10(padj) >= input$table_padj_slider)
    if (input$nes_path_button == 'Positive') {
      table <- filter(table, NES > 0)
      return(table)
    }
    if (input$nes_path_button == 'Negative') {
      table <- filter(table, NES < 0)
      return(table)
    }
    return(table)
  }
  
  
  nes_scatter_plot <- function(dataf) {
    plot <- ggplot(dataf, aes(x = NES, y = -log10(padj), 
                              col = -log10(padj) >= input$scatter_nes_slider)) +
      geom_point() +
      labs(x = "NES", y = "-log10(padj)", title = paste0("Scatter Plot of NES by -log10(padj)",
                                                         color = "Passes P-value Threshold")) +
      scale_color_manual(values = c("grey", "purple")) 
      theme_minimal()
    
    return(plot)
  }
  
  generate_table <- reactive ({
    table <- draw_fgsea_table(load_data())
    return(table)
  })
  
  output$nes_plot <- renderPlot({plot_top_pathways(load_data(), input$num_paths)})
  
  output$nes_table <- renderDT({draw_fgsea_table(load_data())})
  
  output$nes_scatter_plot <- renderPlot({nes_scatter_plot(load_data())})
  
  output$downloadTable <- downloadHandler(
    filename = function() {
      return("fgsea_table.csv")
    },
    content = function(file) {
      write.csv(generate_table(), file)
    }
  )
  
}

# Run the application
shinyApp(ui = ui, server = server)