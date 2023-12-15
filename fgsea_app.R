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
  titlePanel("Samples"),
  
  
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
                 sliderInput(inputId = 'num_paths', min = 10, max = 100,
                             label = 'Top Results by adjusted p-value',
                             value = 60,
                             step = 1),
                 plotOutput('nes_plot')),
        
        tabPanel("Table",
                 sliderInput(inputId = 'table_padj_slider', min = 10, max = 100,
                             label = 'Top Results by adjusted p-value',
                             value = 60,
                             step = 1),
                 radioButtons(inputId = 'nes_path_button', 'Select NES Pathway(s)',
                              choices = c("Positive",
                                          "Negative",
                                          "All"),
                              selected = 'Positive'),
                 downloadButton(inputId = "downloadTable", "Download Table"),
                 DTOutput('nes_table')),
        
        tabPanel("Plots", 
                 sliderInput(inputId = 'scatter_nes_slider', min = 10, max = 1-0,
                            label = 'Top Results by adjusted p-value',
                            value = 60,
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
    df <- read.csv(input$samples_file$datapath)
    
    return(df)
  })
  
  #' 
  
  plot_top_pathways <- function(fgsea_results, num_paths){
    # Sort by NES
    fgsea_results <- arrange(fgsea_results,NES)
    
    # Take top 10
    top_vals <- head(fgsea_results,num_paths)
    
    # Take bottom 10
    bottom_vals <- tail(fgsea_results,num_paths)
    
    # Add top and bottom together
    top_bottom <- add_row(top_vals,bottom_vals)
    
    # Plot top ten and bottom ten NES
    bars <- ggplot(top_bottom, aes(x = reorder(pathway,NES),y= NES, fill = NES > 0)) + 
      geom_bar(stat='identity') +
      labs(title = 'fgsea results for Hallmark MSigDB gene set', 
           y = 'Normalized Enrichment Score (NES)') + coord_flip()
    
    return(bars)
  }
  
  #' Summary Table
  
  draw_table <- function(dataf) {
    # Total number of samples
    table <- as.data.frame(dataf)
    
    return(table)
  }
  
  
  create_plot <- function(dataf, continuous, categorical) {
    plot <- ggplot(dataf, aes(x = !!sym(categorical), y = !!sym(continuous))) +
      geom_violin(fill = "skyblue", color = "black") +
      labs(x = categorical, y = continuous, title = paste0("Violin Plot of ", continuous,"by ", categorical)) +
      theme_minimal()
    
    return(plot)
  }
  output$nes_plot <- renderPlot({plot_top_pathways(load_data(), input$num_paths)})
  
  output$nes_table <- renderDT({draw_table(load_data())})
  
  output$nes_scatter_plot <- renderPlot({create_plot(load_data(),input$continuous,input$categorical)})
  
  
  
}

# Run the application
shinyApp(ui = ui, server = server)