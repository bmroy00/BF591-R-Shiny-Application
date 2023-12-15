library(shiny)
library(bslib)
library(ggplot2)
library(tidyverse)
library(colourpicker)
library(DT)

# Define UI for application that draws a histogram
ui <- fluidPage(
  # Title page
  titlePanel("BF591 Final Project"),
  
  
  sidebarLayout(
    # All in the sidebar
    sidebarPanel(
      # Accept only csv files as input
      fileInput(inputId = "de_file", "Input Differential Expression File:",
                accept=c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv")
      ),
      submitButton(text = 'Submit', icon = NULL)
      ),
    
    # Organize where our plots and table exist
    mainPanel(
      tabsetPanel(
        tabPanel("Table", DTOutput('table')),
        tabPanel("Plot", plotOutput('volcano_plot'))
      )
    )
    
  )
)


# Define server logic required to draw a histogram
server <- function(input, output, session) {
  options(shiny.maxRequestSize = 100*1024^2)
  #' load_Data
  load_data <- reactive({
    df <- read.csv(input$de_file$datapath)

    return(df)
  })
  
  #' Diagnostic plots
  
  volcano_plot <-
    function(dataf) {
      
      dataf <- mutate(
        dataf, 
        volc_plot_status = case_when(log2FoldChange > 0 & padj < .1 ~ "UP",
                                     log2FoldChange < 0 & padj < .1  ~ "DOWN",
                                     padj > .1 ~ "NS"
        ), .after=1)
      
      dataf <- na.omit(dataf)
      
      # Volcano plot of padj by log2FoldChange 
      volcano <- ggplot(dataf, aes(x=log2FoldChange, 
                                             y=-log10(padj),
                                             col = volc_plot_status)) +
        geom_point() +
        labs(title = 'Volcano plot of DESeq2 differential expression results for Huntingtons Data',
             x = "log2FoldChange", y = "-log10(padj)") +
        guides(color = guide_legend(title = "Expressed:"))
        
      return(volcano)
    }
  
  #' Summary Table
  
  draw_table <- function(dataf) {
    table <- as.data.frame(dataf)
    
    return(table)
  }
  #' These outputs aren't really functions, so they don't get a full skeleton, 
  #' but use the renderPlot() and renderTabel() functions to return() a plot 
  #' or table object, and those will be displayed in your application.
  # Same here, just return the table as you want to see it in the web page
  output$table <- renderDT({draw_table(load_data())})
  output$volcano_plot <- renderPlot({volcano_plot(load_data())})
  

  
}

# Run the application
shinyApp(ui = ui, server = server)