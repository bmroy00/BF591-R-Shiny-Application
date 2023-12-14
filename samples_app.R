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
      fileInput(inputId = "samples_file", "Input Samples File:",
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
        tabPanel("Summary", tableOutput('summary')),
        tabPanel("Table",DTOutput('table')),
        tabPanel("Plot", 
                 # Button for category
                 radioButtons(inputId = 'categorical', 'Choose the category to group by',
                              choices = c("Diagnosis",
                                          "Vonsattel Grade"),
                              selected = 'Diagnosis'
                 ),
                 # Button for continuous variable 
                 radioButtons(inputId = 'continuous', 'Choose the column for the y-axis',
                              choices = c("Post Mortem Interval",
                                          "Age of Death",
                                          "RNA Integrity Number",
                                          "mRNA-seq Reads",
                                          "Age of Onset",
                                          "Duration",
                                          "CAG",
                                          "h-v striatal score",
                                          "h-v cortical score"),
                              selected = 'Age of Death'
                 ),
                 plotOutput('histograms'))
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
  
  #' Diagnostic plots
  
  create_summary <-
    function(dataf) {
      column_types <- sapply(dataf, typeof)
      third_col <- c()
      
      
      summary <- data.frame(
        'Column Name' = colnames(dataf),
        'Type' = column_types
        #'Mean (sd) or Distinct Values' = third_col
      )
      
      return(summary)
    }
  
  #' Summary Table
  
  draw_table <- function(dataf) {
    # Total number of samples
    table <- as.data.frame(dataf)
    
    return(table)
  }
  
 
  create_histogram <- function(dataf, continuous, categorical) {
    plot <- ggplot(dataf, aes(x = !!sym(continuous), fill = !!sym(categorical))) +
      geom_bar(position = "dodge", bins = 30, color = "black") +
      labs(title = "Histogram of Continuous Variable Grouped by Category",
           x = continuous, y = categorical) +
      theme_minimal()
    
    return(plot)
  }
  output$summary <- renderTable({create_summary(load_data())})
  
  output$table <- renderDT({draw_table(load_data())})
  
  output$histograms <- renderPlot({create_histogram(load_data(),input$continuous,input$categorical)})
  
  
  
}

# Run the application
shinyApp(ui = ui, server = server)