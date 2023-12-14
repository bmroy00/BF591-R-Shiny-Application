library(shiny)
library(bslib)
library(ggplot2)
library(tidyverse)
library(colourpicker)

# Define UI for application that draws a histogram
ui <- fluidPage(
  # Title page
  titlePanel("BF591 Final Project"),
  
  
  sidebarLayout(
    # All in the sidebar
    sidebarPanel(
      # Accept only csv files as input
      fileInput(inputId = "counts_file", "Input Counts File:",
                accept=c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv")
      ),
      
      
      
      # Slider for padj cutoff threshold
      sliderInput(inputId = 'var_slider', min = 0, max = 100,
                  label = 'Include Genes with X Percentile of Variance:',
                  value = 5,
                  step = 1),
      
      sliderInput(inputId = 'non_zero_slider', min = 0, max = 69,
                  label = 'Include Genes with at least X non-zero samples:',
                  value = 10,
                  step = 1),
      
      # Button to generate plots and table
      submitButton(text = 'Submit', icon = NULL)
    ),
    
    # Organize where our plots and table exist
    mainPanel(
      tabsetPanel(
        tabPanel("Table", tableOutput('table')),
        #tabPanel("Plots", plotOutput('diag_plot')),
        #tabPanel("Heatmap", plotOutput('heatmap'))
        tabPanel("PCA", plotOutput('pca'))
      )
    )
    
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  options(shiny.maxRequestSize = 100*1024^2)
  #' load_Data
  load_data <- reactive({
    df <- read.csv(input$counts_file$datapath)
    dataf <- as.data.frame(df)
    
    
    # Include genes with at least X percentile of variance
    var_slider <- input$var_slider
    non_zero_slider <- input$non_zero_slider
    
    variances <- apply(dataf[,-1], 1, var)
    
    percentile_variance <- sapply(variances, function(x) {
      sum(variances <= x) / length(variances)
    })
    
    var_filter <- which(percentile_variance >= var_slider/100)
    var_gene_names <- rownames(dataf)[var_filter]
    var_filtered_df <- dataf[var_gene_names, ]
    
    # Include genes with at least X non-zero samples
    non_zero_samples <- apply(var_filtered_df[,-1] != 0, 1, sum)
    genes_nonzero <- which(non_zero_samples >= non_zero_slider)
    genes_names_nonzero <- rownames(dataf)[genes_nonzero]
    filtered_counts <- dataf[genes_names_nonzero, ]
    
    return(list(val1 = dataf, val2 = filtered_counts))
  })
  
  #' Diagnostic plots

  diagnostic_plots <-
    function(dataf, filtered_counts) {
      medians <- apply(dataf[,-1], 1, median)
      variances <- apply(dataf[,-1], 1, var)
      zeros_count <- apply(dataf[,-1] != 0, 1, sum)
      
      percentile_variance <- sapply(variances, function(x) {
        sum(variances <= x) / length(variances)
      })
      
      var_filter <- which(percentile_variance >= var_slider/100)
      var_gene_names <- rownames(dataf)[var_filter]
      var_filtered_df <- dataf[var_gene_names, ]
      
      # Include genes with at least X non-zero samples
      non_zero_samples <- apply(var_filtered_df[,-1] != 0, 1, sum)
      genes_nonzero <- which(non_zero_samples >= non_zero_slider)
      genes_names_nonzero <- rownames(dataf)[genes_nonzero]
      filtered_counts <- dataf[genes_names_nonzero, ]
      
      plot <- ggplot(dataf, aes(x = medians, y = zeros_count, 
                        col = zeros_count %in% non_zero_samples)) +
        geom_point() +
        labs(x = "Median Counts", y = "Number of Zeros") +
        title("Median Counts vs. Number of Zeros") +
        scale_color_manual(values=c('purple','white')) +
        theme_dark()
      
      return(plot)
    }
  
  #' Summary Table

  draw_table <- function(dataf, filtered_counts) {
    # Total number of samples
    samp_count <- ncol(dataf[,-1]) 
    
    # Total number of genes
    gene_count <- nrow(dataf)
  
    genes_passing <- nrow(filtered_counts)
    
    table <- tibble(
      "Number of Samples" = samp_count,
      "Total Number of genes" = gene_count,
      "Number of genes passing filter" = genes_passing,
      "Percentage of genes passing filter" = genes_passing/gene_count,
      "Number of genes not passing filter" = gene_count - genes_passing,
      "Percentage of genes not passing filter" = (gene_count - genes_passing) / gene_count
    )
    
    return(table)
  }
  
  make_heatmap <- function(filtered_counts) {
    
    counts_heatmap <- heatmap(as.matrix(filtered_counts[,-1]))
    return(counts_heatmap)
  }
  
  pca_plot <- function(filtered_counts) {
    # Perform PCA on the counts data
    pca_result <- prcomp(filtered_counts[,-1], scale. = TRUE)  # Use scale. = TRUE to scale the data
    
    # Extract the top two principal components
    pca_components <- pca_result$x[, 1:2]
    
    # Convert the PCA components to a data frame for plotting
    pca_df <- as.data.frame(pca_components)
    colnames(pca_df) <- c("PC1", "PC2")
    
    # Scatter plot of the top two principal components
    pca <- ggplot(pca_df, aes(x = PC1, y = PC2)) + 
        geom_point() +
      
         labs(xlab = "Principal Component 1", ylab = "Principal Component 2",
         title = "Scatter plot of Top 2 Principal Components")
    
    return(pca)
  }
  output$table <- renderTable({draw_table(load_data()$val1,load_data()$val2)})
  
  output$diag_plot <- renderPlot({diagnostic_plots(load_data()$val1,load_data()$val2)})
  
  output$heatmap <- renderPlot({make_heatmap(load_data()$val2)})
  
  output$pca <- renderPlot({pca_plot(load_data()$val2)})
  
  
}

# Run the application
shinyApp(ui = ui, server = server)