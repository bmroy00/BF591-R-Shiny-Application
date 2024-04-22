library(shiny)
library(bslib)
library(ggplot2)
library(tidyverse)
library(colourpicker)
library(beeswarm)
library(gplots)

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
      sliderInput(inputId = 'var_slider', min = 0, max = 40,
                  label = 'Include Genes with X Percentile of Variance:',
                  value = 1,
                  step = .5),
      
      sliderInput(inputId = 'non_zero_slider', min = 0, max = 69,
                  label = 'Include Genes with at least X non-zero samples:',
                  value = 50,
                  step = 1),
      
      # Button to generate plots and table
      submitButton(text = 'Submit', icon = NULL)
    ),
    
    # Organize where our plots and table exist
    mainPanel(
      tabsetPanel(
        tabPanel("Table", tableOutput('table')),
        tabPanel("Plots", plotOutput('var_plot'),
                 plotOutput("zeros_plot")),
        tabPanel("Heatmap", plotOutput('heatmap')),
        tabPanel("PCA",
                 sliderInput(inputId = 'pca_slider1', min = 1, max = 10,
                             label = 'Choose N Principal Components to Plot',
                             value = 1,
                             step = 1),
                 sliderInput(inputId = 'pca_slider2', min = 1, max = 10,
                             label = 'Choose N Principal Components to Plot',
                             value = 2,
                             step = 1),
                 plotOutput('pca'))
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
    return(df)
  })
  
  #' Diagnostic plots
  #' Zeros Plot
  zeros_plot <-
    function(dataf, non_zero_slider) {
      medians <- apply(dataf[,-1], 1, median)
      variances <- apply(dataf[,-1], 1, var)
      zeros_count <- apply(dataf[,-1] != 0, 1, sum)
      
      plot <- ggplot(dataf, aes(x = log2(medians), y = zeros_count, 
                        col = zeros_count >= non_zero_slider)) +
        geom_point() +
        labs(x = "log2(Median Counts)", y = "Number of Non-Zeros") +
        title("Median Counts vs. Number of Non-Zeros") +
        scale_color_manual(values=c('white','violet')) +
        theme_dark()
      
      return(plot)
    }
  
  # Variance Plot
  var_plot <-
    function(dataf, var_slider) {
      medians <- apply(dataf[,-1], 1, median)
      variances <- apply(dataf[,-1], 1, var)
      zeros_count <- apply(dataf[,-1] != 0, 1, sum)
    
      plot <- ggplot(dataf, aes(x = log10(medians), y = log10(variances), 
                                col = (variances/sum(variances)) >= (var_slider/100))) +
        geom_point() +
        labs(x = "log10(Median Counts)", y = "log10(Variances)") +
        title("Median Counts vs. Variances") +
        scale_color_manual(values=c('white','limegreen')) +
        theme_dark()
      
      return(plot)
    }
  
  #' Table
  draw_table <- function(dataf, var_slider, non_zero_slider) {
    # Total number of samples
    samp_count <- ncol(dataf[,-1]) 
    
    # Total number of genes
    gene_count <- nrow(dataf)
    
    medians <- apply(dataf[,-1], 1, median)
    variances <- apply(dataf[,-1], 1, var)
    zeros_count <- apply(dataf[,-1] != 0, 1, sum)
    
    proportion_variance <- variances/sum(variances)
    
    var_filter <- which(proportion_variance >= var_slider/100)
    var_gene_names <- rownames(dataf)[var_filter]
    var_filtered_df <- dataf[var_gene_names, ]
    
    # Include genes with at least X non-zero samples
    non_zero_samples <- apply(var_filtered_df[,-1] != 0, 1, sum)
    genes_nonzero <- which(non_zero_samples >= non_zero_slider)
    genes_names_nonzero <- rownames(dataf)[genes_nonzero]
    filtered_counts <- dataf[genes_names_nonzero, ]
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
  
  make_heatmap <- function(dataf, var_slider, non_zero_slider) {
    variances <- apply(dataf[,-1], 1, var)
    zeros_count <- apply(dataf[,-1] != 0, 1, sum)
    
    proportion_variance <- variances/sum(variances)
    
    var_filter <- which(proportion_variance >= var_slider/100)
    var_gene_names <- rownames(dataf)[var_filter]
    var_filtered_df <- dataf[var_gene_names, ]
    
    # Include genes with at least X non-zero samples
    non_zero_samples <- apply(var_filtered_df[,-1] != 0, 1, sum)
    genes_nonzero <- which(non_zero_samples >= non_zero_slider)
    genes_names_nonzero <- rownames(dataf)[genes_nonzero]
    filtered_counts <- dataf[genes_names_nonzero, ]
    
    rownames(filtered_counts) = filtered_counts$Column1
    
    counts_heatmap <- heatmap.2(as.matrix(log10(filtered_counts[,-1] + 1),
                                        Colv = NA, Rowv = NA,
                                        scale = "none",
                                        dendrogram = "none",
                                        trace = "none",
                                        key = TRUE, key.title = "Legend Title",
                                        key.xlab = "Value",
                                        keysize = 1.0,
                                        key.legend = TRUE,
                                        col = colorRampPalette(brewer.pal(9, "Blues"))(100)))  # Adjust color palette as needed))
  
    

    return(counts_heatmap)
  }
  
  pca_plot <- function(dataf, var_slider, non_zero_slider, pca_slider1, pca_slider2) {

    pca_result <- prcomp(dataf[,-1])  # Transpose the count data if necessary
    
    # Number of top principal components to plot
    pc1 <- pca_slider1 
    pc2 <- pca_slider2
    
    # Extract top n principal components
    pca_df <- as.data.frame(pca_result$x[, pc1:pc2])
    
    variance_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
    
    PC1 <- paste0("PC", pc1)
    PC2 <- paste0("PC", pc2)
    
    pca <- ggplot(data = pca_df, aes(x = PC1, y = PC2)) +
      geom_point() +
      labs(x = paste0(PC1, "(", round(variance_explained[pc1] * 100, 2), "%)"), 
           y = paste0(PC2, "(", round(variance_explained[pc2] * 100, 2), "%)")) +
      ggtitle(paste0("PCA: Principal Components ", pc1, " vs. ", pc2))
    
    return(pca)
  }
  output$table <- renderTable({draw_table(load_data(), input$var_slider, input$non_zero_slider)})
  
  output$zeros_plot <- renderPlot({zeros_plot(load_data(), input$non_zero_slider)})
  output$var_plot <- renderPlot({var_plot(load_data(), input$var_slider)})
  
  output$heatmap <- renderPlot({make_heatmap(load_data(),input$var_slider,input$non_zero_slider)})
  
  output$pca <- renderPlot({pca_plot(load_data(),input$var_slider,input$non_zero_slider,input$pca_slider1, input$pca_slider2)})
  
  
}

# Run the application
shinyApp(ui = ui, server = server)