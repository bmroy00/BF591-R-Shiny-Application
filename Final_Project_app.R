library(shiny)
library(bslib)
library(ggplot2)
library(tidyverse)
library(colourpicker)
library(DT)
library(shinythemes)
library(gplots)


ui <- fluidPage(
  theme = shinytheme("simplex"),
  titlePanel("BF 591 Final Project"),
  tabsetPanel(
    
    # Samples UI
    tabPanel("Samples",
             sidebarLayout(
               sidebarPanel(
                 fileInput(inputId = "samples_file", "Input Samples File:",
                           accept=c(
                             "text/csv",
                             "text/comma-separated-values,text/plain",
                             ".csv")
                 ),
                 
                 # Button to generate plots and table
                 submitButton(text = 'Submit', icon = NULL)
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Summary", tableOutput('summary')),
                   tabPanel("Table",DTOutput('samples_table')),
                   tabPanel("Plot", 
                            # Button for category
                            radioButtons(inputId = 'categorical', 'Choose the category to group by',
                                         choices = c("Diagnosis",
                                                     "Vonsattel.Grade"),
                                         selected = 'Diagnosis'
                            ),
                            # Button for continuous variable 
                            radioButtons(inputId = 'continuous', 'Choose the column for the y-axis',
                                         choices = c("Post.Mortem.Interval",
                                                     "Age.of.Death",
                                                     "RNA.Integrity.Number",
                                                     "mRNA.seq.Reads",
                                                     "Age.of.Onset",
                                                     "Duration",
                                                     "CAG",
                                                     "h-v.striatal.score",
                                                     "h-v.cortical.score"),
                                         selected = "Age.of.Death"
                            ),
                            plotOutput('samples_plot'))
                 )
               )
             )
    ),
    
    # Counts UI
    tabPanel("Counts",
             sidebarLayout(
               sidebarPanel(
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
               mainPanel(
                 tabsetPanel(
                   tabPanel("Table", tableOutput('counts_table')),
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
                                        value = 3,
                                        step = 1),
                            plotOutput('pca'))
                 )
               )
             )
    ),
    
    
    # DE UI
    tabPanel("DE",
             sidebarLayout(
               sidebarPanel(
                 fileInput(inputId = "de_file", "Input Differential Expression File:",
                           accept=c(
                             "text/csv",
                             "text/comma-separated-values,text/plain",
                             ".csv")
                 ),
                 submitButton(text = 'Submit', icon = NULL)
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Table", DTOutput('de_table')),
                   tabPanel("Plot", plotOutput('volcano_plot'))
                 )
               )
             )
    ),
    
    # FGSEA UI
    tabPanel("FGSEA",
             sidebarLayout(
               sidebarPanel(
                 fileInput(inputId = "fgsea_file", "Input fgsea File:",
                           accept=c(
                             "text/csv",
                             "text/comma-separated-values,text/plain",
                             ".csv")
                 ),
                 
                 # Button to generate plots and table
                 submitButton(text = 'Submit', icon = NULL)
               ),
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
  )
)

server <- function(input, output) {
  options(shiny.maxRequestSize = 100*1024^2)
  
  
  # Sample Server Portion
  load_sample_data <- reactive({
    df <- read.csv(input$samples_file$datapath)
    
    return(df)
  })
  
  #' 
  
  create_summary <-
    function(dataf) {
      dataf <- data_frame(dataf)
      column_types <- sapply(dataf, typeof)
      third_col <- c()
      
      get_column_summary <- function(column) {
        if (is.numeric(column)) {
          mean_value <- mean(column, na.rm = TRUE)
          sd_value <- sd(column, na.rm = TRUE)
          return(paste(mean_value, " (+/- ", sd_value, ")"))
        } else {
          unique_values <- unique(column)
          return(toString(unique_values))
        }
      }
      
      third_col <- sapply(dataf, get_column_summary)
      
      summary <- data.frame(
        'Column Name' = colnames(dataf),
        'Type' = column_types,
        'Mean (sd) or Distinct Values' = third_col
      )
      
      return(summary)
    }
  
  #' Summary Table
  
  draw_sample_table <- function(dataf) {
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
  output$summary <- renderTable({create_summary(load_sample_data())})
  
  output$samples_table <- renderDT({draw_sample_table(load_sample_data())})
  
  output$samples_plot <- renderPlot({create_plot(load_sample_data(),input$continuous,input$categorical)})
  
  
  # Counts Server Portion
  
  #' load_Data
  load_counts_data <- reactive({
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
  draw_counts_table <- function(dataf, var_slider, non_zero_slider) {
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
  output$counts_table <- renderTable({draw_counts_table(load_counts_data(), input$var_slider, input$non_zero_slider)})
  
  output$zeros_plot <- renderPlot({zeros_plot(load_counts_data(), input$non_zero_slider)})
  output$var_plot <- renderPlot({var_plot(load_counts_data(), input$var_slider)})
  
  output$heatmap <- renderPlot({make_heatmap(load_counts_data(),input$var_slider,input$non_zero_slider)})
  
  output$pca <- renderPlot({pca_plot(load_counts_data(),input$var_slider,input$non_zero_slider,input$pca_slider1, input$pca_slider2)})
  
  # DE Server Portion
  
  load_de_data <- reactive({
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
  
  draw_de_table <- function(dataf) {
    table <- as.data.frame(dataf)
    
    return(table)
  }
  #' These outputs aren't really functions, so they don't get a full skeleton, 
  #' but use the renderPlot() and renderTabel() functions to return() a plot 
  #' or table object, and those will be displayed in your application.
  # Same here, just return the table as you want to see it in the web page
  output$de_table <- renderDT({draw_de_table(load_de_data())})
  output$volcano_plot <- renderPlot({volcano_plot(load_de_data())})
  
  
  # FGSEA Server Portion
  
  load_fgsea_data <- reactive({
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
      labs(x = "NES", y = "-log10(padj)", title = paste0("Scatter Plot of NES by -log10(padj)"),
                                                         color = "Passes P-value Threshold") +
      scale_color_manual(values = c("grey", "purple")) + 
      theme_minimal()
    
    return(plot)
  }
  
  generate_table <- reactive ({
    table <- draw_fgsea_table(load_fgsea_data())
    return(table)
  })
  
  output$nes_plot <- renderPlot({plot_top_pathways(load_fgsea_data(), input$num_paths)})
  
  output$nes_table <- renderDT({draw_fgsea_table(load_fgsea_data())})
  
  output$nes_scatter_plot <- renderPlot({nes_scatter_plot(load_fgsea_data())})
  
  output$downloadTable <- downloadHandler(
    filename = function() {
      return("fgsea_table.csv")
    },
    content = function(file) {
      write.csv(generate_table(), file)
    }
  )
  
}

shinyApp(ui = ui, server = server)