library(shiny)

ui <- fluidPage(
  titlePanel("Tabs with Separate Sidebars"),
  tabsetPanel(
    tabPanel("Tab 1",
             sidebarLayout(
               sidebarPanel(
                 # Input elements for Tab 1 sidebar
                 fileInput("file1", "Choose file for Tab 1")
               ),
               mainPanel(
                 # Main panel content for Tab 1
                 textOutput("text_output1")
               )
             )
    ),
    tabPanel("Tab 2",
             sidebarLayout(
               sidebarPanel(
                 # Input elements for Tab 2 sidebar
                 fileInput("file2", "Choose file for Tab 2")
               ),
               mainPanel(
                 # Main panel content for Tab 2
                 textOutput("text_output2")
               )
             )
    ),
    tabPanel("Tab 3",
             sidebarLayout(
               sidebarPanel(
                 # Input elements for Tab 3 sidebar
                 fileInput("file3", "Choose file for Tab 3")
               ),
               mainPanel(
                 # Main panel content for Tab 3
                 textOutput("text_output3")
               )
             )
    ),
    tabPanel("Tab 4",
             sidebarLayout(
               sidebarPanel(
                 # Input elements for Tab 4 sidebar
                 fileInput("file4", "Choose file for Tab 4")
               ),
               mainPanel(
                 # Main panel content for Tab 4
                 textOutput("text_output4")
               )
             )
    )
  )
)

server <- function(input, output) {
  # Output text based on selected files for each tab
  output$text_output1 <- renderText({
    if (!is.null(input$file1)) {
      paste("File selected for Tab 1:", input$file1$name)
    }
  })
  output$text_output2 <- renderText({
    if (!is.null(input$file2)) {
      paste("File selected for Tab 2:", input$file2$name)
    }
  })
  output$text_output3 <- renderText({
    if (!is.null(input$file3)) {
      paste("File selected for Tab 3:", input$file3$name)
    }
  })
  output$text_output4 <- renderText({
    if (!is.null(input$file4)) {
      paste("File selected for Tab 4:", input$file4$name)
    }
  })
}

shinyApp(ui = ui, server = server)
