library(shiny)

# Define UI
ui <- fluidPage(
  navbarPage("Combined App",
             tabPanel("Samples", uiOutput("samples_app.R")),
             tabPanel("Counts", uiOutput("counts_app.R")),
             tabPanel("DE", uiOutput("de_app.R")),
             tabPanel("FGSEA", uiOutput("fgsea_app.R"))
  )
)

# Define Server
server <- function(input, output) {
  # Source code from separate app files
  source("samples_app.R", local = TRUE)
  source("counts_app.R", local = TRUE)
  source("de_app.R", local = TRUE)
  source("fgsea_app.R", local = TRUE)
}

# Run the app
shinyApp(ui = ui, server = server)
