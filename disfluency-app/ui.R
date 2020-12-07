source("disfluency-app/environment.R")

ui <- fluidPage(
  # Setup title part ---------------------------
  titlePanel(
    windowTitle = "Disfluency app", #this appears in the browser tab 
    title = h1("Typing-disfluency app",
               h4("by Jens Roeser (Nottingham Trent University)"),
               h4("email: jens.roeser@ntu.ac.uk"))),
    # Themes changes app layout
    theme = shinytheme("yeti"),
    navbarPage("Disfluency app",
               # Input for item selector ----------------------
               tabPanel("Example Data", 
                  sidebarPanel(selectInput('N_example',
                                           'Select participants',
                                           choices = unique(d$subj),
                                           selected = unique(d$subj)[sample(seq(5))],
                                           multiple = T),
                               selectInput('scaled',
                                           'IKI scale',
                                           choices = c("Normal", "Log-scaled"),
                                           selected = "Normal"),
                               switchInput(inputId = "group", label = "Grouping",value = TRUE)),
                  # Select number of ppts (later select ppt column)
                  # Select grouping variable
                  mainPanel(
                    tabsetPanel(id = "example",
                    tabPanel("Densities", plotOutput('density', width = "100%")),
                    tabPanel("Time course"),
                    tabPanel("Data preview", dataTableOutput("example")),
                    tabPanel("Mixture modelling", h4("Summary"), dataTableOutput("mog")) # add a by group option
                    # allow download of by ppt results of pause and writing
                  )), 
               )
    )
  )