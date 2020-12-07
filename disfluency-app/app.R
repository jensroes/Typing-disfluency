# Packages
rm(list = ls())
library(shiny)

# get user interface (ui) and server (defines what is happening)
source("disfluency-app/ui.R")
source("disfluency-app/server.R")

# Create app
shinyApp(ui = ui, server = server)
