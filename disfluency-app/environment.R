# Packages
library(shiny)
library(shinythemes)
library(shinyWidgets)
library(tidyverse)
library(grid)
library(DT)
library(magrittr)
library(WVPlots)
source("functions/functions.R")

# Functions
#percentile <- function(x) return((x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm = T)))
percentile <- function(x) return( rank(x)/max(rank(x)) )  

# Plot layout
theme_set(theme_bw(base_size = 14) + 
            theme(strip.text = element_text(hjust = 0, colour = "black"),
                  strip.background = element_rect(fill = "grey97"),
                  panel.grid = element_line(colour = "grey98"),
                  legend.position = "top",
                  legend.justification = "right"))

# Read data
d <- read_csv("data/data.csv")


