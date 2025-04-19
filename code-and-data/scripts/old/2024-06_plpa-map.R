# title: PLPA map
# author: Maddie Wallace
# date: 2024-06-12

# figure 1 for PLPA greenhouse manuscript - a map of PLPA points, table of precip/temp variables, experimental design, pic of PLPA

# packages ----
library(raster)
library(tidyverse)
library(maps)
library(gtsummary)
library(readxl)
library(kableExtra)

# load data
climate_data <- read_excel("code-and-data/data/2021-climate-info.xlsx")

table <- read_excel("writing/table-1.xlsx", sheet = "Sheet2")
