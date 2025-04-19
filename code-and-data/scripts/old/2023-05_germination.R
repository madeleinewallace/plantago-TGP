#load packages
library(tidyverse)
library(readxl)
library(writexl)
library(tidyr)
library(cowplot)
library(ggplot2)
library(janitor)
library(gridExtra)

#set WD
setwd("/Users/madeleinewallace/Library/CloudStorage/OneDrive-UniversityofArizona/1 MS Research/thesis")

#read in data
germination <- read_excel("MADDIE_GreenhouseData.xlsx", 
                          sheet="germination", 
                          col_types = c("text", "text", "text", "text", "text", "date"))
head(germination)

#calcuate total number of plants germinated so far
germination %>% count(`date_germ`)
248/2112
#11.7% germination rate KILL ME

#visualizing
germination_counts <- germination %>% count(pop_id, `date_germ`)
view(germination_counts)


#make df cus i'm dumb
germination_may29 <- data.frame(pop_id = c("1C", "1D", "3C", "3D", "4C", "4D", "5C", "5D", "6C", "6D", "7C", "7D", "8C", "8D", "9C", "9D", "10C", "10D", "11C", "11D", "12C", "12D"),
                                germinated = c(15, 13, 9, 5, 20, 12, 15, 18, 0, 0, 12, 9, 21, 11, 0, 0, 10, 18, 26, 16, 15, 3),
                                not_germinated = c(81, 82, 87, 91, 78, 83, 79, 79, 95, 97, 84, 87, 74, 86, 95, 96, 85, 79, 69, 81, 82, 94))

germination_may29$Pop_ID <- factor(germination_may23$pop_id, levels = c("1C", "1D", "3C", "3D", "4C", "4D", "5C", "5D", "6C", "6D", "7C", "7D", "8C", "8D", "9C", "9D", "10C", "10D", "11C", "11D", "12C", "12D"))

view(germination_may29)

grid.table(germination_may29)

ggplot(germination_may29, aes(x = pop_id, y = germinated)) +
  geom_bar(stat = "identity", width = 0.5) +
  geom_text(aes(label = germinated), vjust = -0.5) +
  scale_y_continuous(limit = c(0, 96)) +
  scale_x_discrete(limits = c("1C", "1D", "3C", "3D", "4C", "4D", "5C", "5D", "6C", "6D", "7C", "7D", "8C", "8D", "9C", "9D", "10C", "10D", "11C", "11D", "12C", "12D")) +
  theme_cowplot()
#each #_C should have 96 individuals

