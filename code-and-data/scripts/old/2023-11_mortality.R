# title: mortality - early analysis
# author: madeleine wallace
# date: 2023-11-1


# PACKAGES ---------------------------------------------------------------------
library(readxl)
library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggpubr)
library(performance)
library(ggplot2)
library(ggthemr)
library(sjPlot)
library(sjlabelled)
library(sjmisc)
library(DHARMa)
library(forcats)
library(plotrix)
library(patchwork)
library(MuMIn)


# CLIMATE DATA -----------------------------------------------------------------

#wrangle
climate_ID_untidy <- read_excel("code-and-data/data/2021-climate-info.xlsx", sheet = 2)

climate_ID_tidy <- climate_ID_untidy |> 
  mutate(population = as.factor(population)) |> 
  rename(temp = ...9) |> 
  dplyr::select(-c(notes)) |> 
  mutate(SAP_level = ifelse(SAP_mm >40, "wet","dry"),
         SAT_level = ifelse(SAT_C >16.6, "hot","cool")) |> 
  filter(!is.na(population))


# MORTALITY ---------------------------------------------------------------------

# wrangle
mort_untidy <- read_excel("code-and-data/data/2023-greenhouse-data_mortality.xlsx")
head(mort_untidy)

mort_tidy <- mort_untidy |> 
  mutate(tray = as.factor(tray)) |>
  mutate(id = as.factor(id)) |> 
  mutate(pop = as.factor(pop)) |> 
  mutate(TGP = as.factor(TGP)) |> 
  mutate(`2021_treat` = as.factor(`2021_treat`)) |> 
  mutate(`2023_treat` = as.factor(`2023_treat`)) |> 
  filter(!is.na(june_8)) |> 
  filter(trial == 2) |>
  pivot_longer(june_8:october_19, 
               names_to = 'date',
               values_to = 'status') |> 
  mutate(status = if_else(status == 'DEAD', 0, 1))  # 0 means dead
#view(mort_tidy)

# july 13th mortality, before most senescence, capturing treatment effects
july_13_mort <- mort_tidy |>  
  filter(date == 'july_13')

# percent dead
percent_dead <- mort_tidy |>  
  group_by(TGP, date) |> 
  summarise(percent_live = mean(status)) |> 
  filter(date == "july_13") |> 
  mutate(percent_dead = 1-percent_live)
#view(percent_dead)

# number dead
number_dead <- mort_tidy |>  
  group_by(TGP, date) |> 
  summarise(total_dead = sum(status)) |> 
  filter(date == "july_13")
#view(number_dead)

# MODELS - not caring about population....
mort_mod_1 <- glm(status ~ TGP, family = binomial(link = "logit"), data = mort_tidy)
summary(mort_mod_1)

mort_mod_2 <- glmer(status ~ TGP + (1|tray/id), family = binomial(link = "logit"), data = mort_tidy)
summary(mort_mod_2)

mort_mod_3 <- glmer(status ~ `2021_treat` + `2023_treat` + (1|tray/id), family = binomial(link = "logit"), data = mort_tidy)
summary(mort_mod_3)

AIC(mort_mod_1)
AIC(mort_mod_2)
AIC(mort_mod_3)
#mort mod 2 is best but no significance

simulationOutput1 <- simulateResiduals(fittedModel = mort_mod_2)
plot(simulationOutput1)
# OVERDISPERSED?

# visualizing

# average dead line graph
average_dead <- mort_tidy |> 
  group_by(TGP, date) |> 
  summarise(mean_status = mean(status))

average_dead$date <- as.Date(paste("2023", average_dead$date, sep = "_"), format = "%Y_%b_%d")
#view(average_dead)

average_dead_early <- average_dead |> 
  filter(date %in% c("2023-06-08", "2023-06-15", "2023-06-22", "2023-06-29", "2023-07-06",
                     "2023-07-13", "2023-07-20", "2023-07-27"))

ggplot(average_dead_early, aes(x = date, y = mean_status, group = TGP, color = TGP)) +
  geom_line() +
  labs(title = "Average Cumulative Dead by Date and TGP",
       x = "Date",
       y = "Average Cumulative Dead")

# mortality rate line graph
mort_rate <- mort_tidy |> 
  group_by(TGP, date) |> 
  summarise(
    total_dead = sum(status),  # Calculate the mortality rate
    n = sum(status == 0)) |>   # Calculate the sample size (n) 
  ungroup() |> 
  mutate(mort_rate = total_dead/n)

mort_rate$date <- as.Date(paste("2023", mort_rate$date, sep = "_"), format = "%Y_%b_%d")
#view(mort_rate)

mort_rate_early <- mort_rate |> 
  filter(date %in% c("2023-06-15", "2023-06-22", "2023-06-29", "2023-07-06",
                     "2023-07-13", "2023-07-20", "2023-07-27"))

ggplot(mort_rate_early, aes(x = date, y = mort_rate, group = TGP, fill = TGP)) +
  geom_bar(stat = "identity", position = position_dodge(width = 5), width = 4)
#scale_fill_manual(values = c("CC" = "black", "CD" = "darkblue", "DC" = "darkgreen", "DD" = "red"))