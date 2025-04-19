# title: sla - early analysis
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


# SLA---------------------------------------------------------------------------

sla_untidy <- read_excel("code-and-data/data/2023-greenhouse-data_sla-ldmc.xlsx")
head(sla_untidy)

sla_tidy <- sla_untidy |> 
  mutate_all(~ifelse(. %in% c("", "NA"), NA, .)) |> 
  filter(!is.na(sla)) |> 
  mutate(tray = as.factor(tray)) |>
  mutate(id = as.factor(id)) |> 
  mutate(pop = as.factor(pop)) |> 
  mutate(TGP = as.factor(TGP)) |>
  mutate(sla = as.numeric(sla))
head(sla_tidy)

hist(sla_tidy$sla)

sla_tidy_1 <- sla_tidy |> 
  filter(sla < 250) |> 
  filter(sla > 45) |> 
  mutate(log_sla = log(sla))

hist(sla_tidy_1$log_sla)

# models
sla_mod_1 <- lmer(log_sla ~ TGP + (1|tray), data = sla_tidy_1)
summary(sla_mod_1)

sla_mod_2 <- lm(log_sla ~ TGP , data = sla_tidy_1)
summary(sla_mod_2)

AIC(sla_mod_1)
AIC(sla_mod_2)
#sla mod 2 is best? interesting

check_model(sla_mod_2)
anova(sla_mod_1, sla_mod_2, test="Chisq")
#model 1 is the best according to chisq. ask rachel about this

sla_mod_2 |> emmeans(pairwise ~ TGP) 

# visualizing
sla_tidy_1 |> 
  #group_by(TGP) |> 
  ggplot(aes(x = TGP, y = sla)) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.2), alpha = 0.5) +
  #geom_signif(comparisons = list(c("CC", "DC")), 
  #map_signif_level=TRUE) +
  xlab(" ") +
  ylab("Height (cm)") +
  ggtitle(" ")