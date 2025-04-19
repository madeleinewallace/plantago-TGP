# title: biomass - early analysis
# author: madeleine wallace
# date: 2024-03-28



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
 



# BIOMASS-----------------------------------------------------------------------

# wrangle
biomass_untidy <- read_excel("code-and-data/data/2023-greenhouse-data_biomass.xlsx")
biomass_untidy <- biomass_untidy |> 
  mutate(total_mass_g = ifelse(total_mass_g == 0, NA, total_mass_g))
head(biomass_untidy)

biomass_tidy <- biomass_untidy |> 
  mutate(tray = as.factor(tray)) |>
  mutate(population = as.factor(population)) |> 
  mutate(id = as.factor(id)) |>
  mutate(TGP = as.factor(TGP)) |> 
  mutate(total_mass_g = as.numeric(total_mass_g)) |> 
  select(-`mass_g_minus_SLA`, -`mass_g_SLA`)

# add R:S ratio
biomass_tidy_1 <- biomass_tidy |> 
  pivot_wider(names_from = plant_part, values_from = total_mass_g) |> 
  mutate(total_biomass = root + shoot) |> 
  mutate(ratio_RS = root/shoot)

# i think i should log the RS ratio
biomass_tidy_1 |>  
  ggplot(aes(x = log(ratio_RS))) + 
  geom_density()

biomass_tidy_2 <- biomass_tidy_1 |> mutate(log_ratio_RS = log(ratio_RS))

#i think i should get rid of the outliers
#outliers <- boxplot(biomass_tidy_2$log_ratio_RS, plot=FALSE)$out
#x <- biomass_tidy_2
#x <- x[-which(x$log_ratio_RS %in% outliers),]



# MODELS------------------------------------------------------------------


#TOTAL BIOMASS------------------------------------------------------------------
total_biomass_mod_1 <- lmer(total_biomass ~ TGP + (1|tray), data = biomass_tidy_2)
summary(total_biomass_mod_1)

total_biomass_mod_2 <- lm(total_biomass ~ TGP , data = biomass_tidy_2)
summary(total_biomass_mod_2)

AIC(total_biomass_mod_1)
AIC(total_biomass_mod_2)
#mod 1 it is

check_model(total_biomass_mod_1)
anova(total_biomass_mod_1, total_biomass_mod_2, test="Chisq")
#model 1 is the best according to chisq also!!

total_biomass_mod_1 |> emmeans(pairwise ~ TGP) 


# ROOT:SHOOT--------------------------------------------------------------------
ratio_biomass_mod_1 <- lmer(log_ratio_RS ~ TGP + (1|tray), data = biomass_tidy_2)
summary(ratio_biomass_mod_1)

ratio_biomass_mod_2 <- lm(log_ratio_RS ~ TGP , data = biomass_tidy_2)
summary(ratio_biomass_mod_2)

AIC(ratio_biomass_mod_1)
AIC(ratio_biomass_mod_2)
#mod 2 it is

check_model(ratio_biomass_mod_2)
anova(ratio_biomass_mod_1, ratio_biomass_mod_2, test="Chisq")
#model 1 is a little better but im sticking with model 2. ask rachel why

ratio_biomass_mod_2 |> emmeans(pairwise ~ TGP) 



# VISUALIZATION------------------------------------------------------------------

# total biomass
biomass_tidy_2  |> 
  ggplot(aes(x = TGP, y = log_ratio_RS)) + geom_violin(aes(color = TGP))

biomass_tidy_2 |> 
  #group_by(TGP) |> 
  ggplot(aes(x = TGP, y = total_biomass)) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.2), alpha = 0.5)

ratio_biomass_mod_2 |> 
  #group_by(TGP) |> 
  ggplot(aes(x = TGP, y = log_ratio_RS)) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.2), alpha = 0.5)

# alvsce poster session graph - R:S ratio

#theme
source("code-and-data/scripts/alvsce-poster-ggplot-theme.R")

heights_x_labels <- setNames(c("control 2021 * control 2023", "control 2021 * drought 2023", "drought 2021 * control 2023", "drought 2021 * drought 2023"),
                             c("CC", "CD", "DC", "DD"))

#graph - violin
biomass_tidy_2$TGP <- factor(biomass_tidy_2$TGP, levels = c("CD", "DD", "DC", "CC"))

biomass_poster_plot <- biomass_tidy_2  |> 
  ggplot(aes(x = TGP, y = log_ratio_RS, fill=TGP)) + 
  geom_violin() +
  alvsce_poster_theme() +
  theme(legend.position="none") +
  scale_color_manual(values=c("CC"="#5CBCD5", "CD"="#FE0000", "DC"="#F98401", "DD"="#01A08A")) +
  scale_fill_manual(values=c("CC"="#5CBCD5", "CD"="#FE0000", "DC"="#F98401", "DD"="#01A08A")) +
  xlab(" ") +
  ylab("root:shoot ratio")

print(biomass_poster_plot)
ggsave("biomass_poster_plot.png", biomass_poster_plot, dpi = 600, path = "code-and-data/output")

# graph - boxplot
biomass_poster_plot_2 <- biomass_tidy_2  |> 
  ggplot(aes(x = TGP, y = log_ratio_RS, fill=TGP)) + 
  geom_boxplot() +
  alvsce_poster_theme() +
  theme(legend.position="none") +
  scale_color_manual(values=c("CC"="#5CBCD5", "CD"="#FE0000", "DC"="#F98401", "DD"="#01A08A")) +
  scale_fill_manual(values=c("CC"="#5CBCD5", "CD"="#FE0000", "DC"="#F98401", "DD"="#01A08A")) +
  xlab(" ") +
  ylab("root:shoot ratio")

print(biomass_poster_plot_2)
ggsave("biomass_poster_plot_2.png", biomass_poster_plot_2, dpi = 600, path = "code-and-data/output")
