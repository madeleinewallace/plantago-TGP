# title: heights - early analysis
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
library(ggpubr)
library(cowplot)


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


# HEIGHTS ----------------------------------------------------------------------

height_untidy <- read_excel("code-and-data/data/2023-greenhouse-data_heights.xlsx")
head(height_untidy)

height_tidy <- height_untidy |> 
  mutate(tray = as.factor(tray)) |>
  mutate(pop = as.factor(pop)) |> 
  mutate(id = as.factor(id)) |> 
  mutate(TGP = as.factor(TGP)) |>
  mutate(`2021_treat` = as.factor(`2021_treat`)) |> 
  mutate(`2023_treat` = as.factor(`2023_treat`)) |> 
  mutate(height_cm_50 = as.numeric(height_cm_50)) |> 
  mutate(height_cm_14 = as.numeric(height_cm_14)) |> 
  filter(!is.na(height_cm_7))

# avg height on measurement day 50 by treatment class - model
height_tidy_50 <- height_tidy |> 
  filter(height_cm_50!= 0 )

hist(height_tidy_50$height_cm_50)
mean(height_tidy_50$height_cm_50)
var(height_tidy_50$height_cm_50)
#poisson distribution with continous data :/
#log did not work :|
#google says forge ahead. okkkkk

height_mod_1 <- lmer(height_cm_50 ~ TGP + (1|tray), data = height_tidy_50)
summary(height_mod_1)

height_mod_2 <- lm(height_cm_50 ~ TGP , data = height_tidy_50)
summary(height_mod_2)

AIC(height_mod_1)
AIC(height_mod_2)
#height mod 1 is best

check_model(height_mod_1)
anova(height_mod_1, height_mod_2)
#model 1 is the best!

height_mod_1 |> emmeans(pairwise ~ TGP) 

# avg height on measurement day 50 by treatment class - visual


# effects plot
plot_model(height_mod_1, 
           vline.color = "black", 
           show.values = TRUE, 
           value.offset = .2, 
           colors = "bw",
           title = "Flowering")

#boxplot
height_tidy_50 |> 
  #group_by(TGP) |> 
  #summarise(mean_height_32 = mean(height_cm_32)) |> 
  ggplot(aes(x = TGP, y = height_cm_50)) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.2), alpha = 0.5) +
  xlab(" ") +
  ylab("Height (cm)") +
  ggtitle(" ")

#barplot wrangle (so messy!)
df0 <- height_tidy_50 |> 
  group_by(TGP) |> 
  summarise(count = n())

df1 <- height_tidy_50 |>  
  group_by(TGP) |> 
  summarise(mean_height_50 = mean(height_cm_50))

df2 <- tibble(TGP = c("CC", "CD", "DC", "DD"),
  se = c(0.25247, 0.36979, 0.16394, 0.37163))

combined_df <- bind_cols(df1, df2 %>% select(-TGP))

DD_CD_heights_rownames <- c("DD", "CD")

DD_CD_heights <- combined_df |> 
  filter(TGP %in% DD_CD_heights_rownames)
DD_CD_heights$TGP <- factor(DD_CD_heights$TGP, levels = c("CD", "DD"))

CC_DC_heights <- combined_df %>%
  filter(!TGP %in% DD_CD_heights_rownames)
CC_DC_heights$TGP <- factor(CC_DC_heights$TGP, levels = c("DC", "CC"))

#theme
source("code-and-data/scripts/alvsce-poster-ggplot-theme.R")

heights_x_labels <- setNames(c("control 2021 * control 2023", "control 2021 * drought 2023", "drought 2021 * control 2023", "drought 2021 * drought 2023"),
                        c("CC", "CD", "DC", "DD"))

#barplots
DD_CD_heights_plot <- DD_CD_heights |>
  ggplot(aes(x = TGP, y = mean_height_50)) +
  geom_bar(stat = "identity", fill="#01A08A", width = 0.5) +
  geom_errorbar(aes(ymin = mean_height_50-se, ymax = mean_height_50+se), 
                width=.1) +
  geom_text(aes(label = " "), position = position_dodge(width = .9), vjust = -3, size = 10) +
  ylim(0, 6) +
  alvsce_poster_theme() +
  xlab(" ") +
  ylab("height at day 50 (cm)")
print(DD_CD_heights_plot)


CC_DC_heights_plot <- CC_DC_heights |>
  ggplot(aes(x = TGP, y = mean_height_50)) +
  geom_bar(stat = "identity", fill="#01A08A", width = 0.5) +
  geom_errorbar(aes(ymin = mean_height_50-se, ymax = mean_height_50+se), 
                width=.1) +
  geom_text(aes(label = "*"), position = position_dodge(width = .9), vjust = -0.5, size = 12) +
  geom_text(aes(label = "p = 0.0265"), position = position_dodge(width = .9), vjust = -4, size = 24) +
  ylim(0, 6) +
  alvsce_poster_theme() +
  xlab(" ") +
  ylab("height at day 50 (cm)")
print(CC_DC_heights_plot)

ggsave("DD_CD_heights_plot.png", DD_CD_heights_plot, dpi = 600, path = "code-and-data/output")
ggsave("CC_DC_heights_plot.png", CC_DC_heights_plot, dpi = 600, path = "code-and-data/output")
