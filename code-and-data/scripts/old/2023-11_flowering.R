# title: flowering - early analysis
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
  filter(!is.na(population))

# FLOWERING --------------------------------------------------------------------

# wrangle
flower_untidy <- read_excel("code-and-data/data/2023-greenhouse-data_flowering.xlsx")
head(flower_untidy)

flower_tidy <- flower_untidy |> 
  mutate(tray = as.factor(tray)) |>
  mutate(id = as.factor(id)) |> 
  mutate(pop = as.factor(pop)) |> 
  mutate(TGP = as.factor(TGP)) |> 
  mutate(flowering = as.factor(flowering)) |> 
  filter(!is.na(date_germ)) |> 
  mutate(flowering = if_else(flowering == 'Y', 1, 0)) |> # 1 means yes flowering
  mutate(flowering = ifelse(is.na(flowering), 0, flowering))
head(flower_tidy)

# join with climate data (if i want it)
flower_climate <- inner_join(flower_tidy, climate_ID_tidy, by = c("pop" = "population")) |> 
  mutate(temp_cv = as.factor(temp_cv), 
         tray = as.factor(tray),
         id = as.factor(id), 
         cv_temp_mon = as.factor(cv_temp_mon), 
         cvlevel = as.factor(cvlevel), 
         temp = as.factor(temp),
         monsoon = as.factor(monsoon))

# percent flowering
percent_flowering <- flower_tidy |>  
  group_by(TGP) |> 
  summarise(percent_flower = mean(flowering),
            sd_flower = sd(flowering), 
            se = std.error(flowering))  |> 
  mutate(percent_noflower = 1-percent_flower)

# number flowering
number_flowering <- flower_tidy |>  
  group_by(TGP) |> 
  summarise(total_flower_num = sum(flowering),
            sd_flower = sd(flowering), 
            se = std.error(flowering))
#view(number_flowering)

# MODELS - not caring about population....
flowering_mod_1 <- glm(flowering ~ TGP, family = binomial(link = "logit"), data = flower_climate)
summary(flowering_mod_1)
anova(flowering_mod_1)

flowering_mod_2 <- glmer(flowering ~ TGP + (1|tray), family = binomial(link = "logit"), data = flower_climate)
summary(flowering_mod_2)

AIC(flowering_mod_1)
AIC(flowering_mod_2)

check_model(flowering_mod_2)

anova(flowering_mod_1, flowering_mod_2, test="Chisq")
# flower_mod_2 it is!

# EMMEANS
emmeans(flowering_mod_2, pairwise ~ TGP, type = "response")
#flower_emmeans_result$contrasts |> 
#summary(infer = TRUE)
#flower_emmeans_result$emmeans |> 
#as.data.frame()
#head(flower_emmeans_result)

# VISUALIZE
# plot 1 - sum of flowered individuals per TGP class
flower_climate |> 
  group_by(TGP) |> 
  summarise(sum_flowered = sum(flowering)) |> 
  ggplot(aes(x = TGP, y = sum_flowered)) +
  geom_bar(stat = "identity", position = "dodge") +
  xlab(" ") +
  ylab("Total Flowered Individuals") +
  ggtitle("Number of Flowering Individuals per Treatment")

# plot 2 - predicted probabilities of flowering based on TGP class
predicted_data <- data.frame(TGP = flower_climate$TGP, 
                             tray = flower_climate$tray)

predicted_data <- data.frame(TGP = flower_climate$TGP, 
                             tray = flower_climate$tray, 
                             Predicted_Probability = predict(flowering_mod_2, newdata = predicted_data, type = "response"))
str(predicted_data)

ggplot(predicted_data, aes(x = TGP, y = Predicted_Probability)) +
  geom_bar(stat = "identity", position = "dodge") +
  xlab("TGP") +
  ylab("Probability of Flowering") +
  ggtitle("Bar Chart of Predicted Probabilities by TGP")

#plot 3 - coefficients plot
plot_model(flowering_mod_2, 
           vline.color = "black", 
           show.values = TRUE, 
           value.offset = .2, 
           colors = "bw",
           title = "Flowering")

#plot 4 - pairwise post hoc multiple comparisons
plot(flower_emmeans_result, comparisons = TRUE)

# flowering (day to)
# flowering (number of seed heads)


## FLOWERING & MORT  -----------------------------------------------------------
flomort <- full_join (percent_flowering, percent_dead)
flomort_2 <- full_join(flower_tidy, mort_tidy)

flomort %>% ggplot(aes(x = percent_flower, y = percent_dead, color = TGP)) + 
  geom_point() + 
  geom_smooth(method = 'lm')

flomort_2 %>% ggplot(aes(x = flowering, y = status, color = TGP)) + 
  geom_point() + 
  geom_smooth(method = 'glm')
# increased flowering increased death duh

# final plot
flower_climate |> 
  group_by(TGP) |> 
  summarise(sum_flowered = sum(flowering)) |> 
  ggplot(aes(x = TGP, y = sum_flowered)) +
  geom_bar(stat = "identity", position = "dodge") +
  xlab(" ") +
  ylab("Total Flowered Individuals") +
  ggtitle("Number of Flowering Individuals per Treatment")

#barplot wrangle (so messy)
df0 <- flower_climate |> 
  group_by(TGP) |> 
  summarise(count = n())

df1 <- flower_climate |>  
  group_by(TGP) |> 
  summarise(total_flowered = sum(flowering))

df2 <- flower_climate |>  
  group_by(TGP) |> 
  summarise(sd = sd(flowering)) |> 
  mutate(se = sd / (sqrt(n())))

combined_df <- bind_cols(df1, df2 %>% select(-TGP))

DD_CD_flowers_rownames <- c("DD", "CD")

DD_CD_flowers <- combined_df |> 
  filter(TGP %in% DD_CD_flowers_rownames)
DD_CD_flowers$TGP <- factor(DD_CD_flowers$TGP, levels = c("CD", "DD"))

CC_DC_flowers <- combined_df %>%
  filter(!TGP %in% DD_CD_flowers_rownames)
CC_DC_flowers$TGP <- factor(CC_DC_flowers$TGP, levels = c("DC", "CC"))

#theme

source("code-and-data/scripts/alvsce-poster-ggplot-theme.R")

heights_x_labels <- setNames(c("control 2021 * control 2023", "control 2021 * drought 2023", "drought 2021 * control 2023", "drought 2021 * drought 2023"),
                             c("CC", "CD", "DC", "DD"))

#barplots
DD_CD_flowers_plot <- DD_CD_flowers |>
  ggplot(aes(x = TGP, y = total_flowered)) +
  geom_bar(stat = "identity", fill="#F1AD00", width = 0.5) +
  #geom_errorbar(aes(ymin = total_flowered-se, ymax = total_flowered+se), 
                #width=.1) +
  geom_text(aes(label = "*"), position = position_dodge(width = .9), vjust = -0.5, size = 12) +
  #geom_text(aes(label = "p = 0.0237"), position = position_dodge(width = .9), vjust = -4, size = 24) +
  ylim(0, 100) +
  alvsce_poster_theme() +
  xlab(" ") +
  ylab("# of flowering plants")
print(DD_CD_flowers_plot)


CC_DC_flowers_plot <- CC_DC_flowers |>
  ggplot(aes(x = TGP, y = total_flowered)) +
  geom_bar(stat = "identity", fill="#F1AD00", width = 0.5) +
  #geom_errorbar(aes(ymin = total_flowered-se, ymax = total_flowered+se), 
                #width=.1) +
  #geom_text(aes(label = "p = 0.0665"), position = position_dodge(width = .9), vjust = -3, size = 24) +
  ylim(0, 100) +
  alvsce_poster_theme() +
  xlab(" ") +
  ylab("# of flowering plants")
print(CC_DC_flowers_plot)

ggsave("DD_CD_flowers_plot.png", DD_CD_flowers_plot, dpi = 600, path = "code-and-data/output")
ggsave("CC_DC_flowers_plot.png", CC_DC_flowers_plot, dpi = 600, path = "code-and-data/output")
