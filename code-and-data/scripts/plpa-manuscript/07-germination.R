# madeleine wallace
# 11-13-2025

# this script is to test if germination rates were different between populations
# as per review 1 from new phytologist!

################################################################################

# packages
library(tidyverse)
library(lmerTest)
library(lme4)
library(emmeans)

################################################################################

# load data and remove trial 1
germ <- read_excel("code-and-data/data/2023-greenhouse-data_germination.xlsx")
germ <- germ |> 
  filter(trial == 2)

# format data
germ$date_germ <- as.Date(germ$date_germ)

# trial 2 started on 5/31/23
sow_date <- as.Date("2023-05-31")

# calculate days to germ
germ$days_to_germ <- as.numeric(germ$date_germ - sow_date)

#checking variables
germ$pop <- as.factor(germ$pop)
germ$`2021_treat` <-as.factor(germ$`2021_treat`)
germ$`2023_treat` <-as.factor(germ$`2023_treat`)

################################################################################

# explore/analyze data
boxplot(days_to_germ ~ pop, data = germ)
boxplot(days_to_germ ~ `2021_treat`, data = germ)
boxplot(days_to_germ ~ `2023_treat`, data = germ)

# fit model
germ_model <- glm(days_to_germ ~ pop*`2021_treat`, family = poisson, data = germ)
summary(germ_model)
anova(germ_model)

# visualize
library(ggplot2)
ggplot(germ, aes(x = pop, y = days_to_germ, color = `2021_treat`)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# emmeans
emm <- emmeans(germ_model, ~ pop * `2021_treat`)
emm
emm_resp <- summary(emm, type = "response")
emm_resp

# plot emmeans output
plot(emm, comparisons = TRUE, by = "pop", type = "response")

emm_df <- as.data.frame(emm_resp)

ggplot(emm_df, aes(x = pop, y = rate, color = `2021_treat`, group = `2021_treat`)) +
  geom_point(position = position_dodge(width = 0.3), size = 2) +
  geom_line(position = position_dodge(width = 0.3)) +
  labs(y = "Estimated days to germination", x = "Population",
       color = "2021 Treatment") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# decided not to add this analysis
# I can’t test for effects of TGP because I had an establishment period of watering
# so offspring treatment doesn’t matter… 
# while this result is interesting, I don’t think it adds to my story?

