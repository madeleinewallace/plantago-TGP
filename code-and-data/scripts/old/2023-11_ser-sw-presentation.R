# title: SER SW Presentation
# author: Maddie Wallace
# date: 2023-11-1

# packages ----
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
view(number_flowering)

# MODELS - not caring about population....
flowering_mod_1 <- glm(flowering ~ TGP, family = binomial(link = "logit"), data = flower_climate)
summary(flowering_mod_1)
anova(flowering_mod_1)

flowering_mod_2 <- glmer(flowering ~ TGP + (1|tray), family = binomial(link = "logit"), data = flower_climate)
summary(flowering_mod_2)

simulationOutput2 <- simulateResiduals(fittedModel = flowering_mod_2)
plot(simulationOutput2)

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


## flowering and mortality! ----------------------------------------------------
flomort <- full_join (percent_flowering, percent_dead)
flomort_2 <- full_join(flower_tidy, mort_tidy)

flomort %>% ggplot(aes(x = percent_flower, y = percent_dead, color = TGP)) + 
  geom_point() + 
  geom_smooth(method = 'lm')

flomort_2 %>% ggplot(aes(x = flowering, y = status, color = TGP)) + 
  geom_point() + 
  geom_smooth(method = 'glm')
# increased flowering increased death duh


# HEIGHTS-----------------------------------------------------------------------
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
  filter(!is.na(height_cm_7))

# avg height on measurement day 50 by treatment class - model
height_tidy_50 <- height_tidy |> 
  filter(height_cm_50!= 0 )

hist(height_tidy_32$height_cm_50)
mean(height_tidy_32$height_cm_50)
var(height_tidy_32$height_cm_50)
#poisson distribution with continous data :/
#log did not work :|
#google says forge ahead. okkkkk

height_mod_1 <- lmer(height_cm_50 ~ TGP + (1|tray), data = height_tidy_50)
summary(height_mod_1)

height_mod_2 <- lm(height_cm_50 ~ TGP , data = height_tidy_50)
summary(height_mod_3)

AIC(height_mod_1)
AIC(height_mod_2)
#height mod 1 is best

check_model(height_mod_3)
anova(height_mod_1, height_mod_3)
#model 1 is the best!

height_mod_1 <- lmer(height_cm_50 ~ TGP + (1|tray) , data = height_tidy_50)
summary(height_mod_1)

height_mod_1 |> emmeans(pairwise ~ TGP) 

flower_emmeans_result$contrasts |> 
  summary(infer = TRUE)
flower_emmeans_result$emmeans |> 
  as.data.frame()
head(flower_emmeans_result)

# avg height on measurement day 32 by treatment class - visual
height_tidy_32 |> 
  #group_by(TGP) |> 
  #summarise(mean_height_32 = mean(height_cm_32)) |> 
  ggplot(aes(x = TGP, y = height_cm_50)) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.2), alpha = 0.5) +
  xlab(" ") +
  ylab("Height (cm)") +
  ggtitle(" ")

plot_model(height_mod_1, 
           vline.color = "black", 
           show.values = TRUE, 
           value.offset = .2, 
           colors = "bw",
           title = "Flowering")

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

# BIOMASS-----------------------------------------------------------------------
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
  select(-`mass g minus 5 lvs`, -`5 lvs mass`)

biomass_tidy_1 <- biomass_tidy |> 
  pivot_wider(names_from = plant_part, values_from = total_mass_g) |> 
  mutate(total_biomass = root + shoot) |> 
  mutate(ratio_RS = root/shoot)

#i think i should log the RS ratio
biomass_tidy_1 |>  
  ggplot(aes(x = log(ratio_RS))) + 
  geom_density()

biomass_tidy_2 <- biomass_tidy_1 |> mutate(log_ratio_RS = log(ratio_RS))

#i think i should get rid of the outliers
boxplot(biomass_tidy_2$log_ratio_RS, plot=FALSE)$out
outliers <- boxplot(biomass_tidy_2$log_ratio_RS, plot=FALSE)$out
x <- biomass_tidy_2
x <- x[-which(x$log_ratio_RS %in% outliers),]

#models!
#total biomass
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

#root to shoot
ratio_biomass_mod_1 <- lmer(log_ratio_RS ~ TGP + (1|tray), data = x)
summary(ratio_biomass_mod_1)

ratio_biomass_mod_2 <- lm(log_ratio_RS ~ TGP , data = x)
summary(ratio_biomass_mod_2)

AIC(ratio_biomass_mod_1)
AIC(ratio_biomass_mod_2)
#mod 2 it is

check_model(ratio_biomass_mod_2)
anova(ratio_biomass_mod_1, ratio_biomass_mod_2, test="Chisq")
#model 1 is a little better but im sticking with model 2. ask rachel why

ratio_biomass_mod_2 |> emmeans(pairwise ~ TGP) 

# visualization

# total biomass
biomass_tidy_2  |> 
  ggplot(aes(x = TGP, y = log_ratio_RS)) + geom_violin(aes(color = TGP))

biomass_tidy_2 |> 
  #group_by(TGP) |> 
  ggplot(aes(x = TGP, y = total_biomass)) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.2), alpha = 0.5)

x |> 
  #group_by(TGP) |> 
  ggplot(aes(x = TGP, y = log_ratio_RS)) +
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.2), alpha = 0.5)
