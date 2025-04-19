# title: 08_plasticity_by_pop

# NEW - calculate a value that captures TGP response (RDPI of CC-DD)
# and see if there is a relationship btwn TGP plasticity and adaptivity
# by regressing it against survival, growth, and fecundity
# also, here, i break it up by POPULATION to see if the correlation btween
# CC-DD and trait is stronger based on different populations

# author: madeleine wallace

# PACKAGES -----------------------------------------------------------------
library(tidyverse)
#devtools::install_github("ameztegui/Plasticity")
library(Plasticity) # calculate rdpi, remember to cite, thx ameztegui
library(agricolae) # dependency
library(psych) # dependency
library(dplyr) # dependency
library(ggplot2) # dependency
library(sciplot) # dependency
library(esquisse) # build plots
source('code-and-data/scripts/plpa-manuscript/01_clean.R') # data


# NEW ANALYSIS - RDPI ADAPTIVE --------------------------------------------
# RDPI of CC-DD total biomass ---------------------------------------------
#first, combine total biomass dataframe and seed num dataframe
total_plastic <- SEED_FINAL |> 
  left_join(BIOMASS_FINAL, by = c("tray", "tray_ab", "id", "pop_id", "pop", "2021_treat", "2023_treat", "TGP", "trial", "ot", "pt"))

mort_day50 <- mort_day50 |> 
  #rename(TGP = TGP.x) |> 
  #rename(trial = trial.x) |> 
  #rename(vpd = vpd.x)

colnames(mort_day50)

total_plastic <- total_plastic %>%
  left_join(mort_day50 %>% select(tray, tray_ab, id, pop_id, pop, `2021_treat`, `2023_treat`, TGP, trial, ot, pt, status), 
            by = c("tray", "tray_ab", "id", "pop_id", "pop", "2021_treat", "2023_treat", "TGP", "trial", "ot", "pt"))

total_plastic <- total_plastic |> 
  #rename(TGP = TGP.x) |> 
  #rename(trial = trial.x) |> 
  rename(vpd = vpd.x)

# calculate avg seed num per VPD group and add it to dataframe
avg_fitness <- total_plastic |> 
  group_by(vpd) |> 
  summarise(avg_seeds = mean(num_total, na.rm = TRUE))

total_plastic <- total_plastic |> 
  group_by(vpd) |>
  mutate(avg_seed_num = mean(num_total, na.rm = TRUE))

# calculate RDPI CC-DD of total biomass per vpd
tgp_biomass <- BIOMASS_FINAL |> 
  filter(TGP %in% c('DD', 'CC')) |> 
  na.omit(total_biomass) |> 
  mutate(vpd = as.factor(vpd)) |> 
  mutate(total_biomass = as.numeric(total_biomass)) |> 
  mutate(TGP = as.factor(TGP))

# rdpi for total biomass
tgp_total_rdpi1 <- rdpi(tgp_biomass, sp = vpd, trait = total_biomass, factor = TGP)

avg_rdpi <- tgp_total_rdpi1 |> 
  group_by(sp) |> 
  summarise(avg_rdpi = mean(rdpi, na.rm = TRUE)) |> 
  rename(vpd = sp)

# filtered biomass_final with only cc and dd
total_plastic <- total_plastic |> 
  filter(TGP %in% c('DD', 'CC'))

# attach avg RDPI to vpd values in filtered biomass dataframe
total_plastic$vpd <- as.numeric(as.character(total_plastic$vpd))
avg_rdpi$vpd <- as.numeric(as.character(avg_rdpi$vpd))

total_plastic <- total_plastic %>%
  left_join(avg_rdpi, by = "vpd")

# now, i think i can regress RDPI and seed num....
cor.test(total_plastic$num_total, total_plastic$avg_rdpi, use = "complete.obs", method = "pearson")
cor.test(total_plastic$mass_total, total_plastic$avg_rdpi, use = "complete.obs", method = "pearson")
cor.test(total_plastic$status, total_plastic$avg_rdpi, use = "complete.obs", method = "pearson")

# summary
#total_plastic_summary <- total_plastic |> 
#group_by(TGP) |> 
#summarise(mean_cv_total = mean(cv_total, na.rm = TRUE),
#se_cv_total = sd(cv_total, na.rm = TRUE) / sqrt(n()),
#mean_fitness_robustness = mean(fitness_robustness, na.rm = TRUE),
#se_fitness_robustness = sd(fitness_robustness, na.rm = TRUE) / sqrt(n()))


# RDPI of CC-DD R:S ratio ---------------------------------------------
#first, combine total biomass dataframe and seed num dataframe
rs_plastic <- SEED_FINAL |> 
  left_join(BIOMASS_FINAL, by = c("tray", "tray_ab", "id", "pop_id", "pop", "2021_treat", "2023_treat", "TGP", "trial", "ot", "pt"))

rs_plastic <- rs_plastic %>%
  rename_with(~ gsub("\\.x$", "", .))

mort_day50 <- mort_day50 |> 
  #rename(TGP = TGP.x) |> 
  rename(trial = trial.x)

colnames(mort_day50)

rs_plastic <- rs_plastic %>%
  left_join(mort_day50 %>% select(tray, tray_ab, id, pop_id, pop, `2021_treat`, `2023_treat`, TGP, trial, ot, pt, status), 
            by = c("tray", "tray_ab", "id", "pop_id", "pop", "2021_treat", "2023_treat", "TGP", "trial", "ot", "pt"))

# calculate RDPI CC-DD of total biomass per vpd
tgp_rs <- BIOMASS_FINAL |> 
  filter(TGP %in% c('DD', 'CC')) |> 
  na.omit(ratio_RS) |> 
  mutate(vpd = as.factor(vpd)) |> 
  mutate(ratio_RS = as.numeric(ratio_RS)) |> 
  mutate(TGP = as.factor(TGP))

# rdpi for total biomass
tgp_rs_rdpi1 <- rdpi(tgp_rs, sp = vpd, trait = ratio_RS, factor = TGP)

avg_rdpi <- tgp_rs_rdpi1 |> 
  group_by(sp) |> 
  summarise(avg_rdpi = mean(rdpi, na.rm = TRUE)) |> 
  rename(vpd = sp)

# filtered biomass_final with only cc and dd
rs_plastic <- rs_plastic |> 
  filter(TGP %in% c('DD', 'CC'))

# attach avg RDPI to vpd values in filtered biomass dataframe
rs_plastic$vpd <- as.numeric(as.character(rs_plastic$vpd))
avg_rdpi$vpd <- as.numeric(as.character(avg_rdpi$vpd))

rs_plastic <- rs_plastic %>%
  left_join(avg_rdpi, by = "vpd")

# now, i think i can regress RDPI and seed num....
cor.test(rs_plastic$num_total, rs_plastic$avg_rdpi, use = "complete.obs", method = "pearson")
cor.test(rs_plastic$mass_total, rs_plastic$avg_rdpi, use = "complete.obs", method = "pearson")
cor.test(rs_plastic$status, rs_plastic$avg_rdpi, use = "complete.obs", method = "pearson")

# summary
rs_plastic_summary <- rs_plastic |> 
  group_by(vpd) |> 
  summarise(mean_mort = mean(status, na.rm = TRUE),
            mean_rdpi = mean(avg_rdpi, na.rm = TRUE))

#rsratio_rdpi <- 
  
  ggplot(rs_plastic_summary, aes(x = mean_mort, y = mean_rdpi)) +
  geom_point(size = 3.5, alpha = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "black", alpha = 0.3) +
  labs(x = "Proportion Alive", y = "RDPI CC-DD R:S ratio") +
  #scale_color_manual(values = c("CC" = "blue3",
  #"CD" = "chartreuse4",
  #"DC" = "orange1",
  #"DD" = "red3")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(color = "black", size = 10, face = "bold"),
        axis.title.y = element_text(color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        legend.title = element_text(color = "black", size = 10, face = "bold"),
        legend.text = element_text(color = "black", size = 10))
view(rsratio_rdpi)

ggsave("rsratio_rdpi.svg", rsratio_rdpi, path = "code-and-data/scripts/plpa-manuscript/figures",
       height = 6, width = 8, dpi = 300)

# RDPI of CC-DD shoot biomass ---------------------------------------------
#first, combine AG biomass dataframe and seed num dataframe
ag_plastic <- SEED_FINAL |> 
  left_join(BIOMASS_FINAL, by = c("tray", "tray_ab", "id", "pop_id", "pop", "2021_treat", "2023_treat", "TGP", "trial", "ot", "pt"))

ag_plastic <- ag_plastic %>%
  rename_with(~ gsub("\\.x$", "", .))

mort_day50 <- mort_day50 |> 
  #rename(TGP = TGP.x) |> 
  rename(trial = trial.x)

colnames(mort_day50)

ag_plastic <- ag_plastic %>%
  left_join(mort_day50 %>% select(tray, tray_ab, id, pop_id, pop, `2021_treat`, `2023_treat`, TGP, trial, ot, pt, status), 
            by = c("tray", "tray_ab", "id", "pop_id", "pop", "2021_treat", "2023_treat", "TGP", "trial", "ot", "pt"))

# calculate avg seed num per VPD group and add it to dataframe
avg_fitness <- ag_plastic |> 
  group_by(vpd) |> 
  summarise(avg_seeds = mean(num_total, na.rm = TRUE))

ag_plastic <- ag_plastic |> 
  group_by(vpd) |>
  mutate(avg_seed_num = mean(num_total, na.rm = TRUE))

# calculate RDPI CC-DD of total biomass per vpd
tgp_ag <- BIOMASS_FINAL |> 
  filter(TGP %in% c('DD', 'CC')) |> 
  na.omit(shoot) |> 
  mutate(vpd = as.factor(vpd)) |> 
  mutate(shoot = as.numeric(shoot)) |> 
  mutate(TGP = as.factor(TGP))

# rdpi for total biomass
tgp_ag_rdpi1 <- rdpi(tgp_rs, sp = vpd, trait = shoot, factor = TGP)

avg_rdpi <- tgp_ag_rdpi1 |> 
  group_by(sp) |> 
  summarise(avg_rdpi = mean(rdpi, na.rm = TRUE)) |> 
  rename(vpd = sp)

# filtered biomass_final with only cc and dd
ag_plastic <- ag_plastic |> 
  filter(TGP %in% c('DD', 'CC'))

# attach avg RDPI to vpd values in filtered biomass dataframe
ag_plastic$vpd <- as.numeric(as.character(ag_plastic$vpd))
avg_rdpi$vpd <- as.numeric(as.character(avg_rdpi$vpd))

ag_plastic <- ag_plastic %>%
  left_join(avg_rdpi, by = "vpd")

# now, i think i can regress RDPI and seed num....
cor.test(ag_plastic$num_total, ag_plastic$avg_rdpi, use = "complete.obs", method = "pearson")
cor.test(ag_plastic$mass_total, ag_plastic$avg_rdpi, use = "complete.obs", method = "pearson")
cor.test(ag_plastic$status, ag_plastic$avg_rdpi, use = "complete.obs", method = "pearson")