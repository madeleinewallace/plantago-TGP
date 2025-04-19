# title: 05_plasticity

# NEW - calculate a value that captures TGP response (RDPI of CC-DD)
# and see if there is a relationship btwn TGP plasticity and adaptivity
# by regressing it against survival, growth, and fecundity


# author: madeleine wallace

# PACKAGES -----------------------------------------------------------------
library(tidyverse)
library(dplyr)
#devtools::install_github("ameztegui/Plasticity")
library(Plasticity) # calculate rdpi, remember to cite, thx ameztegui
library(agricolae) # dependency
library(psych) # dependency
library(dplyr) # dependency
library(ggplot2) # dependency
#library(sciplot) # dependency
#library(esquisse) # build plots
source('code-and-data/scripts/plpa-manuscript/01_clean.R') # data


# ANALYSIS - RDPI ADAPTIVE --------------------------------------------
# RDPI of CC-DD RGR ---------------------------------------------
#first, combine RGR dataframe and seed num dataframe
RGR_FINAL <- RGR_FINAL |> 
  filter(!(tray == 9 & tray_ab == "B" & id == "J2"))
MORT_DAY50 <- MORT_DAY50 |> 
  filter(!(tray == 9 & tray_ab == "B" & id == "J2"))

rgr_plastic <- SEED_FINAL |> 
  left_join(RGR_FINAL, by = c("tray", "tray_ab", "id", "pop_id", "pop", "2021_treat", "2023_treat", "TGP", "trial", "ot", "pt"))

rgr_plastic <- rgr_plastic |> 
  #rename(TGP = TGP.x) |> 
  #rename(trial = trial.x) |> 
  rename(spring_vpd_cv = spring_vpd_cv.x)

colnames(MORT_DAY50)
colnames(rgr_plastic)

#total_plastic <- total_plastic %>%
#dplyr::select(-ends_with(".x"))  # Drop all columns ending with .x

rgr_plastic <- rgr_plastic %>%
  left_join(MORT_DAY50 %>% dplyr::select(tray, tray_ab, id, pop_id, pop, `2021_treat`, `2023_treat`, TGP, trial, ot, pt, status), 
            by = c("tray", "tray_ab", "id", "pop_id", "pop", "2021_treat", "2023_treat", "TGP", "trial", "ot", "pt"))

rgr_plastic <- rgr_plastic |> 
  #rename(TGP = TGP.x)
  #rename(trial = trial.x) |> 
  rename(spring_vpd_cv = spring_vpd_cv.x)

# calculate avg seed num per VPD group and add it to dataframe
avg_fitness <- rgr_plastic |> 
  group_by(spring_vpd_cv) |> 
  summarise(avg_seeds = mean(num_total, na.rm = TRUE))

rgr_plastic <- rgr_plastic |> 
  group_by(spring_vpd_cv) |>
  mutate(avg_seed_num = mean(num_total, na.rm = TRUE))

# calculate RDPI CC-DD of total biomass per vpd
tgp_rgr <- RGR_FINAL |> 
  filter(TGP %in% c('DD', 'CC')) |> 
  filter(!is.na(RGR)) |> 
  mutate(spring_vpd_cv = as.factor(spring_vpd_cv)) |> 
  mutate(RGR = as.numeric(RGR)) |> 
  mutate(TGP = as.factor(TGP))

# rdpi for total biomass
tgp_rgr_rdpi1 <- rdpi(tgp_rgr, sp = spring_vpd_cv, trait = RGR, factor = TGP)

avg_rdpi <- tgp_rgr_rdpi1 |> 
  group_by(sp) |> 
  summarise(avg_rdpi = mean(rdpi, na.rm = TRUE)) |> 
  rename(spring_vpd_cv = sp)

# filtered biomass_final with only cc and dd
rgr_plastic <- rgr_plastic |> 
  filter(TGP %in% c('DD', 'CC'))


# attach avg RDPI to vpd values in filtered biomass dataframe
rgr_plastic$spring_vpd_cv <- as.numeric(as.character(rgr_plastic$spring_vpd_cv))
avg_rdpi$spring_vpd_cv <- as.numeric(as.character(avg_rdpi$spring_vpd_cv))

rgr_plastic <- rgr_plastic %>%
  left_join(avg_rdpi, by = "spring_vpd_cv")

# adding in flowering status!
flower_status <- FLOWER_STATUS |> 
  rename(status_flower = status)

rgr_plastic <- rgr_plastic %>%
  left_join(flower_status %>% select(tray, tray_ab, id, pop_id, pop, `2021_treat`, `2023_treat`, TGP, trial, ot, pt, status_flower), 
            by = c("tray", "tray_ab", "id", "pop_id", "pop", "2021_treat", "2023_treat", "TGP", "trial", "ot", "pt"))

# now, i think i can regress RDPI and seed num....
cor.test(rgr_plastic$num_total, rgr_plastic$avg_rdpi, use = "complete.obs", method = "pearson")
cor.test(rgr_plastic$mass_total, rgr_plastic$avg_rdpi, use = "complete.obs", method = "pearson")
cor.test(rgr_plastic$status, rgr_plastic$avg_rdpi, use = "complete.obs", method = "pearson")
cor.test(rgr_plastic$status_flower, rgr_plastic$avg_rdpi, use = "complete.obs", method = "pearson")
cor.test(rgr_plastic$spring_vpd_cv, rgr_plastic$avg_rdpi, use = "complete.obs", method = "pearson")






# RDPI of CC-DD root biomass ---------------------------------------------
#first, combine root biomass dataframe and seed num dataframe
BIOMASS_FINAL <- BIOMASS_FINAL |> 
  filter(!(tray == 9 & tray_ab == "B" & id == "J2"))
MORT_DAY50 <- MORT_DAY50 |> 
  filter(!(tray == 9 & tray_ab == "B" & id == "J2"))

root_plastic <- SEED_FINAL |> 
  left_join(BIOMASS_FINAL, by = c("tray", "tray_ab", "id", "pop_id", "pop", "2021_treat", "2023_treat", "TGP", "trial", "ot", "pt"))

MORT_DAY50 <- MORT_DAY50 |> 
  #rename(TGP = TGP.x) |> 
  #rename(trial = trial.x) |> 
  rename(spring_vpd_cv = spring_vpd_cv.x)

colnames(MORT_DAY50)
colnames(root_plastic)

#total_plastic <- total_plastic %>%
#dplyr::select(-ends_with(".x"))  # Drop all columns ending with .x

root_plastic <- root_plastic %>%
  left_join(MORT_DAY50 %>% dplyr::select(tray, tray_ab, id, pop_id, pop, `2021_treat`, `2023_treat`, TGP, trial, ot, pt, status), 
            by = c("tray", "tray_ab", "id", "pop_id", "pop", "2021_treat", "2023_treat", "TGP", "trial", "ot", "pt"))

root_plastic <- root_plastic |> 
  #rename(TGP = TGP.x)
  #rename(trial = trial.x) |> 
  rename(spring_vpd_cv = spring_vpd_cv.x)

# calculate avg seed num per VPD group and add it to dataframe
avg_fitness <- root_plastic |> 
  group_by(spring_vpd_cv) |> 
  summarise(avg_seeds = mean(num_total, na.rm = TRUE))

root_plastic <- root_plastic |> 
  group_by(spring_vpd_cv) |>
  mutate(avg_seed_num = mean(num_total, na.rm = TRUE))

# calculate RDPI CC-DD of total biomass per vpd
tgp_root <- BIOMASS_FINAL |> 
  filter(TGP %in% c('DD', 'CC')) |> 
  na.omit(root) |> 
  mutate(spring_vpd_cv = as.factor(spring_vpd_cv)) |> 
  mutate(root = as.numeric(root)) |> 
  mutate(TGP = as.factor(TGP))

# rdpi for total biomass
tgp_root_rdpi1 <- rdpi(tgp_root, sp = spring_vpd_cv, trait = root, factor = TGP)

avg_rdpi <- tgp_root_rdpi1 |> 
  group_by(sp) |> 
  summarise(avg_rdpi = mean(rdpi, na.rm = TRUE)) |> 
  rename(spring_vpd_cv = sp)

# filtered biomass_final with only cc and dd
root_plastic <- root_plastic |> 
  filter(TGP %in% c('DD', 'CC'))


# attach avg RDPI to vpd values in filtered biomass dataframe
root_plastic$spring_vpd_cv <- as.numeric(as.character(root_plastic$spring_vpd_cv))
avg_rdpi$spring_vpd_cv <- as.numeric(as.character(avg_rdpi$spring_vpd_cv))

root_plastic <- root_plastic %>%
  left_join(avg_rdpi, by = "spring_vpd_cv")

# adding in flowering status!
flower_status <- FLOWER_STATUS |> 
  rename(status_flower = status)

root_plastic <- root_plastic %>%
  left_join(flower_status %>% select(tray, tray_ab, id, pop_id, pop, `2021_treat`, `2023_treat`, TGP, trial, ot, pt, status_flower), 
            by = c("tray", "tray_ab", "id", "pop_id", "pop", "2021_treat", "2023_treat", "TGP", "trial", "ot", "pt"))

# now, i think i can regress RDPI and seed num....
cor.test(root_plastic$num_total, root_plastic$avg_rdpi, use = "complete.obs", method = "pearson")
cor.test(root_plastic$mass_total, root_plastic$avg_rdpi, use = "complete.obs", method = "pearson")
cor.test(root_plastic$status, root_plastic$avg_rdpi, use = "complete.obs", method = "pearson")
cor.test(root_plastic$status_flower, root_plastic$avg_rdpi, use = "complete.obs", method = "pearson")
cor.test(root_plastic$spring_vpd_cv, root_plastic$avg_rdpi, use = "complete.obs", method = "pearson")


# FIG 5 RDPI - TRAIT -----------------------------------------------------------

# root biomass - flowering
# summary and viz of rs to flowering
root_plastic_summary2 <- root_plastic |> 
  group_by(spring_vpd_cv) |> 
  summarise(mean_mort = mean(status, na.rm = TRUE),
            mean_flower = mean(status_flower, na.rm = TRUE),
            mean_rdpi = mean(avg_rdpi, na.rm = TRUE))

root_flower_rdpi <- ggplot(root_plastic_summary2, aes(x = mean_flower, y = mean_rdpi)) +
  geom_point(size = 3.5, alpha = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "black", alpha = 0.3) +
  labs(x = "proportion flowering", y = "TGP of root biomass") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(color = "black", size = 10, face = "bold"),
        axis.title.y = element_text(color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        legend.title = element_text(color = "black", size = 10, face = "bold"),
        legend.text = element_text(color = "black", size = 10))

root_flower_rdpi

ggsave("root_flower_rdpi.svg", root_flower_rdpi, path = "code-and-data/scripts/plpa-manuscript/figures",
       dpi = 300, width = 3.25, height = 3.25)


# RGR - mortality
rgr_plastic_summary2 <- rgr_plastic |> 
  group_by(spring_vpd_cv) |> 
  summarise(mean_mort = mean(status, na.rm = TRUE),
            mean_flower = mean(status * n(), na.rm = TRUE),
            mean_rdpi = mean(avg_rdpi, na.rm = TRUE))

rgr_mort_rdpi <- ggplot(rgr_plastic_summary2, aes(x = mean_mort, y = mean_rdpi)) +
  geom_point(size = 3.5, alpha = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "black", alpha = 0.3) +
  labs(x = "proportion alive", y = "TGP of RGR") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(color = "black", size = 10, face = "bold"),
        axis.title.y = element_text(color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        legend.title = element_text(color = "black", size = 10, face = "bold"),
        legend.text = element_text(color = "black", size = 10))

rgr_mort_rdpi

ggsave("rgr_mort_rdpi.svg", rgr_mort_rdpi, path = "code-and-data/scripts/plpa-manuscript/figures",
       dpi = 300, width = 3.25, height = 3.25)


# root biomass - spring VPD CV
root_plastic_summary2 <- root_plastic |> 
  group_by(spring_vpd_cv) |> 
  summarise(mean_mort = mean(status, na.rm = TRUE),
            mean_flower = mean(status * n(), na.rm = TRUE),
            mean_rdpi = mean(avg_rdpi, na.rm = TRUE))

root_springvpd_rdpi <- ggplot(root_plastic_summary2, aes(x = spring_vpd_cv, y = mean_rdpi)) +
  geom_point(size = 3.5, alpha = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "black", alpha = 0.3) +
  labs(x = "spring VPD CV (%)", y = "TGP of root biomass") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(color = "black", size = 10, face = "bold"),
        axis.title.y = element_text(color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        legend.title = element_text(color = "black", size = 10, face = "bold"),
        legend.text = element_text(color = "black", size = 10))

root_springvpd_rdpi

ggsave("root_springvpd_rdpi.svg", root_springvpd_rdpi, path = "code-and-data/scripts/plpa-manuscript/figures",
       dpi = 300, width = 3.25, height = 3.25)
