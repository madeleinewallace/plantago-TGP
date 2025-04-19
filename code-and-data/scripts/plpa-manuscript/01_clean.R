# title: 01_clean

# about: load, clean, and merge data for models in 02/03_models. 
# transformations done in this script when obviously called for

# author: madeleine wallace

# PACKAGES -----------------------------------------------------------------
library(tidyverse)
library(readxl) #read xlsx
library(ggplot2) #visualize data
library(ggpubr) #visualize data
library(rstatix) #check data

# CLIMATE -----------------------------------------------------------------
climate_untidy <- read_excel("code-and-data/data/2021-climate-info-MW.xlsx", sheet = 5)

climate_tidy <- climate_untidy |> 
  mutate(population = as.factor(population))

# BIOMASS ------------------------------------------------------------------
biomass_untidy <- read_excel("code-and-data/data/2023-greenhouse-data_biomass.xlsx")

biomass_untidy <- biomass_untidy |> 
  mutate(total_mass_g = ifelse(total_mass_g == 0, NA, total_mass_g))

biomass_tidy <- biomass_untidy |> 
  mutate(tray = as.factor(tray)) |>
  mutate(population = as.factor(population)) |> 
  rename(pop = population) |> 
  mutate(id = as.factor(id)) |>
  mutate(TGP = as.factor(TGP)) |> 
  mutate(`2021_treat` = as.factor(`2021_treat`)) |> 
  mutate(`2023_treat` = as.factor(`2023_treat`)) |> 
  mutate(total_mass_g = as.numeric(total_mass_g)) |> 
  dplyr::select(-`mass_g_minus_SLA`) |> 
  dplyr::select(-`mass_g_SLA`)
summary(biomass_tidy)

BIOMASS_FINAL <- biomass_tidy |> 
  pivot_wider(names_from = plant_part, values_from = total_mass_g) |> 
  mutate(total_biomass = root + shoot) |> 
  mutate(ratio_RS = root/shoot) |> 
  mutate(log_root = log(root)) |> 
  mutate(log_shoot = log(shoot)) |> 
  mutate(log_total = log(total_biomass)) |>
  mutate(log_RS = log(ratio_RS)) |> 
  mutate(ot = `2023_treat`) |> 
  mutate(pt = `2021_treat`) |> 
  filter(ratio_RS < 150)

BIOMASS_FINAL <- inner_join(BIOMASS_FINAL, climate_tidy, by = c("pop" = "population")) 
head(BIOMASS_FINAL)

BIOMASS_FINAL <- BIOMASS_FINAL |> 
  mutate(pop = as.factor(pop), 
         tray = as.factor(tray),
         id = as.factor(id))

head(BIOMASS_FINAL)

# FLOWERING ----------------------------------------------------------------
flower_untidy <- read_excel("code-and-data/data/2023-greenhouse-data_flowering.xlsx")

flower_tidy <- flower_untidy |> 
  mutate(tray = as.factor(tray)) |>
  mutate(pop = as.factor(pop)) |> 
  mutate(id = as.factor(id)) |>
  mutate(TGP = as.factor(TGP)) |> 
  mutate(`2021_treat` = as.factor(`2021_treat`)) |> 
  mutate(`2023_treat` = as.factor(`2023_treat`)) |> 
  mutate(flowering = as.character(flowering)) |> 
  mutate(`# of structures` = as.numeric(`# of structures`)) |> 
  rename( num_structure = `# of structures`)|> 
  mutate(ot = `2023_treat`) |> 
  mutate(pt = `2021_treat`)

FLOWER_FINAL <- flower_tidy |> 
  mutate(days_to_flower = as.numeric(difftime(ymd(`date flowering`), 
                                          ymd(date_germ), units = "days")))

FLOWER_FINAL <- inner_join(FLOWER_FINAL, climate_tidy, by = c("pop" = "population")) 

FLOWER_FINAL <- FLOWER_FINAL |> 
  mutate(pop = as.factor(pop), 
         tray = as.factor(tray),
         id = as.factor(id))

head(FLOWER_FINAL)

FLOWER_FINAL <- FLOWER_FINAL |> 
  mutate(num_structure = ifelse(!is.na(date_germ) & is.na(`date flowering`), 0, num_structure)) |> 
  mutate(pop = as.factor(pop))

## did the plant flower or not?
FLOWER_STATUS <- flower_tidy |> 
  filter(!is.na(date_germ)) |> 
  mutate(status = ifelse(!is.na(flowering) & flowering == "Y", 1, 0)) |> 
  filter(!(tray == 9 & tray_ab == "B" & id == "J2"))

FLOWER_STATUS <- inner_join(FLOWER_STATUS, climate_tidy, by = c("pop" = "population")) 

FLOWER_FINAL <- FLOWER_FINAL |> 
  filter(!(tray == 9 & tray_ab == "B" & id == "J2"))


# GROWTH ------------------------------------------------------------------
height_untidy <- read_excel("code-and-data/data/2023-greenhouse-data_heights.xlsx")

height_tidy <- height_untidy |> 
  mutate(tray = as.factor(tray)) |>
  mutate(pop = as.factor(pop)) |> 
  mutate(id = as.factor(id)) |> 
  mutate(TGP = as.factor(TGP)) |>
  mutate(`2021_treat` = as.factor(`2021_treat`)) |> 
  mutate(`2023_treat` = as.factor(`2023_treat`)) |> 
  mutate(height_cm_7 = as.numeric(height_cm_7)) |> 
  mutate(height_cm_14 = as.numeric(height_cm_14)) |> 
  mutate(height_cm_32 = as.numeric(height_cm_32)) |> 
  mutate(height_cm_50 = as.numeric(height_cm_50)) |> 
  mutate(height_cm_68 = as.numeric(height_cm_68)) |> 
  mutate(height_cm_86 = as.numeric(height_cm_86)) |> 
  mutate(height_cm_104 = as.numeric(height_cm_104)) |> 
  mutate(height_cm_122 = as.numeric(height_cm_122))|> 
  mutate(ot = `2023_treat`) |> 
  mutate(pt = `2021_treat`)

height_RGR <- height_tidy |> 
  mutate(height_cm_7 = na_if(height_cm_7, 0)) |> 
  mutate(height_cm_14 = na_if(height_cm_14, 0)) |> 
  mutate(height_cm_32 = na_if(height_cm_32, 0)) |> 
  mutate(height_cm_50 = na_if(height_cm_50, 0)) |> 
  mutate(height_cm_68 = na_if(height_cm_68, 0)) |> 
  mutate(height_cm_86 = na_if(height_cm_86, 0)) |> 
  mutate(height_cm_104 = na_if(height_cm_104, 0)) |> 
  mutate(height_cm_122 = na_if(height_cm_122, 0))

height_RGR <- height_RGR |> 
  rowwise() |> 
  mutate(max = max(c_across(height_cm_7:height_cm_122), na.rm = TRUE))

height_RGR <- height_RGR |> 
  rowwise() |> 
  mutate(min = min(c_across(height_cm_7:height_cm_122), na.rm = TRUE))

#calculate final RGR
RGR_FINAL <- height_RGR |> 
  mutate(ratio = max/min) |> 
  mutate(RGR = log(ratio)) |> 
  mutate(across(c(max, min, ratio, RGR), ~ ifelse(is.nan(.), NA, .))) |> 
  mutate(across(c(max, min), ~ ifelse(. == -Inf, NA, .))) |> 
  mutate(across(c(max, min), ~ ifelse(. == Inf, NA, .)))

RGR_FINAL <- inner_join(RGR_FINAL, climate_tidy, by = c("pop" = "population")) 

RGR_FINAL <- RGR_FINAL |> 
  mutate(pop = as.factor(pop), 
         tray = as.factor(tray),
         id = as.factor(id))

head(RGR_FINAL)

# final height and add climate data
HEIGHT_FINAL <- height_tidy |> 
  rowwise() |> 
  mutate(max = max(c_across(height_cm_7:height_cm_122), na.rm = T)) |> 
  mutate(max = na_if(max, -Inf)) |> 
  mutate(max = na_if(max, 0))

HEIGHT_FINAL <- inner_join(HEIGHT_FINAL, climate_tidy, by = c("pop" = "population")) 

HEIGHT_FINAL <- HEIGHT_FINAL |> 
  mutate(pop = as.factor(pop), 
         tray = as.factor(tray),
         id = as.factor(id))

head(HEIGHT_FINAL)

# MORTALITY ---------------------------------------------------------------
mort_untidy <- read_excel("code-and-data/data/2023-greenhouse-data_mortality.xlsx")

mort_tidy <- mort_untidy |> 
  mutate(tray = as.factor(tray)) |>
  mutate(id = as.factor(id)) |> 
  mutate(pop = as.factor(pop)) |> 
  mutate(TGP = as.factor(TGP)) |> 
  mutate(`2021_treat` = as.factor(`2021_treat`)) |> 
  mutate(`2023_treat` = as.factor(`2023_treat`)) |> 
  filter(trial == 2) |>
  pivot_longer(june_8:october_19, 
               names_to = 'date',
               values_to = 'status') |> 
  mutate(status = if_else(status == 'DEAD', 0, 1))|> # 0 means dead, 1 means alive
  mutate(ot = `2023_treat`) |> 
  mutate(pt = `2021_treat`)

# july 13th mortality, before most senescence, capturing treatment effects
MORT_JUNE8_FINAL <- mort_tidy |>  
  filter(date == 'june_8') |> 
  mutate(status = as.numeric(status))

# trying with height data on day 50, before senescence
mort_untidy$tray <- as.factor(mort_untidy$tray)
mort_untidy$pop <- as.factor(mort_untidy$pop)

MORT_DAY50 <- height_tidy |> 
  pivot_longer(height_cm_7:height_cm_122, 
               names_to = 'date',
               values_to = 'height') |> 
  filter(date == 'height_cm_50') |> 
  left_join(mort_untidy, by = c("tray", "tray_ab", "id", "pop_id", "pop", "2021_treat", "2023_treat", "TGP", "trial")) |> 
  dplyr::select(tray, tray_ab, id, pop_id, pop, `2021_treat`, `2023_treat`, TGP, trial, ot, pt, date, height, date_germ) |> 
  filter(!is.na(date_germ)) |> 
  mutate(height = replace_na(height, 0)) |> 
  mutate(status = if_else(height == 0, 0, 1)) |> 
  mutate(status = as.numeric(status)) |> 
  mutate(`2021_treat` = as.factor(`2021_treat`)) |> 
  mutate(`2023_treat` = as.factor(`2023_treat`))

MORT_DAY50 <- inner_join(MORT_DAY50, climate_tidy, by = c("pop" = "population")) 

MORT_DAY50 <- MORT_DAY50 |> 
  mutate(pop = as.factor(pop), 
         tray = as.factor(tray),
         id = as.factor(id))

head(MORT_DAY50)

# SLA & LDMC --------------------------------------------------------------
sla_ldmc_untidy <- read_excel("code-and-data/data/2023-greenhouse-data_sla-ldmc.xlsx")

sla_ldmc_tidy <- sla_ldmc_untidy |> 
  mutate_all(~ifelse(. %in% c("", "NA"), NA, .)) |> 
  #filter(!is.na(sla)) |> 
  mutate(tray = as.factor(tray)) |>
  mutate(id = as.factor(id)) |> 
  mutate(pop = as.factor(pop)) |> 
  mutate(TGP = as.factor(TGP)) |> 
  mutate(`2021_treat` = as.factor(`2021_treat`)) |> 
  mutate(`2023_treat` = as.factor(`2023_treat`)) |> 
  mutate(sla = as.numeric(sla)) |> 
  mutate(ldmc = as.numeric(ldmc)) |> 
  mutate(ot = `2023_treat`) |> 
  mutate(pt = `2021_treat`) |> 
  filter(sla > 12) |> 
  filter(sla < 300) |> 
  filter(ldmc > 0) |> 
  filter(ldmc < 1)

SLA_LDMC_FINAL <- sla_ldmc_tidy

SLA_LDMC_FINAL <- inner_join(SLA_LDMC_FINAL, climate_tidy, by = c("pop" = "population")) 

SLA_LDMC_FINAL <- SLA_LDMC_FINAL |> 
  mutate(pop = as.factor(pop), 
         tray = as.factor(tray),
         id = as.factor(id))

head(SLA_LDMC_FINAL)

# SEED MASS & NUMBER ------------------------------------------------------
seedf2_untidy <- read_excel("code-and-data/data/2023-greenhouse-data_f2maddiehannah-seed-data.xlsx")

seedf2_tidy <- seedf2_untidy |> 
  mutate(tray = as.factor(tray)) |>
  mutate(id = as.factor(id)) |> 
  mutate(pop = as.factor(pop)) |> 
  mutate(TGP = as.factor(TGP)) |> 
  mutate(`2021_treat` = as.factor(`2021_treat`)) |> 
  mutate(`2023_treat` = as.factor(`2023_treat`)) |> 
  mutate(`2023_mass_total_g` = ifelse(`2023_mass_total_g` == 0, NA, `2023_mass_total_g`)) |> 
  mutate(`2023_num_total` = ifelse(`2023_num_total` == "#DIV/0!", NA, `2023_num_total`)) |> 
  mutate(`2023_num_total` = as.numeric(`2023_num_total`))|> 
  mutate(ot = `2023_treat`) |> 
  mutate(pt = `2021_treat`) |> 
  mutate(num_total = `2023_num_total`) |> 
  mutate(mass_total = `2023_mass_total_g`) |> 
  mutate(across(c('num_total'), round, 0))

SEED_FINAL <- seedf2_tidy

SEED_FINAL <- inner_join(SEED_FINAL, climate_tidy, by = c("pop" = "population")) 

SEED_FINAL <- SEED_FINAL |> 
  mutate(pop = as.factor(pop), 
         tray = as.factor(tray),
         id = as.factor(id))

SEED_FINAL <- SEED_FINAL |> 
  filter(!(tray == 9 & tray_ab == "B" & id == "J2"))

head(SEED_FINAL)

# FINAL RESPONSE TRAITS ---------------------------------------------------
head(BIOMASS_FINAL) #AG biomass, BG biomass, total biomass, R:S ratio
head(FLOWER_FINAL) #days to flowering, number of seedheads
head(FLOWER_STATUS) #number flowering, flowering status
head(RGR_FINAL) #relative growth rate
head(HEIGHT_FINAL) #final height
head(MORT_DAY50) #mortality on july 13th, capturing treatment effects before senescense
head(SLA_LDMC_FINAL) #sla and ldmc
head(SEED_FINAL) #total seed mass and number


# climate stuff ----------------------------------------------------------------
climate_untidy1 <- read_csv("code-and-data/data/prism-climate-data/total-clim-data-1.csv")
climate_untidy2 <- read_csv("code-and-data/data/prism-climate-data/total-clim-data-2.csv")
climate_untidy3 <- read_csv("code-and-data/data/prism-climate-data/total-clim-data-3.csv")

total_climate <- bind_rows(climate_untidy1, climate_untidy2, climate_untidy3)
total_climate <- total_climate |> drop_na()
total_climate<- total_climate |> mutate(Date = as.character(Date))

total_climate <- total_climate |> 
                      mutate(date = ym(Date))

total_climate <- total_climate %>%
  mutate(date = as.Date(date))

results <- total_climate |> 
  mutate(
    year = year(date),
    month = month(date)) |> 
  group_by(Name, year) |> 
  summarise(
    annual_ppt = sum(`ppt (mm)`, na.rm = TRUE),
    annual_tmean = mean(`tmean (degrees C)`, na.rm = TRUE),
    spring_ppt = sum(`ppt (mm)`[month %in% c(4, 5, 6)], na.rm = TRUE),
    spring_tmean = mean(`tmean (degrees C)`[month %in% c(4, 5, 6)], na.rm = TRUE),
    spring_vpd = mean(`vpdmax (kPa)`[month %in% c(4, 5, 6)], na.rm = TRUE), .groups = "drop") |> 
  group_by(Name) |> 
  summarise(
    mean_annual_ppt = mean(annual_ppt, na.rm = TRUE),
    sd_annual_ppt = sd(annual_ppt, na.rm = TRUE),
    cv_annual_ppt = sd_annual_ppt / mean_annual_ppt * 100,
    
    mean_annual_tmean = mean(annual_tmean, na.rm = TRUE),
    sd_annual_tmean = sd(annual_tmean, na.rm = TRUE),
    cv_annual_tmean = sd_annual_tmean / mean_annual_tmean * 100,
    
    mean_spring_ppt = mean(spring_ppt, na.rm = TRUE),
    sd_spring_ppt = sd(spring_ppt, na.rm = TRUE),
    cv_spring_ppt = sd_spring_ppt / mean_spring_ppt * 100,
    
    mean_spring_tmean = mean(spring_tmean, na.rm = TRUE),
    sd_spring_tmean = sd(spring_tmean, na.rm = TRUE),
    cv_spring_tmean = sd_spring_tmean / mean_spring_tmean * 100,
    
    mean_spring_vpd = mean(spring_vpd, na.rm = TRUE),
    sd_spring_vpd = sd(spring_vpd, na.rm = TRUE),
    cv_spring_vpd = sd_spring_vpd / mean_spring_vpd * 100)

write.csv(results, "code-and-data/data/prism-climate-data/results.csv", row.names = FALSE)


