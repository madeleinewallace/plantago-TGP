# title: 07_supplemental
# about: tables for manuscript
# author: madeleine wallace



# PACKAGES -----------------------------------------------------------------
library(tidyverse)
library(ggplot2) # graphs
library(patchwork) # collecting graphs
library(gridExtra)
library(grid)
library(gtable)
library(readxl)
source('code-and-data/scripts/plpa-manuscript/01_clean.R')

### trait stats ---------------------------------------------------------------
calculate <- function(df, trait_column) {
  df %>%
    group_by(TGP) %>%
    summarise(
      n = sum(!is.na(.data[[trait_column]])),
      mean = round(mean(.data[[trait_column]], na.rm = TRUE), 6),
      se = round(sd(.data[[trait_column]], na.rm = TRUE) / sqrt(n), 6)
    )
}

calculate(BIOMASS_FINAL, "root")
calculate(BIOMASS_FINAL, "shoot")
calculate(BIOMASS_FINAL, "total_biomass")  
calculate(BIOMASS_FINAL, "ratio_RS")
calculate(HEIGHT_FINAL, "max")
calculate(RGR_FINAL, "RGR")
calculate(SLA_LDMC_FINAL, "sla")
calculate(SLA_LDMC_FINAL, "ldmc")
calculate(mort_day50, "status")
calculate(FLOWER_STATUS, "status")
calculate(SEED_FINAL, "2023_mass_total_g")
calculate(SEED_FINAL, "2023_num_total")
calculate(FLOWER_FINAL, "num_structure")
calculate(FLOWER_FINAL, "days_to_flower")

calculate_vpd <- function(df, trait_column) {
  df %>%
    group_by(TGP, vpd) |> 
    summarise(
      n = sum(!is.na(.data[[trait_column]])),
      mean = round(mean(.data[[trait_column]], na.rm = TRUE), 6),
      se = round(sd(.data[[trait_column]], na.rm = TRUE) / sqrt(n), 6)
    )
}

calculate_flowered <- function(df, trait_column) {
  df %>%
    group_by(TGP) %>%
    summarise(
      total = n(),  # Total individuals in each treatment group
      flowered = sum(.data[[trait_column]], na.rm = TRUE),  # Sum of 1s = number of flowered plants
      proportion_flowered = round(flowered / total, 3)  # Proportion that flowered
    )
}

# Run the function
calculate_flowered(FLOWER_STATUS, "status")


################################# trait by vpd stats----------------------------

calculate_vpd_long <- function(df, trait_column) {
  df %>%
    group_by(spring_vpd_cv, TGP) %>% 
    summarise(
      mean = round(mean(.data[[trait_column]], na.rm = TRUE), 6),
      se = round(sd(.data[[trait_column]], na.rm = TRUE) / sqrt(sum(!is.na(.data[[trait_column]]))), 6),
      .groups = "drop"
    ) %>%
    arrange(spring_vpd_cv, TGP)
}

calculate(BIOMASS_FINAL, "root")
calculate(BIOMASS_FINAL, "shoot")
calculate(BIOMASS_FINAL, "total_biomass")  
calculate(BIOMASS_FINAL, "ratio_RS")
calculate(HEIGHT_FINAL, "max")
calculate(RGR_FINAL, "RGR")
calculate(SLA_LDMC_FINAL, "sla")
calculate(SLA_LDMC_FINAL, "ldmc")
calculate(mort_day50, "status")
calculate(flower_status, "status")
calculate(SEED_FINAL, "2023_mass_total_g")
calculate(SEED_FINAL, "2023_num_total")
calculate(FLOWER_FINAL, "num_structure")
calculate(FLOWER_FINAL, "days_to_flower")

root_output_long <- calculate_vpd_long(BIOMASS_FINAL, "root")
print(root_output_long)
write.table(root_output_long, pipe("pbcopy"), sep = "\t", row.names = FALSE)

shoot_output_long <- calculate_vpd_long(BIOMASS_FINAL, "shoot")
print(shoot_output_long)
write.table(shoot_output_long, pipe("pbcopy"), sep = "\t", row.names = FALSE)

total_output_long <- calculate_vpd_long(BIOMASS_FINAL, "total_biomass")
print(total_output_long)
write.table(total_output_long, pipe("pbcopy"), sep = "\t", row.names = FALSE)

rs_output_long <- calculate_vpd_long(BIOMASS_FINAL, "ratio_RS")
print(rs_output_long)
write.table(rs_output_long, pipe("pbcopy"), sep = "\t", row.names = FALSE)

max_output_long <- calculate_vpd_long(HEIGHT_FINAL, "max")
print(max_output_long)
write.table(max_output_long, pipe("pbcopy"), sep = "\t", row.names = FALSE)

RGR_output_long <- calculate_vpd_long(RGR_FINAL, "RGR")
print(RGR_output_long)
write.table(RGR_output_long, pipe("pbcopy"), sep = "\t", row.names = FALSE)

sla_output_long <- calculate_vpd_long(SLA_LDMC_FINAL, "sla")
print(sla_output_long)
write.table(sla_output_long, pipe("pbcopy"), sep = "\t", row.names = FALSE)

ldmc_output_long <- calculate_vpd_long(SLA_LDMC_FINAL, "ldmc")
print(ldmc_output_long)
write.table(ldmc_output_long, pipe("pbcopy"), sep = "\t", row.names = FALSE)

mort_output_long <- calculate_vpd_long(MORT_DAY50, "status")
print(mort_output_long)
write.table(mort_output_long, pipe("pbcopy"), sep = "\t", row.names = FALSE)

flowerstat_output_long <- calculate_vpd_long(FLOWER_STATUS, "status")
print(flowerstat_output_long)
write.table(flowerstat_output_long, pipe("pbcopy"), sep = "\t", row.names = FALSE)

seednum_output_long <- calculate_vpd_long(SEED_FINAL, "2023_num_total")
print(seednum_output_long)
write.table(seednum_output_long, pipe("pbcopy"), sep = "\t", row.names = FALSE)

seedmass_output_long <- calculate_vpd_long(SEED_FINAL, "2023_mass_total_g")
print(seedmass_output_long)
write.table(seedmass_output_long, pipe("pbcopy"), sep = "\t", row.names = FALSE)

daysflow_output_long <- calculate_vpd_long(FLOWER_FINAL, "days_to_flower")
print(daysflow_output_long)
write.table(daysflow_output_long, pipe("pbcopy"), sep = "\t", row.names = FALSE)

flowstruc_output_long <- calculate_vpd_long(FLOWER_FINAL, "num_structure")
print(flowstruc_output_long)
write.table(flowstruc_output_long, pipe("pbcopy"), sep = "\t", row.names = FALSE)


################################# mappie----------------------------------------

# packages
library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(prettymapr)


# load data
climate <- read_excel("code-and-data/data/2021-climate-info-MW.xlsx")

# convert to sf object
climate_sf <- st_as_sf(climate, coords = c("Long", "Lat"), crs = 4326)

# get U.S. map
us <- ne_states(country = "united states of america", returnclass = "sf")

# filter southwestern states
southwest <- us[us$name %in% c("California", "Nevada", "Arizona", "Texas", "New Mexico", "Utah", "Colorado", "Oklahoma", "Kansas", "Wyoming", "Idaho"), ]

# plot
map <- ggplot() +
  geom_sf(data = southwest, fill = "white", color = "black", size = 0.6) +
  geom_sf(data = climate_sf, color = "gray", alpha = 10, size = 8, shape = 16) +
  geom_sf_text(data = climate_sf,
               aes(label = population),
               color = "black", size = 16 / .pt, nudge_y = 0, nudge_x = 0, family = "Helvetica") +
  coord_sf(xlim = c(-118, -102), ylim = c(30, 43), expand = FALSE) +
  annotation_scale(location = "bl", width_hint = 0.3, text_family = "Helvetica",  line_width = 0.4) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.03, "in"), pad_y = unit(0.2, "in"),
                         style = north_arrow_fancy_orienteering) +
  theme_classic(base_family = "Helvetica", base_size = 16) +
  theme(axis.title = element_text(family = "Helvetica", size = 16),
    axis.text = element_text(family = "Helvetica", size = 16)) +
  xlab("Longitude") + ylab("Latitude")

map

# save

ggsave("map.svg", plot = map, path = "code-and-data/scripts/plpa-manuscript/figures",
       dpi = 600)



################################# seed size difference--------------------------
library(readxl)

seedmass <- read_excel("code-and-data/data/2021-greenhouse-data_f1zoe-seed-data.xlsx", 
                                                    sheet = "weight of 25 seeds")

ggplot(seedmass, aes(x = `2021_treatment`, y = `2021_mass_25_g`, fill = `2021_treatment`)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Seed Mass by Treatment", y = "Mass (g)", x = "Treatment")

t.test(`2021_mass_25_g` ~ `2021_treatment`, data = seedmass)
shapiro.test(seedmass$`2021_mass_25_g`[seedmass$`2021_treatment` == "ambi"])
shapiro.test(seedmass$`2021_mass_25_g`[seedmass$`2021_treatment` == "water"])

library(lme4)
model <- lmer(`2021_mass_25_g` ~ `2021_treatment` + (1|population), data = seedmass)
summary(model)
