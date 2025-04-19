# title: 04_visualization
# about: visualizing results for manuscript
# author: madeleine wallace



# PACKAGES -----------------------------------------------------------------
library(tidyverse)
library(ggplot2) # graphs
library(patchwork) # collecting graphs
library(svglite)
library(ggpattern)
library(cowplot)
source('code-and-data/scripts/plpa-manuscript/01_clean.R')


# FIG 4 GRAPHS WITH POPULATION LEVEL DATA ----------------------------------
# VPD - rootbiomass rxn norm ------------------------------------------------
root_na <- BIOMASS_FINAL |>
  filter(!is.na(ot) & !is.na(pt) & !is.na(root))

root_summary <- root_na |> 
  group_by(TGP, spring_vpd_cv) |> 
  summarise(mean_root = mean(root), 
            se_root = sd(root) / sqrt(n()), .groups = 'drop')



# drought offsrping treatment
root_drought <- ggplot(filter(root_summary, TGP %in% c("CD", "DD")), 
                                aes(x = spring_vpd_cv, y = mean_root, color = TGP)) +
  geom_point(size = 1.5, alpha = 0.8) +
  geom_smooth(aes(color = TGP, linetype = TGP), method = "lm", se = TRUE, size = 1.2, alpha = 0.15) +
  labs(title = "offspring treatment = drought", x = "spring VPD CV (%)", y = "root biomass (g)", color = "treatment group", shape = "treatment group") + 
  scale_color_manual(values = c("DD" = "#A32B6B", "CD" = "#C7B245")) +
  #scale_shape_manual(values = c("DD" = 17, "CD" = 16)) +
  scale_linetype_manual(values = c("DD" = "dotted", "CD" = "dashed")) +
  ylim(0.1, 0.7) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(color = "black", size = 10, face = "bold"),
        axis.title.y = element_text(color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        legend.title = element_text(color = "black", size = 10, face = "bold"),
        legend.text = element_text(color = "black", size = 10),
        plot.title = element_text(color = "black", size = 10, face = "bold", hjust = 0.5),
        legend.position = "right") +
  guides(color = guide_legend(override.aes = list(linetype = c("dotted", "dashed"))),
         linetype = guide_legend(override.aes = list(color = c("#A32B6B", "#C7B245"))))
root_drought

# control offspring treatment
root_control <- ggplot(filter(root_summary, TGP %in% c("DC", "CC")), 
                       aes(x = spring_vpd_cv, y = mean_root, color = TGP)) +
  geom_point(size = 1.5, alpha = 0.8) +
  geom_smooth(aes(color = TGP, linetype = TGP), method = "lm", se = TRUE, size = 1.2, alpha = 0.15) +
  labs(title = "offspring treatment = control",  
       x = "spring VPD CV (%)", 
       y = "root biomass (g)", 
       color = "treatment group", 
       shape = "treatment group") +  
  scale_color_manual(values = c("CC" = "#386C9E", "DC" = "#177D31")) +
  #scale_shape_manual(values = c("CC" = 16, "DC" = 17)) +
  scale_linetype_manual(values = c("CC" = "solid", "DC" = "dotdash")) +
  ylim(0.1, 0.7) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(color = "black", size = 10, face = "bold"),
        axis.title.y = element_text(color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        legend.title = element_text(color = "black", size = 10, face = "bold"),
        legend.text = element_text(color = "black", size = 10),
        plot.title = element_text(color = "black", size = 10, face = "bold", hjust = 0.5),
        legend.position = "right") +
  guides(color = guide_legend(override.aes = list(linetype = c("solid", "dotdash"))),
         linetype = guide_legend(override.aes = list(color = c("#386C9E", "#177D31"))))

root_control

legend_root_drought <- get_legend(root_drought)
legend_root_control <- get_legend(root_control)

ggsave("legend_root_drought.svg", plot = legend_root_drought, path = "code-and-data/scripts/plpa-manuscript/figures",
       dpi = 300, width = 2, height = 2)
ggsave("legend_root_control.svg", plot = legend_root_control, path = "code-and-data/scripts/plpa-manuscript/figures",
       dpi = 300, width = 2, height = 2)

ggsave("root_drought.svg", root_drought, path = "code-and-data/scripts/plpa-manuscript/figures",
       dpi = 300, width = 3, height = 3)
ggsave("root_control.svg", root_control, path = "code-and-data/scripts/plpa-manuscript/figures",
       dpi = 300, width = 3, height = 3)


# VPD - total biomass rxn norm ------------------------------------------------
total_na <- BIOMASS_FINAL |>
  filter(!is.na(ot) & !is.na(pt) & !is.na(total_biomass))

total_summary <- total_na |> 
  group_by(TGP, spring_vpd_cv) |> 
  summarise(mean_total = mean(total_biomass), 
            se_total = sd(total_biomass) / sqrt(n()), .groups = 'drop')

# drought parental treatment
total_drought <- ggplot(filter(total_summary, TGP %in% c("DD", "CD")), 
                       aes(x = spring_vpd_cv, y = mean_total, color = TGP, shape = TGP)) +
  geom_point(size = 3) +
  geom_smooth(aes(color = TGP), method = "lm", se = TRUE, size = 1.2, alpha = 0.15) +
  labs(title = "offspring treatment = drought", 
       x = "spring VPD CV (%)", 
       y = "total biomass (g)", 
       color = "treatment group", 
       shape = "treatment group") + 
  scale_color_manual(values = c("DD" = "#A32B6B", "CD" = "#C7B245")) +
  scale_shape_manual(values = c("DD" = 17, "CD" = 16)) +
  ylim(0.0, 1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(color = "black", size = 10, face = "bold"),
        axis.title.y = element_text(color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        legend.title = element_text(color = "black", size = 10, face = "bold"),
        legend.text = element_text(color = "black", size = 10),
        plot.title = element_text(color = "black", size = 10, face = "bold", hjust = 0.5),
        legend.position = "none")
total_drought

# control parental treatment
total_control <- ggplot(filter(total_summary, TGP %in% c("DC", "CC")), 
                       aes(x = spring_vpd_cv, y = mean_total, color = TGP, shape = TGP)) +
  geom_point(size = 3) +
  geom_smooth(aes(color = TGP), method = "lm", se = TRUE, size = 1.2, alpha = 0.15) +
  labs(title = "offspring treatment = control",  
       x = "spring VPD CV (%)", 
       y = "total biomass (g)", 
       color = "treatment group", 
       shape = "treatment group") +  
  scale_color_manual(values = c("CC" = "#386C9E", "DC" = "#177D31")) +
  scale_shape_manual(values = c("CC" = 16, "DC" = 17)) +
  ylim(0.0, 1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(color = "black", size = 10, face = "bold"),
        axis.title.y = element_text(color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        legend.title = element_text(color = "black", size = 10, face = "bold"),
        legend.text = element_text(color = "black", size = 10),
        plot.title = element_text(color = "black", size = 10, face = "bold", hjust = 0.5),
        legend.position = "none")

total_control

ggsave("total_drought.svg", total_drought, path = "code-and-data/scripts/plpa-manuscript/figures",
       dpi = 300, width = 3, height = 3)
ggsave("total_control.svg", total_control, path = "code-and-data/scripts/plpa-manuscript/figures",
       dpi = 300, width = 3, height = 3)


# VPD - RS ratio rxn norm ------------------------------------------------
rs_na <- BIOMASS_FINAL |>
  filter(!is.na(ot) & !is.na(pt) & !is.na(ratio_RS))

rs_summary <- rs_na |> 
  group_by(TGP, spring_vpd_cv) |> 
  summarise(mean_rs = mean(ratio_RS), 
            se_rs = sd(ratio_RS) / sqrt(n()), .groups = 'drop')

# drought parental treatment
rs_drought <- ggplot(filter(rs_summary, TGP %in% c("DD", "CD")), 
                        aes(x = spring_vpd_cv, y = mean_rs, color = TGP, shape = TGP)) +
  geom_point(size = 3) +
  geom_smooth(aes(color = TGP), method = "lm", se = TRUE, size = 1.2, alpha = 0.15) +
  labs(title = "offspring treatment = drought", 
       x = "spring VPD CV (%)", 
       y = "R:S ratio", 
       color = "treatment group", 
       shape = "treatment group") + 
  scale_color_manual(values = c("DD" = "#A32B6B", "CD" = "#C7B245")) +
  scale_shape_manual(values = c("DD" = 17, "CD" = 16)) +
  ylim(0.0, 10) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(color = "black", size = 10, face = "bold"),
        axis.title.y = element_text(color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        legend.title = element_text(color = "black", size = 10, face = "bold"),
        legend.text = element_text(color = "black", size = 10),
        plot.title = element_text(color = "black", size = 10, face = "bold", hjust = 0.5),
        legend.position = "none")
rs_drought

# control parental treatment
rs_control <- ggplot(filter(rs_summary, TGP %in% c("DC", "CC")), 
                        aes(x = spring_vpd_cv, y = mean_rs, color = TGP, shape = TGP)) +
  geom_point(size = 3) +
  geom_smooth(aes(color = TGP), method = "lm", se = TRUE, size = 1.2, alpha = 0.15) +
  labs(title = "offspring treatment = control",  
       x = "spring VPD CV (%)", 
       y = "R:S ratio", 
       color = "treatment group", 
       shape = "treatment group") +  
  scale_color_manual(values = c("CC" = "#386C9E", "DC" = "#177D31")) +
  scale_shape_manual(values = c("CC" = 16, "DC" = 17)) +
  ylim(0.0, 10) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(color = "black", size = 10, face = "bold"),
        axis.title.y = element_text(color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        legend.title = element_text(color = "black", size = 10, face = "bold"),
        legend.text = element_text(color = "black", size = 10),
        plot.title = element_text(color = "black", size = 10, face = "bold", hjust = 0.5),
        legend.position = "none")

rs_control

ggsave("rs_drought.svg", rs_drought, path = "code-and-data/scripts/plpa-manuscript/figures",
       dpi = 300, width = 3, height = 3)
ggsave("rs_control.svg", rs_control, path = "code-and-data/scripts/plpa-manuscript/figures",
       dpi = 300, width = 3, height = 3)



# VPD - RGR rxn norm ------------------------------------------------
rgr_na <- RGR_FINAL |>
  filter(!is.na(ot) & !is.na(pt) & !is.na(RGR))

rgr_summary <- rgr_na |> 
  group_by(TGP, spring_vpd_cv) |> 
  summarise(mean_rgr = mean(RGR), 
            se_rgr = sd(RGR) / sqrt(n()), .groups = 'drop')

# drought offsrping treatment
rgr_drought <- ggplot(filter(rgr_summary, TGP %in% c("DD", "CD")), 
                     aes(x = spring_vpd_cv, y = mean_rgr, color = TGP)) +
  geom_point(size = 1.5, alpha = 0.8) +
  geom_smooth(aes(color = TGP, linetype = TGP), method = "lm", se = TRUE, size = 1.2, alpha = 0.15) +
  labs(title = "offspring treatment = drought", 
       x = "spring VPD CV (%)", 
       y = "RGR", 
       color = "treatment group", 
       shape = "treatment group") + 
  scale_color_manual(values = c("DD" = "#A32B6B", "CD" = "#C7B245")) +
  #scale_shape_manual(values = c("DD" = 17, "CD" = 16)) +
  scale_linetype_manual(values = c("DD" = "dotted", "CD" = "dashed")) +
  ylim(0.3, 1.8) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(color = "black", size = 10, face = "bold"),
        axis.title.y = element_text(color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        legend.title = element_text(color = "black", size = 10, face = "bold"),
        legend.text = element_text(color = "black", size = 10),
        plot.title = element_text(color = "black", size = 10, face = "bold", hjust = 0.5),
        legend.position = "none")
rgr_drought

# control offspring treatment
rgr_control <- ggplot(filter(rgr_summary, TGP %in% c("DC", "CC")), 
                     aes(x = spring_vpd_cv, y = mean_rgr, color = TGP)) +
  geom_point(size = 1.5, alpha = 0.8) +
  geom_smooth(aes(color = TGP, linetype = TGP), method = "lm", se = TRUE, size = 1.2, alpha = 0.15) +
  labs(title = "offspring treatment = control",  
       x = "spring VPD CV (%)", 
       y = "RGR", 
       color = "treatment group", 
       shape = "treatment group") +  
  scale_color_manual(values = c("CC" = "#386C9E", "DC" = "#177D31")) +
  #scale_shape_manual(values = c("CC" = 16, "DC" = 17)) +
  scale_linetype_manual(values = c("CC" = "solid", "DC" = "dotdash")) +
  ylim(0.3, 1.8) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(color = "black", size = 10, face = "bold"),
        axis.title.y = element_text(color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        legend.title = element_text(color = "black", size = 10, face = "bold"),
        legend.text = element_text(color = "black", size = 10),
        plot.title = element_text(color = "black", size = 10, face = "bold", hjust = 0.5),
        legend.position = "none")

rgr_control

ggsave("rgr_drought.svg", rgr_drought, path = "code-and-data/scripts/plpa-manuscript/figures",
       dpi = 300, width = 3, height = 3)
ggsave("rgr_control.svg", rgr_control, path = "code-and-data/scripts/plpa-manuscript/figures",
       dpi = 300, width = 3, height = 3)

# VPD - number flowered rxn norm ------------------------------------------------
flow_na <- FLOWER_STATUS |>
  filter(!is.na(ot) & !is.na(pt) & !is.na(status))

flow_summary <- flow_na |> 
  group_by(TGP, spring_vpd_cv) |> 
  summarise(mean_flow = mean(status * n()), 
            se_flow = sd(status * n()) / sqrt(n()), .groups = 'drop')

# drought parental treatment
flow_drought <- ggplot(filter(flow_summary, TGP %in% c("DD", "CD")), 
                      aes(x = spring_vpd_cv, y = mean_flow, color = TGP, shape = TGP)) +
  geom_point(size = 3) +
  geom_smooth(aes(color = TGP), method = "lm", se = TRUE, size = 1.2, alpha = 0.15) +
  labs(title = "offspring treatment = drought", 
       x = "spring VPD CV (%)", 
       y = "number flowered ", 
       color = "treatment group", 
       shape = "treatment group") + 
  scale_color_manual(values = c("DD" = "#A32B6B", "CD" = "#C7B245")) +
  scale_shape_manual(values = c("DD" = 17, "CD" = 16)) +
  ylim(-7, 35) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(color = "black", size = 10, face = "bold"),
        axis.title.y = element_text(color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        legend.title = element_text(color = "black", size = 10, face = "bold"),
        legend.text = element_text(color = "black", size = 10),
        plot.title = element_text(color = "black", size = 10, face = "bold", hjust = 0.5),
        legend.position = "none")
flow_drought

# control parental treatment
flow_control <- ggplot(filter(flow_summary, TGP %in% c("DC", "CC")), 
                      aes(x = spring_vpd_cv, y = mean_flow, color = TGP, shape = TGP)) +
  geom_point(size = 3) +
  geom_smooth(aes(color = TGP), method = "lm", se = TRUE, size = 1.2, alpha = 0.15) +
  labs(title = "offspring treatment = control",  
       x = "spring VPD CV (%)", 
       y = "number flowered", 
       color = "treatment group", 
       shape = "treatment group") +  
  scale_color_manual(values = c("CC" = "#386C9E", "DC" = "#177D31")) +
  scale_shape_manual(values = c("CC" = 16, "DC" = 17)) +
  ylim(-7, 35) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(color = "black", size = 10, face = "bold"),
        axis.title.y = element_text(color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        legend.title = element_text(color = "black", size = 10, face = "bold"),
        legend.text = element_text(color = "black", size = 10),
        plot.title = element_text(color = "black", size = 10, face = "bold", hjust = 0.5),
        legend.position = "none")
flow_control

ggsave("flow_drought.svg", flow_drought, path = "code-and-data/scripts/plpa-manuscript/figures",
       dpi = 300, width = 3, height = 3)
ggsave("flow_control.svg", flow_control, path = "code-and-data/scripts/plpa-manuscript/figures",
       dpi = 300, width = 3, height = 3)



# FIG 3 GRAPHS WITHOUT POPULATION LEVEL DATA----------------------------------
# root rxn norm, no pop --------------------------------------------------
root_na <- BIOMASS_FINAL |>
  filter(!is.na(ot) & !is.na(pt) & !is.na(root))

root_summary <- root_na |> 
  group_by(ot, pt) |> 
  summarise(mean_root = mean(root), 
            se_root = sd(root) / sqrt(n()), .groups = 'drop')

root_rxn <- ggplot(root_summary, aes(x = ot, y = mean_root, color = pt, shape = pt, group = pt)) +
  geom_point(aes(shape = pt), size = 2) +
  geom_errorbar(aes(ymin = mean_root - se_root, ymax = mean_root + se_root), width = 0.1, linewidth = 0.5) +
  geom_line(linewidth = 0.5) +
  labs(y = "root biomass (g)", x = " ", color = " ") +
  scale_color_manual(values = c("control" = "#548E90", "drought" = "#7A2218")) +
  scale_shape_manual(values = c("control" = 16, "drought" = 17)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(margin = margin(t = 15), color = "black", size = 14, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 15), color = "black", size = 14, face = "bold"),
        axis.text.x = element_text(color = "black", size = 11),  
        axis.text.y = element_text(color = "black", size = 11),  
        legend.title = element_text(color = "black", size = 14, face = "bold"),  
        legend.text = element_text(color = "black", size = 11), 
        plot.title = element_text(color = "black", size = 14),  
        plot.subtitle = element_text(color = "black", size = 10),
        plot.caption = element_text(color = "black", size = 10),
        legend.position = "none")

root_rxn

ggsave("root_rxn.svg", root_rxn, path = "code-and-data/scripts/plpa-manuscript/figures", dpi = 300, width = 2.83, height = 2.75)


# shoot rxn norm, no pop --------------------------------------------------
shoot_na <- BIOMASS_FINAL |>
  filter(!is.na(ot) & !is.na(pt) & !is.na(shoot))

shoot_summary <- shoot_na |> 
  group_by(ot, pt) |> 
  summarise(mean_shoot = mean(shoot), 
            se_shoot = sd(shoot) / sqrt(n()), .groups = 'drop')

shoot_rxn <- ggplot(shoot_summary, aes(x = ot, y = mean_shoot, color = pt, shape = pt, group = pt)) +
  geom_point(aes(shape = pt), size = 2) +
  geom_errorbar(aes(ymin = mean_shoot - se_shoot, ymax = mean_shoot + se_shoot), width = 0.1, linewidth = 0.5) +
  geom_line(linewidth = 0.5) +
  labs(y = "shoot biomass (g)", x = " ", color = " ") +
  scale_color_manual(values = c("control" = "#548E90", "drought" = "#7A2218")) +
  scale_shape_manual(values = c("control" = 16, "drought" = 17)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(margin = margin(t = 15), color = "black", size = 14, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 15), color = "black", size = 14, face = "bold"),
        axis.text.x = element_text(color = "black", size = 11),  
        axis.text.y = element_text(color = "black", size = 11),  
        legend.title = element_text(color = "black", size = 14, face = "bold"),  
        legend.text = element_text(color = "black", size = 11), 
        plot.title = element_text(color = "black", size = 14),  
        plot.subtitle = element_text(color = "black", size = 10),
        plot.caption = element_text(color = "black", size = 10),
        legend.position = "none")

shoot_rxn

ggsave("shoot_rxn.svg", shoot_rxn, path = "code-and-data/scripts/plpa-manuscript/figures", 
       dpi = 300, width = 2.83, height = 2.75)

# total rxn norm, no pop --------------------------------------------------
total_na <- BIOMASS_FINAL |>
  filter(!is.na(ot) & !is.na(pt) & !is.na(total_biomass))

total_summary <- total_na |> 
  group_by(ot, pt) |> 
  summarise(mean_total = mean(total_biomass), 
            se_total = sd(total_biomass) / sqrt(n()), .groups = 'drop')

total_rxn <- ggplot(total_summary, aes(x = ot, y = mean_total, color = pt, shape = pt, group = pt)) +
  geom_point(aes(shape = pt), size = 2) +
  geom_errorbar(aes(ymin = mean_total - se_total, ymax = mean_total + se_total), width = 0.1, linewidth = 0.5) +
  geom_line(linewidth = 0.5) +
  labs(y = "total biomass (g)", x = " ", color = " ") +
  scale_color_manual(values = c("control" = "#548E90", "drought" = "#7A2218")) +
  scale_shape_manual(values = c("control" = 16, "drought" = 17)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(margin = margin(t = 15), color = "black", size = 14, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 15), color = "black", size = 14, face = "bold"),
        axis.text.x = element_text(color = "black", size = 11),  
        axis.text.y = element_text(color = "black", size = 11),  
        legend.title = element_text(color = "black", size = 14, face = "bold"),  
        legend.text = element_text(color = "black", size = 11),  
        plot.title = element_text(color = "black", size = 14),  
        plot.subtitle = element_text(color = "black", size = 10),
        plot.caption = element_text(color = "black", size = 10),
        legend.position = "none")

total_rxn

ggsave("total_rxn.svg", total_rxn, path = "code-and-data/scripts/plpa-manuscript/figures",
       dpi = 300, width = 2.83, height = 2.75)

# r:s rxn norm, no pop ----------------------------------------------------
rs_na <- BIOMASS_FINAL |>
  filter(!is.na(ot) & !is.na(pt) & !is.na(ratio_RS))

rs_summary <- rs_na |> 
  group_by(ot, pt) |> 
  summarise(mean_rs = mean(ratio_RS), 
            se_rs = sd(ratio_RS) / sqrt(n()), .groups = 'drop')

rs_rxn <- ggplot(rs_summary, aes(x = ot, y = mean_rs, color = pt, shape = pt, group = pt)) +
  geom_point(aes(shape = pt), size = 2) +
  geom_errorbar(aes(ymin = mean_rs - se_rs, ymax = mean_rs + se_rs), width = 0.1, linewidth = 0.5) +
  geom_line(linewidth = 0.5) +
  labs(y = "R:S ratio", x = " ", color = " ") +
  scale_color_manual(values = c("control" = "#548E90", "drought" = "#7A2218")) +
  scale_shape_manual(values = c("control" = 16, "drought" = 17)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(margin = margin(t = 15), color = "black", size = 14, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 15), color = "black", size = 14, face = "bold"),
        axis.text.x = element_text(color = "black", size = 11),  
        axis.text.y = element_text(color = "black", size = 11),  
        legend.title = element_text(color = "black", size = 14, face = "bold"),  
        legend.text = element_text(color = "black", size = 11), 
        plot.title = element_text(color = "black", size = 14),  
        plot.subtitle = element_text(color = "black", size = 10),
        plot.caption = element_text(color = "black", size = 10),
        legend.position = "none")

rs_rxn

ggsave("rs_rxn.svg", rs_rxn, path = "code-and-data/scripts/plpa-manuscript/figures",
       dpi = 300, width = 2.83, height = 2.75)


# rgr rxn norm, no pop ----------------------------------------------------
rgr_na <- RGR_FINAL |>
  filter(!is.na(ot) & !is.na(pt) & !is.na(RGR))

rgr_summary <- rgr_na |> 
  group_by(ot, pt) |> 
  summarise(mean_rgr = mean(RGR), 
            se_rgr = sd(RGR) / sqrt(n()), .groups = 'drop')

rgr_rxn <- ggplot(rgr_summary, aes(x = ot, y = mean_rgr, color = pt, shape = pt, group = pt)) +
  geom_point(aes(shape = pt), size = 2) +
  geom_errorbar(aes(ymin = mean_rgr - se_rgr, ymax = mean_rgr + se_rgr), width = 0.1, linewidth = 0.5) +
  geom_line(linewidth = 0.5) +
  labs(y = "RGR", x = " ", color = " ") +
  scale_color_manual(values = c("control" = "#548E90", "drought" = "#7A2218")) +
  scale_shape_manual(values = c("control" = 16, "drought" = 17)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(margin = margin(t = 15), color = "black", size = 14, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 15), color = "black", size = 14, face = "bold"),
        axis.text.x = element_text(color = "black", size = 11),  
        axis.text.y = element_text(color = "black", size = 11),  
        legend.title = element_text(color = "black", size = 14, face = "bold"),  
        legend.text = element_text(color = "black", size = 11), 
        plot.title = element_text(color = "black", size = 14),  
        plot.subtitle = element_text(color = "black", size = 10),
        plot.caption = element_text(color = "black", size = 10),
        legend.position = "none")

rgr_rxn

ggsave("rgr_rxn.svg", rgr_rxn, path = "code-and-data/scripts/plpa-manuscript/figures",
       dpi = 300, width = 2.83, height = 2.75)

# ldmc rxn norm, no pop ---------------------------------------------------
ldmc_na <- SLA_LDMC_FINAL |>
  filter(!is.na(ot) & !is.na(pt) & !is.na(ldmc))

ldmc_summary <- ldmc_na |> 
  group_by(ot, pt) |> 
  summarise(mean_ldmc = mean(ldmc), 
            se_ldmc = sd(ldmc) / sqrt(n()), .groups = 'drop')

ldmc_rxn <- ggplot(ldmc_summary, aes(x = ot, y = mean_ldmc, color = pt, shape = pt, group = pt)) +
  geom_point(aes(shape = pt), size = 2) +
  geom_errorbar(aes(ymin = mean_ldmc - se_ldmc, ymax = mean_ldmc + se_ldmc), width = 0.1, linewidth = 0.5) +
  geom_line(linewidth = 0.5) +
  labs(y = "LDMC", x = " ", color = " ") +
  scale_color_manual(values = c("control" = "#548E90", "drought" = "#7A2218")) +
  scale_shape_manual(values = c("control" = 16, "drought" = 17)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(margin = margin(t = 15), color = "black", size = 14, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 15), color = "black", size = 14, face = "bold"),
        axis.text.x = element_text(color = "black", size = 11),  
        axis.text.y = element_text(color = "black", size = 11),  
        legend.title = element_text(color = "black", size = 14, face = "bold"),  
        legend.text = element_text(color = "black", size = 11),  
        plot.title = element_text(color = "black", size = 14),  
        plot.subtitle = element_text(color = "black", size = 10),
        plot.caption = element_text(color = "black", size = 10),
        legend.position = "none")

ldmc_rxn

ggsave("ldmc_rxn.svg", ldmc_rxn, path = "code-and-data/scripts/plpa-manuscript/figures",
       dpi = 300, width = 2.83, height = 2.75)

# seed number rxn norm, no pop --------------------------------------------
seed_na <- SEED_FINAL |>
  filter(!is.na(ot) & !is.na(pt) & !is.na(num_total))

seed_summary <- seed_na |> 
  group_by(ot, pt) |> 
  summarise(mean_seed = mean(num_total), 
            se_seed = sd(num_total) / sqrt(n()), .groups = 'drop')

seednum_rxn <- ggplot(seed_summary, aes(x = ot, y = mean_seed, color = pt, shape = pt, group = pt)) +
  geom_point(aes(shape = pt), size = 2) +
  geom_errorbar(aes(ymin = mean_seed - se_seed, ymax = mean_seed + se_seed), width = 0.1, linewidth = 0.5) +
  geom_line(linewidth = 0.5) +
  labs(y = "number of seeds", x = " ", color = " ") +
  scale_color_manual(values = c("control" = "#548E90", "drought" = "#7A2218")) +
  scale_shape_manual(values = c("control" = 16, "drought" = 17)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(margin = margin(t = 15), color = "black", size = 14, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 15), color = "black", size = 14, face = "bold"),
        axis.text.x = element_text(color = "black", size = 11),  
        axis.text.y = element_text(color = "black", size = 11),  
        legend.title = element_text(color = "black", size = 14, face = "bold"),  
        legend.text = element_text(color = "black", size = 11),  
        plot.title = element_text(color = "black", size = 14),  
        plot.subtitle = element_text(color = "black", size = 10),
        plot.caption = element_text(color = "black", size = 10),
        legend.position = "none")

seednum_rxn

ggsave("seednum_rxn.svg", seednum_rxn, path = "code-and-data/scripts/plpa-manuscript/figures",
       dpi = 300, width = 2.83, height = 2.75)

# seed num bar chart TGP
seed_summary_tgp <- seed_na |> 
  group_by(TGP) |>
  summarise(mean_seed = mean(num_total), 
            se_seed = sd(num_total) / sqrt(n()), .groups = 'drop')

pattern_types <- c("CC" = "none", "CD" = "stripe", "DC" = "circle", "DD" = "crosshatch")

seednum_bar <- ggplot(seed_summary_tgp, aes(x = TGP, y = mean_seed, fill = TGP)) +
  geom_col(width = 0.7, color = "black") +
  geom_errorbar(aes(ymin = mean_seed - se_seed, ymax = mean_seed + se_seed), 
                width = 0.2, linewidth = 0.5) + 
  scale_fill_manual(values = c("CC" = "grey80", "CD" = "grey60", "DC" = "grey40", "DD" = "grey20")) +  # Greyscale fills
  theme(legend.key.size = unit(1.5, 'cm')) + 
  labs(y = "number of seeds", x = "treatment", fill = "Treatment") +  # Removed pattern legend
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(margin = margin(t = 15), color = "black", size = 14, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 15), color = "black", size = 14, face = "bold"),
        axis.text.x = element_text(color = "black", size = 11),  
        axis.text.y = element_text(color = "black", size = 11),  
        legend.title = element_text(color = "black", size = 14, face = "bold"),  
        legend.text = element_text(color = "black", size = 11),  
        plot.title = element_text(color = "black", size = 14),  
        plot.subtitle = element_text(color = "black", size = 10),
        plot.caption = element_text(color = "black", size = 10),
        legend.position = "none")

seednum_bar

ggsave("seednum_bar.svg", seednum_bar, path = "code-and-data/scripts/plpa-manuscript/figures",
       dpi = 300, width = 3, height = 3)

# flowering status rxn norm, no pop ----------------------------------------------------
flowerstat_na <- FLOWER_STATUS |>
  filter(!is.na(ot) & !is.na(pt) & !is.na(status))

flowerstat_summary <- flowerstat_na |> 
  group_by(ot, pt) |> 
  summarise(mean_flow = mean(status * n()), 
            se_flow = sd(status * n()) / sqrt(n()), .groups = 'drop')

flowerstat_rxn <- ggplot(flowerstat_summary, aes(x = ot, y = mean_flow, color = pt, shape = pt, group = pt)) +
  geom_point(aes(shape = pt), size = 2) +
  geom_errorbar(aes(ymin = mean_flow - se_flow, ymax = mean_flow + se_flow), width = 0.1, linewidth = 0.5) +
  geom_line(linewidth = 0.5) +
  labs(y = "number flowered", x = " ", color = " ") +
  scale_color_manual(values = c("control" = "#548E90", "drought" = "#7A2218")) +
  scale_shape_manual(values = c("control" = 16, "drought" = 17)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(margin = margin(t = 15), color = "black", size = 14, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 15), color = "black", size = 14, face = "bold"),
        axis.text.x = element_text(color = "black", size = 11),  
        axis.text.y = element_text(color = "black", size = 11),  
        legend.title = element_text(color = "black", size = 14, face = "bold"),  
        legend.text = element_text(color = "black", size = 11),  
        plot.title = element_text(color = "black", size = 14),  
        plot.subtitle = element_text(color = "black", size = 10),
        plot.caption = element_text(color = "black", size = 10),
        legend.position = "none")
flowerstat_rxn

ggsave("flowerstat_rxn.svg", flowerstat_rxn, path = "code-and-data/scripts/plpa-manuscript/figures",
       dpi = 300, width = 2.83, height = 2.75)

# num of flowering structures rxn norm, no pop ----------------------------------------------------
flowernum_na <- FLOWER_FINAL |>
  filter(!is.na(ot) & !is.na(pt) & !is.na(num_structure))

flowernum_summary <- flowernum_na |> 
  group_by(ot, pt) |> 
  summarise(mean_flow = mean(num_structure),  # Mean number of flowering plants
            se_flow = sd(num_structure) / sqrt(n()),  # Standard error
            .groups = 'drop')

flowernum_rxn <- ggplot(flowernum_summary, aes(x = ot, y = mean_flow, color = pt, shape = pt, group = pt)) +
  geom_point(aes(shape = pt), size = 2) +
  geom_errorbar(aes(ymin = mean_flow - se_flow, ymax = mean_flow + se_flow), width = 0.1, linewidth = 0.5) +
  geom_line(linewidth = 0.5) +
  labs(y = "number of flower structures", x = " ", color = " ") +
  scale_color_manual(values = c("control" = "#548E90", "drought" = "#7A2218")) +
  scale_shape_manual(values = c("control" = 16, "drought" = 17)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(margin = margin(t = 15), color = "black", size = 14, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 15), color = "black", size = 14, face = "bold"),
        axis.text.x = element_text(color = "black", size = 11),  
        axis.text.y = element_text(color = "black", size = 11),  
        legend.title = element_text(color = "black", size = 14, face = "bold"),  
        legend.text = element_text(color = "black", size = 11),  
        plot.title = element_text(color = "black", size = 14),  
        plot.subtitle = element_text(color = "black", size = 10),
        plot.caption = element_text(color = "black", size = 10),
        legend.position = "none")
flowernum_rxn

ggsave("flowernum_rxn.svg", flowernum_rxn, path = "code-and-data/scripts/plpa-manuscript/figures",
       dpi = 300, width = 2.83, height = 2.75)


# seed mass rxn norm, no pop --------------------------------------------
seedm_na <- SEED_FINAL |>
  filter(!is.na(ot) & !is.na(pt) & !is.na(mass_total))

seedm_summary <- seedm_na |> 
  group_by(ot, pt) |> 
  summarise(mean_seed = mean(mass_total), 
            se_seed = sd(mass_total) / sqrt(n()), .groups = 'drop')

seedmass_rxn <- ggplot(seedm_summary, aes(x = ot, y = mean_seed, color = pt, shape = pt, group = pt)) +
  geom_point(aes(shape = pt), size = 2) +
  geom_errorbar(aes(ymin = mean_seed - se_seed, ymax = mean_seed + se_seed), width = 0.1, linewidth = 0.5) +
  geom_line(linewidth = 0.5) +
  labs(y = "seed mass (g)", x = " ", color = " ") +
  scale_color_manual(values = c("control" = "#548E90", "drought" = "#7A2218")) +
  scale_shape_manual(values = c("control" = 16, "drought" = 17)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(margin = margin(t = 15), color = "black", size = 14, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 15), color = "black", size = 14, face = "bold"),
        axis.text.x = element_text(color = "black", size = 11),  
        axis.text.y = element_text(color = "black", size = 11),  
        legend.title = element_text(color = "black", size = 14, face = "bold"),  
        legend.text = element_text(color = "black", size = 11),  
        plot.title = element_text(color = "black", size = 14),  
        plot.subtitle = element_text(color = "black", size = 10),
        plot.caption = element_text(color = "black", size = 10),
        legend.position = "none")

seedmass_rxn

ggsave("seedmass_rxn.svg", seedmass_rxn, path = "code-and-data/scripts/plpa-manuscript/figures",
       dpi = 300, width = 2.83, height = 2.75)

# seed mass bar chart TGP
seedm_summary_tgp <- seedm_na |> 
  group_by(TGP) |>
  summarise(mean_seed = mean(mass_total), 
            se_seed = sd(mass_total) / sqrt(n()), .groups = 'drop')

pattern_types <- c("CC" = "none", "CD" = "stripe", "DC" = "circle", "DD" = "crosshatch")

seedmass_bar <- ggplot(seedm_summary_tgp, aes(x = TGP, y = mean_seed, fill = TGP)) +
  geom_col(width = 0.7, color = "black") +
  geom_errorbar(aes(ymin = mean_seed - se_seed, ymax = mean_seed + se_seed), 
                width = 0.2, linewidth = 0.5) + 
  scale_fill_manual(values = c("CC" = "grey80", "CD" = "grey60", "DC" = "grey40", "DD" = "grey20")) +  # Greyscale fills
  theme(legend.key.size = unit(1.5, 'cm')) + 
  labs(y = "seed mass (g)", x = "treatment", fill = "Treatment") +  # Removed pattern legend
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(margin = margin(t = 15), color = "black", size = 14, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 15), color = "black", size = 14, face = "bold"),
        axis.text.x = element_text(color = "black", size = 11),  
        axis.text.y = element_text(color = "black", size = 11),  
        legend.title = element_text(color = "black", size = 14, face = "bold"),  
        legend.text = element_text(color = "black", size = 11),  
        plot.title = element_text(color = "black", size = 14),  
        plot.subtitle = element_text(color = "black", size = 10),
        plot.caption = element_text(color = "black", size = 10),
        legend.position = "none")

seedmass_bar

ggsave("seedmass_bar.svg", seedmass_bar, path = "code-and-data/scripts/plpa-manuscript/figures",
       dpi = 300, width = 3, height = 3)
