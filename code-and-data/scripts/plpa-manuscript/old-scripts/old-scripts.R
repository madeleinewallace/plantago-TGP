# OLD 05_visualization ---------------------------------------------------------

# VPD GRAPHS WITH POPULATION LEVEL DATA, fig 2 ---------------------------------
# color palette for temp, rainfall, seasonality
cv_colors <- scale_color_gradient2(low = "blue3", mid = "purple3", high = "red3", midpoint = 32.5)
vpd_colors <- scale_color_gradient2(low = "blue3", mid ="#FFA425", high = "#DC1C72", midpoint = 26)

# VPD - root rxn norm ------------------------------------------------
root_na <- BIOMASS_FINAL |>
  filter(!is.na(ot) & !is.na(pt) & !is.na(root))

root_droughtpt_summary <- root_na |> 
  group_by(ot, pt, spring_vpd_cv) |> 
  filter(pt == 'drought') |> 
  summarise(mean_root = mean(root), 
            se_root = sd(root) / sqrt(n()), .groups = 'drop')

root_controlpt_summary <- root_na |> 
  group_by(ot, pt, spring_vpd_cv) |> 
  filter(pt == 'control') |> 
  summarise(mean_root = mean(root), 
            se_root = sd(root) / sqrt(n()), .groups = 'drop')

# graph 1
vpdroot_droughtpt_rxn <- ggplot(root_droughtpt_summary, aes(x = ot, y = mean_root, color = spring_vpd_cv)) +
  geom_point(size = 1, alpha = 0.9) +
  #geom_errorbar(aes(ymin = mean_root - se_root, ymax = mean_root + se_root), width = 0.07, alpha = 0.9) +
  geom_line(alpha = 0.9) +
  labs(y = "root biomass (g)", x = "offspring treatment", color = "VPD (Hpa)", title = "drought parental treatment") +
  vpd_colors +
  coord_fixed(ratio = 4) +
  scale_size_continuous(range = c(0.5, 1.5)) +
  scale_x_discrete(expand = c(0.1, 0.1)) +
  ylim(0, 0.85) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(margin = margin(t = 4), color = "black", size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 4), color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10),  
        axis.text.y = element_text(color = "black", size = 10),  
        legend.title = element_text(color = "black", size = 10, face = "bold"),  
        legend.text = element_text(color = "black", size = 10),  
        plot.title = element_text(color = "black", size = 10, face = "bold", hjust = 0.5),  
        plot.subtitle = element_text(color = "black", size = 10),
        plot.caption = element_text(color = "black", size = 10))
vpdroot_droughtpt_rxn

root_controlpt_rxn <- ggplot(root_controlpt_summary, aes(x = ot, y = mean_root, color = SAP_mm, group = pop)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = mean_root - se_root, ymax = mean_root + se_root), width = 0.1) +
  geom_line(size = 1.5) +
  labs(y = "root biomass (g)", x = "offspring treatment", color = "SAP (mm)", title = "control parental treatment") +
  rainfall_colors +
  ylim(0, 0.85) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(margin = margin(t = 4), color = "black", size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 4), color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10),  
        axis.text.y = element_text(color = "black", size = 10),  
        legend.title = element_text(color = "black", size = 10, face = "bold"),  
        legend.text = element_text(color = "black", size = 10),  
        plot.title = element_text(color = "black", size = 10, face = "bold", hjust = 0.5),  
        plot.subtitle = element_text(color = "black", size = 10),
        plot.caption = element_text(color = "black", size = 10))
root_controlpt_rxn

rootbiomass_rxn_pop_fig <- root_controlpt_rxn + root_droughtpt_rxn +
  plot_layout(guides = "collect")
rootbiomass_rxn_pop_fig

ggsave("rootbiomass_rxn_pop_fig2.png", rootbiomass_rxn_pop_fig, path = "code-and-data/scripts/plpa-manuscript/figures",
       height = 4, width = 8, dpi = 300)


# graph 2
root_summary <- root_na |> 
  group_by(TGP, vpd, cv) |> 
  summarise(mean_root = mean(root), 
            se_root = sd(root) / sqrt(n()), .groups = 'drop')

rootgraph2 <- ggplot(root_summary, aes(x = vpd, y = mean_root, linetype = TGP, color = TGP, shape = TGP)) +
  geom_point(aes(color = TGP, shape = TGP), size = 4) +
  geom_smooth(method = "lm", se = TRUE, size = 1.2, color = 'black', alpha = 0.15) +
  labs(x = "VPD (Hpa)", y = "Root Biomass (g)", color = "Treatment Group", 
       linetype = "Treatment Group", shape = "Treatment Group") + 
  scale_linetype_manual(values = c("CC" = "solid", 
                                   "CD" = "dotdash", 
                                   "DC" = "dashed", 
                                   "DD" = "dotted")) +
  scale_shape_manual(values = c("CC" = 15,  
                                "CD" = 16,
                                "DC" = 17,
                                "DD" = 18)) + 
  scale_color_manual(values = c("CC" = "#0068A7",
                                "CD" = "#C66D9B",
                                "DC" = "#E79321",
                                "DD" = "#D25014")) +
  guides(linetype = guide_legend(override.aes = list(size = 5), keywidth = 4)) +
  #annotate('text', x = 30, y = 0.6, label = 'INSERT', size = 4) +
  #annotate('text', x = 30, y = 0.57, label = 'INSERT', size = 4) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(color = "black", size = 10, face = "bold"),
        axis.title.y = element_text(color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        legend.title = element_text(color = "black", size = 10, face = "bold"),
        legend.text = element_text(color = "black", size = 10))
rootgraph2

# VPD - shoot rxn norm ------------------------------------------------
# graph 2
shoot_na <- BIOMASS_FINAL |>
  filter(!is.na(ot) & !is.na(pt) & !is.na(shoot))

shoot_summary <- shoot_na |> 
  group_by(TGP, vpd, cv) |> 
  summarise(mean_shoot = mean(shoot), 
            se_shoot = sd(shoot) / sqrt(n()), .groups = 'drop')

shootgraph2 <- ggplot(shoot_summary, aes(x = vpd, y = mean_shoot, color = cv, linetype = TGP)) +
  geom_point(aes(color = TGP, shape = TGP), size = 4) +
  geom_smooth(method = "lm", se = TRUE, size = 1.2, color = 'black', alpha = 0.15) +
  labs(x = "VPD (Hpa)", y = "Shoot Biomass (g)", color = "Treatment Group", 
       linetype = "Treatment Group", shape = "Treatment Group") + 
  scale_linetype_manual(values = c("CC" = "solid", 
                                   "CD" = "dotdash", 
                                   "DC" = "dashed", 
                                   "DD" = "dotted")) +
  scale_shape_manual(values = c("CC" = 15,  
                                "CD" = 16,
                                "DC" = 17,
                                "DD" = 18)) + 
  scale_color_manual(values = c("CC" = "#0068A7",
                                "CD" = "#C66D9B",
                                "DC" = "#E79321",
                                "DD" = "#D25014")) +
  guides(linetype = guide_legend(override.aes = list(size = 5), keywidth = 4)) +
  #annotate('text', x = 30, y = 0.6, label = 'INSERT', size = 4) +
  #annotate('text', x = 30, y = 0.57, label = 'INSERT', size = 4) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(color = "black", size = 10, face = "bold"),
        axis.title.y = element_text(color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        legend.title = element_text(color = "black", size = 10, face = "bold"),
        legend.text = element_text(color = "black", size = 10))
shootgraph2

# VPD - RS rxn norm ------------------------------------------------
# graph 2
rs_na <- BIOMASS_FINAL |>
  filter(!is.na(ot) & !is.na(pt) & !is.na(ratio_RS))

rs_summary <- rs_na |> 
  group_by(TGP, vpd, cv) |> 
  summarise(mean_rs = mean(ratio_RS), 
            se_rs = sd(ratio_RS) / sqrt(n()), .groups = 'drop')

rsgraph2 <- ggplot(rs_summary, aes(x = vpd, y = mean_rs, color = cv, linetype = TGP)) +
  geom_point(aes(color = TGP, shape = TGP), size = 4) +
  geom_smooth(method = "lm", se = TRUE, size = 1.2, color = 'black', alpha = 0.15) +
  labs(x = "VPD (Hpa)", y = "R:S Ratio", color = "Treatment Group", 
       linetype = "Treatment Group", shape = "Treatment Group") + 
  scale_linetype_manual(values = c("CC" = "solid", 
                                   "CD" = "dotdash", 
                                   "DC" = "dashed", 
                                   "DD" = "dotted")) +
  scale_shape_manual(values = c("CC" = 15,  
                                "CD" = 16,
                                "DC" = 17,
                                "DD" = 18)) + 
  scale_color_manual(values = c("CC" = "#0068A7",
                                "CD" = "#C66D9B",
                                "DC" = "#E79321",
                                "DD" = "#D25014")) +
  guides(linetype = guide_legend(override.aes = list(size = 5), keywidth = 4)) +
  #annotate('text', x = 30, y = 0.6, label = 'INSERT', size = 4) +
  #annotate('text', x = 30, y = 0.57, label = 'INSERT', size = 4) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(color = "black", size = 10, face = "bold"),
        axis.title.y = element_text(color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        legend.title = element_text(color = "black", size = 10, face = "bold"),
        legend.text = element_text(color = "black", size = 10))
rsgraph2

# VPD - height rxn norm ------------------------------------------------
# graph 2
max_na <- HEIGHT_FINAL |>
  filter(!is.na(ot) & !is.na(pt) & !is.na(max))

max_summary <- max_na |> 
  group_by(TGP, vpd, cv) |> 
  summarise(mean_max = mean(max), 
            se_max = sd(max) / sqrt(n()), .groups = 'drop')

maxgraph2 <- ggplot(max_summary, aes(x = vpd, y = mean_max, color = cv, linetype = TGP)) +
  geom_point(aes(color = TGP, shape = TGP), size = 4) +
  geom_smooth(method = "lm", se = TRUE, size = 1.2, color = 'black', alpha = 0.15) +
  labs(x = "VPD (Hpa)", y = "Max Height (cm)", color = "Treatment Group", 
       linetype = "Treatment Group", shape = "Treatment Group") + 
  scale_linetype_manual(values = c("CC" = "solid", 
                                   "CD" = "dotdash", 
                                   "DC" = "dashed", 
                                   "DD" = "dotted")) +
  scale_shape_manual(values = c("CC" = 15,  
                                "CD" = 16,
                                "DC" = 17,
                                "DD" = 18)) + 
  scale_color_manual(values = c("CC" = "#0068A7",
                                "CD" = "#C66D9B",
                                "DC" = "#E79321",
                                "DD" = "#D25014")) +
  guides(linetype = guide_legend(override.aes = list(size = 5), keywidth = 4)) +
  #annotate('text', x = 30, y = 0.6, label = 'INSERT', size = 4) +
  #annotate('text', x = 30, y = 0.57, label = 'INSERT', size = 4) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(color = "black", size = 10, face = "bold"),
        axis.title.y = element_text(color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        legend.title = element_text(color = "black", size = 10, face = "bold"),
        legend.text = element_text(color = "black", size = 10))
maxgraph2

# VPD - RGR rxn norm ------------------------------------------------
# graph 2
rgr_na <- RGR_FINAL |>
  filter(!is.na(ot) & !is.na(pt) & !is.na(RGR))

rgr_summary <- rgr_na |> 
  group_by(TGP, vpd, cv) |> 
  summarise(mean_rgr = mean(RGR), 
            se_rgr = sd(RGR) / sqrt(n()), .groups = 'drop')

rgrgraph2 <- ggplot(rgr_summary, aes(x = vpd, y = mean_rgr, color = cv, linetype = TGP)) +
  geom_point(aes(color = TGP, shape = TGP), size = 4) +
  geom_smooth(method = "lm", se = TRUE, size = 1.2, color = 'black', alpha = 0.15) +
  labs(x = "VPD (Hpa)", y = "RGR", color = "Treatment Group", 
       linetype = "Treatment Group", shape = "Treatment Group") + 
  scale_linetype_manual(values = c("CC" = "solid", 
                                   "CD" = "dotdash", 
                                   "DC" = "dashed", 
                                   "DD" = "dotted")) +
  scale_shape_manual(values = c("CC" = 15,  
                                "CD" = 16,
                                "DC" = 17,
                                "DD" = 18)) + 
  scale_color_manual(values = c("CC" = "#0068A7",
                                "CD" = "#C66D9B",
                                "DC" = "#E79321",
                                "DD" = "#D25014")) +
  guides(linetype = guide_legend(override.aes = list(size = 5), keywidth = 4)) +
  #annotate('text', x = 30, y = 0.6, label = 'INSERT', size = 4) +
  #annotate('text', x = 30, y = 0.57, label = 'INSERT', size = 4) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(color = "black", size = 10, face = "bold"),
        axis.title.y = element_text(color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        legend.title = element_text(color = "black", size = 10, face = "bold"),
        legend.text = element_text(color = "black", size = 10))
rgrgraph2

# VPD - mort rxn norm ------------------------------------------------
# graph 2
mort_na <- mort_day50 |>
  filter(!is.na(ot) & !is.na(pt) & !is.na(status))

mort_summary <- mort_na |> 
  group_by(TGP.x, vpd, cv) |> 
  summarise(mean_mort = mean(status), 
            se_mort = sd(status) / sqrt(n()), .groups = 'drop')

mortgraph2 <- ggplot(mort_summary, aes(x = vpd, y = mean_mort, color = cv, linetype = TGP.x)) +
  geom_point(aes(color = TGP.x, shape = TGP.x), size = 4) +
  geom_smooth(method = "glm", se = TRUE, size = 1.2, color = 'black', alpha = 0.15) +
  labs(x = "VPD (Hpa)", y = "Proportion Alive", color = "Treatment Group", 
       linetype = "Treatment Group", shape = "Treatment Group") + 
  scale_linetype_manual(values = c("CC" = "solid", 
                                   "CD" = "dotdash", 
                                   "DC" = "dashed", 
                                   "DD" = "dotted")) +
  scale_shape_manual(values = c("CC" = 15,  
                                "CD" = 16,
                                "DC" = 17,
                                "DD" = 18)) + 
  scale_color_manual(values = c("CC" = "#0068A7",
                                "CD" = "#C66D9B",
                                "DC" = "#E79321",
                                "DD" = "#D25014")) +
  guides(linetype = guide_legend(override.aes = list(size = 5), keywidth = 4)) +
  #annotate('text', x = 30, y = 0.6, label = 'INSERT', size = 4) +
  #annotate('text', x = 30, y = 0.57, label = 'INSERT', size = 4) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(color = "black", size = 10, face = "bold"),
        axis.title.y = element_text(color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        legend.title = element_text(color = "black", size = 10, face = "bold"),
        legend.text = element_text(color = "black", size = 10))
mortgraph2

# VPD - seed num rxn norm ------------------------------------------------
# graph 2
seednum_na <- SEED_FINAL |>
  filter(!is.na(ot) & !is.na(pt) & !is.na(num_total)) #|> 
#filter(!pop %in% c(3, 5, 11))

seednum_summary <- seednum_na |> 
  group_by(TGP, vpd, cv) |> 
  summarise(mean_num = mean(num_total), 
            se_num = sd(num_total) / sqrt(n()), .groups = 'drop')

seednum_summary <- seednum_summary |> 
  mutate(vpd.k = vpd / 10)

seednumgraph2 <- ggplot(seednum_summary, aes(x = vpd.k, y = mean_num, color = TGP)) +
  geom_point(size = 3, shape = 16) +
  geom_smooth(aes(color = TGP), method = "lm", se = TRUE, size = 1.2, alpha = 0.15) +
  labs(x = "VPD (kPa)", y = "seed number", color = "treatment group") + 
  scale_color_manual(values = c("CC" = "#386C9E",
                                "CD" = "#177D31",
                                "DC" = "#C7B245",
                                "DD" = "#A32B6B")) +
  guides(color = guide_legend(override.aes = list(size = 5), keywidth = 4)) +
  #annotate('text', x = 30, y = 0.6, label = 'INSERT', size = 4) +
  #annotate('text', x = 30, y = 0.57, label = 'INSERT', size = 4) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(color = "black", size = 10, face = "bold"),
        axis.title.y = element_text(color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        legend.title = element_text(color = "black", size = 10, face = "bold"),
        legend.text = element_text(color = "black", size = 10))


seednumgraph2

ggsave("seednumgraph2.svg", seednumgraph2, path = "code-and-data/scripts/plpa-manuscript/figures")

# VPD - days to flower rxn norm ------------------------------------------------
# graph 2
daysflower_na <- FLOWER_FINAL |>
  filter(!is.na(ot) & !is.na(pt) & !is.na(days_to_flower))

daysflow_summary <- daysflower_na |> 
  group_by(TGP, spring_vpd_cv) |> 
  summarise(mean_days = mean(days_to_flower), 
            se_days = sd(days_to_flower) / sqrt(n()), .groups = 'drop')

daysflow_summary <- daysflow_summary |> 
  mutate(vpd.k = vpd / 10)

daysflowgraph2 <- ggplot(daysflow_summary, aes(x = spring_vpd_cv, y = mean_days, color = TGP)) +
  geom_point(size = 3, shape = 16) +
  geom_smooth(aes(color = TGP), method = "lm", se = TRUE, size = 1.2, alpha = 0.15) +
  labs(x = "VPD (kPa)", y = "days to flower", color = "treatment group") + 
  scale_color_manual(values = c("CC" = "#386C9E",
                                "CD" = "#177D31",
                                "DC" = "#C7B245",
                                "DD" = "#A32B6B")) +
  guides(color = guide_legend(override.aes = list(size = 5), keywidth = 4)) +
  #annotate('text', x = 30, y = 0.6, label = 'INSERT', size = 4) +
  #annotate('text', x = 30, y = 0.57, label = 'INSERT', size = 4) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(color = "black", size = 10, face = "bold"),
        axis.title.y = element_text(color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        legend.title = element_text(color = "black", size = 10, face = "bold"),
        legend.text = element_text(color = "black", size = 10))
daysflowgraph2

ggsave("daysflowgraph2.svg", daysflowgraph2, path = "code-and-data/scripts/plpa-manuscript/figures")


vpd_rxn_fig2_1 <- rootgraph2 + shootgraph2
vpd_rxn_fig2_2 <- maxgraph2 + rgrgraph2
vpd_rxn_fig2_3 <- mortgraph2 + seednumgraph2 + daysflowgraph2
vpd_rxn_fig2_1
vpd_rxn_fig2_2

ggsave("root_fig2.png", rootgraph2, path = "code-and-data/scripts/plpa-manuscript/figures",
       height = 6, width = 8, dpi = 300)
ggsave("shoot_fig2.png", shootgraph2, path = "code-and-data/scripts/plpa-manuscript/figures",
       height = 6, width = 8, dpi = 300)
ggsave("rs_fig2.png", rsgraph2, path = "code-and-data/scripts/plpa-manuscript/figures",
       height = 6, width = 8, dpi = 300)
ggsave("max_fig2.png", maxgraph2, path = "code-and-data/scripts/plpa-manuscript/figures",
       height = 6, width = 8, dpi = 300)
ggsave("rgr_fig2.png", rgrgraph2, path = "code-and-data/scripts/plpa-manuscript/figures",
       height = 6, width = 8, dpi = 300)
ggsave("mort_fig2.png", mortgraph2, path = "code-and-data/scripts/plpa-manuscript/figures",
       height = 6, width = 8, dpi = 300)
ggsave("seednum_fig2.png", seednumgraph2, path = "code-and-data/scripts/plpa-manuscript/figures",
       height = 6, width = 8, dpi = 300)
ggsave("daysflow_fig2.png", daysflowgraph2, path = "code-and-data/scripts/plpa-manuscript/figures",
       height = 6, width = 8, dpi = 300)


# OLD GRAPHS --------------------------------------------------------------
# max height rxn norm, with pops ------------------------------------------
max_na <- HEIGHT_FINAL |>
  filter(!is.na(ot) & !is.na(pt) & !is.na(max))

max_droughtpt_summary <- max_na |> 
  group_by(ot, pt, pop, SAP_mm, SAT_C, cv) |> 
  filter(pt == 'drought') |> 
  summarise(mean_max = mean(max), 
            se_max = sd(max) / sqrt(n()), .groups = 'drop')

max_controlpt_summary <- max_na |> 
  group_by(ot, pt, pop, SAP_mm, SAT_C, cv) |> 
  filter(pt == 'control') |> 
  summarise(mean_max = mean(max), 
            se_max = sd(max) / sqrt(n()), .groups = 'drop')

max_droughtpt_rxn <- ggplot(max_droughtpt_summary, aes(x = ot, y = mean_max, color = SAP_mm, group = pop)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = mean_max - se_max, ymax = mean_max + se_max), width = 0.1) +
  geom_line(size=1.5) +
  labs(y = "maximum height (cm)", x = "offspring treatment", color = "SAP (mm)", title = "drought parental treatment") +
  rainfall_colors +
  ylim(0, 9) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(margin = margin(t = 4), color = "black", size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 4), color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10),  
        axis.text.y = element_text(color = "black", size = 10),  
        legend.title = element_text(color = "black", size = 10, face = "bold"),  
        legend.text = element_text(color = "black", size = 10),  
        plot.title = element_text(color = "black", size = 10, face = "bold", hjust = 0.5),  
        plot.subtitle = element_text(color = "black", size = 10),
        plot.caption = element_text(color = "black", size = 10))
max_droughtpt_rxn

max_controlpt_rxn <- ggplot(max_controlpt_summary, aes(x = ot, y = mean_max, color = SAP_mm, group = pop)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = mean_max - se_max, ymax = mean_max + se_max), width = 0.1) +
  geom_line(size=1.5) +
  labs(y = "maximum height (cm)", x = "offspring treatment", color = "SAP (mm)", title = "control parental treatment") +
  rainfall_colors +
  ylim(0, 9) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(margin = margin(t = 4), color = "black", size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 4), color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10),  
        axis.text.y = element_text(color = "black", size = 10),  
        legend.title = element_text(color = "black", size = 10, face = "bold"),  
        legend.text = element_text(color = "black", size = 10),  
        plot.title = element_text(color = "black", size = 10, face = "bold", hjust = 0.5),  
        plot.subtitle = element_text(color = "black", size = 10),
        plot.caption = element_text(color = "black", size = 10))
max_controlpt_rxn

maxheight_rxn_pop_fig <- max_controlpt_rxn + max_droughtpt_rxn +
  plot_layout(guides = "collect")
maxheight_rxn_pop_fig

ggsave("maxheight_rxn_pop_fig2.png", maxheight_rxn_pop_fig, path = "code-and-data/scripts/plpa-manuscript/figures",
       height = 4, width = 8, dpi = 300)


# seed number rxn norm, with pops ------------------------------------------
seednumpop_na <- SEED_FINAL |>
  filter(!is.na(ot) & !is.na(pt) & !is.na(num_total))

seednum_droughtpt_summary <- seednumpop_na |> 
  group_by(ot, pt, pop, SAP_mm, SAT_C, cv) |> 
  filter(pt == 'drought') |> 
  summarise(mean_num = mean(num_total), 
            se_num = sd(num_total) / sqrt(n()), .groups = 'drop') |> 
  filter(!pop %in% c(3, 5, 11))

seednum_controlpt_summary <- seednumpop_na |> 
  group_by(ot, pt, pop, SAP_mm, SAT_C, cv) |> 
  filter(pt == 'control') |> 
  summarise(mean_num = mean(num_total), 
            se_num = sd(num_total) / sqrt(n()), .groups = 'drop') |> 
  filter(!pop %in% c(3, 5, 11))

seednum_droughtpt_rxn <- ggplot(seednum_droughtpt_summary, aes(x = ot, y = mean_num, color = SAP_mm, group = pop)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = mean_num - se_num, ymax = mean_num + se_num), width = 0.1) +
  geom_line(size=1.5) +
  labs(y = "number of seeds", x = "offspring treatment", color = "SAP (mm)", title = "drought parental treatment") +
  rainfall_colors +
  ylim(0, 100) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(margin = margin(t = 4), color = "black", size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 4), color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10),  
        axis.text.y = element_text(color = "black", size = 10),  
        legend.title = element_text(color = "black", size = 10, face = "bold"),  
        legend.text = element_text(color = "black", size = 10),  
        plot.title = element_text(color = "black", size = 10, face = "bold", hjust = 0.5),  
        plot.subtitle = element_text(color = "black", size = 10),
        plot.caption = element_text(color = "black", size = 10))
seednum_droughtpt_rxn

seednum_controlpt_rxn <- ggplot(seednum_controlpt_summary, aes(x = ot, y = mean_num, color = SAP_mm, group = pop)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = mean_num - se_num, ymax = mean_num + se_num), width = 0.1) +
  geom_line(size=1.5) +
  labs(y = "number of seeds", x = "offspring treatment", color = "SAP (mm)", title = "control parental treatment") +
  rainfall_colors +
  ylim(0, 100) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(margin = margin(t = 4), color = "black", size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 4), color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10),  
        axis.text.y = element_text(color = "black", size = 10),  
        legend.title = element_text(color = "black", size = 10, face = "bold"),  
        legend.text = element_text(color = "black", size = 10),  
        plot.title = element_text(color = "black", size = 10, face = "bold", hjust = 0.5),  
        plot.subtitle = element_text(color = "black", size = 10),
        plot.caption = element_text(color = "black", size = 10))
seednum_controlpt_rxn

seednum_rxn_pop_fig <- seednum_controlpt_rxn + seednum_droughtpt_rxn +
  plot_layout(guides = "collect")
seednum_rxn_pop_fig

ggsave("seednum_rxn_pop_fig2.png", seednum_rxn_pop_fig, path = "code-and-data/scripts/plpa-manuscript/figures",
       height = 4, width = 8, dpi = 300)



# mortality rxn norm, with pops ------------------------------------------
mort_na <- mort_day50 |>
  filter(!is.na(ot) & !is.na(pt) & !is.na(status))

mort_droughtpt_summary <- mort_na |> 
  group_by(ot, pt, pop, SAP_mm, SAT_C, cv) |> 
  filter(pt == 'drought') |> 
  summarise(mean_mort = mean(status), 
            se_mort = sd(status) / sqrt(n()), .groups = 'drop')

mort_controlpt_summary <- mort_na |> 
  group_by(ot, pt, pop, SAP_mm, SAT_C, cv) |> 
  filter(pt == 'control') |> 
  summarise(mean_mort = mean(status), 
            se_mort = sd(status) / sqrt(n()), .groups = 'drop')

mort_droughtpt_rxn <- ggplot(mort_droughtpt_summary, aes(x = ot, y = mean_mort, color = SAT_C, group = pop)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = mean_mort - se_mort, ymax = mean_mort + se_mort), width = 0.1) +
  geom_line(size=1.5) +
  labs(y = "proportion alive", x = "offspring treatment", color = "SAT (C)", title = "drought parental treatment") +
  temperature_colors +
  ylim(0, 1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(margin = margin(t = 4), color = "black", size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 4), color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10),  
        axis.text.y = element_text(color = "black", size = 10),  
        legend.title = element_text(color = "black", size = 10, face = "bold"),  
        legend.text = element_text(color = "black", size = 10),  
        plot.title = element_text(color = "black", size = 10, face = "bold", hjust = 0.5),  
        plot.subtitle = element_text(color = "black", size = 10),
        plot.caption = element_text(color = "black", size = 10))
mort_droughtpt_rxn

mort_controlpt_rxn <- ggplot(mort_controlpt_summary, aes(x = ot, y = mean_mort, color = SAT_C, group = pop)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = mean_mort - se_mort, ymax = mean_mort + se_mort), width = 0.1) +
  geom_line(size=1.5) +
  labs(y = "proportion alive", x = "offspring treatment", color = "SAT (C)", title = "control parental treatment") +
  temperature_colors +
  ylim(0, 1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(margin = margin(t = 4), color = "black", size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 4), color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10),  
        axis.text.y = element_text(color = "black", size = 10),  
        legend.title = element_text(color = "black", size = 10, face = "bold"),  
        legend.text = element_text(color = "black", size = 10),  
        plot.title = element_text(color = "black", size = 10, face = "bold", hjust = 0.5),  
        plot.subtitle = element_text(color = "black", size = 10),
        plot.caption = element_text(color = "black", size = 10))
mort_controlpt_rxn

mort_rxn_pop_fig <- mort_controlpt_rxn + mort_droughtpt_rxn +
  plot_layout(guides = "collect")
mort_rxn_pop_fig

ggsave("mort_rxn_pop_fig2.png", mort_rxn_pop_fig, path = "code-and-data/scripts/plpa-manuscript/figures",
       height = 4, width = 8, dpi = 300)


# total rxn norm, with pop ------------------------------------------------
total_na <- BIOMASS_FINAL |>
  filter(!is.na(ot) & !is.na(pt) & !is.na(total_biomass))

total_droughtpt_summary <- total_na |> 
  group_by(ot, pt, pop, SAP_mm, SAT_C, cv) |> 
  filter(pt == 'drought') |> 
  summarise(mean_total = mean(total_biomass), 
            se_total = sd(total_biomass) / sqrt(n()), .groups = 'drop')

total_controlpt_summary <- total_na |> 
  group_by(ot, pt, pop, SAP_mm, SAT_C, cv) |> 
  filter(pt == 'control') |> 
  summarise(mean_total = mean(total_biomass), 
            se_total = sd(total_biomass) / sqrt(n()), .groups = 'drop')

total_droughtpt_rxn <- ggplot(total_droughtpt_summary, aes(x = ot, y = mean_total, color = cv, group = pop)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = mean_total - se_total, ymax = mean_total + se_total), width = 0.1) +
  geom_line(size = 1.5) +
  labs(y = "total biomass (g)", x = "offspring treatment", color = "CV (%)", title = "drought parental treatment") +
  seasonality_colors +
  ylim(0, 1.1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(margin = margin(t = 4), color = "black", size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 4), color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10),  
        axis.text.y = element_text(color = "black", size = 10),  
        legend.title = element_text(color = "black", size = 10, face = "bold"),  
        legend.text = element_text(color = "black", size = 10),  
        plot.title = element_text(color = "black", size = 10, face = "bold", hjust = 0.5),  
        plot.subtitle = element_text(color = "black", size = 10),
        plot.caption = element_text(color = "black", size = 10))
total_droughtpt_rxn

total_controlpt_rxn <- ggplot(total_controlpt_summary, aes(x = ot, y = mean_total, color = cv, group = pop)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = mean_total - se_total, ymax = mean_total + se_total), width = 0.1) +
  geom_line(size = 1.5) +
  labs(y = "total biomass (g)", x = "offspring treatment", color = "CV (%)", title = "control parental treatment") +
  seasonality_colors +
  ylim(0, 1.1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(margin = margin(t = 4), color = "black", size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 4), color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10),  
        axis.text.y = element_text(color = "black", size = 10),  
        legend.title = element_text(color = "black", size = 10, face = "bold"),  
        legend.text = element_text(color = "black", size = 10),  
        plot.title = element_text(color = "black", size = 10, face = "bold", hjust = 0.5),  
        plot.subtitle = element_text(color = "black", size = 10),
        plot.caption = element_text(color = "black", size = 10))
total_controlpt_rxn

totalbiomass_rxn_pop_fig <- total_controlpt_rxn + total_droughtpt_rxn +
  plot_layout(guides = "collect")
totalbiomass_rxn_pop_fig

ggsave("totalbiomass_rxn_pop_fig2.png", totalbiomass_rxn_pop_fig, path = "code-and-data/scripts/plpa-manuscript/figures",
       height = 4, width = 8, dpi = 300)


# R:S ratio rxn norm, with pop ------------------------------------------------
rs_na <- BIOMASS_FINAL |>
  filter(!is.na(ot) & !is.na(pt) & !is.na(ratio_RS))

rs_droughtpt_summary <- rs_na |> 
  group_by(ot, pt, pop, SAP_mm, SAT_C, cv) |> 
  filter(pt == 'drought') |> 
  summarise(mean_rs = mean(ratio_RS), 
            se_rs = sd(ratio_RS) / sqrt(n()), .groups = 'drop')

rs_controlpt_summary <- rs_na |> 
  group_by(ot, pt, pop, SAP_mm, SAT_C, cv) |> 
  filter(pt == 'control') |> 
  summarise(mean_rs = mean(ratio_RS), 
            se_rs = sd(ratio_RS) / sqrt(n()), .groups = 'drop')

rs_droughtpt_rxn <- ggplot(rs_droughtpt_summary, aes(x = ot, y = mean_rs, color = cv, group = pop)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = mean_rs - se_rs, ymax = mean_rs + se_rs), width = 0.1) +
  geom_line(size = 1.5) +
  labs(y = "root:shoot", x = "offspring treatment", color = "CV (%)", title = "drought parental treatment") +
  seasonality_colors +
  ylim(0, 14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(margin = margin(t = 4), color = "black", size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 4), color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10),  
        axis.text.y = element_text(color = "black", size = 10),  
        legend.title = element_text(color = "black", size = 10, face = "bold"),  
        legend.text = element_text(color = "black", size = 10),  
        plot.title = element_text(color = "black", size = 10, face = "bold", hjust = 0.5),  
        plot.subtitle = element_text(color = "black", size = 10),
        plot.caption = element_text(color = "black", size = 10))
rs_droughtpt_rxn

rs_controlpt_rxn <- ggplot(rs_controlpt_summary, aes(x = ot, y = mean_rs, color = cv, group = pop)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = mean_rs - se_rs, ymax = mean_rs + se_rs), width = 0.1) +
  geom_line(size = 1.5) +
  labs(y = "root:shoot", x = "offspring treatment", color = "CV (%)", title = "control parental treatment") +
  seasonality_colors +
  ylim(0, 14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(margin = margin(t = 4), color = "black", size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 4), color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10),  
        axis.text.y = element_text(color = "black", size = 10),  
        legend.title = element_text(color = "black", size = 10, face = "bold"),  
        legend.text = element_text(color = "black", size = 10),  
        plot.title = element_text(color = "black", size = 10, face = "bold", hjust = 0.5),  
        plot.subtitle = element_text(color = "black", size = 10),
        plot.caption = element_text(color = "black", size = 10))
rs_controlpt_rxn


rs_rxn_pop_fig <- rs_controlpt_rxn + rs_droughtpt_rxn +
  plot_layout(guides = "collect")
rs_rxn_pop_fig

ggsave("rs_rxn_pop_fig2.png", rs_rxn_pop_fig, path = "code-and-data/scripts/plpa-manuscript/figures",
       height = 4, width = 8, dpi = 300)


# SLA rxn norm, with pop ------------------------------------------------
sla_na <- SLA_LDMC_FINAL |>
  filter(!is.na(ot) & !is.na(pt) & !is.na(sla))

sla_droughtpt_summary <- sla_na |> 
  group_by(ot, pt, pop, SAP_mm, SAT_C, cv) |> 
  filter(pt == 'drought') |> 
  summarise(mean_sla = mean(sla), 
            se_sla = sd(sla) / sqrt(n()), .groups = 'drop')

sla_controlpt_summary <- sla_na |> 
  group_by(ot, pt, pop, SAP_mm, SAT_C, cv) |> 
  filter(pt == 'control') |> 
  summarise(mean_sla = mean(sla), 
            se_sla = sd(sla) / sqrt(n()), .groups = 'drop')

sla_droughtpt_rxn <- ggplot(sla_droughtpt_summary, aes(x = ot, y = mean_sla, color = SAT_C, group = pop)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = mean_sla - se_sla, ymax = mean_sla + se_sla), width = 0.1) +
  geom_line(size = 1.5) +
  labs(y = "SLA", x = "offspring treatment", color = "SAT (C)", title = "drought parental treatment") +
  temperature_colors +
  ylim(67, 160) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(margin = margin(t = 4), color = "black", size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 4), color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10),  
        axis.text.y = element_text(color = "black", size = 10),  
        legend.title = element_text(color = "black", size = 10, face = "bold"),  
        legend.text = element_text(color = "black", size = 10),  
        plot.title = element_text(color = "black", size = 10, face = "bold", hjust = 0.5),  
        plot.subtitle = element_text(color = "black", size = 10),
        plot.caption = element_text(color = "black", size = 10))
sla_droughtpt_rxn

sla_controlpt_rxn <- ggplot(sla_controlpt_summary, aes(x = ot, y = mean_sla, color = SAT_C, group = pop)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = mean_sla - se_sla, ymax = mean_sla + se_sla), width = 0.1) +
  geom_line(size = 1.5) +
  labs(y = "SLA", x = "offspring treatment", color = "SAT (C)", title = "control parental treatment") +
  temperature_colors +
  ylim(67, 160) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(margin = margin(t = 4), color = "black", size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 4), color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10),  
        axis.text.y = element_text(color = "black", size = 10),  
        legend.title = element_text(color = "black", size = 10, face = "bold"),  
        legend.text = element_text(color = "black", size = 10),  
        plot.title = element_text(color = "black", size = 10, face = "bold", hjust = 0.5),  
        plot.subtitle = element_text(color = "black", size = 10),
        plot.caption = element_text(color = "black", size = 10))
sla_controlpt_rxn

sla_rxn_pop_fig <- sla_controlpt_rxn + sla_droughtpt_rxn +
  plot_layout(guides = "collect")
sla_rxn_pop_fig

ggsave("sla_rxn_pop_fig2.png", sla_rxn_pop_fig, path = "code-and-data/scripts/plpa-manuscript/figures",
       height = 4, width = 8, dpi = 300)


# LDMC rxn norm, with pop ------------------------------------------------
ldmc_na <- SLA_LDMC_FINAL |>
  filter(!is.na(ot) & !is.na(pt) & !is.na(ldmc))

ldmc_droughtpt_summary <- ldmc_na |> 
  group_by(ot, pt, pop, SAP_mm, SAT_C, cv) |> 
  filter(pt == 'drought') |> 
  summarise(mean_ldmc = mean(ldmc), 
            se_ldmc = sd(ldmc) / sqrt(n()), .groups = 'drop')

ldmc_controlpt_summary <- ldmc_na |> 
  group_by(ot, pt, pop, SAP_mm, SAT_C, cv) |> 
  filter(pt == 'control') |> 
  summarise(mean_ldmc = mean(ldmc), 
            se_ldmc = sd(ldmc) / sqrt(n()), .groups = 'drop')

ldmc_droughtpt_rxn <- ggplot(ldmc_droughtpt_summary, aes(x = ot, y = mean_ldmc, color = SAT_C, group = pop)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = mean_ldmc - se_ldmc, ymax = mean_ldmc + se_ldmc), width = 0.1) +
  geom_line(size = 1.5) +
  labs(y = "LDMC", x = "offspring treatment", color = "SAT (C)", title = "drought parental treatment") +
  temperature_colors +
  ylim(0.25, 0.8) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(margin = margin(t = 4), color = "black", size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 4), color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10),  
        axis.text.y = element_text(color = "black", size = 10),  
        legend.title = element_text(color = "black", size = 10, face = "bold"),  
        legend.text = element_text(color = "black", size = 10),  
        plot.title = element_text(color = "black", size = 10, face = "bold", hjust = 0.5),  
        plot.subtitle = element_text(color = "black", size = 10),
        plot.caption = element_text(color = "black", size = 10))
ldmc_droughtpt_rxn

ldmc_controlpt_rxn <- ggplot(ldmc_controlpt_summary, aes(x = ot, y = mean_ldmc, color = SAT_C, group = pop)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = mean_ldmc - se_ldmc, ymax = mean_ldmc + se_ldmc), width = 0.1) +
  geom_line(size = 1.5) +
  labs(y = "LDMC", x = "offspring treatment", color = "SAT (C)", title = "control parental treatment") +
  temperature_colors +
  ylim(0.25, 0.8) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(margin = margin(t = 4), color = "black", size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 4), color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10),  
        axis.text.y = element_text(color = "black", size = 10),  
        legend.title = element_text(color = "black", size = 10, face = "bold"),  
        legend.text = element_text(color = "black", size = 10),  
        plot.title = element_text(color = "black", size = 10, face = "bold", hjust = 0.5),  
        plot.subtitle = element_text(color = "black", size = 10),
        plot.caption = element_text(color = "black", size = 10))
ldmc_controlpt_rxn

ldmc_rxn_pop_fig <- ldmc_controlpt_rxn + ldmc_droughtpt_rxn +
  plot_layout(guides = "collect")
ldmc_rxn_pop_fig

ggsave("ldmc_rxn_pop_fig2.png", ldmc_rxn_pop_fig, path = "code-and-data/scripts/plpa-manuscript/figures",
       height = 4, width = 8, dpi = 300)

# shoot rxn norm, with pop ------------------------------------------------
shoot_na <- BIOMASS_FINAL |>
  filter(!is.na(ot) & !is.na(pt) & !is.na(shoot))

shoot_droughtpt_summary <- shoot_na |> 
  group_by(ot, pt, pop, SAP_mm, SAT_C, cv) |> 
  filter(pt == 'drought') |> 
  summarise(mean_shoot = mean(shoot), 
            se_shoot = sd(shoot) / sqrt(n()), .groups = 'drop')

shoot_controlpt_summary <- shoot_na |> 
  group_by(ot, pt, pop, SAP_mm, SAT_C, cv) |> 
  filter(pt == 'control') |> 
  summarise(mean_shoot = mean(shoot), 
            se_shoot = sd(shoot) / sqrt(n()), .groups = 'drop')

shoot_droughtpt_rxn <- ggplot(shoot_droughtpt_summary, aes(x = ot, y = mean_shoot, color = cv, group = pop)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = mean_shoot - se_shoot, ymax = mean_shoot + se_shoot), width = 0.1) +
  geom_line(size = 1.5) +
  labs(y = "shoot biomass (g)", x = "offspring treatment", color = "CV (%)", title = "drought parental treatment") +
  seasonality_colors +
  ylim(0, 0.7) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(margin = margin(t = 4), color = "black", size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 4), color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10),  
        axis.text.y = element_text(color = "black", size = 10),  
        legend.title = element_text(color = "black", size = 10, face = "bold"),  
        legend.text = element_text(color = "black", size = 10),  
        plot.title = element_text(color = "black", size = 10, face = "bold", hjust = 0.5),  
        plot.subtitle = element_text(color = "black", size = 10),
        plot.caption = element_text(color = "black", size = 10))
shoot_droughtpt_rxn

shoot_controlpt_rxn <- ggplot(shoot_controlpt_summary, aes(x = ot, y = mean_shoot, color = cv, group = pop)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = mean_shoot - se_shoot, ymax = mean_shoot + se_shoot), width = 0.1) +
  geom_line(size = 1.5) +
  labs(y = "shoot biomass (g)", x = "offspring treatment", color = "CV (%)", title = "control parental treatment") +
  seasonality_colors +
  ylim(0, 0.7) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(margin = margin(t = 4), color = "black", size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 4), color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10),  
        axis.text.y = element_text(color = "black", size = 10),  
        legend.title = element_text(color = "black", size = 10, face = "bold"),  
        legend.text = element_text(color = "black", size = 10),  
        plot.title = element_text(color = "black", size = 10, face = "bold", hjust = 0.5),  
        plot.subtitle = element_text(color = "black", size = 10),
        plot.caption = element_text(color = "black", size = 10))
shoot_controlpt_rxn

shootbiomass_rxn_pop_fig <- shoot_controlpt_rxn + shoot_droughtpt_rxn +
  plot_layout(guides = "collect")
shootbiomass_rxn_pop_fig

ggsave("shootbiomass_rxn_pop_fig2.png", shootbiomass_rxn_pop_fig, path = "code-and-data/scripts/plpa-manuscript/figures",
       height = 4, width = 8, dpi = 300)

# root rxn norm, with pop ------------------------------------------------
root_na <- BIOMASS_FINAL |>
  filter(!is.na(ot) & !is.na(pt) & !is.na(root))

root_droughtpt_summary <- root_na |> 
  group_by(ot, pt, pop, SAP_mm, SAT_C, cv) |> 
  filter(pt == 'drought') |> 
  summarise(mean_root = mean(root), 
            se_root = sd(root) / sqrt(n()), .groups = 'drop')

root_controlpt_summary <- root_na |> 
  group_by(ot, pt, pop, SAP_mm, SAT_C, cv) |> 
  filter(pt == 'control') |> 
  summarise(mean_root = mean(root), 
            se_root = sd(root) / sqrt(n()), .groups = 'drop')

root_droughtpt_rxn <- ggplot(root_droughtpt_summary, aes(x = ot, y = mean_root, color = SAP_mm, group = pop)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = mean_root - se_root, ymax = mean_root + se_root), width = 0.1) +
  geom_line(size = 1.5) +
  labs(y = "root biomass (g)", x = "offspring treatment", color = "SAP (mm)", title = "drought parental treatment") +
  rainfall_colors +
  ylim(0, 0.85) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(margin = margin(t = 4), color = "black", size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 4), color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10),  
        axis.text.y = element_text(color = "black", size = 10),  
        legend.title = element_text(color = "black", size = 10, face = "bold"),  
        legend.text = element_text(color = "black", size = 10),  
        plot.title = element_text(color = "black", size = 10, face = "bold", hjust = 0.5),  
        plot.subtitle = element_text(color = "black", size = 10),
        plot.caption = element_text(color = "black", size = 10))
root_droughtpt_rxn

root_controlpt_rxn <- ggplot(root_controlpt_summary, aes(x = ot, y = mean_root, color = SAP_mm, group = pop)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = mean_root - se_root, ymax = mean_root + se_root), width = 0.1) +
  geom_line(size = 1.5) +
  labs(y = "root biomass (g)", x = "offspring treatment", color = "SAP (mm)", title = "control parental treatment") +
  rainfall_colors +
  ylim(0, 0.85) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(margin = margin(t = 4), color = "black", size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 4), color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10),  
        axis.text.y = element_text(color = "black", size = 10),  
        legend.title = element_text(color = "black", size = 10, face = "bold"),  
        legend.text = element_text(color = "black", size = 10),  
        plot.title = element_text(color = "black", size = 10, face = "bold", hjust = 0.5),  
        plot.subtitle = element_text(color = "black", size = 10),
        plot.caption = element_text(color = "black", size = 10))
root_controlpt_rxn

rootbiomass_rxn_pop_fig <- root_controlpt_rxn + root_droughtpt_rxn +
  plot_layout(guides = "collect")
rootbiomass_rxn_pop_fig

ggsave("rootbiomass_rxn_pop_fig2.png", rootbiomass_rxn_pop_fig, path = "code-and-data/scripts/plpa-manuscript/figures",
       height = 4, width = 8, dpi = 300)

# days to flower rxn norm, with pop ------------------------------------------------
flower_na <- FLOWER_FINAL |>
  filter(!is.na(ot) & !is.na(pt) & !is.na(days_to_flower))

flower_droughtpt_summary <- flower_na |> 
  group_by(ot, pt, pop, SAP_mm, SAT_C, cv) |> 
  filter(pt == 'drought') |> 
  summarise(mean_flower = mean(days_to_flower), 
            se_flower = sd(days_to_flower) / sqrt(n()), .groups = 'drop')

flower_controlpt_summary <- flower_na |> 
  group_by(ot, pt, pop, SAP_mm, SAT_C, cv) |> 
  filter(pt == 'control') |> 
  summarise(mean_flower = mean(days_to_flower), 
            se_flower = sd(days_to_flower) / sqrt(n()), .groups = 'drop')

flower_droughtpt_rxn <- ggplot(flower_droughtpt_summary, aes(x = ot, y = mean_flower, color = SAP_mm, group = pop)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = mean_flower - se_flower, ymax = mean_flower + se_flower), width = 0.1) +
  geom_line(size = 1.5) +
  labs(y = "days to flowering", x = "offspring treatment", color = "SAP (mm)", title = "drought parental treatment") +
  rainfall_colors +
  ylim(30, 70) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(margin = margin(t = 4), color = "black", size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 4), color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10),  
        axis.text.y = element_text(color = "black", size = 10),  
        legend.title = element_text(color = "black", size = 10, face = "bold"),  
        legend.text = element_text(color = "black", size = 10),  
        plot.title = element_text(color = "black", size = 10, face = "bold", hjust = 0.5),  
        plot.subtitle = element_text(color = "black", size = 10),
        plot.caption = element_text(color = "black", size = 10))
flower_droughtpt_rxn

flower_controlpt_rxn <- ggplot(flower_controlpt_summary, aes(x = ot, y = mean_flower, color = SAP_mm, group = pop)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = mean_flower - se_flower, ymax = mean_flower + se_flower), width = 0.1) +
  geom_line(size = 1.5) +
  labs(y = "days to flowering", x = "offspring treatment", color = "SAP (mm)", title = "control parental treatment") +
  rainfall_colors +
  ylim(30, 70) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(margin = margin(t = 4), color = "black", size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 4), color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10),  
        axis.text.y = element_text(color = "black", size = 10),  
        legend.title = element_text(color = "black", size = 10, face = "bold"),  
        legend.text = element_text(color = "black", size = 10),  
        plot.title = element_text(color = "black", size = 10, face = "bold", hjust = 0.5),  
        plot.subtitle = element_text(color = "black", size = 10),
        plot.caption = element_text(color = "black", size = 10))
flower_controlpt_rxn

flower_rxn_pop_fig <- flower_controlpt_rxn + flower_droughtpt_rxn +
  plot_layout(guides = "collect")
flower_rxn_pop_fig

ggsave("flower_rxn_pop_fig2.png", flower_rxn_pop_fig, path = "code-and-data/scripts/plpa-manuscript/figures",
       height = 4, width = 8, dpi = 300)



# OLD 04_plasticity ------------------------------------------------------------

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

# calculate RDPI CC-DD of AG biomass per vpd
tgp_ag <- BIOMASS_FINAL |> 
  filter(TGP %in% c('DD', 'CC')) |> 
  na.omit(shoot) |> 
  mutate(vpd = as.factor(vpd)) |> 
  mutate(shoot = as.numeric(shoot)) |> 
  mutate(TGP = as.factor(TGP))

# rdpi for total biomass
tgp_ag_rdpi1 <- rdpi(tgp_ag, sp = vpd, trait = shoot, factor = TGP)

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

# adding in flowering status!
flower_status <- flower_status |> 
  rename(status_flower = status)

ag_plastic <- ag_plastic %>%
  left_join(flower_status %>% select(tray, tray_ab, id, pop_id, pop, `2021_treat`, `2023_treat`, TGP, trial, ot, pt, status_flower), 
            by = c("tray", "tray_ab", "id", "pop_id", "pop", "2021_treat", "2023_treat", "TGP", "trial", "ot", "pt"))

cor.test(ag_plastic$status_flower, ag_plastic$avg_rdpi, use = "complete.obs", method = "pearson")

# quick viz
ag_plastic_summary1 <- ag_plastic |> 
  group_by(vpd) |> 
  summarise(mean_flower = mean(status_flower, na.rm = TRUE),
            mean_rdpi = mean(avg_rdpi, na.rm = TRUE))

ggplot(ag_plastic_summary1, aes(x = mean_flower, y = mean_rdpi, color = vpd)) +
  geom_point(size = 3.5, alpha = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "black", alpha = 0.3) +
  labs(x = "Proportion flowered", y = "RDPI CC-DD total BM")





# RDPI of CC-DD root biomass ---------------------------------------------
#first, combine biomass dataframe and seed num dataframe
bg_plastic <- SEED_FINAL |> 
  left_join(BIOMASS_FINAL, by = c("tray", "tray_ab", "id", "pop_id", "pop", "2021_treat", "2023_treat", "TGP", "trial", "ot", "pt"))

bg_plastic <- bg_plastic %>%
  rename_with(~ gsub("\\.x$", "", .))

mort_day50 <- mort_day50 |> 
  #rename(TGP = TGP.x) |> 
  rename(trial = trial.x)

colnames(mort_day50)

bg_plastic <- bg_plastic %>%
  left_join(mort_day50 %>% select(tray, tray_ab, id, pop_id, pop, `2021_treat`, `2023_treat`, TGP, trial, ot, pt, status), 
            by = c("tray", "tray_ab", "id", "pop_id", "pop", "2021_treat", "2023_treat", "TGP", "trial", "ot", "pt"))

# calculate avg seed num per VPD group and add it to dataframe
avg_fitness <- bg_plastic |> 
  group_by(vpd) |> 
  summarise(avg_seeds = mean(num_total, na.rm = TRUE))

bg_plastic <- bg_plastic |> 
  group_by(vpd) |>
  mutate(avg_seed_num = mean(num_total, na.rm = TRUE))

# calculate RDPI CC-DD of AG biomass per vpd
tgp_bg <- BIOMASS_FINAL |> 
  filter(TGP %in% c('DD', 'CC')) |> 
  na.omit(root) |> 
  mutate(vpd = as.factor(vpd)) |> 
  mutate(root = as.numeric(root)) |> 
  mutate(TGP = as.factor(TGP))

# rdpi for total biomass
tgp_bg_rdpi1 <- rdpi(tgp_bg, sp = vpd, trait = root, factor = TGP)

avg_rdpi <- tgp_bg_rdpi1 |> 
  group_by(sp) |> 
  summarise(avg_rdpi = mean(rdpi, na.rm = TRUE)) |> 
  rename(vpd = sp)

# filtered biomass_final with only cc and dd
bg_plastic <- bg_plastic |> 
  filter(TGP %in% c('DD', 'CC'))

# attach avg RDPI to vpd values in filtered biomass dataframe
bg_plastic$vpd <- as.numeric(as.character(bg_plastic$vpd))
avg_rdpi$vpd <- as.numeric(as.character(avg_rdpi$vpd))

bg_plastic <- bg_plastic %>%
  left_join(avg_rdpi, by = "vpd")

# now, i think i can regress RDPI and seed num....
cor.test(bg_plastic$num_total, bg_plastic$avg_rdpi, use = "complete.obs", method = "pearson")
cor.test(bg_plastic$mass_total, bg_plastic$avg_rdpi, use = "complete.obs", method = "pearson")
cor.test(bg_plastic$status, bg_plastic$avg_rdpi, use = "complete.obs", method = "pearson")

# adding in flowering status!
flower_status <- flower_status |> 
  rename(status_flower = status)

bg_plastic <- bg_plastic %>%
  left_join(flower_status %>% select(tray, tray_ab, id, pop_id, pop, `2021_treat`, `2023_treat`, TGP, trial, ot, pt, status_flower), 
            by = c("tray", "tray_ab", "id", "pop_id", "pop", "2021_treat", "2023_treat", "TGP", "trial", "ot", "pt"))

cor.test(bg_plastic$status_flower, bg_plastic$avg_rdpi, use = "complete.obs", method = "pearson")

# quick viz
bg_plastic_summary1 <- bg_plastic |> 
  group_by(vpd) |> 
  summarise(mean_flower = mean(status_flower, na.rm = TRUE),
            mean_rdpi = mean(avg_rdpi, na.rm = TRUE))

ggplot(bg_plastic_summary1, aes(x = mean_flower, y = mean_rdpi, color = vpd)) +
  geom_point(size = 3.5, alpha = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "black", alpha = 0.3) +
  labs(x = "Proportion flowered", y = "RDPI CC-DD total BM")


# OLD - FITNESS ROBUSTNESS ---------------------------------------
# r:s ratio
#first, combine r:s dataframe and seed num dataframe
rs_plastic <- SEED_FINAL |> 
  left_join(BIOMASS_FINAL, by = c("tray", "tray_ab", "id", "pop_id", "pop", "2021_treat", "2023_treat", "TGP", "trial", "ot", "pt"))

# calculate avg fitness per treatment group
avg_fitness <- rs_plastic |> 
  group_by(TGP) |> 
  summarise(avg_seeds = mean(num_total, na.rm = TRUE))

# find the TGP# find the maximum average fitness
max_avg_fitness <- max(avg_fitness$avg_seeds, na.rm = TRUE)

# calculate fitness robustness for each plant
rs_plastic <- rs_plastic %>%
  mutate(fitness_robustness = num_total / max_avg_fitness)

print(rs_plastic)

# calculate the mean and SD of RS for each treatment group
rs_stats <- rs_plastic |> 
  filter(!is.na(ratio_RS)) |> 
  group_by(TGP) |> 
  summarise(mean_rs = mean(ratio_RS, na.rm = TRUE),
            sd_rs = sd(ratio_RS, na.rm = TRUE))

# calc CV for each treatment group
rs_stats <- rs_stats |> 
  mutate(cv_rs = sd_rs / mean_rs)

# join the cv back to the original df by TGP
rs_plastic <- rs_plastic %>%
  filter(!is.na(ratio_RS)) |>
  left_join(rs_stats %>% select(TGP, cv_rs), by = "TGP")

# correlation test
cor.test(rs_plastic$cv_rs, rs_plastic$fitness_robustness, use = "complete.obs", method = "pearson")

rs_plastic_summary <- rs_plastic |> 
  group_by(TGP) |> 
  summarise(mean_cv_rs = mean(cv_rs, na.rm = TRUE),
            se_cv_rs = sd(cv_rs, na.rm = TRUE) / sqrt(n()),
            mean_fitness_robustness = mean(fitness_robustness, na.rm = TRUE),
            se_fitness_robustness = sd(fitness_robustness, na.rm = TRUE) / sqrt(n()))

# visualize
ggplot(rs_plastic_summary, aes(x = mean_cv_rs, y = mean_fitness_robustness, color = TGP)) +
  geom_point(size = 3.5, alpha = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "black", alpha = 0.3) +
  labs(x = "CV of R:S ratio", y = "fitness robustness") +
  scale_color_manual(values = c("CC" = "blue3",
                                "CD" = "chartreuse4",
                                "DC" = "orange1",
                                "DD" = "red3")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(color = "black", size = 10, face = "bold"),
        axis.title.y = element_text(color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        legend.title = element_text(color = "black", size = 10, face = "bold"),
        legend.text = element_text(color = "black", size = 10))

# total biomass
#first, combine total biomass dataframe and seed num dataframe
total_plastic <- SEED_FINAL |> 
  left_join(BIOMASS_FINAL, by = c("tray", "tray_ab", "id", "pop_id", "pop", "2021_treat", "2023_treat", "TGP", "trial", "ot", "pt"))

# calculate avg fitness per treatment group
avg_fitness <- total_plastic |> 
  group_by(TGP) |> 
  summarise(avg_seeds = mean(num_total, na.rm = TRUE))

# find the TGP# find the maximum average fitness
max_avg_fitness <- max(avg_fitness$avg_seeds, na.rm = TRUE)

# calculate fitness robustness for each plant
total_plastic <- total_plastic %>%
  mutate(fitness_robustness = num_total / max_avg_fitness)

print(total_plastic)

# calculate the mean and SD of RS for each treatment group
total_stats <- total_plastic |> 
  filter(!is.na(total_biomass)) |> 
  group_by(TGP) |> 
  summarise(mean_total = mean(total_biomass, na.rm = TRUE),
            sd_total = sd(total_biomass, na.rm = TRUE))

# calc CV for each treatment group
total_stats <- total_stats |> 
  mutate(cv_total = sd_total / mean_total)

# join the cv back to the original df by TGP
total_plastic <- total_plastic %>%
  filter(!is.na(total_biomass)) |>
  left_join(total_stats %>% select(TGP, cv_total), by = "TGP")

# correlation test
cor.test(total_plastic$cv_total, total_plastic$fitness_robustness, use = "complete.obs", method = "pearson")

total_plastic_summary <- total_plastic |> 
  group_by(TGP) |> 
  summarise(mean_cv_total = mean(cv_total, na.rm = TRUE),
            se_cv_total = sd(cv_total, na.rm = TRUE) / sqrt(n()),
            mean_fitness_robustness = mean(fitness_robustness, na.rm = TRUE),
            se_fitness_robustness = sd(fitness_robustness, na.rm = TRUE) / sqrt(n()))

# visualize
ggplot(total_plastic_summary, aes(x = mean_cv_total, y = mean_fitness_robustness, color = TGP)) +
  geom_point(size = 3.5, alpha = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "black", alpha = 0.3) +
  labs(x = "CV of total biomass", y = "fitness robustness") +
  scale_color_manual(values = c("CC" = "blue3",
                                "CD" = "chartreuse4",
                                "DC" = "orange1",
                                "DD" = "red3")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(color = "black", size = 10, face = "bold"),
        axis.title.y = element_text(color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        legend.title = element_text(color = "black", size = 10, face = "bold"),
        legend.text = element_text(color = "black", size = 10))


# mortality
#first, combine mort dataframe and seed num dataframe
mort_plastic <- SEED_FINAL |> 
  left_join(mort_day50, by = c("tray", "tray_ab", "id", "pop_id", "pop", "2021_treat", "2023_treat", "TGP", "trial", "ot", "pt"))

# calculate avg fitness per treatment group
avg_fitness <- mort_plastic |> 
  group_by(TGP) |> 
  summarise(avg_seeds = mean(num_total, na.rm = TRUE))

# find the TGP# find the maximum average fitness
max_avg_fitness <- max(avg_fitness$avg_seeds, na.rm = TRUE)

# calculate fitness robustness for each plant
mort_plastic <- mort_plastic %>%
  mutate(fitness_robustness = num_total / max_avg_fitness)

print(total_plastic)

# calculate the mean and SD of RS for each treatment group
mort_stats <- mort_plastic |> 
  filter(!is.na(status)) |> 
  group_by(TGP) |> 
  summarise(mean_mort = mean(status, na.rm = TRUE),
            sd_mort = sd(status, na.rm = TRUE))

# calc CV for each treatment group
mort_stats <- mort_stats |> 
  mutate(cv_mort = sd_mort / mean_mort)

# join the cv back to the original df by TGP
mort_plastic <- mort_plastic %>%
  filter(!is.na(status)) |>
  left_join(mort_stats %>% select(TGP, cv_mort), by = "TGP")

# correlation test
cor.test(mort_plastic$cv_mort, mort_plastic$fitness_robustness, use = "complete.obs", method = "pearson")

mort_plastic_summary <- mort_plastic |> 
  group_by(TGP) |> 
  summarise(mean_cv_mort = mean(cv_mort, na.rm = TRUE),
            se_cv_mort = sd(cv_mort, na.rm = TRUE) / sqrt(n()),
            mean_fitness_robustness = mean(fitness_robustness, na.rm = TRUE),
            se_fitness_robustness = sd(fitness_robustness, na.rm = TRUE) / sqrt(n()))

# visualize
ggplot(mort_plastic_summary, aes(x = mean_cv_mort, y = mean_fitness_robustness, color = TGP)) +
  geom_point(size = 3.5, alpha = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "black", alpha = 0.3) +
  labs(x = "CV of mortality", y = "fitness robustness") +
  scale_color_manual(values = c("CC" = "blue3",
                                "CD" = "chartreuse4",
                                "DC" = "orange1",
                                "DD" = "red3")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(color = "black", size = 10, face = "bold"),
        axis.title.y = element_text(color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10),
        axis.text.y = element_text(color = "black", size = 10),
        legend.title = element_text(color = "black", size = 10, face = "bold"),
        legend.text = element_text(color = "black", size = 10))







# OLD ANALYSIS - rdpi of TGP and measure of adaptivity -----------
# here, for traits that are significant for PT or PT x pop

tgp_biomass1 <- BIOMASS_FINAL |> 
  filter(TGP %in% c('DD', 'CC')) |> 
  filter(pop == 1) |> 
  na.omit(shoot)

tgp_height1 <- HEIGHT_FINAL |> 
  filter(TGP %in% c('DD', 'CC')) |> 
  filter(pop == 1) |> 
  na.omit(max)

# rdpi for shoot biomass, rdpi for max height = two plasticity indices for two traits for one pop
tgp_shoot_rdpi1 <- rdpi(tgp_biomass1, sp = pop, trait = shoot, factor = TGP)
tgp_max_rdpi1 <- rdpi(tgp_height1, sp = pop, trait = max, factor = TGP)


# if i combine these dataframes....
total_rdpi1 <- rbind(tgp_shoot_rdpi1, tgp_max_rdpi1)

# i can find the mean and SE/SD and i think that will give me an OVERALL TGP value for pop 1 ?
mean_rdpi <- mean(total_rdpi1$rdpi, na.rm =TRUE)
se_rdpi <- sd(total_rdpi1$rdpi, na.rm =TRUE)/sqrt(length(total_rdpi1$rdpi))
mean_rdpi
se_rdpi

# or, i think i can use the combined RDPI values in the final dataframe to set up
# TGP rdpi of pop. 1 ~ mort + growth + fecundity (of pop 1?)

# but then, if i am just combining the df, some traits will have a stronger influence
# eg, i simply have more data on shoot biomass than i do number of days flowering, because VERY few flowered.
# and also, for my model, pop 1 had more plants germ than pop 3, so the df rows will be different for every pop?

# also, are my x variables the mort, growth, and fecudnity of just pop 1?

# should i include ALL TRAITS measured when calculated total rdpi pop1?













# OLD ANALYSIS -rdpi of WGP/TGP -------------------------------------------

temperature_colors <- scale_color_gradient2(low = "blue3", mid = "purple3", high = "red3", midpoint = 16.5)

rainfall_colors <- scale_color_gradient2(low = "yellow1", mid = "palegreen3", high = "darkgreen", midpoint = 50)

seasonality_colors <- scale_color_gradient2(low = "chartreuse4", mid ="orange2", high = "violetred", midpoint = 32.5)

# shoot biomass rdpi
shoot_rdpi <- rdpi(BIOMASS_FINAL, sp = pop, trait = shoot, factor = TGP)
# i think i need to filter the df to split it up between WGP and TGP pair comparisons!!
# and then run the rdpi fxn
# rn idk what it's comparing. across all four treatments?? idk

# figure where i compare WGP (DD / DC) and TGP (DC / CC) which i think makes more sense to me
tgp_biomass <- BIOMASS_FINAL |> 
  filter(TGP %in% c('DD', 'CC')) |> 
  filter(pop == 1)

tgp_biomass <- BIOMASS_FINAL |> 
  filter(TGP %in% c('DC', 'CC'))

# fancy package
wgp_shoot_rdpi <- rdpi(wgp_biomass, sp = pop, trait = shoot, factor = TGP)

tgp_shoot_rdpi <- rdpi(tgp_biomass, sp = pop, trait = shoot, factor = TGP)


# add climate data for lm later
wgp_shoot_rdpi <- wgp_shoot_rdpi |> 
  rename(pop = sp)
wgp_shoot_rdpi <- full_join(wgp_shoot_rdpi, climate_tidy, by = c("pop" = "population"))

tgp_shoot_rdpi <- tgp_shoot_rdpi |> 
  rename(pop = sp)
tgp_shoot_rdpi <- full_join(tgp_shoot_rdpi, climate_tidy, by = c("pop" = "population"))

# visualize shoot plasticity - clean up
wgp_shoot_summary <- wgp_shoot_rdpi |> 
  group_by(sp) |> 
  summarise(mean_wgp = mean(rdpi), 
            se_wgp = sd(rdpi) / sqrt(n()), .groups = 'drop') |> 
  rename(pop = sp)

tgp_shoot_summary <- tgp_shoot_rdpi |> 
  group_by(sp) |> 
  summarise(mean_tgp = mean(rdpi), 
            se_tgp = sd(rdpi) / sqrt(n()), .groups = 'drop') |> 
  rename(pop = sp)

# join data
shoot_plasticity_summary1 <-full_join(wgp_shoot_summary, tgp_shoot_summary)
shoot_plasticity_summary <-full_join(shoot_plasticity_summary1, climate_tidy, by = c("pop" = "population"))

# reshape data for graph
shoot_plasticity_summary_long <- shoot_plasticity_summary |> 
  pivot_longer(
    cols = c(mean_wgp, mean_tgp),
    names_to = "type",
    values_to = "mean"
  ) |> 
  pivot_longer(
    cols = c(se_wgp, se_tgp),
    names_to = "se_type",
    values_to = "se"
  ) |> 
  filter((type == "mean_wgp" & se_type == "se_wgp") |
           (type == "mean_tgp" & se_type == "se_tgp")) |> 
  mutate(type = ifelse(type == "mean_wgp", "wgp", "tgp"))

shoot_plasticity_summary_long$type <- factor(shoot_plasticity_summary_long$type, levels = c("wgp", "tgp"))

# visualize shoot plasticity - graph
shoot_plasticity <- ggplot(shoot_plasticity_summary_long, aes(x = type, y = mean, group = pop, color = SAT_C)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.1) +
  geom_line(size = 1) +
  labs(y = "RDPI for shoot biomass", x = " ", color = "SAT (C)", title = " ") +
  temperature_colors +
  ylim(0.2, 1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(family = "Helvetica"),
        axis.title.x = element_text(margin = margin(t = 4), color = "black", size = 10, face = "bold"),
        axis.title.y = element_text(margin = margin(r = 4), color = "black", size = 10, face = "bold"),
        axis.text.x = element_text(color = "black", size = 10),  
        axis.text.y = element_text(color = "black", size = 10),  
        legend.title = element_text(color = "black", size = 10, face = "bold"),  
        legend.text = element_text(color = "black", size = 10),  
        plot.title = element_text(color = "black", size = 10, face = "bold", hjust = 0.5),  
        plot.subtitle = element_text(color = "black", size = 10),
        plot.caption = element_text(color = "black", size = 10))
shoot_plasticity

# save it
ggsave("shoot_plasticity_fig3.png", shoot_plasticity, path = "code-and-data/scripts/plpa-manuscript/figures",
       height = 4, width = 8, dpi = 300)

# stats
m1 <- lm(rdpi ~  cv, data = wgp_shoot_rdpi)
m2 <- lm(rdpi ~  SAT_C , data = wgp_shoot_rdpi)
m3 <- lm(rdpi ~  SAP_mm , data = wgp_shoot_rdpi)

m4 <- lm(rdpi ~  cv , data = tgp_shoot_rdpi)
m5 <- lm(rdpi ~  SAT_C , data = tgp_shoot_rdpi)
m6 <- lm(rdpi ~  SAP_mm , data = tgp_shoot_rdpi)

summary(m1)
summary(m2)
summary(m3)

summary(m4)
summary(m5)
summary(m6)

plot(residuals(m))
qqnorm(residuals(m))

plot(residuals(m1))
qqnorm(residuals(m1))

plot(residuals(m2))
qqnorm(residuals(m2))

shoot_plasticity_summary |> ggplot(aes(x = SAT_C, y = mean_wgp)) + geom_point() + geom_smooth( method = 'lm')
shoot_plasticity_summary |> ggplot(aes(x = SAT_C, y = mean_tgp)) + geom_point() + geom_smooth( method = 'lm')

# alternatively, pearson correlation?
cor_data <- wgp_shoot_rdpi %>%
  select(rdpi, SAP_mm, SAT_C, cv)  # Add other bioclimatic variables as needed
cor_matrix <- cor(cor_data, use = "complete.obs", method = "pearson")
print(cor_matrix)

cor_data <- tgp_shoot_rdpi %>%
  select(rdpi, SAP_mm, SAT_C, cv)  # Add other bioclimatic variables as needed
cor_matrix <- cor(cor_data, use = "complete.obs", method = "pearson")
print(cor_matrix)


# RDPI of CC-DD total biomass ---------------------------------------------
#first, combine total biomass dataframe and seed num dataframe
BIOMASS_FINAL <- BIOMASS_FINAL |> 
  filter(!(tray == 9 & tray_ab == "B" & id == "J2"))
MORT_DAY50 <- MORT_DAY50 |> 
  filter(!(tray == 9 & tray_ab == "B" & id == "J2"))

total_plastic <- SEED_FINAL |> 
  left_join(BIOMASS_FINAL, by = c("tray", "tray_ab", "id", "pop_id", "pop", "2021_treat", "2023_treat", "TGP", "trial", "ot", "pt"))

MORT_DAY50 <- MORT_DAY50 |> 
  #rename(TGP = TGP.x) |> 
  #rename(trial = trial.x) |> 
  rename(spring_vpd_cv = spring_vpd_cv.x)

colnames(MORT_DAY50)
colnames(total_plastic)

#total_plastic <- total_plastic %>%
#dplyr::select(-ends_with(".x"))  # Drop all columns ending with .x

total_plastic <- total_plastic %>%
  left_join(MORT_DAY50 %>% dplyr::select(tray, tray_ab, id, pop_id, pop, `2021_treat`, `2023_treat`, TGP, trial, ot, pt, status), 
            by = c("tray", "tray_ab", "id", "pop_id", "pop", "2021_treat", "2023_treat", "TGP", "trial", "ot", "pt"))

total_plastic <- total_plastic |> 
  #rename(TGP = TGP.x)
  #rename(trial = trial.x) |> 
  rename(spring_vpd_cv = spring_vpd_cv.x)

# calculate avg seed num per VPD group and add it to dataframe
avg_fitness <- total_plastic |> 
  group_by(spring_vpd_cv) |> 
  summarise(avg_seeds = mean(num_total, na.rm = TRUE))

total_plastic <- total_plastic |> 
  group_by(spring_vpd_cv) |>
  mutate(avg_seed_num = mean(num_total, na.rm = TRUE))

# calculate RDPI CC-DD of total biomass per vpd
tgp_biomass <- BIOMASS_FINAL |> 
  filter(TGP %in% c('DD', 'CC')) |> 
  na.omit(total_biomass) |> 
  mutate(spring_vpd_cv = as.factor(spring_vpd_cv)) |> 
  mutate(total_biomass = as.numeric(total_biomass)) |> 
  mutate(TGP = as.factor(TGP))

# rdpi for total biomass
tgp_total_rdpi1 <- rdpi(tgp_biomass, sp = spring_vpd_cv, trait = total_biomass, factor = TGP)

avg_rdpi <- tgp_total_rdpi1 |> 
  group_by(sp) |> 
  summarise(avg_rdpi = mean(rdpi, na.rm = TRUE)) |> 
  rename(spring_vpd_cv = sp)

# filtered biomass_final with only cc and dd
total_plastic <- total_plastic |> 
  filter(TGP %in% c('DD', 'CC'))


# attach avg RDPI to vpd values in filtered biomass dataframe
total_plastic$spring_vpd_cv <- as.numeric(as.character(total_plastic$spring_vpd_cv))
avg_rdpi$spring_vpd_cv <- as.numeric(as.character(avg_rdpi$spring_vpd_cv))

total_plastic <- total_plastic %>%
  left_join(avg_rdpi, by = "spring_vpd_cv")

# adding in flowering status!
flower_status <- FLOWER_STATUS |> 
  rename(status_flower = status)

total_plastic <- total_plastic %>%
  left_join(flower_status %>% select(tray, tray_ab, id, pop_id, pop, `2021_treat`, `2023_treat`, TGP, trial, ot, pt, status_flower), 
            by = c("tray", "tray_ab", "id", "pop_id", "pop", "2021_treat", "2023_treat", "TGP", "trial", "ot", "pt"))

# now, i think i can regress RDPI and seed num....
cor.test(total_plastic$num_total, total_plastic$avg_rdpi, use = "complete.obs", method = "pearson")
cor.test(total_plastic$mass_total, total_plastic$avg_rdpi, use = "complete.obs", method = "pearson")
cor.test(total_plastic$status, total_plastic$avg_rdpi, use = "complete.obs", method = "pearson")
cor.test(total_plastic$status_flower, total_plastic$avg_rdpi, use = "complete.obs", method = "pearson")
cor.test(total_plastic$spring_vpd_cv, total_plastic$avg_rdpi, use = "complete.obs", method = "pearson")


# RDPI of CC-DD R:S ratio ---------------------------------------------
#first, combine total biomass dataframe and seed num dataframe
rs_plastic <- SEED_FINAL |> 
  left_join(BIOMASS_FINAL, by = c("tray", "tray_ab", "id", "pop_id", "pop", "2021_treat", "2023_treat", "TGP", "trial", "ot", "pt"))

rs_plastic <- rs_plastic %>%
  rename_with(~ gsub("\\.x$", "", .))

MORT_DAY50 <- MORT_DAY50 |> 
  #rename(TGP = TGP.x) |> 
  rename(trial = trial.x)

colnames(MORT_DAY50)

rs_plastic <- rs_plastic %>%
  left_join(MORT_DAY50 %>% select(tray, tray_ab, id, pop_id, pop, `2021_treat`, `2023_treat`, TGP, trial, ot, pt, status), 
            by = c("tray", "tray_ab", "id", "pop_id", "pop", "2021_treat", "2023_treat", "TGP", "trial", "ot", "pt"))

# calculate RDPI CC-DD of total biomass per vpd
tgp_rs <- BIOMASS_FINAL |> 
  filter(TGP %in% c('DD', 'CC')) |> 
  na.omit(ratio_RS) |> 
  mutate(spring_vpd_cv = as.factor(spring_vpd_cv)) |> 
  mutate(ratio_RS = as.numeric(ratio_RS)) |> 
  mutate(TGP = as.factor(TGP))

# rdpi for root to shoot ratio
tgp_rs_rdpi1 <- rdpi(tgp_rs, sp = spring_vpd_cv, trait = ratio_RS, factor = TGP)

avg_rdpi <- tgp_rs_rdpi1 |> 
  group_by(sp) |> 
  summarise(avg_rdpi = mean(rdpi, na.rm = TRUE)) |> 
  rename(spring_vpd_cv = sp)

# filtered biomass_final with only cc and dd
rs_plastic <- rs_plastic |> 
  filter(TGP %in% c('DD', 'CC'))

# attach avg RDPI to vpd values in filtered biomass dataframe
rs_plastic$spring_vpd_cv <- as.numeric(as.character(rs_plastic$spring_vpd_cv))
avg_rdpi$spring_vpd_cv <- as.numeric(as.character(avg_rdpi$spring_vpd_cv))

rs_plastic <- rs_plastic %>%
  left_join(avg_rdpi, by = "spring_vpd_cv")

# adding in flowering status!
flower_status <- flower_status |> 
  rename(status_flower = status)

rs_plastic <- rs_plastic %>%
  left_join(flower_status %>% select(tray, tray_ab, id, pop_id, pop, `2021_treat`, `2023_treat`, TGP, trial, ot, pt, status_flower), 
            by = c("tray", "tray_ab", "id", "pop_id", "pop", "2021_treat", "2023_treat", "TGP", "trial", "ot", "pt"))


# now, i regress....
cor.test(rs_plastic$num_total, rs_plastic$avg_rdpi, use = "complete.obs", method = "pearson")
cor.test(rs_plastic$mass_total, rs_plastic$avg_rdpi, use = "complete.obs", method = "pearson")
cor.test(rs_plastic$status, rs_plastic$avg_rdpi, use = "complete.obs", method = "pearson")
cor.test(rs_plastic$status_flower, rs_plastic$avg_rdpi, use = "complete.obs", method = "pearson")
cor.test(rs_plastic$spring_vpd_cv, rs_plastic$avg_rdpi, use = "complete.obs", method = "pearson")

