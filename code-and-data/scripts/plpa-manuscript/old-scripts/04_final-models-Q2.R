# title: 04_final-models-population

# about: in the total models, many trait responses differed by population 
# here, i explore how seed source location might affect trait response
# traits: root, shoot, total, RS, height, RGR, SLA, LDMC mort, seed num, days to flower

# author: madeleine wallace

# PACKAGES -----------------------------------------------------------------
library(tidyverse)
library(lme4) #model
library(lmerTest) #model
library(nlme) #model with variance functions if i need it
library(DHARMa) #model diagnostics
library(performance) #model diagnostics
library(car) #testies
library(MuMIn) #r2
library(MASS) #glm.nb
library(emmeans) #contrasts
library(sandwich) #robust standard errors
library(lmtest) #robust standard errors
library(glmmTMB)
library(optimx)
source('code-and-data/scripts/plpa-manuscript/01_clean.R')

add_stars <- function(p) {
  if (p < 0.001) {
    "***"
  } else if (p < 0.01) {
    "**"
  } else if (p < 0.05) {
    "*"
  } else if (p < 0.1) {
    "."
  } else {
    ""
  }
}

# NEW MODELS - incl vpd and cv! -------------------------------------------

# rachel meeting - rather than SAP, SAT, and CV, we can just use VPD and CV
# VPD = combo of SAT and SAP
# high VPD, high CV = more stressful, low VPD and low CV = less stressful
# gets to the driver of TGP; both confirmation and drivers across the geographic range
# redo models with CV and VPD
# remove CV - correlated with VPD
ggplot(climate_tidy, aes(x = vpd, y = cv)) +
  geom_point() +                      
  geom_smooth(method = "lm", se = TRUE, color = "blue") +         
  labs(x = "vpd", y = "cv")
model1 <- lm(vpd ~ cv, data = climate_tidy)

# also, change vpd from hPa to kPa to help with model functioning
# 2-3 kPa vs 20-30 hPa; kPa is alot closer to factor levels with 2 options (0/1)

# CURRENT ANALYSIS ------------------------------------------------------------
# VPD ---------------------------------------------------------------------
# root models -------------------------------------------------------------
BIOMASS_FINAL <- BIOMASS_FINAL %>%
  mutate(vpd.k = vpd / 10)

root_model <- lmer(root ~ ot + pt + vpd.k + 
                     ot*vpd.k + pt*vpd.k + ot*pt + ot*pt*vpd.k + (1 | tray),
                   data = BIOMASS_FINAL,
                   na.action = na.exclude)

# assess
vif(root_model)
plot(BIOMASS_FINAL$ot, resid(root_model), xlab = "ot", ylab = "residuals")
plot(BIOMASS_FINAL$tray, resid(root_model), xlab = "tray", ylab = "residuals")
plot(root_model, add.smooth = FALSE, which = 1)
E <- resid(root_model)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(root_model))
qqline(resid(root_model))

summary(root_model)
Anova(root_model, type="III",  test.statistic= "F")
r.squaredGLMM(root_model)

#contrasts
rootcontrast <- emmeans(root_model, specs = pairwise ~ vpd:ot:pt)
rootcontrast$contrasts

# shoot models ------------------------------------------------------------
#same structure as OG model - tray has a crazy amount of variance btwn groups
shoot_model <- lme(shoot ~ ot + pt + vpd.k + 
                            ot*vpd.k +
                            pt*vpd.k +
                            ot*pt +
                            ot*pt*vpd.k,
                            random = ~ 1 | tray,
                            weights = varIdent(form = ~ 1 | tray),
                            data = BIOMASS_FINAL,
                            na.action = na.exclude,
                            control = lmeControl(opt = "optim", maxIter = 100, msMaxIter = 100))

# assess
vif(shoot_model)
plot(BIOMASS_FINAL$ot, resid(shoot_model), xlab = "ot", ylab = "residuals")
plot(BIOMASS_FINAL$tray, resid(shoot_model), xlab = "tray", ylab = "residuals")
check_model(shoot_model)
plot(shoot_model, add.smooth = FALSE, which = 1)
E <- resid(shoot_model)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(shoot_model))
qqline(resid(shoot_model))

summary(shoot_model)
shoot_anova <- anova(shoot_model, type='marginal')
r.squaredGLMM(shoot_model)

shoot_anova$Significance <- sapply(shoot_anova$`p-value`, add_stars)
print(shoot_anova)

#contrasts
shootcontrast <- emmeans(shoot_model2, specs = pairwise ~ ot:pt)
shootcontrast$contrasts

# total models ------------------------------------------------------------
total_model <- lmer(total_biomass ~ ot + pt + vpd.k + 
                      ot*vpd.k + pt*vpd.k + ot*pt + ot*pt*vpd.k + (1 | tray),
                   data = BIOMASS_FINAL,
                   na.action = na.exclude)

summary(total_model)

# assess
vif(total_model)
plot(BIOMASS_FINAL$ot, resid(total_model), xlab = "ot", ylab = "residuals")
plot(BIOMASS_FINAL$tray, resid(total_model), xlab = "tray", ylab = "residuals")
plot(total_model, add.smooth = FALSE, which = 1)
E <- resid(total_model)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(total_model))
qqline(resid(total_model))              

summary(total_model)
Anova(total_model, type="III",  test.statistic= "F", ddf="Kenward-Roger")
r.squaredGLMM(total_model)

#contrasts
totalcontrast <- emmeans(total_model, specs = pairwise ~ ot:pt)
totalcontrast$contrasts

# RS models -------------------------------------------------------------
rs_model <- lmer(log_RS ~ ot + pt + vpd.k + 
                   ot*vpd.k + pt*vpd.k + ot*pt + ot*pt*vpd.k + (1 | tray),
                 data = BIOMASS_FINAL,
                 na.action = na.exclude)

summary(rs_model)

# assess
vif(rs_model)
plot(BIOMASS_FINAL$ot, resid(rs_model), xlab = "ot", ylab = "residuals")
plot(BIOMASS_FINAL$tray, resid(rs_model), xlab = "tray", ylab = "residuals")
plot(rs_model, add.smooth = FALSE, which = 1)
E <- resid(rs_model)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(rs_model))
qqline(resid(rs_model))              

summary(rs_model)
Anova(rs_model, type="III",  test.statistic= "F", ddf="Kenward-Roger")
r.squaredGLMM(rs_model)

#contrasts
rscontrast <- emmeans(rs_model, specs = pairwise ~ ot:pt)
rscontrast$contrasts

# height models -------------------------------------------------------------
HEIGHT_FINAL <- HEIGHT_FINAL %>%
  mutate(vpd.k = vpd / 10)

height_model <- lme(max ~ ot + pt + vpd.k + 
                      ot*vpd.k +
                      pt*vpd.k +
                      ot*pt +
                      ot*pt*vpd.k,
                      random = ~ 1 | tray,
                      weights = varIdent(form = ~ 1 | tray),
                      data = HEIGHT_FINAL,
                      control = lmeControl(opt = "optim", maxIter = 100, msMaxIter = 100),
                      na.action = na.exclude)

# assess
plot(HEIGHT_FINAL$ot, resid(height_model), xlab = "ot", ylab = "residuals")
plot(HEIGHT_FINAL$tray, resid(height_model), xlab = "tray", ylab = "residuals")
plot(height_model, add.smooth = FALSE, which = 1)
E <- resid(height_model)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(height_model))
qqline(resid(height_model))              

summary(height_model)
height_anova <- anova(height_model, type='marginal')
r.squaredGLMM(height_model)

height_anova$Significance <- sapply(height_anova$`p-value`, add_stars)
print(height_anova)

#contrasts
heightcontrast <- emmeans(height_model, specs = pairwise ~ ot:pt)
heightcontrast$contrasts

# rgr models -------------------------------------------------------------
RGR_FINAL <- RGR_FINAL %>%
  mutate(vpd.k = vpd / 10)

rgr_model <- lme(RGR ~ ot + pt + vpd.k + 
                   ot*vpd.k +
                   pt*vpd.k +
                   ot*pt +
                   ot*pt*vpd.k,
                    random = ~ 1 | tray,
                    weights = varIdent(form = ~ 1 | tray),
                    data = RGR_FINAL,
                    control = lmeControl(opt = "optim", maxIter = 100, msMaxIter = 100),
                    na.action = na.exclude)

# assess
plot(RGR_FINAL$ot, resid(rgr_model), xlab = "ot", ylab = "residuals")
plot(RGR_FINAL$tray, resid(rgr_model), xlab = "tray", ylab = "residuals")
plot(rgr_model, add.smooth = FALSE, which = 1)
E <- resid(rgr_model)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(rgr_model))
qqline(resid(rgr_model))              

summary(rgr_model)
rgr_anova <- anova(rgr_model, type='marginal')
r.squaredGLMM(rgr_model)

rgr_anova$Significance <- sapply(rgr_anova$`p-value`, add_stars)
print(rgr_anova)

#contrasts
rgrcontrast <- emmeans(rgr_model, specs = pairwise ~ ot:pt)
rgrcontrast$contrasts

# sla models --------------------------------------------------------------
SLA_LDMC_FINAL <- SLA_LDMC_FINAL %>%
  mutate(vpd.k = vpd / 10)

sla_model <- lmer(log(sla) ~ ot + pt + vpd.k + 
                    ot*vpd.k + pt*vpd.k + ot*pt + ot*pt*vpd.k + (1 | tray),
                   data = SLA_LDMC_FINAL,
                   na.action = na.exclude)

# assess
plot(SLA_LDMC_FINAL$ot, resid(sla_model), xlab = "ot", ylab = "residuals")
plot(SLA_LDMC_FINAL$tray, resid(sla_model), xlab = "tray", ylab = "residuals")
plot(sla_model, add.smooth = FALSE, which = 1)
E <- resid(sla_model)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(sla_model))
qqline(resid(sla_model))              

summary(sla_model)
Anova(sla_model, type="III",  test.statistic= "F", ddf="Kenward-Roger")
r.squaredGLMM(sla_model)

# ldmc models --------------------------------------------------------------
ldmc_model <- lmer(log(ldmc) ~ ot + pt + vpd.k + 
                     ot*vpd.k + pt*vpd.k + ot*pt + ot*pt*vpd.k + (1 | tray),
                  data = SLA_LDMC_FINAL,
                  na.action = na.exclude)

# assess
plot(SLA_LDMC_FINAL$ot, resid(ldmc_model), xlab = "ot", ylab = "residuals")
plot(SLA_LDMC_FINAL$tray, resid(ldmc_model), xlab = "tray", ylab = "residuals")
plot(ldmc_model, add.smooth = FALSE, which = 1)
E <- resid(ldmc_model)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(ldmc_model))
qqline(resid(ldmc_model))              

summary(ldmc_model)
Anova(ldmc_model, type="III",  test.statistic= "F", ddf="Kenward-Roger")
r.squaredGLMM(ldmc_model)


# mort models -------------------------------------------------------------
mort_day50 <- mort_day50 %>%
  mutate(vpd.k = vpd / 10)

mort_model <- glmer(status ~ ot + pt + vpd.k + 
                      ot*vpd.k + pt*vpd.k + ot*pt + ot*pt*vpd.k + (1 | tray),
               data = mort_day50,
               family = binomial(link = "logit"),
               control = glmerControl(optimizer ="bobyqa",
                                      optCtrl=list(maxfun=1e6)),
               na.action = na.exclude)
# assess
check_model(mort_model)
plot(mort_model, add.smooth = FALSE, which = 1)
E <- resid(mort_model)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(mort_model))
qqline(resid(mort_model))              

summary(mort_model)
r.squaredGLMM(mort_model)
Anova(mort_model, type = 3)

# seed num models ---------------------------------------------------------
SEED_FINAL <- SEED_FINAL %>%
  mutate(vpd.k = vpd / 10)

# Fit the model using the vpd.k
seednum_model1 <- glmer.nb(num_total ~ ot + pt + vpd.k + 
                         ot*vpd.k + pt*vpd.k + ot*pt + ot*pt*vpd.k + 
                         (1 | tray),
                       control = glmerControl(optimizer ="bobyqa", optCtrl=list(maxfun=2000000)),
                       family = poisson(),
                       data = SEED_FINAL,
                       na.action = na.exclude)

# assess
check_model(seednum_model1)
plot(seednum_model1, add.smooth = FALSE, which = 1)
E <- resid(seednum_model1)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(seednum_model1))
qqline(resid(seednum_model1))              

summary(seednum_model1)
r.squaredGLMM(seednum_model1)
Anova(seednum_model1, type = 3)


# days to flower model ----------------------------------------------------
filtered_flower_data <- FLOWER_FINAL[!FLOWER_FINAL$pop %in% c(3, 5, 11), ]

FLOWER_FINAL <- FLOWER_FINAL |> 
  mutate(days_to_flower = as.numeric(days_to_flower))

FLOWER_FINAL <- FLOWER_FINAL %>%
  mutate(vpd.k = vpd / 10)

# full model
daysflower_model <- glmer.nb(days_to_flower ~ ot + pt + vpd.k + 
                               ot*vpd.k + pt*vpd.k + ot*pt + ot*pt*vpd.k + 
                               (1 | tray),
                  control = glmerControl(optimizer ="bobyqa", optCtrl=list(maxfun=2000000)),
                  data = FLOWER_FINAL,
                  na.action = na.exclude)
# assess
check_model(daysflower_model)
plot(daysflower_model, add.smooth = FALSE, which = 1)
E <- resid(daysflower_model)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(daysflower_model))
qqline(resid(daysflower_model))              

summary(daysflower_model)
r.squaredGLMM(daysflower_model)
Anova(daysflower_model, type = 3)


# flower status model ----------------------------------------------------

flower_status <- flower_status %>%
  mutate(vpd.k = vpd / 10)

numflow_model <- glmer(status ~ ot + pt + vpd.k + 
                         ot*vpd.k + pt*vpd.k + ot*pt + ot*pt*vpd.k + (1 | tray),
                       data = flower_status,
                       family = binomial(link = "logit"),
                       control = glmerControl(optimizer ="bobyqa",
                                              optCtrl=list(maxfun=1e6)),
                       na.action = na.exclude)

check_model(numflow_model)
plot(numflow_model, add.smooth = FALSE, which = 1)
E <- resid(numflow_model)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(numflow_model))
qqline(resid(numflow_model))
plot(flower_status$pop, resid(numflow_model), xlab = "pop", ylab = "residuals")
plot(flower_status$`2023_treat`, resid(numflow_model), xlab = "2023 treat", ylab = "residuals")

# final stats
summary(numflow_model)
r.squaredGLMM(numflow_model)
Anova(numflow_model, type = 3)


















# OLD ANALYSIS ------------------------------------------------------------
# CV ----------------------------------------------------------------------
# root models -------------------------------------------------------------
cv_mean <- mean(BIOMASS_FINAL$cv, na.rm = TRUE)
BIOMASS_FINAL$cv.c <- BIOMASS_FINAL$cv - cv_mean
BIOMASS_FINAL <- BIOMASS_FINAL |> 
  mutate(cv.c = as.numeric(cv.c))

root_model <- lmer(root ~ ot + pt + cv.c + 
                     ot*cv.c + pt*cv.c + ot*pt + ot*pt*cv.c + (1 | tray),
                   data = BIOMASS_FINAL,
                   na.action = na.exclude)

# assess
vif(root_model)
plot(BIOMASS_FINAL$ot, resid(root_model), xlab = "ot", ylab = "residuals")
plot(BIOMASS_FINAL$tray, resid(root_model), xlab = "tray", ylab = "residuals")
plot(root_model, add.smooth = FALSE, which = 1)
E <- resid(root_model)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(root_model))
qqline(resid(root_model))

summary(root_model)
Anova(root_model, type="III",  test.statistic= "F")
r.squaredGLMM(root_model)

#contrasts
rootcontrast <- emmeans(root_model, specs = pairwise ~ vpd:ot:pt)
rootcontrast$contrasts

# shoot models ------------------------------------------------------------
#same structure as OG model - tray has a crazy amount of variance btwn groups
shoot_model <- lme(shoot ~ ot + pt + cv.c + 
                     ot*cv.c +
                     pt*cv.c +
                     ot*pt +
                     ot*pt*cv.c,
                   random = ~ 1 | tray,
                   weights = varIdent(form = ~ 1 | tray),
                   data = BIOMASS_FINAL,
                   na.action = na.exclude,
                   control = lmeControl(opt = "optim", maxIter = 100, msMaxIter = 100))

# assess
vif(shoot_model)
plot(BIOMASS_FINAL$ot, resid(shoot_model), xlab = "ot", ylab = "residuals")
plot(BIOMASS_FINAL$tray, resid(shoot_model), xlab = "tray", ylab = "residuals")
check_model(shoot_model)
plot(shoot_model, add.smooth = FALSE, which = 1)
E <- resid(shoot_model)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(shoot_model))
qqline(resid(shoot_model))

summary(shoot_model)
Anova(shoot_model, type="III",  test.statistic= "F", ddf="Kenward-Roger")
r.squaredGLMM(shoot_model)

#contrasts
shootcontrast <- emmeans(shoot_model2, specs = pairwise ~ ot:pt)
shootcontrast$contrasts

# total models ------------------------------------------------------------
total_model <- lmer(total_biomass ~ ot + pt + cv.c + 
                      ot*cv.c + pt*cv.c + ot*pt + ot*pt*cv.c + (1 | tray),
                    data = BIOMASS_FINAL,
                    na.action = na.exclude)

summary(total_model)
vif(total_model)

# assess
plot(BIOMASS_FINAL$ot, resid(total_model), xlab = "ot", ylab = "residuals")
plot(BIOMASS_FINAL$tray, resid(total_model), xlab = "tray", ylab = "residuals")
plot(total_model, add.smooth = FALSE, which = 1)
E <- resid(total_model)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(total_model))
qqline(resid(total_model))              

summary(total_model)
Anova(total_model, type="III",  test.statistic= "F", ddf="Kenward-Roger")
r.squaredGLMM(total_model)

#contrasts
totalcontrast <- emmeans(total_model, specs = pairwise ~ ot:pt)
totalcontrast$contrasts

# RS models -------------------------------------------------------------
rs_model <- lmer(log_RS ~ ot + pt + cv.c + 
                   ot*cv.c + pt*cv.c + ot*pt + ot*pt*cv.c + (1 | tray),
                 data = BIOMASS_FINAL,
                 na.action = na.exclude)

summary(rs_model)
vif(rs_model)

# assess
plot(BIOMASS_FINAL$ot, resid(rs_model), xlab = "ot", ylab = "residuals")
plot(BIOMASS_FINAL$tray, resid(rs_model), xlab = "tray", ylab = "residuals")
plot(rs_model, add.smooth = FALSE, which = 1)
E <- resid(rs_model)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(rs_model))
qqline(resid(rs_model))              

summary(rs_model)
Anova(rs_model, type="III",  test.statistic= "F", ddf="Kenward-Roger")
r.squaredGLMM(rs_model)

#contrasts
rscontrast <- emmeans(rs_model, specs = pairwise ~ ot:pt)
rscontrast$contrasts

# height models -------------------------------------------------------------
cv_mean <- mean(HEIGHT_FINAL$cv, na.rm = TRUE)
HEIGHT_FINAL$cv.c <- HEIGHT_FINAL$cv - cv_mean
HEIGHT_FINAL <- HEIGHT_FINAL |> 
  mutate(cv.c = as.numeric(cv.c))

height_model <- lme(max ~ ot + pt + cv.c + 
                      ot*cv.c +
                      pt*cv.c +
                      ot*pt +
                      ot*pt*cv.c,
                    random = ~ 1 | tray,
                    weights = varIdent(form = ~ 1 | tray),
                    data = HEIGHT_FINAL,
                    control = lmeControl(opt = "optim", maxIter = 100, msMaxIter = 100),
                    na.action = na.exclude)

# assess
plot(HEIGHT_FINAL$ot, resid(height_model), xlab = "ot", ylab = "residuals")
plot(HEIGHT_FINAL$tray, resid(height_model), xlab = "tray", ylab = "residuals")
plot(height_model, add.smooth = FALSE, which = 1)
E <- resid(height_model)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(height_model))
qqline(resid(height_model))              

summary(height_model)
Anova(height_model, type="III",  test.statistic= "F", ddf="Kenward-Roger")
r.squaredGLMM(height_model)

#contrasts
heightcontrast <- emmeans(height_model, specs = pairwise ~ ot:pt)
heightcontrast$contrasts

# rgr models -------------------------------------------------------------
cv_mean <- mean(RGR_FINAL$cv, na.rm = TRUE)
RGR_FINAL$cv.c <- RGR_FINAL$cv - cv_mean
RGR_FINAL <- RGR_FINAL |> 
  mutate(cv.c = as.numeric(cv.c))

rgr_model <- lme(RGR ~ ot + pt + cv.c + 
                   ot*cv.c +
                   pt*cv.c +
                   ot*pt +
                   ot*pt*cv.c,
                 random = ~ 1 | tray,
                 weights = varIdent(form = ~ 1 | tray),
                 data = RGR_FINAL,
                 control = lmeControl(opt = "optim", maxIter = 100, msMaxIter = 100),
                 na.action = na.exclude)

# assess
plot(RGR_FINAL$ot, resid(rgr_model), xlab = "ot", ylab = "residuals")
plot(RGR_FINAL$tray, resid(rgr_model), xlab = "tray", ylab = "residuals")
plot(rgr_model, add.smooth = FALSE, which = 1)
E <- resid(rgr_model)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(rgr_model))
qqline(resid(rgr_model))              

summary(rgr_model)
Anova(rgr_model, type="III",  test.statistic= "F", ddf="Kenward-Roger")
r.squaredGLMM(rgr_model)

#contrasts
rgrcontrast <- emmeans(rgr_model, specs = pairwise ~ ot:pt)
rgrcontrast$contrasts

# sla models --------------------------------------------------------------
cv_mean <- mean(SLA_LDMC_FINAL$cv, na.rm = TRUE)
SLA_LDMC_FINAL$cv.c <- SLA_LDMC_FINAL$cv - cv_mean
SLA_LDMC_FINAL <- SLA_LDMC_FINAL |> 
  mutate(cv.c = as.numeric(cv.c))

sla_model <- lmer(log(sla) ~ ot + pt + cv.c + 
                    ot*cv.c + pt*cv.c + ot*pt + ot*pt*cv.c + (1 | tray),
                  data = SLA_LDMC_FINAL,
                  na.action = na.exclude)

# assess
plot(SLA_LDMC_FINAL$ot, resid(sla_model), xlab = "ot", ylab = "residuals")
plot(SLA_LDMC_FINAL$tray, resid(sla_model), xlab = "tray", ylab = "residuals")
plot(sla_model, add.smooth = FALSE, which = 1)
E <- resid(sla_model)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(sla_model))
qqline(resid(sla_model))              

summary(sla_model)
Anova(sla_model, type="III",  test.statistic= "F", ddf="Kenward-Roger")
r.squaredGLMM(sla_model)

# ldmc models --------------------------------------------------------------
ldmc_model <- lmer(log(ldmc) ~ ot + pt + cv.c + 
                     ot*cv.c + pt*cv.c + ot*pt + ot*pt*cv.c + (1 | tray),
                   data = SLA_LDMC_FINAL,
                   na.action = na.exclude)

# assess
plot(SLA_LDMC_FINAL$ot, resid(ldmc_model), xlab = "ot", ylab = "residuals")
plot(SLA_LDMC_FINAL$tray, resid(ldmc_model), xlab = "tray", ylab = "residuals")
plot(ldmc_model, add.smooth = FALSE, which = 1)
E <- resid(ldmc_model)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(ldmc_model))
qqline(resid(ldmc_model))              

summary(ldmc_model)
Anova(ldmc_model, type="III",  test.statistic= "F", ddf="Kenward-Roger")
r.squaredGLMM(ldmc_model)


# mort models -------------------------------------------------------------
cv_mean <- mean(mort_day50$cv, na.rm = TRUE)
mort_day50$cv.c <- mort_day50$cv - cv_mean
mort_day50 <- mort_day50 |> 
  mutate(cv.c = as.numeric(cv.c))

mort_model <- glmer(status ~ ot + pt + cv.c + 
                      ot*cv.c + pt*cv.c + ot*pt + ot*pt*cv.c + (1 | tray),
                    data = mort_day50,
                    family = binomial(link = "logit"),
                    control = glmerControl(optimizer ="bobyqa",
                                           optCtrl=list(maxfun=1e6)),
                    na.action = na.exclude)
# assess
check_model(mort_model)
plot(mort_model, add.smooth = FALSE, which = 1)
E <- resid(mort_model)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(mort_model))
qqline(resid(mort_model))              

summary(mort_model)
r.squaredGLMM(mort_model)
Anova(mort_model, type = 3)

# seed num models ---------------------------------------------------------
filtered_seed_data <- SEED_FINAL[!SEED_FINAL$pop %in% c(3, 5, 11), ]

cv_mean <- mean(filtered_seed_data$cv, na.rm = TRUE)
filtered_seed_data$cv.c <- filtered_seed_data$cv - cv_mean
filtered_seed_data <- filtered_seed_data |> 
  mutate(cv.c = as.numeric(cv.c))

# Fit the model using the centered vpd
seednum_model <- glmer(num_total ~ ot + pt + cv.c + 
                         ot*cv.c + pt*cv.c + ot*pt + ot*pt*cv.c + 
                         (1 | tray),
                       control = glmerControl(optimizer ="bobyqa", optCtrl=list(maxfun=2000000)),
                       family = poisson(),
                       data = filtered_seed_data,
                       na.action = na.exclude)

# assess
check_model(seednum_model)
plot(seednum_model, add.smooth = FALSE, which = 1)
E <- resid(seednum_model)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(seednum_model))
qqline(resid(seednum_model))              

summary(seednum_model)
r.squaredGLMM(seednum_model)
Anova(seednum_model, type = 3)


# days to flower model ----------------------------------------------------
filtered_flower_data <- FLOWER_FINAL[!FLOWER_FINAL$pop %in% c(3, 5, 11), ]

filtered_flower_data <- filtered_flower_data |> 
  mutate(days_to_flower = as.numeric(days_to_flower))

cv_mean <- mean(filtered_flower_data$cv, na.rm = TRUE)
filtered_flower_data$cv.c <- filtered_flower_data$cv - cv_mean
filtered_flower_data <- filtered_flower_data |> 
  mutate(cv.c = as.numeric(cv.c))

# full model
daysflower_model <- glmer.nb(days_to_flower ~ ot + pt + cv.c + 
                               ot*cv.c + pt*cv.c + ot*pt + ot*pt*cv.c + 
                               (1 | tray),
                             control = glmerControl(optimizer ="bobyqa", optCtrl=list(maxfun=2000000)),
                             data = filtered_flower_data,
                             na.action = na.exclude)
# assess
check_model(daysflower_model)
plot(daysflower_model, add.smooth = FALSE, which = 1)
E <- resid(daysflower_model)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(daysflower_model))
qqline(resid(daysflower_model))              

summary(daysflower_model)
r.squaredGLMM(daysflower_model)
Anova(daysflower_model, type = 3)





























# OLD ANALYSIS ------------------------------------------------------------
# trait ~ OT + PT + SAT_C + SAP_mm + cv +
#         OT*PT +
#         OT*SAT_C + OT*SAP_mm + OT*cv +
#         PT*SAT_C + PT*SAP_mm + PT*cv

# remember you might have to scale SAP, SAT, cv!

# root models -------------------------------------------------------------
# full model
model <- lmer(root ~ ot + pt + SAT_C + SAP_mm + cv +
                    ot*pt + SAT_C*SAP_mm
                    ot*SAT_C + ot*SAP_mm + ot*cv +
                    pt*SAT_C + pt*SAP_mm + pt*cv +
                    (1|tray),
              data = BIOMASS_FINAL,
              na.action = na.exclude)
summary(model)
# assess
plot(model, add.smooth = FALSE, which = 1)
E <- resid(model)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(model))
qqline(resid(model))              

summary(model)
Anova(model, type="III",  test.statistic= "F", ddf="Kenward-Roger")
r.squaredGLMM(model)

#AIC
model1 <- lmer(root ~ ot + pt + SAT_C + SAP_mm + cv +
                ot*pt +
                ot*SAT_C + ot*SAP_mm + ot*cv +
                pt*SAT_C + pt*SAP_mm + cv +
                (1|tray),
              data = BIOMASS_FINAL,
              na.action = na.exclude)
summary(model1)

model2 <- lmer(root ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C + ot*SAP_mm + ot*cv +
                 pt*SAT_C + pt*SAP_mm +
                 (1|tray),
               data = BIOMASS_FINAL,
               na.action = na.exclude)
summary(model2)

model3 <- lmer(root ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C + ot*SAP_mm + ot*cv +
                 pt*SAT_C +
                 (1|tray),
               data = BIOMASS_FINAL,
               na.action = na.exclude)
summary(model3)

model4 <- lmer(root ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C + ot*SAP_mm + ot*cv +
                 (1|tray),
               data = BIOMASS_FINAL,
               na.action = na.exclude)
summary(model4)

model5 <- lmer(root ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C + ot*SAP_mm +
                 (1|tray),
               data = BIOMASS_FINAL,
               na.action = na.exclude)
summary(model5)

model6 <- lmer(root ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C +
                 (1|tray),
               data = BIOMASS_FINAL,
               na.action = na.exclude)
summary(model6)

model7 <- lmer(root ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 (1|tray),
               data = BIOMASS_FINAL,
               na.action = na.exclude)
summary(model7)

model8 <- lmer(root ~ ot + pt + SAT_C + SAP_mm + cv +
                 (1|tray),
               data = BIOMASS_FINAL,
               na.action = na.exclude)
summary(model8)

AIC(model1, model2, model3, model4, model5, model6, model7, model8)

plot(model8, add.smooth = FALSE, which = 1)
E <- resid(model8)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(model8))
qqline(resid(model8))              

summary(model8)
Anova(model8, type="III",  test.statistic= "F", ddf="Kenward-Roger")
r.squaredGLMM(model8)


# shoot models ------------------------------------------------------------
# full model
model <- lmer(log(shoot) ~ ot + pt + SAT_C + SAP_mm + cv +
                ot*pt +
                ot*SAT_C + ot*SAP_mm + ot*cv +
                pt*SAT_C + pt*SAP_mm + pt*cv +
                (1|tray),
              data = BIOMASS_FINAL,
              na.action = na.exclude)
# assess
check_model(model)
plot(model, add.smooth = FALSE, which = 1)
E <- resid(model)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(model))
qqline(resid(model))              

summary(model)
Anova(model, type="III",  test.statistic= "F", ddf="Kenward-Roger")
r.squaredGLMM(model)

#AIC
model1 <- lmer(log(shoot) ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C + ot*SAP_mm + ot*cv +
                 pt*SAT_C + pt*SAP_mm + cv +
                 (1|tray),
               data = BIOMASS_FINAL,
               na.action = na.exclude)
summary(model1)

model2 <- lmer(log(shoot) ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C + ot*SAP_mm + ot*cv +
                 pt*SAT_C + pt*SAP_mm +
                 (1|tray),
               data = BIOMASS_FINAL,
               na.action = na.exclude)
summary(model2)

model3 <- lmer(log(shoot) ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C + ot*SAP_mm + ot*cv +
                 pt*SAT_C +
                 (1|tray),
               data = BIOMASS_FINAL,
               na.action = na.exclude)
summary(model3)

model4 <- lmer(log(shoot) ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C + ot*SAP_mm + ot*cv +
                 (1|tray),
               data = BIOMASS_FINAL,
               na.action = na.exclude)
summary(model4)

model5 <- lmer(log(shoot) ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C + ot*SAP_mm +
                 (1|tray),
               data = BIOMASS_FINAL,
               na.action = na.exclude)
summary(model5)

model6 <- lmer(log(shoot) ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C +
                 (1|tray),
               data = BIOMASS_FINAL,
               na.action = na.exclude)
summary(model6)

model7 <- lmer(log(shoot) ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 (1|tray),
               data = BIOMASS_FINAL,
               na.action = na.exclude)
summary(model7)

model8 <- lmer(log(shoot) ~ ot + pt + SAT_C + SAP_mm + cv +
                 (1|tray),
               data = BIOMASS_FINAL,
               na.action = na.exclude)
summary(model8)

AIC(model1, model2, model3, model4, model5, model6, model7, model8)

plot(model8, add.smooth = FALSE, which = 1)
E <- resid(model8)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(model8))
qqline(resid(model8))              

summary(model8)
Anova(model8, type="III",  test.statistic= "F", ddf="Kenward-Roger")
r.squaredGLMM(model8)


# total models ------------------------------------------------------------
# full model
model <- lmer(total_biomass ~ ot + pt + SAT_C + SAP_mm + cv +
                ot*pt +
                ot*SAT_C + ot*SAP_mm + ot*cv +
                pt*SAT_C + pt*SAP_mm + pt*cv +
                (1|tray),
              data = BIOMASS_FINAL,
              na.action = na.exclude)
# assess
check_model(model)
plot(model, add.smooth = FALSE, which = 1)
E <- resid(model)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(model))
qqline(resid(model))              

summary(model)
Anova(model, type="III",  test.statistic= "F", ddf="Kenward-Roger")
r.squaredGLMM(model)

#AIC
model1 <- lmer(total_biomass ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C + ot*SAP_mm + ot*cv +
                 pt*SAT_C + pt*SAP_mm + cv +
                 (1|tray),
               data = BIOMASS_FINAL,
               na.action = na.exclude)
summary(model1)

model2 <- lmer(total_biomass ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C + ot*SAP_mm + ot*cv +
                 pt*SAT_C + pt*SAP_mm +
                 (1|tray),
               data = BIOMASS_FINAL,
               na.action = na.exclude)
summary(model2)

model3 <- lmer(total_biomass ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C + ot*SAP_mm + ot*cv +
                 pt*SAT_C +
                 (1|tray),
               data = BIOMASS_FINAL,
               na.action = na.exclude)
summary(model3)

model4 <- lmer(total_biomass ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C + ot*SAP_mm + ot*cv +
                 (1|tray),
               data = BIOMASS_FINAL,
               na.action = na.exclude)
summary(model4)

model5 <- lmer(total_biomass ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C + ot*SAP_mm +
                 (1|tray),
               data = BIOMASS_FINAL,
               na.action = na.exclude)
summary(model5)

model6 <- lmer(total_biomass ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C +
                 (1|tray),
               data = BIOMASS_FINAL,
               na.action = na.exclude)
summary(model6)

model7 <- lmer(total_biomass ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 (1|tray),
               data = BIOMASS_FINAL,
               na.action = na.exclude)
summary(model7)

model8 <- lmer(total_biomass ~ ot + pt + SAT_C + SAP_mm + cv +
                 (1|tray),
               data = BIOMASS_FINAL,
               na.action = na.exclude)
summary(model8)

AIC(model1, model2, model3, model4, model5, model6, model7, model8)

plot(model8, add.smooth = FALSE, which = 1)
E <- resid(model8)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(model8))
qqline(resid(model8))              

summary(model8)
Anova(model8, type="III",  test.statistic= "F", ddf="Kenward-Roger")
r.squaredGLMM(model8) 


# RS models ---------------------------------------------------------------
# full model
model <- lmer(log_RS ~ ot + pt + SAT_C + SAP_mm + cv +
                ot*pt +
                ot*SAT_C + ot*SAP_mm + ot*cv +
                pt*SAT_C + pt*SAP_mm + pt*cv +
                (1|tray),
              data = BIOMASS_FINAL,
              na.action = na.exclude)
# assess
check_model(model)
plot(model, add.smooth = FALSE, which = 1)
E <- resid(model)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(model))
qqline(resid(model))              

summary(model)
Anova(model, type="III",  test.statistic= "F", ddf="Kenward-Roger")
r.squaredGLMM(model)

#AIC
model1 <- lmer(log_RS ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C + ot*SAP_mm + ot*cv +
                 pt*SAT_C + pt*SAP_mm + cv +
                 (1|tray),
               data = BIOMASS_FINAL,
               na.action = na.exclude)
summary(model1)

model2 <- lmer(log_RS ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C + ot*SAP_mm + ot*cv +
                 pt*SAT_C + pt*SAP_mm +
                 (1|tray),
               data = BIOMASS_FINAL,
               na.action = na.exclude)
summary(model2)

model3 <- lmer(log_RS ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C + ot*SAP_mm + ot*cv +
                 pt*SAT_C +
                 (1|tray),
               data = BIOMASS_FINAL,
               na.action = na.exclude)
summary(model3)

model4 <- lmer(log_RS ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C + ot*SAP_mm + ot*cv +
                 (1|tray),
               data = BIOMASS_FINAL,
               na.action = na.exclude)
summary(model4)

model5 <- lmer(log_RS ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C + ot*SAP_mm +
                 (1|tray),
               data = BIOMASS_FINAL,
               na.action = na.exclude)
summary(model5)

model6 <- lmer(log_RS ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C +
                 (1|tray),
               data = BIOMASS_FINAL,
               na.action = na.exclude)
summary(model6)

model7 <- lmer(log_RS ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 (1|tray),
               data = BIOMASS_FINAL,
               na.action = na.exclude)
summary(model7)

model8 <- lmer(log_RS ~ ot + pt + SAT_C + SAP_mm + cv +
                 (1|tray),
               data = BIOMASS_FINAL,
               na.action = na.exclude)
summary(model8)

AIC(model1, model2, model3, model4, model5, model6, model7, model8)

plot(model8, add.smooth = FALSE, which = 1)
E <- resid(model8)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(model8))
qqline(resid(model8))              

summary(model8)
Anova(model8, type="III",  test.statistic= "F", ddf="Kenward-Roger")
r.squaredGLMM(model8)


# height models -----------------------------------------------------------
# full model
model <- lmer(max ~ ot + pt + SAT_C + SAP_mm + cv +
                ot*pt +
                ot*SAT_C + ot*SAP_mm + ot*cv +
                pt*SAT_C + pt*SAP_mm + pt*cv +
                (1|tray),
              data = HEIGHT_FINAL,
              na.action = na.exclude)
# assess
check_model(model)
plot(model, add.smooth = FALSE, which = 1)
E <- resid(model)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(model))
qqline(resid(model))              

summary(model)
Anova(model, type="III",  test.statistic= "F", ddf="Kenward-Roger")
r.squaredGLMM(model)

#AIC
model1 <- lmer(max ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C + ot*SAP_mm + ot*cv +
                 pt*SAT_C + pt*SAP_mm + cv +
                 (1|tray),
               data = HEIGHT_FINAL,
               na.action = na.exclude)
summary(model1)

model2 <- lmer(max ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C + ot*SAP_mm + ot*cv +
                 pt*SAT_C + pt*SAP_mm +
                 (1|tray),
               data = HEIGHT_FINAL,
               na.action = na.exclude)
summary(model2)

model3 <- lmer(max ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C + ot*SAP_mm + ot*cv +
                 pt*SAT_C +
                 (1|tray),
               data = HEIGHT_FINAL,
               na.action = na.exclude)
summary(model3)

model4 <- lmer(max ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C + ot*SAP_mm + ot*cv +
                 (1|tray),
               data = HEIGHT_FINAL,
               na.action = na.exclude)
summary(model4)

model5 <- lmer(max ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C + ot*SAP_mm +
                 (1|tray),
               data = HEIGHT_FINAL,
               na.action = na.exclude)
summary(model5)

model6 <- lmer(max ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C +
                 (1|tray),
               data = HEIGHT_FINAL,
               na.action = na.exclude)
summary(model6)

model7 <- lmer(max ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 (1|tray),
               data = HEIGHT_FINAL,
               na.action = na.exclude)
summary(model7)

model8 <- lmer(max ~ ot + pt + SAT_C + SAP_mm + cv +
                 (1|tray),
               data = HEIGHT_FINAL,
               na.action = na.exclude)
summary(model8)

AIC(model1, model2, model3, model4, model5, model6, model7, model8)

plot(model8, add.smooth = FALSE, which = 1)
E <- resid(model8)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(model8))
qqline(resid(model8))              

summary(model8)
Anova(model8, type="III",  test.statistic= "F", ddf="Kenward-Roger")
r.squaredGLMM(model8)


# rgr models --------------------------------------------------------------
# full model
model <- lmer(RGR ~ ot + pt + SAT_C + SAP_mm + cv +
                ot*pt +
                ot*SAT_C + ot*SAP_mm + ot*cv +
                pt*SAT_C + pt*SAP_mm + pt*cv +
                (1|tray),
              data = RGR_FINAL,
              na.action = na.exclude)
# assess
check_model(model)
plot(model, add.smooth = FALSE, which = 1)
E <- resid(model)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(model))
qqline(resid(model))              

summary(model)
Anova(model, type="III",  test.statistic= "F", ddf="Kenward-Roger")
r.squaredGLMM(model)

#AIC
model1 <- lmer(RGR ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C + ot*SAP_mm + ot*cv +
                 pt*SAT_C + pt*SAP_mm + cv +
                 (1|tray),
               data = RGR_FINAL,
               na.action = na.exclude)
summary(model1)

model2 <- lmer(RGR ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C + ot*SAP_mm + ot*cv +
                 pt*SAT_C + pt*SAP_mm +
                 (1|tray),
               data = RGR_FINAL,
               na.action = na.exclude)
summary(model2)

model3 <- lmer(RGR ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C + ot*SAP_mm + ot*cv +
                 pt*SAT_C +
                 (1|tray),
               data = RGR_FINAL,
               na.action = na.exclude)
summary(model3)

model4 <- lmer(RGR ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C + ot*SAP_mm + ot*cv +
                 (1|tray),
               data = RGR_FINAL,
               na.action = na.exclude)
summary(model4)

model5 <- lmer(RGR ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C + ot*SAP_mm +
                 (1|tray),
               data = RGR_FINAL,
               na.action = na.exclude)
summary(model5)

model6 <- lmer(RGR ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C +
                 (1|tray),
               data = RGR_FINAL,
               na.action = na.exclude)
summary(model6)

model7 <- lmer(RGR ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 (1|tray),
               data = RGR_FINAL,
               na.action = na.exclude)
summary(model7)

model8 <- lmer(RGR ~ ot + pt + SAT_C + SAP_mm + cv +
                 (1|tray),
               data = RGR_FINAL,
               na.action = na.exclude)
summary(model8)

AIC(model1, model2, model3, model4, model5, model6, model7, model8)

plot(model8, add.smooth = FALSE, which = 1)
E <- resid(model8)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(model8))
qqline(resid(model8))              

summary(model8)
Anova(model8, type="III",  test.statistic= "F", ddf="Kenward-Roger")
r.squaredGLMM(model8)

# sla models --------------------------------------------------------------
# full model
model <- lmer(log(sla) ~ ot + pt + SAT_C + SAP_mm + cv +
                ot*pt +
                ot*SAT_C + ot*SAP_mm + ot*cv +
                pt*SAT_C + pt*SAP_mm + pt*cv +
                (1|tray),
              data = SLA_LDMC_FINAL,
              na.action = na.exclude)
# assess
check_model(model)
plot(model, add.smooth = FALSE, which = 1)
E <- resid(model)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(model))
qqline(resid(model))              

summary(model)
Anova(model, type="III",  test.statistic= "F", ddf="Kenward-Roger")
r.squaredGLMM(model)

#AIC
model1 <- lmer(log(sla) ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C + ot*SAP_mm + ot*cv +
                 pt*SAT_C + pt*SAP_mm + cv +
                 (1|tray),
               data = SLA_LDMC_FINAL,
               na.action = na.exclude)
summary(model1)

model2 <- lmer(log(sla) ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C + ot*SAP_mm + ot*cv +
                 pt*SAT_C + pt*SAP_mm +
                 (1|tray),
               data = SLA_LDMC_FINAL,
               na.action = na.exclude)
summary(model2)

model3 <- lmer(log(sla) ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C + ot*SAP_mm + ot*cv +
                 pt*SAT_C +
                 (1|tray),
               data = SLA_LDMC_FINAL,
               na.action = na.exclude)
summary(model3)

model4 <- lmer(log(sla) ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C + ot*SAP_mm + ot*cv +
                 (1|tray),
               data = SLA_LDMC_FINAL,
               na.action = na.exclude)
summary(model4)

model5 <- lmer(log(sla) ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C + ot*SAP_mm +
                 (1|tray),
               data = SLA_LDMC_FINAL,
               na.action = na.exclude)
summary(model5)

model6 <- lmer(log(sla) ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C +
                 (1|tray),
               data = SLA_LDMC_FINAL,
               na.action = na.exclude)
summary(model6)

model7 <- lmer(log(sla) ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 (1|tray),
               data = SLA_LDMC_FINAL,
               na.action = na.exclude)
summary(model7)

model8 <- lmer(log(sla) ~ ot + pt + SAT_C + SAP_mm + cv +
                 (1|tray),
               data = SLA_LDMC_FINAL,
               na.action = na.exclude)
summary(model8)

AIC(model1, model2, model3, model4, model5, model6, model7, model8)

plot(model8, add.smooth = FALSE, which = 1)
E <- resid(model8)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(model8))
qqline(resid(model8))              

summary(model8)
Anova(model8, type="III",  test.statistic= "F", ddf="Kenward-Roger")
r.squaredGLMM(model8)

# ldmc models -------------------------------------------------------------
# full model
model <- lmer(log(ldmc) ~ ot + pt + SAT_C + SAP_mm + cv +
                ot*pt +
                ot*SAT_C + ot*SAP_mm + ot*cv +
                pt*SAT_C + pt*SAP_mm + pt*cv +
                (1|tray),
              data = SLA_LDMC_FINAL,
              na.action = na.exclude)
# assess
check_model(model)
plot(model, add.smooth = FALSE, which = 1)
E <- resid(model)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(model))
qqline(resid(model))              

summary(model)
Anova(model, type="III",  test.statistic= "F", ddf="Kenward-Roger")
r.squaredGLMM(model)

#AIC
model1 <- lmer(log(ldmc) ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C + ot*SAP_mm + ot*cv +
                 pt*SAT_C + pt*SAP_mm + cv +
                 (1|tray),
               data = SLA_LDMC_FINAL,
               na.action = na.exclude)
summary(model1)

model2 <- lmer(log(ldmc) ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C + ot*SAP_mm + ot*cv +
                 pt*SAT_C + pt*SAP_mm +
                 (1|tray),
               data = SLA_LDMC_FINAL,
               na.action = na.exclude)
summary(model2)

model3 <- lmer(log(ldmc) ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C + ot*SAP_mm + ot*cv +
                 pt*SAT_C +
                 (1|tray),
               data = SLA_LDMC_FINAL,
               na.action = na.exclude)
summary(model3)

model4 <- lmer(log(ldmc) ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C + ot*SAP_mm + ot*cv +
                 (1|tray),
               data = SLA_LDMC_FINAL,
               na.action = na.exclude)
summary(model4)

model5 <- lmer(log(ldmc) ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C + ot*SAP_mm +
                 (1|tray),
               data = SLA_LDMC_FINAL,
               na.action = na.exclude)
summary(model5)

model6 <- lmer(log(ldmc) ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C +
                 (1|tray),
               data = SLA_LDMC_FINAL,
               na.action = na.exclude)
summary(model6)

model7 <- lmer(log(ldmc) ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 (1|tray),
               data = SLA_LDMC_FINAL,
               na.action = na.exclude)
summary(model7)

model8 <- lmer(log(ldmc) ~ ot + pt + SAT_C + SAP_mm + cv +
                 (1|tray),
               data = SLA_LDMC_FINAL,
               na.action = na.exclude)
summary(model8)

AIC(model1, model2, model3, model4, model5, model6, model7, model8)

plot(model8, add.smooth = FALSE, which = 1)
E <- resid(model8)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(model8))
qqline(resid(model8))              

summary(model8)
Anova(model8, type="III",  test.statistic= "F", ddf="Kenward-Roger")
r.squaredGLMM(model8)

# mort models -------------------------------------------------------------
# full model
scaled <- transform(mort_day50,
                    SAP_mm=scale(SAP_mm),
                    SAT_C=scale(SAT_C),
                    cv = scale(cv))

model <- glmer(status ~ ot + pt + SAT_C + SAP_mm + cv +
                ot*pt +
                ot*SAT_C + ot*SAP_mm + ot*cv +
                pt*SAT_C + pt*SAP_mm + pt*cv +
                (1|tray),
              data = scaled,
              family = binomial(link = "logit"),
              na.action = na.exclude)
# assess
check_model(model)
plot(model, add.smooth = FALSE, which = 1)
E <- resid(model)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(model))
qqline(resid(model))              

summary(model)
r.squaredGLMM(model)
Anova(model, type = 3)

#AIC
model1 <- glmer(status ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C + ot*SAP_mm + ot*cv +
                 pt*SAT_C + pt*SAP_mm + cv +
                 (1|tray),
               data = scaled,
               family = binomial(link = "logit"),
               control = glmerControl(optimizer ="bobyqa", optCtrl=list(maxfun=2000000)),
               na.action = na.exclude)
summary(model1)

model2 <- glmer(status ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C + ot*SAP_mm + ot*cv +
                 pt*SAT_C + pt*SAP_mm +
                 (1|tray),
               data = scaled,
               family = binomial(link = "logit"),
               control = glmerControl(optimizer ="bobyqa", optCtrl=list(maxfun=2000000)),
               na.action = na.exclude)
summary(model2)

model3 <- glmer(status ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C + ot*SAP_mm + ot*cv +
                 pt*SAT_C +
                 (1|tray),
                data = scaled,
                family = binomial(link = "logit"),
                control = glmerControl(optimizer ="bobyqa", optCtrl=list(maxfun=2000000)),
                na.action = na.exclude)
summary(model3)

model4 <- glmer(status ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C + ot*SAP_mm + ot*cv +
                 (1|tray),
                data = scaled,
                family = binomial(link = "logit"),
                control = glmerControl(optimizer ="bobyqa", optCtrl=list(maxfun=2000000)),
                na.action = na.exclude)
summary(model4)

model5 <- glmer(status ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C + ot*SAP_mm +
                 (1|tray),
               data = scaled,
               family = binomial(link = "logit"),
               control = glmerControl(optimizer ="bobyqa", optCtrl=list(maxfun=2000000)),
               na.action = na.exclude)
summary(model5)

model6 <- glmer(status ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C +
                 (1|tray),
                data = scaled,
                family = binomial(link = "logit"),
                control = glmerControl(optimizer ="bobyqa", optCtrl=list(maxfun=2000000)),
                na.action = na.exclude)
summary(model6)

model7 <- glmer(status ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 (1|tray),
                data = scaled,
                family = binomial(link = "logit"),
                control = glmerControl(optimizer ="bobyqa", optCtrl=list(maxfun=2000000)),
                na.action = na.exclude)
summary(model7)

model8 <- glmer(status ~ ot + pt + SAT_C + SAP_mm + cv +
                 (1|tray),
                data = scaled,
                family = binomial(link = "logit"),
                control = glmerControl(optimizer ="bobyqa", optCtrl=list(maxfun=2000000)),
                na.action = na.exclude)
summary(model8)

AIC(model1, model2, model3, model4, model5, model6, model7, model8)

plot(model4, add.smooth = FALSE, which = 1)
E <- resid(model4)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(model4))
qqline(resid(model4))              

summary(model4)
r.squaredGLMM(model4)
Anova(model4, type = 3)

# seed num models ---------------------------------------------------------
filtered_seed_data <- SEED_FINAL[!SEED_FINAL$pop %in% c(3, 5, 11), ]

scaled2 <- transform(filtered_seed_data,
                    SAP_mm=scale(SAP_mm),
                    SAT_C=scale(SAT_C),
                    cv = scale(cv))
# full model
model <- glmer.nb(num_total ~ ot + pt + SAT_C + SAP_mm + cv +
                ot*pt +
                ot*SAT_C + ot*SAP_mm + ot*cv +
                pt*SAT_C + pt*SAP_mm + pt*cv +
                (1|tray),
                control = glmerControl(optimizer ="bobyqa", optCtrl=list(maxfun=2000000)),
                data = scaled2,
                na.action = na.exclude)
# assess
check_model(model)
plot(model, add.smooth = FALSE, which = 1)
E <- resid(model)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(model))
qqline(resid(model))              

summary(model)
r.squaredGLMM(model)
Anova(model, type = 3)

#AIC
model1 <- glmer.nb(num_total ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C + ot*SAP_mm + ot*cv +
                 pt*SAT_C + pt*SAP_mm + cv +
                 (1|tray),
                 control = glmerControl(optimizer ="bobyqa", optCtrl=list(maxfun=2000000)),
                 data = scaled2,
                na.action = na.exclude)
summary(model1)

model2 <- glmer.nb(num_total ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C + ot*SAP_mm + ot*cv +
                 pt*SAT_C + pt*SAP_mm +
                 (1|tray),
                 control = glmerControl(optimizer ="bobyqa", optCtrl=list(maxfun=2000000)),
                 data = scaled2,
               na.action = na.exclude)
summary(model2)

model3 <- glmer.nb(num_total ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C + ot*SAP_mm + ot*cv +
                 pt*SAT_C +
                 (1|tray),
                 control = glmerControl(optimizer ="bobyqa", optCtrl=list(maxfun=2000000)),
                 data = scaled2,
               na.action = na.exclude)
summary(model3)

model4 <- glmer.nb(num_total ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C + ot*SAP_mm + ot*cv +
                 (1|tray),
                 control = glmerControl(optimizer ="bobyqa", optCtrl=list(maxfun=2000000)),
                 data = scaled2,
               na.action = na.exclude)
summary(model4)

model5 <- glmer.nb(num_total ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C + ot*SAP_mm +
                 (1|tray),
                 control = glmerControl(optimizer ="bobyqa", optCtrl=list(maxfun=2000000)),
                 data = scaled2,
               na.action = na.exclude)
summary(model5)

model6 <- glmer.nb(num_total ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 ot*SAT_C +
                 (1|tray),
                 control = glmerControl(optimizer ="bobyqa", optCtrl=list(maxfun=2000000)),
                 data = scaled2,
               na.action = na.exclude)
summary(model6)

model7 <- glmer.nb(num_total ~ ot + pt + SAT_C + SAP_mm + cv +
                 ot*pt +
                 (1|tray),
                 control = glmerControl(optimizer ="bobyqa", optCtrl=list(maxfun=2000000)),
                 data = scaled2,
               na.action = na.exclude)
summary(model7)

model8 <- glmer.nb(num_total ~ ot + pt + SAT_C + SAP_mm + cv +
                 (1|tray),
                 control = glmerControl(optimizer ="bobyqa", optCtrl=list(maxfun=2000000)),
                 data = scaled2,
               na.action = na.exclude)
summary(model8)

AIC(model1, model2, model3, model4, model5, model6, model7, model8)

plot(model7, add.smooth = FALSE, which = 1)
E <- resid(model7)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(model7))
qqline(resid(model7))              

summary(model7)
r.squaredGLMM(model7)
Anova(model7, type = 3)

# days to flower models ---------------------------------------------------
filtered_flower_data <- FLOWER_FINAL[!FLOWER_FINAL$pop %in% c(3, 5, 11), ]

filtered_flower_data <- filtered_flower_data |> 
  mutate(days_to_flower = as.numeric(days_to_flower))

scaled3 <- transform(filtered_flower_data,
                     SAP_mm=scale(SAP_mm),
                     SAT_C=scale(SAT_C),
                     cv = scale(cv))

# full model
model <- glmer.nb(days_to_flower ~ ot + pt + SAT_C + SAP_mm + cv +
                    ot*pt +
                    ot*SAT_C + ot*SAP_mm + ot*cv +
                    pt*SAT_C + pt*SAP_mm + pt*cv +
                    (1|tray),
                  control = glmerControl(optimizer ="bobyqa", optCtrl=list(maxfun=2000000)),
                  data = scaled3,
                  na.action = na.exclude)
# assess
check_model(model)
plot(model, add.smooth = FALSE, which = 1)
E <- resid(model)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(model))
qqline(resid(model))              

summary(model)
r.squaredGLMM(model)
Anova(model, type = 3)

#AIC
model1 <- glmer.nb(days_to_flower ~ ot + pt + SAT_C + SAP_mm + cv +
                     ot*pt +
                     ot*SAT_C + ot*SAP_mm + ot*cv +
                     pt*SAT_C + pt*SAP_mm + cv +
                     (1|tray),
                   control = glmerControl(optimizer ="bobyqa", optCtrl=list(maxfun=2000000)),
                   data = scaled3,
                   na.action = na.exclude)
summary(model1)

model2 <- glmer.nb(days_to_flower ~ ot + pt + SAT_C + SAP_mm + cv +
                     ot*pt +
                     ot*SAT_C + ot*SAP_mm + ot*cv +
                     pt*SAT_C + pt*SAP_mm +
                     (1|tray),
                   control = glmerControl(optimizer ="bobyqa", optCtrl=list(maxfun=2000000)),
                   data = scaled3,
                   na.action = na.exclude)
summary(model2)

model3 <- glmer.nb(days_to_flower ~ ot + pt + SAT_C + SAP_mm + cv +
                     ot*pt +
                     ot*SAT_C + ot*SAP_mm + ot*cv +
                     pt*SAT_C +
                     (1|tray),
                   control = glmerControl(optimizer ="bobyqa", optCtrl=list(maxfun=2000000)),
                   data = scaled3,
                   na.action = na.exclude)
summary(model3)

model4 <- glmer.nb(days_to_flower ~ ot + pt + SAT_C + SAP_mm + cv +
                     ot*pt +
                     ot*SAT_C + ot*SAP_mm + ot*cv +
                     (1|tray),
                   control = glmerControl(optimizer ="bobyqa", optCtrl=list(maxfun=2000000)),
                   data = scaled3,
                   na.action = na.exclude)
summary(model4)

model5 <- glmer.nb(days_to_flower ~ ot + pt + SAT_C + SAP_mm + cv +
                     ot*pt +
                     ot*SAT_C + ot*SAP_mm +
                     (1|tray),
                   control = glmerControl(optimizer ="bobyqa", optCtrl=list(maxfun=2000000)),
                   data = scaled3,
                   na.action = na.exclude)
summary(model5)

model6 <- glmer.nb(days_to_flower ~ ot + pt + SAT_C + SAP_mm + cv +
                     ot*pt +
                     ot*SAT_C +
                     (1|tray),
                   control = glmerControl(optimizer ="bobyqa", optCtrl=list(maxfun=2000000)),
                   data = scaled3,
                   na.action = na.exclude)
summary(model6)

model7 <- glmer.nb(days_to_flower ~ ot + pt + SAT_C + SAP_mm + cv +
                     ot*pt +
                     (1|tray),
                   control = glmerControl(optimizer ="bobyqa", optCtrl=list(maxfun=2000000)),
                   data = scaled3,
                   na.action = na.exclude)
summary(model7)

model8 <- glmer.nb(days_to_flower ~ ot + pt + SAT_C + SAP_mm + cv +
                     (1|tray),
                   control = glmerControl(optimizer ="bobyqa", optCtrl=list(maxfun=2000000)),
                   data = scaled3,
                   na.action = na.exclude)
summary(model8)

AIC(model1, model2, model3, model4, model5, model6, model7, model8)

plot(model7, add.smooth = FALSE, which = 1)
E <- resid(model7)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(model7))
qqline(resid(model7))              

summary(model7)
r.squaredGLMM(model7)
Anova(model7, type = 3)
