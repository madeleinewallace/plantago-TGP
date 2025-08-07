# title: 03.1_final-models-Q1

# about: final models chosen for each trait response. after many years.
# changed the model structure from 03. simplified it per rachel feedback
# just ot, pt, and interaction. removed tray as random effect. added pop as random effect.

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
source('code-and-data/scripts/plpa-manuscript/01_clean.R')
library(pscl) #can't remember :D
library(glmmTMB)
library(multcomp)
library(multcompView)



# root biomass -----------------------------------------------------------------

# model
biomassroot_model1 <- lmer(root ~ ot + pt + ot:pt + (1|pop), 
                  data=BIOMASS_FINAL, na.action = na.exclude)

# check diagnostics
residuals_model1 <- simulateResiduals(biomassroot_model1)
plot(residuals_model1)
testUniformity(residuals_model1)
testDispersion(residuals_model1)
testOutliers(residuals_model1)
check_model(biomassroot_model1)
plot(HEIGHT_FINAL$pop, resid(biomassroot_model1), xlab = "pop", ylab = "residuals")
plot(HEIGHT_FINAL$`2023_treat`, resid(biomassroot_model1), xlab = "2023 treat", ylab = "residuals")

# final stats
summary(biomassroot_model1)
root_anova <- Anova(biomassroot_model1, type="III",  test.statistic= "F")
root_r2 <- r.squaredGLMM(biomassroot_model1)
root_anova
root_r2

# contrasts
rootcontrast <- emmeans(biomassroot_model1, specs = pairwise ~ ot:pt, type = "response")
rootcontrast$contrasts

rootcontrast1 <- emmeans(biomassroot_model1, specs = pairwise ~ ot, type = "response")
rootcontrast1$contrasts



# shoot biomass-----------------------------------------------------------------

# model
biomassshoot_model1 <- lmer(shoot ~ ot + pt + ot:pt + (1|pop), 
                           data=BIOMASS_FINAL, na.action = na.exclude)

# check diagnostics
residuals_model1 <- simulateResiduals(biomassshoot_model1)
plot(residuals_model1)
testUniformity(residuals_model1)
testDispersion(residuals_model1)
testOutliers(residuals_model1)
check_model(biomassshoot_model1)
plot(HEIGHT_FINAL$pop, resid(biomassshoot_model1), xlab = "pop", ylab = "residuals")
plot(HEIGHT_FINAL$`2023_treat`, resid(biomassshoot_model1), xlab = "2023 treat", ylab = "residuals")

# final stats
summary(biomassshoot_model1)
shoot_anova <- Anova(biomassshoot_model1, type="III",  test.statistic= "F")
shoot_r2 <- r.squaredGLMM(biomassshoot_model1)
shoot_anova
shoot_r2

# contrasts
shootcontrast <- emmeans(biomassshoot_model1, specs = pairwise ~ ot:pt, type = "response")
shootcontrast$contrasts

shootcontrast1 <- emmeans(biomassshoot_model1, specs = pairwise ~ ot, type = "response")
shootcontrast1$contrasts



# total biomass ----------------------------------------------------------------

# model
biomasstotal_model1 <- lmer(total_biomass ~ ot + pt + ot:pt + (1|pop), 
                            data=BIOMASS_FINAL, na.action = na.exclude)

# check diagnostics
residuals_model1 <- simulateResiduals(biomasstotal_model1)
plot(residuals_model1)
testUniformity(residuals_model1)
testDispersion(residuals_model1)
testOutliers(residuals_model1)
check_model(biomasstotal_model1)
plot(HEIGHT_FINAL$pop, resid(biomasstotal_model1), xlab = "pop", ylab = "residuals")
plot(HEIGHT_FINAL$`2023_treat`, resid(biomasstotal_model1), xlab = "2023 treat", ylab = "residuals")

# final stats
summary(biomasstotal_model1)
total_anova <- Anova(biomasstotal_model1, type="III",  test.statistic= "F")
total_r2 <-r.squaredGLMM(biomasstotal_model1)
total_anova
total_r2

#contrasts
totalcontrast <- emmeans(biomasstotal_model1, specs = pairwise ~ ot:pt, type = "response")
totalcontrast$contrasts

totalcontrast1 <- emmeans(biomasstotal_model1, specs = pairwise ~ ot, type = "response")
totalcontrast1$contrasts

# root:shoot ratio -------------------------------------------------------------
biomassRS_model1 <- lmer(log_RS ~ ot + pt + ot:pt + (1|pop), 
                            data=BIOMASS_FINAL, na.action = na.exclude)

# check diagnostics
residuals_model1 <- simulateResiduals(biomassRS_model1)
plot(residuals_model1)
testUniformity(residuals_model1)
testDispersion(residuals_model1)
testOutliers(residuals_model1)
check_model(biomassRS_model1)
plot(HEIGHT_FINAL$pop, resid(biomassRS_model1), xlab = "pop", ylab = "residuals")
plot(HEIGHT_FINAL$`2023_treat`, resid(biomassRS_model1), xlab = "2023 treat", ylab = "residuals")

# final stats
summary(biomassRS_model1)
rs_anova <- Anova(biomassRS_model1, type="III",  test.statistic= "F")
rs_r2 <- r.squaredGLMM(biomassRS_model1)
rs_anova
rs_r2

#contrasts
rscontrast <- emmeans(biomassRS_model1, specs = pairwise ~ ot:pt, type = "response")
rscontrast$contrasts

rscontrast1 <- emmeans(biomassRS_model1, specs = pairwise ~ ot, type = "response")
rscontrast1$contrasts

# max height -------------------------------------------------------------------

# model
maxheight_model1 <- lmer(max ~ ot + pt + ot:pt + (1|pop), 
                         data=HEIGHT_FINAL, na.action = na.exclude)

# check diagnostics
residuals_model1 <- simulateResiduals(maxheight_model1)
plot(residuals_model1)
testUniformity(residuals_model1)
testDispersion(residuals_model1)
testOutliers(residuals_model1)
check_model(maxheight_model1)
plot(HEIGHT_FINAL$pop, resid(maxheight_model1), xlab = "pop", ylab = "residuals")
plot(HEIGHT_FINAL$`2023_treat`, resid(maxheight_model1), xlab = "2023 treat", ylab = "residuals")

# final stats
summary(maxheight_model1)
height_anova <- Anova(maxheight_model1, type="III",  test.statistic= "F")
height_r2 <- r.squaredGLMM(maxheight_model1)
height_anova
height_r2


# RGR --------------------------------------------------------------------------

# model
rgr_model1 <- lmer(RGR ~ ot + pt + ot:pt + (1|pop), 
                         data=RGR_FINAL, na.action = na.exclude)

# check diagnostics
residuals_model1 <- simulateResiduals(rgr_model1)
plot(residuals_model1)
testUniformity(residuals_model1)
testDispersion(residuals_model1)
testOutliers(residuals_model1)
check_model(rgr_model1)
plot(HEIGHT_FINAL$pop, resid(rgr_model1), xlab = "pop", ylab = "residuals")
plot(HEIGHT_FINAL$`2023_treat`, resid(rgr_model1), xlab = "2023 treat", ylab = "residuals")

# final stats
summary(rgr_model1)
rgr_anova <- Anova(rgr_model1, type="III",  test.statistic= "F")
rgr_r2 <- r.squaredGLMM(rgr_model1)
rgr_anova
rgr_r2

rgrcontrast1 <- emmeans(rgr_model1, specs = pairwise ~ ot, type = "response")
rgrcontrast1$contrasts


# mortality --------------------------------------------------------------------

# model
mort_model1 <- glmer(status ~ ot + pt + ot:pt + (1|pop),
                     data = MORT_DAY50, family = binomial(link = "logit"), na.action = na.exclude)


# check diagnostics
residuals_model1 <- simulateResiduals(mort_model1)
plot(residuals_model1)
testUniformity(residuals_model1)
testDispersion(residuals_model1)
testOutliers(residuals_model1)
check_model(mort_model1)
plot(HEIGHT_FINAL$pop, resid(mort_model1), xlab = "pop", ylab = "residuals")
plot(HEIGHT_FINAL$`2023_treat`, resid(mort_model1), xlab = "2023 treat", ylab = "residuals")

# final stats
summary(mort_model1)
mort_anova <- Anova(mort_model1, type = "III", test.statistic = "Chisq")
mort_r2 <- r.squaredGLMM(mort_model1)
mort_anova
mort_r2


# sla --------------------------------------------------------------------------

# model
sla_model1 <- lmer(log(sla) ~ ot + pt + ot:pt + (1|pop),
                     data = SLA_LDMC_FINAL, na.action = na.exclude)


# check diagnostics
residuals_model1 <- simulateResiduals(sla_model1)
plot(residuals_model1)
testUniformity(residuals_model1)
testDispersion(residuals_model1)
testOutliers(residuals_model1)
check_model(sla_model1)
plot(HEIGHT_FINAL$pop, resid(sla_model1), xlab = "pop", ylab = "residuals")
plot(HEIGHT_FINAL$`2023_treat`, resid(sla_model1), xlab = "2023 treat", ylab = "residuals")

# final stats
summary(sla_model1)
sla_anova <- Anova(sla_model1, type="III",  test.statistic= "F")
sla_r2 <- r.squaredGLMM(sla_model1)
sla_anova
sla_r2


# ldmc -------------------------------------------------------------------------

# model
ldmc_model1 <- lmer(log(ldmc) ~ ot + pt + ot:pt + (1|pop),
                   data = SLA_LDMC_FINAL, na.action = na.exclude)

# check diagnostics
residuals_model1 <- simulateResiduals(ldmc_model1)
plot(residuals_model1)
testUniformity(residuals_model1)
testDispersion(residuals_model1)
testOutliers(residuals_model1)
check_model(ldmc_model1)
plot(HEIGHT_FINAL$pop, resid(ldmc_model1), xlab = "pop", ylab = "residuals")
plot(HEIGHT_FINAL$`2023_treat`, resid(ldmc_model1), xlab = "2023 treat", ylab = "residuals")

# final stats
summary(ldmc_model1)
ldmc_anova <- Anova(ldmc_model1, type="III",  test.statistic= "F")
ldmc_r2 <- r.squaredGLMM(ldmc_model1)
ldmc_anova
ldmc_r2

ldmccontrast1 <- emmeans(ldmc_model1, specs = pairwise ~ ot, type = "response")
ldmccontrast1$contrasts


# num of plants that flowered --------------------------------------------------
# here, very zero-inflated

# checking for overdispersion
mean_status <- mean(FLOWER_STATUS$status, na.rm = TRUE)
var_status <- var(FLOWER_STATUS$status, na.rm = TRUE)
dispersion_ratio <- var_status / mean_status
dispersion_ratio
# less than one, not overdispersed, no nb needed

# also, gonna use a hurdle model -> response of 0 is part of the ecological process
# hurdle model assumes zero values arise from a distinct structural process
# eg zeros are deterministic (non-flowering plants for biological reasons)

# model
flowerstatus_model1 <- glmmTMB(status ~ ot + pt + ot:pt + (1 | pop), 
                               ziformula = ~ot + pt, 
                               family = binomial, 
                               data = FLOWER_STATUS)

# check diagnostics
residuals_model1 <- simulateResiduals(flowerstatus_model1)
plot(residuals_model1)
testUniformity(residuals_model1)
testDispersion(residuals_model1)
testOutliers(residuals_model1)
check_model(flowerstatus_model1)
plot(HEIGHT_FINAL$pop, resid(flowerstatus_model1), xlab = "pop", ylab = "residuals")
plot(HEIGHT_FINAL$`2023_treat`, resid(flowerstatus_model1), xlab = "2023 treat", ylab = "residuals")

# final stats
summary(flowerstatus_model1)
flowstat_anova <- Anova(flowerstatus_model1, type = "III", test.statistic = "Chisq")
flowstat_r2 <-r.squaredGLMM(flowerstatus_model1)

car::Anova(flowerstatus_model1, type = "III", test.statistic = "Chisq", component="cond")
car::Anova(flowerstatus_model1, type = "III", test.statistic = "Chisq", component="zi")
performance::r2(flowerstatus_model1)

emm <- emmeans(flowerstatus_model1, ~ ot * pt, type = "response")
pairs(emm, adjust = "tukey")

emm <- emmeans(flowerstatus_model1, ~ pt, type = "response")
pairs(emm, adjust = "tukey")



# num of flowering structures --------------------------------------------------
# here, out of ALL the plants that survived - so very zero inflated !

# checking for overdispersion
mean_num_structure <- mean(FLOWER_FINAL$num_structure, na.rm = TRUE)
var_num_structure <- var(FLOWER_FINAL$num_structure, na.rm = TRUE)
dispersion_ratio <- var_num_structure / mean_num_structure
dispersion_ratio
# greater than 1, negative binomial model needed

# also, gonna use a hurdle model -> response of 0 is part of the ecological process
# hurdle model assumes zero values arise from a distinct structural process
# eg zeros are deterministic (non-flowering plants for biological reasons)

# model
flowerstructure_model1 <- glmmTMB(num_structure ~ ot + pt + ot:pt + (1|pop), 
                                        ziformula = ~ot + pt, 
                                        family = truncated_nbinom2, 
                                        data = FLOWER_FINAL)

# check diagnostics
residuals_model1 <- simulateResiduals(flowerstructure_model1)
plot(residuals_model1)
testUniformity(residuals_model1)
testDispersion(residuals_model1)
testOutliers(residuals_model1)
check_model(flowerstructure_model1)
plot(HEIGHT_FINAL$pop, resid(flowerstructure_model1), xlab = "pop", ylab = "residuals")
plot(HEIGHT_FINAL$`2023_treat`, resid(flowerstructure_model1), xlab = "2023 treat", ylab = "residuals")

# final stats
summary(flowerstructure_model1)
flowstruc_anova <- Anova(flowerstructure_model1, type = "III", test.statistic = "Chisq")
flowstruc_r2 <- r.squaredGLMM(flowerstructure_model1)

car::Anova(flowerstructure_model1, type = "III", test.statistic = "Chisq", component="cond")
car::Anova(flowerstructure_model1, type = "III", test.statistic = "Chisq", component="zi")

# conditional model r2
performance::r2(flowerstatus_model1)

# WOWWWW offspring drought treatments reduce the number of structures, but parental drought has no effect
# parental drought reduces the odds of having a structural zero - parental drought reduces the odds of NOT flowering
# amazing

flowercontrast1 <- emmeans(flowerstructure_model1, specs = pairwise ~ ot, type = "response")
flowercontrast1$contrasts


# days to flower ---------------------------------------------------------------

# checking data
hist(FLOWER_FINAL$days_to_flower)
# checking for overdispersion
mean_days <- mean(FLOWER_FINAL$days_to_flower, na.rm = TRUE)
var_days <- var(FLOWER_FINAL$days_to_flower, na.rm = TRUE)
dispersion_ratio <- var_days / mean_days
dispersion_ratio
# dispersion ratio over 1, which means data is overdispersed

# model
flowerdays_model1 <- glmmTMB(days_to_flower ~ ot + pt + ot:pt + (1 | pop),
                              family = nbinom2, 
                              control = glmmTMBControl(optimizer = "nlminb"),
                              data = FLOWER_FINAL)

# check diagnostics
residuals_model1 <- simulateResiduals(flowerdays_model1)
plot(residuals_model1)
testUniformity(residuals_model1)
testDispersion(residuals_model1)
testOutliers(residuals_model1)
check_model(flowerdays_model1)
plot(HEIGHT_FINAL$pop, resid(flowerdays_model1), xlab = "pop", ylab = "residuals")
plot(HEIGHT_FINAL$`2023_treat`, resid(flowerdays_model1), xlab = "2023 treat", ylab = "residuals")

# final stats
summary(flowerdays_model1)
flowdays_anova <- Anova(flowerdays_model1, type = "III", test.statistic = "Chisq")
flowdays_r2 <- r.squaredGLMM(flowerdays_model1)
flowdays_anova
flowdays_r2


# seed mass --------------------------------------------------------------------

# checking data
hist(SEED_FINAL$`2023_mass_total_g`)
summary(SEED_FINAL$`2023_mass_total_g`)

# model
# data is continuous, positive, and skewed = gamma dist
seedmass_model1 <- glmer(mass_total ~ ot + pt + ot:pt + (1 | pop),
                         data = SEED_FINAL,
                         family = Gamma(link = "log"),
                         control = glmerControl(optimizer = "bobyqa"))



# check diagnostics
residuals_model1 <- simulateResiduals(seedmass_model1)
plot(residuals_model1)
testUniformity(residuals_model1)
testDispersion(residuals_model1)
testOutliers(residuals_model1)
check_model(seedmass_model1)
plot(HEIGHT_FINAL$pop, resid(seedmass_model1), xlab = "pop", ylab = "residuals")
plot(HEIGHT_FINAL$`2023_treat`, resid(seedmass_model1), xlab = "2023 treat", ylab = "residuals")

# final stats
summary(seedmass_model1)
seedmass_anova <- Anova(seedmass_model1, type = "III", test.statistic = "Chisq")
seedmass_r2 <- r.squaredGLMM(seedmass_model1)
seedmass_anova
seedmass_r2

emm <- emmeans(seedmass_model1, ~ ot * pt, type = "response")
pairs(emm, adjust = "tukey")

# seed num ---------------------------------------------------------------------

# checking for overdispersion
mean_status <- mean(SEED_FINAL$num_total, na.rm = TRUE)
var_status <- var(SEED_FINAL$num_total, na.rm = TRUE)
dispersion_ratio <- var_status / mean_status
dispersion_ratio
# VERY overdispersed lolz

# checking data
hist(SEED_FINAL$num_total)
summary(SEED_FINAL$num_total)

# model
seednum_model1 <- glmmTMB(num_total ~ ot + pt + ot:pt + (1 | pop),
                             family = nbinom2, 
                             control = glmmTMBControl(optimizer = "nlminb"),
                             data = SEED_FINAL)

# check diagnostics
residuals_model1 <- simulateResiduals(seednum_model1)
plot(residuals_model1)
testUniformity(residuals_model1)
testDispersion(residuals_model1)
testOutliers(residuals_model1)
check_model(seednum_model1)
plot(HEIGHT_FINAL$pop, resid(seednum_model1), xlab = "pop", ylab = "residuals")
plot(HEIGHT_FINAL$`2023_treat`, resid(seednum_model1), xlab = "2023 treat", ylab = "residuals")

# final stats
summary(seednum_model1)
seednum_anova <- Anova(seednum_model1, type = "III", test.statistic = "Chisq")
seednum_r2 <- r.squaredGLMM(seednum_model1)
seednum_anova
seednum_r2

emm <- emmeans(seednum_model1, ~ ot * pt, type = "response")
pairs <- pairs(emm, adjust = "tukey")
pairs




# model summaries --------------------------------------------------------------
summary(biomassroot_model1)
root_anova <- Anova(biomassroot_model1, type="III",  test.statistic= "F")
root_r2 <- r.squaredGLMM(biomassroot_model1)
root_anova
root_r2

summary(biomassshoot_model1)
shoot_anova <- Anova(biomassshoot_model1, type="III",  test.statistic= "F")
shoot_r2 <- r.squaredGLMM(biomassshoot_model1)
shoot_anova
shoot_r2

summary(biomasstotal_model1)
total_anova <- Anova(biomasstotal_model1, type="III",  test.statistic= "F")
total_r2 <-r.squaredGLMM(biomasstotal_model1)
total_anova
total_r2

summary(biomassRS_model1)
rs_anova <- Anova(biomassRS_model1, type="III",  test.statistic= "F")
rs_r2 <- r.squaredGLMM(biomassRS_model1)
rs_anova
rs_r2

summary(maxheight_model1)
height_anova <- Anova(maxheight_model1, type="III",  test.statistic= "F")
height_r2 <- r.squaredGLMM(maxheight_model1)
height_anova
height_r2

summary(rgr_model1)
rgr_anova <- Anova(rgr_model1, type="III",  test.statistic= "F")
rgr_r2 <- r.squaredGLMM(rgr_model1)
rgr_anova
rgr_r2

summary(mort_model1)
mort_anova <- Anova(mort_model1, type = "III", test.statistic = "Chisq")
mort_r2 <- r.squaredGLMM(mort_model1)
mort_anova
mort_r2

summary(sla_model1)
sla_anova <- Anova(sla_model1, type="III",  test.statistic= "F")
sla_r2 <- r.squaredGLMM(sla_model1)
sla_anova
sla_r2

summary(ldmc_model1)
ldmc_anova <- Anova(ldmc_model1, type="III",  test.statistic= "F")
ldmc_r2 <- r.squaredGLMM(ldmc_model1)
ldmc_anova
ldmc_r2

summary(flowerstatus_model1)
flowstat_anova <- Anova(flowerstatus_model1, type = "III", test.statistic = "Chisq")
flowstat_r2 <-r.squaredGLMM(flowerstatus_model1)

car::Anova(flowerstatus_model1, type = "III", test.statistic = "Chisq", component="cond")
car::Anova(flowerstatus_model1, type = "III", test.statistic = "Chisq", component="zi")
performance::r2(flowerstatus_model1)

summary(flowerstructure_model1)
flowstruc_anova <- Anova(flowerstructure_model1, type = "III", test.statistic = "Chisq")
flowstruc_r2 <- r.squaredGLMM(flowerstructure_model1)

car::Anova(flowerstructure_model1, type = "III", test.statistic = "Chisq", component="cond")
car::Anova(flowerstructure_model1, type = "III", test.statistic = "Chisq", component="zi")

# conditional model r2
performance::r2(flowerstatus_model1)

summary(flowerdays_model1)
flowdays_anova <- Anova(flowerdays_model1, type = "III", test.statistic = "Chisq")
flowdays_r2 <- r.squaredGLMM(flowerdays_model1)
flowdays_anova
flowdays_r2

summary(seedmass_model1)
seedmass_anova <- Anova(seedmass_model1, type = "III", test.statistic = "Chisq")
seedmass_r2 <- r.squaredGLMM(seedmass_model1)
seedmass_anova
seedmass_r2

summary(seednum_model1)
seednum_anova <- Anova(seednum_model1, type = "III", test.statistic = "Chisq")
seednum_r2 <- r.squaredGLMM(seednum_model1)
seednum_anova
seednum_r2
