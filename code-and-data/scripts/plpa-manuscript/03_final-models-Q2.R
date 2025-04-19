# title: 04.1_final-models-population

# here, i explore how seed source location might affect trait response
# specifically spring cv!

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

# i think i am going to use spring VPD CV
# combines SAP and SAT into one term
# per cocchlico & herman

# The coefficient of variation (CV) for spring VPD (spring_vpd) 
# was calculated over a 30-year period (1989–2019) by dividing the standard deviation 
# of mean monthly vapor pressure deficit values for the spring months (April, May, June) 
# by the corresponding mean and expressing the result as a percentage, to quantify relative variability.

# root biomass -----------------------------------------------------------------

# scale cv term
BIOMASS_FINAL$scv_scaled <- scale(BIOMASS_FINAL$cv_term_spring)
BIOMASS_FINAL$vpd_scaled <- scale(BIOMASS_FINAL$spring_vpd_cv)

# model
biomassroot_1 <- lmer(root ~ ot * pt * vpd_scaled + (1|pop), 
                           data=BIOMASS_FINAL, na.action = na.exclude)
vif(biomassroot_1)

# check diagnostics
residuals_model1 <- simulateResiduals(biomassroot_1)
plot(residuals_model1)
testUniformity(residuals_model1)
testDispersion(residuals_model1)
testOutliers(residuals_model1)
check_model(biomassroot_1)

# final stats
summary(biomassroot_1)
root_anova <- Anova(biomassroot_1, type="III",  test.statistic= "F")
root_r2 <- r.squaredGLMM(biomassroot_1)
root_anova
root_r2

# contrasts
rootcontrast <- emmeans(biomassroot_1, specs = pairwise ~ ot:pt, type = "response")
rootcontrast$contrasts

root_slopes <- emtrends(biomassroot_1, pairwise ~ ot:pt, var="vpd_scaled")
root_slopes

# shoot biomass -----------------------------------------------------------------

# scale cv term
BIOMASS_FINAL$vpd_scaled <- scale(BIOMASS_FINAL$spring_vpd_cv)

# model
biomassshoot_1 <- lmer(shoot ~ ot * pt * vpd_scaled + (1|pop), 
                      data=BIOMASS_FINAL, na.action = na.exclude)
vif(biomassshoot_1)

# check diagnostics
residuals_model1 <- simulateResiduals(biomassshoot_1)
plot(residuals_model1)
testUniformity(residuals_model1)
testDispersion(residuals_model1)
testOutliers(residuals_model1)
check_model(biomassshoot_1)

# final stats
summary(biomassshoot_1)
shoot_anova <- Anova(biomassshoot_1, type="III",  test.statistic= "F")
shoot_r2 <- r.squaredGLMM(biomassshoot_1)
shoot_anova
shoot_r2

# contrasts
shootcontrast <- emmeans(biomassshoot_1, specs = pairwise ~ ot:pt, type = "response")
shootcontrast$contrasts



# total biomass -----------------------------------------------------------------

# scale cv term
BIOMASS_FINAL$vpd_scaled <- scale(BIOMASS_FINAL$spring_vpd_cv)

# model
biomasstotal_1 <- lmer(total_biomass ~ ot * pt * vpd_scaled + (1|pop), 
                       data=BIOMASS_FINAL, na.action = na.exclude)
vif(biomasstotal_1)

# check diagnostics
residuals_model1 <- simulateResiduals(biomasstotal_1)
plot(residuals_model1)
testUniformity(residuals_model1)
testDispersion(residuals_model1)
testOutliers(residuals_model1)
check_model(biomasstotal_1)

# final stats
summary(biomasstotal_1)
total_anova <- Anova(biomasstotal_1, type="III",  test.statistic= "F")
total_r2 <- r.squaredGLMM(biomasstotal_1)
total_anova
total_r2

# contrasts
totalcontrast <- emmeans(biomasstotal_1, specs = pairwise ~ ot:pt, type = "response")
totalcontrast$contrasts

total_slopes <- emtrends(biomasstotal_1, pairwise ~ ot:pt, var="vpd_scaled")
total_slopes



# root:shoot ratio -------------------------------------------------------------

# scale cv term
BIOMASS_FINAL$vpd_scaled <- scale(BIOMASS_FINAL$spring_vpd_cv)

# model
biomassrs_1 <- lmer(log_RS ~ ot * pt * vpd_scaled + (1|pop), 
                       data=BIOMASS_FINAL, na.action = na.exclude)
vif(biomasstotal_1)

# check diagnostics
residuals_model1 <- simulateResiduals(biomassrs_1)
plot(residuals_model1)
testUniformity(residuals_model1)
testDispersion(residuals_model1)
testOutliers(residuals_model1)
check_model(biomassrs_1)

# final stats
summary(biomassrs_1)
rs_anova <- Anova(biomassrs_1, type="III",  test.statistic= "F")
rs_r2 <- r.squaredGLMM(biomassrs_1)
rs_anova
rs_r2

rs_slopes <- emtrends(biomassrs_1, pairwise ~ ot:pt, var="vpd_scaled")
rs_slopes



# max height -------------------------------------------------------------------

# scale cv term
HEIGHT_FINAL$vpd_scaled <- scale(HEIGHT_FINAL$spring_vpd_cv)

# model
height_1 <- lmer(max ~ ot * pt * vpd_scaled + (1|pop), 
                    data=HEIGHT_FINAL, na.action = na.exclude)
vif(height_1)

# check diagnostics
residuals_model1 <- simulateResiduals(height_1)
plot(residuals_model1)
testUniformity(residuals_model1)
testDispersion(residuals_model1)
testOutliers(residuals_model1)
check_model(height_1)

# final stats
summary(height_1)
height_anova <- Anova(height_1, type="III",  test.statistic= "F")
height_r2 <- r.squaredGLMM(height_1)
height_anova
height_r2


height_slopes <- emtrends(height_1, pairwise ~ ot:pt, var="vpd_scaled")
height_slopes


# contrasts
heightcontrast <- emmeans(height_1, specs = pairwise ~ ot:pt, type = "response")
heightcontrast$contrasts



# RGR --------------------------------------------------------------------------

# scale cv term
RGR_FINAL$vpd_scaled <- scale(RGR_FINAL$spring_vpd_cv)

# model
rgr_1 <- lmer(RGR ~ ot * pt * vpd_scaled + (1|pop), 
                 data=RGR_FINAL, na.action = na.exclude)
vif(rgr_1)

# check diagnostics
residuals_model1 <- simulateResiduals(rgr_1)
plot(residuals_model1)
testUniformity(residuals_model1)
testDispersion(residuals_model1)
testOutliers(residuals_model1)
check_model(rgr_1)

# final stats
summary(rgr_1)
rgr_anova <- Anova(rgr_1, type="III",  test.statistic= "F")
rgr_r2 <- r.squaredGLMM(rgr_1)
rgr_anova
rgr_r2

# contrasts
rgrcontrast <- emmeans(rgr_1, specs = pairwise ~ ot:pt, type = "response")
rgrcontrast$contrasts

rgr_slopes <- emtrends(rgr_1, pairwise ~ ot:pt, var="vpd_scaled")
rgr_slopes



# mortality --------------------------------------------------------------------

# scale cv term
MORT_DAY50$vpd_scaled <- scale(MORT_DAY50$spring_vpd_cv)

# model
mort_1 <- glmer(status ~ ot * pt * vpd_scaled + (1|pop), 
               data = MORT_DAY50, family = binomial(link = "logit"), na.action = na.exclude)
vif(mort_1)

# check diagnostics
residuals_model1 <- simulateResiduals(mort_1)
plot(residuals_model1)
testUniformity(residuals_model1)
testDispersion(residuals_model1)
testOutliers(residuals_model1)
check_model(mort_1)

# final stats
summary(mort_1)
mort_anova <- Anova(mort_1, type = "III", test.statistic = "Chisq")
mort_r2 <- r.squaredGLMM(mort_1)
mort_anova
mort_r2

# contrasts
mortcontrast <- emmeans(mort_1, specs = pairwise ~ ot:pt, type = "response")
mortcontrast$contrasts



# SLA --------------------------------------------------------------------------

# scale cv term
SLA_LDMC_FINAL$vpd_scaled <- scale(SLA_LDMC_FINAL$spring_vpd_cv)

# model
sla_1 <- lmer(log(sla) ~ ot * pt * vpd_scaled + (1|pop), 
              data=SLA_LDMC_FINAL, na.action = na.exclude)
vif(sla_1)

# check diagnostics
residuals_model1 <- simulateResiduals(sla_1)
plot(residuals_model1)
testUniformity(residuals_model1)
testDispersion(residuals_model1)
testOutliers(residuals_model1)
check_model(sla_1)

# final stats
summary(sla_1)
sla_anova <- Anova(sla_1, type="III",  test.statistic= "F")
sla_r2 <- r.squaredGLMM(sla_1)
sla_anova
sla_r2

# contrasts
slacontrast <- emmeans(sla_1, specs = pairwise ~ ot:pt, type = "response")
slacontrast$contrasts



# LDMC -------------------------------------------------------------------------

# scale cv term
SLA_LDMC_FINAL$vpd_scaled <- scale(SLA_LDMC_FINAL$spring_vpd_cv)

# model
ldmc_1 <- lmer(log(ldmc) ~ ot * pt * vpd_scaled + (1|pop), 
              data=SLA_LDMC_FINAL, na.action = na.exclude)
vif(ldmc_1)

# check diagnostics
residuals_model1 <- simulateResiduals(ldmc_1)
plot(residuals_model1)
testUniformity(residuals_model1)
testDispersion(residuals_model1)
testOutliers(residuals_model1)
check_model(ldmc_1)

# final stats
summary(ldmc_1)
ldmc_anova <- Anova(ldmc_1, type="III",  test.statistic= "F")
ldmc_r2 <- r.squaredGLMM(ldmc_1)
ldmc_anova
ldmc_r2

# contrasts
ldmccontrast <- emmeans(ldmc_1, specs = pairwise ~ ot:pt, type = "response")
ldmccontrast$contrasts



# num of plants that flowered --------------------------------------------------

# scale cv term
FLOWER_STATUS$vpd_scaled <- scale(FLOWER_STATUS$spring_vpd_cv)

# model
flowerstatus_1 <- glmmTMB(status ~ ot * pt * vpd_scaled + (1|pop), 
                               ziformula = ~ot + pt + vpd_scaled, 
                               family = binomial, 
                               data = FLOWER_STATUS)
vif(flowerstatus_1)

# check diagnostics
residuals_model1 <- simulateResiduals(flowerstatus_1)
plot(residuals_model1)
testUniformity(residuals_model1)
testDispersion(residuals_model1)
testOutliers(residuals_model1)
check_model(flowerstatus_1)

# final stats
summary(flowerstatus_1)
car::Anova(flowerstatus_1, type = "III", test.statistic = "Chisq", component="cond")
car::Anova(flowerstatus_1, type = "III", test.statistic = "Chisq", component="zi")
performance::r2(flowerstatus_1)

# contrasts
flowerstatuscontrast <- emmeans(flowerstatus_1, specs = pairwise ~ ot:pt, type = "response", var="vpd_scaled")
flowerstatuscontrast$contrasts

flow_slopes <- emtrends(flowerstatus_1, pairwise ~ ot:pt, var="vpd_scaled")
flow_slopes



# num of flowering structures --------------------------------------------------

# scale cv term
FLOWER_FINAL$vpd_scaled <- scale(FLOWER_FINAL$spring_vpd_cv)

# model
flowerstructure_1 <- glmmTMB(num_structure ~ ot * pt * vpd_scaled + (1|pop), 
                                  ziformula = ~ot + pt + vpd_scaled, 
                                  family = truncated_nbinom2, 
                                  data = FLOWER_FINAL)
vif(flowerstructure_1)

# check diagnostics
residuals_model1 <- simulateResiduals(flowerstructure_1)
plot(residuals_model1)
testUniformity(residuals_model1)
testDispersion(residuals_model1)
testOutliers(residuals_model1)
check_model(flowerstructure_1)

# final stats
summary(flowerstructure_1)
car::Anova(flowerstructure_1, type = "III", test.statistic = "Chisq", component="cond")
car::Anova(flowerstructure_1, type = "III", test.statistic = "Chisq", component="zi")
performance::r2(flowerstructure_1)

# contrasts
flowerstructurecontrast <- emmeans(flowerstructure_1, specs = pairwise ~ ot:pt, type = "response")
flowerstructurecontrast$contrasts



# days to flower ---------------------------------------------------------------

# scale cv term
FLOWER_FINAL$vpd_scaled <- scale(FLOWER_FINAL$spring_vpd_cv)

# model
flowerdays_1 <- glmmTMB(days_to_flower ~ ot * pt * vpd_scaled + (1|pop),
                             family = nbinom2, 
                             control = glmmTMBControl(optimizer = "nlminb"),
                             data = FLOWER_FINAL)
vif(flowerdays_1)

# check diagnostics
residuals_model1 <- simulateResiduals(flowerdays_1)
plot(residuals_model1)
testUniformity(residuals_model1)
testDispersion(residuals_model1)
testOutliers(residuals_model1)
check_model(flowerdays_1)

# final stats
summary(flowerdays_1)
flowdays_anova <- Anova(flowerdays_1, type = "III", test.statistic = "Chisq")
flowdays_r2 <- r.squaredGLMM(flowerdays_1)
flowdays_anova
flowdays_r2

# contrasts
mortcontrast <- emmeans(mort_1, specs = pairwise ~ ot:pt, type = "response")
mortcontrast$contrasts



# seed mass --------------------------------------------------------------------

# scale cv term
SEED_FINAL$vpd_scaled <- scale(SEED_FINAL$spring_vpd_cv)

# model
seedmass_1 <- glmmTMB(mass_total ~ ot * pt * vpd_scaled + (1 | pop),
                        data = SEED_FINAL,
                        family = Gamma(link = "log"))


vif(seedmass_1)

# check diagnostics
residuals_model1 <- simulateResiduals(seedmass_1)
plot(residuals_model1)
testUniformity(residuals_model1)
testDispersion(residuals_model1)
testOutliers(residuals_model1)
check_model(seedmass_1)

# final stats
summary(seedmass_1)
seedmass_anova <- Anova(seedmass_1, type = "III", test.statistic = "Chisq")
seedmass_r2 <- r.squaredGLMM(seedmass_1)
seedmass_anova
seedmass_r2

# contrasts
seedmasscontrast <- emmeans(seedmass_1, specs = pairwise ~ ot:pt, type = "response")
seedmasscontrast$contrasts



# seed num ---------------------------------------------------------------------

# scale cv term
SEED_FINAL$vpd_scaled <- scale(SEED_FINAL$spring_vpd_cv)

# model
seednum_1 <- glmmTMB(num_total ~ ot * pt * vpd_scaled + (1 | pop),
                          family = nbinom2, 
                          control = glmmTMBControl(optimizer = "nlminb"),
                          data = SEED_FINAL)

# check diagnostics
residuals_model1 <- simulateResiduals(seednum_1)
plot(residuals_model1)
testUniformity(residuals_model1)
testDispersion(residuals_model1)
testOutliers(residuals_model1)
check_model(seednum_1)

# final stats
summary(seednum_1)
seednum_anova <- Anova(seednum_1, type = "III", test.statistic = "Chisq")
seednum_r2 <- r.squaredGLMM(seednum_1)
seednum_anova
seednum_r2

# contrasts
seednumcontrast <- emmeans(seednum_1, specs = pairwise ~ ot:pt, type = "response")
seednumcontrast$contrasts




# model summaries --------------------------------------------------------------
summary(biomassroot_1)
root_anova <- Anova(biomassroot_1, type="III",  test.statistic= "F")
root_r2 <- r.squaredGLMM(biomassroot_1)
root_anova
root_r2


summary(biomassshoot_1)
shoot_anova <- Anova(biomassshoot_1, type="III",  test.statistic= "F")
shoot_r2 <- r.squaredGLMM(biomassshoot_1)
shoot_anova
shoot_r2

summary(biomasstotal_1)
total_anova <- Anova(biomasstotal_1, type="III",  test.statistic= "F")
total_r2 <- r.squaredGLMM(biomasstotal_1)
total_anova
total_r2

summary(biomassrs_1)
rs_anova <- Anova(biomassrs_1, type="III",  test.statistic= "F")
rs_r2 <- r.squaredGLMM(biomassrs_1)
rs_anova
rs_r2

summary(height_1)
height_anova <- Anova(height_1, type="III",  test.statistic= "F")
height_r2 <- r.squaredGLMM(height_1)
height_anova
height_r2

summary(rgr_1)
rgr_anova <- Anova(rgr_1, type="III",  test.statistic= "F")
rgr_r2 <- r.squaredGLMM(rgr_1)
rgr_anova
rgr_r2

summary(mort_1)
mort_anova <- Anova(mort_1, type = "III", test.statistic = "Chisq")
mort_r2 <- r.squaredGLMM(mort_1)
mort_anova
mort_r2

summary(sla_1)
sla_anova <- Anova(sla_1, type="III",  test.statistic= "F")
sla_r2 <- r.squaredGLMM(sla_1)
sla_anova
sla_r2

summary(ldmc_1)
ldmc_anova <- Anova(ldmc_1, type="III",  test.statistic= "F")
ldmc_r2 <- r.squaredGLMM(ldmc_1)
ldmc_anova
ldmc_r2

summary(flowerstatus_1)
car::Anova(flowerstatus_1, type = "III", test.statistic = "Chisq", component="cond")
car::Anova(flowerstatus_1, type = "III", test.statistic = "Chisq", component="zi")
performance::r2(flowerstatus_1)

summary(flowerstructure_1)
car::Anova(flowerstructure_1, type = "III", test.statistic = "Chisq", component="cond")
car::Anova(flowerstructure_1, type = "III", test.statistic = "Chisq", component="zi")
performance::r2(flowerstructure_1)

summary(flowerdays_1)
flowdays_anova <- Anova(flowerdays_1, type = "III", test.statistic = "Chisq")
flowdays_r2 <- r.squaredGLMM(flowerdays_1)
flowdays_anova
flowdays_r2

summary(seedmass_1)
seedmass_anova <- Anova(seedmass_1, type = "III", test.statistic = "Chisq")
seedmass_r2 <- r.squaredGLMM(seedmass_1)
seedmass_anova
seedmass_r2

summary(seednum_1)
seednum_anova <- Anova(seednum_1, type = "III", test.statistic = "Chisq")
seednum_r2 <- r.squaredGLMM(seednum_1)
seednum_anova
seednum_r2