# title: 02_models

# about: run initial models as proposed in manuscript.
# many times they don't work, so there are many models here and it is awful
# check for normality, heteroscedacsity, and outliers on each piece of data
# generally figure out wtf to do, lots of level fixed effects rather than continuous means my models are Hard

# author: madeleine wallace

# PACKAGES -----------------------------------------------------------------
library(tidyverse)
library(lme4) #model
library(lmerTest) #model
library(nlme) #model with variance functions if i need it
library(DHARMa) #model diagnostics
library(performance) #model diagnostics
library(MuMIn) #r2
library(parameters) # fix heteroscedascity?
library(car) #anove for glmers
source('code-and-data/scripts/plpa-funeco-paper/01_clean.R')

# trait ~ parental treatment + 
        # offspring treatment + 
        # population + 
        # PT x OT + 
        # PT x population + 
        # OT x population + 
        # PT x OT x population +
        # (1|tray)



# BIOMASS MODELS ----------------------------------------------------------

### root 7/25
biomassroot_model <- lmer(root ~ `2023_treat` + `2021_treat` + pop + `2023_treat`:`2021_treat` + pop:`2023_treat` + pop:`2021_treat` + `2021_treat`:`2023_treat`:pop + (1|tray), data=BIOMASS_FINAL, na.action = na.omit)



### root - good
# root with nlme
biomassroot_model_nlme <- lme(root ~ ot + pt + pop + ot:pt + pop:ot + pop:pt + pt:ot:pop,
                              random = ~ 1 | tray,
                              data = BIOMASS_FINAL,
                              weights = varIdent(form = ~ 1 | pop),  # different variances for different levels of factors
                              na.action = na.omit)

biomassroot_model_nlme_2 <- lme(log_root ~ ot + pt + pop + ot:pt + pop:ot + pop:pt + pt:ot:pop,
                              random = ~ 1 | tray,
                              data = BIOMASS_FINAL,
                              weights = varIdent(form = ~ 1 | pop),  # different variances for different levels of factors
                              na.action = na.omit)

qqnorm(resid(biomassroot_model_nlme))
qqline(resid(biomassroot_model_nlme))
plot(biomassroot_model_nlme)

# not nlme
biomassroot_model <- lmer(root ~ `2023_treat` + `2021_treat` + pop + `2023_treat`:`2021_treat` + pop:`2023_treat` + pop:`2021_treat` + `2021_treat`:`2023_treat`:pop + (1|tray), data=BIOMASS_FINAL, na.action = na.omit)
biomassroot_model_2 <- lmer(log_root ~ `2023_treat` + `2021_treat` + pop + `2023_treat`:`2021_treat` + pop:`2023_treat` + pop:`2021_treat` + `2021_treat`:`2023_treat`:pop + (1|tray), data=BIOMASS_FINAL, na.action = na.omit)

qqnorm(resid(biomassroot_model))
qqline(resid(biomassroot_model))
plot(biomassroot_model)
simulateResiduals(biomassroot_model, plot = TRUE)
check_model(biomassroot_model)
check_heteroscedasticity(biomassroot_model)

# root - data is heteroskedastic, but logging makes it worse
# and using nlme with variance parameters didn't seem to have much of a difference in final models
# i think i can just ignore the heteroskedascity
# also there's a paper where it says lmer is robust enough to deal with heteroskedascitiy esp when the slopes and intercepts aren't fixed
# i think the thing that matters the most here is that the residuals are normally distributed - and the unlogged model with lmer looks the best

# final stats
summary(biomassroot_model)
Anova(biomassroot_model, type="III",  test.statistic= "F", ddf="Kenward-Roger")
r.squaredGLMM(biomassroot_model)


### shoot
biomassshoot_model <- lmer(shoot ~ `2023_treat` + `2021_treat` + pop + `2023_treat`:`2021_treat` + pop:`2023_treat` + pop:`2021_treat` + `2021_treat`:`2023_treat`:pop + (1|tray), data=BIOMASS_FINAL, na.action = na.omit)

biomassshoot_model_2 <- lmer(log_shoot ~ `2023_treat` + `2021_treat` + pop + `2023_treat`:`2021_treat` + pop:`2023_treat` + pop:`2021_treat` + `2021_treat`:`2023_treat`:pop + (1|tray), data=BIOMASS_FINAL, na.action = na.omit)

biomassshoot_model_nlme <- lme(shoot ~ ot + pt + pop + ot:pt + pop:ot + pop:pt + pt:ot:pop,
                               random = ~ 1 | tray,
                               data = BIOMASS_FINAL,
                               weights = varIdent(form = ~ 1 | pop),  # different variances for different levels of factors
                               na.action = na.omit)

# check
simulateResiduals(biomassshoot_model, plot = TRUE)
check_model(biomassshoot_model)
check_heteroscedasticity(biomassshoot_model)
plot(biomassshoot_model)
qqnorm(resid(biomassshoot_model))
qqline(resid(biomassshoot_model))

simulateResiduals(biomassshoot_model_2, plot = TRUE)
check_model(biomassshoot_model_2)
qqnorm(resid(biomassshoot_model_2))
qqline(resid(biomassshoot_model_2))
plot(biomassshoot_model_2)

check_model(biomassshoot_model_nlme)
qqnorm(resid(biomassshoot_model_nlme))
qqline(resid(biomassshoot_model_nlme))
plot(biomassshoot_model_nlme)

# again, just looking at the model assumptions and making sure they fit and ignoring heteroscedasticity, 
# the lmer model with no log seems to have the best normality and linearity of residuals

# final stats
summary(biomassshoot_model)
Anova(biomassshoot_model, type="III",  test.statistic= "F", ddf="Kenward-Roger")
r.squaredGLMM(biomassshoot_model)


### total biomass
biomasstotal_model <- lmer(total_biomass ~ `2023_treat` + `2021_treat` + pop + `2023_treat`:`2021_treat` + pop:`2023_treat` + pop:`2021_treat` + `2021_treat`:`2023_treat`:pop + (1|tray), data=BIOMASS_FINAL, na.action = na.omit)
biomasstotal_model_2 <- lmer(log_total ~ `2023_treat` + `2021_treat` + pop + `2023_treat`:`2021_treat` + pop:`2023_treat` + pop:`2021_treat` + `2021_treat`:`2023_treat`:pop + (1|tray), data=BIOMASS_FINAL, na.action = na.omit)


# check
plot(biomasstotal_model) # slight heteroscedascity here
plot(biomasstotal_model_2)

simulateResiduals(biomasstotal_model, plot = TRUE) # better here
simulateResiduals(biomasstotal_model_2, plot = TRUE)

check_model(biomasstotal_model) # linearity and normality of residuals is better here but homoegeity of variance is worse
check_model(biomasstotal_model_2) 

qqnorm(resid(biomasstotal_model))
qqline(resid(biomasstotal_model))

# heteroscedasticity detected - logging made it better
# but made the residuals and linearity worse, which i think matters more than the two things above if i am following my rule

# stats
summary(biomasstotal_model)
Anova(biomasstotal_model, type="III",  test.statistic= "F", ddf="Kenward-Roger")
r.squaredGLMM(biomasstotal_model)


### RS ratio - used the transformed values (logged)!
biomassratio_model <- lmer(ratio_RS ~ `2023_treat` + `2021_treat` + pop + `2023_treat`:`2021_treat` + pop:`2023_treat` + pop:`2021_treat` + `2021_treat`:`2023_treat`:pop + (1|tray), data=BIOMASS_FINAL, na.action = na.omit)
biomassratio_model_2 <- lmer(log_RS ~ `2023_treat` + `2021_treat` + pop + `2023_treat`:`2021_treat` + pop:`2023_treat` + pop:`2021_treat` + `2021_treat`:`2023_treat`:pop + (1|tray), data=BIOMASS_FINAL, na.action = na.omit)


# check
plot(biomassratio_model)
plot(biomassratio_model_2) #WAYYYYY better

simulateResiduals(biomassratio_model, plot = TRUE) 
simulateResiduals(biomassratio_model_2, plot = TRUE) #WAYYYY better

check_model(biomassratio_model) 
check_model(biomassratio_model_2) # literally the best

# stats
summary(biomassratio_model_2)
Anova(biomassratio_model_2, type="III",  test.statistic= "F", ddf="Kenward-Roger")
r.squaredGLMM(biomassratio_model_2)




# FLOWER MODELS -----------------------------------------------------------
# i think i need to include an analysis of what plants flowered and what plants did not flower!!
# oh lord this one weird


# number of structures
flowerstructure_model <- glmer(num_structure ~ `2023_treat` + 
                                 `2021_treat` + 
                                 pop + 
                                 `2023_treat`:`2021_treat` + 
                                 pop:`2023_treat` + 
                                 pop:`2021_treat` + 
                                 `2021_treat`:`2023_treat`:pop + 
                                 (1|tray), family = poisson, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)), na.action = na.omit, data=FLOWER_FINAL)
# check
simulateResiduals(flowerstructure_model, plot = TRUE)
check_model(flowerstructure_model)

plot(flowerstructure_model)
qqnorm(resid(flowerstructure_model))
qqline(resid(flowerstructure_model))

# stats
summary(flowerstructure_model)
car::Anova(flowerstructure_model, type=3)
r.squaredGLMM(flowerstructure_model) 


# days to flower
hist(log(FLOWER_FINAL$days_to_flower))
flowerdays_model <- lmer(log(days_to_flower) ~ `2023_treat` + `2021_treat` + pop + `2023_treat`:`2021_treat` + pop:`2023_treat` + pop:`2021_treat` + `2021_treat`:`2023_treat`:pop + (1|tray), data=FLOWER_FINAL, na.action = na.omit)

simulateResiduals(flowerdays_model, plot = TRUE)
check_model(flowerdays_model)
check_heteroscedasticity(flowerdays_model)

summary(flowerdays_model)
Anova(flowerdays_model, type="III",  test.statistic= "F", ddf="Kenward-Roger")
r.squaredGLMM(flowerdays_model)




# GROWTH MODELS -----------------------------------------------------------


# max height
hist(HEIGHT_FINAL$max)
maxheight_model <- lmer(max ~ `2023_treat` + `2021_treat` + pop + `2023_treat`:`2021_treat` + pop:`2023_treat` + pop:`2021_treat` + `2021_treat`:`2023_treat`:pop + (1|tray), data=HEIGHT_FINAL, na.action = na.omit)

simulateResiduals(maxheight_model, plot = TRUE)
check_model(maxheight_model)
check_heteroscedasticity(maxheight_model)

summary(maxheight_model)
Anova(maxheight_model, type="III",  test.statistic= "F", ddf="Kenward-Roger")
r.squaredGLMM(maxheight_model)


# RGR
hist(sqrt(RGR_FINAL$RGR))
RGR_model <- lmer(RGR ~ `2023_treat` + `2021_treat` + pop + `2023_treat`:`2021_treat` + pop:`2023_treat` + pop:`2021_treat` + `2021_treat`:`2023_treat`:pop + (1|tray), data=RGR_FINAL, na.action = na.omit)

simulateResiduals(RGR_model, plot = TRUE)
check_model(RGR_model)
check_heteroscedasticity(RGR_model)

qqnorm(resid(RGR_model))

summary(RGR_model)
Anova(RGR_model, type="III",  test.statistic= "F", ddf="Kenward-Roger")
r.squaredGLMM(RGR_model)




# MORTALITY MODEL --------------------------------------------------------
# this is binary data! 0 or 1 on july 13th
# removed the random effect here, overspecification of random effects (Matuschek et al. 2017)

mort_model <- glm(status ~ `2023_treat` + 
                    `2021_treat` + 
                    pop + 
                    `2023_treat`:`2021_treat` + 
                    pop:`2023_treat` + 
                    pop:`2021_treat` + 
                    `2021_treat`:`2023_treat`:pop, 
                  data=MORT_JULY13_FINAL, family = binomial(link = "logit"), na.action = na.omit)

simulateResiduals(mort_model, plot = TRUE)
check_model(mort_model)

summary(mort_model)
Anova(mort_model, type="III",  test.statistic= "F", ddf="Kenward-Roger")
r.squaredGLMM(mort_model)




# SLA AND LDMC MODEL ------------------------------------------------------


# sla
hist(sqrt(SLA_LDMC_FINAL$sla))
sla_model <- lmer(sla ~ `2023_treat` + `2021_treat` + pop + `2023_treat`:`2021_treat` + pop:`2023_treat` + pop:`2021_treat` + `2021_treat`:`2023_treat`:pop + (1|tray), data=SLA_LDMC_FINAL, na.action = na.omit)

qqnorm(resid(sla_model))
qqline(resid(sla_model))

qqnorm(resid(sla_model_2))
qqline(resid(sla_model_2))

simulateResiduals(sla_model, plot = TRUE)
check_model(sla_model)
check_heteroscedasticity(sla_model)

summary(sla_model)
Anova(sla_model, type="III",  test.statistic= "F", ddf="Kenward-Roger")
r.squaredGLMM(sla_model)


# ldmc
hist(log(SLA_LDMC_FINAL$ldmc))
ldmc_model <- lmer(ldmc ~ `2023_treat` + `2021_treat` + pop + `2023_treat`:`2021_treat` + pop:`2023_treat` + pop:`2021_treat` + `2021_treat`:`2023_treat`:pop + (1|tray), data=SLA_LDMC_FINAL, na.action = na.omit)
ldmc_model2 <- lmer(log(ldmc) ~ `2023_treat` + `2021_treat` + pop + `2023_treat`:`2021_treat` + pop:`2023_treat` + pop:`2021_treat` + `2021_treat`:`2023_treat`:pop + (1|tray), data=SLA_LDMC_FINAL, na.action = na.omit)


simulateResiduals(ldmc_model, plot = TRUE)
simulateResiduals(ldmc_model2, plot = TRUE)
check_model(ldmc_model)
check_model(ldmc_model2)
check_heteroscedasticity(ldmc_model)

summary(ldmc_model2)
Anova(ldmc_model2, type="III",  test.statistic= "F", ddf="Kenward-Roger")
r.squaredGLMM(ldmc_model2)




# SEED MODELS -------------------------------------------------------------


# mass
hist(SEED_FINAL$`2023_mass_total_g`)
seedmass_model <- lmer(`2023_mass_total_g` ~ `2023_treat` + `2021_treat` + pop + `2023_treat`:`2021_treat` + pop:`2023_treat` + pop:`2021_treat` + `2021_treat`:`2023_treat`:pop + (1|tray), data=SEED_FINAL, na.action = na.omit)
seedmass_model2 <- lmer(sqrt(`2023_mass_total_g`) ~ `2023_treat` + `2021_treat` + pop + `2023_treat`:`2021_treat` + pop:`2023_treat` + pop:`2021_treat` + `2021_treat`:`2023_treat`:pop + (1|tray), data=SEED_FINAL, na.action = na.omit)


simulateResiduals(seedmass_model, plot = TRUE)
simulateResiduals(seedmass_model2, plot = TRUE)
check_model(seedmass_model)
check_model(seedmass_model2)

qqnorm(resid(seedmass_model))
qqline(resid(seedmass_model))

qqnorm(resid(seedmass_model2))
qqline(resid(seedmass_model2))

summary(seedmass_model2)
Anova(seedmass_model2, type="III",  test.statistic= "F", ddf="Kenward-Roger")
r.squaredGLMM(seedmass_model2)


# number
hist((SEED_FINAL$`2023_num_total`))
seednum_model <- lmer(`2023_num_total` ~ `2023_treat` + `2021_treat` + pop + `2023_treat`:`2021_treat` + pop:`2023_treat` + pop:`2021_treat` + `2021_treat`:`2023_treat`:pop + (1|tray), data=SEED_FINAL, na.action = na.omit)
seednum_model2 <- lmer(sqrt(`2023_num_total`) ~ `2023_treat` + `2021_treat` + pop + `2023_treat`:`2021_treat` + pop:`2023_treat` + pop:`2021_treat` + `2021_treat`:`2023_treat`:pop + (1|tray), data=SEED_FINAL, na.action = na.omit)


simulateResiduals(seednum_model, plot = TRUE)
simulateResiduals(seednum_model2, plot = TRUE)
check_model(seednum_model)
check_model(seednum_model2)
check_heteroscedasticity(seednum_model)

summary(seednum_model2)
Anova(seednum_model2, type="III",  test.statistic= "F", ddf="Kenward-Roger")
r.squaredGLMM(seedmass_model2)


#### FINAL MODELS..... I HOPE

# ROOT
#straight up lmer, random effect, no variance
biomassroot_model1 <- lmer(root ~ ot + pt + pop + ot:pt + pop:ot + pop:pt + pt:ot:pop + (1|tray), 
                           data=BIOMASS_FINAL, 
                           na.action = na.exclude)

# lme, random effect plus variance pop
biomassroot_model2 <- lme(root ~ ot + pt + pop + ot:pt + pop:ot + pop:pt + pt:ot:pop,
                          random = ~ 1 | tray,
                          data = BIOMASS_FINAL,
                          weights = varIdent(form = ~ 1 | pop),  # different variances for different levels of factors
                          na.action = na.exclude)

# lme, random effect plus variance ot
biomassroot_model3 <- lme(root ~ ot + pt + pop + ot:pt + pop:ot + pop:pt + pt:ot:pop,
                          random = ~ 1 | tray,
                          data = BIOMASS_FINAL,
                          weights = varIdent(form = ~ 1 | ot),  # different variances for different levels of factors
                          na.action = na.exclude)

# gls, just variance pop, no random effect
biomassroot_model4 <- gls(root ~ ot + pt + pop + ot:pt + pop:ot + pop:pt + pt:ot:pop,
                          data = BIOMASS_FINAL,
                          weights = varIdent(form = ~ 1 | pop),  # different variances for different levels of factors
                          na.action = na.exclude)

# gls, just variance ot, no random effect
biomassroot_model5 <- gls(root ~ ot + pt + pop + ot:pt + pop:ot + pop:pt + pt:ot:pop,
                          data = BIOMASS_FINAL,
                          weights = varIdent(form = ~ 1 | ot),  # different variances for different levels of factors
                          na.action = na.exclude)

AIC(biomassroot_model1, biomassroot_model2, biomassroot_model3, biomassroot_model4, biomassroot_model5)

# model 4 has lowest AIC

# SHOOT
biomassshoot_model1 <- lmer(shoot ~ ot + pt + pop + ot:pt + pop:ot + pop:pt + pt:ot:pop + (1|tray), 
                           data=BIOMASS_FINAL, 
                           na.action = na.exclude)

# lme, random effect plus variance pop
biomassshoot_model2 <- lme(shoot ~ ot + pt + pop + ot:pt + pop:ot + pop:pt + pt:ot:pop,
                          random = ~ 1 | tray,
                          data = BIOMASS_FINAL,
                          weights = varIdent(form = ~ 1 | pop),  # different variances for different levels of factors
                          na.action = na.exclude)

# lme, random effect plus variance ot
biomassshoot_model3 <- lme(shoot ~ ot + pt + pop + ot:pt + pop:ot + pop:pt + pt:ot:pop,
                          random = ~ 1 | tray,
                          data = BIOMASS_FINAL,
                          weights = varIdent(form = ~ 1 | ot),  # different variances for different levels of factors
                          na.action = na.exclude)

# gls, just variance pop, no random effect
biomassshoot_model4 <- gls(shoot ~ ot + pt + pop + ot:pt + pop:ot + pop:pt + pt:ot:pop,
                          data = BIOMASS_FINAL,
                          weights = varIdent(form = ~ 1 | pop),  # different variances for different levels of factors
                          na.action = na.exclude)

# gls, just variance ot, no random effect
biomassshoot_model5 <- gls(shoot ~ ot + pt + pop + ot:pt + pop:ot + pop:pt + pt:ot:pop,
                          data = BIOMASS_FINAL,
                          weights = varIdent(form = ~ 1 | ot),  # different variances for different levels of factors
                          na.action = na.exclude)

AIC(biomassshoot_model1, biomassshoot_model2, biomassshoot_model3, biomassshoot_model4, biomassshoot_model5)

# model 4 has lowest AIC

# TOTAL
biomasstotal_model1 <- lmer(total_biomass ~ ot + pt + pop + ot:pt + pop:ot + pop:pt + pt:ot:pop + (1|tray), 
                            data=BIOMASS_FINAL, 
                            na.action = na.exclude)

# lme, random effect plus variance pop
biomasstotal_model2 <- lme(total_biomass ~ ot + pt + pop + ot:pt + pop:ot + pop:pt + pt:ot:pop,
                           random = ~ 1 | tray,
                           data = BIOMASS_FINAL,
                           weights = varIdent(form = ~ 1 | pop),  # different variances for different levels of factors
                           na.action = na.exclude)

# lme, random effect plus variance ot
biomasstotal_model3 <- lme(total_biomass ~ ot + pt + pop + ot:pt + pop:ot + pop:pt + pt:ot:pop,
                           random = ~ 1 | tray,
                           data = BIOMASS_FINAL,
                           weights = varIdent(form = ~ 1 | ot),  # different variances for different levels of factors
                           na.action = na.exclude)

# gls, just variance pop, no random effect
biomasstotal_model4 <- gls(total_biomass ~ ot + pt + pop + ot:pt + pop:ot + pop:pt + pt:ot:pop,
                           data = BIOMASS_FINAL,
                           weights = varIdent(form = ~ 1 | pop),  # different variances for different levels of factors
                           na.action = na.exclude)

# gls, just variance ot, no random effect
biomasstotal_model5 <- gls(total_biomass ~ ot + pt + pop + ot:pt + pop:ot + pop:pt + pt:ot:pop,
                           data = BIOMASS_FINAL,
                           weights = varIdent(form = ~ 1 | ot),  # different variances for different levels of factors
                           na.action = na.exclude)

AIC(biomasstotal_model1, biomasstotal_model2, biomasstotal_model3, biomasstotal_model4, biomasstotal_model5)

# model 3 has lowest AIC


# RS RATIO
biomassRS_model1 <- lmer(log_RS ~ ot + pt + pop + ot:pt + pop:ot + pop:pt + pt:ot:pop + (1|tray), 
                            data=BIOMASS_FINAL, 
                            na.action = na.exclude)

# lme, random effect plus variance pop
biomassRS_model2 <- lme(log_RS ~ ot + pt + pop + ot:pt + pop:ot + pop:pt + pt:ot:pop,
                           random = ~ 1 | tray,
                           data = BIOMASS_FINAL,
                           weights = varIdent(form = ~ 1 | pop),  # different variances for different levels of factors
                           na.action = na.exclude)

# lme, random effect plus variance ot
biomassRS_model3 <- lme(log_RS ~ ot + pt + pop + ot:pt + pop:ot + pop:pt + pt:ot:pop,
                           random = ~ 1 | tray,
                           data = BIOMASS_FINAL,
                           weights = varIdent(form = ~ 1 | ot),  # different variances for different levels of factors
                           na.action = na.exclude)

# gls, just variance pop, no random effect
biomassRS_model4 <- gls(log_RS ~ ot + pt + pop + ot:pt + pop:ot + pop:pt + pt:ot:pop,
                           data = BIOMASS_FINAL,
                           weights = varIdent(form = ~ 1 | pop),  # different variances for different levels of factors
                           na.action = na.exclude)

# gls, just variance ot, no random effect
biomassRS_model5 <- gls(log_RS ~ ot + pt + pop + ot:pt + pop:ot + pop:pt + pt:ot:pop,
                           data = BIOMASS_FINAL,
                           weights = varIdent(form = ~ 1 | ot),  # different variances for different levels of factors
                           na.action = na.exclude)

AIC(biomassRS_model1, biomassRS_model2, biomassRS_model3, biomassRS_model4, biomassRS_model5)

# model 2 has lowest AIC