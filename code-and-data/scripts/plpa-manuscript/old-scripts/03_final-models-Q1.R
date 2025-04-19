# title: 03_final-models

# about: final models chosen for each trait response. after many years

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
library(pscl)

# PROBLEM SOLVING ---------------------------------------------------------

# i think the problem is that populations 3, 5, and 11 had very few flowering individuals when the offspring treatment was drought
# so if i just remove pop 3, 5, and 11 from the dataset before i fill in the model, maybve that would be better?

# for all flowering data: iu removed population 3, 5, and 11 because of missing factor combinations, because of low flowering rates across populations
# for all flowering data: i removed the interaction pop:pt:ot because of missing factor combinations, because of low flowering rates across populations
# i left in random effect for flowering data

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

# root biomass - medium amount of heterogeneity ---------------------------
biomassroot_model <- lme(root ~ ot + pt + pop + ot:pt + pop:ot + pop:pt + pt:ot:pop,
                          random = ~ 1 | tray,
                          data = BIOMASS_FINAL,
                          weights = varIdent(form = ~ 1 | pop),
                          na.action = na.exclude)

# check
plot(biomassroot_model, add.smooth = FALSE, which = 1)
E <- resid(biomassroot_model)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(biomassroot_model))
qqline(resid(biomassroot_model))
plot(BIOMASS_FINAL$pop, resid(biomassroot_model), xlab = "pop", ylab = "residuals")
plot(BIOMASS_FINAL$`2023_treat`, resid(biomassroot_model), xlab = "2023 treat", ylab = "residuals")

# final stats
summary(biomassroot_model)
root_anova <- anova(biomassroot_model, type='marginal')
r.squaredGLMM(biomassroot_model)

root_anova$Significance <- sapply(root_anova$`p-value`, add_stars)
print(root_anova)



# ## shoot biomass - lots of heterogeneity! -------------------------------
biomassshoot_model <- lme(shoot ~ ot + pt + pop + ot:pt + pop:ot + pop:pt + pt:ot:pop,
                           random = ~ 1 | tray,
                           data = BIOMASS_FINAL,
                           weights = varIdent(form = ~ 1 | ot),
                           na.action = na.exclude)


# check
plot(biomassshoot_model, add.smooth = FALSE, which = 1)
E <- resid(biomassshoot_model)
hist(E, xlab = "residuals", main = " ") #normal!
qqnorm(resid(biomassshoot_model))
qqline(resid(biomassshoot_model))
plot(BIOMASS_FINAL$pop, resid(biomassshoot_model), xlab = "pop", ylab = "residuals")
plot(BIOMASS_FINAL$`2023_treat`, resid(biomassshoot_model), xlab = "2023 treat", ylab = "residuals")


# final stats
summary(biomassshoot_model)
shoot_anova <- anova(biomassshoot_model, type='marginal')
r.squaredGLMM(biomassshoot_model)

shoot_anova$Significance <- sapply(shoot_anova$`p-value`, add_stars)
print(shoot_anova)

#contrasts
shootcontrast <- emmeans(biomassshoot_model, specs = pairwise ~ ot:pt)
shootcontrast$contrasts



# total biomass -----------------------------------------------------------
biomasstotal_model <- lme(total_biomass ~ ot + pt + pop + ot:pt + pop:ot + pop:pt + pt:ot:pop,
                           random = ~ 1 | tray,
                           data = BIOMASS_FINAL,
                           weights = varIdent(form = ~ 1 | pop),
                           na.action = na.exclude)

# check
plot(biomasstotal_model, add.smooth = FALSE, which = 1)
E <- resid(biomasstotal_model)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(biomasstotal_model))
qqline(resid(biomasstotal_model))
plot(BIOMASS_FINAL$pop, resid(biomasstotal_model), xlab = "pop", ylab = "residuals")
plot(BIOMASS_FINAL$`2023_treat`, resid(biomasstotal_model), xlab = "2023 treat", ylab = "residuals")

# final stats
summary(biomasstotal_model)
total_anova <- anova(biomasstotal_model, type='marginal')
r.squaredGLMM(biomasstotal_model)

total_anova$Significance <- sapply(total_anova$`p-value`, add_stars)
print(total_anova)

#contrasts
totalcontrast <- emmeans(biomasstotal_model, specs = pairwise ~ ot:pt)
totalcontrast$contrasts



# ratio RS ----------------------------------------------------------------
biomassRS_model1 <- lme(log_RS ~ ot + pt + pop + ot:pt + pop:ot + pop:pt + pt:ot:pop,
                        random = ~ 1 | tray,
                        data = BIOMASS_FINAL,
                        weights = varIdent(form = ~ 1 | pop), 
                        na.action = na.exclude)

# check
plot(biomassRS_model1, add.smooth = FALSE, which = 1)
E <- resid(biomassRS_model1)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(biomassRS_model1))
qqline(resid(biomassRS_model1))
plot(BIOMASS_FINAL$pop, resid(biomassRS_model1), xlab = "pop", ylab = "residuals")
plot(BIOMASS_FINAL$`2023_treat`, resid(biomassRS_model1), xlab = "2023 treat", ylab = "residuals")

# final stats
summary(biomassRS_model1)
rs_anova <- anova(biomassRS_model1, type='marginal')
r.squaredGLMM(biomassRS_model1)

rs_anova$Significance <- sapply(rs_anova$`p-value`, add_stars)
print(rs_anova)

#contrasts
rscontrast <- emmeans(biomassRS_model1, specs = pairwise ~ ot:pt, type = "response")
rscontrast$contrasts



# max height --------------------------------------------------------------
hist(HEIGHT_FINAL$max)
maxheight_model2 <- lme(max ~ ot + pt + pop + ot:pt + pop:ot + pop:pt + pt:ot:pop,
                          random = ~ 1 | tray,
                          data = HEIGHT_FINAL,
                          weights = varIdent(form = ~ 1 | pop),  # different variances for different levels of factors
                          na.action = na.exclude)

check_model(maxheight_model2)
plot(maxheight_model2, add.smooth = FALSE, which = 1)
E <- resid(maxheight_model3)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(maxheight_model3))
qqline(resid(maxheight_model3))
plot(HEIGHT_FINAL$pop, resid(maxheight_model3), xlab = "pop", ylab = "residuals")
plot(HEIGHT_FINAL$`2023_treat`, resid(maxheight_model3), xlab = "2023 treat", ylab = "residuals")

# final stats
summary(maxheight_model2)
maxheight_anova <- anova(maxheight_model2, type='marginal')
anova(maxheight_model2, type='marginal')
r.squaredGLMM(maxheight_model2)

maxheight_anova$Significance <- sapply(maxheight_anova$`p-value`, add_stars)
print(maxheight_anova)



# RGR ---------------------------------------------------------------------
head(RGR_FINAL)
RGR_model2 <- lme(RGR ~ ot + pt + pop + ot:pt + pop:ot + pop:pt + pt:ot:pop,
                        random = ~ 1 | tray,
                        data = RGR_FINAL,
                        weights = varIdent(form = ~ 1 | pop),  # different variances for different levels of factors
                        na.action = na.exclude)

check_model(RGR_model2)
plot(RGR_model3, add.smooth = FALSE, which = 1)
E <- resid(RGR_model3)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(RGR_model2))
qqline(resid(RGR_model2))
plot(RGR_FINAL$pop, resid(RGR_model2), xlab = "pop", ylab = "residuals")
plot(RGR_FINAL$`2023_treat`, resid(RGR_model3), xlab = "2023 treat", ylab = "residuals")

# final stats
summary(RGR_model2)
rgr_anova <- anova(RGR_model2, type='marginal')
anova(RGR_model2, type='marginal')
r.squaredGLMM(RGR_model2)

rgr_anova$Significance <- sapply(rgr_anova$`p-value`, add_stars)
print(rgr_anova)



# mortality --------------------------------------------------------------------
# i removed random effect because there wasn't enough variance in it
head(mort_day50)
mort_model <- glm(status ~ ot + pt + pop + ot:pt + pop:ot + pop:pt + pt:ot:pop,
                  data=mort_day50, family = binomial(link = "logit"), na.action = na.exclude)

check_model(mort_model)
plot(mort_model, add.smooth = FALSE, which = 1)
E <- resid(mort_model)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(mort_model))
qqline(resid(mort_model))
plot(mort_day50$pop, resid(mort_model), xlab = "pop", ylab = "residuals")
plot(mort_day50$`2023_treat`, resid(mort_model), xlab = "2023 treat", ylab = "residuals")

# final stats
summary(mort_model)
1-(mort_model$deviance/mort_model$null.deviance) #PseudoR^2 (MacFayden)
Anova(mort_model, type = 3, test.statistic = "F")



# sla ---------------------------------------------------------------------
## sla - no clear variance pattern between populations, so i used lmer and logged SLA
# to help with heteroscedascity
hist(SLA_LDMC_FINAL$sla)
sla_model <- lmer(log(sla) ~ ot + pt + pop + ot:pt + pop:ot + pop:pt + pt:ot:pop + (1|tray), 
                  data=SLA_LDMC_FINAL, na.action = na.exclude)

check_model(sla_model)
plot(sla_model, add.smooth = FALSE, which = 1)
E <- resid(sla_model)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(sla_model))
qqline(resid(sla_model))
plot(SLA_LDMC_FINAL$pop, resid(sla_model), xlab = "pop", ylab = "residuals")
plot(SLA_LDMC_FINAL$`2023_treat`, resid(sla_model), xlab = "2023 treat", ylab = "residuals")

# final stats
summary(sla_model)
sla_anova <- Anova(sla_model, type="III",  test.statistic= "F")
Anova(sla_model, type="III",  test.statistic= "F")
r.squaredGLMM(sla_model)

sla_anova$Significance <- sapply(sla_anova$`p-value`, add_stars)
print(sla_anova)

# contrasts
slacontrast <- emmeans(sla_model, specs = pairwise ~ ot:pt, type = "response")
slacontrast$contrasts



# ldmc --------------------------------------------------------------------
## ldmc- no clear variance pattern between populations, so i used lmer and logged SLA
# to help with heteroscedascity
ldmc_model2 <- lmer(log(ldmc) ~ ot + pt + pop + ot:pt + pop:ot + pop:pt + pt:ot:pop + 
                      (1|tray), data=SLA_LDMC_FINAL, na.action = na.exclude)

check_model(ldmc_model2)
plot(ldmc_model2, add.smooth = FALSE, which = 1)
E <- resid(ldmc_model2)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(ldmc_model2))
qqline(resid(ldmc_model2))
plot(SLA_LDMC_FINAL$pop, resid(ldmc_model2), xlab = "pop", ylab = "residuals")
plot(SLA_LDMC_FINAL$`2023_treat`, resid(ldmc_model2), xlab = "2023 treat", ylab = "residuals")

# final stats
summary(ldmc_model2)
ldmc_anova <- Anova(ldmc_model2, type="III",  test.statistic= "F")
Anova(ldmc_model2, type="III",  test.statistic= "F")
r.squaredGLMM(ldmc_model2)

ldmc_anova$Significance <- sapply(ldmc_anova$`p-value`, add_stars)
print(ldmc_anova)

# contrasts
ldmccontrast <- emmeans(ldmc_model2, specs = pairwise ~ ot:pt, type = "response")
ldmccontrast$contrasts



# number of structures ----------------------------------------------------
#filtered_flower_data <- FLOWER_FINAL[!FLOWER_FINAL$pop %in% c(3, 5, 11), ]
head(FLOWER_FINAL)
mean(FLOWER_FINAL$num_structure, na.rm=T)
var(FLOWER_FINAL$num_structure, na.rm=T)
hist(FLOWER_FINAL$num_structure)

flowerstructure_model3 <- glmer(num_structure ~ ot + pt + pop + ot:pt + pop:ot + pop:pt + pt:ot:pop + (1|tray),
                                family = "poisson",
                                control = glmerControl(optimizer ="bobyqa", optCtrl=list(maxfun=2000000)),
                                na.action = na.exclude, data=FLOWER_FINAL)

# check
check_model(flowerstructure_model3)
plot(flowerstructure_model3, add.smooth = FALSE, which = 1)
E <- resid(flowerstructure_model3)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(flowerstructure_model3))
qqline(resid(flowerstructure_model3))
plot(filtered_flower_data$pop, resid(flowerstructure_model3), xlab = "pop", ylab = "residuals")
plot(filtered_flower_data$`2023_treat`, resid(flowerstructure_model3), xlab = "2023 treat", ylab = "residuals")


# final stats
summary(flowerstructure_model3)
r.squaredGLMM(flowerstructure_model3)
Anova(flowerstructure_model3, type = 3)



# days to flower ----------------------------------------------------------
hist(FLOWER_FINAL$days_to_flower)
mean(FLOWER_FINAL$days_to_flower, na.rm=T)
var(FLOWER_FINAL$days_to_flower, na.rm=T) #overdispersion i think!
filtered_flower_data <- filtered_flower_data |> 
  mutate(days_to_flower = as.numeric(days_to_flower))
head(filtered_flower_data)

flowerdays_model4 <- glmer.nb(days_to_flower ~ ot + pt + pop + ot:pt + pop:ot + pop:pt + pop:pt:ot + (1|tray),
                              control = glmerControl(optimizer ="bobyqa", optCtrl=list(maxfun=2000000)),
                              na.action = na.exclude, data=FLOWER_FINAL)

# check
check_model(flowerdays_model5)
plot(flowerdays_model4, add.smooth = FALSE, which = 1) 
E <- resid(flowerdays_model4)
hist(E, xlab = "residuals", main = " ") 
qqnorm(resid(flowerdays_model4))
qqline(resid(flowerdays_model4))
plot(FLOWER_FINAL$pop, resid(flowerdays_model4), xlab = "pop", ylab = "residuals")
plot(FLOWER_FINAL$`2023_treat`, resid(flowerdays_model4), xlab = "2023 treat", ylab = "residuals")


# final stats
summary(flowerdays_model4)
r.squaredGLMM(flowerdays_model4)
Anova(flowerdays_model4, type = 3)



# seed mass ---------------------------------------------------------------
## seed mass- no clear variance pattern between populations, so i used lmer and sqrted total mass
# to help with heteroscedascity
hist(SEED_FINAL$mass_total)
head(SEED_FINAL)
seedmass_model1 <- lmer(sqrt(`2023_mass_total_g`) ~ ot + pt + pop + ot:pt + pop:ot + pop:pt + pt:ot:pop + (1|tray), 
                        data=SEED_FINAL, na.action = na.exclude)

check_model(seedmass_model1)
plot(seedmass_model1, add.smooth = FALSE, which = 1)
E <- resid(seedmass_model1)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(seedmass_model1))
qqline(resid(seedmass_model1))
plot(SEED_FINAL$pop, resid(seedmass_model1), xlab = "pop", ylab = "residuals")
plot(SEED_FINAL$`2023_treat`, resid(seedmass_model1), xlab = "2023 treat", ylab = "residuals")

# final stats
summary(seedmass_model1)
seedmass_anova <- Anova(seedmass_model1, type="III",  test.statistic= "F")
Anova(seedmass_model1, type="III",  test.statistic= "F")
r.squaredGLMM(seedmass_model1)

seedmass_anova$Significance <- sapply(seedmass_anova$`p-value`, add_stars)
print(seedmass_anova)




# seed num ----------------------------------------------------------------
hist(SEED_FINAL$num_total)
seednum_model1 <- glmer.nb(num_total ~ ot + pt + ot:pt + pop + pop:ot + pop:pt + pop:pt:ot + (1|tray),
                              control = glmerControl(optimizer ="bobyqa", optCtrl=list(maxfun=2000000)),
                              na.action = na.exclude, data=SEED_FINAL)

check_model(seednum_model1)
plot(seednum_model1, add.smooth = FALSE, which = 1)
E <- resid(seednum_model1)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(seednum_model1))
qqline(resid(seednum_model1))
plot(filtered_seed_data$pop, resid(seednum_model1), xlab = "pop", ylab = "residuals")
plot(filtered_seed_data$`2023_treat`, resid(seednum_model1), xlab = "2023 treat", ylab = "residuals")

# final stats
summary(seednum_model1)
r.squaredGLMM(seednum_model1)
Anova(seednum_model1, type = 3)

# contrasts
seednumcontrast <- emmeans(seednum_model1, specs = pairwise ~ ot:pt, type = "response")
seednumcontrast$contrasts

# number flowering ----------------------------------------------------------------
flower_status$ot <- factor(flower_status$ot)  # Ensure the correct order
flower_status$pt <- factor(flower_status$pt)  # Ensure the correct order
flower_status$status_flower <- numeric(flower_status$status_flower)

numflow_model <- glm(status_flower ~ ot + pt + pop + ot:pt + pop:ot + pop:pt + pt:ot:pop,
                       data = flower_status, 
                       family = binomial(link = "logit"), 
                       na.action = na.exclude)

zero_inflated_model <- zeroinfl(status_flower ~ ot + pt + pop + ot:pt + pop:ot + pop:pt + pt:ot:pop + (1|tray) |
                                  1,
                                data = flower_status, 
                                dist = "negbin")

mixed_model <- glmer(status_flower ~ ot + pt + pop + ot:pt + pop:ot + pop:pt + pt:ot:pop + 
                       (1|tray), 
                     data = flower_status, 
                     family = binomial(link = "logit"))

summary(zero_inflated_model)


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
1-(numflow_model$deviance/numflow_model$null.deviance) #PseudoR^2 (MacFayden)
Anova(numflow_model, type = 3, test.statistic = "F")

# neg binomial
check_model(zero_inflated_model)
plot(zero_inflated_model, add.smooth = FALSE, which = 1)
E <- resid(zero_inflated_model)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(zero_inflated_model))
qqline(resid(zero_inflated_model))
plot(flower_status$pop, resid(zero_inflated_model), xlab = "pop", ylab = "residuals")
plot(flower_status$`2023_treat`, resid(zero_inflated_model), xlab = "2023 treat", ylab = "residuals")

# final stats
summary(zero_inflated_model)
1-(numflow_model$deviance/numflow_model$null.deviance) #PseudoR^2 (MacFayden)
Anova(zero_inflated_model, type = 3, test.statistic = "F")


emm_results_zero_inflated <- emmeans(zero_inflated_model, ~ ot * pt)

# Perform pairwise comparisons
pairwise_results_zero_inflated <- pairs(emm_results_zero_inflated, adjust = "tukey")

# Display the results
summary(pairwise_results_zero_inflated)

## comparisons
# height
emmeans_results <- emmeans(maxheight_model2, ~ ot * pt)
print(emmeans_results)

pairwise_comparisons <- contrast(emmeans_results, method = "pairwise")
print(pairwise_comparisons)

emmeans_results_by_pop <- emmeans(maxheight_model2, ~ ot * pt | pop)
print(emmeans_results_by_pop)
contrast(emmeans_results_by_pop, method = "pairwise")

# mort
emmeans_results <- emmeans(mort_model, ~ ot * pt)
contrast(emmeans_results, method = "pairwise")

emmeans_results_by_pop <- emmeans(mort_model, ~ ot * pt | pop)
print(emmeans_results_by_pop)
contrast(emmeans_results_by_pop, method = "pairwise")



# summary of differences --------------------------------------------------


# structural models: biomassroot_model, biomassshoot_model, biomasstotal_model,
# biomassRS_model1, maxheight_model2, RGR_model2, sla_model, ldmc_model2

# survival: mort_model

# reproductive models: flowerstructure_model3, flowerdays_model4
# seedmass_model1, seednum_model1, numflow_model


# structural
emmeans_results <- emmeans(biomassroot_model, ~ ot * pt)
pairs(emmeans_results, adjust = "tukey")

emmeans_results <- emmeans(biomassshoot_model, ~ ot * pt)
pairs(emmeans_results, adjust = "tukey")

emmeans_results <- emmeans(biomasstotal_model, ~ ot * pt)
pairs(emmeans_results, adjust = "tukey")

emmeans_results <- emmeans(biomassRS_model1, ~ ot * pt)
pairs(emmeans_results, adjust = "tukey")

emmeans_results <- emmeans(maxheight_model2, ~ ot * pt)
contrast(emmeans_results, method = "pairwise")

emmeans_results <- emmeans(RGR_model2, ~ ot * pt)
contrast(emmeans_results, method = "pairwise")

emmeans_results <- emmeans(sla_model, ~ ot + pt)
pairs(emmeans_results, adjust = "tukey")
# see below - had to filter some pops out!


emmeans_results <- emmeans(ldmc_model2, ~ ot * pt)
pairs(emmeans_results, adjust = "tukey")
# see below - had to filter some pops out!


emmeans_results <- emmeans(mort_model, ~ ot * pt)
contrast(emmeans_results, method = "pairwise")

# reproductive
# reproductive models: flowerstructure_model3, flowerdays_model4
# seedmass_model1, seednum_model1, numflow_model
emmeans_results <- emmeans(flowerstructure_model3, ~ ot * pt | pop)
contrast(emmeans_results, method = "pairwise")

emmeans_results <- emmeans(flowerdays_model4, ~ ot * pt | pop)
contrast(emmeans_results, method = "pairwise")

emmeans_results <- emmeans(seedmass_model1, ~ ot * pt | pop)
contrast(emmeans_results, method = "pairwise")

emmeans_results <- emmeans(seednum_model1, ~ ot * pt | pop)
contrast(emmeans_results, method = "pairwise")

emmeans_results <- emmeans(numflow_model, ~ ot * pt)
contrast(emmeans_results, method = "pairwise")
custom_contrasts <- list(
  "CC_vs_DC" = c(1, -1, 0, 0),
  "DD_vs_CD" = c(0, 0, -1, 1))
contrast_results <- contrast(emmeans_results, custom_contrasts)
summary(contrast_results)


############# sla filtered #####
table(SLA_LDMC_FINAL$ot, SLA_LDMC_FINAL$pt, SLA_LDMC_FINAL$pop)

SLA_LDMC_FINAL_filtered <- SLA_LDMC_FINAL %>%
  filter(!pop %in% c(6, 2))  # Replace 6 and 2 with the actual identifiers for your populations

sla_model_filtered <- lmer(log(sla) ~ ot + pt + pop + ot:pt + pop:ot + pop:pt + pt:ot:pop + (1|tray), 
                           data = SLA_LDMC_FINAL_filtered, na.action = na.exclude)

emm_filtered <- emmeans(sla_model_filtered, ~ ot * pt)

pairs(emm_filtered, adjust = "tukey")

summary(sla_model_filtered)
sla_anova <- Anova(sla_model_filtered, type="III",  test.statistic= "F")
Anova(sla_model_filtered, type="III",  test.statistic= "F")
r.squaredGLMM(sla_model_filtered)




############# ldmc filtered #####
table(SLA_LDMC_FINAL$ot, SLA_LDMC_FINAL$pt, SLA_LDMC_FINAL$pop)

SLA_LDMC_FINAL_filtered <- SLA_LDMC_FINAL %>%
  filter(!pop %in% c(6, 2))  # Replace 6 and 2 with the actual identifiers for your populations

ldmc_model_filtered <- lmer(log(ldmc) ~ ot + pt + pop + ot:pt + pop:ot + pop:pt + pt:ot:pop + (1|tray), 
                           data = SLA_LDMC_FINAL_filtered, na.action = na.exclude)

emm_filtered <- emmeans(ldmc_model_filtered, ~ ot * pt)

pairs(emm_filtered, adjust = "tukey")

summary(ldmc_model_filtered)
sla_anova <- Anova(ldmc_model_filtered, type="III",  test.statistic= "F")
Anova(ldmc_model_filtered, type="III",  test.statistic= "F")
r.squaredGLMM(ldmc_model_filtered)


############# flowering rate filtered #####
table(flower_status$ot, flower_status$pt, flower_status$pop)

flower_status_filtered <- flower_status %>%
  filter(!pop %in% c(2))  # Replace 6 and 2 with the actual identifiers for your populations

numflow_model_filtered <- glm(status_flower ~ ot + pt + pop + ot:pt + pop:ot + pop:pt + pt:ot:pop,
                     data = flower_status_filtered, 
                     family = binomial(link = "logit"), 
                     na.action = na.exclude)

emm_filtered <- emmeans(numflow_model_filtered, ~ ot * pt)

pairs(emm_filtered, adjust = "tukey")
pairs(emm_filtered, adjust = "bonferroni")

summary(ldmc_model_filtered)
sla_anova <- Anova(ldmc_model_filtered, type="III",  test.statistic= "F")
Anova(ldmc_model_filtered, type="III",  test.statistic= "F")
r.squaredGLMM(ldmc_model_filtered)