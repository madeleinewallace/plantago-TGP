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

# CLEANING -------

#seedf2_untidy <- read_excel("code-and-data/data/2023-greenhouse-data_f2maddiehannah-seed-data.xlsx")
#climate_untidy <- read_excel("code-and-data/data/2021-climate-info.xlsx", sheet = 2)

SEED_FINAL <- seedf2_untidy |> 
  mutate(tray = as.factor(tray)) |>
  mutate(id = as.factor(id)) |> 
  mutate(pop = as.factor(pop)) |> 
  mutate(TGP = as.factor(TGP)) |> 
  mutate(`2021_treat` = as.factor(`2021_treat`)) |> 
  mutate(`2023_treat` = as.factor(`2023_treat`)) |> 
  mutate(`2023_mass_total_g` = ifelse(`2023_mass_total_g` == 0, NA, `2023_mass_total_g`)) |> 
  mutate(`2023_num_total` = ifelse(`2023_num_total` == "#DIV/0!", NA, `2023_num_total`)) |> 
  mutate(`2023_num_total` = as.numeric(`2023_num_total`))|> 
  mutate(`2023_mass_total_g` = as.numeric(`2023_mass_total_g`))|> 
  mutate(ot = `2023_treat`) |> 
  mutate(pt = `2021_treat`) |> 
  mutate(num_total = `2023_num_total`) |> 
  mutate(mass_total = `2023_mass_total_g`) |> 
  mutate(across(c('num_total'), round, 0))

SEED_FINAL <- inner_join(SEED_FINAL, climate_tidy, by = c("pop" = "population")) 

SEED_FINAL <- SEED_FINAL |> 
  mutate(temp_cv = as.factor(temp_cv), 
         tray = as.factor(tray),
         id = as.factor(id), 
         cv_temp_mon = as.factor(cv_temp_mon), 
         cvlevel = as.factor(cvlevel), 
         temp = as.factor(temp),
         monsoon = as.factor(monsoon),
         MAP_level = as.factor(MAP_level),
         ot = as.factor(ot),
         pt = as.factor(pt),
         vpd = as.numeric(vpd),
         num_total = as.numeric(num_total),
         mass_total = as.numeric(mass_total))

filtered_seed_data <- SEED_FINAL[!SEED_FINAL$pop %in% c(3, 5, 11), ]
# pop 3, 5, and 11 only produced seed under DC condition, not any other the other 3 conditions
# Additionally, due to low flowering rates across populations in the greenhouse, 
# the final interaction term (OT:PT:pop) was removed from the model, as most factor combinations weren't present. 


# MODEL - NOT CENTERED
# ot and pt are FACTORS here, vpd is NUMERIC
seednum_model <- glmer(num_total ~ ot + pt + vpd + 
                         ot*vpd + pt*vpd + ot*pt + ot*pt*vpd + 
                         (1 | tray),
                       control = glmerControl(optimizer ="bobyqa", optCtrl=list(maxfun=2000000)),
                       family = poisson(),
                       data = filtered_seed_data,
                       na.action = na.exclude)
# wouldn't converge - "Model is nearly unidentifiable: very large eigenvalue Rescale variables?"


# here, ot and pt are CHARACTERS, vpd is NUMERIC
filtered_seed_data <- filtered_seed_data |> 
  mutate(ot = as.character(ot), 
         pt = as.character(pt))

seednum_model <- glmer(num_total ~ ot + pt + vpd + 
                         ot*vpd + pt*vpd + ot*pt + ot*pt*vpd + 
                         (1 | tray),
                       control = glmerControl(optimizer ="bobyqa", optCtrl=list(maxfun=2000000)),
                       family = poisson(),
                       data = filtered_seed_data,
                       na.action = na.exclude)
# wouldn't converge - "Model is nearly unidentifiable: very large eigenvalue Rescale variables?"

# here, i convert vpd (hPa) into (kPa) and then use the kPa values in the model
filtered_seed_data <- filtered_seed_data %>%
  mutate(vpd.k = vpd / 10)

seednum_model <- glmer(num_total ~ ot + pt + vpd.k + 
                         ot*vpd.k + pt*vpd.k + ot*pt + ot*pt*vpd.k + 
                         (1 | tray),
                       control = glmerControl(optimizer ="bobyqa", optCtrl=list(maxfun=2000000)),
                       family = poisson(),
                       data = filtered_seed_data,
                       na.action = na.exclude)

plot(seednum_model, add.smooth = FALSE, which = 1)
E <- resid(seednum_model)
hist(E, xlab = "residuals", main = " ")
qqnorm(resid(seednum_model))
qqline(resid(seednum_model))              

summary(seednum_model)
r.squaredGLMM(seednum_model)
Anova(seednum_model, type = 3)

# that did work without centering it, but results are nearly the same, which is good!!!!