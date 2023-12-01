# evaluate pH - recruitment  and temperarture-recruitment relationships

library(tidyverse)
library(mgcv)

library(rstan)
library(brms)
library(bayesplot)
source("./code/stan_utils.R")

theme_set(theme_bw())

cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# load data
rdat <- read.csv("./data/_23_0a_recruit_mfem_out.csv", row.names = 1)

ph_dat <- read.csv("./data/pH_annual_values.csv") %>%
  rename(year = Year,
         BB_pH = Bristol.Bay.mean) %>%
  select(year, BB_pH)

ph_dat <- ph_dat %>%
  mutate(BB_pH_lag4 = lag(BB_pH, 4),
         BB_pH_lag3 = lag(BB_pH, 3),
         BB_pH_lag2 = lag(BB_pH, 2),
         BB_pH_lag1 = lag(BB_pH, 1))


ph_dat$BB_ph_lag4_3 <- ph_dat$BB_ph_lag4_3_2 <- ph_dat$BB_ph_lag4_3_2_1 <- ph_dat$BB_ph_lag4_3_2_1_0 <- NA

for(i in 1:nrow(ph_dat)){
  
  ph_dat$BB_ph_lag4_3_2_1_0[i] <- mean(as.vector(c(ph_dat[i,2], ph_dat[i,3], ph_dat[i,4], ph_dat[i,5], ph_dat[i,6]))) 
  
  ph_dat$BB_ph_lag4_3_2_1[i] <- mean(as.vector(c(ph_dat[i,3], ph_dat[i,4], ph_dat[i,5], ph_dat[i,6])))
  
  ph_dat$BB_ph_lag4_3_2[i] <- mean(as.vector(c(ph_dat[i,3], ph_dat[i,4], ph_dat[i,5])))
  
  ph_dat$BB_ph_lag4_3[i] <- mean(as.vector(c(ph_dat[i,3], ph_dat[i,4])))
  
}

ph_dat <- ph_dat %>%
  select(1,3,10,9,8,7)

# load temperature DFA trend
temp_dat <- read.csv("./output/dfa_trend.csv") %>%
  rename(year = t, 
         temp_index = estimate) %>%
  select(year, temp_index)

temp_dat <- temp_dat %>%
  mutate(temp_index_lag4 = lag(temp_index, 4),
         temp_index_lag3 = lag(temp_index, 3),
         temp_index_lag2 = lag(temp_index, 2),
         temp_index_lag1 = lag(temp_index, 1))


temp_dat$temp_index_lag4_3 <- temp_dat$temp_index_lag4_3_2 <- temp_dat$temp_index_lag4_3_2_1 <- temp_dat$temp_index_lag4_3_2_1_0 <- NA

for(i in 1:nrow(temp_dat)){
  
  temp_dat$temp_index_lag4_3_2_1_0[i] <- mean(as.vector(c(temp_dat[i,2], temp_dat[i,3], temp_dat[i,4], temp_dat[i,5], temp_dat[i,6]))) 
  
  temp_dat$temp_index_lag4_3_2_1[i] <- mean(as.vector(c(temp_dat[i,3], temp_dat[i,4], temp_dat[i,5], temp_dat[i,6])))
  
  temp_dat$temp_index_lag4_3_2[i] <- mean(as.vector(c(temp_dat[i,3], temp_dat[i,4], temp_dat[i,5])))
  
  temp_dat$temp_index_lag4_3[i] <- mean(as.vector(c(temp_dat[i,3], temp_dat[i,4])))
  
}

temp_dat <- temp_dat %>%
  select(1,3,10,9,8,7)
# age distributions of recruits in assessment model
# these are the proportions assigned to
# 65-70mm, 70-75mm, 75-80mm, 80-85mm, 85-90mm, 90-95mm, 95-100mm

# Males proportions:
#   
#   0.26821432   0.32387433   0.24116113   0.11710765   0.03887032   0.00917743   0.00159482
# 
# Female proportions:
#   
#   0.27245822   0.38053177   0.24900488   0.08292921   0.01507593

# start with assumed age of 5 for model recruits

dat <- rdat %>%
  mutate(lag5_S = lag(mat_fem_GE90, 5),
         log_R_S = log(rec/lag5_S)) %>%
  select(year, lag5_S, log_R_S)

# add ph and temperature data
dat <- left_join(dat, ph_dat) %>%
  left_join(., temp_dat)

# remove 2022 for now (no ph estimates)
dat <- dat %>%
  filter(year <= 2021)

ggplot(dat, aes(lag5_S, log_R_S)) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4)) # overfit


ggplot(dat, aes(lag5_S, log_R_S)) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 3))

# evaluate S-R relationship before recent era of low R
mod <- gam(log_R_S ~ s(lag5_S, k = 3), data = dat[dat$year <= 2010,])
summary(mod)
plot(mod, resid = T, pch = 19)

# and for the full time series
mod <- gam(log_R_S ~ s(lag5_S, k = 3), data = dat)
summary(mod)
plot(mod, resid = T, pch = 19)

# no S-R relationship to account for

## fit brms models with ar and sigma ~ X-----------------------------------------------------------
priors <- c(set_prior("student_t(3, 0, 3)", class = "Intercept"),
            set_prior("student_t(3, 0, 3)", class = "b"),
            set_prior("student_t(3, 0, 3)", class = "sds"),
            set_prior("normal(0, 0.5)", class = "ar"))

## lag 4 only

# ph version
ph_form1 <- bf(log_R_S ~ 1 + s(BB_pH_lag4, k = 4) + ar(time = year, p = 1, cov = TRUE), sigma ~ BB_pH_lag4)

## fit ph
brms_ph_model1 <- brm(ph_form1,
                      data = dat,
                      prior = priors,
                      cores = 4, chains = 4, iter = 2000,
                      save_pars = save_pars(all = TRUE),
                      control = list(adapt_delta = 0.999, max_treedepth = 10))

saveRDS(brms_ph_model1, file = "./output/brms_ph_model1.rds")

brms_ph_model1 <- readRDS("./output/brms_ph_model1.rds")

check_hmc_diagnostics(brms_ph_model1$fit)

neff_lowest(brms_ph_model1$fit)

rhat_highest(brms_ph_model1$fit)

summary(brms_ph_model1)

bayes_R2(brms_ph_model1)

plot(conditional_smooths(brms_ph_model1), ask = FALSE)


## fit temp
temp_form1 <- bf(log_R_S ~ 1 + s(temp_index_lag4, k = 4) + ar(time = year, p = 1, cov = TRUE), sigma ~ temp_index_lag4)

brms_temp_model1 <- brm(temp_form1,
                      data = dat,
                      prior = priors,
                      cores = 4, chains = 4, iter = 2000,
                      save_pars = save_pars(all = TRUE),
                      control = list(adapt_delta = 0.999, max_treedepth = 10))

saveRDS(brms_temp_model1, file = "./output/brms_temp_model1.rds")

brms_temp_model1 <- readRDS("./output/brms_temp_model1.rds")

check_hmc_diagnostics(brms_temp_model1$fit)

neff_lowest(brms_temp_model1$fit)

rhat_highest(brms_temp_model1$fit)

summary(brms_temp_model1)

bayes_R2(brms_temp_model1)

plot(conditional_smooths(brms_temp_model1), ask = FALSE)

## lag 4_3 ##

# ph version
ph_form2 <- bf(log_R_S ~ 1 + s(BB_ph_lag4_3, k = 4) + ar(time = year, p = 1, cov = TRUE), sigma ~ BB_ph_lag4_3)

## fit ph
brms_ph_model2 <- brm(ph_form2,
                      data = dat,
                      prior = priors,
                      cores = 4, chains = 4, iter = 2000,
                      save_pars = save_pars(all = TRUE),
                      control = list(adapt_delta = 0.999, max_treedepth = 10))

saveRDS(brms_ph_model2, file = "./output/brms_ph_model2.rds")

brms_ph_model2 <- readRDS("./output/brms_ph_model2.rds")

check_hmc_diagnostics(brms_ph_model2$fit)

neff_lowest(brms_ph_model2$fit)

rhat_highest(brms_ph_model2$fit)

summary(brms_ph_model2)

bayes_R2(brms_ph_model2)

plot(conditional_smooths(brms_ph_model2), ask = FALSE)


## fit temp
temp_form2 <- bf(log_R_S ~ 1 + s(temp_index_lag4_3, k = 4) + ar(time = year, p = 1, cov = TRUE), sigma ~ temp_index_lag4_3)

brms_temp_model2 <- brm(temp_form2,
                        data = dat,
                        prior = priors,
                        cores = 4, chains = 4, iter = 2000,
                        save_pars = save_pars(all = TRUE),
                        control = list(adapt_delta = 0.999, max_treedepth = 10))

saveRDS(brms_temp_model2, file = "./output/brms_temp_model2.rds")

brms_temp_model2 <- readRDS("./output/brms_temp_model2.rds")

check_hmc_diagnostics(brms_temp_model2$fit)

neff_lowest(brms_temp_model2$fit)

rhat_highest(brms_temp_model2$fit)

summary(brms_temp_model2)

bayes_R2(brms_temp_model2)

plot(conditional_smooths(brms_temp_model2), ask = FALSE)

## lag 4_3_2
# ph version
ph_form3 <- bf(log_R_S ~ 1 + s(BB_ph_lag4_3_2, k = 4) + ar(time = year, p = 1, cov = TRUE), sigma ~ BB_ph_lag4_3_2)

## fit ph
brms_ph_model3 <- brm(ph_form3,
                      data = dat,
                      prior = priors,
                      cores = 4, chains = 4, iter = 2000,
                      save_pars = save_pars(all = TRUE),
                      control = list(adapt_delta = 0.999, max_treedepth = 10))

saveRDS(brms_ph_model3, file = "./output/brms_ph_model3.rds")

brms_ph_model3 <- readRDS("./output/brms_ph_model3.rds")

check_hmc_diagnostics(brms_ph_model3$fit)

neff_lowest(brms_ph_model3$fit)

rhat_highest(brms_ph_model3$fit)

summary(brms_ph_model3)

bayes_R2(brms_ph_model3)

plot(conditional_smooths(brms_ph_model3), ask = FALSE)


## fit temp
temp_form3 <- bf(log_R_S ~ 1 + s(temp_index_lag4_3_2, k = 4) + ar(time = year, p = 1, cov = TRUE), sigma ~ temp_index_lag4_3_2)

brms_temp_model3 <- brm(temp_form3,
                        data = dat,
                        prior = priors,
                        cores = 4, chains = 4, iter = 2000,
                        save_pars = save_pars(all = TRUE),
                        control = list(adapt_delta = 0.999, max_treedepth = 10))

saveRDS(brms_temp_model3, file = "./output/brms_temp_model3.rds")

brms_temp_model3 <- readRDS("./output/brms_temp_model3.rds")

check_hmc_diagnostics(brms_temp_model3$fit)

neff_lowest(brms_temp_model3$fit)

rhat_highest(brms_temp_model3$fit)

summary(brms_temp_model3)

bayes_R2(brms_temp_model3)

plot(conditional_smooths(brms_temp_model3), ask = FALSE)

## lag 4_3_2_1
# ph version
ph_form4 <- bf(log_R_S ~ 1 + s(BB_ph_lag4_3_2_1, k = 4) + ar(time = year, p = 1, cov = TRUE), sigma ~ BB_ph_lag4_3_2_1)

## fit ph
brms_ph_model4 <- brm(ph_form4,
                      data = dat,
                      prior = priors,
                      cores = 4, chains = 4, iter = 2000,
                      save_pars = save_pars(all = TRUE),
                      control = list(adapt_delta = 0.999, max_treedepth = 10))

saveRDS(brms_ph_model4, file = "./output/brms_ph_model4.rds")

brms_ph_model4 <- readRDS("./output/brms_ph_model4.rds")

check_hmc_diagnostics(brms_ph_model4$fit)

neff_lowest(brms_ph_model4$fit)

rhat_highest(brms_ph_model4$fit)

summary(brms_ph_model4)

bayes_R2(brms_ph_model4)

plot(conditional_smooths(brms_ph_model4), ask = FALSE)


## fit temp
temp_form4 <- bf(log_R_S ~ 1 + s(temp_index_lag4_3_2_1, k = 4) + ar(time = year, p = 1, cov = TRUE), sigma ~ temp_index_lag4_3_2_1)

brms_temp_model4 <- brm(temp_form4,
                        data = dat,
                        prior = priors,
                        cores = 4, chains = 4, iter = 2000,
                        save_pars = save_pars(all = TRUE),
                        control = list(adapt_delta = 0.999, max_treedepth = 10))

saveRDS(brms_temp_model4, file = "./output/brms_temp_model4.rds")

brms_temp_model4 <- readRDS("./output/brms_temp_model4.rds")

check_hmc_diagnostics(brms_temp_model4$fit)

neff_lowest(brms_temp_model4$fit)

rhat_highest(brms_temp_model4$fit)

summary(brms_temp_model4)

bayes_R2(brms_temp_model4)

plot(conditional_smooths(brms_temp_model4), ask = FALSE)

## lag 4_3_2_1_0
# ph version
ph_form5 <- bf(log_R_S ~ 1 + s(BB_ph_lag4_3_2_1_0, k = 4) + ar(time = year, p = 1, cov = TRUE), sigma ~ BB_ph_lag4_3_2_1_0)

## fit ph
brms_ph_model5 <- brm(ph_form5,
                      data = dat,
                      prior = priors,
                      cores = 4, chains = 4, iter = 2000,
                      save_pars = save_pars(all = TRUE),
                      control = list(adapt_delta = 0.999, max_treedepth = 10))

saveRDS(brms_ph_model5, file = "./output/brms_ph_model5.rds")

brms_ph_model5 <- readRDS("./output/brms_ph_model5.rds")

check_hmc_diagnostics(brms_ph_model5$fit)

neff_lowest(brms_ph_model5$fit)

rhat_highest(brms_ph_model5$fit)

summary(brms_ph_model5)

bayes_R2(brms_ph_model5)

plot(conditional_smooths(brms_ph_model5), ask = FALSE)


## fit temp
temp_form5 <- bf(log_R_S ~ 1 + s(temp_index_lag4_3_2_1_0, k = 4) + ar(time = year, p = 1, cov = TRUE), sigma ~ temp_index_lag4_3_2_1_0)

brms_temp_model5 <- brm(temp_form5,
                        data = dat,
                        prior = priors,
                        cores = 4, chains = 4, iter = 2000,
                        save_pars = save_pars(all = TRUE),
                        control = list(adapt_delta = 0.999, max_treedepth = 10))

saveRDS(brms_temp_model5, file = "./output/brms_temp_model5.rds")

brms_temp_model5 <- readRDS("./output/brms_temp_model5.rds")

check_hmc_diagnostics(brms_temp_model5$fit)

neff_lowest(brms_temp_model5$fit)

rhat_highest(brms_temp_model5$fit)

summary(brms_temp_model5)

bayes_R2(brms_temp_model5)

plot(conditional_smooths(brms_temp_model5), ask = FALSE)

#####################################
# load all the model objects

temp_brms1 <- readRDS("./output/brms_temp_model1.rds")
temp_brms2 <- readRDS("./output/brms_temp_model2.rds")
temp_brms3 <- readRDS("./output/brms_temp_model3.rds")
temp_brms4 <- readRDS("./output/brms_temp_model4.rds")
temp_brms5<- readRDS("./output/brms_temp_model5.rds")

ph_brms1 <- readRDS("./output/brms_ph_model1.rds")
ph_brms2 <- readRDS("./output/brms_ph_model2.rds")
ph_brms3 <- readRDS("./output/brms_ph_model3.rds")
ph_brms4 <- readRDS("./output/brms_ph_model4.rds")
ph_brms5<- readRDS("./output/brms_ph_model5.rds")

loo_compare <- loo(temp_brms1, temp_brms2, temp_brms3, temp_brms4, temp_brms5,
                   ph_brms1, ph_brms2, ph_brms3, ph_brms4, ph_brms5, moment_match = T)

# save 
saveRDS(loo_compare, "./output/ph_temp_brms_model_comparison.rds")

loo_compare <- readRDS("./output/ph_temp_brms_model_comparison.rds")

plot <- as.data.frame(loo_compare$diffs) %>%
  mutate(Covariate = rep(c("pH", "Warming index"), each = 5), 
           age_order = c(1,2,4,3,5,
                       5,4,1,2,3),
         Age = c("1-5", "1-4", "1-2", "1-3", "1",
                 "1", "1-2", "1-5", "1-4", "1-3")) 

# change se_diff = 0 to NA
change <- plot$se_diff == 0
plot$se_diff[change] <- NA

# order ages for plot
plot$Age <- reorder(plot$Age, plot$age_order)

# add difference from smallest looic value
plot <- plot %>%
  mutate(looic.diff = looic - min(looic))

compare_plot <- ggplot(plot, aes(Age, looic.diff, color = Covariate, fill = Covariate)) +
  geom_col(position = "dodge") +
  labs(x = "Age of modeled effect",
       y = "LOOIC difference") +
  scale_fill_manual(values = cb[c(2,6)]) +
  scale_color_manual(values = cb[c(2,6)])

compare_plot

# and temp v. pH compare for each set of lags
age_diff <- plot %>%
  select(Covariate, Age, looic) %>%
  pivot_wider(names_from = Covariate, 
              values_from = looic) %>%
  mutate(diff = `Warming index` - pH)

age_diff

# plot best models
bayes_R2(ph_brms5)

conditional_effects(ph_brms5)

summary(ph_brms5)

## plot pH effect

## 95% CI
ce1s_1 <- conditional_effects(ph_brms5, effect = "BB_ph_lag4_3_2_1_0", re_formula = NA,
                              prob = 0.95)
## 90% CI
ce1s_2 <- conditional_effects(ph_brms5, effect = "BB_ph_lag4_3_2_1_0", re_formula = NA,
                              prob = 0.9)
## 80% CI
ce1s_3 <- conditional_effects(ph_brms5, effect = "BB_ph_lag4_3_2_1_0", re_formula = NA,
                              prob = 0.8)
dat_ce <- ce1s_1$BB_ph_lag4_3_2_1_0

#########################

dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$BB_ph_lag4_3_2_1_0[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$BB_ph_lag4_3_2_1_0[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$BB_ph_lag4_3_2_1_0[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$BB_ph_lag4_3_2_1_0[["lower__"]]

plot_dat <- dat %>%
  mutate(year = str_sub(year, start = 3))

ph_model_plot <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1.5, color = "red3") +
  labs(x = "Mean pH, ages 1-5", y = "ln(R/S)") +
  theme_bw() + 
  geom_text(data = plot_dat, aes(x=BB_ph_lag4_3_2_1_0, y=log_R_S, label = year)) + 
  scale_x_reverse()

ph_model_plot

## plot temp effect
bayes_R2(temp_brms1)
## 95% CI
ce1s_1 <- conditional_effects(temp_brms1, effect = "temp_index_lag4", re_formula = NA,
                              prob = 0.95)
## 90% CI
ce1s_2 <- conditional_effects(temp_brms1, effect = "temp_index_lag4", re_formula = NA,
                              prob = 0.9)
## 80% CI
ce1s_3 <- conditional_effects(temp_brms1, effect = "temp_index_lag4", re_formula = NA,
                              prob = 0.8)
dat_ce <- ce1s_1$temp_index_lag4

#########################

dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$temp_index_lag4[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$temp_index_lag4[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$temp_index_lag4[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$temp_index_lag4[["lower__"]]


temp_model_plot <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1.5, color = "red3") +
  labs(x = "Temperature index, age 1", y = "ln(R/S)") +
  theme_bw() + 
  geom_text(data = plot_dat, aes(x=temp_index_lag4, y=log_R_S, label = year))

temp_model_plot


##
png("./figs/fig2.png", width = 5, height = 10, units = 'in', res = 300)

ggpubr::ggarrange(compare_plot, ph_model_plot, temp_model_plot, 
                  ncol = 1,
                  labels = "auto")

dev.off()



## plot time series ------------------

plot_rkc <- dat %>%
  select(year, mat_fem_GE90, rec) %>%
  filter(year >= 1975)

names(plot_rkc)[2:3] <- c("Mature females", "Recruits")

plot_rkc <- plot_rkc %>%
  pivot_longer(cols = -year)

crab_plot <- ggplot(plot_rkc, aes(year, value, color = name)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = cb[c(2,6)]) +
  theme(axis.title.x = element_blank(),
        legend.position = c(0.8, 0.8),
        legend.title = element_blank()) +
  ylab("Abundance (millions)")

crab_plot

plot_ph <- ggplot(filter(dat, year >= 1975), aes(year, BB_pH)) +
  geom_point() +
  geom_line() +
  theme(axis.title.x = element_blank()) +
  ylab("Mean pH")

plot_ph

plot_R.S <- ggplot(filter(dat, year >= 1975), aes(year, log_R_S)) +
  geom_point() +
  geom_line() +
  theme(axis.title.x = element_blank()) +
  ylab("ln(recruits/spawner)")

plot_R.S

# combine

png("./figs/time_series_plot.png", width = 4.5, height  = 8, units = 'in', res = 300)

ggpubr::ggarrange(plot_ph, crab_plot, plot_R.S, ncol = 1, labels = c("b", "c", "d"))

dev.off()
