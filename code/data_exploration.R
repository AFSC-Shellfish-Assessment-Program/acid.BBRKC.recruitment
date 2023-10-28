# evaluate pH - recruitment relationship

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
phdat <- read.csv("./data/pH_annual_values.csv") %>%
  rename(year = Year,
         BB_pH = Bristol.Bay.mean,
         BB_pH50m = Bristol.Bay...50m.mean) %>%
  select(year, BB_pH, BB_pH50m)

dat <- left_join(phdat, rdat) %>%
  mutate(log_R = log(rec),
         log_S = log(mat_fem_GE90),
         lag_BB_ph = lag(BB_pH),
         lag_BB_ph50m = lag(BB_pH50m))

ggplot(dat, aes(lag_BB_ph50m, log_R)) +
  geom_point()

ggplot(dat, aes(lag_BB_ph, log_R)) +
  geom_point()

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

dat <- dat %>%
  mutate(lag5_S = lag(mat_fem_GE90, 5),
         log_R_S = log(rec/lag5_S))


ggplot(dat, aes(lag5_S, rec)) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4))

ggplot(dat, aes(lag5_S, log_R_S)) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4)) # overfit


ggplot(dat, aes(lag5_S, log_R_S)) +
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 3))


mod1_truncated <- gam(log_R_S ~ s(lag5_S, k = 3), data = dat[dat$year %in% 1977:2010,])
summary(mod1_truncated)
plot(mod1_truncated, resid = T, pch = 19)

mod1 <- gam(log_R_S ~ s(lag5_S, k = 3), data = dat)
summary(mod1)
plot(mod1, resid = T, pch = 19)

# no S-R relationship to account for

# add in pH data at appropriate lags
temp <- dat %>%
  filter(year >= 1980)

ccf(temp$BB_pH, temp$log_R_S)


dat <- dat %>%
  mutate(BB_pH_lag4 = lag(BB_pH, 4),
         BB_pH_lag3 = lag(BB_pH, 3),
         BB_pH_lag2 = lag(BB_pH, 2),
         BB_pH_lag1 = lag(BB_pH, 1))


dat$BB_ph_lag4_3 <- dat$BB_ph_lag4_3_2 <- dat$BB_ph_lag4_3_2_1 <- dat$BB_ph_lag4_3_2_1_0 <- NA

for(i in 5:nrow(dat)){
  
  dat$BB_ph_lag4_3_2_1_0[i] <- mean(as.vector(c(dat[i,2], dat[i,12], dat[i,13], dat[i,14], dat[i,15]))) 

  dat$BB_ph_lag4_3_2_1[i] <- mean(as.vector(c(dat[i,12], dat[i,13], dat[i,14], dat[i,15])))

  dat$BB_ph_lag4_3_2[i] <- mean(as.vector(c(dat[i,12], dat[i,13], dat[i,14])))

  dat$BB_ph_lag4_3[i] <- mean(as.vector(c(dat[i,12], dat[i,13])))
  
}

# start with the expectation that we observe a pH effect at lag 4 (1-yr-olds) and build out with 
# more lags to account for possible effects at age 2, 3, 4, 5

mod2 <- gam(log_R_S ~ s(BB_pH_lag4, k = 4), data = dat)
summary(mod2)

mod3 <- gam(log_R_S ~ s(BB_ph_lag4_3, k = 4), data = dat)
summary(mod3)

mod4 <- gam(log_R_S ~ s(BB_ph_lag4_3_2, k = 4), data = dat)
summary(mod4)

mod5 <- gam(log_R_S ~ s(BB_ph_lag4_3_2_1, k = 4), data = dat)
summary(mod5)

mod6 <- gam(log_R_S ~ s(BB_ph_lag4_3_2_1_0, k = 4), data = dat)
summary(mod6)


MuMIn::AICc(mod2, mod3, mod4, mod5, mod6) # very similar, mod6 is *slightly* better than mod2

gam.check(mod2)

plot(mod2, resid = T, se = F, pch = 19)

gam.check(mod6)

plot(mod6, resid = T, se = F, pch = 19)



## fit brms version of model 2
form <- bf(log_R_S ~ 1 + s(BB_pH_lag4, k = 4) + ar(time = year, p = 1, cov = TRUE))

priors <- c(set_prior("student_t(3, 0, 3)", class = "Intercept"),
            set_prior("student_t(3, 0, 3)", class = "b"),
            set_prior("student_t(3, 0, 3)", class = "sds"),
            set_prior("student_t(3, 0, 3)", class = "sigma"),
            set_prior("normal(0, 0.5)", class = "ar"))



## fit model with covariance in ar term --------------------------------------
brms_model2 <- brm(form,
                      data = dat,
                      prior = priors,
                      cores = 4, chains = 4, iter = 3000,
                      save_pars = save_pars(all = TRUE),
                      control = list(adapt_delta = 0.999, max_treedepth = 10))

saveRDS(brms_model2, file = "./output/brms_model2.rds")

brms_model2 <- readRDS("./output/brms_model2.rds")

check_hmc_diagnostics(brms_model2$fit)

neff_lowest(brms_model2$fit)

rhat_highest(brms_model2$fit)

summary(brms_model2)

bayes_R2(brms_model2)

plot(conditional_smooths(brms_model2), ask = FALSE)

######
## fit brms version of model 3
form <- bf(log_R_S ~ 1 + s(BB_ph_lag4_3, k = 4) + ar(time = year, p = 1, cov = TRUE))

priors <- c(set_prior("student_t(3, 0, 3)", class = "Intercept"),
            set_prior("student_t(3, 0, 3)", class = "b"),
            set_prior("student_t(3, 0, 3)", class = "sds"),
            set_prior("student_t(3, 0, 3)", class = "sigma"),
            set_prior("normal(0, 0.5)", class = "ar"))



## fit model with covariance in ar term --------------------------------------
brms_model3 <- brm(form,
                   data = dat,
                   prior = priors,
                   cores = 4, chains = 4, iter = 3000,
                   save_pars = save_pars(all = TRUE),
                   control = list(adapt_delta = 0.999, max_treedepth = 10))

saveRDS(brms_model3, file = "./output/brms_model3.rds")

brms_model3 <- readRDS("./output/brms_model3.rds")

check_hmc_diagnostics(brms_model3$fit)

neff_lowest(brms_model3$fit)

rhat_highest(brms_model3$fit)

summary(brms_model3)

bayes_R2(brms_model3)

plot(conditional_smooths(brms_model3), ask = FALSE)


######

######
## fit brms version of model 4
form <- bf(log_R_S ~ 1 + s(BB_ph_lag4_3_2, k = 4) + ar(time = year, p = 1, cov = TRUE))

priors <- c(set_prior("student_t(3, 0, 3)", class = "Intercept"),
            set_prior("student_t(3, 0, 3)", class = "b"),
            set_prior("student_t(3, 0, 3)", class = "sds"),
            set_prior("student_t(3, 0, 3)", class = "sigma"),
            set_prior("normal(0, 0.5)", class = "ar"))



## fit model with covariance in ar term --------------------------------------
brms_model4 <- brm(form,
                   data = dat,
                   prior = priors,
                   cores = 4, chains = 4, iter = 3000,
                   save_pars = save_pars(all = TRUE),
                   control = list(adapt_delta = 0.999, max_treedepth = 10))

saveRDS(brms_model4, file = "./output/brms_model4.rds")

brms_model4 <- readRDS("./output/brms_model4.rds")

check_hmc_diagnostics(brms_model4$fit)

neff_lowest(brms_model4$fit)

rhat_highest(brms_model4$fit)

summary(brms_model4)

bayes_R2(brms_model4)

plot(conditional_smooths(brms_model4), ask = FALSE)


######
######
## fit brms version of model 5
form <- bf(log_R_S ~ 1 + s(BB_ph_lag4_3_2_1, k = 4) + ar(time = year, p = 1, cov = TRUE))

priors <- c(set_prior("student_t(3, 0, 3)", class = "Intercept"),
            set_prior("student_t(3, 0, 3)", class = "b"),
            set_prior("student_t(3, 0, 3)", class = "sds"),
            set_prior("student_t(3, 0, 3)", class = "sigma"),
            set_prior("normal(0, 0.5)", class = "ar"))



## fit model with covariance in ar term --------------------------------------
brms_model5 <- brm(form,
                   data = dat,
                   prior = priors,
                   cores = 4, chains = 4, iter = 3000,
                   save_pars = save_pars(all = TRUE),
                   control = list(adapt_delta = 0.999, max_treedepth = 10))

saveRDS(brms_model5, file = "./output/brms_model5.rds")

brms_model5 <- readRDS("./output/brms_model5.rds")

check_hmc_diagnostics(brms_model5$fit)

neff_lowest(brms_model5$fit)

rhat_highest(brms_model5$fit)

summary(brms_model5)

bayes_R2(brms_model5)

plot(conditional_smooths(brms_model5), ask = FALSE)


######



## fit brms version of model 6
form <- bf(log_R_S ~ 1 + s(BB_ph_lag4_3_2_1_0, k = 4) + ar(time = year, p = 1, cov = TRUE))

priors <- c(set_prior("student_t(3, 0, 3)", class = "Intercept"),
            set_prior("student_t(3, 0, 3)", class = "b"),
            set_prior("student_t(3, 0, 3)", class = "sds"),
            set_prior("student_t(3, 0, 3)", class = "sigma"),
            set_prior("normal(0, 0.5)", class = "ar"))



## fit model with covariance in ar term --------------------------------------
brms_model6 <- brm(form,
                   data = dat,
                   seed = 99,
                   prior = priors,
                   cores = 4, chains = 4, iter = 3000,
                   save_pars = save_pars(all = TRUE),
                   control = list(adapt_delta = 0.999, max_treedepth = 10))

saveRDS(brms_model6, file = "./output/brms_model6.rds")

brms_model6 <- readRDS("./output/brms_model6.rds")

check_hmc_diagnostics(brms_model6$fit)

neff_lowest(brms_model6$fit)

rhat_highest(brms_model6$fit)

summary(brms_model6)

bayes_R2(brms_model6)

plot(conditional_smooths(brms_model6), ask = FALSE)

###############
## fit version with sigma modeled on ph
form <- bf(log_R_S ~ 1 + s(BB_ph_lag4_3_2_1_0, k = 4) + ar(time = year, p = 1, cov = TRUE), sigma ~ BB_ph_lag4_3_2_1_0)

priors <- c(set_prior("student_t(3, 0, 3)", class = "Intercept"),
            set_prior("student_t(3, 0, 3)", class = "b"),
            set_prior("student_t(3, 0, 3)", class = "sds"),
            set_prior("student_t(3, 0, 3)", class = "sigma"),
            set_prior("normal(0, 0.5)", class = "ar"))



## fit model with covariance in ar term --------------------------------------
brms_model6b <- brm(form,
                   data = dat,
                   seed = 999,
                   # prior = priors,
                   cores = 4, chains = 4, iter = 3000,
                   save_pars = save_pars(all = TRUE),
                   control = list(adapt_delta = 0.999, max_treedepth = 10))

saveRDS(brms_model6b, file = "./output/brms_model6b.rds")

brms_model6b <- readRDS("./output/brms_model6b.rds")

check_hmc_diagnostics(brms_model6b$fit)

neff_lowest(brms_model6b$fit)

rhat_highest(brms_model6b$fit)

summary(brms_model6b)

bayes_R2(brms_model6b)

plot(conditional_smooths(brms_model6b), ask = FALSE)


###############

####
# load all the model objects

brms2 <- readRDS("./output/brms_model2.rds")
brms3 <- readRDS("./output/brms_model3.rds")
brms4 <- readRDS("./output/brms_model4.rds")
brms5 <- readRDS("./output/brms_model5.rds")
brms6 <- readRDS("./output/brms_model6.rds")
brms6b <- readRDS("./output/brms_model6b.rds")

loo_compare <- loo(brms2, brms3, brms4, brms5, brms6)

str(loo_compare)

# save 
saveRDS(loo_compare, "./output/brms_model_comparison.rds")

loo_compare <- readRDS("./output/brms_model_comparison.rds")

plot <- as.data.frame(loo_compare$diffs) %>%
  mutate(model_number = 1:5,
         Age = c("1-5", "1", "1-4", "1-2", "1-3")) %>%
  arrange(model_number)

# change se_diff = 0 to NA
change <- plot$se_diff == 0
plot$se_diff[change] <- NA


compare_plot <- ggplot(plot, aes(model_number, elpd_diff)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = elpd_diff - 1.96*se_diff,
                    ymax = elpd_diff + 1.96*se_diff)) +
  scale_x_continuous(breaks = 1:5, labels = plot$Age) +
  labs(x = "Age of modeled effect",
       y = "ELPD difference")

# really no difference, but model 6 nominally better



loo(brms6, brms6b)

## plot pH effect

## 95% CI
ce1s_1 <- conditional_effects(brms6b, effect = "BB_ph_lag4_3_2_1_0", re_formula = NA,
                              prob = 0.95)
## 90% CI
ce1s_2 <- conditional_effects(brms6b, effect = "BB_ph_lag4_3_2_1_0", re_formula = NA,
                              prob = 0.9)
## 80% CI
ce1s_3 <- conditional_effects(brms6b, effect = "BB_ph_lag4_3_2_1_0", re_formula = NA,
                              prob = 0.8)
dat_ce <- ce1s_1$BB_ph_lag4_3_2_1_0

#########################

dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$BB_ph_lag4_3_2_1_0[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$BB_ph_lag4_3_2_1_0[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$BB_ph_lag4_3_2_1_0[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$BB_ph_lag4_3_2_1_0[["lower__"]]



g2 <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1.5, color = "red3") +
  labs(x = "Bristol Bay mean pH, ages 1-5", y = "ln(R/S)") +
  theme_bw() + 
  geom_text(data = dat, aes(x=BB_ph_lag4_3_2_1_0, y=log_R_S, label = year)) + 
  scale_x_reverse()

print(g2)

ggsave("./figs/brms_model6b_pH_effect.png", width = 4.5, height = 3, units = 'in')

## model selection in brms------------------

# using the same formulation as for mgcv model selection above
mod3 <- gam(log_R_S ~ s(BB_ph_lag4_3, k = 4), data = dat)

mod4 <- gam(log_R_S ~ s(BB_ph_lag4_3_2, k = 4), data = dat)

mod5 <- gam(log_R_S ~ s(BB_ph_lag4_3_2_1, k = 4), data = dat)

mod6 <- gam(log_R_S ~ s(BB_ph_lag4_3_2_1_0, k = 4), data = dat)


## fit brms version of model 3
form <- bf(log_R_S ~ 1 + s(BB_ph_lag4_3, k = 4) + ar(time = year, p = 1, cov = TRUE))

priors <- c(set_prior("student_t(3, 0, 3)", class = "Intercept"),
            set_prior("student_t(3, 0, 3)", class = "b"),
            set_prior("student_t(3, 0, 3)", class = "sds"),
            set_prior("student_t(3, 0, 3)", class = "sigma"),
            set_prior("normal(0, 0.5)", class = "ar"))



## fit model with covariance in ar term 
brms_model3 <- brm(form,
                   data = dat,
                   prior = priors,
                   cores = 4, chains = 4, iter = 3000,
                   save_pars = save_pars(all = TRUE),
                   control = list(adapt_delta = 0.999, max_treedepth = 10))

saveRDS(brms_model3, file = "./output/brms_model3.rds")

## fit brms version of model 4
form <- bf(log_R_S ~ 1 + s(BB_ph_lag4_3_2, k = 4) + ar(time = year, p = 1, cov = TRUE))

priors <- c(set_prior("student_t(3, 0, 3)", class = "Intercept"),
            set_prior("student_t(3, 0, 3)", class = "b"),
            set_prior("student_t(3, 0, 3)", class = "sds"),
            set_prior("student_t(3, 0, 3)", class = "sigma"),
            set_prior("normal(0, 0.5)", class = "ar"))



## fit model with covariance in ar term 
brms_model4 <- brm(form,
                   data = dat,
                   prior = priors,
                   cores = 4, chains = 4, iter = 3000,
                   save_pars = save_pars(all = TRUE),
                   control = list(adapt_delta = 0.999, max_treedepth = 10))

saveRDS(brms_model4, file = "./output/brms_model4.rds")

## fit brms version of model 5
form <- bf(log_R_S ~ 1 + s(BB_ph_lag4_3_2_1, k = 4) + ar(time = year, p = 1, cov = TRUE))

priors <- c(set_prior("student_t(3, 0, 3)", class = "Intercept"),
            set_prior("student_t(3, 0, 3)", class = "b"),
            set_prior("student_t(3, 0, 3)", class = "sds"),
            set_prior("student_t(3, 0, 3)", class = "sigma"),
            set_prior("normal(0, 0.5)", class = "ar"))



## fit model with covariance in ar term 
brms_model5 <- brm(form,
                   data = dat,
                   prior = priors,
                   cores = 4, chains = 4, iter = 3000,
                   save_pars = save_pars(all = TRUE),
                   control = list(adapt_delta = 0.999, max_treedepth = 10))

saveRDS(brms_model5, file = "./output/brms_model5.rds")

## brms model selection ---------------------------------

# load all the model objects

brms2 <- readRDS("./output/brms_model2.rds")
brms3 <- readRDS("./output/brms_model3.rds")
brms4 <- readRDS("./output/brms_model4.rds")
brms5 <- readRDS("./output/brms_model5.rds")
brms6 <- readRDS("./output/brms_model6.rds")


loo_compare <- loo(brms2, brms3, brms4, brms5, brms6)

str(loo_compare)

# save 
saveRDS(loo_compare, "./output/brms_model_comparison.rds")

loo_compare <- readRDS("./output/brms_model_comparison.rds")

plot <- as.data.frame(loo_compare$diffs) %>%
  mutate(model_number = 1:5,
         Age = c("1-5", "1", "1-4", "1-2", "1-3")) %>%
  arrange(model_number)

# change se_diff = 0 to NA
change <- plot$se_diff == 0
plot$se_diff[change] <- NA


compare_plot <- ggplot(plot, aes(model_number, elpd_diff)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = elpd_diff - 1.96*se_diff,
                    ymax = elpd_diff + 1.96*se_diff)) +
  scale_x_continuous(breaks = 1:5, labels = plot$Age) +
  labs(x = "Age of modeled effect",
       y = "ELPD difference")

# combine with g2 (best model)

png("./figs/model_selection_best_model_plot.png", width = 9.5, height = 4, units = 'in', res = 300)

ggpubr::ggarrange(compare_plot, g2, ncol = 2, widths = c(0.45, 0.55), labels = "auto")

dev.off()

## plot least-supported model -------------------------

brms_model4 <- readRDS("./output/brms_model4.rds")

check_hmc_diagnostics(brms_model4$fit)

neff_lowest(brms_model4$fit)

rhat_highest(brms_model4$fit)

summary(brms_model4)

bayes_R2(brms_model6)

## plot 

# pH effect

## 95% CI
ce1s_1 <- conditional_effects(brms_model4, effect = "BB_ph_lag4_3_2", re_formula = NA,
                              prob = 0.95)
## 90% CI
ce1s_2 <- conditional_effects(brms_model4, effect = "BB_ph_lag4_3_2", re_formula = NA,
                              prob = 0.9)
## 80% CI
ce1s_3 <- conditional_effects(brms_model4, effect = "BB_ph_lag4_3_2", re_formula = NA,
                              prob = 0.8)
dat_ce <- ce1s_1$BB_ph_lag4_3_2

#########################

dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$BB_ph_lag4_3_2[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$BB_ph_lag4_3_2[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$BB_ph_lag4_3_2[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$BB_ph_lag4_3_2[["lower__"]]



g_mod4 <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1.5, color = "red3") +
  labs(x = "Bristol Bay mean pH, ages 1-3", y = "ln(R/S)") +
  theme_bw() + 
  geom_text(data = dat, aes(x=BB_ph_lag4_3_2, y=log_R_S, label = year)) + 
  scale_x_reverse()

print(g_mod4)

ggsave("./figs/brms_model6_pH_effect.png", width = 4.5, height = 3, units = 'in')

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
