# evaluate pH - recruitment relationship

library(tidyverse)
library(mgcv)

library(rstan)
library(brms)
library(bayesplot)
source("./code/stan_utils.R")

theme_set(theme_bw())

# load data
rdat <- read.csv("./data/recruit_mfem_out.csv", row.names = 1)
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

# cross-correlations for S and R
check <- dat %>%
  filter(year %in% 1976:2010)

png("./figs/R_S_ccf.png", width = 5, height = 5, units = 'in', res = 300)
ccf(check$log_S, check$log_R, main = "log(S), log(R) cross-correlation", ) # looks like lag 2 is strongest
dev.off()


check.ccf <- ccf(check$log_S, check$log_R)

check.ccf

dat <- dat %>%
  mutate(lag2_S = lag(mat_fem_GE90, 2),
         log_R_S = log(rec/lag2_S),
         BB_pH_lag2 = lag(BB_pH, 2),
         BB_pH50m_lag2 = lag(BB_pH50m, 2))

dat$BB_ph_lag0_1 <- dat$BB_ph_lag1_2 <- dat$BB_ph_50m_lag1_2 <- NA

for(i in 3:nrow(dat)){
 
dat$BB_ph_lag0_1[i] <- mean(dat[i,2], dat[i,8]) 
dat$BB_ph_lag1_2[i] <- mean(dat[i,8], dat[i,12])
dat$BB_ph_50m_lag1_2[i] <- mean(dat[i,9], dat[i,13])
}

mod1_truncated <- gam(log_R_S ~ s(lag2_S, k = 4), data = dat[dat$year %in% 1977:2010,])
summary(mod1_truncated)
plot(mod1_truncated)

mod1 <- gam(log_R_S ~ s(lag2_S, k = 4), data = dat)
summary(mod1)

mod2 <- gam(log_R_S ~ s(lag2_S, k = 4) + s(BB_ph_lag1_2, k = 4), data = dat)
summary(mod2)

mod3 <- gam(log_R_S ~ s(lag2_S, k = 4) + s(BB_ph_50m_lag1_2, k = 4), data = dat)
summary(mod3)

MuMIn::AICc(mod1, mod2, mod3) # mod 2 and 3 are nearly identical

gam.check(mod2)

plot(mod2, resid = T, se = F, pch = 19)

mod4 <- gam(log_R_S ~ s(lag2_S, k = 4) + s(lag_BB_ph, k = 4), data = dat)
summary(mod4)

mod5 <- gam(log_R_S ~ s(lag2_S, k = 4) + s(BB_pH_lag2, k = 4), data = dat)
summary(mod5)

mod6 <- gam(log_R_S ~ s(lag2_S, k = 4) + s(BB_pH, k = 4), data = dat)
summary(mod6)

MuMIn::AICc(mod1, mod2, mod3, mod4, mod5, mod6)

mod7 <- gam(log_R_S ~ s(lag2_S, k = 4) + s(BB_ph_lag0_1, k = 4), data = dat)
summary(mod7)

MuMIn::AICc(mod1, mod2, mod3, mod4, mod5, mod6, mod7)

## fit brms version of model 2
form <- bf(log_R_S ~ 1 + s(lag2_S) + s(BB_ph_lag1_2) + ar(time = year, p = 1, cov = TRUE))

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

## plot 


# first, spawner effect

## 95% CI
ce1s_1 <- conditional_effects(brms_model2, effect = "lag2_S", re_formula = NA,
                              prob = 0.95)
## 90% CI
ce1s_2 <- conditional_effects(brms_model2, effect = "lag2_S", re_formula = NA,
                              prob = 0.9)
## 80% CI
ce1s_3 <- conditional_effects(brms_model2, effect = "lag2_S", re_formula = NA,
                              prob = 0.8)
dat_ce <- ce1s_1$lag2_S

#########################

dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$lag2_S[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$lag2_S[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$lag2_S[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$lag2_S[["lower__"]]
dat_ce[["rug.anom"]] <- c(jitter(unique(dat$lag2_S), amount = 0.0051),
                          rep(NA, 100-length(unique(dat$lag2_S))))


g1 <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1.5, color = "red3") +
  labs(x = "Mature female abundance, lag 2", y = "log(R/S)") +
  theme_bw() +
  geom_hline(yintercept = 0, lty = 2) +
  geom_rug(aes(x=rug.anom, y=NULL)) 

print(g1)

# now, pH effect

## 95% CI
ce1s_1 <- conditional_effects(brms_model2, effect = "BB_ph_lag1_2", re_formula = NA,
                              prob = 0.95)
## 90% CI
ce1s_2 <- conditional_effects(brms_model2, effect = "BB_ph_lag1_2", re_formula = NA,
                              prob = 0.9)
## 80% CI
ce1s_3 <- conditional_effects(brms_model2, effect = "BB_ph_lag1_2", re_formula = NA,
                              prob = 0.8)
dat_ce <- ce1s_1$BB_ph_lag1_2

#########################

dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$BB_ph_lag1_2[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$BB_ph_lag1_2[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$BB_ph_lag1_2[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$BB_ph_lag1_2[["lower__"]]
dat_ce[["rug.anom"]] <- c(jitter(unique(dat$BB_ph_lag1_2), amount = 0.0051),
                          rep(NA, 100-length(unique(dat$BB_ph_lag1_2))))


g2 <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1.5, color = "red3") +
  labs(x = "Bristol Bay mean pH, lag 1-2", y = "log(R/S)") +
  theme_bw() + 
  geom_rug(aes(x=rug.anom, y=NULL)) + 
  scale_x_reverse()

print(g2)

ggsave("./figs/brms_model2_pH_effect.png", width = 4.5, height = 3, units = 'in')
