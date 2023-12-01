# evaluate pH - recruitment  and temperarture-recruitment relationships

library(tidyverse)
library(mgcv)
library(nlme)
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



# start with assumed age of 5 for model recruits

dat <- rdat %>%
  mutate(lag5_S = lag(mat_fem_GE90, 5),
         log_R_S = log(rec/lag5_S)) %>%
  select(year, lag5_S, log_R_S)

# add ph data
dat <- left_join(dat, ph_dat) 

# remove 2022 for now (no ph estimates)
dat <- dat %>%
  filter(year <= 2021)

# and select only the columns we need, and remove years with NA
dat <- dat %>%
  select(year, log_R_S, BB_ph_lag4_3_2_1_0) %>% 
  filter(year >= 1980)

## loop through candidate thresholds--------------------

# minimum era length = 20% of years
0.2*nrow(dat) # 8 years is minimum window length

breaks <- 1987:2013
output <- data.frame()

# and loop through each
for(i in 1:length(breaks)){
  
  # i <- 1
  
  dat$era <- ifelse(dat$year <= breaks[i], "early", "late")
  
  mod <- gls(log_R_S ~ BB_ph_lag4_3_2_1_0*era, data = dat, correlation = corAR1())
    
  output <- rbind(output,
                  data.frame(breakpoint = breaks[i],
                             AICc = MuMIn::AICc(mod)))
  
}


ggplot(output, aes(breakpoint, AICc)) +
  geom_point() +
  geom_line()
