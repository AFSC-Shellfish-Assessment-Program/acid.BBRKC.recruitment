# plot Fig1

library(tidyverse)
library(ggpubr)

theme_set(theme_bw())

cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

rdat <- read.csv("./data/_23_0a_recruit_mfem_out.csv", row.names = 1)

plot_rkc <- rdat %>%
  select(year, mat_fem_GE90, rec) %>%
  filter(year < 2023)

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

ph_dat <- read.csv("./data/pH_annual_values.csv") %>%
  rename(year = Year,
         BB_pH = Bristol.Bay.mean) %>%
  select(year, BB_pH)

ph_plot <- ggplot(filter(ph_dat, year >= 1975), aes(year, BB_pH)) +
  geom_point() +
  geom_line() +
  theme(axis.title.x = element_blank()) +
  ylab("Mean pH")

ph_plot

plot_recr <- plot_rkc %>%
  pivot_wider(names_from = name, values_from = value) %>%
  mutate(`Mature females` = lag(`Mature females`, 5),
         log_R_S = log(Recruits/`Mature females`))
  

R.S_plot <- ggplot(filter(plot_recr, year >= 1980), aes(year, log_R_S)) +
  geom_point() +
  geom_line() +
  theme(axis.title.x = element_blank()) +
  ylab("ln(recruits/spawner)")

R.S_plot


temp_ice <- read.csv("./data/sst_temp_original_units.csv")

ice <- temp_ice %>%
  select(1:3) %>%
  rename(`January-February` = Jan_Feb_plot_ice,
         `March-April` = Mar_Apr_plot_ice) %>%
  pivot_longer(cols = -year,
               values_to = "Ice concentration")

ice_plot <- ggplot(ice, aes(year, `Ice concentration`, color = name)) +
  geom_point() + 
  geom_line() +
  scale_color_manual(values = cb[c(2,4)]) +
  theme(legend.title = element_blank(),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.margin = margin(t = -1, r = 0, b = -1, l = 0, unit = "pt"),
        axis.title.x = element_blank())

ice_plot

temp <- temp_ice %>%
  select(c(1,4,5)) %>%
  rename(`January-June surface` = Jan_Jun_SST,
         `June bottom` = June_bottom_temp) %>%
  pivot_longer(cols = -year,
               values_to = "Temperature (°C)")

temp_plot <- ggplot(temp, aes(year, `Temperature (°C)`, color = name)) +
  geom_point() + 
  geom_line() +
  scale_color_manual(values = cb[c(6,7)]) +
  theme(legend.title = element_blank(),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.margin = margin(t = -1, r = 0, b = -1, l = 0, unit = "pt"),
        axis.title.x = element_blank())


temp_plot

plot.CI <- read.csv("./Output/dfa_loadings.csv")
trend <- read.csv("./Output/dfa_trend.csv")

# plot loadings and trend

dodge <- position_dodge(width=0.9)

# clean up names
plot.CI$names <- c("June bot. temp.", 
                   "Jan-Feb ice",
                   "Mar-Apr ice",
                   "Jan-Jun SST")

plot.CI$names <- reorder(plot.CI$names, plot.CI$mean)



loadings.plot <- ggplot(plot.CI, aes(x=names, y=mean)) +
  geom_bar(position=dodge, stat="identity", fill=cb[6]) +
  geom_errorbar(aes(ymax=upCI, ymin=lowCI), position=dodge, width=0.5) +
  ylab("Loading") +
  xlab("") +
  theme(axis.text.x  = element_text(angle=60, hjust=1,  size=9), legend.title = element_blank(), legend.position = 'top') +
  geom_hline(yintercept = 0)

trend.plot <- ggplot(trend, aes(t, estimate)) +
  theme_bw() +
  geom_line(color=cb[6]) +
  geom_hline(yintercept = 0) +
  geom_point(color=cb[6]) +
  geom_ribbon(aes(x=t, ymin=conf.low, ymax=conf.high), linetype=2, alpha=0.1, fill=cb[6]) + xlab("") + ylab("Temperature index")

# and combine
void.plot <- ggplot() + theme_void()

png("./figs/fig1.png", width = 8, height = 10, units = 'in', res = 300)

ggarrange(void.plot, crab_plot, R.S_plot, ph_plot, ice_plot, temp_plot, loadings.plot, trend.plot,
          ncol = 2, nrow = 4,
          labels = "auto")

dev.off()

