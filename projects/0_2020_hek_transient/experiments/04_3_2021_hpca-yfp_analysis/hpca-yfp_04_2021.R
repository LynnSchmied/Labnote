# HPCA-EYFP data analysis, uncaging with different stimulation power
# Copyright © 2021 Borys Olifirov

require(dplyr)
require(tidyr)
require(purrr)

require(TTR)

require(ggplot2)
require(ggpubr)
require(ggsci)

setwd('/home/astria/bio/note/projects/0_2020_hek_transient/experiments/04_3_2021_hpca-yfp_analysis')

df.fluo <- read.csv('results_fluo.csv') %>%
           subset(select = c(power, time, int)) %>%
           group_by(time, power) %>%
           summarise_all(list(mean, sd), na.rm = TRUE) %>%
           rename(mean = fn1, sd = fn2, frame = time) %>%
           mutate(time = frame - 3) %>%
           mutate(power = as.factor(power))

df.yfp <- rbind(read.csv('results_3_04.csv')) %>%  # ,read.csv('results_17_03.csv')
          mutate(power = as.factor(power)) %>%
          mutate(cell = as.factor(cell)) %>%
          mutate(date = as.factor(date)) %>%
          mutate(rep = as.factor(rep)) %>%
          subset(select = c(ID, cell, date, rep, power, time, mask, delta, rel))

df.yfp.03 <- read.csv('results_17_03.csv') %>%
             mutate(power = as.factor(power)) %>%
             mutate(cell = as.factor(cell)) %>%
             mutate(date = as.factor(date)) %>%
             mutate(rep = as.factor(rep)) %>%
             subset(select = c(ID, cell, date, rep, power, time, mask, delta, rel))

##### RAW #####
# REPEATED STIMULATIONS
selected.power <- '50'
selected.mask <- 'up'

df.power.norm <- df.yfp %>%
                 filter(mask == selected.mask,
                        power == selected.power) %>%
                 subset(select = c(ID, cell, rep, time, delta)) %>%
                 group_by(cell, rep) %>%
                 mutate(delta_norm = delta / max(delta))

# 17.03.2021 data
df.power.norm.03 <- df.yfp.03 %>%
  filter(mask == selected.mask,
         power == selected.power) %>%
  subset(select = c(ID, cell, rep, time, delta)) %>%
  group_by(cell, rep) %>%
  mutate(delta_norm = delta / max(delta))

ggplot(data = df.power.norm.03) +
  geom_vline(xintercept = 0) +
  geom_line(aes(x=time, y= delta_norm, colour=rep),
            size=1.2) +
  facet_grid(rows = vars(cell)) +
  scale_x_continuous(name = 'Time (s)',
                    limits = c(-20, 90),
                   breaks = seq(-100, 100, 5)) +
  scale_y_continuous(name = 'Normalized ΔF/F0',
                     breaks = seq(-100, 100, 10)) +
  labs(title = sprintf("Mask %s, 405 nm power %s", selected.mask, selected.power),
       colour = 'Stimulation') +
  scale_color_jco() +
  scale_fill_jco() +
  theme_light()

# DROP BAD CELLS
df.yfp.drop <- df.yfp %>%
               filter(ID == 'cell1_01_3_04' |  # 405 nm 50%
                      ID == 'cell1_02_3_04' |
                      ID == 'cell1_03_3_04' |
                      ID == 'cell10_01_3_04'|
                      ID == 'cell4_01_3_04' |
                      ID == 'cell7_01_3_04' |
                      ID == 'cell2_01_3_04' |  # 405 nm 75%
                      ID == 'cell2_02_3_04' |
                      ID == 'cell2_03_3_04' |
                      ID == 'cell8_01_3_04' |
                      ID == 'cell8_02_3_04' |
                      ID == 'cell8_03_3_04' |
                      ID == 'cell3_01_3_04' |  # 405 nm 100%
                      ID == 'cell3_02_3_04' |
                      ID == 'cell6_01_3_04' |
                      ID == 'cell9_01_3_04')

df.yfp.drop.03 <- df.yfp.03 %>%
                  filter(ID == 'cell7_01_17_03'  |  # 405 nm 20%
                         ID == 'cell9_01_17_03'  |
                         ID == 'cell10_01_17_03' |  # 405 nm 50%
                         ID == 'cell12_01_17_03' |
                         ID == 'cell13_01_17_03' |
                         ID == 'cell20_01_17_03') 
# ALL GOOD CELLS               
df.yfp.good <- rbind(df.yfp.drop, df.yfp.drop.03)        

# DOSE DEP
start.time <- 5
end.time <- 95
decim.num <- 10

w.fun <- function(x, a, b){(exp(-x/a)+b)/max(exp(-x/a)+b)}
a <- 30
b <- 5

ggscatter(data.frame(y = w.fun(seq(0, 100, 1), a, b), x = seq(0, 100, 1)),
          x = 'x', y = 'y')

time.slice.yfp <- df.yfp.good %>%
                  filter(power == selected.power,
                         mask == selected.mask,
                         time>=start.time,
                         time<=end.time) %>%
                  subset(select = c(ID, time, delta)) %>%
                  group_by(ID) %>%
                  mutate(norm = delta/max(delta)) %>%
                  mutate(dec = HMA(norm, n = 10)) %>%  # w = w.fun(seq(0, length(norm)-1), a, b))
                  rename(yfp = delta)
time.slice.fluo <- df.fluo %>%
                   filter(power == selected.power,
                   time>=start.time,
                   time<=end.time) %>%
                   subset(select = c(mean, time)) %>%
                   rename(fluo = mean)

dose.dep.norm <- left_join(time.slice.yfp, time.slice.fluo,
                           by = c('time')) %>%
                 subset(select = c(ID, yfp, norm, dec, fluo))

ggplot() +
  geom_line(data = dose.dep.norm, 
            aes(x = fluo, y = dec, colour = ID),
            size = 1) +
  geom_line(data = dose.dep.norm, 
            aes(x = fluo, y = norm, colour = ID),
            size = 0.2,
            linetype = 5) +
  scale_x_continuous(name = 'Fluo-4 ΔF/F0',
                     limits = c(0.5, 2),
                     breaks = seq(0, 10, .1)) +
  scale_y_continuous(name = 'Normalized HPCA-YFP ΔF/F0',
                     limits = c(-0.5, 1),
                     breaks = seq(-1, 1, .1)) +
  labs(title = sprintf("Dose dep, 405 nm power %s", selected.power),
       colour = 'Cell') +
  scale_color_jco() +
  scale_fill_jco() +
  theme_minimal()

# ggsave(sprintf('dose_dep_raw_power_%s.png', selected.power),
#        width = 30, height = 15, units = 'cm')

##### MEAN #####
# TIME SERIES
selected.mask <- 'up'
coeff <- 7

df.mean.yfp <- subset(df.yfp, select = c(power, time, mask, rel)) %>%
               group_by(time, power, mask) %>%
               summarise_all(list(mean, sd), na.rm = TRUE) %>%
               rename(mean = fn1, sd = fn2)

  ggplot() +
  geom_line(data = filter(df.mean.yfp, mask == selected.mask),
            aes(x = time, y = mean, colour = power, fill = power),
            size = 1) +
  geom_ribbon(data = filter(df.mean.yfp, mask == selected.mask),
              aes(x = time,
                  ymin = mean-sd,
                  ymax = mean+sd,
                  fill = as.factor(power),
                  colour = as.factor(power)),
              alpha=.2,
              size=0) +
  geom_line(data = subset(df.fluo, power == c('20', '50')),
            aes(x = time, y = mean/coeff, colour = power, fill = power),
            size = 1,
            linetype = 5,
            alpha = 0.6) +
  scale_x_continuous(name = 'Time (s)',
                     limits = c(-20, 90),
                     breaks = seq(-100, 100, 5)) +
  scale_y_continuous(name = 'HPCA-YFP ΔF/F0',
                     breaks = seq(-100, 100, 0.05),
                     sec.axis = sec_axis(trans = ~.*coeff,
                                       name = 'Fluo-4 ΔF/F0',
                                       breaks = seq(-100, 100, 0.5))) +
  labs(title = sprintf("Mask %s mean", selected.mask),
       colour = '405 nm power (%)',
       fill = '405 nm power (%)') +
  scale_color_jco() +
  scale_fill_jco() +
  theme_minimal()

# ggsave(sprintf('mean_mask_%s.png', selected.mask, selected.power),
#        width = 30, height = 15, units = 'cm')


# DOSE DEP
start.time <- 4
end.time <- 100

time.slice.yfp <- df.mean.yfp %>%
                  filter(mask == selected.mask,
                         time>=start.time,
                         time<=end.time) %>%
                  subset(select = c(power, mean, time)) %>%
                  rename(mean_yfp = mean)
time.slice.fluo <- df.fluo %>%
                   filter(time>=start.time & time<=end.time) %>%
                   subset(select = c(power, mean, time)) %>%
                   rename(mean_fluo = mean)
dose.dep.mean <- left_join(time.slice.yfp, time.slice.fluo,
                            by = c('power', 'time')) %>%
                 subset(select = c(power, mean_yfp, mean_fluo))
dose.dep.rel <- dose.dep.mean

write.csv(file = 'DR_50_17_03.csv', dose.dep.rel)

ggplot(dose.dep.rel, aes(x = mean_fluo, y = mean_yfp, colour = as.factor(power))) +
  geom_line(size = 1) +
  scale_x_continuous(name = 'Fluo-4 ΔF/F0',
                     limits = c(0.2, 2.1),
                     breaks = seq(-100, 100, .1)) +
  scale_y_continuous(name = 'HPCA-YFP ΔF/F0',
                     limits = c(0, 0.35),
                     breaks = seq(-100, 100, 0.05)) +
  labs(title = "Dose dep mean",
       color = '405 nm power (%)') + 
  scale_color_jco() +
  scale_fill_jco() +
  theme_minimal()

# ggsave(sprintf('dose_dep_mean.png', selected.power),
#        width = 30, height = 15, units = 'cm')


