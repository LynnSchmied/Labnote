# HPCA-EYFP data analysis, uncaging with different stimulation power
# Copyright © 2021 Borys Olifirov

require(dplyr)
require(tidyr)
require(purrr)
require(ggplot2)

# require(magrittr)

setwd('/home/astria/bio/note/projects/0_2020_hek_transient/experiments/1_03_2021_hpca-yfp_analysis/data')

df.fluo <- read.csv('results_fluo.csv') %>%  # 20.10.2020, results of test uncaging, stimulation after 5th frame
           subset(select = c(power, time, int)) %>%
           group_by(time, power) %>%
           summarise_all(list(mean, sd), na.rm = TRUE) %>%
           rename(mean = fn1, sd = fn2, frame = time)
df.fluo$time <- as.numeric(df.fluo$frame) - 4

df.yfp <- rbind(read.csv('results_10.csv'), read.csv('results_20.csv'))
df.mask <- subset(df.yfp,
                  ((mask == 'up') & (cell == 'cell2_01')) |
                  ((mask == 'up') & (cell == 'cell7_01')))


coeff <- 10
ggplot() +
  geom_line(data=df.mask,
            aes(x=time, y=delta, colour=as.factor(cell)),
            size=1) +
  geom_line(data=subset(df.fluo, (power == 20)),
            aes(x=time, y=mean/coeff, colour=as.factor(power)),
            size=1,
            linetype=5) +
  scale_x_continuous(name='Time (s)',
                     limits=c(-20, 120),
                     breaks=seq(-100, 100, 10)) +
  scale_y_continuous(name='HPCA-YFP ΔF/F0',
                     sec.axis=sec_axis(trans=~.*coeff,
                                       name='Fluo-4 ΔF/F0',
                                       breaks=seq(-100, 100, 0.5))) +
  theme_minimal()


ggplot(df.sum, aes(x=time)) +
  geom_ribbon(aes(ymin=fn1-fn2,
                  ymax=fn1+fn2,
                  fill=as.factor(power),
                  colour=as.factor(power)),
              alpha=.3,
              size=.2) +
  geom_line(aes(y=fn1, colour=as.factor(power)),
            size=1) +
  geom_point(aes(y=fn1, colour=as.factor(power)),
             size=1.5) +
  geom_vline(xintercept = 4, size=.5, alpha = .5) +
  scale_x_continuous(limits = c(0, 120),
                     breaks = seq(0, 120, 20)) +
  labs(y ='ΔF/F0', x = 'Time (s)', title = 'Fluo-4 transient') +
  guides(fill=guide_legend(title='405 nm power (%)'),
         colour=guide_legend(title='405 nm power (%)')) +
  theme_minimal(base_size = 18)
