# Experiments with Hill equation.
# Copyright © 2021 Borys Olifirov
#
#       [C]^n
# Y = -----------
#    Ka^n + [C]^n
#
# Y - fraction of occupied target protein
# [C] - ligand (Ca2+) concentration
# n - Hill coeficient
# Ka - ligand concentration producing half occupation
#
# Kd = Ka^n, Kd - dissociation constant

require(ggplot2)
require(ggsci)


Ka <- 0.5
n <- 2
obs <- 5  # number of observations
hill.sd <- 0.01  # additive gaussian noise sd
hill <- function(x){((x^n)/(Ka^n+x^n))}  # Hill equation function

model.hill <- data.frame(c = seq(0, 1, 1/obs)) %>%
  mutate(y = hill(c)) %>%
  mutate(c_noise = c + rnorm(obs+1, mean = 0, sd = hill.sd)) %>%
  mutate(y_noise = y + rnorm(obs+1, mean = 0, sd = hill.sd))

# Dose-responce plot building
ggplot(model.hill) +
  stat_function(data = data.frame(x = c(0, 1)),
                aes(x),
                fun = function(x){((x^n)/(Ka^n+x^n))},
                color = 'red') +
  geom_point(aes(x = c_noise, y = y_noise)) +
  scale_x_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, 0.1)) +
  labs(title = sprintf('Dose-responce plot, Ka=%s, n=%s', Ka, n),
       subtitle = sprintf('Red - theoretical curve, dots - model data + gaussian noise (sd=%s)', hill.sd), 
       x = 'C',
       y = 'Y') +
  theme_minimal()

# Hill plot building
ggplot(model.hill) + 
  geom_point(aes(x = log10(c_noise), y = log10(y_noise/(1-y_noise)))) +
  geom_line(aes(x = log10(c), y = log10(y/(1-y))), color = 'red') +
  stat_smooth(aes(x = log10(c_noise), y = log10(y_noise/(1-y_noise))),
              method = 'lm') +
  scale_x_continuous(limits = c(-2, 0),
                     breaks = seq(-10, 10, 0.1)) +
  scale_y_continuous(limits = c(-3.5, 1),
                     breaks = seq(-10, 10, 0.5)) +
  labs(title = sprintf('Hill plot, Ka = %s, n = %s', Ka, n),
       subtitle = 'Red - model data without noise, blue - linear model', 
       x = 'log(C)',
       y = 'log(Y/1-Y)') +
  theme_minimal()

ggplot() +
  stat_function(data = data.frame(x = c(0, c)), aes(x),
                fun = function(x){(x^n)/(Ka^n+x^n)}) +
  scale_x_continuous(limits = c(0, c),
                     breaks = seq(0, c, c/10)) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, c, c/10)) +
  labs(title = sprintf('Ka = %s, n = %s', Ka, n),
       x = 'Ligand concentration',
       y = 'Binded target') +
  theme_minimal()

# linear model building
hill.mod <- lm(log10(y_noise/(1-y_noise)) ~ log10(c_noise), data = model.hill)
summary(hill.mod)

