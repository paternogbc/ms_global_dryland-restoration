library(tidyverse)
library(glmmTMB)
library(ggeffects)
library(effects)
library(patchwork)
library(cowplot)
library(extrafont)
library(ggnewscale)
loadfonts()
source("R/ZZZ_functions.R")

q3 <- readRDS("outputs/models/Q3_fitted_model.RDs")
dc <- readRDS(file = "outputs/data/Q3_used_data.Rds")

mnc1  <- ggeffect(model = q3, terms = c("tsrLZ[all]", "lifeform2"))
mnc2  <- ggeffect(model = q3, terms = "weed_control")
mnc3  <- ggeffect(model = q3, terms = c("seededrateLZ[all]", "aridityLZ[-2, -1, 0, 1, 2]"))
mnc4  <- ggeffect(model = q3, terms = "seededrateLZ[all]")
mnc5  <- ggeffect(model = q3, terms = "aridityLZ[all]")

# Time vs lifeform-----
# Raw data
preds <- as.data.frame(mnc1)
Xlabs <- unzt(formula = log(tsr) ~ tsrLZ, data = dc, xu = c(10, 50, 100, 200, 300))

# Color-blind friendly palette
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

gq3a <-
  ggplot() +
  geom_jitter(data =  dc %>% select(predicted = Ps, x = tsrLZ, n = n, group = lifeform2),
              aes(y = predicted, x = x, size = n),
              alpha = .1, show.legend = FALSE, color = gray(.7)) +
  scale_size_area() +
  geom_ribbon(data = preds, 
              aes(x = x, y = predicted, ymax = conf.high, ymin = conf.low,
                  fill = group), alpha = .25) +
  geom_line(data = preds, 
            aes(x = x, y = predicted, color = group,), size = .8) +
  scale_color_manual(values = cbp1, 
                     name = "Lifeform", aesthetics = c("color", "fill")) +
  scale_x_continuous(breaks = Xlabs$xz, 
                     labels = Xlabs$xu) +
  scale_y_continuous(labels = seq(0, 100, 20),
                     breaks = seq(0, 1, .2)) +
  theme_classic(base_size = 14) +
  theme(legend.position = c(0.8, 0.7),
        legend.background = element_blank()) +
  theme(axis.title = element_text(size = 16)) +
  theme(text=element_text(family="Garamond")) + 
  labs(x = expression(paste("Time since seeding (weeks)")),
       y = "Success (%)"); gq3a

# Seed rate ---------------------------------------------------------------
preds <- as.data.frame(mnc4)

# Seed rate plot Q2------------
Xlabs <- unzt(formula = log(seededrate) ~ seededrateLZ, data = dc, xu = c(1, 10, 100, 1000, 10000))

gq3b <-
  ggplot() + 
  geom_jitter(data =  dc %>% select(predicted = Ps, x = seededrateLZ, n = n),
              aes(y = predicted, x = x, size = n, color = predicted),
              alpha = .5, color = gray(.75), show.legend = FALSE) +
  geom_ribbon(data = preds, 
              aes(x = x, y = predicted, ymax = conf.high, ymin = conf.low),
              fill = "#b8e186", alpha = .7) +
  geom_line(data = preds, 
            aes(x = x, y = predicted), color = "darkgreen", size = 1.2) +
  scale_x_continuous(
    breaks = Xlabs$xz,
    labels = Xlabs$xu,
    limits = c(-2.7, 3)
  ) +  
  scale_y_continuous(labels = seq(0, 100, 20),
                     breaks = seq(0, 1, .2)) +
  scale_size_area() +
  theme_classic(base_size = 14) +
  theme(axis.title = element_text(size = 16)) +
  theme(text=element_text(family="Garamond")) + 
  labs(x = expression(paste("Seed rate (", m ^ 2, ")")),
       y = "Success (%)"); gq3b


# Aridity -----------------------------------------------------------------
preds <- as.data.frame(mnc5)
Xlabs <- unzt(formula = log(aridity) ~ aridityLZ, data = dc, xu = c(0.05 ,0.1, .2, .4, .7))

gq3c <- 
  ggplot() + 
  geom_point(data =  dc %>% select(predicted = Ps, x = aridityLZ, n = n),
             aes(y = predicted, x = x, size = n, color = predicted),
             alpha = .5, color = gray(.75), show.legend = FALSE) +
  geom_ribbon(data = preds, 
              aes(x = x, y = predicted, ymax = conf.high, ymin = conf.low),
              fill = "orange", alpha = .5) +
  geom_line(data = preds, 
            aes(x = x, y = predicted), color = "tomato", size = .8) +
  scale_x_continuous(
    breaks = Xlabs$xz,
    labels = Xlabs$xu,
    limits = c(-2, 2)
  ) +  
  scale_y_continuous(labels = seq(0, 100, 20),
                     breaks = seq(0, 1, .2)) +
  scale_size_area() +
  theme_classic(base_size = 14) +
  theme(axis.title = element_text(size = 16)) +
  theme(text=element_text(family="Garamond")) + 
  labs(x = expression(paste("Aridity index")),
       y = "Success (%)"); gq3c

# Weed control-----
# Raw data
preds <- as.data.frame(mnc2)

gq3d <-
  ggplot() +
  geom_jitter(data =  dc %>% select(predicted = Ps, x = weed_control, n = n),
              aes(y = predicted, x = x, size = n),
              alpha = .05, show.legend = FALSE, width = .05, color = gray(.7)) +
  geom_pointrange(data = preds, size = 1.5, color = gray(.2),
                  aes(x = x, y = predicted, ymax = conf.high, ymin = conf.low)) +
  scale_size_area() +
  scale_y_continuous(breaks = seq(0, 1, .2), 
                     labels = 100 * seq(0, 1, .2)) +
  
  theme_classic(base_size = 14) +
  theme(axis.title = element_text(size = 16)) +
  theme(text=element_text(family="Garamond")) + 
  labs(x = "Weed control",
       y = "Success (%)"); gq3d

# Interaction aridity seed rate
preds <- as.data.frame(mnc3)
Xlabs1 <- unzt(formula = log(seededrate) ~ seededrateLZ, data = dc, xu = c(1, 10, 100, 1000, 10000))
Xlabs2 <- unzt(formula = log(aridity) ~ aridityLZ, data = dc, xu = c(0.05, 0.1, .2, .4, .7))


dc %>% 
  mutate(group = cut(aridity, breaks = c(-Inf, .1, .2, .4, .7, Inf))) -> dc2
table(dc2$group)

gq3b2 <- 
  ggplot(data =  dc %>% select(predicted = Ps, x = seededrateLZ, n = n, group = aridity)) +
  geom_point(aes(y = predicted, x = x, size = n),
              alpha = .2, show.legend = FALSE, color = "gray") +
  scale_size_area(guide = 'none') +
  scale_color_gradient(low = "red", high = "blue") +
  new_scale_color() +
  geom_ribbon(data = preds, 
              aes(x = x, y = predicted, ymax = conf.high, ymin = conf.low,
                  fill = group), alpha = .25, show.legend = TRUE) +
  geom_line(data = preds, show.legend = TRUE, 
            aes(x = x, y = predicted, color = group,), size = .8) +
  scale_color_manual(values = c("#FFE526", "#A0D93C", "#1DA287", "#365C8D", "#430052"), 
                     labels = c("0.05", "0.10", "0.20", "0.40", "0.70"),
                     aesthetics = c("color", "fill"),
                     name = "Aridity index") +
  scale_x_continuous(breaks = Xlabs1$xz, 
                     labels = Xlabs1$xu, 
                     limits = c(-2.7, 3)) +
  scale_y_continuous(labels = seq(0, 100, 20),
                     breaks = seq(0, 1, .2)) +
  theme_classic(base_size = 14) +
  theme(legend.position = c(0.15, 0.7)) +
  theme(axis.title = element_text(size = 16)) +
  theme(text=element_text(family="Garamond")) + 
  theme(legend.title = element_text(size = 12),
        legend.background = element_blank()) +
  labs(x = expression(paste("Seed rate (seed / ", m ^ 2, ")")),
       y = "Success (%)"); gq3b2

gq3 <- gq3b2 + theme(text=element_text(family="Garamond")) +
  gq3c + theme(text=element_text(family="Garamond")) +
  gq3a + theme(text=element_text(family="Garamond")) + 
  gq3d + theme(text=element_text(family="Garamond")) + plot_annotation(tag_levels = "A")
gq3

# Save plot
ggsave(
  plot = gq3,
  filename = "outputs/figures/Figure_3_ABCD.png",
  width = 11,
  height = 8,
)