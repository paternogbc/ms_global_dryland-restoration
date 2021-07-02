library(tidyverse)
library(glmmTMB)
library(ggeffects)
library(effects)
library(patchwork)
library(cowplot)
library(extrafont)
library(ggnewscale)
loadfonts(device = "win")
source("R/ZZZ_functions.R")

q3 <- readRDS("outputs/models/revision/Q3_fitted_model.RDs")
dc <- readRDS(file = "outputs/data/revision/Q3_used_data.Rds")
dc_sp <- read_csv("outputs/tables/revision/Supp_species_average_success_full.csv")

# Predicted marginal effects------
mnc1  <- ggeffect(model = q3, terms = c("tsrLZ[all]", "lifeform2"))
mnc2  <- ggeffect(model = q3, terms = "weed_control")
mnc3  <- ggeffect(model = q3, terms = c("seededrateLZ[all]", "aridityLZ[-2, -1, 0, 1, 2]"))
mnc4  <- ggeffect(model = q3, terms = "seededrateLZ[all]")
mnc5  <- ggeffect(model = q3, terms = "aridityLZ[all]")
mnc6  <- ggeffect(model = q3, terms = "seedmassLZ[all]")


# Time : lifeform-----
# Raw data
preds <- as.data.frame(mnc1)
Xlabs <- unzt(formula = log(tsr) ~ tsrLZ, data = dc, xu = c(10, 50, 100, 200, 300))

# Color-blind friendly palette
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

gq3life <-
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
  scale_y_continuous(labels = seq(0, 1, .2),
                     breaks = seq(0, 1, .2)) +
  theme_classic(base_size = 14) +
  theme(legend.position = c(0.8, 0.7),
        legend.background = element_blank()) +
  theme(axis.title = element_text(size = 16)) +
  theme(text=element_text(family="Garamond")) + 
  labs(x = expression(paste("Time since seeding (weeks)")),
       y = "Species success"); gq3life

# Seed rate ---------------------------------------------------------------
preds <- as.data.frame(mnc4)
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
  scale_y_continuous(labels = seq(0, 1, .2),
                     breaks = seq(0, 1, .2)) +
  scale_size_area() +
  theme_classic(base_size = 14) +
  theme(axis.title = element_text(size = 16)) +
  theme(text=element_text(family="Garamond")) + 
  labs(x = expression(paste("Seed rate (", m ^ 2, ")")),
       y = "Species success"); gq3b

# Aridity -----------------------------------------------------------------
preds <- as.data.frame(mnc5)
Xlabs <- unzt(formula = log(aridity) ~ aridityLZ, data = dc, xu = c(0.05 ,0.1, .2, .4, .65))

gq3aridity <- 
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
  scale_y_continuous(labels = seq(0, 1, .2),
                     breaks = seq(0, 1, .2)) +
  scale_size_area() +
  theme_classic(base_size = 14) +
  theme(axis.title = element_text(size = 16)) +
  theme(text=element_text(family="Garamond")) + 
  labs(x = expression(paste("Aridity index")),
       y = "Species success"); gq3aridity

# Seed mass -----------------------------------------------------------------
preds <- as.data.frame(mnc6)
Xlabs <- unzt(formula = log(seedmass) ~ seedmassLZ, data = dc, xu = c(0.1, 1, 10, 100, 1000))

gq3seed <- 
  ggplot() + 
  geom_point(data =  dc %>% select(predicted = Ps, x = seedmassLZ, n = n),
             aes(y = predicted, x = x, size = n, color = predicted),
             alpha = .5, color = gray(.75), show.legend = FALSE) +
  geom_ribbon(data = preds, 
              aes(x = x, y = predicted, ymax = conf.high, ymin = conf.low),
              fill = "chartreuse4", alpha = .5) +
  geom_line(data = preds, 
            aes(x = x, y = predicted), color = "darkgreen", size = .8) +
  scale_x_continuous(
    breaks = Xlabs$xz,
    labels = Xlabs$xu,
    limits = c(-3, 5)
  ) +  
  scale_y_continuous(labels = seq(0, 1, .2),
                     breaks = seq(0, 1, .2)) +
  scale_size_area() +
  theme_classic(base_size = 14) +
  theme(axis.title = element_text(size = 16)) +
  theme(text=element_text(family="Garamond")) + 
  labs(x = expression(paste("Seed mass (1000 seeds / gram)")),
       y = "Species success"); gq3seed

# Weed control-----
# Raw data
preds <- as.data.frame(mnc2)

gq3weed <-
  ggplot() +
  geom_jitter(data =  dc %>% select(predicted = Ps, x = weed_control, n = n),
              aes(y = predicted, x = x, size = n),
              alpha = .05, show.legend = FALSE, width = .05, color = gray(.7)) +
  geom_pointrange(data = preds, size = 1.5, color = gray(.2),
                  aes(x = x, y = predicted, ymax = conf.high, ymin = conf.low)) +
  scale_size_area() +
  scale_y_continuous(breaks = seq(0, 1, .2), 
                     labels = seq(0, 1, .2)) +
  
  theme_classic(base_size = 14) +
  theme(axis.title = element_text(size = 16)) +
  theme(text=element_text(family="Garamond")) + 
  labs(x = "Weed control",
       y = "Species success"); gq3weed

# aridity : seed rate-----
preds <- as.data.frame(mnc3)
Xlabs1 <- unzt(formula = log(seededrate) ~ seededrateLZ, data = dc, xu = c(1, 10, 100, 1000, 10000))
Xlabs2 <- unzt(formula = log(aridity) ~ aridityLZ, data = dc, xu = c(0.05, 0.1, .2, .4, .65))


dc %>% 
  mutate(group = cut(aridity, breaks = c(-Inf, .1, .2, .4, .65, Inf))) -> dc2
table(dc2$group)

gq3seed_aridity <- 
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
                     labels = c("0.05", "0.10", "0.20", "0.40", "0.65"),
                     aesthetics = c("color", "fill"),
                     name = "Aridity index",
                     guide = guide_legend(reverse = TRUE)) +
  scale_x_continuous(breaks = Xlabs1$xz, 
                     labels = Xlabs1$xu, 
                     limits = c(-2.7, 3)) +
  scale_y_continuous(labels = seq(0, 1, .2),
                     breaks = seq(0, 1, .2)) +
  theme_classic(base_size = 14) +
  theme(legend.position = c(0.2, 0.7)) +
  theme(axis.title = element_text(size = 16)) +
  theme(text=element_text(family="Garamond")) + 
  theme(legend.title = element_text(size = 12),
        legend.background = element_blank()) +
  labs(x = expression(paste("Seed rate (seed / ", m ^ 2, ")")),
       y = "Species success"); gq3seed_aridity


# Genus level -------------------------------------------------------------
# Top and Lower 15 genus with at least 5 species per genus
dc_sp %>% 
  group_by(genus) %>% 
  summarise(n = n()) %>% 
  filter(n > 3) %>% 
  pull(genus) -> top_genus

dc_sp %>% 
  filter(genus %in% top_genus) %>% 
  group_by(genus) %>% 
  summarise(ps_gr = mean(Success),
            sd_gr = plotrix::std.error(Success),
            n = n()) %>% 
  arrange(desc(ps_gr)) -> dc_genus

# Plot Genus
gggenus <-
  ggplot(dc_sp %>% filter(genus %in% top_genus), 
         aes(y = Success, x = reorder(genus, -Success))) +
  geom_jitter(show.legend = F, width = .2, color = "lightgreen") +
  geom_boxplot(show.legend = F, alpha = .15, fill = "lightgreen" ) +
  stat_summary(color = "black", show.legend = FALSE, size = .2) +
  coord_flip() +
  scale_y_continuous(labels = seq(0, 1, .2),
                     breaks = seq(0, 1, .2)) +
  theme_classic(base_size = 14) +
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 9.5)) +
  theme(text=element_text(family="Arial")) + 
  theme(legend.title = element_text(size = 12),
        legend.background = element_blank()) +
  labs(x = "Genus",
       y = "Species success"); gggenus

# Fig3: Merge all plots---------------------------------------------------------------

gq3 <- 
  gq3seed_aridity + theme(text=element_text(family="Garamond")) +
  gq3aridity + theme(text=element_text(family="Garamond")) +
  gq3weed + theme(text=element_text(family="Garamond")) + 
  gq3life + theme(text=element_text(family="Garamond")) + 
  gggenus + theme(text=element_text(family="Garamond")) +
  gq3seed + theme(text=element_text(family="Garamond")) + 
  plot_annotation(tag_levels = "a")

# Save plot
ggsave(
  plot = gq3,
  filename = "outputs/figures/revision/Figure_3.png",
  width = 12,
  height = 8,
)
