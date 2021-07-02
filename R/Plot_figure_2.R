# Packages-----------------------------------------------------------
library(tidyverse)
library(here)
library(patchwork)
library(ggridges)
library(extrafont)
library(glmmTMB)
library(ggeffects)
loadfonts()

q2_pairs  <- readRDS(here("outputs","data/revision", "Q2_pairs_seeded_vs_unseeded.RDs"))
q2_model  <- readRDS(here("outputs","models/revision", "Q2_fitted_model.RDs"))

rd  <- readRDS(here("data","temp","merged_transformed_gazp_data.RDS"))
stu <- openxlsx::read.xlsx(xlsxFile = here("data", "raw", "GAZP_Database-master", "projectdata.xlsx"), sheet = "study")
rd <- left_join(rd, stu)

# Raw Success-----------------------------------------------------------
rd %>% 
  filter(seededrate > 0)  %>%
  filter(!is.na(genus)) %>%   # Only species identified into at least the genus level
  group_by(region, projectid, siteid,  disturbance, treatmentid, speciesid, tsr) %>% 
  summarise(n = n(),
            succ  = sum(occ == 1, na.rm = TRUE),
            fail  = sum(occ == 0, na.rm = TRUE),
            ps    = succ / n) %>% 
  ungroup() -> d_raw

sum(d_raw$ps == 1) / nrow(d_raw) # Percentage of species that failed completely
sum(d_raw$ps == 0) / nrow(d_raw) # Percentage of species that failed completely
mean(d_raw$ps)                   # Average success 

# Figure 2. A-------------------------------------------------------------------
g2_A <- 
  ggplot(d_raw, aes(x = ps*100)) +
  geom_histogram(fill = "#addd8e", color = "darkgreen", alpha = .75) +
  theme_classic(base_size = 16) +
  labs(y = "Frequency", x = "Present (%)",
       subtitle = "") +
  theme(text=element_text(family="Garamond")); g2_A

# Figure 2. B-------------------------------------------------------------------
# Predictions
pmd <- ggpredict(q2_model, terms = c("seed"))
pmd
plot(pmd)

# Paired plot
g2_B <-
  ggplot(q2_pairs, aes(y = Ps*100, x = seed, group = id, color = sign)) +
  geom_line(show.legend = F, alpha = .3, size = 1) +
  geom_point(alpha = .5, show.legend = FALSE, size = 2) +
  scale_color_manual(values = c("#ca0020", "gray", "#0571b0")) +
  scale_x_discrete(limits = c("unseeded", "seeded")) +
  labs(x = "Treatment",
       y = "Present (%)",
       title = "") +
  theme_classic(base_size = 16) +
  theme(axis.title = element_text(size = 16)) +
  theme(text=element_text(family="Garamond")); g2_B

# both together
g2 <- 
  g2_A + theme(plot.subtitle = element_blank()) + 
  g2_B + 
  plot_annotation(tag_levels = "A")

ggsave(
  plot = g2,
  filename = "outputs/figures/revision/Figure_2_AB.png",
  width = 10, height = 6
)
