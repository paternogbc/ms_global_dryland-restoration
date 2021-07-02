library(tidyverse)
library(here)
library(patchwork)
library(ggridges)
library(extrafont)
loadfonts()
rd  <- readRDS(here("data","temp","merged_transformed_gazp_data.RDS"))
stu <- openxlsx::read.xlsx(xlsxFile = here("data", "raw", "GAZP_Database-master", "projectdata.xlsx"), sheet = "study")
rd <- left_join(rd, stu)

# Plot 1- Raw Success-----------------------------------------------------------
rd %>% 
  filter(seededrate > 0 
         & aridity <= 0.65)  %>%
  filter(!is.na(genus)) %>%   # Only species identified into at least the genus level
  group_by(region, projectid, siteid,  disturbance, treatmentid, speciesid, tsr) %>% 
  summarise(n = n(),
            succ  = sum(occ == 1, na.rm = TRUE),
            fail  = sum(occ == 0, na.rm = TRUE),
            ps    = succ / n) %>% 
  ungroup() -> d_raw

# Success of projects in terms of % of species-------------------------
d_raw %>%
  group_by(region, projectid, siteid, disturbance, treatmentid, speciesid) %>%
  summarise(obs       = sum(n),
            n_succ    = sum(succ),
            n_fail    = sum(fail),
            ps        = n_succ / obs) %>%
  ungroup() -> d_raw_tsr

# Success by treatment (number of species with ps > 0)
d_trt <-
  d_raw_tsr %>%
  group_by(region, projectid, siteid, disturbance, treatmentid) %>%
  summarise(
    obs         = sum(obs),
    n_succ      = sum(n_succ),
    n_fail      = sum(n_fail),
    n_sp_succ   = sum(ps > 0),
    n_sp_fail   = sum(ps == 0),
    n_sp        = n_distinct(speciesid),
    ps_avg      = mean(ps),
    prop_sp     = n_sp_succ / n_sp
  ) %>% 
  ungroup()

# Add categories for number of seeded species
d_trt <- 
  d_trt %>% 
  mutate(
    cat_n = cut(n_sp, breaks = c(-Inf, 10,25,40, Inf), 
                labels = c("1-10 \n n = 1112","10-25 \n n = 450",
                           "25-40 \n n = 46", "> 40 \n n = 42")),
    cat = cut(n_sp, breaks = c(-Inf, 10,25,40, Inf), 
              labels = c("1-10","10-25",
                         "25-40", "> 40")))

sum(d_trt$prop_sp == 1)/nrow(d_trt)
sum(d_trt$prop_sp == 0)/nrow(d_trt)

# Figure 1. B-------------------------------------------------------------------
# proportion of successful species in a treatment --------------------------
d_trt$disturbance <- d_trt$disturbance %>%  
  fct_collapse(grazing = c("grazing", "grazing|fire"))

dist_colors <- c(Cropping = "#984EA3", Invasion = "#4DAF4A", Clearing = "#FF7E00",
                 Erosion	= "#999999", Fire	= "#377EB8", Oil_and_gas =	"#FBFB30",
                 Overgrazing	= "#A65627", Mining	= "#E4191A", Roads = "#F881BF")

g2b <-
  ggplot(d_trt %>% drop_na(disturbance),
         aes(x = prop_sp * 100, fill = fct_reorder(disturbance, prop_sp))) +
  geom_histogram(alpha = .75, bins = 20, 
                 color = gray(.2), size = .1,
                 show.legend = FALSE) +
  scale_y_continuous(limits = c(0,500)) +
  scale_fill_manual(values = as.character(dist_colors[c(4, 6, 2, 5, 8, 1, 9, 3, 7)]),
                    name = "") +
  theme_classic(base_size = 16) +
  labs(y = "Frequency", x = "Present species (%)") +
  theme(text=element_text(family="Garamond")); g2b

# Figure 1. C-------------------------------------------------------------------
# Average success against species richness v2
g2c <-
  ggplot(d_trt, aes(x =n_sp_succ, y = cat)) +
  geom_density_ridges(show.legend = FALSE, rel_min_height = 0.005,
                      fill = "#addd8e", 
                      color = "darkgreen", 
                      alpha = .75,
                      quantile_lines = TRUE, quantiles = 2,
                      jittered_points = FALSE, scale = .9) +
  geom_jitter(width = .2, height = .02, size = .05, color = gray(.2), show.legend = F) +
  annotate(geom = "text", x = 42.5, y = 4.3, label = "N = 42", hjust = 0) +
  annotate(geom = "text", x = 42.5, y = 3.3, label = "N = 46", hjust = 0) +
  annotate(geom = "text", x = 42.5, y = 2.3, label = "N = 355", hjust = 0) +
  annotate(geom = "text", x = 42.5, y = 1.3, label = "N = 839", hjust = 0) +
  stat_summary(show.legend = FALSE, geom = "pointrange", fun = median, color = "red") +
  scale_y_discrete(expand = expand_scale(add = c(0.2, .75))) +
  #scale_x_continuous(expand = c(0.05, 0.05)) +
  scale_x_continuous(limits = c(0,50), breaks = seq(0, 50, 10)) +
  scale_fill_distiller(palette = 15, direction = 1) +
  theme_bw(base_size = 16) +
  labs(y = "Number of seeded species", 
       x = "Number of present species",
       title = "") +
  theme(text=element_text(family="Garamond"));g2c

# Merge both plots
ggout <- g2b+labs(tag = "B") + g2c +  labs(tag = "C") +   plot_layout(widths = c(1.5, 1))
ggout
ggsave(ggout, filename = "outputs/figures/revision/Figure_1_BC.png", width = 11.8, height = 5)
