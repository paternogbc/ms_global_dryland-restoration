# Packages ---------------------------------------------------------------------
library(tidyverse)
library(glmmTMB)
library(sjPlot)
library(DHARMa)
library(car)
library(inspectdf)
library(patchwork)

# Own functions-----------------------------------------------------------------
source(file = "R/ZZZ_functions.R")

# Load Data --------------------------------------------------------------------
# Raw data
# Raw data
da <- readRDS(file = "data/processed/revision/clean_PS_data-no-pooling.RDs")

# Pool across replicates
# Should mirror major decisions in other analyses-------------------------------
## Density
dc <- da %>%
  filter(tsr < 312 & tsr > 8) %>%
  filter(seededrate < 10000 & seededrate > 0)  %>%
  group_by(projectid,
           siteid,
           tsr,
           treatmentid,
           speciesid) %>%
  summarise(
    seededrate        = unique(seededrate),
    aridity           = unique(aridity),
    lifeform2         = unique(lifeform2),
    weed_control      = unique(weed_control),
    seedmass          = unique(seedmass),
    region            = unique(region),
    n = n(),
    n_fail   = sum(occ == 0),
    n_succ   = sum(occ == 1),
    Ps       = n_succ / n
  ) %>%
  ungroup()

dc <- na.omit(dc)

# Transform and standardize variables
dc <- dc %>%
  mutate(
    seededrateL = log(seededrate),
    seededrateLZ = scale(seededrateL),
    seedmassL = log(seedmass),
    seedmassLZ = round(scale(seedmassL),1),
    tsrL = log(tsr),
    tsrLZ = scale(tsrL),
    aridityL = log(aridity),
    aridityLZ = scale(aridityL)
  )
glimpse(dc)
inspect_na(dc)

# Save data used in models---------------
saveRDS(object = dc, file = "outputs/data/revision/Q3_used_data.Rds")
saveRDS(object = dc[dc$region == "North America", ], file = "outputs/data/revision/Q3_NA_used_data.Rds")

# Variables description---------------------------------------------------------
# Life form
g1 <-
  ggplot(dc, aes(x = lifeform2)) +
  geom_bar(fill = "steelblue", alpha = .8) +
  scale_x_discrete(labels = c("annual", "p.forb", "p.grass", "woody")) +
  theme_bw(base_size = 16) +
  labs(y = "Frequency",
       x = "Life form")
g1

# Weed control
g2 <-
  ggplot(dc, aes(x = weed_control)) +
  geom_bar(fill = "steelblue", alpha = .8) +
  theme_bw(base_size = 16) +
  labs(y = "Frequency",
       x = "Weed control")
g2

# Aridity
g3 <-
  ggplot(dc, aes(x = aridity)) +
  geom_histogram(fill = "steelblue", alpha = .8) +
  theme_bw(base_size = 16) +
  # scale_x_continuous(trans = "log10") +
  labs(y = "Frequency",
       x = "Aridity")
g3

# Seeded rate
g4 <-
  ggplot(dc, aes(x = seededrate)) +
  geom_histogram(fill = "steelblue", alpha = .8) +
  theme_bw(base_size = 16) +
  scale_x_continuous(trans = "log10", breaks = c(1, 10, 100, 1000, 10000)) +
  labs(y = "Frequency",
       x = "Seed rate [log10]")
g4

# Seed mass
g7 <-
  ggplot(dc, aes(x = seedmass)) +
  geom_histogram(fill = "steelblue", alpha = .8) +
  theme_bw(base_size = 16) +
  scale_x_continuous(trans = "log10", 
                     breaks = c(0.01, 0.1, 0, 1, 10, 100, 1000),
                     labels = c(0.01,0.1, 0, 1, 10, 100, 1000)) +
  labs(y = "Frequency",
       x = "Seed mass [log10]")
g7

# Time since restoration
g5 <-
  ggplot(dc, aes(x = tsr)) +
  geom_histogram(fill = "steelblue", alpha = .8) +
  theme_bw(base_size = 16) +
  labs(y = "Frequency",
       x = "Weeks since restoration")
g5

# Ps
g6 <-
  ggplot(dc, aes(x = Ps)) +
  geom_histogram(alpha = 1, fill = "orange") +
  theme_bw(base_size = 16) +
  labs(y = "Frequency",
       x = "Species Success")

gall <- g6 + g1 + g2 + g3 + g4 + g5 + (g7) +
  plot_annotation(tag_levels = 'a') &
  theme(plot.title = element_text(size = 18),
        plot.subtitle = element_text(size = 16))
gall

ggsave(gall,
       filename = "outputs/figures/revision/Q3_supp_Fig_variables_description.png",
       width = 11,
       height = 9)

# 2. Fit models----------------------------------------------------------
# 2.1 Full model scaled------
full_mod_Z <-
  glmmTMB(
    Ps ~ aridityLZ * seededrateLZ + 
      weed_control  + seedmassLZ +
      tsrLZ * lifeform2 + 
      (1 |
         projectid / siteid / treatmentid) + (1 | speciesid),
    family = binomial(link = "logit"),
    ziformula = ~ 1,
    weights = n,
    data = dc
  )
summary(full_mod_Z)

# Model table------
tab_model(full_mod_Z, transform = NULL)
tab_model(full_mod_Z, transform = NULL, show.reflvl = TRUE, 
          file = "outputs/tables/revision/Q3_table_fitted_model.html")

# 4.2. Model residuals
res <- simulateResiduals(full_mod_Z)
plot(res)

plotResiduals(res)
testUniformity(res)     # test uniformity
testZeroInflation(res)  # Zero-inflation test
testDispersion(res)     # Overdispersion

# 4.3. Model residuals against predictors
plotResiduals(simulationOutput = res,  form = dc$seededrateLZ)
plotResiduals(simulationOutput = res,  form = dc$aridityLZ)
plotResiduals(simulationOutput = res,  form = dc$lifeform2)
plotResiduals(simulationOutput = res,  form = dc$weed_control)
plotResiduals(simulationOutput = res,  form = factor(dc$seedmassLZ))
plotResiduals(simulationOutput = res,  form = dc$tsrLZ)

# Random effects
plot_model(full_mod_Z, type = "diag")

# Save models-------------------------------------------------------------------
# Save model object-------------------------------------------------------------
saveRDS(object = full_mod_Z, file = "outputs/models/revision/Q3_fitted_model.RDs")
saveRDS(object = res, file = "outputs/models/revision/Q3_simulated_residuals.Rds")

# 3. Fit North America only models----------------------------------------------
# Variables description---------------------------------------------------------
# Life form
g1 <-
  ggplot(dc[dc$region == "North America", ], aes(x = lifeform2)) +
  geom_bar(fill = "steelblue", alpha = .8) +
  scale_x_discrete(labels = c("annual", "p.forb", "p.grass", "woody")) +
  theme_bw(base_size = 16) +
  labs(y = "Frequency",
       x = "Life form")
g1

# Weed control
g2 <-
  ggplot(dc[dc$region == "North America", ], aes(x = weed_control)) +
  geom_bar(fill = "steelblue", alpha = .8) +
  theme_bw(base_size = 16) +
  labs(y = "Frequency",
       x = "Weed control")
g2

# Aridity
g3 <-
  ggplot(dc[dc$region == "North America", ], aes(x = aridity)) +
  geom_histogram(fill = "steelblue", alpha = .8) +
  theme_bw(base_size = 16) +
  # scale_x_continuous(trans = "log10") +
  labs(y = "Frequency",
       x = "Aridity")
g3

# Seeded rate
g4 <-
  ggplot(dc[dc$region == "North America", ], aes(x = seededrate)) +
  geom_histogram(fill = "steelblue", alpha = .8) +
  theme_bw(base_size = 16) +
  scale_x_continuous(trans = "log10", breaks = c(1, 10, 100, 1000, 10000)) +
  labs(y = "Frequency",
       x = "Seed rate [log10]")
g4

# Seed mass
g7 <-
  ggplot(dc[dc$region == "North America", ], aes(x = seedmass)) +
  geom_histogram(fill = "steelblue", alpha = .8) +
  theme_bw(base_size = 16) +
  scale_x_continuous(trans = "log10", 
                     breaks = c(0.01, 0.1, 0, 1, 10, 100, 1000),
                     labels = c(0.01,0.1, 0, 1, 10, 100, 1000)) +
  labs(y = "Frequency",
       x = "Seed mass [log10]")
g7

# Time since restoration
g5 <-
  ggplot(dc[dc$region == "North America", ], aes(x = tsr)) +
  geom_histogram(fill = "steelblue", alpha = .8) +
  theme_bw(base_size = 16) +
  labs(y = "Frequency",
       x = "Weeks since restoration")
g5

# Ps
g6 <-
  ggplot(dc[dc$region == "North America", ], aes(x = Ps)) +
  geom_histogram(alpha = 1, fill = "orange") +
  theme_bw(base_size = 16) +
  labs(y = "Frequency",
       x = "Species Success")

gall <- g6 + g1 + g2 + g3 + g4 + g5 + (g7) +
  plot_annotation(tag_levels = 'a') &
  theme(plot.title = element_text(size = 18),
        plot.subtitle = element_text(size = 16))
gall

ggsave(gall,
       filename = "outputs/figures/revision/Q3_NA_supp_Fig_variables_description.png",
       width = 11,
       height = 9)

# 3.1 Full model scaled------
full_mod_na_Z <-
  glmmTMB(
    Ps ~ aridityLZ * seededrateLZ + 
      weed_control  + seedmassLZ +
      tsrLZ * lifeform2 + 
      (1 |
         projectid / siteid / treatmentid) + (1 | speciesid),
    family = binomial(link = "logit"),
    ziformula = ~ 1,
    weights = n,
    data = dc[dc$region == "North America", ]
  )
summary(full_mod_na_Z)

# Model table------
tab_model(full_mod_na_Z, transform = NULL)
tab_model(full_mod_na_Z, transform = NULL, show.reflvl = TRUE, 
          file = "outputs/tables/revision/Q3_NA_table_fitted_model.html")

# 4.2. Model residuals
res_na <- simulateResiduals(full_mod_na_Z)
plot(res_na)

plotResiduals(res_na)
testUniformity(res_na)     # test uniformity
testZeroInflation(res_na)  # Zero-inflation test
testDispersion(res_na)     # Overdispersion

# 4.3. Model residuals against predictors
plotResiduals(simulationOutput = res_na,  form = dc$seededrateLZ[dc$region == "North America"])
plotResiduals(simulationOutput = res_na,  form = dc$aridityLZ[dc$region == "North America"])
plotResiduals(simulationOutput = res_na,  form = dc$lifeform2[dc$region == "North America"])
plotResiduals(simulationOutput = res_na,  form = dc$weed_control[dc$region == "North America"])
plotResiduals(simulationOutput = res_na,  form = factor(dc$seedmassLZ[dc$region == "North America"]))
plotResiduals(simulationOutput = res_na,  form = dc$tsrLZ[dc$region == "North America"])

# Random effects
plot_model(full_mod_na_Z, type = "diag")

# Save models-------------------------------------------------------------------
# Save model object-------------------------------------------------------------
saveRDS(object = full_mod_na_Z, file = "outputs/models/revision/Q3_NA_fitted_model.RDs")
saveRDS(object = res_na, file = "outputs/models/revision/Q3_NA_simulated_residuals.Rds")
# END----