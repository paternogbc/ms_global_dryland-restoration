# Packages ---------------------------------------------------------------------
library(tidyverse)
library(glmmTMB)
library(sjPlot)
library(DHARMa)
library(car)
library(cowplot)
library(inspectdf)
library(patchwork)
library(ggeffects)

# Load data ---------------------------------------------------------------
dat <- readRDS("data/processed/revision/clean_PS_data-pooling_rep_early_time.RDs")

# Ensure factors and rm undesired levels + group_by-----------------------------
dat <- dat %>% 
  mutate(projectid         = as.factor(projectid),
         region            = as.factor(region),
         siteid            = as.factor(siteid),
         treatmentid       = as.factor(treatmentid),
         speciesid         = as.factor(speciesid),
         seed              = ifelse(test = seededrate > 0, 
                                    yes = "seeded", 
                                    no = "unseeded")) %>% 
  droplevels() %>% 
  select(projectid, 
         region,
         siteid, 
         treatmentid, 
         tsr, 
         speciesid,
         n, 
         seed, 
         seededrate, 
         n_fail, 
         n_succ, 
         Ps)

inspect_na(dat)
glimpse(dat)

# Remove NAs--------------------------------------------------------------------
dat <- na.omit(dat) %>% droplevels()
nrow(dat)

# Remove species only with life form identiy------------------------------------
# Remove unseeded species communities-------------------------------------------
# Remove species unseeded but with community data-------------------------------
dat %>% 
  group_by(projectid, speciesid) %>% 
  summarise(n = n(),
            n_unseed = sum(seededrate == 0),
            comm     = n == n_unseed) %>% 
  filter(comm == TRUE) %>% 
  ungroup() %>% 
  droplevels() -> comm_spp 
remo <- split(comm_spp, comm_spp$projectid)

# Remove species from dataset
dc <- dat
for (i in 1:length(remo)) {
  dc <- 
    dc[!(dc$projectid %in% as.character(remo[[i]]$projectid) &
           dc$speciesid %in% as.character(remo[[i]]$speciesid)), ]
}

# tsr < 150  & seed rate < 10 000
dc <- 
  dat %>% 
  filter(tsr < 156) %>% 
  filter(seededrate < 10000) %>% 
  droplevels()

inspect_cat(dc) %>% show_plot()
inspect_num(dc) %>% show_plot()

# Check control treatments within a project
dc %>% 
  group_by(projectid, treatmentid, siteid, speciesid) %>%
  filter(seed == "unseeded") -> controls

dc %>% 
  group_by(projectid, treatmentid, siteid, speciesid) %>%
  filter(seed == "seeded") -> seeded

# Which project/sites have seed and controls together?
dc %>% group_by(projectid, siteid, speciesid) %>% 
  summarise(cont = any(seed == "unseeded"),
            seed = any(seed == "seeded"),
            both = cont == TRUE & seed == TRUE) %>% 
  filter(both == TRUE) -> projects_with_controls

# Create categorical seed rate
dc$seed_c = cut(dc$seededrate, breaks = c(-Inf, 1, 10, 100, 1000, Inf), 
                labels = c("unseeded", "1-10", "10-100", "100-1000", "> 1000"))

# Filter data into these projects
dc %>% 
  filter(projectid %in% projects_with_controls$projectid & 
           siteid %in% projects_with_controls$siteid & 
           speciesid %in% projects_with_controls$speciesid) %>% 
  droplevels() -> dc_clean

# get pairs
pairs_ps <- dc %>% filter(seededrate < 0)
projs <- unique(projects_with_controls$projectid)

set.seed(543210)
for (j in projs) {
  projects_with_controls %>% filter(projectid == j) %>% 
    pull(siteid) %>% 
    unique() -> sites
  
  for(i in sites) {
    projects_with_controls %>% 
      filter(projectid == j) %>% 
      filter(siteid == i) %>% 
      pull(speciesid) %>% 
      unique() -> species
    
    for(s in species){ 
      filter(dc, projectid == j & siteid == i & speciesid == s) %>% 
        filter(seed == "unseeded") %>% 
        pull(treatmentid) %>% 
        as.character() %>% 
        sample(size = 1) -> control_trt
      
      filter(dc, projectid == j & siteid == i & speciesid == s) %>% 
        filter(seed == "seeded") %>% 
        pull(treatmentid) %>% 
        as.character() %>% 
        sample(size = 1) -> trt
      
      pair <- filter(dc, projectid == j & siteid == i & speciesid == s &
                       treatmentid %in% c(control_trt, trt)) %>% mutate(
                         id = paste(projectid, siteid, speciesid, sep = "_")
                       )
      
      pairs_ps <- rbind(pairs_ps, pair)
    }
  }
}


# Split pairs
gp_unseed <- pairs_ps %>% filter(seed == "unseeded") %>% arrange(projectid, siteid, speciesid, tsr, n)
gp_seeded <- pairs_ps %>% filter(seed == "seeded") %>% arrange(projectid, siteid, speciesid, tsr, n)
sign <- rep(gp_seeded$Ps - gp_unseed$Ps, each = 2)

pairs_ps$diff <- sign
pairs_ps$sign <- "neutral"

pairs_ps[pairs_ps$diff > 0, ]$sign = "positive"
pairs_ps[pairs_ps$diff < 0, ]$sign = "negative"

pairs_ps$seed <- factor(pairs_ps$seed, levels = c("unseeded", "seeded")) 
pairs_ps <-pairs_ps %>% droplevels()

# Data description--------------------------------------------------------------
# Variables description---------------------------------------------------------
# Life form
g1 <-
  ggplot(dc_clean, aes(x = seed)) +
  geom_bar(fill = "steelblue", alpha = .8) +
  #scale_x_discrete(labels = c("annual", "p.forb", "p.grass", "woody")) +
  theme_bw(base_size = 16) +
  labs(y = "Frequency",
       x = "Treatment")
g1

# Ps
g2 <-
  ggplot(dc_clean, aes(x = Ps)) +
  geom_histogram(alpha = 1, fill = "orange") +
  theme_bw(base_size = 16) +
  labs(y = "Frequency",
       x = "Species success")
g2


gall <- g2 + g1  +
  plot_annotation(tag_levels = 'A') &
  theme(plot.title = element_text(size = 18),
        plot.subtitle = element_text(size = 16))
gall

ggsave(gall,
       filename = "outputs/figures/revision/Q2_supp_Fig_variables_description.png",
       width = 9,
       height = 5)

# Saved used data---------------------------------------------------------------
saveRDS(dc_clean, "outputs/data/revision/Q2_used_data.Rds")
saveRDS(pairs_ps, "outputs/data/revision/Q2_pairs_seeded_vs_unseeded.RDs")

# Fit model ------------------------------------------------------------
t.test(gp_seeded$Ps, gp_unseed$Ps, paired = TRUE) # paired t-test
mod <- glmmTMB(Ps ~ seed + 
                 (1 | projectid / siteid ),
               family = binomial(link = "logit"), 
               weights = n,
               ziformula = ~ 1,
               data = dc_clean)
summary(mod)

# Model table------
tab_model(mod, transform = NULL)
tab_model(mod, transform = NULL, show.reflvl = TRUE, 
          file = "outputs/tables/revision/Q2_table_fitted_model.html")

# Model diagnostic------------
res <- simulateResiduals(mod)
plot(res)
table(pairs_ps$sign)
testUniformity(res)
testZeroInflation(res)
testDispersion(res)

# Model residuals against predictors & other variables
plotResiduals(simulationOutput = res,  form = dc_clean$seed)       # Seed x unseed
plotResiduals(simulationOutput = res,  form = log10(dc_clean$tsr)) # time

# Random effects
plot_model(mod, type = "diag")

# Plot predicted effect---------------------------------------------------------
mnc1  <- ggeffect(model = mod, terms = c("seed"))
preds <- as.data.frame(mnc1)

ggplot() +
  geom_jitter(data =  dc %>% select(predicted = Ps, x = seed, n = n),
              aes(y = predicted, x = x, size = n),
              alpha = .45, show.legend = FALSE, width = .07, color = gray(.7)) +
  geom_pointrange(data = preds, size = 1.5, color = gray(.2),
                  aes(x = x, y = predicted, ymax = conf.high, ymin = conf.low)) +
  scale_size_area() +
  scale_y_continuous(breaks = seq(0, 1, .2), 
                     labels = 100 * seq(0, 1, .2)) +
  scale_x_discrete(limits = c("unseeded", "seeded")) +
  
  theme_classic(base_size = 14) +
  theme(axis.title = element_text(size = 16)) +
  theme(text=element_text(family="Garamond")) + 
  labs(x = "Treatment",
       y = "Present (%)")

# Save model object-------------------------------------------------------------
saveRDS(object = mod, file = "outputs/models/revision/Q2_fitted_model.RDs")
saveRDS(object = res, file = "outputs/models/revision/Q2_simulated_residuals.Rds")
