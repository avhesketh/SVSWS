# Littorine and limpet abundance
# January 2023 AH

#######################
# Cleaning grazer data##

# load packages
pkgs <- c("tidyverse","glmmTMB","DHARMa", "plotrix", "lubridate", "emmeans","car")
lapply(pkgs, library, character.only = TRUE)
rm(pkgs)

# names of grazers we want to model
limpets <- c("Lottia_digitalis","Lottia_scutum","Lottia_paradigitalis",
             "Lottia_pelta","Lottia_sp_recruits")
littorines <- c("Littorina_sitkana","Littorina_scutulata")

# clean data for analyses and plotting
grazers_y2 <- surveys %>% 
  left_join(survey_numbers) %>% 
  # retain only species data collected in year 2 of study and just limpets and littorines
  filter(species %in% c(limpets, littorines) & date >= "2020-03-15") %>% 
  mutate(grazer = if_else(species %in% limpets, "Lottia spp.", "Littorina spp."),
      timesincestart = difftime(date, ymd("2019-04-12"), units = c("weeks")),
         count = as.integer(count)) %>% 
  group_by(date, treatment, grazer, timesincestart,block, tile_id, trt_y1, trt_y2, original_herb_trt) %>% 
  # calculate total grazer abundance for limpets and littorines (across species)
  summarise(total_grazers = sum(count)) %>% ungroup() %>% 
  # only retain tiles that continued into year 2 (some were damaged and excluded early in year 2)
  filter((treatment %in% c("C","W")) == F) %>% 
  # flag tiles that had herbivore treatments originally vs. those with herbivores
  mutate(herb = if_else(is.na(original_herb_trt), "control","herb")) %>%
  mutate(day = yday(date), treatment = factor(treatment), block = factor(block),
         tile_id = factor(tile_id)) %>% 
  mutate(treatment = factor(treatment, levels = c("CC", "CW", "WC", "WW")))


#######################
# Plotting

# define colour palette for plotting
pal.trt.y2 <- c("#014779", "#7985CB","#9C0098", "#EE4B2B")

# summarize mean grazer abundance across all tiles in treatment
grazers_plot <- grazers_y2 %>% group_by(date, grazer, treatment) %>% 
  summarize(mean_grazers = mean(total_grazers),
            se_grazers = std.error(total_grazers)) %>% ungroup()

# limpet and littorine abundance through time during year two
FigA7 <- ggplot(grazers_plot, aes(x=date, y=mean_grazers, col = treatment)) +
  geom_point()+
  geom_line(lwd = 0.8) +
  theme_classic() +
  facet_grid(rows = vars(grazer), scales = "free") +
  theme(plot.tag = element_text(face = "bold"),
        strip.text = element_text(face = "italic")) +
  labs(x = "Date", y = "Grazer abundance", col = "Treatment") +
  scale_color_manual(values = pal.trt.y2) +
  geom_errorbar(aes(ymax = mean_grazers + se_grazers, 
                    ymin = mean_grazers - se_grazers), width = 5) +
  scale_x_date(date_labels = "%Y-%m")

png("./figures/FigA7.png", res = 700, width = 8.5, height = 5, units = "in")
FigA7
dev.off()

# boxplot of abundances at key timepoints (post-summer, post-winter) only

grazers_keytimes <- grazers_y2 %>% 
  filter(date %in% c("2020-09-14", "2021-02-24")) %>% 
  mutate(date = if_else(date == "2020-09-14", "Post-summer", "Winter"),
         grazer = factor(grazer, levels = c("Lottia spp.", "Littorina spp.")))

Fig3 <- ggplot(grazers_keytimes, aes(x = treatment, y = total_grazers, col = treatment)) +
  geom_boxplot(size = 0.4, outlier.color = NA) +
  geom_jitter(size = 0.7) +
  labs(x = "Treatment", y = "Total abundance", col = "Treatment") +
  facet_grid(cols = vars(date), rows = vars(grazer), scales = "free_y") +
  scale_color_manual(values = pal.trt.y2) +
  theme_classic() +
  theme(strip.text.y = element_text(face = "italic"))

png("./figures/Fig3.png", res = 700, width = 7, height = 5, units = "in")
Fig3
dev.off()

#######################
# Modeling limpet & littorine abundance after summer heat and after winter recovery


# Model for limpet abundance after summer heat stress
limpet.abund.t1 <- glmmTMB(total_grazers ~ trt_y1*trt_y2 +
                               (1|block), data = grazers_y2 %>% filter(date == "2020-09-14" &
                                                                         grazer == "Lottia spp."),
                             family = nbinom1())
limpet.abund.herb.t1 <- glmmTMB(total_grazers ~ trt_y1*trt_y2 + herb +
                               (1|block), data = grazers_y2 %>% filter(date == "2020-09-14" &
                                                                        grazer == "Lottia spp."),
                             family = nbinom1())

# AIC(limpet.abund.t1, limpet.abund.herb.t1) # herbivore term does not improve fit

plot(simulateResiduals(limpet.abund.t1))
summary(limpet.abund.t1)
Anova(limpet.abund.t1, type = 3, contrasts=list(topic=contr.sum, sys=contr.sum))
Anova(limpet.abund.t1, type = 2)

# Model for limpet abundance after winter recovery
limpet.abund.t2 <- glmmTMB(total_grazers ~ trt_y1*trt_y2 +
                             (1|block), data = grazers_y2 %>% filter(date == "2021-02-24" &
                                                                       grazer == "Lottia spp."),
                           family = nbinom1())
limpet.abund.herb.t2 <- glmmTMB(total_grazers ~ trt_y1*trt_y2 + herb +
                                  (1|block),data = grazers_y2 %>% filter(date == "2021-02-24" &
                                                                           grazer == "Lottia spp."),
                                family = nbinom1())

# AIC(limpet.abund.t2, limpet.abund.herb.t2) # does not improve fit

plot(simulateResiduals(limpet.abund.t2))
summary(limpet.abund.t2)
Anova(limpet.abund, type = 3, contrasts=list(topic=contr.sum, sys=contr.sum))
Anova(limpet.abund, type = 2)



# Model for littorine abundance after summer heat stress
litt.abund.t1 <- glmmTMB(total_grazers ~ trt_y1*trt_y2 
                      + (1|block),
                      data = grazers_y2 %>% filter(date == "2020-09-14" &
                                                   grazer == "Littorina spp."),
                      family = nbinom1())
litt.abund.herb.t1 <- glmmTMB(total_grazers ~ trt_y1*trt_y2 + herb
                           + (1|block), 
                           data = grazers_y2 %>% filter(date == "2020-09-14" &
                                                          grazer == "Littorina spp."),
                      family = nbinom1())
AIC(litt.abund.t1, litt.abund.herb.t1) # no difference

plot(simulateResiduals(litt.abund.t1))
summary(litt.abund.t1)
Anova(litt.abund.t1, type = 3, contrasts=list(topic=contr.sum, sys=contr.sum))
Anova(litt.abund.t1, type = 2)


# Model for littorine abundance after winter recovery
litt.abund.t2 <- glmmTMB(total_grazers ~ trt_y1*trt_y2 
                         + (1|block),
                         data = grazers_y2 %>% filter(date == "2021-02-24" &
                                                        grazer == "Littorina spp."),
                         family = nbinom1())
litt.abund.herb.t2 <- glmmTMB(total_grazers ~ trt_y1*trt_y2 + herb
                              + (1|block), 
                              data = grazers_y2 %>% filter(date == "2021-02-24" &
                                                             grazer == "Littorina spp."),
                              family = nbinom1())
AIC(litt.abund.t2, litt.abund.herb.t2) # no difference

plot(simulateResiduals(litt.abund.t2))
summary(litt.abund.t2)
Anova(litt.abund.t2, type = 3, contrasts=list(topic=contr.sum, sys=contr.sum))
Anova(litt.abund.t2, type = 2)

