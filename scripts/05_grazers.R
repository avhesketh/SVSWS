# Littorine and limpet abundance
# January 2023 AH

#######################
# Cleaning grazer data #

# load packages
pkgs <- c("tidyverse","glmmTMB","DHARMa", "plotrix", "lubridate", "emmeans","car")
lapply(pkgs, library, character.only = TRUE)
rm(pkgs)

# load in tidy survey data
surveys <- read_csv("./clean_data/SVSWS_survey_clean.csv")


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
  mutate(original_herb_trt = if_else(is.na(original_herb_trt), "control", original_herb_trt), 
         day = yday(date), treatment = factor(treatment), block = factor(block),
         tile_id = factor(tile_id)) %>% 
  mutate(treatment = factor(treatment, levels = c("CC", "CW", "WC", "WW")))

#######################
# Modeling limpet & littorine abundance after summer heat and after winter recovery

grazers_y2 <- grazers_y2 %>% mutate(date = factor(date))

# Model for limpet abundance
limpet.abund <- glmmTMB(total_grazers ~ trt_y1*trt_y2 + date +
                               (1|block), data = grazers_y2 %>% filter(date %in% c("2020-09-14","2021-02-24") &
                                                                         grazer == "Lottia spp."),
                             family = nbinom1())
limpet.abund.herb <- glmmTMB(total_grazers ~ trt_y1*trt_y2 + date + original_herb_trt +
                               (1|block), data = grazers_y2 %>% filter(date %in% c("2020-09-14","2021-02-24") &
                                                                        grazer == "Lottia spp."),
                             family = nbinom1())

AIC(limpet.abund, limpet.abund.herb) # herbivore term does not improve fit

plot(simulateResiduals(limpet.abund))
summary(limpet.abund)
Anova(limpet.abund, type = 3)


emm.lott <- emmeans(limpet.abund, specs = c("trt_y1","trt_y2"))
contrast(emm.lott, method = "pairwise", adjust = "tukey")


# Model for littorine abundance
litt.abund <- glmmTMB(total_grazers ~ trt_y1*trt_y2 + date  
                      + (1|block),
                      data = grazers_y2 %>% filter(date %in% c("2020-09-14","2021-02-24") &
                                                   grazer == "Littorina spp."),
                      family = nbinom1())
litt.abund.herb <- glmmTMB(total_grazers ~ trt_y1*trt_y2 + date + original_herb_trt
                           + (1|block), 
                           data = grazers_y2 %>% filter(date %in% c("2020-09-14","2021-02-24") &
                                                          grazer == "Littorina spp."),
                      family = nbinom1())
AIC(litt.abund, litt.abund.herb) # herbivore term does not improve fit

plot(simulateResiduals(litt.abund))
summary(litt.abund)
Anova(litt.abund, type = 3)

emm.litt <- emmeans(litt.abund, specs = c("trt_y1","trt_y2"))
contrast(emm.litt, method = "pairwise", adjust = "tukey")


# create dataframe with significance labels from emmeans comparisons
labels.3 <- as.data.frame(cbind(
  c(rep(c("CC","CW","WC","WW"),times = 2)),
  c(rep("Lottia spp.",times = 4), rep("Littorina spp.", times = 4)),
  c(33,20,20,20, 205,148,148,148),
  c("a","bc","ac","b","d","e","e","e"))
)
colnames(labels.3) <- c("treatment","grazer","total_grazers","label")
labels.3 <- labels.3 %>% mutate(total_grazers = as.numeric(total_grazers),
                                grazer = factor(grazer, levels = c("Lottia spp.", "Littorina spp.")))

#######################
# Plotting

# define colour palette for plotting
pal.trt.y2 <- c("#014779", "#7985CB","#9C0098", "#EE4B2B")

# boxplot of abundances at key timepoints (post-summer, post-winter) only

grazers_keytimes <- grazers_y2 %>% 
  filter(date %in% c("2020-09-14", "2021-02-24")) %>% 
  mutate(date = if_else(date == "2020-09-14", "Post-summer", "Winter"),
         grazer = factor(grazer, levels = c("Lottia spp.", "Littorina spp.")))

Fig3 <- ggplot(grazers_keytimes, aes(x = treatment, y = total_grazers,
                                     col = treatment,
                                     fill = treatment,
                                     alpha = date)) +
  geom_boxplot(size = 0.4, outlier.color = NA) +
  labs(x = "Treatment", y = "Total abundance", col = "Treatment", pch = "Time of year") +
  facet_grid(rows = vars(grazer), scales = "free_y") +
  scale_color_manual(values = pal.trt.y2, guide = "none") +
  scale_fill_manual(values = pal.trt.y2, guide = "none") +
  scale_alpha_manual(values = c(0.4,0)) +
  scale_shape_manual(values = c(16,17)) +
  guides(fill = "none", alpha = "none", shape = guide_legend(override.aes = list(size = 3))) +
  geom_jitter(size = 0.7, alpha = 1,
              position = position_jitterdodge(jitter.width = 0.5, 
                                              jitter.height = 0, 
                                              dodge.width = 0.7), aes(shape = date)) +
  theme_classic() +
  theme(plot.tag = element_text(face = "bold", size = 10),
        strip.text = element_text(size = 10, face = "italic"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.position = "top",
        legend.box.margin=margin(-5,-5,-5,-5)) +
  geom_text(data = labels.3, 
            aes(x = treatment, y = total_grazers, label = label), 
            col = "black", fontface = "bold", size = 2.5, inherit.aes = FALSE) +
  guides(pch = guide_legend(override.aes = list(size = 2) ) )


ggsave(Fig3, filename = "./figures/Fig3.pdf", device = cairo_pdf, 
       width = 8.5, height = 10, units = "cm")


#########################
# Appendix figure for grazer abundance through time #

# summarize mean grazer abundance across all tiles in treatment
grazers_plot <- grazers_y2 %>% group_by(date, grazer, treatment) %>% 
  summarize(mean_grazers = mean(total_grazers),
            se_grazers = std.error(total_grazers)) %>% ungroup() %>% 
  mutate(date = ymd(date))

# limpet and littorine abundance through time during year two
FigS7 <- ggplot(grazers_plot, aes(x=date, y=mean_grazers, col = treatment)) +
  geom_point()+
  geom_line(lwd = 0.8) +
  theme_classic() +
  facet_grid(rows = vars(grazer), scales = "free") +
  theme(plot.tag = element_text(face = "bold"),
        strip.text = element_text(face = "italic", size = 12)) +
  labs(x = "Date", y = "Grazer abundance", col = "Treatment") +
  scale_color_manual(values = pal.trt.y2) +
  geom_errorbar(aes(ymax = mean_grazers + se_grazers, 
                    ymin = mean_grazers - se_grazers), width = 5) +
  scale_x_date(date_labels = "%Y-%m", breaks=c(ymd("2020-04-01"), 
                                               ymd("2020-06-01"),
                                               ymd("2020-08-01"),
                                               ymd("2020-10-01"),
                                               ymd("2020-12-01"),
                                               ymd("2021-02-01")))

png("./figures/FigS7.png", res = 700, width = 8, height = 4, units = "in",
    pointsize = 22)
FigS7
dev.off()
