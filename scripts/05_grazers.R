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
  group_by(date, treatment, grazer, timesincestart,block, tile_id, trt_y1, trt_y2, second_herb_trt) %>% 
  # calculate total grazer abundance for limpets and littorines (across species)
  summarise(total_grazers = sum(count)) %>% ungroup() %>% 
  # only retain tiles that continued into year 2 (some were damaged and excluded early in year 2)
  filter((treatment %in% c("C","W")) == F) %>% 
  mutate(herb_trt = if_else(is.na(second_herb_trt), "control", "grazer"), 
         day = yday(date), treatment = factor(treatment), block = factor(block),
         tile_id = factor(tile_id)) %>% 
  mutate(treatment = factor(treatment, levels = c("CC", "CW", "WC", "WW")))

#######################
# Modeling limpet & littorine abundance after summer heat and after winter recovery

grazers_y2 <- grazers_y2 %>% mutate(date = factor(date))

# Model for limpet abundance post-summer

limpet.abund.ps.herb <- glmmTMB(total_grazers ~ trt_y1*trt_y2 + herb_trt + 
                        (1|block), data = grazers_y2 %>% filter(date == "2020-09-14" &
                                                                  grazer == "Lottia spp."),
                        family = poisson())
plot(simulateResiduals(limpet.abund.ps.herb))
summary(limpet.abund.ps.herb)
Anova(limpet.abund.ps.herb, type = 3, contrasts=list(trt_y1 = "contr.sum", trt_y2 = "contr.sum")) 
# effect of grazers in y1 not significant; drop this term


limpet.abund.ps <- glmmTMB(total_grazers ~ trt_y1*trt_y2 +
                             (1|block), data = grazers_y2 %>% filter(date %in% c("2020-09-14") &
                                                                       grazer == "Lottia spp."),
                           family = poisson())
plot(simulateResiduals(limpet.abund.ps))
summary(limpet.abund.ps)
Anova(limpet.abund.ps, type = 3, contrasts=list(trt_y1 = "contr.sum", trt_y2 = "contr.sum"))


emm.lott.ps <- emmeans(limpet.abund.ps, specs = c("trt_y1","trt_y2"))
contrast(emm.lott.ps, method = "pairwise", adjust = "tukey")


# limpet abundance in winter

limpet.abund.w <- glmmTMB(total_grazers ~ trt_y1*trt_y2 +
                          (1|block), data = grazers_y2 %>% filter(date == "2021-02-24" &
                                                                    grazer == "Lottia spp."),
                        family = nbinom1())

plot(simulateResiduals(limpet.abund.w))
summary(limpet.abund.w)
Anova(limpet.abund.w, type = 3, contrasts=list(trt_y1 = "contr.sum", trt_y2 = "contr.sum"))

emm.lott.w <- emmeans(limpet.abund.w, specs = c("trt_y1","trt_y2"))
contrast(emm.lott.w, method = "pairwise", adjust = "tukey")


# Littorine abundance post-summer
litt.abund.ps.herb <- glmmTMB(total_grazers ~ trt_y1*trt_y2 + herb_trt
                              + (1|block), 
                              data = grazers_y2 %>% filter(date == "2020-09-14" &
                                                             grazer == "Littorina spp."),
                              family = nbinom1())
plot(simulateResiduals(litt.abund.ps.herb))
summary(litt.abund.ps.herb)
Anova(litt.abund.ps.herb, type = 3, contrasts=list(trt_y1 = "contr.sum", trt_y2 = "contr.sum"))
# grazer manipulation in y1 not significant; drop this term

litt.abund.ps <- glmmTMB(total_grazers ~ trt_y1*trt_y2 
                      + (1|block),
                      data = grazers_y2 %>% filter(date == "2020-09-14" &
                                                   grazer == "Littorina spp."),
                      family = nbinom1())
plot(simulateResiduals(litt.abund.ps))
summary(litt.abund.ps)
Anova(litt.abund.ps, type = 3, contrasts=list(trt_y1 = "contr.sum", trt_y2 = "contr.sum"))

emm.litt.ps <- emmeans(litt.abund.ps, specs = c("trt_y1","trt_y2"))
contrast(emm.litt.ps, method = "pairwise", adjust = "tukey")

# Littorine abundance in winter
litt.abund.w <- glmmTMB(total_grazers ~ trt_y1*trt_y2 
                      + (1|block),
                      data = grazers_y2 %>% filter(date == "2021-02-24" &
                                                     grazer == "Littorina spp."),
                      family = nbinom1())
plot(simulateResiduals(litt.abund.w))
summary(litt.abund.w)
Anova(litt.abund.w, type = 3, contrasts=list(trt_y1 = "contr.sum", trt_y2 = "contr.sum"))

emm.litt.w <- emmeans(litt.abund.w, specs = c("trt_y1","trt_y2"))
contrast(emm.litt.w, method = "pairwise", adjust = "tukey")


# create dataframe with significance labels from emmeans comparisons
labels.4a <- as.data.frame(cbind(
  c(rep("Post-summer", times = 4), rep("Winter", times = 4)),
  c(rep(c("CC","CW","WC","WW"),times = 2)),
  c(20,20,20,20,34,20,20,20),
  c("a","bc","ab","c","d","de","d","e"))
)
colnames(labels.4a) <- c("season","treatment","total_grazers","label")
labels.4a <- labels.4a %>% mutate(total_grazers = as.numeric(total_grazers),
                                season = factor(season, levels = c("Post-summer", "Winter")))

labels.4b <- as.data.frame(cbind(
  c(rep("Post-summer", times = 4), rep("Winter", times = 4)),
  c(rep(c("CC","CW","WC","WW"),times = 2)),
  c(215,150,150,150, 150, 150,150,150),
  c("a","b","bc","c","d","de","de","e"))
)
colnames(labels.4b) <- c("season","treatment","total_grazers","label")
labels.4b <- labels.4b %>% mutate(total_grazers = as.numeric(total_grazers),
                                  season = factor(season, levels = c("Post-summer", "Winter")))

#######################
# Plotting

# define colour palette for plotting
pal.trt.y2 <- c("#014779", "#7985CB","#9C0098", "#EE4B2B")

# boxplot of abundances at key timepoints (post-summer, post-winter) only

grazers_keytimes <- grazers_y2 %>% 
  filter(date %in% c("2020-09-14", "2021-02-24")) %>% 
  mutate(season = if_else(date == "2020-09-14", "Post-summer", "Winter"),
         grazer = factor(grazer, levels = c("Lottia spp.", "Littorina spp.")))

Fig4a <- ggplot(grazers_keytimes %>% filter(grazer == "Lottia spp."), 
                aes(x = treatment, y = total_grazers,
                                     col = treatment,
                                     alpha = season)) +
  geom_boxplot(size = 0.4, outlier.color = NA) +
  labs(x = "Treatment", y = expression("Abundance of"~italic("Lottia")~"spp."), col = "Treatment", pch = "Season") +
  facet_grid(cols = vars(season), scales = "free_y") +
  scale_color_manual(values = pal.trt.y2, guide = "none") +
  scale_fill_manual(values = pal.trt.y2, guide = "none") +
  scale_alpha_manual(values = c(0.4,0)) +
  scale_shape_manual(values = c(16,17)) +
  guides(fill = "none", alpha = "none", shape = guide_legend(override.aes = list(size = 3))) +
  geom_jitter(size = 1, alpha = 0.7,
              position = position_jitterdodge(jitter.width = 1, 
                                              jitter.height = 0, 
                                              dodge.width = 1), aes(shape = season)) +
  theme_classic() +
  theme(plot.tag = element_text(face = "bold", size = 10),
        strip.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.box.margin=margin(-5,-5,-5,-5)) +
  geom_text(data = labels.4a, 
            aes(x = treatment, y = total_grazers, label = label), 
            col = "black", fontface = "bold", size = 3, inherit.aes = FALSE) +
  lims(y = c(0,35)) +
  guides(pch = guide_legend(override.aes = list(size = 2) ) )
Fig4a

Fig4b <- ggplot(grazers_keytimes %>% filter(grazer == "Littorina spp."), 
                aes(x = treatment, y = total_grazers,
                    col = treatment,
                    alpha = season)) +
  geom_boxplot(size = 0.4, outlier.color = NA) +
  labs(x = "Treatment", y = expression("Abundance of"~italic("Littorina")~"spp."), col = "Treatment", pch = "Season") +
  facet_grid(cols = vars(season), scales = "free_y") +
  scale_color_manual(values = pal.trt.y2, guide = "none") +
  scale_fill_manual(values = pal.trt.y2, guide = "none") +
  scale_alpha_manual(values = c(0.4,0)) +
  scale_shape_manual(values = c(16,17)) +
  guides(fill = "none", alpha = "none", shape = guide_legend(override.aes = list(size = 3))) +
  geom_jitter(size = 1, alpha = 0.7,
              position = position_jitterdodge(jitter.width = 1, 
                                              jitter.height = 0, 
                                              dodge.width = 1), aes(shape = season)) +
  theme_classic() +
  theme(plot.tag = element_text(face = "bold", size = 10),
        strip.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.box.margin=margin(-5,-5,-5,-5)) +
  lims(y = c(0,230)) +
  geom_text(data = labels.4b, 
            aes(x = treatment, y = total_grazers, label = label), 
            col = "black", fontface = "bold", size = 3, inherit.aes = FALSE) +
  guides(pch = guide_legend(override.aes = list(size = 2) ) )
Fig4b

Fig4 <- (Fig4a / Fig4b) + plot_annotation(tag_levels = "a", tag_prefix = "(",
                                          tag_suffix = ")") + plot_layout(guides = "collect") & theme(legend.position = "bottom")
Fig4

ggsave(Fig4, filename = "./figures/Fig4.pdf", device = cairo_pdf, 
       width = 8.5, height = 15, units = "cm")


#########################
# Appendix figure for grazer abundance through time #

# summarize mean grazer abundance across all tiles in treatment
grazers_plot <- grazers_y2 %>% group_by(date, grazer, treatment) %>% 
  summarize(mean_grazers = mean(total_grazers),
            se_grazers = std.error(total_grazers)) %>% ungroup() %>% 
  mutate(date = ymd(date))

# limpet and littorine abundance through time during year two
FigS9 <- ggplot(grazers_plot, aes(x=date, y=mean_grazers, col = treatment)) +
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

png("./figures/FigS9.png", res = 700, width = 8, height = 4, units = "in",
    pointsize = 22)
FigS9
dev.off()
