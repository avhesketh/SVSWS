# Algal cover shifts & composition changes through time
# February 2024

##################
# Cleaning data for analyses & plotting

# load packages
pkgs <- c("tidyverse","glmmTMB","DHARMa", "plotrix", "lubridate", "emmeans","car",
          "vegan")
lapply(pkgs, library, character.only = TRUE)
rm(pkgs)

# read in algal cover data
algae <- read_csv("./clean_data/SVSWS_survey_clean.csv") %>% 
  filter(species %in% c("Ulothrix_sp", "Ulva_sp", "Savoiea_robusta", "Pyropia_sp",
                       "Leathesia_marina", "Fucus_distichus","Endocladia_muricata",
                       "Petalonia_fascia", "Scunge","Mastocarpus_sp_crust")) %>% 
  left_join(survey_numbers) %>% select(-count)

algae_y1 <- algae %>% 
  filter(date <= as.Date("2020-03-15")) %>% 
  mutate(treatment = trt_y1) %>% 
  group_by(date, block, tile_id, trt_y1, trt_y2, treatment) %>% 
  summarize(algal_cover = sum(percent_cover)) %>% 
  mutate(
         day = yday(date),
         block = factor(block),
         treatment = factor(treatment, levels = c("C","W")),
         tile_id = factor(tile_id),
         period = "Year 1") %>% ungroup()

algae_y2 <- algae %>% 
  filter(date > "2020-03-15" & (treatment %in% c("W","C"))==F) %>% 
  group_by(date, block, tile_id, trt_y1, trt_y2, treatment) %>% 
  summarize(algal_cover = sum(percent_cover)) %>%
  mutate(
         block = factor(block),
         trt_y1 = factor(trt_y1, levels = c("C","W")),
         trt_y2 = factor(trt_y2, levels = c("C","W")),
         tile_id = factor(tile_id),
         day = yday(date),
         treatment = factor(treatment),
         period = "Year 2") %>% ungroup()

##############
# Visualizing data

pal.trt <- c("#014779", "#EE4B2B", "#014779","#7985CB", "#9C0098", "#EE4B2B")
pch.trt <- c(1,1,16,16,16,16)
line.trt <- c("dotdash","dotdash","solid","solid","solid","solid")

# plot algal cover over time
algae_summ <- algae_all %>% 
  group_by(date, treatment) %>% 
  summarize(mean_cover = mean(algal_cover), se_cover = std.error(algal_cover)) %>% 
  mutate(treatment = factor(treatment, levels = c("C","W","CC","CW","WC","WW"))) %>% 
  ungroup()

FigA8 <- ggplot(aes(x = date, y = mean_cover, col = treatment, pch = treatment, lty=treatment), data = algae_summ) +
  plot_theme +
  geom_point() +
  geom_line(lwd = 0.8 ,aes(group = treatment)) +
  geom_errorbar(aes(ymax = mean_cover + se_cover, ymin = mean_cover - se_cover, lty = NULL)) +
  labs(x = "Date", y ="Algal cover (%)", col = "Treatment",
       lty = "Treatment", pch = "Treatment") +
  scale_color_manual(values = pal.trt)+
  scale_linetype_manual(values = line.trt) +
  scale_shape_manual(values = pch.trt) +
  theme(legend.key.width = unit(1, "cm"),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10)) +
  geom_vline(aes(xintercept = ymd("2020-04-07")), linetype = "dotted", col = "grey30", lwd = 0.8)

png("./figures/FigA8.png", res = 700, width = 8, height = 4, units = "in")
FigA8
dev.off()

## plotting algal cover at key timepoints

algae_keytimes <- algae_y1 %>% 
  full_join(algae_y2) %>% 
  filter(date %in% c("2019-10-20", "2020-03-15","2020-09-14", "2021-02-24")) %>% 
  mutate(timept = if_else(date %in% c("2019-10-20", "2020-09-14"), "Post-summer",
                          "Winter"),
         treatment = factor(treatment, levels = c("C","W","CC","CW","WC","WW")))

Fig4 <- ggplot(aes(x = timept, y = algal_cover, col = treatment, pch = treatment), 
       data = algae_keytimes) +
  geom_boxplot(outlier.color = NA) +
  geom_point(position = position_jitterdodge(), size = 0.8)+
  facet_wrap(~period, scales ="free_x") +
  theme_classic() +
  scale_color_manual(values = pal.trt) +
  scale_shape_manual(values = pch.trt) +
  labs(x = "Time", pch = "Treatment", col = "Treatment", y = "Algal cover (%)") 
#+ theme(legend.position = "none",plot.tag = element_text(face = "bold"))

png("./figures/Fig4.png", res = 700, width = 6, height = 4, units = "in")
Fig4
dev.off()

##############
# modeling data at each timepoint
algae_models <- algae_keytimes %>% 
  mutate(algal_cover = if_else(algal_cover > 100, 100, algal_cover)) %>% 
  mutate(algae.tr = ((algal_cover*(nrow(.)-1)+0.5)/nrow(.))/100)

## model year 1: post-summer

algae.model.stress.y1 <- glmmTMB(algae.tr ~ trt_y1 + (1|block),
                            family = beta_family(),
                            data = algae_models %>% filter(date =="2019-10-20"))
res.ams1 <- simulateResiduals(algae.model.stress.y1)
plot(res.ams1)
summary(algae.model.stress.y1)
Anova(algae.model.stress.y1, type = 2)

## model year 1: winter recovery

algae.model.recovery.y1 <- glmmTMB(algal_cover ~ trt_y1 + (1|block),
                                 family = tweedie(),
                                 data = algae_models %>% filter(date =="2020-03-15"))
res.amr1 <- simulateResiduals(algae.model.recovery.y1)
plot(res.amr1)
summary(algae.model.recovery.y1)
Anova(algae.model.recovery.y1, type = 2)

## model year 2: can't model post-summer since there are too few non-zero data points
## model year 2: winter recovery


algae.model.recovery.y2 <- glmmTMB(algae.tr ~ trt_y1*trt_y2 + (1|block),
                          family = beta_family(), data = algae_models %>% filter(date=="2021-02-24"))

res.amr2 <- simulateResiduals(algae.model.recovery.y2)
plot(res.amr2)
summary(algae.model.recovery.y2)
Anova(algae.model.recovery.y2, type = 3)

comparison.y2 <- emmeans(algae.model.recovery.y2, specs = pairwise ~ trt_y1:trt_y2)
comparison.y2$emmeans
comparison.y2$contrasts
