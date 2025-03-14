# Algal cover shifts & composition changes through time
# February 2024

##################
# Cleaning data for analyses & plotting

# load packages
pkgs <- c("tidyverse","glmmTMB","DHARMa", "plotrix", "lubridate", "emmeans","car", "mgcv")
lapply(pkgs, library, character.only = TRUE)
rm(pkgs)

# read in algal cover data
algae <- read_csv("./clean_data/SVSWS_survey_clean.csv") %>% 
  filter(species %in% c("Ulothrix_sp", "Ulva_sp", "Savoiea_robusta", "Pyropia_sp",
                       "Leathesia_marina", "Fucus_distichus","Endocladia_muricata",
                       "Petalonia_fascia", "Scunge","Mastocarpus_sp_crust")) %>% 
  left_join(survey_numbers) %>% select(-count)

# Retain Year 1 data only
algae_y1 <- algae %>% 
  filter(date <= as.Date("2020-03-15")) %>% 
  mutate(treatment = trt_y1) %>% 
  group_by(date, block, tile_id, trt_y1, trt_y2, treatment, second_herb_trt) %>% 
  summarize(algal_cover = sum(percent_cover)) %>% 
  mutate(
         day = yday(date),
         block = factor(block),
         treatment = factor(treatment, levels = c("C","W")),
         tile_id = factor(tile_id),
         period = "Year 1",
         herb_trt = if_else(is.na(second_herb_trt), "control", "grazer")) %>% ungroup()

algae_y1_reduced <- algae_y1 %>% 
  filter(date > as.Date("2019-08-27"))

# Retain Year 2 data only
algae_y2 <- algae %>% 
  filter(date >= "2020-03-15" & (treatment %in% c("W","C"))==F) %>% 
  group_by(date, block, tile_id, trt_y1, trt_y2, treatment ,second_herb_trt) %>% 
  summarize(algal_cover = sum(percent_cover)) %>%
  mutate(
         block = factor(block),
         trt_y1 = factor(trt_y1, levels = c("C","W")),
         trt_y2 = factor(trt_y2, levels = c("C","W")),
         tile_id = factor(tile_id),
         day = yday(date),
         treatment = factor(treatment),
         period = "Year 2",
         herb_trt = if_else(is.na(second_herb_trt), "control", "grazer")) %>% ungroup()

algae_all <- algae_y1 %>% full_join(algae_y2) %>% 
  mutate(original_herb_trt = if_else(is.na(original_herb_trt), "control",original_herb_trt))

##############
# Visualizing algal cover over time

# Define palettes
pal.trt <- c("#014779", "#EE4B2B", "#014779","#7985CB", "#9C0098", "#EE4B2B")
pch.trt <- c(1,1,16,16,16,16)
line.trt <- c("dotdash","dotdash","solid","solid","solid","solid")

# Create summary dataframe with one point per date per treatment
algae_summ <- algae_all %>% 
  group_by(date, treatment) %>% 
  summarize(mean_cover = mean(algal_cover), se_cover = std.error(algal_cover)) %>% 
  mutate(treatment = factor(treatment, levels = c("C","W","CC","CW","WC","WW"))) %>% 
  ungroup()

# Define axis breaks where labels will apear
breaks <- c(ymd("2019-05-01"), ymd("2019-09-01"),ymd("2020-01-01"),
            ymd("2020-05-01"),ymd("2020-09-01"),ymd("2021-01-01"))

# Panel 4a: Algal cover in each treatment over the course of the experiment
Fig5a <- ggplot(aes(x = date, y = mean_cover, col = treatment, pch = treatment, lty=treatment), data = algae_summ) +
  theme_classic()+
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
  scale_x_date(date_labels = "%Y-%m", breaks=breaks, 
               limits = c(ymd("2019-05-03"),ymd("2021-02-24"))) +
  geom_vline(aes(xintercept = ymd("2020-04-03")), col = "grey60", lwd = 0.8) +
  geom_vline(aes(xintercept = ymd("2019-08-27")), col = "grey70", lty = "dashed", lwd = 0.8) +
  theme(plot.tag = element_text(face = "bold", size = 10),
        strip.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))
Fig5a

#########################
# GAM of algal cover over time in each treatment #

# generate dataframe from which to develop models
algae_gam_y1 <- algae_y1 %>% mutate(days.since.start = as.numeric(difftime(date, ymd("2019-04-12"), units = "days")),
                                    oTreatment = ordered(treatment, levels = c("C","W")))

# reduced dataframe = remove data from period during which herbivore community was manipulated
algae_gam_y1_red <- algae_y1_reduced %>% mutate(days.since.start = as.numeric(difftime(date, ymd("2019-04-12"), units = "days")),
                                    oTreatment = ordered(treatment, levels = c("C","W")))

# model of algal cover in year 1: include smooth terms for time, time x treatment, and a random effect
# of block. Include AR1 correlation structure to account for autocorrelation through time.

# herbivore treatment is significant: do not use this model, only model of reduced data
algae.y1.gam.herb <- gamm(algal_cover ~ treatment + herb_trt + s(days.since.start) +
                            s(days.since.start, by = oTreatment) + s(block, bs = "re"),
                          correlation = corAR1(),
                          data = algae_gam_y1, method = "REML")

# this is the final y1 model of algal cover
algae.y1.gam.red <- gamm(algal_cover ~ treatment + s(days.since.start, k =3) +
                       s(days.since.start, by = oTreatment, k =3) + s(block, bs = "re"),
                     correlation = corAR1(),
                     data = algae_gam_y1_red, method = "REML")

summary(algae.y1.gam.herb$gam)
plot(algae.y1.gam$gam, pages = 1, scale = 0 , seWithMean = TRUE)

# look at autocorrelation
acf(residuals(algae.y1.gam.red$gam)) # raw residuals
acf(residuals(algae.y1.gam.red$lme, type = "normalized")) # standardized residuals with AR1 process

summary(algae.y1.gam.red$gam) 
summary(algae.y1.gam.red$lme) 
plot(algae.y1.gam.red$gam)

# Create dataframe for generating model predictions
pdat.red <- expand.grid(days.since.start = seq(137, 338, by = 1),
                    treatment = c("C","W"),
                    block = c("A","B","C","D","E","F")) %>% mutate(oTreatment = treatment)

# Find differences between smoothers from reduced y1 model
comp1.red <- smooth_diff(algae.y1.gam.red$gam, pdat, 'W', 'C', 'treatment')
# Add columns for the days since experiment start and actual date for plotting purposes
comp.red <- cbind(days.since.start = seq(137, 338, by = 1),
              rbind(comp1.red)) %>% 
  mutate(date = ymd("2019-04-12") + days.since.start)

# full gamm for putting on plot

# Create dataframe for generating model predictions
pdat <- expand.grid(days.since.start = seq(27, 338, by = 1),
                        treatment = c("C","W"),
                        block = c("A","B","C","D","E","F")) %>% mutate(oTreatment = treatment)

# Find differences between smoothers from model
comp1 <- smooth_diff(algae.y1.gam$gam, pdat, 'W', 'C', 'treatment')
# Add columns for the days since experiment start and actual date for plotting purposes
comp <- cbind(days.since.start = seq(27, 338, by = 1),
                  rbind(comp1)) %>% 
  mutate(date = ymd("2019-04-12") + days.since.start,
         pair = "W-C")

# GAM for algal cover in year 2
algae_gam_y2 <- algae_y2 %>% mutate(days.since.start = as.numeric(difftime(date, ymd("2019-04-12"), units = "days")),
                                    oTreatment = ordered(treatment, levels = c("CC","CW","WC","WW")),
                                    treatment = factor(treatment, levels = c("CC","CW","WC","WW")))


algae.y2.gam.herb <- gamm(algal_cover ~ trt_y1*trt_y2 + herb_trt +
                       s(days.since.start,k=5) + 
                       s(days.since.start, by = oTreatment, k=5) + 
                       s(block, bs = "re"),
                       correlation = corAR1(),
                     data = algae_gam_y2, method = "REML")
summary(algae.y2.gam.herb$gam) # herbivore term is not significant; drop this

algae.y2.gam <- gamm(algal_cover ~ trt_y1*trt_y2 + s(days.since.start,k=5) + 
                       s(days.since.start, by = oTreatment, k=5) + 
                       s(block, bs = "re"),
                     correlation = corAR1(),
                     data = algae_gam_y2, method = "REML")

# Autocorrelation is corrected
acf(residuals(algae.y2.gam$gam)) # raw residuals
acf(residuals(algae.y2.gam$lme, type = "normalized")) # standardized residuals with AR1 process

summary(algae.y2.gam$gam)
plot(algae.y2.gam$gam)

# Dataframe for generating model predictions for Year 2 GAM
pdat2 <- expand.grid(days.since.start = seq(338, 684, by = 1),
                    treatment = c("CC","CW", "WC","WW"),
                    block = c("A","B","C","D","E","F")) %>% mutate(oTreatment = treatment)

# Generate three sets of pairwise comparisons between control (CC) and all others
y2.comp1 <- cbind(days.since.start = seq(338, 684, by = 1), 
                  smooth_diff(algae.y2.gam$gam, pdat2, 'WC', 'CC', 'treatment'))

y2.comp2 <- cbind(days.since.start = seq(338, 684, by = 1), 
                  smooth_diff(algae.y2.gam$gam, pdat2, 'CW', 'CC', 'treatment'))

y2.comp3 <- cbind(days.since.start = seq(338, 684, by = 1), 
                  smooth_diff(algae.y2.gam$gam, pdat2, 'WW', 'CC', 'treatment'))

# Join comparisons between smoothers into single dataframe
y2.comp <- y2.comp1 %>% full_join(y2.comp2) %>% 
                full_join(y2.comp3) %>% 
                mutate(date = ymd("2019-04-12") + days.since.start)
# Join with Year 1 also
all.comp <- comp %>% full_join(y2.comp) %>%
  mutate(pair = factor(pair, levels = c("W-C","CW-CC","WC-CC","WW-CC")))

# Define palettes
pal.gam <- c("#EE4B2B","#7985CB", "#9C0098", "#EE4B2B")
lty.gam <- c("dotdash","solid","solid","solid")

# Panel B: plot of differences between gam smooths for algal cover in Y1 and Y2
Fig5b <- ggplot(all.comp, aes(x = date, y = diff, group = pair, lty = pair, col = pair, fill = pair)) +
  geom_vline(aes(xintercept = ymd("2020-04-03")), col = "grey60", lwd = 0.8) +
  geom_vline(aes(xintercept = ymd("2019-08-27")), col = "grey70", lty = "dashed", lwd = 0.8) +
  geom_line(data = comp, col = "grey70")+
  geom_ribbon(data = comp, aes(ymin = lower, ymax = upper), col = NA, fill = "grey80",alpha = 0.2)+
  geom_ribbon(aes(ymin = lower, ymax = upper), col = NA, alpha = 0.2) +
  scale_color_manual(values = pal.gam) +
  scale_fill_manual(values = pal.gam) +
  scale_linetype_manual(values = lty.gam)+
  geom_line() +
  labs(x = "Date", y = 'Difference in algal cover (%)', fill = "Comparison",
       lty = "Comparison", col = "Comparison") +
  theme_classic() +
  geom_hline(aes(yintercept = 0)) +
  scale_x_date(date_labels = "%Y-%m", breaks=breaks, 
               limits = c(ymd("2019-05-03"),ymd("2021-02-24"))) +
  theme(plot.tag = element_text(face = "bold", size = 10),
        strip.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))
Fig5b

# Assemble multipanel figure with raw cover and difference smoothers
Fig5 <- (Fig5a / Fig5b) + plot_annotation(tag_levels = "a", tag_prefix = "(",
                                          tag_suffix = ")")
Fig5

# save figure
ggsave(Fig5, filename = "./figures/Fig5.pdf", device = cairo_pdf, 
       width = 16, height = 13, units = "cm")
