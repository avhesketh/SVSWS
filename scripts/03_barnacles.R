# Barnacle abundance and recruitment
# Amelia Hesketh January 2024

# load in packages
pkgs <- c("tidyverse","glmmTMB","DHARMa", "plotrix", 
          "lubridate", "emmeans","car", "patchwork")
lapply(pkgs, library, character.only = TRUE)
rm(pkgs)

# load in experimental design information
block_design <- read_csv("./raw_data/design/SVSHS_tilesetup.csv") %>% 
  select(tile_id, original_block, original_no, new_block, new_no)

# define palettes for use in later plotting
pal.trt <- c("#014779", "#EE4B2B", "#014779", "#7985CB", "#9C0098", "#EE4B2B", "grey50")
pch.trt <- c(1,1,16,16,16,16,16)
lty.trt <- c("longdash","longdash","solid","solid","solid","solid","solid","solid")


# read in the cleaned survey data
surveys <- read_csv("./clean_data/SVSWS_survey_clean.csv")

# Some surveys occurred over multiple days in a tide survey.
# Need to unify survey dates for ease of plotting... use last day of survey.
survey_numbers <- read_csv("./raw_data/design/SVSWS_survey_times.csv") %>% 
  group_by(survey_no) %>% 
  summarize(date = max(date))

# Year 1 barnacle data ... 
barnacles_y1 <- surveys %>% 
  left_join(survey_numbers) %>% 
  filter(species %in% c("Balanus_glandula", "Chthamalus_dalli") & date <= "2020-03-15") %>% 
  mutate(timesincestart = difftime(date, ymd("2019-04-12"), units = c("weeks")),
         count = as.integer(count)) %>% 
  mutate(treatment = trt_y1)

# Year 2 barnacle data ...
barnacles_y2 <- surveys %>% 
  left_join(survey_numbers) %>% 
  filter(species %in% c("Balanus_glandula", "Chthamalus_dalli") & date >= "2020-03-15" & (treatment %in% c("W","C"))==F) %>% 
  mutate(timesincestart = difftime(date, ymd("2019-04-12"), units = c("weeks")),
         count = as.integer(count))

# Join the two data frames into one ... 
barnacles_all <- barnacles_y1 %>% full_join(barnacles_y2) %>% 
  mutate(treatment = factor(treatment, levels = c("C","W","CC","CW","WC","WW"))) %>% 
  mutate(month = month(date)) %>% select(treatment, block, tile_id, date, survey_no, species, count, original_herb_trt) %>% 
  left_join(block_design) %>% 
  # add in information about tile number and block at the time of surveys from tile id
  mutate(tile_number = ifelse(date <= "2019-06-06", original_no, new_no),
         block = ifelse(date <= "2019-06-06", original_block, new_block)) %>% 
  select(!c(original_block, original_no, new_block, new_no))

###############
# Plotting barnacle recruitment and adult abundance

pal.trt.y2 <- c("#014779","#7985CB", "#9C0098", "#EE4B2B")

# Add recruitment data

# In year 1, all barnacles were recruits. Balanus peaked in recruiment in May 2019
# While C. dalli peaked in survey 2 ... filter these data for analyses & plots
recruitment_y1 <- barnacles_all %>% 
  filter((species == "Balanus_glandula" & survey_no == 1) |
           (species == "Chthamalus_dalli" & survey_no == 2)) %>% 
  mutate(period = "Year 1") %>% 
  mutate(herb = if_else(is.na(original_herb_trt), "control","herb"),
         species = if_else(species == "Balanus_glandula", "Balanus glandula", "Chthamalus dalli"))
           
# In year 2, recruits and adults were explicitly counted for most surveys
# These data are in a separate dataframe
recruitment <- read_csv("./raw_data/tile_surveys/SVSWS_barnacle_recruit.csv") %>%
  rename(new_block = block, new_no = number) %>% 
  left_join(block_design) %>% 
  mutate(species = if_else(species == "balanus", "Balanus glandula", "Chthamalus dalli")) %>% 
  select(new_block, new_no, tile_id, date, species, size, count, treatment,
         original_herb_trt) %>% 
  mutate(herb = if_else(is.na(original_herb_trt), "control","herb"))

recruitment_y2 <- recruitment %>% filter(size == "recruit" & date == "2020-06-04") %>% 
  mutate(period = "Year 2")

peak_recruitment <- recruitment_y1 %>% full_join(recruitment_y2) %>% 
  mutate(treatment = factor(treatment, levels = c("C", "W", "CC","CW","WC","WW", "Rock")),
         period = if_else(period == "Year 1", "2019", "2020"))

# Plot of peak recruitment for Balanus in 2019 and 2020
Fig2A <- ggplot(aes(x = treatment, y = count, col = treatment), 
                data = peak_recruitment %>% filter(species == "Balanus glandula")) +
  geom_boxplot(lwd = 0.4, outlier.color = NA) +
  geom_jitter(aes(pch = treatment), size = 0.7) +
  labs(y = expression(~italic("Balanus glandula")~"recruit abundance"), 
       x = "Treatment", col = "Treatment",
       pch = "Treatment") +
  scale_color_manual(values = pal.trt) +
  scale_shape_manual(values = pch.trt) +
  theme_classic() +
  facet_wrap(~period, scales = "free") +
  theme(plot.tag = element_text(face = "bold")) +
  ylim(c(-1,600))
Fig2A

# Same for C. dalli
Fig2B <- ggplot(aes(x = treatment, y = count, col = treatment), 
                data = peak_recruitment %>% filter(species == "Chthamalus dalli")) +
  geom_boxplot(lwd = 0.4, outlier.color = NA) +
  geom_jitter(aes(pch = treatment), size = 0.7) +
  labs(y = expression(~italic("Chthamalus dalli")~"recruit abundance"),
       x = "Treatment", col = "Treatment",
       pch = "Treatment") +
  scale_color_manual(values = pal.trt) +
  scale_shape_manual(values = pch.trt) +
  theme_classic() +
  facet_wrap(~period, scales = "free") +
  theme(plot.tag = element_text(face = "bold")) +
  ylim(c(-1,85))
Fig2B


# Now plot the abundance of adult barnacles remaining @ end of experiment

adult_barnacles <- recruitment %>% filter(size == "adult" & date == "2021-02-24") %>% 
  mutate(trt_y1 = substr(treatment ,1,1), trt_y2 = substr(treatment, 2,2))
  

Fig2C <- ggplot(adult_barnacles, aes(x = treatment, y = count, col = treatment)) +
  geom_boxplot(outlier.color = NA, lwd = 0.4) +
  geom_jitter(width = 0.25, size = 0.7) +
  scale_color_manual(values = pal.trt.y2) +
  facet_wrap(~species) +
  theme_classic() +
  theme(strip.text= element_text(face = "italic"),
        plot.tag = element_text(face = "bold"),
        legend.position = "none") +
  labs(y = "Number of adult barnacles", x = "Treatment", col = "Treatment",
       pch = "Treatment") +
  theme(plot.tag = element_text(face = "bold"))
Fig2C



Fig2 <- ((Fig2A / Fig2B) | Fig2C) + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") +
  plot_layout(widths = c(0.6,0.4))


png("./figures/Fig2.png", res = 700, width = 9, height = 6, units = "in")
Fig2
dev.off()


############### 
# Models of barnacle abundance and recruitment

# recruitment models
recruits <- peak_recruitment %>% mutate(trt_y1 = substr(treatment ,1,1), trt_y2 = substr(treatment, 2,2))

bal.rec.y1 <- glmmTMB(count ~ trt_y1 + (1|block), 
                   family = nbinom1(), 
                   data = recruits %>% filter(survey_no == 1 & species == "Balanus glandula"))
bal.rec.herb <- glmmTMB(count ~ trt_y1 + herb + (1|block), 
                        family = nbinom1(), 
                        data = recruits %>% filter(survey_no == 1 & species == "Balanus glandula"))

plot(simulateResiduals(bal.rec.y1)) # slight patterning indicating overdispersion, but doesn't violate any assumptions
summary(bal.rec.y1)
Anova(bal.rec.y1)

# AIC(bal.rec.y1, bal.rec.herb) # no substantial herbivore effect, exclude from model

bal.rec.y2 <- glmmTMB(count ~ trt_y1*trt_y2 + (1|new_block), 
                      family = nbinom1(), 
                      data = recruits %>% filter(date == "2020-06-04" & species == "Balanus glandula"))
bal.rec.herb <- glmmTMB(count ~ trt_y1*trt_y2 + herb + (1|new_block), 
                        family = nbinom1(), 
                        data = recruits %>% filter(date == "2020-06-04" & species == "Balanus glandula"))

plot(simulateResiduals(bal.rec.y2)) # slight patterning indicating overdispersion, but doesn't violate any assumptions
summary(bal.rec.y2)
Anova(bal.rec.y2, type = 2)

# AIC(bal.rec.y2, bal.rec.herb) # no substantial herbivore effect, exclude from model

# now for Chthamalus

cd.rec.y1 <- glmmTMB(count ~ trt_y1 + (1|block), 
                      family = nbinom1(), 
                      data = recruits %>% filter(survey_no == 2 & species == "Chthamalus dalli"))
cd.rec.herb <- glmmTMB(count ~ trt_y1 + herb + (1|block), 
                        family = nbinom1(), 
                        data = recruits %>% filter(survey_no == 2 & species == "Chthamalus dalli"))

plot(simulateResiduals(cd.rec.y1)) # slight patterning indicating overdispersion, but doesn't violate any assumptions
summary(cd.rec.y1)
Anova(cd.rec.y1)

# AIC(cd.rec.y1, cd.rec.herb) # no substantial herbivore effect, exclude from model

cd.rec.y2 <- glmmTMB(count ~ trt_y1*trt_y2 + (1|new_block), 
                      family = nbinom1(), 
                      data = recruits %>% filter(date == "2020-06-04" & species == "Chthamalus dalli"))
cd.rec.herb <- glmmTMB(count ~ trt_y1*trt_y2 + herb + (1|new_block), 
                        family = nbinom1(), 
                        data = recruits %>% filter(date == "2020-06-04" & species == "Chthamalus dalli"))

plot(simulateResiduals(cd.rec.y2)) # slight patterning indicating overdispersion, but doesn't violate any assumptions
summary(cd.rec.y2)
Anova(cd.rec.y2, type = 2)

# AIC(cd.rec.y2, cd.rec.herb) # no substantial herbivore effect, exclude from model

# now create model for final number of adults on tiles

# Balanus glandula

bg.adults.y2 <- glmmTMB(count ~ trt_y1*trt_y2 + (1|new_block), 
                        family = nbinom1(), 
                        data = adult_barnacles %>% filter(species == "Balanus glandula"))

plot(simulateResiduals(bg.adults.y2)) # slight patterning indicating overdispersion, but doesn't violate any assumptions
summary(bg.adults.y2)
Anova(bg.adults.y2, type = 2)

cd.adults.y2 <- glmmTMB(count ~ trt_y1*trt_y2 + (1|new_block), 
                        family = nbinom1(), 
                        data = adult_barnacles %>% filter(species == "Chthamalus dalli"))

plot(simulateResiduals(cd.adults.y2)) # slight patterning indicating overdispersion, but doesn't violate any assumptions
summary(cd.adults.y2)
Anova(cd.adults.y2, type = 2)

####################
# Plotting barnacle abundance through time (all surveys)
barnacles_summarized <- barnacles_all %>% group_by(date, species, survey_no, treatment) %>% 
  summarize(mean_abund = mean(count), se_abund = std.error(count)) %>% ungroup()

# add in zero data from when experiment was deployed
zero_point <- data.frame(date = as.Date("2019-04-12"))
zero_point$mean_abund <- 0
zero_point$se_abund <- 0
treatment <- c("C","W")
species <- c("Balanus_glandula", "Chthamalus_dalli")
# expand across all treatments and species
zero_point <- zero_point %>% expand_grid(treatment) %>% expand_grid(species)

# this dataframe will be used for plot of abundance over time
barnacles_summarized <- barnacles_summarized %>% full_join(zero_point) %>% 
  mutate(treatment = factor(treatment, levels = c("C","W","CC","CW","WC","WW")),
         species = if_else(species == "Balanus_glandula", "Balanus glandula",
                           "Chthamalus dalli"))

FigA6 <- ggplot(aes(x = date, y = mean_abund, col = treatment,
                      lty = treatment, shape = treatment), data = barnacles_summarized) +
  geom_line(lwd = 0.7) +
  geom_point(size = 1.2)+
  scale_color_manual(values = pal.trt) +
  geom_errorbar(aes(ymin = mean_abund - se_abund, ymax = mean_abund + se_abund, lty = NULL), width = 7) +
  scale_linetype_manual(values = c("twodash","twodash","solid","solid","solid","solid")) +
  scale_shape_manual(values = c(1,1,16,16,16,16)) +
  labs(x = "Date", 
       y = "Barnacle abundance",
       col = "Treatment", lty = "Treatment", pch = "Treatment") +
  geom_vline(aes(xintercept = as.Date("2020-04-07")), linetype = "dotted", col = "grey30", lwd = 0.8) +
  theme_classic() +
  theme(strip.text = element_text(face = "italic")) +
  facet_grid(rows = vars(species), scales = "free")#+
FigA6

png("./figures/FigA6.png", res = 700, width = 8.5, height = 5, units = "in")
FigA6
dev.off()
