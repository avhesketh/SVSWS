# Barnacle abundance and recruitment
# Amelia Hesketh January 2024

# load in packages
pkgs <- c("tidyverse","glmmTMB","DHARMa", "plotrix", 
          "lubridate", "emmeans","car", "patchwork")
lapply(pkgs, library, character.only = TRUE)
rm(pkgs)

# load in experimental design information
block_design <- read_csv("./raw_data/design/SVSWS_tilesetup.csv") %>% 
  select(tile_id, original_block, original_no, new_block, new_no, treatment, first_herb_trt, 
         second_herb_trt)


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

block_design2 <- block_design %>% 
  select(-treatment)

# Join the two data frames into one ... 
barnacles_all <- barnacles_y1 %>% full_join(barnacles_y2) %>% 
  mutate(treatment = factor(treatment, levels = c("C","W","CC","CW","WC","WW"))) %>% 
  mutate(month = month(date)) %>% select(treatment, block, tile_id, date, 
                                         survey_no, species, count, 
                                         first_herb_trt, second_herb_trt) %>% 
  left_join(block_design2, by = c("tile_id", "first_herb_trt","second_herb_trt")) %>% 
  # add in information about tile number and block at the time of surveys from tile id
  mutate(tile_number = ifelse(date <= "2019-06-06", original_no, new_no),
         block = ifelse(date <= "2019-06-06", original_block, new_block)) %>% 
  select(!c(original_block, original_no, new_block, new_no))

###############
# Plotting barnacle recruitment and adult abundance

pal.trt.y2 <- c("#014779","#7985CB", "#9C0098", "#EE4B2B")
pal.trt.y1 <- c("#014779", "#EE4B2B")

# Add recruitment data

# In year 1, all barnacles were recruits. Balanus peaked in recruiment in May 2019
# While C. dalli peaked in survey 2 ... filter these data for analyses & plots
recruitment_y1 <- barnacles_all %>% 
  filter((species == "Balanus_glandula" & survey_no == 1) |
           (species == "Chthamalus_dalli" & survey_no == 2)) %>% 
  mutate(period = "Year 1") %>% 
  mutate(herb_trt = if_else(is.na(first_herb_trt), "control","grazers"),
         species = if_else(species == "Balanus_glandula", "Balanus glandula", "Chthamalus dalli"))
           
# In year 2, recruits and adults were explicitly counted for most surveys
# These data are in a separate dataframe
recruitment <- read_csv("./raw_data/tile_surveys/SVSWS_barnacle_recruit.csv") %>%
  rename(new_block = block, new_no = number) %>% 
  left_join(block_design) %>% 
  mutate(species = if_else(species == "balanus", "Balanus glandula", "Chthamalus dalli")) %>% 
  select(new_block, new_no, tile_id, date, species, size, count, treatment,
         second_herb_trt) %>% 
  mutate(herb_trt = if_else(is.na(second_herb_trt), "control",
                                     "grazers"))

recruitment_y2 <- recruitment %>% filter(size == "recruit" & date == "2020-06-04") %>% 
  mutate(period = "Year 2")

peak_recruitment <- recruitment_y1 %>% full_join(recruitment_y2) %>% 
  mutate(treatment = factor(treatment, levels = c("C", "W", "CC","CW","WC","WW", "Rock")),
         period = if_else(period == "Year 1", "2019", "2020"))

############### 
# Models of barnacle abundance and recruitment

# recruitment models
recruits <- peak_recruitment %>% mutate(trt_y1 = substr(treatment ,1,1), trt_y2 = substr(treatment, 2,2))

# since herbivore treatments were not added by this point, no need to test its inclusion within the model
bal.rec.y1 <- glmmTMB(count ~ trt_y1 + (1|block), 
                      family = nbinom1(), 
                      data = recruits %>% filter(survey_no == 1 & species == "Balanus glandula"))

plot(simulateResiduals(bal.rec.y1)) # slight patterning indicating overdispersion, but doesn't violate any assumptions
summary(bal.rec.y1)
Anova(bal.rec.y1, type = 2)


# Year 2
bal.rec.y2 <- glmmTMB(count ~ trt_y1*trt_y2 + (1|new_block), 
                      family = nbinom1(), 
                      data = recruits %>% filter(date == "2020-06-04" & species == "Balanus glandula"))
plot(simulateResiduals(bal.rec.y2)) # no violations of assumptions
summary(bal.rec.y2)
Anova(bal.rec.y2, type = 3, contrasts=list(trt_y1 = "contr.sum", trt_y2 = "contr.sum"))

emm2.bg <- emmeans(bal.rec.y2, specs = c("trt_y1","trt_y2"))
contrast(emm2.bg, method = "pairwise", adjust = "tukey")

# here, just creating a dataframe with labels indicating significance 
# of emmeans comparisons
labels.3a <- as.data.frame(cbind(
  c("C","W","CC","CW","WC","WW"),
  c("2019","2019","2020","2020","2020","2020"),
  c(625,625,300,300,300,300),
  c("a","a","c","d","c","d"))
)
colnames(labels.3a) <- c("treatment","period", "count","label")
labels.3a <- labels.3a %>% mutate(count = as.numeric(count))

# now for Chthamalus

# Year 1: may be affected by early herbivore manipulations -- need to rule this effect out
cd.rec.herb <- glmmTMB(count ~ trt_y1 + herb_trt + (1|block), 
                       family = nbinom1(), 
                       data = recruits %>% filter(survey_no == 2 & species == "Chthamalus dalli"))
plot(simulateResiduals(cd.rec.herb)) 
summary(cd.rec.herb)
Anova(cd.rec.herb, type = 2) 
# herbivore treatment did not exert a significant effect; drop this

cd.rec.y1 <- glmmTMB(count ~ trt_y1 + (1|block), 
                     family = nbinom1(), 
                     data = recruits %>% filter(survey_no == 2 & species == "Chthamalus dalli"))

plot(simulateResiduals(cd.rec.y1)) # slight patterning indicating overdispersion, but doesn't violate any assumptions
summary(cd.rec.y1)
Anova(cd.rec.y1, type = 2)


# Year 2
cd.rec.y2 <- glmmTMB(count ~ trt_y1*trt_y2 + (1|new_block), 
                     family = nbinom1(), 
                     data = recruits %>% filter(date == "2020-06-04" & species == "Chthamalus dalli"))

plot(simulateResiduals(cd.rec.y2)) # slight patterning indicating overdispersion, but doesn't violate any assumptions
summary(cd.rec.y2)
Anova(cd.rec.y2, type = 3, contrasts=list(trt_y1 = "contr.sum", trt_y2 = "contr.sum"))

emm2.cd <- emmeans(cd.rec.y2, specs = c("trt_y1","trt_y2"))
contrast(emm2.cd, method = "pairwise", adjust = "tukey")

# comparison labels for chthamalus recruitment panel
labels.3b <- as.data.frame(cbind(
  c("C","W","CC","CW","WC","WW"),
  c("2019","2019","2020","2020","2020","2020"),
  c(75,40,89,51,73,51),
  c("a","b","c","d","cd","d"))
)
colnames(labels.3b) <- c("treatment","period", "count","label")
labels.3b <- labels.2b %>% mutate(count = as.numeric(count))

############### 
# Visualizing recruitment data

# Plot of peak recruitment for Balanus in 2019 and 2020
Fig3A <- ggplot(aes(x = treatment, y = count, col = treatment), 
                data = peak_recruitment %>% filter(species == "Balanus glandula")) +
  geom_boxplot(lwd = 0.4, outlier.color = NA) +
  geom_jitter(aes(pch = treatment), height = 0, size = 0.7) +
  labs(y = expression(~italic("Balanus glandula")~"recruit abundance"), 
       x = "Treatment", col = "Treatment",
       pch = "Treatment") +
  scale_color_manual(values = pal.trt) +
  scale_shape_manual(values = pch.trt) +
  theme_classic() +
  facet_wrap(~period, scales = "free") +
  theme(plot.tag = element_text(face = "bold", size = 10),
        strip.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8),
        legend.position = "none") +
  ylim(c(-1,625)) +
  geom_text(data = labels.3a, aes(label = label), size = 3, col = "black", fontface = "bold")
Fig3A

# Same for C. dalli
Fig3B <- ggplot(aes(x = treatment, y = count, col = treatment), 
                data = peak_recruitment %>% filter(species == "Chthamalus dalli")) +
  geom_boxplot(lwd = 0.4, outlier.color = NA) +
  geom_jitter(aes(pch = treatment), height = 0, size = 0.7) +
  labs(y = expression(~italic("Chthamalus dalli")~"recruit abundance"),
       x = "Treatment", col = "Treatment",
       pch = "Treatment") +
  scale_color_manual(values = pal.trt) +
  scale_shape_manual(values = pch.trt) +
  theme_classic() +
  facet_wrap(~period, scales = "free") +
  theme(plot.tag = element_text(face = "bold", size = 10),
        strip.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8),
        legend.position = "none") +
  ylim(c(0,90)) +
  geom_text(data = labels.3b, aes(label = label), size = 3, col = "black", fontface = "bold")
Fig3B

###############
# Modeling adult barnacle data

# create dataframe for summer
adult_barnacles_presummer <- barnacles_all %>% 
  filter(date == "2020-03-15") %>% 
  mutate(treatment = substr(treatment,1,1),
         herb_trt = if_else(is.na(second_herb_trt), "control", "grazers"),
         species = if_else(species == "Balanus_glandula", "Balanus glandula", "Chthamalus dalli"),
         period = "Winter 2020") %>% 
  mutate(treatment = factor(treatment, levels = c("C","W","CC","CW","WC","WW")))

adult_barnacles_postsummer <- recruitment %>% filter(size == "adult" & date == "2021-02-24") %>% 
  mutate(trt_y1 = substr(treatment ,1,1), trt_y2 = substr(treatment, 2,2),
         herb_trt = if_else(is.na(second_herb_trt), "control", "grazers"),
         period = "Winter 2021") %>% 
  mutate(treatment = factor(treatment, levels = c("C","W","CC","CW","WC","WW"))) %>% 
  rename(block = new_block, number = new_no)

balanus_adults <- adult_barnacles_presummer %>% 
  full_join(adult_barnacles_postsummer) %>% filter(species == "Balanus glandula")

chthamalus_adults <- adult_barnacles_presummer %>% 
  full_join(adult_barnacles_postsummer) %>% filter(species == "Chthamalus dalli")


# Balanus glandula

# End of year 1: test effect of grazer manipulations first

bg.adults.y1.herb <- glmmTMB(count ~ treatment + herb_trt + (1|block), 
                             family = nbinom1(), 
                             data = balanus_adults %>% filter(period == "Winter 2020"))
plot(simulateResiduals(bg.adults.y1.herb)) # slight patterning indicating overdispersion, but doesn't violate any assumptions
summary(bg.adults.y1.herb)
print(Anova(bg.adults.y1.herb, type = 2, digits =7))
# grazer inclusion not significant

bg.adults.y1 <- glmmTMB(count ~ treatment + (1|block), 
                        family = nbinom1(), 
                        data = balanus_adults %>% filter(period == "Winter 2020"))
plot(simulateResiduals(bg.adults.y1)) # slight patterning indicating overdispersion, but doesn't violate any assumptions
summary(bg.adults.y1)
print(Anova(bg.adults.y1, type = 2, digits =7))

# End of year 2

bg.adults.y2 <- glmmTMB(count ~ trt_y1*trt_y2 + (1|block), 
                        family = nbinom1(), 
                        data = balanus_adults %>% filter(period == "Winter 2021"))
plot(simulateResiduals(bg.adults.y2)) # slight patterning indicating overdispersion, but doesn't violate any assumptions
summary(bg.adults.y2)
Anova(bg.adults.y2, type = 3, contrasts=list(trt_y1 = "contr.sum", trt_y2 = "contr.sum"))

emm2.bga <- emmeans(bg.adults.y2, specs = c("trt_y1","trt_y2"))
contrast(emm2.bga, method = "pairwise", adjust = "tukey")

# significance labels for adult Bg abundance plot
labels.3c <- as.data.frame(cbind(
  c("C","W","CC","CW","WC","WW"),
  c("Winter 2020","Winter 2020","Winter 2021","Winter 2021","Winter 2021","Winter 2021"),
  c(118,33,132,75,75,32),
  c("a","b","c","cd","cd","d"))
)
colnames(labels.3c) <- c("treatment","period", "count","label")
labels.3c <- labels.3c %>% mutate(count = as.numeric(count))


# Adult Chthamalus dalli abundance

# Year 1: test the effect of grazer treatments first
cd.adults.y1.herb <- glmmTMB(count ~ treatment + herb_trt + (1|block), 
                        family = nbinom1(), 
                        data = chthamalus_adults %>% filter(period == "Winter 2020"))
plot(simulateResiduals(cd.adults.y1.herb)) # slight patterning indicating overdispersion, but doesn't violate any assumptions
summary(cd.adults.y1.herb)
Anova(cd.adults.y1.herb, type = 2)
# grazer manipulation term not significant; drop this

cd.adults.y1 <- glmmTMB(count ~ treatment + (1|block), 
                             family = nbinom1(), 
                             data = chthamalus_adults %>% filter(period == "Winter 2020"))
plot(simulateResiduals(cd.adults.y1)) # slight patterning indicating overdispersion, but doesn't violate any assumptions
summary(cd.adults.y1)
Anova(cd.adults.y1, type = 2)

# Year 2
cd.adults.y2 <- glmmTMB(count ~ trt_y1*trt_y2 + (1|block), 
                        family = nbinom1(), 
                        data = chthamalus_adults %>% filter(period == "Winter 2021"))

plot(simulateResiduals(cd.adults.y2)) # slight patterning indicating overdispersion, but doesn't violate any assumptions
summary(cd.adults.y2)
Anova(cd.adults.y2, type = 3, contrasts=list(trt_y1 = "contr.sum", trt_y2 = "contr.sum"))

emm2.cda <- emmeans(cd.adults.y2, specs = c("trt_y1", "trt_y2"))
contrast(emm2.cda, method = "pairwise", adjust = "tukey")

# labels for C. dalli adult abundance plot
labels.3d <- as.data.frame(cbind(
  c("C","W","CC","CW","WC","WW"),
  c("Winter 2020","Winter 2020","Winter 2021","Winter 2021","Winter 2021","Winter 2021"),
  c(15,15,110,110,110,110),
  c("a","a","c","c","c","c"))
)
colnames(labels.3d) <- c("treatment","period", "count","label")
labels.3d <- labels.3d %>% mutate(count = as.numeric(count))

###############
# Figure panels for adult barnacle abundance

# Plotting abundance of adult B. glandula before y2 summer and at the end of the experiment

Fig3C <- ggplot(balanus_adults, aes(x = treatment, y = count, col = treatment)) +
  geom_boxplot(outlier.color = NA, lwd = 0.4) +
  geom_jitter(aes(pch = treatment), height = 0, size = 0.7) +
  scale_color_manual(values = pal.trt) +
  scale_shape_manual(values = pch.trt) +
  facet_wrap(~period, scales = "free_x") +
  theme_classic() +
  theme(
        plot.tag = element_text(face = "bold"),
        legend.position = "none") +
  labs(y = expression(~italic("Balanus glandula")~"adult abundance"), x = "Treatment", col = "Treatment",
       pch = "Treatment") +
  theme(plot.tag = element_text(face = "bold", size = 10),
        strip.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8),
        legend.position = "none") +
  scale_y_continuous(breaks = c(0,25,50,75,100,125)) +
  geom_text(data = labels.3c, aes(label = label), size = 3, col = "black", fontface = "bold")
Fig3C

# Plotting abundance of adult C. dalli before y2 summer and at the end of the experiment
Fig3D <- ggplot(chthamalus_adults, aes(x = treatment, y = count, col = treatment)) +
  geom_boxplot(outlier.color = NA, lwd = 0.4) +
  geom_jitter(aes(pch = treatment), height = 0, size = 0.7) +
  scale_color_manual(values = pal.trt) +
  scale_shape_manual(values = pch.trt) +
  facet_wrap(~period, scales = "free_x") +
  theme_classic() +
  theme(
        plot.tag = element_text(face = "bold"),
        legend.position = "none") +
  labs(y = expression(~italic("Chthamalus dalli")~"adult abundance"), x = "Treatment", col = "Treatment",
       pch = "Treatment") +
  theme(plot.tag = element_text(face = "bold", size = 10),
        strip.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8),
        legend.position = "none") +
  scale_y_continuous(breaks = c(0,25,50,75,100)) +
  geom_text(data = labels.3d, aes(label = label), size = 3, col = "black", fontface = "bold")
Fig3D

# Assemble multiplanel figure
Fig3 <- ((Fig3A / Fig3B) | (Fig3C / Fig3D)) + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")

# save figure
ggsave(Fig3, filename = "./figures/Fig3.pdf", device = cairo_pdf, 
       width = 16, height = 14, units = "cm")


#################### 
# Appendix plots of barnacle data

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

# setting dates where labels appear on the x axis
breaks <- c(ymd("2019-04-01"), ymd("2019-08-01"),ymd("2019-12-01"),
            ymd("2020-04-01"),ymd("2020-08-01"),ymd("2020-12-01"),
            ymd("2021-04-01"))

# Plot of barnacle (B. glandula and C. dalli) abundance over time in each treatment 
FigS6 <- ggplot(aes(x = date, y = mean_abund, col = treatment,
                      lty = treatment, shape = treatment), data = barnacles_summarized) +
  geom_vline(aes(xintercept = ymd("2020-04-03")), col = "grey60", lwd = 0.8) +
  geom_vline(aes(xintercept = ymd("2019-08-27")), col = "grey70", lty = "dashed", lwd = 0.8) +
  geom_line(lwd = 0.7) +
  geom_point(size = 1.2) +
  scale_color_manual(values = pal.trt) +
  geom_errorbar(aes(ymin = mean_abund - se_abund, ymax = mean_abund + se_abund, lty = NULL), width = 7) +
  scale_linetype_manual(values = c("twodash","twodash","solid","solid","solid","solid")) +
  scale_shape_manual(values = c(1,1,16,16,16,16)) +
  labs(x = "Date", 
       y = "Barnacle abundance",
       col = "Treatment", lty = "Treatment", pch = "Treatment") +
  theme_classic() +
  theme(strip.text = element_text(face = "italic")) +
  facet_grid(rows = vars(species), scales = "free") +
  scale_x_date(date_labels = "%Y-%m", breaks=breaks)
FigS6

png("./figures/FigS6.png", res = 700, width = 8.5, height = 5, units = "in")
FigS6
dev.off()

############# 
# Supplemental size & mortality analysis using data extracted from photos

# 1: Mortality

# Isolate number of live barnacles on tiles (already counted) from Balanus dataset
# Photos are from August 27, but the closest live counts we have are from two weeks earlier (August 14). 
# Mortality should be minimal between low tides.
live_balanus <- barnacles_all %>% 
  filter(species == "Balanus_glandula") %>% 
  mutate(date = if_else(date == ymd("2019-08-14"),ymd("2019-08-27"),date)) %>% 
  select(-treatment)

# Load the mortality data

mortality <- read_csv("./raw_data/tile_surveys/SVSWS_barnacle_mortality.csv") %>% 
  rename(new_block = block, new_no = number) %>% 
  # join with other pertinent information (e.g. tile id)
  left_join(block_design) %>% 
  mutate(treatment = factor(treatment, levels = c("C","W","CC", "CW","WC","WW"))) %>% 
  select(date, treatment, new_block, new_no, number_dead) %>% 
  rename(block = new_block, tile_number = new_no) %>% 
  # join to live Balanus counts
  left_join(live_balanus, by = c("block","tile_number","date")) %>% 
  select(date, treatment, block, tile_number, number_dead, count) %>% 
  # calculate a proportional mortality
  mutate(mortality_prop = number_dead/(number_dead + count),
         treatment = if_else(date < "2020-03-15", substr(treatment, 1,1), treatment))

# summarize for plotting: calculate mean & std error for treatment groups on each date
mort_summary <- mortality %>% 
  group_by(date, treatment) %>% 
  summarize(mean_mort = mean(mortality_prop,na.rm =T), se_mort = std.error(mortality_prop, na.rm = T),
            mean_dead = mean(number_dead, na.rm = T), se_dead = std.error(number_dead, na.rm = T)) %>% 
  mutate(treatment = factor(treatment, levels = c("C","W","CC",
                                                  "CW","WC","WW")))

# Mortality supplemental figure: proportional mortality
FigS7a <-ggplot(mortality, aes(x = treatment, y = mortality_prop, col = treatment, pch = treatment)) +
  facet_grid(cols = vars(date), scales = "free_x") +
  geom_jitter(alpha = 0.2, height=0, width = 0.1) +
  geom_point(data = mort_summary, aes(y = mean_mort), size = 2, stroke = 1.5) +
  geom_errorbar(data = mort_summary, aes(y = mean_mort,
                                         ymax = mean_mort + se_mort,
                                         ymin = mean_mort - se_mort),
                col = "grey30", width = 0.2, linewidth = 0.5) +
  scale_color_manual(values = pal.trt) +
  scale_shape_manual(values = pch.trt)+
  theme_classic() +
  labs(x = "Date", y = "Proportion mortality",col = "Treatment", pch = "Treatment") +
  theme(axis.text = element_text(size =8),
        axis.title = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        plot.tag = element_text(size = 10, face= "bold"),
        legend.position = "none")

# Additional figure showing absolute number of dead barnacles (density of empty tests)
FigS7b <-ggplot(mortality, aes(x = treatment, y = number_dead, col = treatment, pch = treatment)) +
  facet_grid(cols = vars(date), scales = "free_x") +
  geom_jitter(alpha = 0.2, height=0, width = 0.1) +
  geom_point(data = mort_summary, aes(y = mean_dead), size = 2, stroke = 1.5) +
  geom_errorbar(data = mort_summary, aes(y = mean_dead,
                                         ymax = mean_dead + se_dead,
                                         ymin = mean_dead - se_dead),
                col = "grey30", width = 0.2, linewidth = 0.5) +
  scale_color_manual(values = pal.trt) +
  scale_shape_manual(values = pch.trt)+
  theme_classic() +
  labs(x = "Date", y = expression("Number of dead"~italic("Balanus glandula")),
       col = "Treatment", pch = "Treatment") +
  theme(axis.text = element_text(size =8),
        axis.title = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        plot.tag = element_text(size = 10, face= "bold"))

# Put it together in patchwork
FigS7 <- (FigS7a / FigS7b) + plot_layout(guides = "collect") + 
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")

# save output
ggsave(FigS7, filename = "./figures/FigS7.png",
       width = 12, height = 7, units = "cm", scale = 2)

# modeling mortality

# need one dataframe for first year, one for second
mort_y1 <- mortality %>% filter(date == "2019-08-27")
mort_y2 <- mortality %>% filter(date != "2019-08-27") %>% 
  mutate(trt_y1 = factor(substr(treatment, 1,1)),
         trt_y2 = factor(substr(treatment,2,2)),
         date = as.factor(date))

# in y1, only one data point -- simple analysis, using Tweedie distribution to help with zero inflation
mort.y1 <- glmmTMB(mortality_prop*100 ~ treatment + (1|block), 
                   family = tweedie(),
                   data = mort_y1)
plot(simulateResiduals(mort.y1)) #assumptions met
summary(mort.y1)
Anova(mort.y1, type = 2)

# in y2, multiple data points, so including date important (just as an additive effect)
mort.y2 <- glmmTMB(mortality_prop*100 ~ trt_y1*trt_y2 + date + (1|block/tile_number), 
                   family = tweedie(),
                   data = mort_y2)
plot(simulateResiduals(mort.y2)) #assumptions met
acf(residuals(mort.y2)) # autocorrelation not an issue despite repeated measures
summary(mort.y2)
Anova(mort.y2, type = 3, contrasts = list(trt_y1 = "contr.sum", trt_y2 = "contr.sum",
                                          date = "contr.sum"))

emm2.mort <- emmeans(mort.y2, specs = c("trt_y1", "trt_y2"))
contrast(emm2.mort, method = "pairwise", adjust = "tukey")


## 2: Size data

# Read in size data extracted from photographs
size <- read_csv("./raw_data/tile_surveys/SVSWS_barnacle_size.csv") %>% 
  mutate(treatment = factor(treatment, levels = c("C","W","CC", "CW","WC","WW")),
         time_since_start = difftime(date, "2019-04-12"),
         date = factor(date))  %>%  filter(status=="alive")

# summarize size data for plots
size_summary <- size %>% 
  group_by(date, treatment) %>% 
  summarize(mean_size = mean(basal_diameter_mm), se_size = std.error(basal_diameter_mm))

# plot density of barnacles at each time point separately
FigS8 <- ggplot(size, aes(x = treatment, y = basal_diameter_mm, col = treatment, pch = treatment)) +
  facet_grid(cols = vars(date), scales = "free_x") +
  geom_jitter(alpha = 0.07, height=0, width = 0.2) +
  geom_point(data = size_summary, aes(y = mean_size), alpha = 1, size = 2, stroke = 1.5) +
  geom_errorbar(data = size_summary, aes(y = mean_size,
                                         ymax = mean_size + se_size,
                                         ymin = mean_size - se_size),
                col = "grey20", width = 0.2, linewidth = 0.5) +
  scale_color_manual(values = pal.trt) +
  scale_shape_manual(values = pch.trt)+
  geom_hline(aes(yintercept = 5), lty = "dashed") +
  theme_classic() +
  labs(x = "Date", y = "Basal diameter (mm)",
       col = "Treatment", pch = "Treatment") +
  theme(axis.text = element_text(size =8),
        axis.title = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        plot.tag = element_text(size = 10, face= "bold"))

# save figure file
ggsave(FigS8, filename = "./figures/FigS8.png",
       width = 8, height = 5, units = "cm", scale = 2)


## Analysis of barnacle size

# break into two dataframes, one for each year to analyse separately
size_y1 <- size %>% filter(date == "2019-08-27")
size_y2 <- size %>% filter(date != "2019-08-27") %>% 
  mutate(trt_y1 = factor(substr(treatment, 1,1)),
         trt_y2 = factor(substr(treatment,2,2)),
         date = as.factor(date),
         time_since_start = difftime( date, "2019-04-12"))

# Simple mixed effects linear model, including random effect for individual tile
# log transform to better meet assumptions of normality for highly dispersed response variable
size.y1 <- lmer(log(basal_diameter_mm) ~ treatment + (1|block/number), data = size_y1)

plot(simulateResiduals(size.y1)) # everything looks great re:assumptions of normality

summary(size.y1)
Anova(size.y1, type = 2, test.statistic = "F")

# tried a model with date as an additive effect, but AIC drops a lot if it's an interactive effect
size.y2 <- lmer(log(basal_diameter_mm) ~ trt_y1*trt_y2*date + (1|block/number),
                   data = size_y2)

plot(simulateResiduals(size.y2)) # fails the KS test (distribution); however, lots of sample points and 
# looks like this is failing because of underdispersion -- proceed
plot(residuals(size.y2)) # more clear here: no substantial patterns in these residuals
acf(residuals(size.y2)) # autocorrelation looks ok

summary(size.y2)
Anova(size.y2, type = 3, contrasts=list(trt_y1 = "contr.sum", trt_y2 = "contr.sum",
                                        date = "contr.sum"), test.statistic = "F")

