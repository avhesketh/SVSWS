# Barnacle abundance and recruitment
# Amelia Hesketh January 2024

# load in packages
pkgs <- c("tidyverse","glmmTMB","DHARMa", "plotrix", 
          "lubridate", "emmeans","car", "patchwork")
lapply(pkgs, library, character.only = TRUE)
rm(pkgs)

# load in experimental design information
block_design <- read_csv("./raw_data/design/SVSWS_tilesetup.csv") %>% 
  select(tile_id, original_block, original_no, new_block, new_no, treatment, original_herb_trt)

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
  mutate(month = month(date)) %>% select(treatment, block, tile_id, date, survey_no, species, count, original_herb_trt) %>% 
  left_join(block_design2, by = c("tile_id", "original_herb_trt")) %>% 
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
  mutate(original_herb_trt = if_else(is.na(original_herb_trt), "control",
                                     original_herb_trt))

recruitment_y2 <- recruitment %>% filter(size == "recruit" & date == "2020-06-04") %>% 
  mutate(period = "Year 2")

peak_recruitment <- recruitment_y1 %>% full_join(recruitment_y2) %>% 
  mutate(treatment = factor(treatment, levels = c("C", "W", "CC","CW","WC","WW", "Rock")),
         period = if_else(period == "Year 1", "2019", "2020"))

############### 
# Models of barnacle abundance and recruitment

# recruitment models
recruits <- peak_recruitment %>% mutate(trt_y1 = substr(treatment ,1,1), trt_y2 = substr(treatment, 2,2))

bal.rec.y1 <- glmmTMB(count ~ trt_y1 + (1|block), 
                      family = nbinom1(), 
                      data = recruits %>% filter(survey_no == 1 & species == "Balanus glandula"))
bal.rec.herb <- glmmTMB(count ~ trt_y1 + original_herb_trt + (1|block), 
                        family = nbinom1(), 
                        data = recruits %>% filter(survey_no == 1 & species == "Balanus glandula"))

AIC(bal.rec.y1, bal.rec.herb) # no substantial herbivore effect, exclude from model

plot(simulateResiduals(bal.rec.y1)) # slight patterning indicating overdispersion, but doesn't violate any assumptions
summary(bal.rec.y1)
Anova(bal.rec.y1)

bal.rec.y2 <- glmmTMB(count ~ trt_y1*trt_y2 + (1|new_block), 
                      family = nbinom1(), 
                      data = recruits %>% filter(date == "2020-06-04" & species == "Balanus glandula"))
bal.rec.herb <- glmmTMB(count ~ trt_y1*trt_y2 + original_herb_trt + (1|new_block), 
                        family = nbinom1(), 
                        data = recruits %>% filter(date == "2020-06-04" & species == "Balanus glandula"))

AIC(bal.rec.y2, bal.rec.herb) # no substantial herbivore effect, exclude from model

plot(simulateResiduals(bal.rec.y2)) # slight patterning indicating overdispersion, but doesn't violate any assumptions
summary(bal.rec.y2)
Anova(bal.rec.y2, type = 3)

emm2.bg <- emmeans(bal.rec.y2, specs = c("trt_y1","trt_y2"))
contrast(emm2.bg, method = "pairwise", adjust = "tukey")

# here, just creating a dataframe with labels indicating significance 
# of emmeans comparisons
labels.2a <- as.data.frame(cbind(
  c("C","W","CC","CW","WC","WW"),
  c("2019","2019","2020","2020","2020","2020"),
  c(625,625,300,300,300,300),
  c("a","a","c","d","c","d"))
)
colnames(labels.2a) <- c("treatment","period", "count","label")
labels.2a <- labels.2a %>% mutate(count = as.numeric(count))

# now for Chthamalus

cd.rec.y1 <- glmmTMB(count ~ trt_y1 + (1|block), 
                     family = nbinom1(), 
                     data = recruits %>% filter(survey_no == 2 & species == "Chthamalus dalli"))
cd.rec.herb <- glmmTMB(count ~ trt_y1 + original_herb_trt + (1|block), 
                       family = nbinom1(), 
                       data = recruits %>% filter(survey_no == 2 & species == "Chthamalus dalli"))

AIC(cd.rec.y1, cd.rec.herb) # no substantial herbivore effect, exclude from model

plot(simulateResiduals(cd.rec.y1)) # slight patterning indicating overdispersion, but doesn't violate any assumptions
summary(cd.rec.y1)
Anova(cd.rec.y1)


cd.rec.y2 <- glmmTMB(count ~ trt_y1*trt_y2 + (1|new_block), 
                     family = nbinom1(), 
                     data = recruits %>% filter(date == "2020-06-04" & species == "Chthamalus dalli"))
cd.rec.herb <- glmmTMB(count ~ trt_y1*trt_y2 + original_herb_trt + (1|new_block), 
                       family = nbinom1(), 
                       data = recruits %>% filter(date == "2020-06-04" & species == "Chthamalus dalli"))

AIC(cd.rec.y2, cd.rec.herb) # no substantial herbivore effect, exclude from model

plot(simulateResiduals(cd.rec.y2)) # slight patterning indicating overdispersion, but doesn't violate any assumptions
summary(cd.rec.y2)
Anova(cd.rec.y2, type = 3)

emm2.cd <- emmeans(cd.rec.y2, specs = c("trt_y1","trt_y2"))
contrast(emm2.cd, method = "pairwise", adjust = "tukey")

# comparison labels for chthamalus recruitment panel
labels.2b <- as.data.frame(cbind(
  c("C","W","CC","CW","WC","WW"),
  c("2019","2019","2020","2020","2020","2020"),
  c(75,40,89,51,73,51),
  c("a","b","c","d","cd","d"))
)
colnames(labels.2b) <- c("treatment","period", "count","label")
labels.2b <- labels.2b %>% mutate(count = as.numeric(count))

############### 
# Visualizing recruitment data

# Plot of peak recruitment for Balanus in 2019 and 2020
Fig2A <- ggplot(aes(x = treatment, y = count, col = treatment), 
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
  geom_text(data = labels.2a, aes(label = label), size = 3, col = "black", fontface = "bold")
Fig2A

# Same for C. dalli
Fig2B <- ggplot(aes(x = treatment, y = count, col = treatment), 
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
  geom_text(data = labels.2b, aes(label = label), size = 3, col = "black", fontface = "bold")
Fig2B

###############
# Modeling adult barnacle data

# create dataframe for summe=
adult_barnacles_presummer <- barnacles_all %>% 
  filter(date == "2020-03-15") %>% 
  mutate(treatment = substr(treatment,1,1),
         original_herb_trt = if_else(is.na(original_herb_trt), "control", original_herb_trt),
         species = if_else(species == "Balanus_glandula", "Balanus glandula", "Chthamalus dalli"),
         period = "Winter 2020") %>% 
  mutate(treatment = factor(treatment, levels = c("C","W","CC","CW","WC","WW")))

adult_barnacles_postsummer <- recruitment %>% filter(size == "adult" & date == "2021-02-24") %>% 
  mutate(trt_y1 = substr(treatment ,1,1), trt_y2 = substr(treatment, 2,2),
         original_herb_trt = if_else(is.na(original_herb_trt), "control",original_herb_trt),
         period = "Winter 2021") %>% 
  mutate(treatment = factor(treatment, levels = c("C","W","CC","CW","WC","WW"))) %>% 
  rename(block = new_block, number = new_no)

balanus_adults <- adult_barnacles_presummer %>% 
  full_join(adult_barnacles_postsummer) %>% filter(species == "Balanus glandula")

chthamalus_adults <- adult_barnacles_presummer %>% 
  full_join(adult_barnacles_postsummer) %>% filter(species == "Chthamalus dalli")


# Balanus glandula

# Year 1
bg.adults.y1 <- glmmTMB(count ~ treatment + (1|block), 
                        family = nbinom1(), 
                        data = balanus_adults %>% filter(period == "Winter 2020"))
bg.adults.y1.herb <- glmmTMB(count ~ treatment + original_herb_trt + (1|block), 
                        family = nbinom1(), 
                        data = balanus_adults %>% filter(period == "Winter 2020"))

AIC(bg.adults.y1, bg.adults.y1.herb) # herbivore effect not significant

plot(simulateResiduals(bg.adults.y1)) # slight patterning indicating overdispersion, but doesn't violate any assumptions
summary(bg.adults.y1)
print(Anova(bg.adults.y1, type = 2), digits =7)

# Year 2
bg.adults.y2 <- glmmTMB(count ~ trt_y1*trt_y2 + (1|block), 
                        family = nbinom1(), 
                        data = balanus_adults %>% filter(period == "Winter 2021"))
bg.adults.y2.herb <- glmmTMB(count ~ trt_y1*trt_y2 + original_herb_trt + (1|block), 
                        family = nbinom1(), 
                        data = balanus_adults %>% filter(period == "Winter 2021"))

AIC(bg.adults.y2, bg.adults.y2.herb) # herbivore effect not significant

plot(simulateResiduals(bg.adults.y2)) # slight patterning indicating overdispersion, but doesn't violate any assumptions
summary(bg.adults.y2)
Anova(bg.adults.y2, type = 3)

emm2.bga <- emmeans(bg.adults.y2, specs = c("trt_y1","trt_y2"))
contrast(emm2.bga, method = "pairwise", adjust = "tukey")

# significance labels for adult Bg abundance plot
labels.2c <- as.data.frame(cbind(
  c("C","W","CC","CW","WC","WW"),
  c("Winter 2020","Winter 2020","Winter 2021","Winter 2021","Winter 2021","Winter 2021"),
  c(118,33,132,75,75,32),
  c("a","b","c","cd","cd","d"))
)
colnames(labels.2c) <- c("treatment","period", "count","label")
labels.2c <- labels.2c %>% mutate(count = as.numeric(count))


# Adult Chthamalus dalli abundance

# Year 1
cd.adults.y1 <- glmmTMB(count ~ treatment + (1|block), 
                        family = nbinom1(), 
                        data = chthamalus_adults %>% filter(period == "Winter 2020"))

# Model converges with original herbivore treatment included -- probably not enough observations
#cd.adults.y1.herb <- glmmTMB(count ~ treatment + original_herb_trt + (1|block), 
                        #family = nbinom1(), 
                        #data = chthamalus_adults %>% filter(period == "Winter 2020"))
#AIC(cd.adults.y1, cd.adults.y1.herb)

plot(simulateResiduals(cd.adults.y1)) # slight patterning indicating overdispersion, but doesn't violate any assumptions
summary(cd.adults.y1)
Anova(cd.adults.y1, type = 2)

# Year 2
cd.adults.y2 <- glmmTMB(count ~ trt_y1*trt_y2 + (1|block), 
                        family = nbinom1(), 
                        data = chthamalus_adults %>% filter(period == "Winter 2021"))
cd.adults.y2.herb <- glmmTMB(count ~ trt_y1*trt_y2 + original_herb_trt + (1|block), 
                        family = nbinom1(), 
                        data = chthamalus_adults %>% filter(period == "Winter 2021"))
AIC(cd.adults.y2, cd.adults.y2.herb) # not significant; exclude herbivore term from model

plot(simulateResiduals(cd.adults.y2)) # slight patterning indicating overdispersion, but doesn't violate any assumptions
summary(cd.adults.y2)
Anova(cd.adults.y2, type = 3)

emm2.cda <- emmeans(cd.adults.y2, specs = c("trt_y1", "trt_y2"))
contrast(emm2.cda, method = "pairwise", adjust = "tukey")

# labels for C. dalli adult abundance plot
labels.2d <- as.data.frame(cbind(
  c("C","W","CC","CW","WC","WW"),
  c("Winter 2020","Winter 2020","Winter 2021","Winter 2021","Winter 2021","Winter 2021"),
  c(15,15,110,110,110,110),
  c("a","a","c","c","c","c"))
)
colnames(labels.2d) <- c("treatment","period", "count","label")
labels.2d <- labels.2d %>% mutate(count = as.numeric(count))

###############
# Figure panels for adult barnacle abundance

# Plotting abundance of adult B. glandula before y2 summer and at the end of the experiment

Fig2C <- ggplot(balanus_adults, aes(x = treatment, y = count, col = treatment)) +
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
  geom_text(data = labels.2c, aes(label = label), size = 3, col = "black", fontface = "bold")
Fig2C

# Plotting abundance of adult C. dalli before y2 summer and at the end of the experiment
Fig2D <- ggplot(chthamalus_adults, aes(x = treatment, y = count, col = treatment)) +
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
  geom_text(data = labels.2d, aes(label = label), size = 3, col = "black", fontface = "bold")
Fig2D

# Assemble multiplanel figure
Fig2 <- ((Fig2A / Fig2B) | (Fig2C / Fig2D)) + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")

# save figure
ggsave(Fig2, filename = "./figures/Fig2.pdf", device = cairo_pdf, 
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
  facet_grid(rows = vars(species), scales = "free") +
  scale_x_date(date_labels = "%Y-%m", breaks=breaks)
FigS6

png("./figures/FigS6.png", res = 700, width = 8.5, height = 5, units = "in")
FigS6
dev.off()
