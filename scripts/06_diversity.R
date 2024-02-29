## Alpha & beta diversity
## Amelia Hesketh February 2024

# load packages
pkgs <- c("vegan", "tidyverse","glmmTMB","DHARMa", 
          "plotrix", "lubridate", "emmeans","car", "patchwork",
          "glmer", "biodiversityR")
lapply(pkgs, library, character.only = TRUE)
rm(pkgs)

################
# loading and formatting data for analyses and plotting

# load information about surveys
visits_info <- read_csv("./raw_data/design/SVSWS_survey_times.csv") %>% 
  group_by(survey_no) %>% summarize(date = max(date))

# join survey data and metadata
diversity <- read_csv("./clean_data/SVSWS_survey_clean.csv") %>% 
  left_join(visits_info) %>% 
  mutate(timesincestart = difftime(date, ymd("2019-04-12"), units = c("weeks")),
         count = as.integer(count)) %>% 
  group_by(survey_no) %>% 
  mutate(timesincestart = max(timesincestart),
         date = max(date)) %>% 
  mutate(abund = if_else(is.na(count), percent_cover, as.numeric(count))) %>% 
  select(species, abund, survey_no, tile_id, original_herb_trt) %>% 
  unique() %>% 
  pivot_wider(names_from = species, values_from = abund,values_fill = 0)

# calculate species richness and Shannon diversity
diversity.factors <- diversity %>% 
  select(survey_no, tile_id, original_herb_trt)

# note that, because Shannon diversity from count & cover data cannot be combined
# invertebrate and algal diversity metrics must be calculated separately
diversity.invert.abund <- diversity %>% 
  ungroup() %>% 
  select(any_of(invert_list$species))

diversity.algae.abund <- diversity %>% 
  ungroup() %>% 
  select(any_of(alga_list$species))

richness <- specnumber(diversity %>% ungroup() %>% select(-survey_no, -tile_id, -original_herb_trt))
shannon.algae <- diversity(diversity.algae.abund)
shannon.invert <- diversity(diversity.invert.abund)

# join diversity response metrics that were just calculated
# with tile-specific information (treatment, date, etc.)
diversity.responses <- diversity.factors %>% cbind(richness) %>% 
  cbind(shannon.algae) %>% cbind(shannon.invert) %>% 
  rename(richness = 4, shannon.algae = 5, shannon.invert = 6) %>% 
  left_join(tile_info) %>% left_join(survey_numbers) %>% 
  mutate(timesincestart = difftime(date, ymd("2019-04-12"), units = c("weeks")))

# Make sure treatments are coded correctly....
# in year 1, treatment is described by trt_y1 (C or W) but in year 2, treatment is
# described by the combination of trt_y1 and trt_y2.
diversity <- diversity.responses %>% 
  mutate(treatment = if_else(date < "2020-03-15", as.factor(trt_y1), as.factor(treatment)),
         block = as.factor(block),
         date = as.factor(date),
         herb = factor(if_else(is.na(original_herb_trt), "control","herb")),
         tile_id = as.factor(tile_id))

# prepare data for plotting
diversity.plot <- diversity %>% group_by(timesincestart, date, treatment) %>% 
  # summarize the mean and SE for richness and Shannon diversity
  summarize(mean.rich = mean(richness),se.rich = std.error(richness),
            mean.shann.i = mean(shannon.invert), se.shann.i = std.error(shannon.invert),
            mean.shann.a = mean(shannon.algae), se.shann.a = std.error(shannon.algae)) %>% 
  mutate(treatment = factor(treatment, levels = c("C","W","CC","CW","WC","WW")),
         date = ymd(date)) %>% 
  ungroup()

################
# plots of alpha diversity metrics over time for Appendix 1

# define palettes for plotting
pal.trt <- c("#014779", "#EE4B2B", "#014779", "#7985CB", "#9C0098", "#EE4B2B", "grey50")
pch.trt <- c(1,1,16,16,16,16,16)
lty.trt <- c("dotdash","dotdash","solid","solid","solid","solid","solid","solid")

# richness over time
FigA9a <- ggplot(diversity.plot, aes(y = mean.rich, x =date, 
                                       col = treatment, shape = treatment,
                                     lty = treatment)) +
  geom_point() +
  geom_line(lwd = 0.8) +
  scale_color_manual(values = pal.trt) +
  scale_linetype_manual(values = lty.trt) +
  scale_shape_manual(values = pch.trt) +
  geom_errorbar(aes(ymax = mean.rich + se.rich, ymin = mean.rich - se.rich, lty = NULL)) +
  labs(y = "Species richness", x = "Date", pch = "Treatment", col = "Treatment", lty = "Treatment") +
  geom_vline(aes(xintercept = as.Date("2020-04-07")), linetype = "dotted", col = "grey30", lwd = 0.8) +
  theme_classic() +
  theme(legend.key.width = unit(1,"cm"), plot.tag = element_text(face = "bold"))
FigA9a

# Shannon diversity of invertebrate community over time
FigA9b <- ggplot(diversity.plot, aes(y = mean.shann.i, x =date, col = treatment,
            theme_classic() +                            lty = treatment, pch=treatment)) +
  geom_point() +
  geom_line(lwd = 0.8) +
  scale_color_manual(values = pal.trt) +
  scale_linetype_manual(values = lty.trt) +
  scale_shape_manual(values = pch.trt) +
  theme_classic() + theme(legend.key.width = unit(1,"cm"), plot.tag = element_text(face = "bold")) +
  geom_errorbar(aes(ymax = mean.shann.i + se.shann.i, ymin = mean.shann.i - se.shann.i, lty = NULL)) +
  labs(y = "Invert. Shannon diversity", x = "Date", pch = "Treatment", col = "Treatment", lty = "Treatment") +
  geom_vline(aes(xintercept = as.Date("2020-04-07")), linetype = "dotted", col = "grey30", lwd = 0.8)
FigA9b

# Shannon diversity of algal community over time
FigA9c <- ggplot(diversity.plot, aes(y = mean.shann.a, x =date, 
                                        col = treatment, pch = treatment, lty = treatment)) +
  geom_point() +
  geom_line(lwd = 0.8) +
  scale_color_manual(values = pal.trt) +
  scale_linetype_manual(values = lty.trt) +
  scale_shape_manual(values = pch.trt) +
  theme_classic() + theme(legend.key.width = unit(1,"cm"),plot.tag = element_text(face = "bold")) +
  geom_errorbar(aes(ymax = mean.shann.a + se.shann.a, ymin = mean.shann.a - se.shann.a, lty = NULL)) +
  labs(y = "Algal Shannon diversity", x = "Date", pch = "Treatment", col = "Treatment", lty = "Treatment") +
  geom_vline(aes(xintercept = as.Date("2020-04-07")), linetype = "dotted", col = "grey30", lwd = 0.8)
FigA9c

# stitch them all together
FigA9 <- FigA9a / FigA9b / FigA9c +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") +
  plot_layout(guides = "collect") & theme(plot.tag = element_text(size = 14, face = "bold")) 
FigA9

# save plot
png("./figures/FigA9.png", res = 700, width = 9, height = 9, units = "in")
FigA9
dev.off()


## plotting instead diversity metrics at key timepoints

diversity_keytimes <- diversity %>% 
  filter(date %in% c("2019-10-20", "2020-03-15", "2020-09-14","2021-02-24")) %>% 
  mutate(period = if_else(date %in% c("2019-10-20", "2020-03-15"), "Year 1", "Year 2"),
         timept = if_else(date %in% c("2019-10-20", "2020-09-14"), "Post-summer", "Winter")) %>% 
  mutate(treatment = if_else(period == "Year 1", trt_y1, treatment),
         treatment = factor(treatment, levels = c("C","W","CC","CW","WC","WW")))

diversity_plots <- diversity_keytimes %>% 
  group_by(treatment, timept, period) %>% 
  summarize(mean.rich = mean(richness), se.rich = std.error(richness),
            mean.isd = mean(shannon.invert), se.isd = std.error(shannon.invert),
            mean.asd = mean(shannon.algae), se.asd = std.error(shannon.algae)) %>% ungroup()

# mean richness post-summer and wintertime in year 1 and year 2
Fig5a <- ggplot(diversity_plots, aes(x = timept, y = mean.rich, col = treatment, pch = treatment)) +
  geom_point(position = position_dodge(width = 1)) +
  geom_errorbar(aes(ymax = mean.rich + se.rich, ymin = mean.rich - se.rich),
                position = position_dodge(width = 1), width = 0.3) +
  facet_wrap(~period, scales ="free_x") +
  theme_classic() +
  scale_color_manual(values = pal.trt) +
  scale_shape_manual(values = pch.trt) +
  labs(x = "Time", pch = "Treatment", col = "Treatment", y = "Species richness") +
  theme(plot.tag = element_text(face = "bold"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

# mean Shannon diversity of inverts post-summer and wintertime in year 1 and year 2
Fig5b <- ggplot(diversity_plots, aes(x = timept, y = mean.isd, col = treatment, pch = treatment)) +
  geom_point(position = position_dodge(width = 1)) +
  geom_errorbar(aes(ymax = mean.isd + se.isd, ymin = mean.isd - se.isd),
                position = position_dodge(width = 1), width = 0.3) +
  facet_wrap(~period, scales ="free_x") +
  theme_classic() +
  scale_color_manual(values = pal.trt) +
  scale_shape_manual(values = pch.trt) +
  labs(x = "Time", pch = "Treatment", col = "Treatment", y = "Invertebrate Shannon diversity") +
  theme(plot.tag = element_text(face = "bold"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

################
# Modeling alpha diversity

# Species richness

# Year 1: post-summer
richness.stress.y1 <- glmmTMB(richness ~ trt_y1 + (1|block), 
                       family = poisson(),
                       data = diversity_keytimes %>% filter(date == "2019-10-20"))
richness.stress.y1.herb <- glmmTMB(richness ~ trt_y1 + herb+ (1|block), 
                         family = poisson(), 
                         data = diversity_y1 %>% filter(date == "2019-10-20"))

plot(simulateResiduals(richness.stress.y1))
compare_performance(richness.stress.y1, richness.stress.y1.herb) # better without herbivores 
summary(richness.stress.y1)
Anova(richness.stress.y1, type = 2)

# Year 1: winter

richness.recovery.y1 <- glmmTMB(richness ~ trt_y1 + (1|block), 
                              family = poisson(),
                              data = diversity_keytimes %>% filter(date == "2020-03-15"))

richness.recovery.y1.herb <- glmmTMB(richness ~ trt_y1 + herb + (1|block), 
                                family = poisson(),
                                data = diversity_y1 %>% filter(date == "2020-03-15"))

plot(simulateResiduals(richness.recovery.y1))
compare_performance(richness.recovery.y1, richness.recovery.y1.herb) # better without herbivores
summary(richness.recovery.y1)
Anova(richness.recovery.y1, type = 2)

# Year 2: post-summer
# including random effect causes singularity, so dropped from model here 

richness.stress.y2 <- glmmTMB(richness ~ trt_y1*trt_y2 + (1|block),
                       family = poisson(),
                       data = diversity_keytimes %>% filter(date == "2020-09-14"))
richness.stress.y2.herb <- glmmTMB(richness ~ trt_y1*trt_y2 + herb + (1|block), 
                       family = poisson(),
                       data = diversity_y2 %>% filter(date == "2020-09-14"))

compare_performance(richness.stress.y2, richness.stress.y2.herb) # better without herbivores
plot(simulateResiduals(richness.stress.y2)) # model is under-dispersed, which will yield more conservative p values
# no other distribution is appropriate, so we will keep this as-is despite violations
summary(richness.stress.y2)
Anova(richness.stress.y2, type = 3)
Anova(richness.stress.y2, type = 2)

# Year 2: winter

richness.recovery.y2 <- glmmTMB(richness ~ trt_y1*trt_y2 + (1|block), 
                              family = poisson(),
                              data = diversity_keytimes %>% filter(date == "2021-02-24"))
richness.recovery.y2.herb <- glmmTMB(richness ~ trt_y1*trt_y2 + herb +(1|block), 
                                   family = poisson(),
                                   data = diversity_y2 %>% filter(date == "2021-02-24"))

compare_performance(richness.recovery.y2, richness.recovery.y2.herb) # better without herbivores
plot(simulateResiduals(richness.recovery.y2))
summary(richness.recovery.y2)
Anova(richness.recovery.y2, type = 3)
Anova(richness.recovery.y2, type = 2)

# Invertebrate Shannon diversity

# Year 1: post-summer
isd.stress.y1 <- glmmTMB(shannon.invert ~ trt_y1 + (1|block), 
                         family = tweedie(), dispformula = ~trt_y1,
                            data = diversity_keytimes %>% filter(date == "2019-10-20"))
isd.stress.y1.herb <- glmmTMB(shannon.invert ~ trt_y1 + herb+ (1|block), 
                              dispformula = ~trt_y1,
                              family = tweedie(), data = diversity_y1 %>% filter(date == "2019-10-20"))

plot(simulateResiduals(isd.stress.y1))
compare_performance(isd.stress.y1, isd.stress.y1.herb) # better without herbivores 
summary(isd.stress.y1)
Anova(isd.stress.y1, type = 2)

# Year 1: winter

isd.recovery.y1 <- glmmTMB(shannon.invert ~ trt_y1 + (1|block), 
                              family = tweedie(),
                              data = diversity_keytimes %>% filter(date == "2020-03-15"))

isd.recovery.y1.herb <- glmmTMB(shannon.invert ~ trt_y1 + herb + (1|block), 
                                   family = tweedie(),
                                   data = diversity_y1 %>% filter(date == "2020-03-15"))

plot(simulateResiduals(isd.recovery.y1))
compare_performance(isd.recovery.y1, isd.recovery.y1.herb) # better without herbivores
summary(isd.recovery.y1)
Anova(isd.recovery.y1, type = 2)

# Year 2: post-summer
# including random effect causes singularity, so dropped from model here 

isd.stress.y2 <- glmmTMB(shannon.invert ~ trt_y1*trt_y2 + (1|block),
                          data = diversity_keytimes %>% filter(date == "2020-09-14"))
isd.stress.y2.herb <- glmmTMB(shannon.invert ~ trt_y1*trt_y2 + herb + (1|block), 
                               data = diversity_y2 %>% filter(date == "2020-09-14"))

compare_performance(isd.stress.y2, isd.stress.y2.herb) # better without herbivores
plot(simulateResiduals(isd.stress.y2))
summary(isd.stress.y2)
Anova(isd.stress.y2, type = 3)
Anova(isd.stress.y2, type = 2)

# Year 2: winter

isd.recovery.y2 <- glmmTMB(shannon.invert ~ trt_y1*trt_y2 + (1|block), 
                              data = diversity_keytimes %>% filter(date == "2021-02-24"))
isd.recovery.y2.herb <- glmmTMB(shannon.invert ~ trt_y1*trt_y2 + herb +(1|block), 
                                   data = diversity_y2 %>% filter(date == "2021-02-24"))

compare_performance(isd.recovery.y2, isd.recovery.y2.herb) # better without herbivores
plot(simulateResiduals(isd.recovery.y2))
summary(isd.recovery.y2)
Anova(isd.recovery.y2, type = 3)
Anova(isd.recovery.y2, type = 2)

##################
# beta analysis of invertebrate community living within foundation spp

set.seed(26)

# load community data
epifauna_sept20 <- read_csv("./raw_data/epifauna/SVSWS_20200914_epifauna.csv")
epifauna_feb21 <- read_csv("./raw_data/epifauna/SVSWS_20210224_epifauna.csv")


# read in relevant metadata
block_design <- read_csv("./raw_data/design/SVSWS_tilesetup.csv")
codes <- read_csv("./raw_data/epifauna/taxonomic_codes.csv")

# all tiles were collected late in the experiment, so block & number of tiles corresponds
# to new block and new number
treatments_epifauna <- block_design %>% select(new_block, new_no, treatment) %>% 
  rename(block = new_block, number = new_no)

# since we ended up lumping some of the taxa, ensure all like taxa (e.g. amphipod adults & larvae)
# are united when summing abundance
epifauna_summer <- epifauna_sept20 %>%
  left_join(treatments_epifauna) %>% 
  left_join(codes) %>% select(-taxon, -taxon_repaired, -sp_code) %>% 
  group_by(block, number, treatment, unified_code) %>% 
  summarize(abund = sum(abund)) %>% ungroup() %>% 
  # pivot for dbRDA (wide format required)
  pivot_wider(id_cols = c(block, number, treatment), names_from = unified_code, values_from = abund,
              values_fill = list(abund =0))

# extract factors for analyses
epifauna_factors.summer <- epifauna_summer %>% 
  select(block, number, treatment) 

# convert wide-format data to abundance matrix
epifauna_matrix.summer <- epifauna_summer %>% 
  select(-block, -number, - treatment) %>% 
  as.matrix()


epi_summer <- dbrda(epifauna_matrix.summer ~ treatment, epifauna_factors.summer, dist = "bray")
summary(epi_summer)
screeplot(epi_summer)

# Create biplot

# extract coordinates of tiles
data.scores.summer <- as.data.frame(scores(epi_summer)$sites)
# associate these with experimental treatments
data.scores.summer$treatment <- epifauna_factors.summer$treatment

# calculate the average position of tiles in ordination space for plotting centroids
data.scores.summer <- data.scores.summer %>% group_by(treatment) %>% 
  mutate(mean.x = mean(dbRDA1), mean.y = mean(dbRDA2)) %>% ungroup()

## adding ellipses to plot

veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

df_ell_summer <- data.frame()
for(g in levels(as.factor(data.scores.summer$treatment))){
  df_ell_summer <- rbind(df_ell_summer, cbind(as.data.frame(with(data.scores.summer[data.scores.summer$treatment==g,],
                                                                 veganCovEllipse(cov.wt(cbind(dbRDA1,dbRDA2),
                                                                                        wt=rep(1/length(dbRDA1),
                                                                                               length(dbRDA1)))$cov,center=c(mean(dbRDA1),mean(dbRDA2)))))
                                              ,group=g))
}

df_ell_summer$treatment <- df_ell_summer$group

pal.trt.y2 <- c("#014779", "#7985CB", "#9C0098", "#EE4B2B")

# ordination in ggplot2
Fig5c <- ggplot() + 
  geom_path(data=df_ell_summer, aes(x=dbRDA1,y=dbRDA2,colour=treatment)) +
  geom_point(data=data.scores.summer, aes(x=dbRDA1,y=dbRDA2,colour=treatment),
             size=3, alpha = 0.8) + # add the point markers
  scale_color_manual(values = pal.trt.y2) +
  labs(color = "Treatment", x = "dbRDA1", y = "dbRDA2") +
  plot_theme +
  theme(plot.tag = element_text(face = "bold"),
        legend.position = "none")
Fig5c
# PerMANOVA
# need to constrain permutations within experimental blocks.
# also increase number of permutations to make up for constrained permutation.
h <- how(blocks = epifauna_factors.summer$block, nperm = 9999) #constrain permutation
anova.cca(epi_summer, permutations = h, by = "term")
multiconstrained(method = "dbrda", epifauna_matrix.summer ~ treatment, data = epifauna_factors.summer,
                 distance = "bray")

# Permdisp
dist.epi.summer <- vegdist(epifauna_matrix.summer, method = "bray")

disp.epi.summer <- betadisper(dist.epi.summer, type = "centroid",  bias.adjust = T,
                              group = epifauna_factors.summer$treatment)
anova(disp.epi.summer)
TukeyHSD(disp.epi.summer)


## and second data point for spring samples

epifauna_winter <- epifauna_feb21 %>%
  left_join(treatments_epifauna) %>% 
  left_join(codes) %>% select(-taxon, -taxon_repaired, -sp_code) %>% 
  group_by(block, number, treatment, unified_code) %>% 
  summarize(abund = sum(abund)) %>% ungroup() %>% 
  filter(is.na(unified_code) == F) %>% 
  pivot_wider(id_cols = c(block, number, treatment), names_from = unified_code, values_from = abund,
              values_fill = list(abund =0))

epifauna_factors.winter <- epifauna_winter %>% 
  select(block, number, treatment) %>% 
  mutate(trty1 = substring(treatment, 1, 1),
         trty2 = substring(treatment, 2,2))

epifauna_matrix.winter <- epifauna_winter %>% 
  select(-block, -number, - treatment) %>% 
  as.matrix() 

epi_winter <- dbrda(epifauna_matrix.winter ~ treatment, epifauna_factors.winter, dist = "bray")
summary(epi_winter)
screeplot(epi_winter)

# extract coordinates of tiles
data.scores.winter <- as.data.frame(scores(epi_winter)$sites)
# associate these with experimental treatments
data.scores.winter$treatment <- epifauna_factors.winter$treatment

# calculate the average position of tiles in ordination space for plotting centroids
data.scores.winter <- data.scores.winter %>% group_by(treatment) %>% 
  mutate(mean.x = mean(dbRDA1), mean.y = mean(dbRDA2)) %>% ungroup()

## adding ellipses to plot

veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

df_ell.winter <- data.frame()
for(g in levels(as.factor(data.scores.winter$treatment))){
  df_ell.winter <- rbind(df_ell.winter, cbind(as.data.frame(with(data.scores.winter[data.scores.winter$treatment==g,],
                                                                 veganCovEllipse(cov.wt(cbind(dbRDA1,dbRDA2),
                                                                                        wt=rep(1/length(dbRDA1),
                                                                                               length(dbRDA1)))$cov,center=c(mean(dbRDA1),mean(dbRDA2)))))
                                              ,group=g))
}

df_ell.winter$treatment <- df_ell.winter$group

# ordination in ggplot2
Fig5d <- ggplot() + 
  geom_path(data=df_ell.winter, aes(x=dbRDA1,y=dbRDA2,colour=treatment)) +
  geom_point(data=data.scores.winter, aes(x=dbRDA1,y=dbRDA2,colour=treatment),
             size=3, alpha = 0.8) + # add the point markers
  scale_color_manual(values = pal.trt.y2) +
  scale_shape_manual(values = pch.trt.y2) +
  labs(color = "Treatment", x = "dbRDA1", y = "dbRDA2") +
  plot_theme  +
  theme(plot.tag = element_text(face = "bold"),
        legend.position = "none")
Fig5d
# PerMANOVA
# need to constrain permutations within experimental blocks.
# also increase number of permutations to make up for constrained permutation.
h <- how(blocks = epifauna_factors.winter$block, nperm = 9999) #constrain permutation
anova.cca(epi_winter, permutations = h, by = "term")
multiconstrained(method = "dbrda", epifauna_matrix.winter ~ treatment, data = epifauna_factors.winter,
                 distance = "bray")


# Permdisp
dist.epi.winter <- vegdist(epifauna_matrix.winter, method = "bray")

disp.epi.winter <- betadisper(dist.epi.winter, type = "centroid",  bias.adjust = T,
                              group = epifauna_factors.winter$treatment)
anova(disp.epi.winter)


Fig5 <- ((Fig5a / Fig5b) | (Fig5c / Fig5d)) +
  plot_layout(guides = "collect", widths = c(0.6, 0.4)) +
  plot_annotation(tag_levels = "a", tag_prefix = "(",
                  tag_suffix = ")") 
Fig5

png(filename = "./figures/Fig5.png", width = 10, height = 7, 
    res = 700, units = "in")
Fig5
dev.off()



## Algal Shannon diversity (Appendix 1 content)

# mean Shannon diversity of algae post-summer and wintertime in year 1 and year 2
FigA10 <- ggplot(diversity_plots, aes(x = timept, y = mean.asd, col = treatment, pch = treatment)) +
  geom_point(position = position_dodge(width = 1)) +
  geom_errorbar(aes(ymax = mean.asd + se.asd, ymin = mean.asd - se.asd),
                position = position_dodge(width = 1), width = 0.5) +
  facet_wrap(~period, scales ="free_x") +
  theme_classic() +
  scale_color_manual(values = pal.trt) +
  scale_shape_manual(values = pch.trt) +
  labs(x = "Time", pch = "Treatment", col = "Treatment", y = "Algal Shannon diversity") +
  theme(plot.tag = element_text(face = "bold"))

png("./figures/FigA10.png", res = 700, width = 6, height = 5,
    units = "in")
FigA10
dev.off()

# try combining with algal cover figure to see how this looks

Fig4.1 <- Fig4/FigA10 + plot_layout(guides = "collect") + plot_annotation(tag_levels = "a",
                                                                          tag_prefix = "(",
                                                                          tag_suffix = ")")
png("./figures/Fig4.1.png", res = 700, width = 5.5, height = 6,
    units = "in")
Fig4.1
dev.off()

# Year 1: post-summer
asd.stress.y1 <- glmmTMB(shannon.algae ~ trt_y1 + (1|block), 
                         family = tweedie(), dispformula = ~trt_y1,
                         data = diversity_keytimes %>% filter(date == "2019-10-20"))
asd.stress.y1.herb <- glmmTMB(shannon.algae ~ trt_y1 + herb+ (1|block), 
                              dispformula = ~trt_y1,
                              family = tweedie(), data = diversity_y1 %>% filter(date == "2019-10-20"))

plot(simulateResiduals(asd.stress.y1))
compare_performance(asd.stress.y1, asd.stress.y1.herb) # better without herbivores 
summary(asd.stress.y1)
Anova(asd.stress.y1, type = 2)

# Year 1: winter
asd.recovery.y1 <- glmmTMB(shannon.algae ~ trt_y1 + (1|block), 
                           family = tweedie(),
                           data = diversity_keytimes %>% filter(date == "2020-03-15"))

asd.recovery.y1.herb <- glmmTMB(shannon.algae ~ trt_y1 + herb + (1|block), 
                                family = tweedie(),
                                data = diversity_y1 %>% filter(date == "2020-03-15"))

plot(simulateResiduals(asd.recovery.y1))
compare_performance(asd.recovery.y1, asd.recovery.y1.herb) # better without herbivores
summary(asd.recovery.y1)
Anova(asd.recovery.y1, type = 2)

# Year 2: post-summer has too many zeroes, can't model data

# Year 2: winter

asd.recovery.y2 <- glmmTMB(shannon.algae ~ trt_y1*trt_y2 + (1|block), 
                           family = tweedie(),
                           data = diversity_keytimes %>% filter(date == "2021-02-24"))
asd.recovery.y2.herb <- glmmTMB(shannon.algae ~ trt_y1*trt_y2 + herb +(1|block), 
                                data = diversity_y2 %>% filter(date == "2021-02-24"))

compare_performance(asd.recovery.y2, asd.recovery.y2.herb) # better without herbivores
plot(simulateResiduals(asd.recovery.y2))
summary(asd.recovery.y2)
Anova(asd.recovery.y2, type = 3)
Anova(asd.recovery.y2, type = 2)


