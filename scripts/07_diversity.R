## Alpha & beta diversity
## Amelia Hesketh February 2024

# load packages
pkgs <- c("vegan", "tidyverse","glmmTMB","DHARMa", 
          "plotrix", "lubridate", "emmeans","car", "patchwork",
          "lme4", "BiodiversityR")
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
  select(species, abund, survey_no, tile_id, second_herb_trt) %>% 
  unique() %>% 
  pivot_wider(names_from = species, values_from = abund,values_fill = 0)

# calculate species richness and Shannon diversity
diversity.factors <- diversity %>% 
  select(survey_no, tile_id, second_herb_trt)

# note that, because Shannon diversity from count & cover data cannot be combined
# invertebrate and algal diversity metrics must be calculated separately

spp_list <- read_csv("./clean_data/SVSWS_species_list.csv")
alga_list <- spp_list %>% filter(type == "alga")
invert_list <- spp_list %>% filter(type == "invert")

diversity.invert.abund <- diversity %>% 
  ungroup() %>% 
  select(any_of(invert_list$species))

diversity.algae.abund <- diversity %>% 
  ungroup() %>% 
  select(any_of(alga_list$species))

richness <- specnumber(diversity %>% ungroup() %>% select(-survey_no, -tile_id, -second_herb_trt))
shannon.algae <- diversity(diversity.algae.abund)
shannon.invert <- diversity(diversity.invert.abund)

# join diversity response metrics that were just calculated
# with tile-specific information (treatment, date, etc.)

tile_info <- read_csv("./clean_data/SVSWS_tile_treatments.csv")

diversity.responses <- diversity.factors %>% cbind(richness) %>% 
  cbind(shannon.algae) %>% cbind(shannon.invert) %>% 
  rename(richness = 4, shannon.algae = 5, shannon.invert = 6) %>% 
  left_join(tile_info) %>% left_join(survey_numbers) %>% 
  mutate(timesincestart = difftime(date, ymd("2019-04-12"), units = c("weeks")))

# Make sure treatments are coded correctly....
# in year 1, treatment is described by trt_y1 (C or W) but in year 2, treatment is
# described by the combination of trt_y1 and trt_y2.
diversity_y1 <- diversity.responses %>% 
  filter(date <= "2020-03-15") %>% 
  mutate(treatment = as.factor(trt_y1),
         block = as.factor(block),
         date = as.factor(date),
         tile_id = as.factor(tile_id))

diversity_y2 <- diversity.responses %>% 
  filter(date >= "2020-03-15") %>% 
  mutate(treatment = as.factor(treatment),
         block = as.factor(block),
         date = as.factor(date),
         tile_id = as.factor(tile_id))

diversity <- diversity_y1 %>% full_join(diversity_y2) %>% 
  mutate(herb_trt = if_else(is.na(second_herb_trt), "control","grazer"))

# prepare data for plotting
diversity.plot <- diversity %>% group_by(timesincestart, date, treatment) %>% 
  # summarize the mean and SE for richness and Shannon diversity
  summarize(mean.rich = mean(richness),se.rich = std.error(richness),
            mean.shann.i = mean(shannon.invert), se.shann.i = std.error(shannon.invert),
            mean.shann.a = mean(shannon.algae), se.shann.a = std.error(shannon.algae)) %>% 
  mutate(treatment = factor(treatment, levels = c("C","W","CC","CW","WC","WW")),
         date = ymd(date)) %>% 
  ungroup() %>% 
  filter((treatment %in% c("C","W") & date > "2020-04-15") == F)

################
# Modeling alpha diversity

# Species richness

# Year 1 post summer
richness.y1.ps.herb <- glmmTMB(richness ~ trt_y1+ herb_trt + (1|block), 
                         family = poisson(), 
                         data = diversity %>% filter(date == "2019-10-20"))
plot(simulateResiduals(richness.y1.ps.herb))
summary(richness.y1.ps.herb)
Anova(richness.y1.ps.herb, type = 2)
# grazer treatment not significant; drop this from the model

richness.y1.ps <- glmmTMB(richness ~ trt_y1 + (1|block),
                          family = poisson(), 
                          data = diversity %>% filter(date == "2019-10-20"))
plot(simulateResiduals(richness.y1.ps))
plot(residuals(richness.y1.ps)) # slightly underdispersed residuals, but overall not too bad
summary(richness.y1.ps)
Anova(richness.y1.ps, type = 2)


emm1.rich.ps <- emmeans(richness.y1.ps, specs = c("trt_y1"))
contrast(emm1.rich.ps, method = "pairwise", adjust = "tukey")

# Year 1 winter
richness.y1.w.herb <- glmmTMB(log(richness+1) ~ trt_y1+ herb_trt + (1|block), 
                               data = diversity %>% filter(date == "2020-03-15"))

plot(simulateResiduals(richness.y1.w.herb)) # just over significance for distribution test; should be OK for anova
summary(richness.y1.w.herb)
Anova(richness.y1.w.herb, type = 2)

# grazer treatment not significant; drop this from the model

richness.y1.w <- glmmTMB(richness ~ trt_y1 + (1|block), 
                          family = poisson(),
                          data = diversity %>% filter(date == "2020-03-15"))
plot(simulateResiduals(richness.y1.w))
plot(residuals(richness.y1.w)) # KS test failed (p = 0.04); this is fairly marginal, 
# and the distribution is overall appropriate to the data (count)
summary(richness.y1.w)
Anova(richness.y1.w, type = 2)


emm1.rich.w <- emmeans(richness.y1.w, specs = c("trt_y1"))
contrast(emm1.rich.w, method = "pairwise", adjust = "tukey")

# labels for plot of richness (emmeans results)
labels.6a <- as.data.frame(cbind(
  c(rep("Post-summer", times = 2), rep("Winter", times = 2)),
  c("C","W", "C","W"),
  c(7.8,5.9,9,5),
  c("a","b","a","b"))
)
colnames(labels.6a) <- c("season", "treatment","richness","label")
labels.6a <- labels.6a %>% mutate(richness = as.numeric(richness))


# Year 2 post-summer

richness.y2.ps <- glmmTMB(log(richness+1) ~ trt_y1*trt_y2 + (1|block),
                       data = diversity %>% filter(date == "2020-09-14"))

plot(simulateResiduals(richness.y2.ps)) 
summary(richness.y2.ps)
Anova(richness.y2.ps, type = 3, contrasts=list(trt_y1 = "contr.sum", trt_y2 = "contr.sum"))

emm2.rich.ps <- emmeans(richness.y2.ps, specs = c("trt_y1","trt_y2"))
contrast(emm2.rich.ps, method = "pairwise", adjust = "tukey")

richness.y2.w <- glmmTMB(richness ~ trt_y1*trt_y2 + (1|block),
                       family = poisson(),
                       data = diversity %>% filter(date == "2021-02-24"))

plot(simulateResiduals(richness.y2.w))
summary(richness.y2.w)
Anova(richness.y2.w, type = 3, contrasts=list(trt_y1 = "contr.sum", trt_y2 = "contr.sum"))

emm2.rich.w <- emmeans(richness.y2.w, specs = c("trt_y1","trt_y2"))
contrast(emm2.rich.w, method = "pairwise", adjust = "tukey")

# create dataframe with significance labels for eventual plot
labels.6b <- as.data.frame(cbind(
  c(rep("Post-summer", times = 4), rep("Winter", times = 4)),
  c(rep(c("CC","CW","WC","WW"), times = 2)),
  c(7.8,7.8,7.8,7.8,9.9,7.8,9.9,7.8),
  c("a","ab","a","b", "c","cd","c","d"))
)
colnames(labels.6b) <- c("season","treatment","richness","label")
labels.6b <- labels.6b %>% mutate(richness = as.numeric(richness))


################
# plots of alpha diversity metrics at key times (post summer and in winter)

# define palettes for plotting
pal.trt <- c("#014779", "#EE4B2B", "#014779", "#7985CB", "#9C0098", "#EE4B2B")
pch.trt <- c(1,1,16,16,16,16)

# Create summary dataframe
diversity_keytimes <- diversity %>% 
  filter(date %in% c("2019-10-20", "2020-03-15", "2020-09-14","2021-02-24")) %>% 
  mutate(period = if_else(date %in% c("2019-10-20", "2020-03-15"), "Year 1", "Year 2"),
         season = if_else(date %in% c("2019-10-20", "2020-09-14"), "Post-summer", "Winter")) %>% 
  mutate(treatment = if_else(period == "Year 1", trt_y1, treatment),
         treatment = factor(treatment, levels = c("C","W","CC","CW","WC","WW")),
         date = factor(date))

# Plot as single point ± one standard error, so summarize mean and standard error
diversity_plots <- diversity_keytimes %>% 
  group_by(treatment, season, period) %>% 
  summarize(mean.rich = mean(richness), se.rich = std.error(richness),
            mean.isd = mean(shannon.invert), se.isd = std.error(shannon.invert),
            mean.asd = mean(shannon.algae), se.asd = std.error(shannon.algae)) %>% ungroup()

# Extract year 1 and year 2 data for separate plots
diversity_plot_y1_summary <- diversity_plots %>% filter(period == "Year 1") 
diversity_plot_y2_summary <- diversity_plots %>% filter(period == "Year 2") %>% mutate(factor(treatment))

diversity_plot_y1 <- diversity_keytimes %>% filter(period == "Year 1") 
diversity_plot_y2 <- diversity_keytimes %>% filter(period == "Year 2")

# mean richness post-summer and wintertime in year 1
Fig6a <- ggplot(diversity_plot_y1, aes(x = treatment, y = richness, 
                                       col = treatment, pch = season)) +
  geom_point(data = diversity_plot_y1_summary, aes(x = treatment, y = mean.rich),
             size = 2, stroke = 1, alpha = 0.9) +
  geom_errorbar(data = diversity_plot_y1_summary, 
                aes(y = mean.rich, 
                    ymax = mean.rich + se.rich, 
                    ymin = mean.rich - se.rich),
                width = 0.25, color = "grey20") +
  geom_jitter(alpha = 0.1, height = 0, width = 0.2) +
  facet_wrap(~season) +
  theme_classic() +
  scale_color_manual(values = pal.trt, guide = "none") +
  scale_shape_manual(values = c(1,2), guide = "none") +
  labs(x = "Treatment", pch = "Season", col = "Treatment", y = "Species richness") +
  theme(plot.tag = element_text(face = "bold", size = 10),
        strip.text = element_text(size = 9.5),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8)) +
  scale_y_continuous(limits = c(0,10), breaks = c(0,2,4,6,8,10)) +
  geom_text(data = labels.6a, aes(x = treatment, y = richness, label = label),
           col = "black", size = 3, fontface = "bold")
Fig6a

# mean richness post-summer and wintertime in year 2

Fig6b <- ggplot(diversity_plot_y2_summary, aes(x = treatment, y = mean.rich, 
                                       col = treatment, pch = season)) +
  geom_jitter(data = diversity_plot_y2, aes(y = richness),
              alpha = 0.2, height = 0, width = 0.2, show.legend = F) +
  geom_point(size = 3, alpha = 0.8, show.legend = T) +
  geom_errorbar(aes(ymax = mean.rich + se.rich, 
                    ymin = mean.rich - se.rich),
                width = 0.25, color = "grey20") +
  theme_classic() +
  facet_wrap(~season) +
  labs(x = "Treatment", pch = "Season", col = "Treatment", y = "Species richness") +
  geom_text(data = labels.6b, aes(x = treatment, y = richness, label = label),
            col = "black", size = 3, fontface = "bold") +
  scale_color_manual(values = pal.trt, drop = FALSE) +
  scale_y_continuous(limits = c(0,10), breaks = c(0,2,4,6,8,10)) +
  guides(color = guide_legend(override.aes = list(shape = pch.trt, size = 3))) +
  theme(plot.tag = element_text(face = "bold", size = 10),
        strip.text = element_text(size = 9.5),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8),
        legend.box = "horizontal")
Fig6b


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
              values_fill = list(abund = 0))

# extract factors for analyses
epifauna_factors.summer <- epifauna_summer %>% 
  select(block, number, treatment) %>% 
  mutate(trty1 = substring(treatment, 1, 1),
         trty2 = substring(treatment, 2,2))

# convert wide-format data to abundance matrix
epifauna_matrix.summer <- epifauna_summer %>% 
  select(-block, -number, - treatment) %>% 
  as.matrix() 


epi_summer <- dbrda(epifauna_matrix.summer ~ trty1*trty2, epifauna_factors.summer, dist = "bray")
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


## Ellipses for each treatment
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

df_ell_summer <- df_ell_summer %>% mutate(treatment = factor(case_when(treatment == "CC" ~ "C | CC",
                                                                treatment == "CW" ~ "CW",
                                                                treatment == "WC" ~ "WC",
                                                                treatment == "WW" ~ "W | WW"),
                                                             levels = c("C | CC", "CW","WC","W | WW")))

data.scores.summer <- data.scores.summer %>% mutate(treatment = factor(case_when(treatment == "CC" ~ "C | CC",
                                                                          treatment == "CW" ~ "CW",
                                                                          treatment == "WC" ~ "WC",
                                                                          treatment == "WW" ~ "W | WW"),
                                                    levels = c("C | CC", "CW","WC","W | WW")))  


# ordination in ggplot2
Fig6c <- ggplot() + 
  # add ellipses
  geom_path(data=df_ell_summer, aes(x=dbRDA1,y=dbRDA2,colour=treatment),
            show.legend = F) +
  # add points for each tile
  geom_point(data=data.scores.summer, aes(x=dbRDA1,y=dbRDA2,colour=treatment),
             size=2, alpha = 0.8, show.legend = F) + # add the point markers
  scale_color_manual(values = pal.trt.y2) +
  theme_classic()+
  labs(color = "Treatment", x = "dbRDA1", y = "dbRDA2") +
  theme(plot.tag = element_text(face = "bold", size = 10),
        strip.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8), 
        legend.position = "none") +
  coord_equal()+
  # add labels showing results of pairwise comparisons
  annotate(geom = "text", x = -0.9, y = 1, label = "a", col ="#014779", size = 3, fontface = "bold") +
  annotate(geom = "text", x = -0.5, y = -0.53, label = "ab", col = "#7985CB", size = 3, fontface = "bold") +
  annotate(geom = "text", x = -0.08, y = -0.7, label = "ab", col = "#9C0098", size = 3, fontface = "bold") +
  annotate(geom = "text", x = 1.25, y = 0.8, label = "b", col = "#EE4B2B", size = 3, fontface = "bold") 
Fig6c

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


## Repeat for samples taken at the end of the experiment (February 2021)

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

epi_winter <- dbrda(epifauna_matrix.winter ~ trty1*trty2, epifauna_factors.winter, dist = "bray")
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
Fig6d <- ggplot() + 
  geom_path(data=df_ell.winter, aes(x=dbRDA1,y=dbRDA2,colour=treatment),show.legend = F) +
  geom_point(data=data.scores.winter, aes(x=dbRDA1,y=dbRDA2,colour=treatment),
             size=2, alpha = 0.8, shape = 17, show.legend = F) + # add the point markers
  scale_color_manual(values = pal.trt.y2) +
  labs(color = "Treatment", x = "dbRDA1", y = "dbRDA2") +
  theme_classic()  + 
  theme(plot.tag = element_text(face = "bold", size = 10),
        legend.position = "none") +
  coord_equal()+
  annotate(geom = "text", x = 1, y = -0.7, label = "c", col ="#014779", size = 3, fontface = "bold") +
  annotate(geom = "text", x = -0.6, y = 0.8, label = "ab", col = "#7985CB", size = 3, fontface = "bold") +
  annotate(geom = "text", x = 0.52, y = 0.5, label = "b", col = "#9C0098", size = 3, fontface = "bold") +
  annotate(geom = "text", x = -0.9, y = -0.6, label = "a", col = "#EE4B2B", size = 3, fontface = "bold") 
Fig6d

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

# define layout
layout <- 
"AAAAAABBBBBBBBBB
 CCCCCCCDDDDDDDDD
"

# assemble multipanel figure
Fig6 <- Fig6a + Fig6b + Fig6c + Fig6d +
  plot_layout(guides = "collect", 
              design = layout, heights = c(0.4,0.6)) +
  plot_annotation(tag_levels = "a", tag_prefix = "(",
                  tag_suffix = ")") &
  theme(legend.position = "bottom",
        legend.box = "horizontal",
        legend.box.just = "left") 
Fig6


# save figure
ggsave(Fig6, filename = "./figures/Fig6.pdf", device = cairo_pdf, 
       width = 18, height = 18, units = "cm")

################ 
# Appendix content: Figures & analyses

# Plot of species richness of treatments through time
FigS10a <- ggplot(diversity.plot, aes(y = mean.rich, x =date, col = treatment,lty = treatment, pch=treatment)) +
  geom_vline(aes(xintercept = ymd("2020-04-03")), col = "grey60", lwd = 0.8) +
  geom_vline(aes(xintercept = ymd("2019-08-27")), col = "grey70", lty = "dashed", lwd = 0.8) +
  theme_classic() +                            
  geom_point() +
  geom_line(lwd = 0.8) +
  scale_color_manual(values = pal.trt) +
  scale_linetype_manual(values = lty.trt) +
  scale_shape_manual(values = pch.trt) +
  theme_classic() + theme(legend.key.width = unit(1,"cm"), plot.tag = element_text(face = "bold")) +
  geom_errorbar(aes(ymax = mean.rich + se.rich, ymin = mean.rich - se.rich, lty = NULL)) +
  labs(y = "Species richness", x = "Date", pch = "Treatment", col = "Treatment", lty = "Treatment") +
  scale_x_date(date_labels = "%Y-%m", breaks=breaks)
FigS10a

# Shannon diversity of invertebrate community over time
FigS10b <- ggplot(diversity.plot, aes(y = mean.shann.i, x =date, col = treatment,lty = treatment, pch=treatment)) +
  geom_vline(aes(xintercept = ymd("2020-04-03")), col = "grey60", lwd = 0.8) +
  geom_vline(aes(xintercept = ymd("2019-08-27")), col = "grey70", lty = "dashed", lwd = 0.8) +
  theme_classic() +                            
  geom_point() +
  geom_line(lwd = 0.8) +
  scale_color_manual(values = pal.trt) +
  scale_linetype_manual(values = lty.trt) +
  scale_shape_manual(values = pch.trt) +
  theme_classic() + theme(legend.key.width = unit(1,"cm"), plot.tag = element_text(face = "bold")) +
  geom_errorbar(aes(ymax = mean.shann.i + se.shann.i, ymin = mean.shann.i - se.shann.i, lty = NULL)) +
  labs(y = "Invert. Shannon diversity", x = "Date", pch = "Treatment", col = "Treatment", lty = "Treatment") +
  scale_x_date(date_labels = "%Y-%m", breaks=breaks)
FigS10b

# Shannon diversity of algal community over time
FigS10c <- ggplot(diversity.plot, aes(y = mean.shann.a, x =date, 
                                     col = treatment, pch = treatment, lty = treatment)) +
  geom_vline(aes(xintercept = ymd("2020-04-03")), col = "grey60", lwd = 0.8) +
  geom_vline(aes(xintercept = ymd("2019-08-27")), col = "grey70", lty = "dashed", lwd = 0.8) +
  geom_point() +
  geom_line(lwd = 0.8) +
  scale_color_manual(values = pal.trt) +
  scale_linetype_manual(values = lty.trt) +
  scale_shape_manual(values = pch.trt) +
  theme_classic() + theme(legend.key.width = unit(1,"cm"),plot.tag = element_text(face = "bold")) +
  geom_errorbar(aes(ymax = mean.shann.a + se.shann.a, ymin = mean.shann.a - se.shann.a, lty = NULL)) +
  labs(y = "Algal Shannon diversity", x = "Date", pch = "Treatment", col = "Treatment", lty = "Treatment") +
  scale_x_date(date_labels = "%Y-%m", breaks=breaks)
FigS10c

# stitch them all together
FigS10 <- FigS10a / FigS10b / FigS10c +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") +
  plot_layout(guides = "collect") & theme(plot.tag = element_text(size = 14, face = "bold")) 
FigS10

# save plot
png("./figures/FigS10.png", res = 700, width = 9, height = 9, units = "in")
FigS10
dev.off()


## Shannon diversity (Appendix 1 content)

# Invertebrate Shannon diversity

# Year 1 post-summer

isd.y1.ps.herb <- glmmTMB(shannon.invert ~ trt_y1 + herb_trt+ (1|block), 
                          dispformula = ~trt_y1,
                          family = tweedie(), data = diversity %>% filter(date == "2019-10-20"))
plot(simulateResiduals(isd.y1.ps.herb))
summary(isd.y1.ps.herb)
Anova(isd.y1.ps.herb, type = 2) # herbivore treatments are not significant; drop this term

isd.y1.ps <- glmmTMB(shannon.invert ~ trt_y1 + (1|block), 
                         family = tweedie(), dispformula = ~trt_y1,
                         data = diversity %>% filter(date == "2019-10-20"))
plot(simulateResiduals(isd.y1.ps))
summary(isd.y1.ps)
Anova(isd.y1.ps, type = 2) 

# Year 1 winter

isd.y1.w.herb <- glmmTMB(shannon.invert ~ trt_y1 + herb_trt+ (1|block), 
                          family = tweedie(), data = diversity %>% filter(date == "2020-03-15"))
plot(simulateResiduals(isd.y1.w.herb))
summary(isd.y1.w.herb)
Anova(isd.y1.w.herb, type = 2) # herbivore treatments are not significant; drop this term

isd.y1.w <- glmmTMB(shannon.invert ~ trt_y1 + (1|block), 
                     family = tweedie(), dispformula = ~trt_y1,
                     data = diversity %>% filter(date == "2020-03-15"))
plot(simulateResiduals(isd.y1.w))
summary(isd.y1.w)
Anova(isd.y1.w, type = 2)

labels.S11a <- as.data.frame(cbind(
  c("C","W","C","W"),
  c("Post-summer", "Post-summer",
    "Winter", "Winter"),
  c(1.45,1.45,1.15,1.15),
  c("a","a","b","c"))
)
colnames(labels.S11a) <- c("treatment","season", "mean.isd","label")
labels.S11a <- labels.S11a %>% mutate(mean.isd = as.numeric(mean.isd))


# Year 2 post-summer
isd.y2.ps <- glmmTMB(shannon.invert ~ trt_y1*trt_y2 + (1|block), 
                       data = diversity %>% filter(date == "2020-09-14"))

plot(simulateResiduals(isd.y2.ps))
summary(isd.y2.ps)
Anova(isd.y2.ps, type = 3,contrasts=list(trt_y1 = "contr.sum", trt_y2 = "contr.sum"))

emm2.isd.ps <- emmeans(isd.y2.ps, specs = c("trt_y1","trt_y2"))
contrast(emm2.isd.ps, method = "pairwise", adjust = "tukey")

# Year 2 winter
isd.y2.w <- glmmTMB(shannon.invert ~ trt_y1*trt_y2 + (1|block), 
                     data = diversity %>% filter(date == "2021-02-24"))

plot(simulateResiduals(isd.y2.w))
summary(isd.y2.w)
Anova(isd.y2.w, type = 3,contrasts=list(trt_y1 = "contr.sum", trt_y2 = "contr.sum"))

emm2.isd.w <- emmeans(isd.y2.w, specs = c("trt_y1","trt_y2"))
contrast(emm2.isd.w, method = "pairwise", adjust = "tukey")

labels.S11b <- as.data.frame(cbind(
  c(rep(c("CC","CW","WC","WW"), times = 2)),
  c(rep("Post-summer", times = 4), rep("Winter", times = 4)),
  c(1.6,1.6,1.5,1.5,1.55,1.55,1.55,1.4),
  c("a","ab","bc","c", "d","de","de","e"))
)
colnames(labels.S11b) <- c("treatment","season","mean.isd","label")
labels.S11b <- labels.S11b %>% mutate(mean.isd = as.numeric(mean.isd))


# Invert Shannon diversity post-summer and wintertime in year 1 and year 2
FigS11a <- ggplot(diversity_plot_y1, aes(x = treatment, y = shannon.invert, 
                                        col = treatment, pch = season)) +
  geom_point(data = diversity_plot_y1_summary, aes(x = treatment, y = mean.isd),
             size = 2, stroke = 1, alpha = 0.9) +
  geom_errorbar(data = diversity_plot_y1_summary, 
                aes(y = mean.isd, 
                    ymax = mean.isd + se.isd, 
                    ymin = mean.isd - se.isd),
                width = 0.15, color = "grey20") +
  geom_jitter(alpha = 0.2, height = 0, width = 0.2) +
  facet_wrap(~season) +
  theme_classic() +
  scale_color_manual(values = pal.trt, guide = "none") +
  scale_shape_manual(values = c(1,2), guide = "none") +
  labs(x = "Treatment", pch = "Season", col = "Treatment", y = "Invertebrate Shannon diversity") +
  theme(plot.tag = element_text(face = "bold", size = 10),
        strip.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8)) +
  scale_y_continuous(limits = c(0,1.5), breaks = c(0,0.3,0.6,0.9,1.2,1.5)) +
  geom_text(data = labels.S11a, aes(x = treatment, y = mean.isd, label = label),
            col = "black", size = 3, fontface = "bold")
FigS11a

FigS11b <- ggplot(diversity_plot_y2, aes(x = treatment, y = shannon.invert, 
                                                col = treatment, pch = season)) +
  geom_jitter(alpha = 0.2, height = 0, width = 0.2, show.legend = F) +
  geom_errorbar(data = diversity_plot_y2_summary, aes(y = mean.isd, 
                                                      ymax = mean.isd + se.isd, 
                                                      ymin = mean.isd - se.isd),
                width = 0.15, color = "grey20") +
  geom_point(data = diversity_plot_y2_summary, aes(y = mean.isd), 
             size = 3, alpha = 0.8, show.legend = T) +
  theme_classic() +
  facet_wrap(~season) +
  labs(x = "Treatment", pch = "Season", col = "Treatment", y = "Invertebrate Shannon diversity") +
  geom_text(data = labels.S11b, aes(x = treatment, y = mean.isd, label = label),
            col = "black", size = 3, fontface = "bold") +
  scale_color_manual(values = pal.trt, drop = FALSE) +
  scale_y_continuous(limits = c(0,1.6), breaks = c(0,0.3,0.6,0.9,1.2,1.5)) +
  guides(color = guide_legend(override.aes = list(shape = pch.trt, size = 3))) +
  theme(plot.tag = element_text(face = "bold", size = 10),
        strip.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8),
        legend.box = "horizontal")
FigS11b


# Year 1: post-summer
asd.y1.ps.herb <- glmmTMB(shannon.algae ~ trt_y1 + herb_trt + (1|block), 
                         family = tweedie(), 
                         data = diversity %>% filter(date == "2019-10-20"))
plot(simulateResiduals(asd.y1.ps.herb))
summary(asd.y1.ps.herb)
Anova(asd.y1.ps.herb, type = 2) # significant! do not include this test or these data

asd.y1.w.herb

asd.y1.w.herb <- glmmTMB(shannon.algae ~ trt_y1 + herb_trt + (1|block), 
                              family = tweedie(), data = diversity %>% filter(date == "2020-03-15"))
plot(simulateResiduals(asd.y1.w.herb))
summary(asd.y1.w.herb)
Anova(asd.y1.w.herb, type = 2) # herbivory not significant; drop this term

asd.y1.w <- glmmTMB(shannon.algae ~ trt_y1 + (1|block), 
                         family = tweedie(), data = diversity %>% filter(date == "2020-03-15"))
plot(simulateResiduals(asd.y1.w))
summary(asd.y1.w)
Anova(asd.y1.w, type = 2) # herbivory not significant; drop this term

# Year 2: post-summer has too many zeroes, can't model data

# Year 2: winter

asd.y2.w <- glmmTMB(shannon.algae ~ trt_y1*trt_y2 + (1|block), 
                           family = tweedie(),
                           data = diversity %>% filter(date == "2021-02-24"))

plot(simulateResiduals(asd.y2.w))
summary(asd.y2.w)
Anova(asd.y2.w, type = 3,contrasts=list(trt_y1 = "contr.sum", trt_y2 = "contr.sum"))

emm2.asd.w <- emmeans(asd.y2.w, specs = c("trt_y1","trt_y2"))
contrast(emm2.asd.w, method = "pairwise", adjust = "tukey")

labels.S11c <- as.data.frame(cbind(
  c("Year 1", "Year 1", "Year 2","Year 2", "Year 2", "Year 2"),
  c("C","W","CC","CW","WC","WW"),
  c(0.8, 0.65,1.1,1.1,1.1,1.1),
  c("a","b","c","c","c","c"))
)
colnames(labels.S11c) <- c("period", "treatment", "mean.asd","label")
labels.S11c <- labels.S11c %>% mutate(mean.asd = as.numeric(mean.asd))

# mean Shannon diversity of algae post-summer and wintertime in year 1 and year 2

algae_plot <- diversity_keytimes %>% filter(date %in% c("2020-03-15", "2021-02-24"))
algae_plot_summary <- diversity_plots %>% filter(season == "Winter")

FigS11c <- ggplot(algae_plot, aes(x = treatment, y = shannon.algae, shape = treatment,
                                         col = treatment)) +
  geom_point(data = algae_plot_summary, aes(x = treatment, y = mean.asd),
             size = 2, stroke = 1, alpha = 0.9) +
  geom_errorbar(data = algae_plot_summary, 
                aes(y = mean.asd, 
                    ymax = mean.asd + se.asd, 
                    ymin = mean.asd - se.asd),
                width = 0.15, color = "grey20") +
  geom_jitter(alpha = 0.2, height = 0, width = 0.2) +
  facet_grid(cols = vars(period), scales = "free_x") +
  theme_classic() +
  scale_color_manual(values = pal.trt, guide = "none") +
  scale_shape_manual(values = c(2,2,17,17,17,17), guide = "none") +
  labs(x = "Treatment", pch = "Season", col = "Treatment", y = "Algal Shannon diversity") +
  theme(plot.tag = element_text(face = "bold", size = 10),
        strip.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8)) +
  scale_y_continuous(limits = c(0,1.1), breaks = c(0,0.2,0.4,0.6,0.8,1)) +
  geom_text(data = labels.S11c, aes(x = treatment, y = mean.asd, label = label),
            col = "black", size = 3, fontface = "bold")
FigS11c


layout <- "
AABBB
CCCCC
"

FigS11 <- (FigS11a + FigS11b + FigS11c) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") +
  plot_layout(guides = "collect", design = layout) & theme(plot.tag = element_text(size = 14, face = "bold")) 
FigS11

# save plot
png("./figures/FigS11.png", res = 700, width = 9, height = 7, units = "in")
FigS11
dev.off()


## Species richness and shannon diversity of epifauna from destructive sampling

# extract factor variables
epi_summer_factors <- epifauna_summer %>% select(block, number, treatment)

# calculate richness and shannon diversity of invertebrate epifauna
richness <- specnumber(epifauna_summer)
shannon <- diversity(epifauna_summer %>% ungroup() %>%  select(-block,-number,-treatment))

# put factors and calculated factors together
summer.epi.rich <- epi_summer_factors %>% cbind(richness) %>% cbind(shannon) %>% 
  mutate(period = "Post-summer")

# repeat for winter sampling point
epi_winter_factors <- epifauna_winter %>% select(block, number, treatment)

richness.winter <- specnumber(epifauna_winter)
shannon.winter <- diversity(epifauna_winter %>% ungroup() %>%  select(-block,-number,-treatment))

winter.epi.rich <- epi_winter_factors %>% cbind(richness.winter) %>% cbind(shannon.winter) %>% 
  rename(richness = richness.winter, shannon = shannon.winter) %>% 
  mutate(period = "Winter")

# join all data together
epi.rich <- summer.epi.rich %>% full_join(winter.epi.rich) %>% 
  mutate(trt_y1 = as.factor(substr(treatment, 1,1)), 
         trt_y2=as.factor(substr(treatment, 2,2)),
         period = as.factor(period))



## Generalized linear models of richness (Poisson distribution) with treatments
# and date of sampling (period; post-summer or winter)

# Block is included as a random intercept effect
destructive.rich.ps <- glmmTMB(richness ~ trt_y1*trt_y2 + (1|block),
                           family = poisson(),
                             data = epi.rich %>% filter(period == "Post-summer"))
plot(simulateResiduals(destructive.rich.ps))
summary(destructive.rich.ps)
Anova(destructive.rich.ps,type=3,contrasts=list(trt_y1 = "contr.sum", trt_y2 = "contr.sum"))

# Post-hoc comparisons between treatments
emm.rich.ps <- emmeans(destructive.rich.ps, specs = c("trt_y1","trt_y2"))
contrast(emm.rich.ps, method = "pairwise", adjust = "tukey")

destructive.rich.w <- glmmTMB(richness ~ trt_y1*trt_y2 + (1|block),
                               family = poisson(),
                               data = epi.rich %>% filter(period == "Winter"))
plot(simulateResiduals(destructive.rich.w))
summary(destructive.rich.w)
Anova(destructive.rich.w,type=3, contrasts=list(trt_y1 = "contr.sum", trt_y2 = "contr.sum"))

# Post-hoc comparisons between treatments
emm.rich.w <- emmeans(destructive.rich.w, specs = c("trt_y1","trt_y2"))
contrast(emm.rich.w, method = "pairwise", adjust = "tukey")

# Create dataframe with label values
labels.S12a <- as.data.frame(cbind(
  c(rep("Post-summer", times = 4), rep("Winter", times = 4)),
  c(rep(c("CC", "CW", "WC","WW"),times = 2)),
  c(16,16,16,10,19,19,19,11),
  c("a","ab","ab","b","a","ab","a","b"))
)

colnames(labels.S12a) <- c("period", "treatment", "mean_rich", "label")
labels.S12a <- labels.S12a %>% mutate(mean_rich = as.numeric(mean_rich))

# Generalized linear model of Shannon diversity of epifauna with same model structure
destructive.sd.ps <- glmmTMB(shannon ~ trt_y1*trt_y2 + (1|block),
                            data = epi.rich %>%  filter(period == "Post-summer"))
plot(simulateResiduals(destructive.sd.ps))
summary(destructive.sd.ps)
Anova(destructive.sd.ps,type=3, contrasts=list(trt_y1 = "contr.sum", trt_y2 = "contr.sum"))

emm.sd.ps <- emmeans(destructive.sd.ps, specs = c("trt_y1","trt_y2"))
contrast(emm.sd.ps, method = "pairwise", adjust = "tukey")

destructive.sd.w <- glmmTMB(shannon ~ trt_y1*trt_y2 + (1|block),
                             data = epi.rich %>%  filter(period == "Winter"))
plot(simulateResiduals(destructive.sd.w))
summary(destructive.sd.w)
Anova(destructive.sd.w,type=3, contrasts=list(trt_y1 = "contr.sum", trt_y2 = "contr.sum"))

emm.sd.w <- emmeans(destructive.sd.w, specs = c("trt_y1","trt_y2"))
contrast(emm.sd.w, method = "pairwise", adjust = "tukey")

labels.S12b <- as.data.frame(cbind(
  c(rep("Post-summer", times = 4), rep("Winter", times = 4)),
  c(rep(c("CC", "CW", "WC","WW"),times = 2)),
  c(2,2,2,1.5,2.45,2.45,2.45,1.95),
  c("a","ab","a","b","a","a","a","b"))
)

colnames(labels.S12b) <- c("period","treatment", "mean_sd","label")
labels.S12b <- labels.S12b %>% mutate(mean_sd = as.numeric(mean_sd))

# Create a plotting dataframe
epi.plot <- epi.rich %>% 
  group_by(period, treatment) %>% 
  summarize(mean_rich = mean(richness), se_rich = std.error(richness),
            mean_shannon = mean(shannon), se_shannon = std.error(shannon))

# Species richness of epifauna community
FigS12a <- ggplot(aes(x = treatment, y = mean_rich, col = treatment, pch = period), 
                  data = epi.plot) + 
  geom_errorbar(aes(ymin = mean_rich - se_rich, 
                    ymax = mean_rich + se_rich), width = 0.2,
                color = "grey20",
                position = position_dodge(width = 1)) +
  facet_wrap(~period) +
  geom_point(position = position_dodge(width = 1), size = 3, alpha = 0.8) +
  geom_jitter(data = epi.rich, aes(y = richness), alpha = 0.2, height = 0, width = 0.2) +
  scale_colour_manual(values = pal.trt.y2) +
  theme_classic() +
  theme(plot.tag = element_text(face = "bold"))+
  labs(x = "Treatment", y = "Species richness", col = "Treatment", pch = "Season") +
  scale_y_continuous(limits = c(4,19), breaks = c(4,6,8,10,12,14,16,18)) +
  geom_text(data = labels.S12a, 
            aes(x = treatment, y = mean_rich, label = label), 
            col = "black", fontface = "bold", inherit.aes = FALSE)
FigS12a

# Shannon diversity of epifauna community
FigS12b <- ggplot(aes(x = treatment, y = mean_shannon, col = treatment, pch = period,), data = epi.plot) + 
  facet_wrap(~period)+
  geom_errorbar(aes(ymin = mean_shannon - se_shannon, 
                    ymax = mean_shannon + se_shannon), width = 0.2,
                color = "grey20",
                position = position_dodge(width = 1)) +
  geom_point(position = position_dodge(width = 1), size = 3, alpha = 0.8) +
  geom_jitter(data = epi.rich, aes(y = shannon), alpha = 0.2, height = 0, width = 0.2) +
  scale_colour_manual(values = pal.trt.y2) +
  theme_classic() +
  theme(plot.tag = element_text(face = "bold"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))+
  labs(x = "Treatment", y = "Shannon diversity", pch = "Season", 
       color = "Treatment") +
  scale_y_continuous(limits = c(0,2.5), breaks = c(0,0.5,1,1.5,2,2.5)) +
  geom_text(data = labels.S12b, 
            aes(x = treatment, y = mean_sd, label = label), 
            col = "black", fontface = "bold", inherit.aes = FALSE)
FigS12b

# stitch together plot
FigS12 <- (FigS12a | FigS12b) + 
  plot_layout(guides = "collect") + 
  plot_annotation(tag_levels = "a", tag_prefix = "(",
                  tag_suffix = ")")
FigS12

# save plot
png("./figures/FigS12.png", res = 700, width = 8, height = 4.5,
    units = "in")
FigS12
dev.off()
