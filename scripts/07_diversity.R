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
  select(species, abund, survey_no, tile_id, original_herb_trt) %>% 
  unique() %>% 
  pivot_wider(names_from = species, values_from = abund,values_fill = 0)

# calculate species richness and Shannon diversity
diversity.factors <- diversity %>% 
  select(survey_no, tile_id, original_herb_trt)

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

richness <- specnumber(diversity %>% ungroup() %>% select(-survey_no, -tile_id, -original_herb_trt))
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
  mutate(original_herb_trt = if_else(is.na(original_herb_trt), "control",original_herb_trt))

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

# Year 1
richness.y1 <- glmmTMB(richness ~ trt_y1*date + (1|block), 
                       family = poisson(),
                       data = diversity %>% filter(date %in% c("2019-10-20","2020-03-15")))
richness.y1.herb <- glmmTMB(richness ~ trt_y1*date + original_herb_trt + (1|block), 
                         family = poisson(), 
                         data = diversity %>% filter(date %in% c("2019-10-20","2020-03-15")))
AIC(richness.y1, richness.y1.herb) # better without herbivores 

plot(simulateResiduals(richness.y1))
summary(richness.y1)
Anova(richness.y1, type = 3)

emm1.rich <- emmeans(richness.y1, specs = "trt_y1", by = "date")
contrast(emm1.rich, method = "pairwise", adjust = "tukey")


# Year 2
# including random effect causes singularity, so dropped from model here 

richness.y2 <- glmmTMB(richness ~ trt_y1*trt_y2 + date + (1|block),
                       family = poisson(),
                       data = diversity %>% filter(date %in% c("2020-09-14","2021-02-24")))
richness.y2.herb <- glmmTMB(richness ~ trt_y1*trt_y2 + date + original_herb_trt + (1|block), 
                       family = poisson(),
                       data = diversity %>% filter(date %in% c("2020-09-14", "2021-02-24")))

AIC(richness.y2, richness.y2.herb) # better without herbivores
plot(simulateResiduals(richness.y2)) # model is under-dispersed, which will yield more conservative p values
# no other distribution is appropriate, so we will keep this as-is despite violations
summary(richness.y2)
Anova(richness.y2, type = 3)

emm2.rich <- emmeans(richness.y2, specs = c("trt_y1","trt_y2"))
contrast(emm2.rich, method = "pairwise", adjust = "tukey")

# create dataframe with significance labels for eventual plot
labels.5b <- as.data.frame(cbind(
  c("CC","CW","WC","WW"),
  c(6.7,5.05,5.85,4.3),
  c("a","ab","a","b"))
)
colnames(labels.5b) <- c("treatment","richness","label")
labels.5b <- labels.5b %>% mutate(richness = as.numeric(richness))


################
# plots of alpha diversity metrics at key times (post summer and in winter)

# define palettes for plotting
pal.trt <- c("#014779", "#EE4B2B", "#014779", "#7985CB", "#9C0098", "#EE4B2B", "grey50")
pch.trt <- c(1,1,16,16,16,16,16)
lty.trt <- c("dotdash","dotdash","solid","solid","solid","solid","solid","solid")

# Create summary dataframe
diversity_keytimes <- diversity %>% 
  filter(date %in% c("2019-10-20", "2020-03-15", "2020-09-14","2021-02-24")) %>% 
  mutate(period = if_else(date %in% c("2019-10-20", "2020-03-15"), "Year 1", "Year 2"),
         timept = if_else(date %in% c("2019-10-20", "2020-09-14"), "Post-summer", "Winter")) %>% 
  mutate(treatment = if_else(period == "Year 1", trt_y1, treatment),
         treatment = factor(treatment, levels = c("C","W","CC","CW","WC","WW")),
         date = factor(date))

# Plot as single point ± one standard error, so summarize mean and standard error
diversity_plots <- diversity_keytimes %>% 
  group_by(treatment, timept, period) %>% 
  summarize(mean.rich = mean(richness), se.rich = std.error(richness),
            mean.isd = mean(shannon.invert), se.isd = std.error(shannon.invert),
            mean.asd = mean(shannon.algae), se.asd = std.error(shannon.algae)) %>% ungroup()

# Extract year 1 and year 2 data for separate plots
diversity_plot_y1 <- diversity_plots %>% filter(period == "Year 1")
diversity_plot_y2 <- diversity_plots %>% filter(period == "Year 2")

# mean richness post-summer and wintertime in year 1
Fig5a <- ggplot(diversity_plot_y1, aes(x = treatment, y = mean.rich, 
                                       col = treatment, pch = timept)) +
  geom_point(position = position_dodge(width = 1), size = 2.5) +
  geom_errorbar(aes(ymax = mean.rich + se.rich, ymin = mean.rich - se.rich),
                position = position_dodge(width = 1), width = 0.3) +
  theme_classic() +
  scale_color_manual(values = pal.trt, guide = "none") +
  scale_shape_manual(values = c(1,2), guide = "none") +
  labs(x = "Treatment", pch = "Time of year", col = "Treatment", y = "Species richness") +
  theme(plot.tag = element_text(face = "bold", size = 10),
        strip.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8)) +
  scale_y_continuous(limits = c(0,7), breaks = c(0,1,2,3,4,5,6,7)) +
  annotate("segment", x = 0.75, xend = 0.75, y = 5.3, yend = 5.7, col = "black") +
  annotate("segment", x = 1.75, xend = 1.75, y = 5.3, yend = 5.7, col = "black") +
  annotate("segment", x = 0.75, xend = 1.75, y = 5.5, yend = 5.5, col = "black") +
  annotate("segment", x = 1.25, xend = 1.25, y = 4.1, yend = 4.5, col = "black") +
  annotate("segment", x = 2.25, xend = 2.25, y = 4.1, yend = 4.5, col = "black") +
  annotate("segment", x = 1.25, xend = 2.25, y = 4.3, yend = 4.3, col = "black") +
  annotate("text", x = 1.25, y = 5.95, label = "p < 0.001", size = 2.5) +
  annotate("text", x = 1.75, y = 4.75, label = "p < 0.001", size = 2.5)
Fig5a

# mean richness post-summer and wintertime in year 2

Fig5b <- ggplot(diversity_plot_y2, aes(x = treatment, y = mean.rich, 
                                       col = treatment, pch = timept)) +
  geom_point(position = position_dodge(width = 1), size = 3, alpha = 0.8) +
  geom_errorbar(aes(ymax = mean.rich + se.rich, ymin = mean.rich - se.rich),
                position = position_dodge(width = 1), width = 0.3) +
  theme_classic() +
  scale_color_manual(values = pal.trt.y2, guide = "none") +
  scale_shape_manual(values = c(16,17)) +
  labs(x = "Treatment", pch = "Time of year", col = "Treatment", y = "Species richness") +
  scale_y_continuous(limits = c(0,7), breaks = c(0,1,2,3,4,5,6,7)) +
  geom_text(data = labels.5b, 
            aes(x = treatment, y = richness, label = label), 
            col = "black", fontface = "bold", inherit.aes = FALSE, size = 2.5) +
  theme(plot.tag = element_text(face = "bold", size = 10),
        strip.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))
  
Fig5b

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
Fig5c <- ggplot() + 
  # add ellipses
  geom_path(data=df_ell_summer, aes(x=dbRDA1,y=dbRDA2,colour=treatment)) +
  # add points for each tile
  geom_point(data=data.scores.summer, aes(x=dbRDA1,y=dbRDA2,colour=treatment),
             size=2, alpha = 0.8) + # add the point markers
  scale_color_manual(values = pal.trt.y2) +
  theme_classic()+
  labs(color = "Treatment", x = "dbRDA1", y = "dbRDA2") +
  theme(plot.tag = element_text(face = "bold", size = 10),
        strip.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8)) +
  # add labels showing results of pairwise comparisons
  annotate(geom = "text", x = -0.9, y = 1, label = "a", col ="#014779", size = 2.5, fontface = "bold") +
  annotate(geom = "text", x = -0.5, y = -0.53, label = "ab", col = "#7985CB", size = 2.5, fontface = "bold") +
  annotate(geom = "text", x = -0.08, y = -0.7, label = "ab", col = "#9C0098", size = 2.5, fontface = "bold") +
  annotate(geom = "text", x = 1.25, y = 0.8, label = "b", col = "#EE4B2B", size = 2.5, fontface = "bold") 

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
             size=2, alpha = 0.8, shape = 17) + # add the point markers
  scale_color_manual(values = pal.trt.y2) +
  labs(color = "Treatment", x = "dbRDA1", y = "dbRDA2") +
  theme_classic()  + 
  theme(plot.tag = element_text(face = "bold"),
        legend.position = "none") +
  annotate(geom = "text", x = 1, y = 0.7, label = "c", col ="#014779", size = 2.5, fontface = "bold") +
  annotate(geom = "text", x = -0.6, y = -0.8, label = "ab", col = "#7985CB", size = 2.5, fontface = "bold") +
  annotate(geom = "text", x = 0.52, y = -0.5, label = "b", col = "#9C0098", size = 2.5, fontface = "bold") +
  annotate(geom = "text", x = -0.9, y = 0.6, label = "a", col = "#EE4B2B", size = 2.5, fontface = "bold") 
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

# assemble multipanel figure
Fig5 <- ((Fig5a | Fig5b) / ((Fig5c) | (Fig5d))) +
  plot_layout(guides = "collect", widths = c(0.6, 0.4)) +
  plot_annotation(tag_levels = "a", tag_prefix = "(",
                  tag_suffix = ")") 
Fig5

# save figure
ggsave(Fig5, filename = "./figures/Fig5.pdf", device = cairo_pdf, 
       width = 16, height = 12, units = "cm")

################ 
# Appendix content: Figures & analyses

# Plot of species richness of treatments through time
FigS8a <- ggplot(diversity.plot, aes(y = mean.rich, x =date, 
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
  theme(legend.key.width = unit(1,"cm"), plot.tag = element_text(face = "bold")) +
  scale_x_date(date_labels = "%Y-%m", breaks=breaks)
FigS8a

# Shannon diversity of invertebrate community over time
FigS8b <- ggplot(diversity.plot, aes(y = mean.shann.i, x =date, col = treatment,lty = treatment, pch=treatment)) +
  theme_classic() +                            
  geom_point() +
  geom_line(lwd = 0.8) +
  scale_color_manual(values = pal.trt) +
  scale_linetype_manual(values = lty.trt) +
  scale_shape_manual(values = pch.trt) +
  theme_classic() + theme(legend.key.width = unit(1,"cm"), plot.tag = element_text(face = "bold")) +
  geom_errorbar(aes(ymax = mean.shann.i + se.shann.i, ymin = mean.shann.i - se.shann.i, lty = NULL)) +
  labs(y = "Invert. Shannon diversity", x = "Date", pch = "Treatment", col = "Treatment", lty = "Treatment") +
  geom_vline(aes(xintercept = as.Date("2020-04-07")), linetype = "dotted", col = "grey30", lwd = 0.8) +
  scale_x_date(date_labels = "%Y-%m", breaks=breaks)
FigS8b

# Shannon diversity of algal community over time
FigS8c <- ggplot(diversity.plot, aes(y = mean.shann.a, x =date, 
                                     col = treatment, pch = treatment, lty = treatment)) +
  geom_point() +
  geom_line(lwd = 0.8) +
  scale_color_manual(values = pal.trt) +
  scale_linetype_manual(values = lty.trt) +
  scale_shape_manual(values = pch.trt) +
  theme_classic() + theme(legend.key.width = unit(1,"cm"),plot.tag = element_text(face = "bold")) +
  geom_errorbar(aes(ymax = mean.shann.a + se.shann.a, ymin = mean.shann.a - se.shann.a, lty = NULL)) +
  labs(y = "Algal Shannon diversity", x = "Date", pch = "Treatment", col = "Treatment", lty = "Treatment") +
  geom_vline(aes(xintercept = as.Date("2020-04-07")), linetype = "dotted", col = "grey30", lwd = 0.8) +
  scale_x_date(date_labels = "%Y-%m", breaks=breaks)
FigS8c

# stitch them all together
FigS8 <- FigS8a / FigS8b / FigS8c +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") +
  plot_layout(guides = "collect") & theme(plot.tag = element_text(size = 14, face = "bold")) 
FigS8

# save plot
png("./figures/FigS8.png", res = 700, width = 9, height = 9, units = "in")
FigS8
dev.off()


## Shannon diversity (Appendix 1 content)

# Invertebrate Shannon diversity

# Year 1
isd.y1 <- glmmTMB(shannon.invert ~ trt_y1*date + (1|block), 
                         family = tweedie(), dispformula = ~trt_y1,
                         data = diversity %>% filter(date %in% c("2019-10-20", "2020-03-15")))
isd.y1.herb <- glmmTMB(shannon.invert ~ trt_y1*date + original_herb_trt+ (1|block), 
                              dispformula = ~trt_y1,
                              family = tweedie(), data = diversity %>% filter(date %in% c("2019-10-20", "2020-03-15")))

AIC(isd.y1, isd.y1.herb) # better without herbivores 

plot(simulateResiduals(isd.y1))
summary(isd.y1)
Anova(isd.y1, type = 3)

emm1.isd <- emmeans(isd.y1, specs = "trt_y1", by = "date")
contrast(emm1.isd, method = "pairwise", adjust = "tukey")

# Year 2

isd.y2 <- glmmTMB(shannon.invert ~ trt_y1*trt_y2 + date + (1|block),
                         data = diversity %>% filter(date %in% c("2020-09-14","2021-02-24")))
isd.y2.herb <- glmmTMB(shannon.invert ~ trt_y1*trt_y2 + date + original_herb_trt + (1|block), 
                              data = diversity %>% filter(date %in% c("2020-09-14","2021-02-24")))

AIC(isd.y2, isd.y2.herb) # better without herbivores

plot(simulateResiduals(isd.y2))
summary(isd.y2)
Anova(isd.y2, type = 3)

emm2.isd <- emmeans(isd.y2, specs = c("trt_y1","trt_y2"))
contrast(emm2.isd, method = "pairwise", adjust = "tukey")

labels.S9b <- as.data.frame(cbind(
  c("CC","CW","WC","WW"),
  c(1.35,1.25,1.15,1.05),
  c("a","ab","bc","c"))
)
colnames(labels.S9b) <- c("treatment","mean.isd","label")
labels.S9b <- labels.S9b %>% mutate(mean.isd = as.numeric(mean.isd))


# Invert Shannon diversity post-summer and wintertime in year 1 and year 2
FigS9a <- ggplot(diversity_plot_y1, aes(x = treatment, y = mean.isd, 
                                        col = treatment, pch = timept)) +
  geom_point(position = position_dodge(width = 1), size = 2.5) +
  geom_errorbar(aes(ymax = mean.isd + se.isd, ymin = mean.isd - se.isd),
                position = position_dodge(width = 1), width = 0.3) +
  theme_classic() +
  scale_color_manual(values = pal.trt, guide = "none") +
  scale_shape_manual(values = c(1,2), guide = "none") +
  labs(x = "Treatment", pch = "Time of year", col = "Treatment", y = "Invertebrate Shannon diversity") +
  theme(plot.tag = element_text(face = "bold"), 
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 11)) +
  scale_y_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5), limits = c(0,0.55)) +
  annotate("segment", x = 0.75, xend = 0.75, y = 0.51, yend = 0.53, col = "black") +
  annotate("segment", x = 1.75, xend = 1.75, y = 0.51, yend = 0.53, col = "black") +
  annotate("segment", x = 0.75, xend = 1.75, y = 0.52, yend = 0.52, col = "black") +
  annotate("segment", x = 1.25, xend = 1.25, y = 0.43, yend = 0.45, col = "black") +
  annotate("segment", x = 2.25, xend = 2.25, y = 0.43, yend = 0.45, col = "black") +
  annotate("segment", x = 1.25, xend = 2.25, y = 0.44, yend = 0.44, col = "black") +
  annotate("text", x = 1.25, y = 0.55, label = "ns") +
  annotate("text", x = 1.75, y = 0.475, label = "p < 0.001")
FigS9a

FigS9b <- ggplot(diversity_plot_y2, aes(x = treatment, y = mean.isd, 
                                       col = treatment, pch = timept)) +
  geom_point(position = position_dodge(width = 1), size = 3, alpha = 0.8) +
  geom_errorbar(aes(ymax = mean.isd + se.isd, ymin = mean.isd - se.isd),
                position = position_dodge(width = 1), width = 0.3) +
  theme_classic() +
  scale_color_manual(values = pal.trt.y2, guide = "none") +
  scale_shape_manual(values = c(16,17)) +
  labs(x = "Treatment", pch = "Time of year", col = "Treatment", y = "Invertebrate Shannon diversity") +
  theme(plot.tag = element_text(face = "bold")) +
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1,1.25,1.5), limits = c(0,1.35)) +
  geom_text(data = labels.S9b, 
            aes(x = treatment, y = mean.isd, label = label), 
            col = "black", fontface = "bold", inherit.aes = FALSE)

FigS9b

pal.trt.time <- c("transparent", "transparent", "#014779", "#7985CB" ,"#9C0098", "#EE4B2B")

# Year 1: post-summer
asd.y1 <- glmmTMB(shannon.algae ~ trt_y1 + date + (1|block), 
                         family = tweedie(), 
                         data = diversity %>% filter(date %in% c("2019-10-20","2020-03-15")))
asd.y1.herb <- glmmTMB(shannon.algae ~ trt_y1 + date + original_herb_trt + (1|block), 
                              family = tweedie(), data = diversity %>% filter(date %in% c("2019-10-20","2020-03-15")))

AIC(asd.y1, asd.y1.herb) # better without herbivores 

plot(simulateResiduals(asd.stress.y1))
summary(asd.y1)
Anova(asd.y1, type = 2)

emm1.asd <- emmeans(asd.y1, specs = c("trt_y1"))
contrast(emm1.asd, method = "pairwise", adjust = "tukey")


# Year 2: post-summer has too many zeroes, can't model data

# Year 2: winter

asd.recovery.y2 <- glmmTMB(shannon.algae ~ trt_y1*trt_y2 + (1|block), 
                           family = tweedie(),
                           data = diversity %>% filter(date == "2021-02-24"))

# adding herbivore term causes convergence. Based on other models, will exclude from model
#asd.recovery.y2.herb <- glmmTMB(shannon.algae ~ trt_y1*trt_y2 + original_herb_trt +(1|block), 
                             #   family = tweedie(),
                             #   data = diversity %>% filter(date == "2021-02-24"))

plot(simulateResiduals(asd.recovery.y2))
summary(asd.recovery.y2)
Anova(asd.recovery.y2, type = 3)

emm2.asd <- emmeans(asd.recovery.y2, specs = c("trt_y1","trt_y2"))
contrast(emm2.asd, method = "pairwise", adjust = "tukey")

labels.S9c <- as.data.frame(cbind(
  c(1,2,1.25,2.25,3.25,4.25),
  c("Year 1", "Year 1", "Year 2","Year 2", "Year 2", "Year 2"),
  c(0.21, 0.11,0.36, 0.265,0.61,0.22),
  c("a","b","c","c","c","c"))
)
colnames(labels.S9c) <- c("x","period", "mean.asd","label")
labels.S9c <- labels.A9c %>% mutate(mean.asd = as.numeric(mean.asd),
                                    x = as.numeric(x))

# mean Shannon diversity of algae post-summer and wintertime in year 1 and year 2

FigS9c <- ggplot(diversity_plots, aes(x = treatment, y = mean.asd, fill = treatment, col = treatment, pch = timept)) +
  geom_point(position = position_dodge(width = 1), size = 3, alpha = 0.8) +
  geom_errorbar(aes(ymax = mean.asd + se.asd, ymin = mean.asd - se.asd),
                position = position_dodge(width = 1), width = 0.3) +
  facet_wrap(~period, scales ="free_x") +
  theme_classic() +
  scale_color_manual(values = pal.trt) +
  scale_fill_manual(values = pal.trt.time, guide = "none") +
  scale_alpha_manual(values = c(0.1,1), guide = "none") +
  scale_shape_manual(values = c(21,24), guide = "none") +
  labs(x = "Treatment", pch = "Time of year", 
       col = "Treatment", y = "Algal Shannon diversity") +
  theme(plot.tag = element_text(face = "bold")) +
  geom_text(data = labels.S9c, 
            aes(x = x, y = mean.asd, label = label), 
            col = "black", fontface = "bold", inherit.aes = FALSE)
FigS9c

FigS9 <- ((FigS9a | FigS9b)/FigS9c) + plot_annotation(tag_levels = "a",
                                                tag_prefix = "(",
                                                tag_suffix = ")") + 
  plot_layout(guides = "collect")

png("./figures/FigS9.png", res = 700, width = 7.5, height = 6,
    units = "in")
FigS9
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
destructive.rich <- glmmTMB(richness ~ trt_y1*trt_y2 + period + (1|block),
                           family = poisson(),
                             data = epi.rich)
plot(simulateResiduals(destructive.rich))
summary(destructive.rich)
Anova(destructive.rich,type=3)

# Post-hoc comparisons between treatments
emm.rich <- emmeans(destructive.rich, specs = c("trt_y1","trt_y2"))
contrast(emm.rich, method = "pairwise", adjust = "tukey")

# Create dataframe with label values
labels.S10a <- as.data.frame(cbind(
  c("CC", "CW", "WC","WW"),
  c(12.65,9,10,6.9),
  c("a","bc","ab","c"))
)

colnames(labels.S10a) <- c("treatment", "mean_rich","label")
labels.S10a <- labels.S10a %>% mutate(mean_rich = as.numeric(mean_rich))

# Generalized linear model of Shannon diversity of epifauna with same model structure
destructive.sd <- glmmTMB(shannon ~ trt_y1*trt_y2 + period + (1|block),
                            data = epi.rich)
plot(simulateResiduals(destructive.sd))
summary(destructive.sd)
Anova(destructive.sd,type=3)

emm.sd <- emmeans(destructive.sd, specs = c("trt_y1","trt_y2"))
contrast(emm.sd, method = "pairwise", adjust = "tukey")

labels.S10b <- as.data.frame(cbind(
  c("CC", "CW", "WC","WW"),
  c(1.8, 1.5,1.5,0.98),
  c("a","a","a","b"))
)

colnames(labels.S10b) <- c("treatment", "mean_sd","label")
labels.S10b <- labels.S10b %>% mutate(mean_sd = as.numeric(mean_sd))

# Create a plotting dataframe
epi.plot <- epi.rich %>% 
  group_by(period, treatment) %>% 
  summarize(mean_rich = mean(richness), se_rich = std.error(richness),
            mean_shannon = mean(shannon), se_shannon = std.error(shannon))

# Species richness of epifauna community
FigS10a <- ggplot(aes(x = treatment, y = mean_rich, col = treatment, pch = period), 
                  data = epi.plot) + 
  geom_point(position = position_dodge(width = 1), size = 3, alpha = 0.8) +
  geom_errorbar(aes(ymin = mean_rich - se_rich, 
                    ymax = mean_rich + se_rich), width = 0.5,
                position = position_dodge(width = 1)) +
  scale_colour_manual(values = pal.trt.y2) +
  theme_classic() +
  theme(plot.tag = element_text(face = "bold"))+
  labs(x = "Treatment", y = "Species richness", col = "Treatment", pch = "Time of year") +
  geom_text(data = labels.S10a, 
            aes(x = treatment, y = mean_rich, label = label), 
            col = "black", fontface = "bold", inherit.aes = FALSE)
FigS10a

# Shannon diversity of epifauna community
FigS10b <- ggplot(aes(x = treatment, y = mean_shannon, col = treatment, pch = period,), data = epi.plot) + 
  geom_point(position = position_dodge(width = 1), size = 3, alpha = 0.8) +
  geom_errorbar(aes(ymin = mean_shannon - se_shannon, 
                    ymax = mean_shannon + se_shannon), width = 0.5,
                position = position_dodge(width = 1)) +
  scale_colour_manual(values = pal.trt.y2) +
  theme_classic() +
  theme(plot.tag = element_text(face = "bold"))+
  labs(x = "Treatment", y = "Shannon diversity", pch = "Time of year", 
       color = "Treatment") +
  geom_text(data = labels.S10b, 
            aes(x = treatment, y = mean_sd, label = label), 
            col = "black", fontface = "bold", inherit.aes = FALSE)
FigS10b

# stitch together plot
FigS10 <- (FigS10a | FigS10b) + 
  plot_layout(guides = "collect") + 
  plot_annotation(tag_levels = "a", tag_prefix = "(",
                  tag_suffix = ")")
FigS10

# save plot
png("./figures/FigS10.png", res = 700, width = 8, height = 4.5,
    units = "in")
FigS10
dev.off()
