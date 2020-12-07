pkgs <- c("tidyverse", "vegan", "pair")
lapply(pkgs, library, character.only = TRUE)


# need to get community data in useful format
surveys <- read_csv("./clean_data/SVSHW_survey_clean.csv")
infauna_sept20 <- read_csv("./clean_data/SVSHW_infauna_20200914.csv")


# just pick out end of summer 2019 and end of summer 2020

comm_oct19 <- surveys %>%
  filter(date == "2019-10-20") %>%
  unite(abund, c(count, percent_cover), sep = "") %>% 
  mutate(abund = as.numeric(str_remove_all(abund, "NA"))) %>% 
  select(-date) %>% 
  pivot_wider(id_cols = c(block, number, treatment), names_from = species, values_from = abund,
              values_fill = list(abund = 0)) %>% 
  select(-mytilus, -number) %>% 
  mutate(annelida = ifelse(is.na(annelida), 0, annelida)) %>% 
  mutate(ulva = ifelse(is.na(ulva), 0, ulva))


comm_oct_factors <- comm_oct19 %>% 
  select(treatment, block) %>% 
  mutate(treatment = str_replace_all(treatment, c("WW" = "W", "WA" = "W",
                                                  "AA" = "A", "AW" = "A")))

#get rid of everything and convert to matrix
comm_matrix_oct19 <- comm_oct19 %>% 
  select(-block, - treatment) 
  
# convert to matrix
comm_matrix_oct19 <- as.matrix(comm_matrix_oct19)

# ready to analyze
comm_oct <- metaMDS(comm_matrix_oct19, k = 2, try = 400)

tiff("fall2019.tiff", units = "in", width = 6, height = 6, res = 450)

oct_ord <- ordiplot(comm_oct, type = "none", xlim = c(-2,2), ylim = c(-2,2))
points(oct_ord, "sites", cex = 0.7, col = "black")
orditorp(oct_ord,display="species",col="grey30", cex = 0.8, air=1)
ordiellipse(oct_ord, groups = comm_oct_factors$treatment, 
            col = c("blue", "red"), label = TRUE)
dev.off()

# PERMANOVA & PERMDISP

perm.y1 <- adonis(comm_matrix_oct19 ~ treatment, strata = comm_oct_factors$block, data = comm_oct_factors,
        perm = 99)
perm.y1

dist.y1 <- vegdist(comm_matrix_oct19, method = "bray")

disp.y1 <- betadisper(dist.y1, type = "centroid", group = comm_oct_factors$treatment)
anova(disp.y1)
boxplot(disp.y1)


# After year 2 of warming
foundation_sept20 <- surveys %>%
  filter(date == "2020-09-14") %>% 
  filter(species %in% c("balanus","chthamalus","ulva","ulothrix","pyropia","fucus")) %>% 
  rename(taxon = species) %>% 
  unite(abund, c(count, percent_cover), sep = "") %>% 
  mutate(abund = as.numeric(str_remove_all(abund, "NA")))

key_list <- paste(infauna_sept20$block, infauna_sept20$number, sep = "_")

comm_sept20 <- foundation_sept20 %>%
  full_join(infauna_sept20) %>% 
  unite(key, c(block, number), sep = "_", remove = FALSE) %>% 
  filter(key %in% key_list) %>% 
  select(-key) %>% 
  pivot_wider(id_cols = c(block, number, treatment), names_from = taxon, values_from = abund,
              values_fill = list(abund =0))

comm_sept_factors <- comm_sept20 %>% 
  select(treatment, block) %>% 
  mutate(trty1 = substring(treatment, 1, 1),
           trty2 = substring(treatment, 2,2))

comm_matrix_sept20 <- comm_sept20 %>% 
  dplyr::select(-block, - treatment, -number) %>% 
  as.matrix()

comm_sept <- metaMDS(comm_matrix_sept20, k = 2, try = 400)

tiff("fall2020.tiff", units = "in", width = 6, height = 6, res = 450)

sept_ord <- ordiplot(comm_sept, type = "none", xlim=c(-1,1),ylim = c(-1.5,1.5))
points(sept_ord, what = "sites", col = "black")
orditorp(sept_ord, display="species",col="grey30", cex = 0.8, air=1)
ordiellipse(comm_sept, groups = comm_sept_factors$treatment, 
            col = c("blue", "purple","orange","red"), label = TRUE)

dev.off()

# PERMANOVA & PERMDISP

perm.y2 <- adonis(comm_matrix_sept20 ~ trty1*trty2, data = comm_sept_factors,
                  perm = 99)
perm.y2

dist.y2 <- vegdist(comm_matrix_sept20, method = "bray")

disp.y2 <- betadisper(dist.y2, group = comm_sept_factors$treatment)
anova(disp.y2)
TukeyHSD(disp.y2)
boxplot(disp.y2)


# lump taxa and remove spurious things that aren't actually using barnacles long term/couldn't be IDed

comm_sept20_2 <- comm_sept20 %>% 
  mutate(amphipoda = amphipoda_1 + amphipoda_2,
         polychaeta = polychaeta_1 + polychaeta_2 +sabellidae,
         nemertea = nemertean_1 + emplectonema_gracile) %>% 
  select(-c(amphipoda_1, amphipoda_2, polychaeta_1,polychaeta_2,
            nemertean_1, emplectonema_gracile, insecta_2, insecta,
            limpet_recruit, copepoda, worm_thing, sabellidae,
            onchidoris_bliamellata))

comm_sept_factors_2 <- comm_sept20 %>% 
  select(block, treatment) %>% 
  mutate(trty1 = substring(treatment, 1, 1),
         trty2 = substring(treatment, 2,2))

comm_matrix_sept20_2 <- comm_sept20_2 %>% 
  dplyr::select(-block, - treatment, -number) %>% 
  as.matrix() 

comm_sept_2 <- metaMDS(comm_matrix_sept20_2, k = 2, try = 400)

tiff("fall2020_2.tiff", units = "in", width = 6, height = 6, res = 450)

sept_ord <- ordiplot(comm_sept_2, type = "none", xlim = c(-2.5,2.5))
points(sept_ord, what = "sites", col = "black")
orditorp(sept_ord, display="species",col="grey30", cex = 0.8, air=1)
ordiellipse(comm_sept_2, groups = comm_sept_factors$treatment, 
            col = c("blue", "purple","orange","red"), label = TRUE)

dev.off()

# PERMANOVA & PERMDISP

perm.y2.2 <- adonis(comm_matrix_sept20_2 ~ trty1*trty2, data = comm_sept_factors_2,
                  perm = 99)
perm.y2.2

dist.y2.2 <- vegdist(comm_matrix_sept20_2, method = "bray")

disp.y2.2 <- betadisper(dist.y2.2, type = "centroid", group = comm_sept_factors_2$treatment)
anova(disp.y2.2)
TukeyHSD(disp.y2.2)

tiff("fall2020_centroids.tiff", units = "in", width = 6, height = 4, res = 450)

boxplot(disp.y2.2, xlab = "Treatment", col = c("blue", "purple","orange","red"))

dev.off()

# spring 2020?

comm_mar20 <- surveys %>%
  filter(date == "2020-03-15") %>%
  unite(abund, c(count, percent_cover), sep = "") %>% 
  mutate(abund = as.numeric(str_remove_all(abund, "NA"))) %>% 
  select(-date) %>% 
  pivot_wider(id_cols = c(block, number, treatment), names_from = species, values_from = abund,
              values_fill = list(abund = 0)) %>% 
  select(-number)

# remove rows where nothing is present

comm_abund_mar20 <- comm_mar20 %>% 
select(-block, -treatment)

totals <- rowSums(comm_abund_mar20, na.rm = TRUE)

bare <- which(as.numeric(totals) == 0)

bare_comm <- comm_abund_mar20 %>% 
  slice(bare)

full_comm <- comm_mar20 %>% 
  anti_join(bare_comm)

#extract treatment information

comm_mar_factors <- full_comm %>% 
  select(treatment, block) %>% 
  mutate(treatment = str_replace_all(treatment, c("WW" = "W", "WA" = "W",
                                                  "AA" = "A", "AW" = "A")))

comm_matrix_mar20 <- full_comm %>% 
  anti_join(bare_comm) %>% 
  select(-block, -treatment) %>% 
  as.matrix()

# ready to analyze
comm_mar <- metaMDS(comm_matrix_mar20, k = 2, try = 400)

tiff("spring20.tiff", units = "in", width = 6, height = 6, res = 450)

mar_ord <- ordiplot(comm_mar, type = "none", xlim = c(-3,3), ylim = c(-2,2))
points(mar_ord, "sites", cex = 0.7, col = "black")
orditorp(mar_ord, display="species",col="grey30", cex = 0.8, air=1)
ordiellipse(mar_ord, groups = comm_mar_factors$treatment, 
            col = c("blue", "red"), label = TRUE)
dev.off()

# PERMANOVA & PERMDISP

perm.y1s <- adonis(comm_matrix_mar20 ~ treatment, strata = comm_mar_factors$block, data = comm_mar_factors,
                  perm = 99)
perm.y1s

dist.y1s <- vegdist(comm_matrix_mar20, method = "bray")

disp.y1s <- betadisper(dist.y1s, type = "centroid", group = comm_mar_factors$treatment)
anova(disp.y1s)
boxplot(disp.y1s)

# just infauna for fall 2020

infauna_only <- infauna_sept20 %>%
  pivot_wider(id_cols = c(block, number, treatment), names_from = taxon, values_from = abund,
            values_fill = list(abund =0)) %>% 
  mutate(amphipoda = amphipoda_1 + amphipoda_2,
         polychaeta = polychaeta_1 + polychaeta_2 +sabellidae,
         nemertea = nemertean_1 + emplectonema_gracile) %>% 
  select(-c(amphipoda_1, amphipoda_2, polychaeta_1,polychaeta_2,
            nemertean_1, emplectonema_gracile, insecta_2, insecta,
            limpet_recruit, copepoda, worm_thing, sabellidae,
            onchidoris_bliamellata))

infauna_factors <- infauna_only %>% 
  select(block, number, treatment) %>% 
  mutate(trty1 = substring(treatment, 1, 1),
         trty2 = substring(treatment, 2,2))

infauna_matrix <- infauna_only %>% 
  select(-block, -number, - treatment) %>% 
  as.matrix()

infaun_sept_1 <- metaMDS(infauna_matrix, k = 2, try = 400)

tiff("infauna.tiff", units = "in", width = 6, height = 6, res = 450)

inf_ord <- ordiplot(infaun_sept_1, type = "none")
points(inf_ord, what = "sites", col = "black")
orditorp(inf_ord, display="species",col="grey30", cex = 0.8, air=1)
ordiellipse(infaun_sept_1, groups = infauna_factors$treatment, 
            col = c("blue", "purple","orange","red"), label = TRUE)

dev.off()

perm.inf <- adonis(infauna_matrix ~ trty1*trty2, data = infauna_factors,
                    perm = 99)
perm.inf

dist.inf <- vegdist(infauna_matrix, method = "bray")

disp.inf <- betadisper(dist.inf, type = "centroid", group = comm_sept_factors_2$treatment)
anova(disp.inf)
TukeyHSD(disp.inf)
