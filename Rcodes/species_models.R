#abundance data stuff
pkgs <- c("tidyverse","glmmTMB","DHARMa", "plotrix", "lubridate","gamm4")
lapply(pkgs, library, character.only = TRUE)
rm(pkgs)

## Balanus glandula

surveys <- read_csv("./clean_data/SVSHW_survey_clean.csv")
survey_numbers <- read_csv("./raw_data/design/SVSHW_survey_times.csv") %>% 
  group_by(survey_no) %>% 
  summarize(date = max(date))

balanus_y1 <- surveys %>% 
  left_join(survey_numbers) %>% 
  filter(species == "Balanus_glandula" & date <= "2020-03-15") %>% 
  mutate(timesincestart = difftime(date, ymd("2019-04-12"), units = c("weeks")),
         count = as.integer(count)) %>% 
  mutate(treatment = trt_y1)

balanus_y2 <- surveys %>% 
  left_join(survey_numbers) %>% 
  filter(species == "Balanus_glandula" & date >= "2020-03-15" & (treatment %in% c("W","C"))==F) %>% 
  mutate(timesincestart = difftime(date, ymd("2019-04-12"), units = c("weeks")),
         count = as.integer(count))
         
balanus_all <- balanus_y1 %>% full_join(balanus_y2) %>% 
  mutate(treatment = factor(treatment, levels = c("C","W","CC","CW","WC","WW"))) %>% 
  mutate(month = month(date))

time_continuous <- balanus_all %>% select(survey_no, block, timesincestart) %>% 
  group_by(survey_no) %>% mutate(timesincestart = max(timesincestart)) %>% ungroup() %>% unique()

balanus_summarized <- balanus_all %>% group_by(survey_no, treatment) %>% 
  summarize(mean_abund = mean(count), se_abund = std.error(count)) %>% 
  full_join(time_continuous)

write_csv(balanus_summarized, "./plotting_df/balanus_time.csv")

treatment_colours <- c("cornflowerblue","tomato3","blue","purple","orange","darkred")

ggplot(aes(x = timesincestart, y = mean_abund, col = treatment), data = balanus_summarized) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = treatment_colours) +
  geom_errorbar(aes(ymin = mean_abund - se_abund, ymax = mean_abund + se_abund), width = 1) +
  labs(x = "Time since experiment start (weeks)", 
       y = expression("Mean abundance of"~italic("Balanus glandula")),
       col = "Treatment") +
  theme_classic()

balanus_time_m1 <- glmmTMB(count ~ trt_y1*trt_y2 + poly(month, 4) + (1|block),
                           data = balanus_all, family = nbinom1)
summary(balanus_time_m1)

plot(simulateResiduals(balanus_time_m1))
hist(residuals(balanus_time_m1)) # super overdispersed - pattern of barnacles through time is NOT linear
# Try using a gamm to capture recruitment dynamics

balanus_gamm <- balanus_all %>% mutate(trt_y1 = if_else(trt_y1 == "C", 0 , 1),
                                       trt_y2 = if_else(trt_y2 == "C", 0 , 1),
                                       month = month(date),
                                       year = year(date),
                                       block = case_when(block == "A" ~ 1,
                                                         block == "B" ~ 2,
                                                         block == "C" ~ 3,
                                                         block == "D" ~ 4,
                                                         block == "E" ~ 5))
balanus_gamm_y1 <- balanus_all %>% filter(date <= "2020-03-15")
balanus_gamm_y2 <- balanus_all %>% filter(date >= "2020-03-15")


bt_gamm.y1 <- gamm4(count ~ s(month, bs = "cc", k = 4) + trt_y1, random = ~(1|tile_id),
                         data = balanus_gamm_y1, family = negative.binomial(1))

summary(bt_gamm.y1$mer)

plot.gam(bt_gamm.y1$gam)
plot(bt_gamm.y1$mer)
anova(bt_gamm.y1$gam)

# doesn't look bad

bt_gamm.y2 <- gamm4(count ~ s(month, bs = "cc", k = 4) + trt_y2*trt_y1, random = ~(1|tile_id),
                    data = balanus_gamm_y2, family = negative.binomial(1))

summary(bt_gamm.y2$gam)
anova(bt_gamm.y2$gam)

plot.gam(bt_gamm.y2$gam)
plot(bt_gamm.y2$mer)


# final timepoint

balanus_final <- balanus_all %>% filter(survey_no == max(survey_no))

bal.final <- glmmTMB(count ~ trt_y1*trt_y2 + (1|tile_id), 
                  data = balanus_final,
                  family = nbinom1)

plot(simulateResiduals(bal.final))
summary(bal.final)
Anova(bal.final)


## now for ulothrix

ulo_y1 <- surveys %>% 
  left_join(survey_numbers) %>% 
  filter(species == "Ulothrix_sp" & date <= "2020-03-15") %>% 
  mutate(timesincestart = difftime(date, ymd("2019-04-12"), units = c("weeks"))) %>% 
  mutate(treatment = trt_y1)

ulo_y2 <- surveys %>% 
  left_join(survey_numbers) %>% 
  filter(species == "Ulothrix_sp" & date >= "2020-03-15") %>% 
  mutate(timesincestart = difftime(date, ymd("2019-04-12"), units = c("weeks")))

ulo_all <- ulo_y1 %>% full_join(ulo_y2) %>% 
  mutate(treatment = factor(treatment, levels = c("C","W","CC","CW","WC","WW"))) %>% 
  mutate(month = month(date))

time_continuous <- ulo_all %>% select(survey_no, block, date, timesincestart) %>% 
  group_by(survey_no) %>% mutate(date = max(date), timesincestart = max(timesincestart)) %>% ungroup() %>% unique()

ulo_summarized <- ulo_all %>% group_by(survey_no, treatment) %>% 
  summarize(mean_abund = mean(percent_cover), se_abund = std.error(percent_cover)) %>% 
  full_join(time_continuous)

write_csv(ulo_summarized, "./plotting_df/ulo_time.csv")


ulo_gamm <- ulo_all %>% mutate(trt_y1 = if_else(trt_y1 == "C", 0 , 1),
                                       trt_y2 = if_else(trt_y2 == "C", 0 , 1),
                                       month = month(date),
                                       year = year(date),
                                       block = case_when(block == "A" ~ 1,
                                                         block == "B" ~ 2,
                                                         block == "C" ~ 3,
                                                         block == "D" ~ 4,
                                                         block == "E" ~ 5))
ulo_y1 <- ulo_all %>% filter(date <= "2020-03-15")
ulo_y2 <- ulo_all %>% filter(date >= "2020-03-15")

M1_ulo <- gamm4(percent_cover ~ s(month, bs = "cc", k = 4) + trt_y1*trt_y2, random = ~(1|tile_id),
                data = ulo_y2)
plot.gam(M1_ulo$gam)
plot.gam(M1_ulo$mer)
anova(M1_ulo$gam)


peak_timept <- glmmTMB(percent_cover/100 ~ trt_y1,
                       family = beta_family(), dispformula = ~trt_y1,
                       data = ulo_y1 %>% filter(date == "2019-08-14"))
plot(simulateResiduals(peak_timept))
Anova(peak_timept)

ulo_y1_oct <- ulo_y1 %>%
  filter(date == "2019-10-20") %>% 
  mutate(percent_cover = (percent_cover*(length(percent_cover) -1) + 0.5)/length(percent_cover)) %>% 
  mutate(percent_cover = as.numeric(percent_cover), trty1 = as.factor(trty1))

M2_ulo <- glmmTMB(percent_cover/100 ~ trty1 + (1|block),
                  data = ulo_y1_oct,
                  family = beta_family())

plot(simulateResiduals(M1_ulo))

Anova(M1_ulo)

### chthamalus

chtham_y1 <- surveys %>% 
  left_join(survey_numbers) %>% 
  filter(species == "Chthamalus_dalli" & date <= "2020-03-15") %>% 
  mutate(timesincestart = difftime(date, ymd("2019-04-12"), units = c("weeks")),
         count = as.integer(count)) %>% 
  mutate(treatment = trt_y1)

chtham_y2 <- surveys %>% 
  left_join(survey_numbers) %>% 
  filter(species == "Chthamalus_dalli" & date >= "2020-03-15") %>% 
  mutate(timesincestart = difftime(date, ymd("2019-04-12"), units = c("weeks")),
         count = as.integer(count)) 

chtham_all <- chtham_y1 %>% full_join(chtham_y2) %>% 
  mutate(treatment = factor(treatment, levels = c("C","W","CC","CW","WC","WW")))

time_continuous <- chtham_all %>% select(survey_no, block, timesincestart) %>% 
  group_by(survey_no) %>% mutate(timesincestart = max(timesincestart)) %>% ungroup() %>% unique()

chtham_summarized <- chtham_all %>% group_by(survey_no, treatment) %>% 
  summarize(mean_abund = mean(count), se_abund = std.error(count)) %>% full_join(time_continuous)

write_csv(chtham_summarized, "./plotting_df/chthamalus_time.csv")

treatment_colours <- c("cornflowerblue","tomato3","blue","purple","orange","darkred")

ggplot(aes(x = timesincestart, y = mean_abund, col = treatment), data = chtham_summarized) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = treatment_colours) +
  geom_errorbar(aes(ymin = mean_abund - se_abund, ymax = mean_abund + se_abund), width = 1) +
  labs(x = "Time since experiment start (weeks)", 
       y = expression("Mean abundance of"~italic("Chthamalus dalli")),
       col = "Treatment") +
  theme_classic()

## models

chtham_gamm <- chtham_all %>% mutate(trt_y1 = if_else(trt_y1 == "C", 0 , 1),
                                       trt_y2 = if_else(trt_y2 == "C", 0 , 1),
                                       month = month(date),
                                       year = year(date),
                                       block = case_when(block == "A" ~ 1,
                                                         block == "B" ~ 2,
                                                         block == "C" ~ 3,
                                                         block == "D" ~ 4,
                                                         block == "E" ~ 5))
chtham_gamm_y2 <- chtham_gamm %>% filter(date >= "2020-03-15")


ct_gamm.y2 <- gamm4(count ~ s(month, bs = "cc", k =4) + trt_y1*trt_y2, random = ~(1|tile_id),
                    data = chtham_gamm_y2, family = negative.binomial(1))

summary(ct_gamm.y2$gam)

plot(ct_gamm.y2$gam)
plot(ct_gamm.y2$mer)
anova(ct_gamm.y2$gam)

# doesn't look bad

chtham_final <- chtham_all %>% filter(survey_no == max(survey_no))

chtham.final <- glmmTMB(count ~ trt_y1*trt_y2 + (1|tile_id), 
                     data = chtham_final,
                     family = nbinom1)
plot(simulateResiduals(chtham.final))
summary(chtham.final)
Anova(chtham.final)


## all the grazers!

limpets <- c("Lottia_digitalis","Lottia_scutum","Lottia_paradigitalis",
             "Lottia_pelta","Lottia_sp_recruits")
littorines <- c("Littorina_sitkana","Littorina_scutulata")

limpets_y2 <- surveys %>% 
  left_join(survey_numbers) %>% 
  filter(species %in% limpets & date >= "2020-03-15") %>% 
  mutate(timesincestart = difftime(date, ymd("2019-04-12"), units = c("weeks")),
         count = as.integer(count)) %>% 
  group_by(date, treatment, timesincestart,block, tile_id, trt_y1, trt_y2) %>% 
  summarise(total_limpets = sum(count)) %>% 
  filter((treatment %in% c("C","W")) == F)

limpet.abund <- glmmTMB(total_limpets ~ trt_y1*trt_y2 + timesincestart
                        + (1|date/tile_id) + (1|block), data = limpets_y2,
                        family = nbinom2())
plot(simulateResiduals(limpet.abund))
summary(limpet.abund)
Anova(limpet.abund, type = 3)

litts_y2 <- surveys %>% 
  left_join(survey_numbers) %>% 
  filter(species %in% littorines & date >= "2020-03-15") %>% 
  mutate(timesincestart = difftime(date, ymd("2019-04-12"), units = c("weeks")),
         count = as.integer(count)) %>% 
  group_by(date, treatment, timesincestart,block, tile_id, trt_y1, trt_y2) %>% 
  summarise(total_litts = sum(count)) %>% 
  filter((treatment %in% c("C","W")) == F)

litt.abund <- glmmTMB(total_litts ~ trt_y1*trt_y2 + timesincestart
                        + (1|date/tile_id) + (1|block), dispformula = ~treatment,
                        data = litts_y2,
                        family = nbinom2())
plot(simulateResiduals(litt.abund))
summary(litt.abund)
Anova(litt.abund, type = 3)

## Alpha diversity metrics

diversity <- surveys %>% 
  left_join(survey_numbers) %>% 
  mutate(timesincestart = difftime(date, ymd("2019-04-12"), units = c("weeks")),
         count = as.integer(count)) %>% 
  group_by(survey_no) %>% 
  mutate(timesincestart = max(timesincestart),
         date = max(date)) %>% 
  mutate(abund = if_else(is.na(count), percent_cover, as.numeric(count))) %>% 
  select(species, abund, survey_no, tile_id) %>% 
  unique() %>% 
  pivot_wider(names_from = species, values_from = abund,values_fill = 0)

diversity.factors <- diversity %>% 
  select(survey_no, tile_id)

diversity.invert.abund <- diversity %>% 
  ungroup() %>% 
  select(any_of(invert_list$species))

diversity.algae.abund <- diversity %>% 
  ungroup() %>% 
  select(any_of(alga_list$species))
  

richness <- specnumber(diversity)
shannon.algae <- diversity(diversity.algae.abund)
shannon.invert <- diversity(diversity.invert.abund)

diversity.responses <- diversity.factors %>% cbind(richness) %>% 
  cbind(shannon.algae) %>% cbind(shannon.invert) %>% 
  rename(richness = 3, shannon.algae = 4, shannon.invert = 5) %>% 
  left_join(tile_info) %>% left_join(survey_numbers) %>% 
  mutate(timesincestart = difftime(date, ymd("2019-04-12"), units = c("weeks")))

diversity_y1 <- diversity.responses %>% 
  filter(date <= "2020-03-15") %>% 
  mutate(treatment = trt_y1)

diversity_y2 <- diversity.responses %>% 
  filter(date >= "2020-03-15" & (treatment %in% c("W","C"))==F)
  
diversity <- diversity_y1 %>% full_join(diversity_y2) %>% 
  mutate(treatment = factor(treatment, levels = c("C","W","CC","CW","WC","WW"))) 

diversity.plot <- diversity %>% group_by(timesincestart, date, treatment) %>% 
  summarize(mean.rich = mean(richness),se.rich = std.error(richness),
            mean.shann.i = mean(shannon.invert), se.shann.i = std.error(shannon.invert),
            mean.shann.a = mean(shannon.algae), se.shann.a = std.error(shannon.algae))

write_csv(diversity.plot, "./plotting_df/diversity.csv")

richness.time <- glmmTMB(richness~ as.numeric(timesincestart) + trt_y1*trt_y2 + (1|block),
                         family = genpois(),
                         dispformula = ~treatment,
                      data = diversity_y2)
plot(simulateResiduals(richness.time))
summary(richness.time)
Anova(richness.time, type = 3)

# peak diversity only

richness.peak <- glmmTMB(richness~ trt_y1*trt_y2 + (1|block),
                         family = genpois(),
                         data = diversity_y2 %>% filter(date == "2020-12-11"))
plot(simulateResiduals(richness.peak))
summary(richness.peak)
Anova(richness.peak, type = 3)

# Shannon, inverts

si.time <- glmmTMB(shannon.invert ~ as.numeric(timesincestart) + trt_y1*trt_y2 + (1|block),
                         data = diversity_y2)
plot(simulateResiduals(si.time))
summary(si.time)
Anova(si.time, type = 3)

# Shannon, algae

sa.time <- glmmTMB(shannon.algae ~ as.numeric(timesincestart) + trt_y1*trt_y2 + (1|block),
                   ziformula = ~treatment,
                   data = diversity_y2)
plot(simulateResiduals(sa.time))
summary(sa.time)
Anova(si.time, type = 3)

# algal cover

algae <- read_csv("./clean_data/SVSHW_survey_clean.csv") %>% 
  filter(species %in% alga_list$species) %>% 
  left_join(survey_numbers)

algae_y1 <- algae %>% 
  filter(date <= "2020-03-15") %>% 
  mutate(treatment = trt_y1) %>% 
  group_by(date, block, tile_id, trt_y1, trt_y2, treatment) %>% 
  summarize(algal_cover = sum(percent_cover))

algae_y2 <- algae %>% 
  filter(date >= "2020-03-15" & (treatment %in% c("W","C"))==F) %>% 
  group_by(date, block, tile_id, trt_y1, trt_y2, treatment) %>% 
  summarize(algal_cover = sum(percent_cover))

algae_all <- algae_y1 %>% full_join(algae_y2)

algae_summ <- algae_all %>% 
  group_by(date, treatment) %>% 
  summarize(mean_cover = mean(algal_cover), se_cover = std.error(algal_cover))

write_csv(algae_summ, "./plotting_df/algal_cover.csv")

# beta family with transform doesn't work. Tweedie w zero inflation?

cover.time <- glmmTMB(algal_cover ~ date + trt_y1*trt_y2, family = tweedie(),
                      ziformula = ~date,
                      data = algae_y2)
plot(simulateResiduals(cover.time))
Anova(cover.time, type = 3)
