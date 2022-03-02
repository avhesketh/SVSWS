# barnacle recruitment spring 2020

se <- function(x) {
  sd(x)/sqrt(length(x))
}

library(tidyverse)
library(glmmTMB)
library(car)
library(DHARMa)

recruitment <- read_csv("./raw_data/tile_surveys/SVSHW_bncle_recruit.csv") %>%
  rename(new_block = block, new_no = number) %>% 
  left_join(block_design) %>% 
  mutate(species = if_else(species == "balanus", "Balanus glandula", "Chthamalus dalli")) %>% 
  select(new_block, new_no, tile_id, date, species, size, count, treatment) %>% 
  left_join(survey_join) 


# negative binomial for a first approximation
# 
# but let's plot it first!

barnacle_summary <- recruitment %>% group_by(species, size,treatment, survey_no) %>% 
  summarize(mean_abund = mean(count), se_abund = std.error(count))

time_continuous <- barnacle_summary %>% left_join(survey_join) %>% 
  mutate(timesincestart = difftime(date, ymd("2019-04-12"), units = c("weeks"))) %>% 
  group_by(survey_no) %>% mutate(timesincestart = max(timesincestart)) %>% ungroup() %>% unique()

barn_fig <- ggplot(aes(x = timesincestart, y = mean_abund, col = treatment), data = time_continuous) +
  facet_grid(species~size) +
  theme_bw() +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymax = mean_abund + se_abund, ymin = mean_abund - se_abund), 
                width = 1) +
  labs(y = "Mean abundance", x = "Time since experiment start (weeks)", col = "Treatment") +
  scale_color_manual(values = c("blue","purple","orange","darkred"))
barn_fig

## recruits only

recruitment_only <- recruitment %>% filter(size == "recruit" & survey_no == min(survey_no,na.rm=T))

#write_csv(recruitment_only, "./plotting_df/recruitment.csv")

only_rec <- ggplot(aes(x = species, y = count, col = treatment), data = recruitment_only) +
  geom_boxplot() +
  labs(y = "Barnacle abundance", x = "Temperature treatment", col = "Treatment") +
  scale_color_manual(values = c("blue","purple","orange","darkred")) +
  theme_classic()
only_rec

# now for models

model_df <- recruitment_only %>% mutate(trt_y1 = substr(treatment ,1,1), trt_y2 = substr(treatment, 2,2))

bal.rec <- glmmTMB(count ~ trt_y1*trt_y2 + (1|new_block), family = nbinom1(), data = model_df %>% filter(species == "Balanus glandula"))
chtham.rec <- glmmTMB(count ~ trt_y1*trt_y2 + (1|new_block), family = nbinom1(), data = model_df %>% filter(species == "Chthamalus dalli"))

plot(simulateResiduals(bal.rec))
summary(bal.rec)
Anova(bal.rec, type = 3)

plot(simulateResiduals(chtham.rec))
summary(chtham.rec)
Anova(chtham.rec)
                   