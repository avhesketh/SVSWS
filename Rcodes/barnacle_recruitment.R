# barnacle recruitment spring 2020

se <- function(x) {
  sd(x)/sqrt(length(x))
}

library(tidyverse)
library(glmmTMB)
library(car)
library(DHARMa)

recruitment <- read_csv("./raw_data/tile_surveys/SVSHW_bncle_recruit.csv") %>% 
  mutate(treatment = ifelse(trt == "A" & consecutive == "N", "AW",
                            ifelse(trt == "A" & consecutive == "Y", "AA",
                                   ifelse(trt == "W" & consecutive == "N", "WA", "WW")))) 

# negative binomial for a first approximation
# 
# but let's plot it first!

bal_rec <- recruitment %>% 
  filter(species == "balanus" & size == "recruit")

bal_adult <- recruitment %>% 
  filter(species == "balanus" & size == "adult")

chtham_rec <- recruitment %>% 
  filter(species == "chthamalus" & size == "recruit")

chtham_adult <- recruitment %>% 
  filter(species == "chthamalus" & size == "adult")

barn_summ <- recruitment %>% 
  group_by(treatment, species, size) %>% 
  summarize(mean_abund = mean(count), se_abund = se(count))
barn_summ
barn_summ$treatment <- factor(barn_summ$treatment, levels = c("AA","WA","AW","WW"))
?as.factor
bal_fig <- ggplot(aes(x = trt, y = mean_abund, fill = treatment), data = barn_summ) +
  facet_grid(species ~ size) +
  theme_bw() +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymax = mean_abund + se_abund, ymin = mean_abund - se_abund), 
                position = position_dodge(width = 0.9),
                width = 0.5) +
  ylab("Abundance") +
  xlab("Temperature treatment") +
  theme(plot.title = element_text(hjust = 0.5, size = 16)) +
  theme(axis.title = element_text(size = 16)) +
  theme(axis.text = element_text(size = 12)) +
  theme(axis.text.x = element_blank()) +
  theme(strip.text = element_text(size = 14)) +
  theme(legend.title = element_text(size = 14)) +
  theme(legend.text = element_text(size = 12)) +
  scale_fill_manual(values = c("blue","purple","orange","red"))
bal_fig


chtham_fig <- ggplot(aes(x = trt, y = mean_abund, fill = consecutive), data = chtham_summ) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymax = mean_abund + se_abund, ymin = mean_abund - se_abund), 
                position = position_dodge(width = 0.9),
                width = 0.5) +
  facet_wrap(~size) +
  theme_classic() +
  ggtitle("Size class") +
  ylab("Abundance") +
  xlab("Temperature treatment") +
  theme(plot.title = element_text(hjust = 0.5, size = 16)) +
  theme(axis.title = element_text(size = 16)) +
  theme(axis.text = element_text(size = 12)) +
  theme(strip.text = element_text(size = 14)) +
  theme(legend.title = element_text(size = 14)) +
  theme(legend.text = element_text(size = 12)) +
  scale_fill_manual(values = c("grey30", "grey70"))
chtham_fig

# now for a model

brec.glmm.1 <- glmmTMB(count ~ treatment
                      + (1 | block), data = bal_rec,
                      dispformula = ~block,
                      family = nbinom1())
dh.brec.1 <- simulateResiduals(brec.glmm.1)
plot(dh.brec.1)

summary(brec.glmm.1)
Anova(brec.glmm.1)

bad.glmm.1 <- glmmTMB(count ~ treatment
                      + (1 | block), data = bal_adult,
                      family = nbinom1())
dh.bad.1 <- plot(simulateResiduals(bad.glmm.1))

summary(bad.glmm.1)
Anova(bad.glmm.1)


crec.glmm.1 <- glmmTMB(count ~ treatment
                       + (1 | block), data = chtham_rec,
                       family = nbinom1())
dh.crec.1 <- simulateResiduals(crec.glmm.1)
plot(dh.crec.1)

summary(crec.glmm.1)
Anova(crec.glmm.1)

cad.glmm.1 <- glmmTMB(count ~ treatment
                        + (1 | block), data = chtham_adult,
                      family = nbinom1())
dh.cad.1 <- simulateResiduals(cad.glmm.1)
plot(dh.cad.1)

summary(cad.glmm.1)
Anova(cad.glmm.1)