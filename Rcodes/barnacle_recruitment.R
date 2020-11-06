# barnacle recruitment spring 2020

se <- function(x) {
  sd(x)/sqrt(length(x))
}

library(tidyverse)
library(glmmTMB)
library(car)
library(DHARMa)

bar_rec <- read_csv("barnacle_recruitment.csv")

bar_rec$species <- as.factor(bar_rec$species)
bar_rec$size <- as.factor(bar_rec$size)

hist(bar_rec$abund)

# negative binomial for a first approximation
# 
# but let's plot it first!

bal_rec <- bar_rec %>% 
  filter(species == "balanus")

chtham_rec <- bar_rec %>% 
  filter(species == "chthamalus")

bal_summ <- bal_rec %>% 
  group_by(trt, size, consecutive) %>% 
  summarize(mean_abund = mean(abund), se_abund = se(abund))

chtham_summ <- chtham_rec %>% 
  group_by(trt, size, consecutive) %>% 
  summarize(mean_abund = mean(abund), se_abund = se(abund))

bal_summ$size <- gsub("adult", "> 2 mm", bal_summ$size)
bal_summ$size <- gsub("recruit", "< 2 mm", bal_summ$size)
chtham_summ$size <- gsub("adult", "> 2 mm", chtham_summ$size)
chtham_summ$size <- gsub("recruit", "< 2 mm", chtham_summ$size)

bal_fig <- ggplot(aes(x = trt, y = mean_abund, fill = consecutive), data = bal_summ) +
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

brec.glmm.1 <- glmmTMB(abund ~ (trt + consecutive + size)^3
                      + (1 | block), data = bal_rec,
                      family = nbinom1())
dh.brec.1 <- simulateResiduals(brec.glmm.1)
plot(dh.brec.1)

summary(brec.glmm.1)
Anova(brec.glmm.1)

crec.glmm.1 <- glmmTMB(abund ~ (trt + consecutive + size)^3
                       + (1 | block), data = chtham_rec,
                       family = nbinom1())
dh.crec.1 <- simulateResiduals(crec.glmm.1)
plot(dh.crec.1)

drop1(crec.glmm.1, test = "Chisq")

crec.glmm.2 <- update(crec.glmm.1, ~. - trt:consecutive:size)

drop1(crec.glmm.2, test = "Chisq")

crec.glmm.3 <- update(crec.glmm.2, ~. - trt:size)

drop1(crec.glmm.3, test = "Chisq")

crec.glmm.4 <- update(crec.glmm.3, ~. -consecutive:size)

drop1(crec.glmm.4, test = "Chisq")

dh.crec.4 <- simulateResiduals(crec.glmm.4)
plot(dh.crec.4)
print(crec.glmm.4)
summary(crec.glmm.4)
Anova(crec.glmm.4)

## July 2020

# get the barnacle data

barn_jul <- read.csv("barnacle_recruitment_jul.csv")

bal_summ <- barn_jul %>% 
  filter(species == "balanus")

chtham_summ <- barn_jul %>% 
  filter(species == "chthamalus")

bal_jul_summ <- bal_summ %>% 
  group_by(trt, size, species, consecutive) %>% 
  summarize(mean_count = mean(count), se_count = se(count))

chtham_jul_summ <- chtham_summ %>% 
  group_by(trt, size, species, consecutive) %>% 
  summarize(mean_count = mean(count), se_count = se(count))


bal_fig <- ggplot(aes(x = trt, y = mean_count, fill = consecutive), data = bal_jul_summ) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymax = mean_count + se_count, ymin = mean_count - se_count), 
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
bal_fig

chtham_fig <- ggplot(aes(x = trt, y = mean_count, fill = consecutive), data = chtham_jul_summ) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymax = mean_count + se_count, ymin = mean_count - se_count), 
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

hist(bal_summ$count, breaks = 100)

bal.jul.1 <- glmmTMB(count ~ trt*consecutive*size + (1|block),
                     family = nbinom1(),
                     data = bal_summ)

plot(simulateResiduals(bal.jul.1))

summary(bal.jul.1)

drop1(bal.jul.1, test = "Chisq")

bal.jul.2 <- update(bal.jul.1, ~. -trt:consecutive:size)

drop1(bal.jul.2, test = "Chisq")

bal.jul.3 <- update(bal.jul.2, ~. -consecutive:size)

Anova(bal.jul.3)

## now chthamalus

chtham.jul.1 <- glmmTMB(count ~ trt*consecutive*size + (1|block),
                        family = nbinom1(),
                        data = chtham_summ)

summary(chtham.jul.1)

drop1(chtham.jul.1, test = "Chisq")

chtham.jul.2 <- update(chtham.jul.1, ~. -trt:consecutive:size)

drop1(chtham.jul.2, test = "Chisq")

chtham.jul.3 <- update(chtham.jul.2, ~. -trt:consecutive)

drop1(chtham.jul.3, test = "Chisq")

chtham.jul.4 <- update(chtham.jul.3, ~. - trt:size)

drop1(chtham.jul.4, test = "Chisq")

chtham.jul.5 <- update(chtham.jul.4, ~. -consecutive:size)

plot(simulateResiduals(chtham.jul.5))

Anova(chtham.jul.5)
