# look for effect of herbivores

library(tidyverse)

cover <- read_csv("limpets_august.csv")

cover <- cover %>% 
  rename(abundance = number_1)

cover$temp <- gsub("ambeint", "ambient", cover$temp)

# exploratory stuff

balanus <- cover %>% 
  filter(species == "balanus") %>% 
  select(-cover)

balanus$treatment <- as.factor(balanus$treatment)
balanus$temp <- as.factor(balanus$temp)

ulothrix <- cover %>% 
  filter(species == "ulothrix") %>% 
  select(-abundance)

# plot balanus

balanus <- balanus %>% 
  mutate(log_abundance = log(abundance+1))

balanus_summary <- balanus %>% 
  group_by(treatment, temp) %>% 
  summarize(av_abundance = mean(log_abundance), se_abundance = sd(log_abundance)/length(abundance))

bal_abund <- ggplot(aes(x = treatment, y = av_abundance, fill = temp), data = balanus_summary) +
  geom_bar(stat = "identity", position = "dodge") +
  xlab("limpet treatment") +
  ylab("Balanus abundance (log scale)") +
  geom_errorbar(aes(ymax = av_abundance + se_abundance, ymin = av_abundance - se_abundance), position = "dodge") +
  theme_classic()

bal_abund

# greater variance in abundance depending on the treatment.temp combination
# use varIdent variance structure?

M1_bal <- lm(abundance ~ treatment*temp, data = balanus)
plot(resid(M1_bal)~balanus$treatment)
plot(resid(M1_bal)~balanus$temp)

#normalize abundance data

M2_bal <- gls(log_abundance~treatment +temp, data=balanus, method = "REML")
anova(M2_bal)
summary(M2_bal)
plot(M2_bal)
resid.M2bal <- residuals(M2_bal, type = "normalized")
plot(resid.M2bal~treatment, data = balanus)

# add block effect
library(nlme)
M3_bal <- lme(log_abundance ~ treatment*temp, random = ~1|block, data=balanus)
summary(M3_bal)
anova(M2_bal, M3_bal)


# try Poisson since it's count data
library(lme4)
M4_bal <- glmer(abundance ~ treatment*temp + (1 + number| block), family = "poisson", data = balanus)
anova(M4_bal)
summary(M4_bal)
plot(M4_bal, pch = balanus$temp)

# zero-inflated data
hist(balanus$abundance, breaks = 20)

#download glmmTMB and load... package that does mixed models with zero inflation
library(glmmTMB)
M5_bal <- glmmTMB(abundance ~ treatment*temp + (1|block),
                  data = balanus,
                  family = poisson,
                  ziformula=~1)
summary(M5_bal)

M6_bal <- glmmTMB(abundance ~ treatment*temp + (1|block),
                  data = balanus,
                  family = poisson,
                  ziformula=~0)
summary(M6_bal)

M7_bal <- glmmTMB(abundance ~ treatment*temp + (1|block),
                  data = balanus,
                  family = nbinom2,
                  ziformula=~treatment*temp)
summary(M7_bal)

#by AIC M7_bal is best. now plot assumptions
resid.M7 <- residuals(M7_bal)
plot(resid.M7 ~ balanus$abundance) # trend of higher resid at higher abundances
plot(resid.M7 ~ balanus$temp) #greater variance at higher temp
plot(resid.M7 ~ balanus$treatment) #dig pdig has higher variation
pred.M7 <- balanus[,1:6] 
abund_pred <- predict(M7_bal, newdata = pred.M7, exclude = "s(block)", type = "link")
head(abund_pred)
plot(balanus$abundance ~ transform.pred)
transform.pred <- exp(abund_pred)


## same with ulothrix
ulothrix <- cover %>% 
  filter(species == "ulothrix")

ulothrix$treatment <- as.factor(ulothrix$treatment)
ulothrix$temp <- as.factor(ulothrix$temp)

# plot ulothrix

hist(asin(sqrt(ulothrix$cover/100)))

ulothrix <- ulothrix %>% 
  mutate(asin_cover = asin(sqrt(cover/100)))

ulo_summary <- ulothrix %>% 
  group_by(treatment, temp) %>% 
  summarize(av_cover = mean(asin_cover), se_cover = sd(asin_cover)/length(cover))

ulo_cover <- ggplot(aes(x = treatment, y = av_cover, fill = temp), data = ulo_summary) +
  geom_bar(stat = "identity", position = "dodge") +
  xlab("limpet treatment") +
  ylab("Ulothrix cover (arcsine sqrt transformed)") +
  geom_errorbar(aes(ymax = av_cover + se_cover, ymin = av_cover - se_cover), position = "dodge") +
  theme_classic()

ulo_cover

M1_ulo <- gls(asin_cover ~ treatment*temp, data = ulothrix)
anova(M1_ulo)

M2_ulo <- lme(asin_cover ~ temp, random = ~1|block, data = ulothrix)
anova(M2_ulo)

resid <- residuals(M2_ulo)
plot(resid~temp, data = ulothrix)
