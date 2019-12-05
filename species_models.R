#abundance data stuff

library(tidyverse)
library(glmmTMB)

surveys <- read_csv("surveys.csv")

glimpse(surveys)

block_design <- read_csv("moving_guide.csv")

# need to first assign treatments to each tile, but 
# first get old tile assignments translated to new.

surveys$date <- as.Date(surveys$date, format = "%y-%m-%d")

survey_postmove <- surveys %>% 
  filter(date >= "2019-06-06") %>% 
  unite(b_no, c(block, number), sep = "_", remove = TRUE)

survey_premove <- surveys %>% 
  filter(date < "2019-06-06") %>% 
  unite(b_no, c(block, number), sep = "_", remove = TRUE)

joint_legend <- block_design %>% 
  unite(b_no_original, c(original_block, original_no), sep = "_", remove = TRUE) %>% 
  unite(b_no_new, c(new_block, new_no), sep = "_", remove = TRUE)

# now compare each entry against block design guide. if a tile
# was moved, then replace the old tile number/block with the 
# one after the move
for (i in 1:length(survey_premove$b_no)){
  for (j in 1:length(joint_legend$b_no_original)){
    if (survey_premove$b_no[i] == joint_legend$b_no_original[j]){
      survey_premove$b_no[i] = joint_legend$b_no_new[j]
    }
  }
}

block_new <- block_design %>%
  select(-original_block, -original_no) %>% 
  rename(block = new_block, number = new_no)

surveys_clean <- survey_premove %>%
  full_join(survey_postmove) %>% 
  select(-X7, -X8) %>% 
  separate(b_no, c("block","number"), sep = "_")

surveys_clean$date <- gsub("2019-05-08", "2019-05-09", surveys_clean$date)
surveys_clean$date <- gsub("2019-06-04", "2019-06-05", surveys_clean$date)
surveys_clean$date <- gsub("2019-07-03", "2019-07-04", surveys_clean$date)
surveys_clean$date <- gsub("2019-07-17", "2019-07-18", surveys_clean$date)
surveys_clean$date <- gsub("2019-07-30", "2019-07-31", surveys_clean$date)
surveys_clean$date <- gsub("2019-10-18", "2019-10-20", surveys_clean$date)

surveys_clean$number <- as.numeric(surveys_clean$number)
surveys_clean$date <- as.Date(surveys_clean$date)

surveys_clean <- surveys_clean %>% 
  full_join(block_new) %>% 
  mutate(colour = if_else(colour == "black", "warm", "ambient")) %>% 
  mutate(timediff = difftime(date, "2019-04-12", unit = "weeks"))

write.csv(surveys_clean, "clean_survey_data.csv")

# barnacles in diff treatments

surveys_clean <- read_csv("clean_survey_data.csv")

balanus <- surveys_clean %>% 
  select(-X1, -percent_cover) %>% 
  na.omit() %>% 
  filter(species == "balanus")

write.csv(balanus, "balanus.csv")

# ulothrix in diff treatments

ulothrix <- surveys_clean %>% 
  select(-count) %>% 
  na.omit() %>% 
  filter(species == "ulothrix")

write.csv(ulothrix, "ulothrix.csv")

library(lme4)

glimpse(balanus)

balanus$block <- as.factor(balanus$block)
balanus$colour <- as.factor(balanus$colour)
balanus$timediff <- as.numeric(balanus$timediff)
balanus$count <- as.integer(balanus$count)

hist(balanus$count, breaks = 100)

M1_bal <- glmmTMB(count ~ timediff*colour + (1|block), data = balanus, family = nbinom2)

summary(M1_bal)

library(DHARMa)

M1_bal_resdh <- simulateResiduals(M1_bal)
plot(M1_bal_resdh)

# residuals normal but patterned.

M1_bal_res <- residuals(M1_bal)

plot(M1_bal_res ~ balanus$timediff) # more variance at the start
plot(M1_bal_res ~ balanus$colour) # more variance for white tiles

# add time error structure??

M2_bal <- glmmTMB(count ~ timediff*colour + (1|block), 
                  data = balanus, family = nbinom2,
                  dispformula = ~timediff)
AIC(M1_bal, M2_bal) # better AIC with dispformula
anova(M1_bal, M2_bal) # and M2 significatnly better by anova comparison

summary(M2_bal)

M2_bal_resdh <- simulateResiduals(M2_bal)
plot(M2_bal_resdh)

M2_bal_res <- residuals(M2_bal)

# pattern in residuals is a little better with time variance accounted for

# try linearizing timediff

balanus_exp <- balanus %>% 
  mutate(exp_time = exp(-1/timediff))

balanus_exp$colour <- as.factor(balanus_exp$colour)

M3_bal <- glmmTMB(count ~ exp_time*colour + (1|block), 
                  data = balanus_exp, family = nbinom2,
                  dispformula = ~exp_time)
AIC(M1_bal, M3_bal)

M3_bal_resdh <- simulateResiduals(M3_bal)
plot(M3_bal_resdh)

M3_bal_res <- residuals(M3_bal)

plot(M3_bal_res ~ balanus_exp$timediff)
plot(M3_bal_res ~ balanus_exp$colour)

summary(M3_bal)
library(car)
Anova(M3_bal)

## now for ulothrix

hist(ulothrix$percent_cover, breaks = 100)

ulothrix <- ulothrix %>% 
  mutate(prop_cover = percent_cover/100)

# sort of a zero-inflated beta situation here gamlss?
library(gamlss)
M1_ulo <- gamlss(prop_cover ~ timediff*colour,
                 data = ulothrix, family = BEINF())
summary(M1_ulo)
plot(M1_ulo)
glimpse(ulothrix)
M1_ulo_res <- residuals(M1_ulo)
plot(M1_ulo_res ~ ulothrix$timediff)
plot(M1_ulo_res ~ as.factor(ulothrix$colour))

M6_ulo <- gamlss(prop_cover ~ poly(timediff,2)*colour + re(random = ~timediff|block),
                 data = ulothrix, family = BEINF())

summary(M6_ulo)

M6_ulo_res <- residuals(M6_ulo)
plot(M6_ulo_res ~ ulothrix$timediff)

plot(M6_ulo)

# want an interaction between timediff and colour, could use beta in glmmTMB

a <- length(ulothrix$prop_cover)
ulothrix$pcp <- (ulothrix$prop_cover*(a-1)+0.5)/a

M2_ulo <- glmmTMB(pcp ~ timediff*colour, 
                  data = ulothrix, family = beta_family())

M2_ulo_resdh <- simulateResiduals(M2_ulo)
plot(M2_ulo_resdh)

M2_ulo_res <- residuals(M2_ulo)
plot(M2_ulo_res ~ ulothrix$timediff)
plot(M2_ulo_res ~ as.factor(ulothrix$colour))

# put in block term

M3_ulo <- glmmTMB(pcp~ timediff*colour + (1|block), data = ulothrix,
                  family = beta_family())

summary(M3_ulo)
AIC(M2_ulo, M3_ulo)

M3_ulo_res_dh <- simulateResiduals(M3_ulo)
plot(M3_ulo_res_dh)

M3_ulo_res <- residuals(M3_ulo)
plot(M3_ulo_res ~ ulothrix$timediff)
plot(M3_ulo_res ~ as.factor(ulothrix$colour))

# still a pattern in the residuals. maybe try disp formula?

M4_ulo <- glmmTMB(pcp~ timediff*colour + (1|block), data = ulothrix,
                  family = beta_family(), dispformula = ~timediff)

AIC(M3_ulo, M4_ulo) 
anova(M3_ulo, M4_ulo) # M4 is better

M4_ulo_res_dh <- simulateResiduals(M4_ulo)
plot(M4_ulo_res_dh) # makes really weird patterns though.

# best model is M4 so far, so M3 is best
# correct for time?

acf(M3_ulo_res)
#Yes - positive autocorrelation throughout.

ulothrix$b_no <- paste(ulothrix$block, ulothrix$number, sep = "_")

M5_ulo <- glmmTMB(pcp~ timediff*colour + (1|block) +
                  diag(timediff|number/block), data = ulothrix,
                  family = beta_family())

M5_ulo_res <- residuals(M5_ulo)
acf(M5_ulo_res)
plot(M5_ulo_res ~ ulothrix$timediff)

library(lmtest)

M5_ulo_resdh <- simulateResiduals(M5_ulo)
plot(M5_ulo_resdh)

AIC(M3_ulo, M5_ulo)
anova(M3_ulo, M5_ulo)

plot(ulothrix$prop_cover ~ exp(-1/(ulothrix$timediff)))
plot(ulothrix$prop_cover ~ log(ulothrix$timediff))

ulothrix$exp_time <- exp(-1/ulothrix$timediff)

M7_ulo <-  glmmTMB(pcp~ exp_time*colour + (1|block) +
                     diag(exp_time|number/block), data = ulothrix,
                   family = beta_family())
AIC(M3_ulo, M7_ulo)     
anova(M3_ulo, M7_ulo)

M7_ulo_resdh <- simulateResiduals(M7_ulo)
plot(M7_ulo_resdh)

M7_ulo_res <- residuals(M7_ulo)
plot(M7_ulo_res ~ ulothrix$exp_time)
plot(M7_ulo_res ~ as.factor(ulothrix$colour))
# better by AIC, anova, but still a linear trend in the residuals. leave it here for now?
acf(M7_ulo_res)
# autocorrelation still present, but massively attenuated.

## looked at model predictions, and gamlss M6 is way better at fitting the data.
# need to talk to Devin about whether it's ok to have interactions with smoothed terms

ulothrix$pred <- predict(M6_ulo, type = "response")*100

AIC(M1_ulo, M6_ulo)
