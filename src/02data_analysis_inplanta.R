#!/usr/bin/env Rscript

###   Analysis in planta
#     Libraries
library(tidyverse)
library(nlme)
library(broom)
library(multcomp)
library(multcompView)
library(forcats)
library(patchwork)
library(ggrepel)
library(emmeans)

#     Dependencies
source(here("src", "00data_clean.R"))

##    Correlation between CFU and qPCR data
#     Extract the estimated coefficient and its standard error
model_corr <- lmer(logcopies ~ logcfu + (logcfu|type), data = corr_cfu_biofilm)
summary_model <- summary(model_corr)

#     Extract the random effects for each group
random_effects <- ranef(model_corr)$type
random_effects

#     Extract the fixed effects
fixed_effects <- fixef(model_corr)
fixed_effects

#     Extract the random slopes and their standard errors
var_cor <- as.data.frame(VarCorr(model_corr))
std_errs <- sqrt(var_cor[var_cor$grp == "type" & var_cor$var1 == "logcfu", "vcov"])

#     Perform hypothesis tests
results <- data.frame(type = rownames(random_effects), 
                      slope = random_effects[, "logcfu"], 
                      std_err = std_errs) %>% 
      mutate(
            intercept = fixed_effects["(Intercept)"] + random_effects[,"(Intercept)"],
            m = slope + fixed_effects["logcfu"],
            z_value = (m - 1)/std_err,
            p_value = 2 * pnorm(-abs(z_value)))

# Output the results
results
results$type <- factor(results$type, levels = c('none/weak', 'moderate', 'strong', 'extreme'))

# Correlations
corr_cfu_biofilm %>% 
      group_by(type) %>% 
      tally

correlation_all <- corr_cfu_biofilm %>% 
      cor_test(logcfu, logcopies, alternative = "greater", method = "pearson")

correlation_qpcr_cfu <- corr_cfu_biofilm %>% 
      group_by(type) %>% 
      cor_test(logcfu, logcopies, alternative = "greater", method = "pearson")

##    GLS
#     Variance structures
vf1Id = varIdent(form = ~1|dpi)
vf2Id = varIdent(form = ~1|type)
vf3Id = varIdent(form = ~1|strain)
vf4Id = varIdent(form = ~1|exp)

####  TYPE OF BIOFILM   ####
#     GLS
M1    <- gls(logcopies ~ type * as.factor(dpi), data = corr_cfu_biofilm)
M1id  <- gls(logcopies ~ type * as.factor(dpi), weights = vf1Id, data = corr_cfu_biofilm)
M2id  <- gls(logcopies ~ type * as.factor(dpi), weights = vf2Id, data = corr_cfu_biofilm)
M3id  <- gls(logcopies ~ type * as.factor(dpi), weights = vf3Id, data = corr_cfu_biofilm)
M4id  <- gls(logcopies ~ type * as.factor(dpi), weights = vf4Id, data = corr_cfu_biofilm)

#     MODEL SELECTION
anova(M1, M1id, M2id, M3id, M4id)
anova(M1, M4id) # best model

#     CHECK RESIDUALS AND VARIANCE
E2 <- resid(M4id, type = "normalized")
coplot(E2 ~ type | dpi, ylab = "Ordinary residuals", data = corr_cfu_biofilm)
qqnorm(E2)

#     ANOVA ON BEST MODEL
anova(M4id)

#     MEAN COMPARISON
em_M4id = emmeans(M4id, ~ type * dpi, data = corr_cfu_biofilm)
contrast(em_M4id, 'pairwise', type = 'response', adjust = "BH") %>% 
      tidy %>% 
      write.csv(., 'output/data/pairwise_copies_M4id_type.csv')

#     CREATE DATA FRAME
df_M4id = cld(em_M4id, Letters = c("a", "b", "c", "d")) %>% 
      tidy() %>% 
      data.frame() %>% 
      arrange(dpi,type) %>% 
      mutate_if(is.character, str_trim)
df_M4id$type <- factor(df_M4id$type, 
                       levels = c("none/weak", "moderate", "strong", "extreme"))


####  STRAIN   ####
#     GLS
M1s   <- gls(logcopies ~ strain * as.factor(dpi), data = corr_cfu_biofilm)
M1sid <- gls(logcopies ~ strain * as.factor(dpi), weights = vf1Id, data = corr_cfu_biofilm)
M2sid <- gls(logcopies ~ strain * as.factor(dpi), weights = vf2Id, data = corr_cfu_biofilm)
M3sid <- gls(logcopies ~ strain * as.factor(dpi), weights = vf3Id, data = corr_cfu_biofilm)
M4sid <- gls(logcopies ~ strain * as.factor(dpi), weights = vf4Id, data = corr_cfu_biofilm)

#     MODEL SELECTION
anova(M1s, M1sid, M2sid, M3sid, M4sid)
anova(M1s, M4sid)

#     CHECK RESIDUALS AND VARIANCE
E2 <- resid(M4sid, type = "normalized")
coplot(E2 ~ strain | dpi, ylab = "Ordinary residuals", data = corr_cfu_biofilm)
qqnorm(E2)
summary(M4sid)

#     ANOVA ON BEST MODEL
anova(M4sid)

#     MEAN COMPARISON
em_M4sid_strain = emmeans(M4sid, ~ strain * as.factor(dpi), data = corr_cfu_biofilm)
contrast(em_M4sid_strain, 'pairwise', type = 'response', adjust = "bonferroni") %>% 
      tidy %>% 
      write.csv(., 'output/data/pairwise_copies_M4sid_strain.csv')

#     CREATE DATA FRAME
df_M4sid = cld(em_M4sid_strain, Letters = c("a", "b", "c", "d")) %>% 
      tidy() %>% 
      left_join(., unique(cfu_biofilm[,c(1,3,6)]), by = c("strain", "dpi")) %>% 
      data.frame() %>% 
      mutate_if(is.character, str_trim)
df_M4sid$type = factor(df_M4sid$type, levels = c("none/weak", "moderate", "strong", "extreme"))
df_M4sid$dpi = factor(df_M4sid$dpi, levels = c("0", "3", "7", "14", "21"))
df_M4sid$strain = factor(df_M4sid$strain, levels = c("H2", "B456", "C1", "C13", "B471", "B545", "C30", "C160", "B368", "B466", "C15"))

### MAR vs PERSISTENCE
mod_mar <- lm(logcopies ~ mar*type*dpi, data = biofilm_plant_mar[biofilm_plant_mar$type != 'none/weak',])
summary(mod_mar)
resid(mod_mar) %>% qqnorm

anova(mod_mar)
aov(logcopies ~ mar*type*dpi, data = biofilm_plant_mar[biofilm_plant_mar$type != 'none/weak',]) %>% 
      tidy() %>% 
      mutate(partialR2 = 100*sumsq/sum(sumsq))