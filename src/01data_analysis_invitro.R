#!/usr/bin/env Rscript

##    In vitro Biofilm Data
#     Dependencies
library(tidyverse)
library(here)
library(rstatix)
library(lme4)
library(emmeans)

#     Dependencies
source(here("src", "00data_clean.R"))

##    Data analysis
biofilm_filter <- biofilm %>% filter(category!='biofilm')
biofilm_abtcaa <- biofilm_filter %>% filter(medium == 'ABTCAA')
biofilm_lbnonacl <- biofilm_filter %>% filter(medium == 'LBNaCl')

##    Paired t-test NaCl vs ABTCAA
ttest_result_nacl_abtcaa <- t.test(biofilm$logOD[biofilm$medium == 'LBNaCl'],
       biofilm$logOD[biofilm$medium == 'ABTCAA'],
       paired = TRUE)
#     Output
ttest_result_nacl_abtcaa$statistic
ttest_result_nacl_abtcaa$parameter
ttest_result_nacl_abtcaa$estimate
ttest_result_nacl_abtcaa$p.value

##    Mixed effect model for type and medium in biofilm formation
mod_type_medium <- lmer(logOD ~ medium * type + (1|strain), data = biofilm)
summary(mod_type_medium)
resid(mod_type_medium) %>% qqnorm

Anova(mod_type_medium, type = "III")

em_type_medium <- emmeans(mod_type_medium, specs = ~ medium | type)
con_type_medium <- contrast(em_type_medium, 'pairwise', adjust = 'BH')


##    Mixed effect model for category and medium in biofilm formation
mod_cat_medium <- lmer(logOD ~ medium * category + (1|strain), data = biofilm_filter)
summary(mod_cat_medium)
resid(mod_cat_medium) %>% qqnorm

Anova(mod_cat_medium, type = "III")

em_cat_medium1 <- emmeans(mod_cat_medium, specs = ~ medium | category)
con_cat_medium1 <- contrast(em_cat_medium1, 'pairwise', adjust = 'BH')

em_cat_medium2 <- emmeans(mod_cat_medium, specs = ~ category | medium)
con_cat_medium2 <- contrast(em_cat_medium2, 'pairwise', adjust = 'BH')


##    Mixed effect model for phylogroup and medium in biofilm formation
mod_phylo_medium <- lmer(logOD ~ medium * phylogroup + (1|strain), data = biofilm[biofilm$phylogroup!="not determined",])
summary(mod_phylo_medium)
resid(mod_phylo_medium) %>% qqnorm

Anova(mod_phylo_medium, type = "III")

em_phylo_medium <- emmeans(mod_phylo_medium, specs = ~ phylogroup | medium)
con_phylo_medium <- contrast(em_phylo_medium, 'pairwise', adjust = 'BH')


##    Biofilm formation in ABTCAA
#     Assumptions
#     Normality
hist(biofilm_abtcaa$logOD, breaks = 20)
shapiro.test(biofilm_abtcaa$logOD)

#     Homoscedasticity
bartlett.test(biofilm_abtcaa$logOD, g = biofilm_abtcaa$type)

#     Kruskal-Wallis
biofilm_abtcaa %>% kruskal_test(logOD ~ type)
biofilm_abtcaa %>% kruskal_test(logOD ~ category)

#     Dunn test
biofilm_abtcaa %>% dunn_test(logOD ~ type, p.adjust.method = 'BH')


##    Biofilm formation in LBnoNaCl
#     Assumptions
#     Normality
hist(biofilm_lbnonacl$logOD, breaks = 20)
shapiro.test(biofilm_lbnonacl$logOD)

#     Homoscedasticity
bartlett.test(biofilm_lbnonacl$logOD, g = biofilm_lbnonacl$type)

#     Kruskal-Wallis
biofilm_lbnonacl %>% kruskal_test(logOD ~ type)
biofilm_lbnonacl %>% kruskal_test(logOD ~ category)

#     Dunn test
biofilm_lbnonacl %>% dunn_test(logOD ~ type, p.adjust.method = 'BH')
biofilm_lbnonacl %>% dunn_test(logOD ~ category, p.adjust.method = 'BH')

#     MAR vs biofilm type
#     Linear model ABTCAA medium
medium <- c('ABTCAA', 'LBNaCl')
list_mar_mod <- list()
df_mod_mar <- data.frame()

for (i in 1:length(medium)) {
      temp_mod <- lm(logOD ~ mar, data = biofilm_mar[biofilm_mar$medium==medium[i],])
      list_mar_mod[[i]] <- temp_mod
      temp_df <- temp_mod %>% 
            summary %>% 
            tidy %>% 
            mutate(medium = medium[i]) %>% 
            filter(term == "mar")
      df_mod_mar <- rbind(df_mod_mar, temp_df)
      
}

#     Check residuals
list_mar_mod[[1]] %>% resid %>% qqnorm
list_mar_mod[[2]] %>% resid %>% qqnorm

#     Correlation
df_cor_mar <-biofilm_mar %>% 
      group_by(medium) %>% 
      cor_test(logOD, mar, method = 'pearson')
cor.test(biofilm_mar$logOD[biofilm_mar$medium=="ABTCAA"], biofilm_mar$mar[biofilm_mar$medium=="ABTCAA"])
cor.test(biofilm_mar$logOD[biofilm_mar$medium=="LBNaCl"], biofilm_mar$mar[biofilm_mar$medium=="LBNaCl"])

#     Combine data frames
df_mar <- left_join(df_mod_mar, df_cor_mar, by = 'medium', suffix = c('.mod', '.cor'))
