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

#     Dependencies
source(here("src", "theme_rs.R"))
source(here("src", "01data_analysis_invitro.R"))

#     Data
cfu_biofilm <- cfu_biofilm %>% filter(dpi != 0)
cfu_biofilm$strain <- factor(cfu_biofilm$strain, levels = c("H2", "B456", "C1", "C13", "B471", "B545", "C30", "C160", "B368", "B466", "C15"))

#     Summary
cfu_biofilm_summaryAll = cfu_biofilm %>% 
      group_by(dpi, type, strain, OD) %>% 
      summarise(mean_cfu = mean(logcfu),
                sd_cfu = sd(logcfu),
                cv_cfu = 100*sd_cfu/mean_cfu,
                mean_copies = mean(logcopies),
                sd_copies = sd(logcopies),
                cv_copies = 100*sd_copies/mean_copies,
                n = length(logcfu),
                .groups = "drop")
cfu_biofilm_summaryAll$strain = factor(cfu_biofilm_summaryAll$strain, # for ABTCAA
                                       levels = c("H2", "B456", "C1", "C13", "B471", "B545", "C30", "C160", "B368", "B466", "C15"))

cfu_biofilm_summary_type = cfu_biofilm %>% 
      group_by(dpi, type) %>% 
      summarise(mean_cfu = mean(logcfu),
                sd_cfu = sd(logcfu),
                cv_cfu = 100*sd_cfu/mean_cfu,
                mean_copies = mean(logcopies),
                sd_copies = sd(logcopies),
                cv_copies = 100*sd_copies/mean_copies,
                n = length(logcfu))

cfu_biofilm_summary_exp = cfu_biofilm %>% 
      group_by(dpi, exp) %>% 
      summarise(mean_cfu = mean(logcfu),
                sd_cfu = sd(logcfu),
                cv_cfu = 100*sd_cfu/mean_cfu,
                mean_copies = mean(logcopies),
                sd_copies = sd(logcopies),
                cv_copies = 100*sd_copies/mean_copies,
                n = length(logcfu)) %>% na.omit
lab_exp = cfu_biofilm_summary_exp %>% filter(dpi == "21")
lab_exp$exp = factor(lab_exp$exp)

##    GLS
#     Variance structures
vf1Id = varIdent(form = ~1|dpi)
vf2Id = varIdent(form = ~1|type)
vf3Id = varIdent(form = ~1|strain)
vf4Id = varIdent(form = ~1|exp)
vf5Id = varIdent(form = ~1|dpi*type)
vf6Id = varIdent(form = ~1|dpi*strain)
vf7Id = varIdent(form = ~1|dpi*exp)

####  TYPE OF BIOFILM   ####
#     GLS
M1    <- gls(logcopies ~ type * as.factor(dpi), data = cfu_biofilm)
M1id  <- gls(logcopies ~ type * as.factor(dpi), weights = vf1Id, data = cfu_biofilm)
M2id  <- gls(logcopies ~ type * as.factor(dpi), weights = vf2Id, data = cfu_biofilm)
M3id  <- gls(logcopies ~ type * as.factor(dpi), weights = vf3Id, data = cfu_biofilm)
M4id  <- gls(logcopies ~ type * as.factor(dpi), weights = vf4Id, data = cfu_biofilm)
M5id  <- gls(logcopies ~ type * as.factor(dpi), weights = vf5Id, data = cfu_biofilm)
M6id  <- gls(logcopies ~ type * as.factor(dpi), weights = vf6Id, data = cfu_biofilm)
M7id  <- gls(logcopies ~ type * as.factor(dpi), weights = vf7Id, data = cfu_biofilm)

#     MODEL SELECTION
anova(M1, M1id, M2id, M3id, M4id, M5id, M6id, M7id)
anova(M1, M7id) # best model

#     CHECK RESIDUALS AND VARIANCE
E2 <- resid(M7id, type = "normalized")
coplot(E2 ~ type | dpi, ylab = "Ordinary residuals", data = cfu_biofilm)
qqnorm(E2)

#     ANOVA ON BEST MODEL
anova(M7id)

#     MEAN COMPARISON
em_M7id = emmeans(M7id, ~ type * dpi, data = cfu_biofilm)
contrast(em_M7id, 'pairwise', type = 'response', adjust = "bonferroni") %>% tidy %>% write.csv(., 'output/data/pairwise_copies_M7id_type.csv')

#     CREATE DATA FRAME
df_M7id = cld(em_M7id) %>% 
      tidy() %>% 
      data.frame() %>% 
      arrange(dpi,type) %>% 
      mutate_if(is.character, str_trim)
df_M7id$type <- factor(df_M7id$type, levels = c("none/weak", "moderate", "strong", "extreme"))



####  STRAIN   ####
#     GLS
M1s   <- gls(logcopies ~ strain * as.factor(dpi), data = cfu_biofilm)
M1sid <- gls(logcopies ~ strain * as.factor(dpi), weights = vf1Id, data = cfu_biofilm)
M2sid <- gls(logcopies ~ strain * as.factor(dpi), weights = vf2Id, data = cfu_biofilm)
M3sid <- gls(logcopies ~ strain * as.factor(dpi), weights = vf3Id, data = cfu_biofilm)
M4sid <- gls(logcopies ~ strain * as.factor(dpi), weights = vf4Id, data = cfu_biofilm)
M5sid <- gls(logcopies ~ strain * as.factor(dpi), weights = vf5Id, data = cfu_biofilm)
M6sid <- gls(logcopies ~ strain * as.factor(dpi), weights = vf6Id, data = cfu_biofilm)
M7sid <- gls(logcopies ~ strain * as.factor(dpi), weights = vf7Id, data = cfu_biofilm)

#     MODEL SELECTION
anova(M1s, M1sid, M2sid, M3sid, M4sid, M5sid, M6sid, M7sid)
anova(M1s, M7sid)

#     CHECK RESIDUALS AND VARIANCE
E2 <- resid(M7sid, type = "normalized")
coplot(E2 ~ strain | dpi, ylab = "Ordinary residuals", data = cfu_biofilm)
qqnorm(E2)
summary(M7sid)

#     ANOVA ON BEST MODEL
anova(M7sid)

#     MEAN COMPARISON
em_M7sid_strain = emmeans(M7sid, ~ strain * as.factor(dpi), data = cfu_biofilm)
contrast(em_M7sid_strain, 'pairwise', type = 'response', adjust = "bonferroni") %>% tidy %>% write.csv(., 'output/data/pairwise_copies_M7sid_strain.csv')

#     CREATE DATA FRAME
df_M7sid = cld(em_M7sid_strain) %>% 
      tidy() %>% 
      left_join(., unique(cfu_biofilm[,c(1,3,6)]), by = c("strain", "dpi")) %>% 
      data.frame() %>% 
      mutate_if(is.character, str_trim)
df_M7sid$type = factor(df_M7sid$type, levels = c("none/weak", "moderate", "strong", "extreme"))
df_M7sid$dpi = factor(df_M7sid$dpi, levels = c("3", "7", "14", "21"))
df_M7sid$strain = factor(df_M7sid$strain, levels = c("H2", "B456", "C1", "C13", "B471", "B545", "C30", "C160", "B368", "B466", "C15"))


####  FIGURE 5 #####
f5a <- df_M7id %>% 
      ggplot(aes(x = dpi, y = estimate, fill = type))+
      facet_wrap(~type, ncol = 4)+
      geom_jitter(data = cfu_biofilm, aes(x = dpi, y = logcopies, color = type),
                  width = 0.9, alpha = 0.8, size = 2, stroke = 0)+
      geom_point(size = 2, stroke = 0.5, fill = "black", color = "grey", pch=21, position = position_dodge(width = 2))+
      geom_text(aes(label = .group, y = 11.5), position = position_dodge(width = 2))+
      geom_line(alpha = 0.5, linetype = "dashed")+
      theme_rs()+
      theme(aspect.ratio = 1)+
      guides(color = guide_legend(title = "Biofilm type", override.aes = list(size = 4, alpha = 1)), fill = "none")+
      scale_y_continuous(name = "Log10 copies gFW-1", limits = c(2,12), breaks = seq(2,12,2))+
      scale_x_continuous(name = "Time [dpi]", limits = c(0, 25), breaks = c(3,7,14,21))+
      scale_color_manual(values = palette_biofilm, labels = lab_biofilm)

f5b <- df_M7sid %>% 
      ggplot(aes(x = strain, y = estimate, fill = type))+
      facet_wrap(~as.factor(dpi), ncol = 4)+
      geom_jitter(data = cfu_biofilm, aes(x = strain, y = logcopies, color = type),
                  width = 0.2, alpha = 0.8, size = 2, stroke = 0)+
      geom_point(size = 2, stroke = 0.5, fill = "black", color = "grey", pch=21, position = position_dodge(width = 2))+
      geom_text(aes(label = .group, y=11.5), position = position_dodge(width = 0.9))+
      theme_rs()+
      theme(aspect.ratio = 1,
            axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1))+
      guides(color = guide_legend(title = "Biofilm type", override.aes = list(size = 4, alpha = 1)),
             fill = "none")+
      scale_y_continuous(name = "Log10 copies gFW-1", limits = c(2,12), breaks = seq(2,12,2))+
      scale_x_discrete(name = "Strain")+
      scale_color_manual(values = palette_biofilm, labels = lab_biofilm)

f5a/f5b+
      plot_annotation(tag_levels = "A")+
      plot_layout(guides = "collect")
ggsave(here("output", "fig5.pdf"), width = 7.2, dpi = 300)


####  FIGURE S2   ####
fs2a <- cfu_biofilm_summary_type %>% 
      ggplot(aes(x = dpi, y = cv_copies, fill = type, color = type))+
      geom_line(linewidth = 1, alpha = 0.8)+
      geom_point(pch = 21, size = 2.5, stroke = 0.5, color = "black")+
      theme_rs()+
      theme(aspect.ratio = 0.75)+
      scale_shape(solid = TRUE)+
      scale_y_continuous(name = "Coef. of Variation [%]", limits = c(0,35), expand = c(0,0))+
      scale_x_continuous(name = "Time [dpi]", limits = c(-2,23), expand = c(0,0))+
      scale_color_manual(values = palette_biofilm, labels = lab_biofilm)+
      scale_fill_manual(values = palette_biofilm, labels = lab_biofilm)+
      guides(color = guide_legend(title = "Biofilm type", override.aes = list(size = 4, alpha = 1)), fill = "none")

fs2b <- cfu_biofilm_summaryAll %>% 
      ggplot(aes(x = dpi, y = cv_copies, fill = type, color = type, group = strain))+
      geom_line(aes(group = interaction(type,strain)), linewidth = 1, alpha = 0.8)+
      geom_point(pch = 21, size = 2.5, stroke = 0.5, color = "black", alpha = 1)+
      theme_rs()+
      theme(aspect.ratio = 0.75)+
      scale_shape(solid = TRUE)+
      scale_y_continuous(name = "Coef. of Variation [%]", limits = c(0,35), expand = c(0,0))+
      scale_x_continuous(name = "Time [dpi]", limits = c(-2,23), expand = c(0,0))+
      scale_color_manual(values = palette_biofilm, labels = lab_biofilm)+
      scale_fill_manual(values = palette_biofilm, labels = lab_biofilm)+
      guides(color = guide_legend(title = "Biofilm type", override.aes = list(size = 4, alpha = 1)), fill = "none")

fs2c <- cfu_biofilm_summary_exp %>% 
      ggplot(., aes(x = dpi, y = cv_copies, group = as.factor(exp)))+
      geom_line(linewidth = 0.7)+
      geom_point(size = 2.5, stroke = 0 , color = "black", alpha = 1)+
      geom_text_repel(data = lab_exp, aes(label = exp),
                      color = "red",
                      force             = 0.1,
                      nudge_x           = 1,
                      direction         = "y",
                      hjust             = -1,
                      segment.size      = 0.5)+
      theme_rs()+
      theme(aspect.ratio = 0.75)+
      scale_shape(solid = TRUE)+
      scale_y_continuous(name = "Coef. of Variation [%]", limits = c(0,38), expand = c(0,0), breaks = seq(0,30,10))+
      scale_x_continuous(name = "Time [dpi]", limits = c(-2,24), expand = c(0,0))

fs2a + fs2b + fs2c +
      plot_annotation(tag_levels = "A")+
      plot_layout(guides = "collect")
ggsave(here("output", "figs2.pdf"), width = 7.2, dpi = 300)
