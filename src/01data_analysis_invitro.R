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
biofilm$logOD <- log2(biofilm$OD)
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


####  FIGURE 2 ####
corBiofilmOD %>% 
      filter(category != "biofilm") %>% 
      ggplot(., aes(ABTCAA, LBNaCl))+
      geom_point(aes(color = category), alpha = 0.75, stroke = 0, size = 2)+
      geom_text_repel(data = labels_strains, aes(label=strain), color = 'black',
                      force_pull   = 0.5, # do not pull toward data points
                      nudge_y      = 6,
                      direction    = "x",
                      angle        = 90,
                      hjust        = 0,
                      segment.color = "grey80",
                      segment.size = 0.2)+
      geom_abline(aes(slope=1, intercept=0))+
      geom_point(data = labels_strains, aes(fill=category), pch = 21, color = 'black', size = 2)+
      theme_rs()+
      theme(aspect.ratio = 0.75)+
      guides(color = guide_legend(title = "Source", override.aes = list(alpha = 1, size = 3)), fill = "none")+
      scale_fill_manual(values = c("#F8766D", '#00BFC4'))+
      scale_color_manual(values = c("#F8766D", '#7CAE00', '#00BFC4'), labels = c("Fresh produce", "Soil", "Water"))+
      scale_y_continuous(limits=c(0,10), breaks = seq(0,10,2))+
      labs(y = "Biofilm in LB-NaCl(0)\n(OD600)", 
           x = "Biofilm in ABTCAA\n(OD600)")
ggsave(here("output", "fig2.pdf"), width = 3.5, dpi = 300)



####  FIGURE 3 ####
#   Annotation
pal <- palette()[c(2:4,6)]
names(pal) <- c("none/weak", "moderate", "strong", "extreme")

lb.pal <- palette_biofilm[match(abr$biofilm.in.LB.NaCl.0., names(palette_biofilm))]
abt.pal <- palette_biofilm[match(abr$biofilm.in.ABTCAA, names(palette_biofilm))]
col_fun = colorRamp2(c(0, 0.25, 0.5, 0.75, 1), brewer.pal(name = "YlGnBu", n = 9)[c(1,2,4,7,9)])

annotation_biofilm <- rowAnnotation(
      MAR = abr_long$mar,
      LBNaCl = abr$biofilm.in.LB.NaCl.0.,
      ABTCAA = abr$biofilm.in.ABTCAA,
      col = list(MAR = col_fun,
                 ABTCAA = abt.pal,
                 LBNaCl = lb.pal),
      gp = gpar(col = "black"),
      simple_anno_size = unit(4, "mm"),
      border = TRUE,
      annotation_name_gp= gpar(fontsize = 12))

heatmap_abr <- abr_long %>% 
      column_to_rownames(var = "strain") %>% 
      dplyr::select("AM10":"CTX") %>%
      as.matrix %>% 
      Heatmap(
            right_annotation = annotation_biofilm,
            cluster_rows = FALSE,
            border = TRUE,
            rect_gp = gpar(col = "gray", lwd = 1),
            column_dend_height = unit(0, "cm"),
            col = c("grey90", "black"),
            heatmap_legend_param = list(title = "ABR",
                                        labels = c("Yes", "No")))

pdf(here("output", "fig3.pdf"), width = 7.2, height = 3)
heatmap_abr
dev.off()

####  FIGURE 4 ####



corr_cfu_biofilm %>% 
      mutate(type = factor(type, levels = names(palette_biofilm))) %>% 
      ggplot(aes(x = logcfu, y = logcopies, group = type, fill = type))+
      facet_wrap(~ type, ncol = 2, labeller = labeller(type = lab_biofilm))+
      geom_point(pch = 21, alpha = 0.6, size = 1.5)+
      geom_abline(aes(intercept = intercept, slope = m), data = results)+
      geom_abline(slope = 1, intercept = 0, linetype = "dashed")+
      theme_rs()+
      theme(aspect.ratio = 1)+
      scale_y_continuous(name = "Bacterial density\n[log10 CFU gFW-1]", limits = c(2,11), breaks = seq(2,11,2))+
      scale_x_continuous(name = "Gene copy number \n[log10 yccT copies gFW-1]", limits = c(1,10), breaks = seq(2,10,2))+
      scale_fill_manual(labels = lab_biofilm, values = palette_biofilm)+
      guides(fill = "none")


cfu_biofilm %>% 
      ggplot(aes(logcfu, logcopies, fill = type))+
      facet_wrap(~strain)+
      geom_point(pch = 21, alpha = 0.6, size = 1.5)+
      geom_abline(slope = 1)+
      theme_rs()+
      theme(aspect.ratio = 1)+
      scale_y_continuous(name = "Bacterial density\n[log10 CFU gFW-1]", limits = c(2,11), breaks = seq(2,11,2))+
      scale_x_continuous(name = "qPCR data\n[log10 yccT copies gFW-1]", limits = c(1,10), breaks = seq(2,10,2))+
      guides(fill = guide_legend(title = "Strain", override.aes = list(size = 2, alpha = 0.8)))
ggsave(here("output", "fig4.pdf"), width = 3.5, dpi = 300)

####  FIGURE S1 ####
#     Biofilm OD values for every strains
vec = biofilm %>% 
      filter(medium =="ABTCAA") %>% 
      arrange(-OD) %>% 
      dplyr::select(strain) %>% 
      t %>% 
      as.vector()
biofilm2$strain = factor(biofilm2$strain, levels = vec)

fS1a <- corBiofilmOD %>% 
      arrange(-LBNaCl) %>% 
      slice(1:57) %>% 
      ggplot()+
      geom_segment(aes(x = log2(ABTCAA), xend = log2(LBNaCl), y = reorder(strain, LBNaCl), yend = strain), color="black") +
      geom_point(aes(x = log2(ABTCAA), y = strain), color = rgb(0.2, 0.7, 0.1, 0.5), size = 3, alpha = .9) +
      geom_point(aes(x = log2(LBNaCl), y = strain), color = rgb(0.7, 0.2, 0.1, 0.5), size = 3, alpha = .9)+
      scale_x_continuous(limits = c(-4,4))+
      labs(x = "Biofilm (log2 OD600)", y = "Strain")+
      theme_rs()
fS1b <- corBiofilmOD %>% 
      arrange(-LBNaCl) %>% 
      slice(58:115) %>% 
      ggplot()+
      geom_segment(aes(x = log2(ABTCAA), xend = log2(LBNaCl), y = reorder(strain, LBNaCl), yend = strain), color="black") +
      geom_point(aes(x = log2(ABTCAA), y = strain), color = rgb(0.2, 0.7, 0.1, 0.5), size = 3, alpha = .9) +
      geom_point(aes(x = log2(LBNaCl), y = strain), color = rgb(0.7, 0.2, 0.1, 0.5), size = 3, alpha = .9)+
      scale_x_continuous(limits = c(-4,4))+
      labs(x = "Biofilm (log2 OD600)", y = "")+
      theme_rs()
fS1c <- corBiofilmOD %>% 
      arrange(-LBNaCl) %>% 
      slice(116:174) %>% 
      ggplot()+
      geom_segment(aes(x = log2(ABTCAA), xend = log2(LBNaCl), y = reorder(strain, LBNaCl), yend = strain), color="black") +
      geom_point(aes(x = log2(ABTCAA), y = strain), color = rgb(0.2, 0.7, 0.1, 0.5), size = 3, alpha = .9) +
      geom_point(aes(x = log2(LBNaCl), y = strain), color = rgb(0.7, 0.2, 0.1, 0.5), size = 3, alpha = .9)+
      scale_x_continuous(limits = c(-4,4))+
      labs(x = "Biofilm (log2 OD600)", y = "")+
      theme_rs()

fS1a + fS1b + fS1c
ggsave(here("output", "figS1.pdf"), width = 7.2, height = 8, dpi = 300)

#### Data analysis ####
