#!/usr/bin/env Rscript

##    In vitro Biofilm Data
#     Dependencies
library(tidyverse)
library(here)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)

#     Dependencies
source(here("src", "00data_clean.R"))
source(here("src", "theme_rs.R"))

#     Labels
lab_biofilm <- c("None/Weak", "Moderate", "Strong", "Extreme")
lab_medium <- c("ABTCAA", "LB-NaCl(0)")

#     Color palette biofilm type
palette_biofilm <- brewer.pal(name = "BuPu", 9)[c(3,5,7,9)]
names(palette_biofilm) <- c("none/weak", "moderate", "strong", "extreme")

####  FIGURE 1   ####
f1.a <- biofilm_strength %>% 
      ggplot(aes(x = medium, y = log2(OD)))+
      geom_violin(size = 0.5, trim = TRUE, scale = 'width', fill = 'grey90', adjust = 1.2)+
      geom_boxplot(size = 0.5, outlier.alpha = 0, width=0.2)+
      geom_line(aes(group = strain, color = is_LBhigher), alpha = 0.3)+
      geom_point(aes(fill = type), pch = 21, size = 1.5, alpha = 0.5, stroke = 0.2)+
      theme_rs()+
      theme(aspect.ratio = 0.75)+
      labs(y = "Biofilm (log2 OD)", x = "Medium")+
      scale_fill_manual(values = palette_biofilm, labels = lab_biofilm)+
      scale_color_manual(values = c("#009900", "#cc0000", "black"))+
      scale_y_continuous(limits = c(-4,4), expand = c(0,0))+
      scale_x_discrete(labels = lab_medium)+
      stat_compare_means(method = 't.test', aes(label=..p.signif..), size = 6, label.x = 1.5, label.y = 3, comparisons = list(c("ABTCAA","LBNaCl")))+
      guides(color = "none", 
             fill = guide_legend(title = "Biofilm type", override.aes = list(size = 4, alpha = 1)))

f1.b <- biofilm %>% 
      ggplot(aes(x = type, log2(OD), fill = medium, group = interaction(medium,type)))+
      geom_violin(size=0.25, trim = TRUE, scale = 'width', adjust = 1.2)+
      geom_jitter(aes(color = type), alpha = 0.5, size = 1.5, stroke = 0, position = position_jitterdodge(jitter.width = 0.8, dodge.width = 0.9))+
      geom_boxplot(fill = "white", size = 0.25, outlier.alpha = 0, width = 0.1, position = position_dodge(width=0.9))+
      theme_rs()+
      theme(aspect.ratio = 0.75)+
      scale_x_discrete(name = "Biofilm type", labels = lab_biofilm)+
      scale_y_continuous(name = "Biofilm (log2 OD)", limits = c(-4,4), expand = c(0,0))+
      scale_fill_manual(name = "Medium", values = c("grey50", 'grey80'), labels = lab_medium)+
      scale_color_manual(values = palette_biofilm, labels = lab_biofilm)+
      guides(fill = guide_legend(override.aes = list(shape = 1, size = 4)),
             color = guide_legend(title = "Biofilm type", override.aes = list(size = 4, alpha = 1)))

f1.c <- biofilm %>% 
      filter(category != "biofilm") %>% 
      ggplot(aes(category, log2(OD), fill = medium, group = interaction(medium, category)))+
      geom_violin(size = 0.25, trim = TRUE, scale = 'width', adjust = 1.2)+
      geom_jitter(aes(color = type), alpha = 0.5, size = 1.5, stroke = 0, position = position_jitterdodge(jitter.width = 0.8, dodge.width = 0.9))+
      geom_boxplot(fill = "white", size=0.25, outlier.alpha = 0, width = 0.1, position = position_dodge(width=0.9))+
      theme_rs()+
      theme(aspect.ratio = 0.75)+
      scale_x_discrete(name = "Source", labels = c("Fresh produce", "Soil", "Water"))+
      scale_fill_manual(name = "Medium", values = c("grey50", 'grey80'), labels = lab_medium)+
      scale_y_continuous(name = "Biofilm (log2 OD)", limits = c(-4,4), expand = c(0,0))+
      scale_color_manual(values = palette_biofilm, labels = lab_biofilm)+
      guides(fill = guide_legend(override.aes = list(shape = 1, size = 4)),
             color = guide_legend(title = "Biofilm type", override.aes = list(size = 4, alpha = 1)))

f1.d <- biofilm %>% 
      ggplot(aes(phylogroup, log2(OD)))+
      facet_wrap(~ medium, ncol = 1)+
      geom_violin(size = 0.25, trim = TRUE, scale = 'width', adjust = 1.2, fill = 'grey90')+
      geom_jitter(aes(color = type), alpha = 0.5, size = 1.5, stroke = 0, width = 0.2)+
      geom_boxplot(fill = "white", size = 0.25, outlier.alpha = 0, width = 0.1, position = position_dodge(width = 0.9))+
      theme_rs()+
      theme(aspect.ratio = 0.25)+
      scale_x_discrete(name = "Phylogroup", labels = c("A", "B1", "B2", "C", "D", "E", "F", "n.d."))+
      scale_y_continuous(name = "Biofilm (log2 OD)", limits = c(-4,4), expand = c(0,0), breaks = seq(-4,4,2))+
      scale_color_manual(values = palette_biofilm, labels = lab_biofilm)+
      guides(color = guide_legend(title = "Biofilm type", override.aes = list(size = 4, alpha = 1)))


((f1.a + f1.b)/(f1.c + f1.d) +
            plot_annotation(tag_levels = "A") +
            plot_layout(guides = "collect"))
ggsave(here("output", "fig1.pdf"), width = 7.2, dpi = 300)