#!/usr/bin/env Rscript

##    In vitro Biofilm Data
#     Dependencies
library(tidyverse)
library(here)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(ggpubr)
library(patchwork)

#     Dependencies
source(here("src", "00data_clean.R"))
source(here("src", "theme_rs.R"))

#     Labels
lab_biofilm <- c("None/Weak", "Moderate", "Strong", "Extreme")
names(lab_biofilm) <- c("none/weak", "moderate", "strong", "extreme")

lab_medium <- c("ABTCAA", "LB-NaCl(0)")
names(lab_medium) <- c("ABTCAA", "LBNaCl")

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
      labs(y = bquote("Biofilm (" ~ log[2] ~ OD[600] ~ ")"), x = "Medium")+
      scale_fill_manual(values = palette_biofilm, labels = lab_biofilm)+
      scale_color_manual(values = c("#009900", "#cc0000", "black"))+
      scale_y_continuous(limits = c(-4,4), expand = c(0,0))+
      scale_x_discrete(labels = lab_medium)+
      stat_compare_means(method = 't.test', aes(label=..p.signif..), size = 6, label.x = 1.5, label.y = 3, comparisons = list(c("ABTCAA","LBNaCl")))+
      guides(color = "none", 
             fill = guide_legend(title = "Biofilm type", override.aes = list(size = 4, alpha = 1)))

f1.b <- biofilm %>% 
      ggplot(aes(x = type, logOD, fill = medium, group = interaction(medium,type)))+
      geom_violin(size=0.25, trim = TRUE, scale = 'width', adjust = 1.2)+
      geom_jitter(aes(color = type), alpha = 0.5, size = 1.5, stroke = 0, position = position_jitterdodge(jitter.width = 0.8, dodge.width = 0.9))+
      geom_boxplot(fill = "white", size = 0.25, outlier.alpha = 0, width = 0.1, position = position_dodge(width=0.9))+
      theme_rs()+
      theme(aspect.ratio = 0.75)+
      scale_x_discrete(name = "Biofilm type", labels = lab_biofilm)+
      scale_y_continuous(name = bquote("Biofilm (" ~ log[2] ~ OD[600] ~ ")"), limits = c(-4,4), expand = c(0,0))+
      scale_fill_manual(name = "Medium", values = c("grey50", 'grey80'), labels = lab_medium)+
      scale_color_manual(values = palette_biofilm, labels = lab_biofilm)+
      guides(fill = guide_legend(override.aes = list(shape = 1, size = 4)),
             color = guide_legend(title = "Biofilm type", override.aes = list(size = 4, alpha = 1)))

f1.c <- biofilm %>% 
      filter(category != "biofilm") %>% 
      ggplot(aes(category, logOD, fill = medium, group = interaction(medium, category)))+
      geom_violin(size = 0.25, trim = TRUE, scale = 'width', adjust = 1.2)+
      geom_jitter(aes(color = type), alpha = 0.5, size = 1.5, stroke = 0, position = position_jitterdodge(jitter.width = 0.8, dodge.width = 0.9))+
      geom_boxplot(fill = "white", size=0.25, outlier.alpha = 0, width = 0.1, position = position_dodge(width=0.9))+
      theme_rs()+
      theme(aspect.ratio = 0.75)+
      scale_x_discrete(name = "Source", labels = c("Fresh produce", "Soil", "Water"))+
      scale_fill_manual(name = "Medium", values = c("grey50", 'grey80'), labels = lab_medium)+
      scale_y_continuous(name = bquote("Biofilm (" ~ log[2] ~ OD[600] ~ ")"), limits = c(-4,4), expand = c(0,0))+
      scale_color_manual(values = palette_biofilm, labels = lab_biofilm)+
      guides(fill = guide_legend(override.aes = list(shape = 1, size = 4)),
             color = guide_legend(title = "Biofilm type", override.aes = list(size = 4, alpha = 1)))

f1.d <- biofilm %>% 
      ggplot(aes(phylogroup, logOD))+
      facet_wrap(~ medium, ncol = 1, labeller = labeller(medium = lab_medium))+
      geom_violin(size = 0.25, trim = TRUE, scale = 'width', adjust = 1.2, fill = 'grey90')+
      geom_jitter(aes(color = type), alpha = 0.5, size = 1.5, stroke = 0, width = 0.2)+
      geom_boxplot(fill = "white", size = 0.25, outlier.alpha = 0, width = 0.1, position = position_dodge(width = 0.9))+
      theme_rs()+
      theme(aspect.ratio = 0.25)+
      scale_x_discrete(name = "Phylogroup", labels = c("A", "B1", "B2", "C", "D", "E", "F", "n.d."))+
      scale_y_continuous(name = bquote("Biofilm (" ~ log[2] ~ OD[600] ~ ")"), limits = c(-4,4), expand = c(0,0), breaks = seq(-4,4,2))+
      scale_color_manual(values = palette_biofilm, labels = lab_biofilm)+
      guides(color = guide_legend(title = "Biofilm type", override.aes = list(size = 4, alpha = 1)))


((f1.a + f1.b)/(f1.c + f1.d) +
            plot_annotation(tag_levels = "A") +
            plot_layout(guides = "collect"))
ggsave(here("output", "fig1.pdf"), width = 7.2, dpi = 300)

####  FIGURE 2 ####
corBiofilmOD %>% 
      filter(category != "biofilm") %>% 
      ggplot(., aes(log2(ABTCAA), log2(LBNaCl)))+
      geom_point(aes(color = category), alpha = 0.75, stroke = 0, size = 2)+
      geom_text_repel(data = labels_strains, aes(label=strain), color = 'black',
                      force_pull   = 0.5, # do not pull toward data points
                      nudge_y      = 2,
                      direction    = "x",
                      angle        = 90,
                      hjust        = -1,
                      segment.color = "grey80",
                      segment.size = 0.4)+
      geom_abline(aes(slope=1, intercept=0))+
      geom_point(data = labels_strains, aes(fill=category), pch = 21, color = 'black', size = 2)+
      theme_rs()+
      theme(aspect.ratio = 0.75)+
      guides(color = guide_legend(title = "Source", override.aes = list(alpha = 1, size = 3)), fill = "none")+
      scale_fill_manual(values = c("#F8766D", '#00BFC4'))+
      scale_color_manual(values = c("#F8766D", '#7CAE00', '#00BFC4'), labels = c("Fresh produce", "Soil", "Water"))+
      scale_x_continuous(limits=c(-4,2), breaks = seq(-4,2,2))+
      scale_y_continuous(limits=c(-4,4), breaks = seq(-4,4,2))+
      labs(y = bquote("Biofilm in LB-NaCl(0) (" ~ log[2] ~ OD[600] ~ ")"), 
           x = bquote("Biofilm in ABTCAA (" ~ log[2] ~ OD[600] ~ ")"))
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

####  FIGURE 4 #####
corr_cfu_biofilm %>% 
      mutate(type = factor(type, levels = names(palette_biofilm))) %>% 
      ggplot(aes(x = logcfu, y = logcopies, group = type, fill = type))+
      facet_wrap(~ type, ncol = 4, labeller = labeller(type = lab_biofilm))+
      geom_point(pch = 21, alpha = 0.6, size = 1.5)+
      geom_abline(aes(intercept = intercept, slope = m), data = results)+
      geom_abline(slope = 1, intercept = 0, linetype = "dashed")+
      geom_text(data = results, aes(x = 1, y = 10.5, label = paste('b = ',sprintf('%.2f', m))), hjust = 'inward')+
      geom_text(data = correlation_qpcr_cfu, aes(x = 1, y = 9.5, label = paste('r = ',sprintf('%.2f', cor))), hjust = 'inward')+
      geom_text(data = correlation_qpcr_cfu, aes(x = 1, y = 8.5, label = paste('p < 0.05')), hjust = 'inward')+
      theme_rs()+
      theme(aspect.ratio = 1)+
      scale_y_continuous(name = bquote("Bacterial density ["~log[10] ~ "CFU" ~ gFW^-1~"]"), limits = c(2,11), breaks = seq(2,11,2))+
      scale_x_continuous(name = bquote("Gene copy number ["~log[10] ~ "yccT copies" ~ gFW^-1~"]"), limits = c(1,10), breaks = seq(2,10,2))+
      scale_fill_manual(values = palette_biofilm)+
      guides(fill = "none")
ggsave(here("output", "fig4.pdf"), width = 7, dpi = 300)

####  FIGURE 5 #####
f5a <- df_M4id %>% 
      ggplot(aes(x = dpi, y = estimate, fill = type))+
      facet_wrap(~type, ncol = 4, labeller = labeller(type = lab_biofilm))+
      geom_jitter(data = cfu_biofilm, aes(x = dpi, y = logcopies, color = type),
                  width = 0.9, alpha = 0.8, size = 2, stroke = 0)+
      geom_point(size = 2, stroke = 0.5, fill = "black", color = "grey", pch=21, position = position_dodge(width = 2))+
      geom_text(aes(label = .group, y = 13), hjust = 1, size = 3, 
                position = position_dodge(width = 2), angle = 90)+
      geom_line(alpha = 0.5, linetype = "dashed")+
      theme_rs()+
      theme(aspect.ratio = 1)+
      guides(color = guide_legend(title = "Biofilm type", override.aes = list(size = 4, alpha = 1)), fill = "none")+
      scale_y_continuous(name = bquote(log[10] ~ "yccT copies" ~ gFW^-1), limits = c(2,13), breaks = seq(2,12,2))+
      scale_x_continuous(name = "Time [dpi]", limits = c(0, 25), breaks = c(0, 3, 7, 14, 21))+
      scale_color_manual(values = palette_biofilm, labels = lab_biofilm)

f5b <- df_M4sid %>% 
      ggplot(aes(x = strain, y = estimate, fill = type))+
      facet_wrap(~as.factor(dpi), ncol = 5)+
      geom_jitter(data = cfu_biofilm, aes(x = strain, y = logcopies, color = type),
                  width = 0.2, alpha = 0.8, size = 2, stroke = 0)+
      geom_point(size = 2, stroke = 0.5, fill = "black", color = "grey", pch=21, position = position_dodge(width = 2))+
      geom_text(aes(label = .group, y = 13), hjust = 1, size = 3, 
                position = position_dodge(width = 2), angle = 90)+
      theme_rs()+
      theme(aspect.ratio = 1,
            axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1, size = 8))+
      guides(color = guide_legend(title = "Biofilm type", override.aes = list(size = 4, alpha = 1)),
             fill = "none")+
      scale_y_continuous(name = bquote(log[10] ~ "yccT copies" ~ gFW^-1), limits = c(2,13), breaks = seq(2,12,2))+
      scale_x_discrete(name = "Strain")+
      scale_color_manual(values = palette_biofilm, labels = lab_biofilm)

f5a/f5b+
      plot_annotation(tag_levels = "A")+
      plot_layout(guides = "collect")
ggsave(here("output", "fig5.pdf"), width = 9, dpi = 300)
ggsave(here("output", "fig5.png"), width = 9, dpi = 300)


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
      labs(x = bquote("Biofilm (" ~ log[2] ~ OD[600] ~ ")"), y = "Strain")+
      theme_rs()
fS1b <- corBiofilmOD %>% 
      arrange(-LBNaCl) %>% 
      slice(58:115) %>% 
      ggplot()+
      geom_segment(aes(x = log2(ABTCAA), xend = log2(LBNaCl), y = reorder(strain, LBNaCl), yend = strain), color="black") +
      geom_point(aes(x = log2(ABTCAA), y = strain), color = rgb(0.2, 0.7, 0.1, 0.5), size = 3, alpha = .9) +
      geom_point(aes(x = log2(LBNaCl), y = strain), color = rgb(0.7, 0.2, 0.1, 0.5), size = 3, alpha = .9)+
      scale_x_continuous(limits = c(-4,4))+
      labs(x = bquote("Biofilm (" ~ log[2] ~ OD[600] ~ ")"), y = "")+
      theme_rs()
fS1c <- corBiofilmOD %>% 
      arrange(-LBNaCl) %>% 
      slice(116:174) %>% 
      ggplot()+
      geom_segment(aes(x = log2(ABTCAA), xend = log2(LBNaCl), y = reorder(strain, LBNaCl), yend = strain), color="black") +
      geom_point(aes(x = log2(ABTCAA), y = strain), color = rgb(0.2, 0.7, 0.1, 0.5), size = 3, alpha = .9) +
      geom_point(aes(x = log2(LBNaCl), y = strain), color = rgb(0.7, 0.2, 0.1, 0.5), size = 3, alpha = .9)+
      scale_x_continuous(limits = c(-4,4))+
      labs(x = bquote("Biofilm (" ~ log[2] ~ OD[600] ~ ")"), y = "")+
      theme_rs()

fS1a + fS1b + fS1c
ggsave(here("output", "figs1.pdf"), width = 7.2, height = 8, dpi = 300)

####  FIGURE S2   ####
biofilm_mar %>% 
      ggplot(aes(x = mar, y = logOD))+
      facet_grid(cols = vars(medium))+
      geom_point(aes(color = type), size = 3)+
      geom_text(aes(x = 0.1, y = 3, label = paste("p =", sprintf("%.3f", p.value))), 
                hjust = "inward", data = df_mar)+
      geom_text(aes(x = 0.1, y = 2.5, label = paste("r =", sprintf("%.2f", cor))), 
                hjust = "inward", data = df_mar)+
      theme_rs()+
      theme(aspect.ratio = 1)+
      scale_color_manual(name = "Biofilm type", values = palette_biofilm)+
      labs(y = bquote("Biofilm (" ~ log[2] ~ OD[600] ~ ")"),
           x = "MAR")
ggsave(here("output", "figs2.pdf"), width = 6, dpi = 300)

####  FIGURE S3   ####
fs3a <- cfu_biofilm_summary_type %>% 
      ggplot(aes(x = dpi, y = cv_copies, fill = type, color = type))+
      geom_line(linewidth = 1, alpha = 0.8)+
      geom_point(pch = 21, size = 2.5, stroke = 0.5, color = "black")+
      theme_rs()+
      theme(aspect.ratio = 1)+
      scale_shape(solid = TRUE)+
      scale_y_continuous(name = "CV [%]", limits = c(0,35), expand = c(0,0))+
      scale_x_continuous(name = "Time [dpi]", limits = c(-2,23), expand = c(0,0))+
      scale_color_manual(values = palette_biofilm, labels = lab_biofilm)+
      scale_fill_manual(values = palette_biofilm, labels = lab_biofilm)+
      guides(color = guide_legend(title = "Biofilm type", override.aes = list(size = 4, alpha = 1)), fill = "none")

fs3b <- cfu_biofilm_summaryAll %>% 
      ggplot(aes(x = dpi, y = cv_copies, fill = type, color = type, group = strain))+
      geom_line(aes(group = interaction(type,strain)), linewidth = 1, alpha = 0.8)+
      geom_point(pch = 21, size = 2.5, stroke = 0.5, color = "black", alpha = 1)+
      theme_rs()+
      theme(aspect.ratio = 1)+
      scale_shape(solid = TRUE)+
      scale_y_continuous(name = "CV [%]", limits = c(0,35), expand = c(0,0))+
      scale_x_continuous(name = "Time [dpi]", limits = c(-2,23), expand = c(0,0))+
      scale_color_manual(values = palette_biofilm, labels = lab_biofilm)+
      scale_fill_manual(values = palette_biofilm, labels = lab_biofilm)+
      guides(color = guide_legend(title = "Biofilm type", override.aes = list(size = 4, alpha = 1)), fill = "none")

fs3c <- cfu_biofilm_summary_exp %>% 
      ggplot(., aes(x = dpi, y = cv_copies, group = as.factor(exp)))+
      geom_line(linewidth = 0.7)+
      geom_point(size = 2.5, stroke = 0 , color = "black", alpha = 1)+
      geom_text_repel(data = lab_exp, aes(label = exp),
                      color = "black",
                      force             = 0.1,
                      nudge_x           = 3,
                      direction         = "y",
                      hjust             = 0,
                      segment.size      = 0.5)+
      theme_rs()+
      theme(aspect.ratio = 1)+
      scale_shape(solid = TRUE)+
      scale_y_continuous(name = "CV [%]", limits = c(0,38), expand = c(0,0), breaks = seq(0,30,10))+
      scale_x_continuous(name = "Time [dpi]", limits = c(-2,27), expand = c(0,0))

fs3a + fs3b + fs3c +
      plot_annotation(tag_levels = "A")+
      plot_layout(guides = "collect")
ggsave(here("output", "figs3.pdf"), width = 9, dpi = 300)

#### FIGURE S4 ####
biofilm_plant_mar %>% 
      ggplot(aes(x = mar, y = logcopies))+
      facet_grid(cols = vars(type), labeller = labeller(type = lab_biofilm))+
      geom_point(aes(color = type, shape = as.factor(dpi)), alpha = 0.8)+
      geom_smooth(method = "lm", se = FALSE, fullrange = TRUE, color = "black", 
                  linewidth = 0.3, alpha = 0.1)+
      theme_rs()+
      theme(aspect.ratio = 1)+
      scale_color_manual(name = "Biofilm type", values = palette_biofilm, labels = lab_biofilm)+
      scale_x_continuous(limits = c(0.05, 0.7), breaks = seq(0.2, 0.6, 0.2))+
      labs(y = bquote(log[10] ~ "yccT copies" ~ gFW^-1),
           x = "MAR")+
      guides(shape = guide_legend(title = "Time [dpi]"))
ggsave(here("output", "figs4.pdf"), width = 9, dpi = 300)
