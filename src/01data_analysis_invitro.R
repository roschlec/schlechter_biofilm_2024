#!/usr/bin/env Rscript

##    In vitro Biofilm Data
#     Dependencies
library(tidyverse)
library(here)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)

#     Dependencies
source(here("src", "theme_rs.R"))

#     Labels
lab_biofilm <- c("None/Weak", "Moderate", "Strong", "Extreme")
lab_medium <- c("ABTCAA", "LB-NaCl(0)")

#     Color palette biofilm type
palette_biofilm <- brewer.pal(name = "BuPu", 9)[c(3,5,7,9)]
names(palette_biofilm) <- c("none/weak", "moderate", "strong", "extreme")

#     Import Metadata
strain = read.csv(here("data", "strain.csv"), header = T)

#     Import in vitro biofilm dataset
biofilm = read.csv(here("data", "biofilm.csv"), header = T) %>% 
      left_join(., strain, by = c("year", "strain")) %>% 
      mutate(type = case_when(
                  type == "weak" ~ "none/weak",
                  type == "none" ~ "none/weak",
                  TRUE ~ type),
             type = factor(type, levels = c("none/weak", "moderate", "strong", "ultra"), 
                           labels = c("none/weak", "moderate", "strong", "extreme")),
             year = factor(year),
             strain = factor(strain, levels = unique(strain)),
             medium = factor(medium, labels = c("ABTCAA", "LBNaCl")),
             source = factor(source),
             category = factor(category),
             phylogroup = factor(phylogroup))

#     Selected strains for ABR and in planta experiments
selected_strains <- c("H2", "B456", "C1", "C13", "B368", "B545", "C30", "B466", "B471", "C15", "C160")

#     Data set for correlations between media
corBiofilmOD = biofilm %>% 
      pivot_wider(id_cols = c(strain, phylogroup, category), names_from = medium, values_from = OD)

labels_strains = corBiofilmOD %>% 
      filter(strain %in% selected_strains)

#     In vitro biofilm strength
biofilm_strength = biofilm %>% 
      mutate(dummy_type = case_when(
            type=="none/weak" ~ 0,
            type=="moderate" ~ 1,
            type=="strong" ~ 2,
            type=="extreme" ~ 3)) %>% 
      pivot_wider(id_cols = strain, names_from = medium, values_from = dummy_type) %>% 
      mutate(
            is_LBhigher = case_when(ABTCAA == LBNaCl ~ "same", 
                                    ABTCAA < LBNaCl ~ "high",
                                    ABTCAA > LBNaCl ~ "low"),
            mean_type = (ABTCAA + LBNaCl)/2) %>% 
      pivot_longer(cols=c(ABTCAA, LBNaCl), names_to = "medium", values_to = "dummy_type") %>% 
      left_join(., biofilm, by = c("strain", "medium"))

#     In vitro antibiotic resistance
abr <- read.csv(here("data", "strain_abr.csv"), header = TRUE)

abr_long <- abr %>% 
      separate_rows(antibiotic.resistance, sep = ", ") %>% 
      mutate(resistance = 1) %>% 
      pivot_wider(id_cols = c(strain, source, biofilm.in.LB.NaCl.0., biofilm.in.ABTCAA, PG),
                  names_from = antibiotic.resistance,
                  values_from = resistance,
                  values_fill = 0) %>% 
      mutate(sum_ar = rowSums(across(where(is.numeric))),
             mar = sum_ar/25)

#     In planta
#     CFU dataset
cfu = read.csv(here("data", "raw", "cfu.csv"), header = T) %>% 
      drop_na %>% 
      filter(copies > 0)

#     Combining in planta and in vitro data sets
cfu_biofilm = cfu %>% 
      left_join(., biofilm[biofilm$medium=="ABTCAA",], by = "strain") %>%  
      dplyr::select(strain, exp, dpi, cfu, copies, type, OD, category, phylogroup) %>% 
      mutate(logcfu = log10(cfu), logcopies = log10(copies),
             strain = factor(strain, levels = selected_strains))

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
cfu_biofilm %>% 
      ggplot(aes(logcfu, logcopies, fill = strain))+
      geom_point(pch = 21, alpha = 0.6, size = 1.5)+
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
##    Data analysis
#     GLS
biofilm_filter = biofilm %>% 
      filter(category!='biofilm' & type != 'none/weak')

#     Regression models in ABTCAA
lm_abtcaa   = lm(log2(OD) ~ type * category, data = biofilm_filter[biofilm_filter$medium=='ABTCAA',])
summary(lm_abtcaa)
lm_abtcaa %>% resid %>% qqnorm
bartlett.test(log2(biofilm_filter$OD[biofilm_filter$medium == 'ABTCAA']), 
              g = biofilm_filter$type[biofilm_filter$medium == 'ABTCAA'])
anova(lm_abtcaa)
emmeans(lm_abtcaa, ~ category | type, data = biofilm_filter[biofilm_filter$medium=='ABTCAA',]) %>% contrast(., 'pairwise', adjust = "bonferroni", type = "response")
emmeans(lm_abtcaa, ~ type | category, data = biofilm_filter[biofilm_filter$medium=='ABTCAA',]) %>% contrast(., 'pairwise', adjust = "bonferroni", type = "response")

#     Regression models in LB-NaCl(0)
lm_lbnonacl = lm(log2(OD) ~ type * category, data = biofilm_filter[biofilm_filter$medium=='LBNaCl',])
summary(lm_lbnonacl)
lm_lbnonacl %>% resid %>% qqnorm
anova(lm_lbnonacl)
emmeans(lm_lbnonacl, ~ category | type, data = biofilm_filter[biofilm_filter$medium=='LBNaCl',]) %>% contrast(., "pairwise", adjust = "bonferroni", type = "response")
emmeans(lm_lbnonacl, ~ type | category, data = biofilm_filter[biofilm_filter$medium=='LBNaCl',]) %>% contrast(., "pairwise", adjust = "bonferroni", type = "response")
