#!/usr/bin/env Rscript

##    Clean up data biofilm in vitro and in planta
#     Dependencies
library(tidyverse)
library(here)

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
             mar = sum_ar/32)

#     In planta
#     CFU dataset
cfu = read.csv(here("data", "cfu.csv"), header = T) %>% 
      drop_na %>% 
      filter(copies > 0)

#     Combining in planta and in vitro data sets
cfu_biofilm = cfu %>% 
      left_join(., biofilm[biofilm$medium=="ABTCAA",], by = "strain") %>%  
      dplyr::select(strain, exp, dpi, cfu, copies, type, OD, category, phylogroup) %>% 
      mutate(logcfu = log10(cfu), logcopies = log10(copies),
             strain = factor(strain, levels = selected_strains))

#     Remove outliers CFU and qPCR data
lm_cfu_qpcr <- lm(logcopies ~ logcfu, data = cfu_biofilm)
cds_cfu_qpcr <- cooks.distance(lm_cfu_qpcr)

corr_cfu_biofilm <- cds_cfu_qpcr %>% 
      as_tibble() %>% 
      rename(cookD = value) %>% 
      cbind(cfu_biofilm, .) %>% 
      filter(cookD < 4*mean(cookD))

#     Data
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