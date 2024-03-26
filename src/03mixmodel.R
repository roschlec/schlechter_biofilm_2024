#!/usr/bin/env Rscript

###   Mixture models of bacterial populations in planta
#    Libraries
library(here)
library(mixtools)
library(magrittr)
library(patchwork)
library(tidyverse)

#     Dependencies
source(here("src", "01data_analysis_invitro.R"))

#     Set seed
set.seed(1)

#     Define variables
dpi = c("7", "14", "21") # Days post inoculation
population = c("population1", "population2", "populationtotal") # (sub)populations
biofilm_type = c("none/weak", "moderate", "strong", "extreme") # in vitro type of biofilm

#     Plotting function to overlay mixture models
plot_mix_comps <- function(x, mu, sigma, lam) {
      lam * dnorm(x, mu, sigma)
}

##    MIX MODELS
#     Create an empty list to populate with distribution plots
list_plot <- list()

#     Create an empty list to with mixture models
list_mixmdl <- list()
l <- 1

#     Loop over each biofilm type and time
for(i in 1:length(time)){
  for(j in 1:length(biofilm_type)){
        x = cfu_biofilm$logcopies[cfu_biofilm$type==biofilm_type[j] & cfu_biofilm$dpi==dpi[i]]
        m = mean(cfu_biofilm$logcopies[cfu_biofilm$type==biofilm_type[j] & cfu_biofilm$dpi==dpi[i]])
        s = sd(cfu_biofilm$logcopies[cfu_biofilm$type==biofilm_type[j] & cfu_biofilm$dpi==dpi[i]])
        mixmdl = normalmixEM(cfu_biofilm$logcopies[cfu_biofilm$type==biofilm_type[j] & cfu_biofilm$dpi==dpi[i]])
        list_mixmdl[[l]] = c(mixmdl, mean = m, sd = s)
        list_plot[[l]] = data.frame(x = mixmdl$x) %>%
              ggplot() +
              geom_histogram(aes(x, ..density..), binwidth = 0.5, colour = "black", fill = "grey", alpha = 0.8, size=0.5) +
              stat_function(geom = "line", fun = plot_mix_comps,
                        args = list(m, s, lam = 1),
                        colour = "black", lwd = 1) +
              stat_function(geom = "line", fun = plot_mix_comps,
                        args = list(mixmdl$mu[1], mixmdl$sigma[1], lam = mixmdl$lambda[1]),
                        colour = "green", lwd = 1) +
              stat_function(geom = "line", fun = plot_mix_comps,
                        args = list(mixmdl$mu[2], mixmdl$sigma[2], lam = mixmdl$lambda[2]),
                        colour = "magenta", lwd = 1) +
              labs(title=paste(biofilm_type[j], dpi[i], sep="_"))+
              theme_rs()+
              theme(aspect.ratio = 1, 
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank())+
              scale_x_continuous(limits = c(2,12), breaks = seq(2,12,2), name = "")+
              scale_y_continuous(limits = c(0,0.7), name = "")
        names(list_plot)[[l]] <- paste(biofilm_type[j], dpi[i], sep="_")
        names(list_mixmdl)[[l]] <- paste(biofilm_type[j], dpi[i], sep="_")
        l <- l+1
        }
      }

##    Data frame with mixture model parameters
df = data.frame()
for (i in 1:length(list_mixmdl)){
  tmp1 = data.frame(
    label = rep(names(list_mixmdl)[[i]],3),
    success = population,
    lambda = c(list_mixmdl[[i]][[2]], 1),
    mu = c(list_mixmdl[[i]][[3]], list_mixmdl[[i]][[10]]),
    sigma = c(list_mixmdl[[i]][[4]],list_mixmdl[[i]][[11]]))
  
  tmp2 = list_mixmdl[[i]][[6]] %>% 
    data.frame() %>% 
    mutate(success = ifelse(comp.1 > 0.5, "population1", "population2")) %>% 
    group_by(success) %>% 
    tally()
  
  new_row = data.frame(success = "populationtotal", n = sum(tmp2$n))
  
  tmp2 = rbind(tmp2, new_row)
  
  tmp3 = left_join(tmp1,tmp2, by=c("success"))
  
  df = rbind(df, tmp3)
}

df %<>% separate(label, into = c('type', 'dpi'), sep = "_") %>% na.omit
df$type = factor(df$type, levels = c("none/weak", "moderate", "strong", "extreme"))
df$dpi = factor(df$dpi, levels = c("7", "14", "21"))
df$success = factor(df$success, levels = population)

#### FIGURE 6 ####
lapply(names(list_plot), function(.x) {
  ggsave(
    path = here("output"),
    filename = paste0(.x, ".pdf"),
    scale = 1,
    plot = list_plot[[.x]]
  )
})
