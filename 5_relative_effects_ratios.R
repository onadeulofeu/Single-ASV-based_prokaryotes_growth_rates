##En aquest script comparem els diferents tractaments i veiem quins efectes controlen més els diferents grups taxonomics
sessionInfo()

library(tidyverse)
library(magrittr)
library(ggplot2)
library(multipanelfigure) ##per unir diferents gràfics
library(gridExtra)
library(ggpmisc) ##stat_poli_line and equation

#functions needed
# source(src/filter.season.treatment.time.R) ###filtrar les pseudoabund per treatment season i time per calcular regressions
# source(src/multiple.linear.regressions.R) ##per fer multiples regressions d'un dataframe
# source(src/comparing.reg.R) ##per comparar les taxes de creiement entre 4 i 3 punts (les regressions)

#palettes------
palette_large <- c("#a6cee3","#1f78b4", "#009E73","#943132", "#fb9a99", "#FFE867", "#fdbf6f","#ff7f00","#cab2d6",
                   "#6a3d9a", "#d33f6a", "#5562bd", "#c6dec7", "#ef9708", "#c398d9", "#CC6666", "#CABFAB", 
                   "#637b88", "#e3a6ce", "#cee3a6", "#a6b0e3", "#a6e3da", "#e3bba6", "#e0c643", "#84c680", 
                   "#e49273", "#004643", "#80ada0", "#5f5566", "#773344", "#F2BB05", "#cbf3f0", "#ffbc42", 
                   "#cd5334", "#005c69", "#2c1a1d", "#d496a7", "#a6a15e", "#b0413e", "#1d7874", "#fffffc", 
                   "#247ba0", "#fcca46", "#183059", "#ffa737", "#669bbc", "#bcc4db", "#000000")
new_palette_large <- c("#bba76c","#c8811c","#fa8775","#ea5f94","#cd34b5","#4827cc","#a6cee3","#a6a15e","#b0413e",
                       "#1f78b4","#009e73","#61000c","#fb9a99","#ffe867","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a",
                       "#d33f6a","#5562bd","#c6dec7","#992800","#c398d9","#cc6666","#cabfab","#637b88","#e3a6ce",
                       "#cee3a6","#a6b0e3","#75b1a9","#e3bba6","#e0c643","#84c680","#e49273","#003029","#3a5048",
                       "#5f5566","#773344","#ba8b00","#9ac0be","#752800","#cd5334","#005c69","#2c1a1d","#d496a7",
                       "#1d7874","#cccccc","#247ba0","#fcca46","#183059","#ffa737","#669bbc","#bcc4db","#000000")

palette_treatments_remei <- c("#cab2d6","#6a3d9a", "#a6cee3","#1f78b4", "#009E73", "#F2BB05")
palette_seasons <- c( "#005c69",  "#80ada0", "#ffbc42", "#cd5334", "#2c1a1d")
#palette_seasons_4 <- c( "#1f78b4", "#005c69", "#ffbc42", "#cd5334")
palette_seasons_4 <- c("#002562", "#519741", "#ffb900", "#96220a") #verd-blau no es confonen ara en principi

#this function allows us to create many colors from our palette
palf_large <- colorRampPalette(palette_large) 
##to use it scale_color_manual(values=palf_large(x))+ x=number of colors you want
palette_phylums <- c("#fcca46", "#b0413e", "#669bbc", "#009e73", "#ffa737", "#cccccc",
                     "#69267e", "#fb9a99", "#1f78b4", "#8c789d", "#637b88", "#003029",
                     "#e3a6ce", "#002960", "#ba9864", "#005c69")
#fixing colors by phylum ----
# reg_all_slopes_chosen_silva_tax_filt <- reg_all_slopes_chosen_silva_tax %>%
#   filter(slope_chosen_days > 0,
#          pvalue_slope_chosen < 0.05) %>%
#   group_by(phylum_f) %>%
#   #summarize(n = n())
#   filter(n() > 2)
# 
# names(palette_phylums) = unique(reg_all_slopes_chosen_silva_tax_filt$phylum_f)

palette_phylums_assigned <- c('Proteobacteria' = "#fcca46","Bacteroidota" = "#669bbc" , 'Crenarchaeota' = "#b0413e",
                              'Cyanobacteria' = "#009e73",'Actinobacteriota' = "#ffa737", 'Verrucomicrobiota' = "#cccccc",
                              'Planctomycetota' = "#69267e", 'Acidobacteriota' = "#1f78b4",'Bdellovibrionota' = "#8c789d",
                              'Firmicutes' = "#637b88", 'Myxococcota'= "#003029", 'Nitrospinota'= "#e3a6ce",
                              'Campilobacterota' = "#002960", 'Deinococcota'= "#ba9864",'Fusobacteriota' ="#fb9a99",
                              'Desulfobacterota' = "#005c69", 'Methylomirabilota' = "#000000" ,
                              'Gemmatimonadota' = "#c55e5c", 'Chloroflexi' = "#00d198", 'Spirochaetota' = "#5cb1d6",
                              'Calditrichota' = "#8a007a", 'Halobacterota' = "#b79c64", 'Nitrospirota' = "#41815f",
                              'Dependentiae' = "#5b95e5", 'Patescibacteria' = "#33af9c",'Cloacimonadota' = "#fbed5c",
                              'Synergistota' = "#ce7800", 'Abditibacteriota' = "#87878b", 'Deferribacterota' = "#4dbaa9")

##Facet_wrap/ grid labels----
##two rows
effects_labels <-  as_labeller(c(effect_top_down_virus = 'Top-down effect \nviruses (VL/DL)',
                                 effect_top_down_grazers_dark = 'Top-down effect \ngrazers dark (PD/CD)',
                                 effect_top_down_grazers_light = 'Top-down effect \ngrazers light (PL/CL)',
                                 effect_bottom_up = 'Bottom-up effect \n(DL/PL)',
                                 effect_light_C = 'Light effect \n(CL/CD)',
                                 effect_light_P = 'Light effect \n(PL/PD)'))
##one row
effects_labels2 <-  as_labeller(c(effect_top_down_virus = 'Top-down effect viruses (VL-DL)',
                                  effect_top_down_grazers_dark = 'Top-down effect grazers dark (PD-CD)',
                                  effect_top_down_grazers_light = 'Top-down effect grazers light (PL-CL)',
                                  effect_bottom_up = 'Bottom-up effect (DL-PL)',
                                  effect_light_C = 'Light effect (CL-CD)',
                                  effect_light_P = 'Light effect (PL-PD)'))

#Import data go to line 2640 for COMPLETE DATASET WITH IN SITU DATA
reg_all_slopes_chosen_silva_tax <- read.csv("data/intermediate_files/reg_all_slopes_chosen_silva_tax.csv", sep=",") %>%
  filter(season != "Early_fall")

#OPCIÓ 1: Calculating relative effects (divisions) HO DESCARTEM MASSES UNCALCULATED-----
##Relative  effects at ASV level calculation
relative_effects <- reg_all_slopes_chosen_silva_tax %>% 
  distinct(treatment, season, asv_num, domain, phylum, class, order, family, genus, species, tax_ed, .keep_all = TRUE) %>%
  filter(pvalue_slope_chosen < 0.05) %>%
  pivot_wider(names_from = treatment, values_from = c(slope_chosen_days, slope_chosen, pvalue_slope_chosen), 
              id_cols = c(asv_num, season, domain, phylum, class, order, family, genus, species, tax_ed)) %>%
  as_tibble %>%
  mutate(across(!c(asv_num, domain, phylum, class, order, family, genus, 
                   species, tax_ed, season), as.numeric)) %>%
  mutate(effect_top_down_virus = slope_chosen_days_VL/slope_chosen_days_DL,
         effect_top_down_grazers_dark = slope_chosen_days_PD/slope_chosen_days_CD,
         effect_top_down_grazers_light = slope_chosen_days_PL/slope_chosen_days_CL,
         effect_bottom_up = slope_chosen_days_DL/slope_chosen_days_PL,
         effect_light_C = slope_chosen_days_CL/slope_chosen_days_CD,
         effect_light_P = slope_chosen_days_PL/slope_chosen_days_PD)

relative_effects %>%
  colnames()
relative_effects_l <- relative_effects %>%
  pivot_longer(cols = starts_with('effect'),
               names_to = 'effects',
               values_to = 'ratio') %>%
  as_tibble()

relative_effects_l <- relative_effects_l %>%
  mutate(phylum_f = as_factor(phylum),
         family_f = as_factor(family),
         class_f = as_factor(class))

### % of relative effects uncalculated at asv level ----
total_rows <- relative_effects %>%
  summarise(n = n()) %>%
  as.numeric()

relative_effects_nas <- relative_effects %>%
  select(starts_with("effect")) %>% 
  summarise_all(funs(sumNA = sum(is.na(.))))

relative_effects_nas/total_rows #total uncalculated relative effects

create_summary <- function(inputdata, ...) {
  inputdata %>%
    summarise_at(vars(...), list(na= ~sum(is.na(.)), nonna= ~sum(!is.na(.))))
}

asv_ratios_nas <- create_summary(relative_effects, effect_top_down_virus,
               effect_top_down_grazers_dark,
               effect_top_down_grazers_light,
               effect_bottom_up, effect_light_C, effect_light_P)

percentage <- asv_ratios_nas %>%
  select(contains('_na')) %>% 
  summarize(sum())
  mutate(percentage = select(contains('_na'))/select(((contains('_na'))+select((contains('nonna'))))))

## Relative effects at Phylum level calculation-----
relative_effects_phylum <- reg_all_slopes_chosen_silva_tax %>% 
  distinct(treatment, season, asv_num, domain, phylum, phylum, order, family, genus, species, tax_ed, .keep_all = TRUE) %>%
  filter(pvalue_slope_chosen < 0.05) %>%
  group_by(phylum,treatment,season) %>%
  summarize(slope_chosen_days_mean = mean(slope_chosen_days)) %>%
  pivot_wider(names_from = treatment, values_from = c(slope_chosen_days_mean), 
              id_cols = c(season, phylum)) %>%
  as_tibble %>%
  mutate(across(!c(phylum, season), as.numeric)) %>%
  mutate(effect_top_down_virus= VL/DL,
         effect_top_down_grazers_dark = PD/CD,
         effect_top_down_grazers_light = PL/CL,
         effect_bottom_up = DL/PL,
         effect_light_C = CL/CD,
         effect_light_P = PL/PD)

relative_effects_phylum_l <- relative_effects_phylum %>%
  pivot_longer(cols = starts_with('effect'),
               names_to = 'effects',
               values_to = 'ratio') %>%
  as_tibble()

relative_effects_phylum_l <- relative_effects_phylum_l %>%
  mutate(phylum_f = as_factor(phylum))

## Relative effects at Class level calculation
relative_effects_class <- reg_all_slopes_chosen_silva_tax %>% 
  distinct(treatment, season, asv_num, domain, phylum, class, order, genus, species, tax_ed, asv_num, order, .keep_all = TRUE) %>%
  filter(pvalue_slope_chosen < 0.05) %>%
  group_by(treatment, season, domain, phylum, class) %>%
  summarise(slope_chosen_days_mean = mean(slope_chosen_days),
            slope_chosen_mean = mean(slope_chosen),
            na.rm = TRUE) %>%
  pivot_wider(names_from = treatment, values_from = c(slope_chosen_days_mean, slope_chosen_mean), 
              id_cols = c(season, domain, phylum, class)) %>%
  as_tibble() %>%
  mutate(across(!c(domain, phylum, class, season), as.numeric)) %>%
  mutate(effect_top_down_virus = slope_chosen_days_mean_VL/slope_chosen_days_mean_DL,
         effect_top_down_grazers_dark = slope_chosen_days_mean_PD/slope_chosen_days_mean_CD,
         effect_top_down_grazers_light = slope_chosen_days_mean_PL/slope_chosen_days_mean_CL,
         effect_bottom_up = slope_chosen_days_mean_DL/slope_chosen_days_mean_PL,
         effect_light_C = slope_chosen_days_mean_CL/slope_chosen_days_mean_CD,
         effect_light_P = slope_chosen_days_mean_PL/slope_chosen_days_mean_PD)

relative_effects_class_l <- relative_effects_class %>%
  pivot_longer(cols = starts_with('effect'),
               names_to = 'effects',
               values_to = 'ratio') %>%
  as_tibble()

## Relative effects at Order level calculation
relative_effects_order <- reg_all_slopes_chosen_silva_tax %>% 
  distinct(treatment, season, asv_num, domain, phylum, class, order, genus, species, tax_ed, asv_num, order, .keep_all = TRUE) %>%
  filter(pvalue_slope_chosen < 0.05) %>%
  group_by(treatment, season, domain, phylum, class, order) %>%
  summarise(slope_chosen_days_mean = mean(slope_chosen_days),
            slope_chosen_mean = mean(slope_chosen),
            na.rm = TRUE) %>%
  pivot_wider(names_from = treatment, values_from = c(slope_chosen_days_mean, slope_chosen_mean), 
              id_cols = c(season, domain, phylum, class, order)) %>%
  as_tibble() %>%
  mutate(across(!c(domain, phylum, class, order, season), as.numeric)) %>%
  mutate(effect_top_down_virus = slope_chosen_days_mean_VL/slope_chosen_days_mean_DL,
         effect_top_down_grazers_dark = slope_chosen_days_mean_PD/slope_chosen_days_mean_CD,
         effect_top_down_grazers_light = slope_chosen_days_mean_PL/slope_chosen_days_mean_CL,
         effect_bottom_up = slope_chosen_days_mean_DL/slope_chosen_days_mean_PL,
         effect_light_C = slope_chosen_days_mean_CL/slope_chosen_days_mean_CD,
         effect_light_P = slope_chosen_days_mean_PL/slope_chosen_days_mean_PD)

relative_effects_order_l <- relative_effects_order %>%
  pivot_longer(cols = starts_with('effect'),
               names_to = 'effects',
               values_to = 'ratio') %>%
  as_tibble()

### % of relative effects uncalculated at order level----
total_rows <- relative_effects_order %>%
  summarise(n = n()) %>%
  as.numeric()

relative_effects_nas_order <- relative_effects_order %>%
  select(starts_with("effect")) %>% 
  summarise_all(funs(sumNA = sum(is.na(.))))

relative_effects_nas_order/total_rows #total uncalculated relative effects

## Relative  effects at family level calculation ------
relative_effects_family <- reg_all_slopes_chosen_silva_tax %>% 
  distinct(treatment, season, asv_num, domain, phylum, class, order, genus, species, tax_ed, asv_num, family, .keep_all = TRUE) %>%
  filter(pvalue_slope_chosen < 0.05) %>%
  group_by(treatment, season, domain, phylum, class, order, family) %>%
  summarise(slope_chosen_days_mean = mean(slope_chosen_days),
            slope_chosen_mean = mean(slope_chosen),
            na.rm = TRUE) %>%
  pivot_wider(names_from = treatment, values_from = c(slope_chosen_days_mean, slope_chosen_mean), 
              id_cols = c(season, domain, phylum, class, order, family)) %>%
  as_tibble() %>%
  mutate(across(!c(domain, phylum, class, order, family, season), as.numeric)) %>%
  mutate(effect_top_down_virus = slope_chosen_days_mean_VL/slope_chosen_days_mean_DL,
         effect_top_down_grazers_dark = slope_chosen_days_mean_PD/slope_chosen_days_mean_CD,
         effect_top_down_grazers_light = slope_chosen_days_mean_PL/slope_chosen_days_mean_CL,
         effect_bottom_up = slope_chosen_days_mean_DL/slope_chosen_days_mean_PL,
         effect_light_C = slope_chosen_days_mean_CL/slope_chosen_days_mean_CD,
         effect_light_P = slope_chosen_days_mean_PL/slope_chosen_days_mean_PD)

relative_effects_family_l <- relative_effects_family %>%
  pivot_longer(cols = starts_with('effect'),
               names_to = 'effects',
               values_to = 'ratio') %>%
  as_tibble()

### % of relative effects uncalculated at family level----
total_rows <- relative_effects_family %>%
  summarise(n = n()) %>%
  as.numeric()

relative_effects_nas_family <- relative_effects_family %>%
  select(starts_with("effect")) %>% 
  summarise_all(funs(sumNA = sum(is.na(.))))

relative_effects_nas_family/total_rows #total uncalculated relative effects


##a nivell de diferències
total_rows <- effects_difference_family %>%
  summarise(n = n()) %>%
  as.numeric()

differences_effects_nas_family <- effects_difference_family %>%
  select(starts_with("effect")) %>% 
  summarise_all(funs(sumNA = sum(is.na(.))))

differences_effects_nas_family/total_rows #total uncalculated relative effects

## Reordering factors for ploting all datasets -----
relative_effects_phylum_l$season <-relative_effects_phylum_l$season %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))

relative_effects_class_l$season <-relative_effects_class_l$season %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))

relative_effects_order_l$season <-relative_effects_order_l$season %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))

relative_effects_family_l$season <-relative_effects_family_l$season %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))

relative_effects_l$season <-relative_effects_l$season %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))

relative_effects_phylum_l$effects <- relative_effects_phylum_l$effects %>% 
  factor(levels=(c("effect_bottom_up", "effect_top_down_grazers_dark", 
                   "effect_top_down_grazers_light", "effect_top_down_virus",
                   "effect_light_C", "effect_light_P")))

relative_effects_class_l$effects <- relative_effects_class_l$effects %>% 
  factor(levels=(c("effect_bottom_up", "effect_top_down_grazers_dark", 
                   "effect_top_down_grazers_light", "effect_top_down_virus",
                   "effect_light_C", "effect_light_P")))

relative_effects_order_l$effects <- relative_effects_order_l$effects %>% 
  factor(levels=(c("effect_bottom_up", "effect_top_down_grazers_dark", 
                   "effect_top_down_grazers_light", "effect_top_down_virus",
                   "effect_light_C", "effect_light_P")))

relative_effects_family_l$effects <- relative_effects_family_l$effects %>% 
  factor(levels=(c("effect_bottom_up", "effect_top_down_grazers_dark", 
                   "effect_top_down_grazers_light", "effect_top_down_virus",
                   "effect_light_C", "effect_light_P")))

relative_effects_l$effects <- relative_effects_l$effects %>% 
  factor(levels=(c("effect_bottom_up", "effect_top_down_grazers_dark", 
                   "effect_top_down_grazers_light", "effect_top_down_virus",
                   "effect_light_C", "effect_light_P")))

#creating taxonomic factors
relative_effects_phylum_l <- relative_effects_phylum_l %>%
  mutate(phylum_f = as_factor(phylum))

relative_effects_class_l <- relative_effects_class_l %>%
  mutate(phylum_f = as_factor(phylum),
         class_f = as_factor(class))

relative_effects_order_l <- relative_effects_order_l %>%
  mutate(phylum_f = as_factor(phylum),
         order_f = as_factor(order),
         class_f = as_factor(class))

relative_effects_family_l <- relative_effects_family_l %>%
  mutate(phylum_f = as_factor(phylum),
         class_f = as_factor(class),
         order_f = as_factor(order),
         family_f = as_factor(family))

relative_effects_l <- relative_effects_l %>%
  mutate(domain_f = as_factor(domain),
         phylum_f = as_factor(phylum),
         class_f = as_factor(class), 
         order_f = as_factor(order),
         family_f = as_factor(family),
         asv_f = as_factor(asv_num))

#els ordeno en funció del phylum
relative_effects_class_l$class_f <-  factor(relative_effects_class_l$class_f, 
                                            levels=unique(relative_effects_class_l$class_f[order(relative_effects_class_l$phylum_f)]), 
                                            ordered=TRUE)

relative_effects_order_l$order_f <-  factor(relative_effects_order_l$order_f, 
                                            levels=unique(relative_effects_order_l$order_f[order(relative_effects_order_l$phylum_f,
                                                                                                 relative_effects_order_l$class_f)]), 
                                            ordered=TRUE)

relative_effects_family_l$family_f <-  factor(relative_effects_family_l$family_f, 
                                              levels=unique(relative_effects_family_l$family_f[order(relative_effects_family_l$phylum_f,
                                                                                                     relative_effects_family_l$class_f,
                                                                                                     relative_effects_family_l$order_f)]), 
                                              ordered=TRUE)

relative_effects_l$asv_f <-  factor(relative_effects_l$asv_f, 
                                              levels=unique(relative_effects_l$asv_f[order(relative_effects_l$phylum_f,
                                                                                                     relative_effects_l$class_f,
                                                                                                     relative_effects_l$order_f)]), 
                                              ordered=TRUE)

#se'm ordena malament provo d'ordenar la classe bé
relative_effects_l$class_f <-  factor(relative_effects_l$class_f, 
                                    levels=unique(relative_effects_l$class_f[order(relative_effects_l$domain_f,
                                                                                 relative_effects_l$phylum_f,
                                                                                 relative_effects_l$class_f,
                                                                                 relative_effects_l$order_f,
                                                                                 relative_effects_l$family_f)]), 
                                    ordered=TRUE)



## Plotting relative effects at ASV level----
#Subseting by effects
relative_effects_l %$%
  phylum %>%
  unique()
relative_effects_l %>%
  subset(effects == "effect_light_P") %>%
  ggplot(aes(class_f, ratio, color = phylum_f, na.rm = TRUE))+ #, na.rm = TRUE
  geom_boxplot()+
  #facet_wrap(~effects, labeller =  effects)+
  guides(color = guide_legend(ncol = 2))+
  scale_color_manual(values = c(new_palette_large))+
  theme_bw()+
  theme(strip.text.x = element_text(size = 5), 
        axis.text.x = element_text(angle = 90, size = 10), 
        legend.position = "right")

relative_effects_l %>%
  subset(effects == "effect_light_C") %>%
  ggplot(aes(class_f, ratio, color = phylum_f, na.rm = TRUE))+ #, na.rm = TRUE
  geom_boxplot()+
  #facet_wrap(~effects, labeller =  effects)+
  guides(color = guide_legend(ncol = 2))+
  scale_color_manual(values = c(new_palette_large))+
  theme_bw()+
  theme(strip.text.x = element_text(size = 5), 
        axis.text.x = element_text(angle = 90, size = 10), 
        legend.position = "right")

relative_effects_l %>%
  subset(effects == "effect_top_down_virus") %>%
  ggplot(aes(class_f, ratio, color = phylum_f, na.rm = TRUE))+ #, na.rm = TRUE
  geom_boxplot()+
  #facet_wrap(~effects, labeller =  effects)+
  guides(color = guide_legend(ncol = 2))+
  scale_color_manual(values = c(new_palette_large))+
  theme_bw()+
  theme(strip.text.x = element_text(size = 5), 
        axis.text.x = element_text(angle = 90, size = 10), 
        legend.position = "right")

relative_effects_l %>%
  subset(effects == "effect_top_down_grazers_dark"  | 
           effects == "effect_top_down_grazers_light") %>%
  ggplot(aes(class_f, ratio, color = phylum_f, na.rm = TRUE))+ #, na.rm = TRUE
  geom_boxplot()+
  #facet_wrap(~effects, labeller =  effects)+
  guides(color = guide_legend(ncol = 2))+
  scale_color_manual(values = c(new_palette_large))+
  theme_bw()+
  theme(strip.text.x = element_text(size = 5), 
        axis.text.x = element_text(angle = 90, size = 10), 
        legend.position = "right")

relative_effects_l %>%
  subset(effects == "effect_bottom_up") %>%
  #filter(ratio > 10) %>%
  ggplot(aes(class_f, ratio, color = phylum_f, na.rm = TRUE))+ #, na.rm = TRUE
  geom_boxplot()+
  #facet_wrap(~effects, labeller =  effects)+
  guides(color = guide_legend(ncol = 2))+
  scale_color_manual(values = c(new_palette_large))+
  theme_bw()+
  theme(strip.text.x = element_text(size = 5), 
        axis.text.x = element_text(angle = 90, size = 10), 
        legend.position = "right")

relative_effects_l %>%
  #subset(family == "Alteromonadaceae") %>%
  #filter(ratio > 50 | ratio < -50) %>%
  ggplot(aes(effects, family, fill=ratio, na.rm = TRUE))+
  #geom_tile(aes(fill = ratio), color = 'white', lwd = 1.5, linetype = 0.5)+
  scale_fill_gradient2(low = '#BD4247', high =  '#42BDB8', na.value = '#fdfdfd', midpoint = 0)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

relative_effects_l %>%
  #subset(family == "Alteromonadaceae") %>%
  filter(ratio > 1 | ratio < -1) %>%
  ggplot(aes(order, ratio, na.rm = TRUE))+
  #geom_tile(aes(fill = ratio), color = 'white', lwd = 1.5, linetype = 0.5)+
  #scale_fill_gradient2(low = '#BD4247', high =  '#42BDB8', na.value = '#fdfdfd', midpoint = 0)+
  geom_boxplot()+
  geom_point(aes(color = season))+
  facet_wrap(vars(effects), scales = 'free')+
  scale_color_manual(values = palette_seasons_4)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.text.x = element_text(size = 5), axis.text.x = element_text(angle = 90, size = 10), 
        legend.position = "none")
# 
# relative_effects_l$treatment <-relative_effects_l$treatment %>% 
#   factor(levels=(c("CL", "CD", "PL", "PD", "DL", "VL")))
relative_effects_l$season <- relative_effects_l$season %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))

##Flavobacteriaceae (veiem top-down)
relative_effects_l %>% 
  subset(family == "Flavobacteriaceae") %>%
  #filter(ratio > 50 | ratio < -50) %>%
  ggplot(aes(x = asv_num, y = ratio, na.rm = TRUE))+
  #geom_tile(aes(fill = ratio), color = 'white', lwd = 1.5, linetype = 0.5)+
  #scale_fill_gradient2(low = '#BD4247', high =  '#42BDB8', na.value = '#fdfdfd', midpoint = 0)+
  #geom_boxplot()+
  geom_point(aes(color = season))+
  facet_wrap(vars(effects), scales = 'fixed')+
  scale_color_manual(values = palette_seasons_4)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.text.x = element_text(size = 5), axis.text.x = element_text(angle = 90, size = 10), 
        legend.position = "none")

relative_effects_phylum_l %>% 
  #filter(ratio > 1 | ratio < -1) %>%
  ggplot(aes(x = phylum, y = ratio, na.rm = TRUE))+
  #geom_tile(aes(fill = ratio), color = 'white', lwd = 1.5, linetype = 0.5)+
  #scale_fill_gradient2(low = '#BD4247', high =  '#42BDB8', na.value = '#fdfdfd', midpoint = 0)+
  #geom_boxplot()+
  geom_hline(yintercept = 1)+
  geom_point(aes(color = season))+
  geom_boxplot()+
  labs(y='Ratio', x='Phylum', color = 'Season')+
  facet_wrap(vars(effects), scales = 'fixed', labeller =  effects_labels)+
  scale_color_manual(values = palette_seasons_4)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.text.x = element_text(size = ), axis.text.x = element_text(angle = 90, size = 10), 
        legend.position = "right")

## Plotting relative effects at Class level----
relative_effects_class_l %>% 
  filter(ratio > 1 | ratio < -1) %>%
  ggplot(aes(x = class_f, y = ratio, na.rm = TRUE))+
  #geom_tile(aes(fill = ratio), color = 'white', lwd = 1.5, linetype = 0.5)+
  #scale_fill_gradient2(low = '#BD4247', high =  '#42BDB8', na.value = '#fdfdfd', midpoint = 0)+
  #geom_boxplot()+
  geom_hline(yintercept = 1)+
  geom_point(aes(color = season))+
  geom_boxplot()+
  labs(y='Response ratio', x='Phylum', color = 'Season')+
  facet_wrap(vars(effects), scales = 'free', labeller =  effects_labels)+
  scale_color_manual(values = palette_seasons_4)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.text.x = element_text(size = 10), axis.text.x = element_text(angle = 90, size = 10), 
        legend.position = "right")

## Plotting relative effects at Order level------
relative_effects_order_l %>% 
  filter(ratio > 2 | ratio < -2) %>%
  ggplot(aes(x = order_f, y = ratio, na.rm = TRUE))+
  #geom_tile(aes(fill = ratio), color = 'white', lwd = 1.5, linetype = 0.5)+
  #scale_fill_gradient2(low = '#BD4247', high =  '#42BDB8', na.value = '#fdfdfd', midpoint = 0)+
  #geom_boxplot()+
  geom_point(aes(color = season))+
  #geom_boxplot()+
  geom_hline(yintercept = 1)+
  #facet_grid(vars(effects), scales = 'free')+
  labs(y='Response ratio', x='Phylum', color = 'Season')+
  facet_wrap(vars(effects), scales = 'fixed', labeller =  effects_labels)+
  scale_color_manual(values = palette_seasons_4)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.text.x = element_text(size = 14), axis.text.x = element_text(angle = 90, size = 10), 
        legend.position = "none")

## Ploting subsets by interesting taxonomical levels-----
relative_effects_f <- reg_all_slopes_chosen_silva_tax %>% 
  subset(class != "Proteobacteria") %>%
  distinct(treatment, season, asv_num, domain, phylum, class, order, family, genus, species, tax_ed, .keep_all = TRUE) %>%
  filter(pvalue_slope_chosen < 0.05) %>%
  group_by(family, treatment, season) %>%
  summarize(slope_chosen_days_mean = mean(slope_chosen_days)) %>%
  pivot_wider(names_from = treatment, values_from = c(slope_chosen_days_mean), 
              id_cols = c(season, family)) %>%
  as_tibble %>%
  mutate(across(!c(family, season), as.numeric)) %>%
  mutate(effect_top_down_virus= VL/DL,
         effect_top_down_grazers_dark = PD/CD,
         effect_top_down_grazers_light = PL/CL,
         effect_bottom_up = DL/PL,
         effect_light_C = CL/CD,
         effect_light_P = PL/PD)

relative_effects_f_l <- relative_effects_f %>%
  pivot_longer(cols = starts_with('effect'),
               names_to = 'effects',
               values_to = 'ratio') %>%
  as_tibble()

relative_effects_f_l %>% 
  #filter(ratio > 10 | ratio < -10) %>%
  ggplot(aes(x = family, y = ratio, na.rm = TRUE))+
  #geom_tile(aes(fill = ratio), color = 'white', lwd = 1.5, linetype = 0.5)+
  #scale_fill_gradient2(low = '#BD4247', high =  '#42BDB8', na.value = '#fdfdfd', midpoint = 0)+
  #geom_boxplot()+
  geom_point(aes(color = season))+
  #geom_boxplot()+
  geom_hline(yintercept = 1)+
  labs(y='Response ratio', x='Family', color = 'Season')+
  facet_wrap(vars(effects), scales = 'free', labeller =  effects_labels)+
  scale_color_manual(values = palette_seasons_4)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.text.x = element_text(size = 10), axis.text.x = element_text(angle = 90, size = 10), 
        legend.position = "right")

relative_effects_f_l %>% 
  #filter(ratio > 10 | ratio < -10) %>%
  ggplot(aes(x = season, y = ratio, na.rm = TRUE))+
  #geom_tile(aes(fill = ratio), color = 'white', lwd = 1.5, linetype = 0.5)+
  #scale_fill_gradient2(low = '#BD4247', high =  '#42BDB8', na.value = '#fdfdfd', midpoint = 0)+
  geom_boxplot()+
  geom_point(aes(color = season))+
  #geom_boxplot()+
  labs( x = 'Season', y = 'Response ratio', color = 'Season') +
  #scale_y_continuous(breaks = c(-500, -200, -300,-50,  1, 50, 100,  300, 500))+
  geom_hline(yintercept = 1)+
  facet_wrap(vars(effects), scales = 'fixed', labeller =  effects_labels)+
  scale_color_manual(values = palette_seasons_4)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.text.x = element_text(size = 8), axis.text.x = element_text(angle = 0, size = 10), 
        legend.position = "right")

### Alphaproteobacteria----
relative_effects_alpha <- reg_all_slopes_chosen_silva_tax %>% 
  subset(season != "Early_fall") %>%
  subset(order != "Alphaproteobacteria") %>%
  distinct(treatment, season, asv_num, domain, phylum, class, order, family, genus, species, tax_ed, .keep_all = TRUE) %>%
  filter(pvalue_slope_chosen < 0.05) %>%
  group_by(family, treatment, season) %>%
  summarize(slope_chosen_days_mean = mean(slope_chosen_days)) %>%
  pivot_wider(names_from = treatment, values_from = c(slope_chosen_days_mean), 
              id_cols = c(season, family)) %>%
  as_tibble %>%
  mutate(across(!c(family, season), as.numeric)) %>%
  mutate(effect_top_down_virus= VL/DL,
         effect_top_down_grazers_dark = PD/CD,
         effect_top_down_grazers_light = PL/CL,
         effectbottom_up = DL/PL,
         effect_light_C = CL/CD,
         effect_light_P = PL/PD)

relative_effects_alpha_l <- relative_effects_f %>%
  pivot_longer(cols = starts_with('effect'),
               names_to = 'effects',
               values_to = 'ratio') %>%
  as_tibble()

relative_effects_alpha_l$season <-relative_effects_alpha_l$season %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))

relative_effects_alpha_l %>% 
  #filter(ratio > 10 | ratio < -10) %>%
  ggplot(aes(x = family, y = ratio, na.rm = TRUE))+
  #geom_tile(aes(fill = ratio), color = 'white', lwd = 1.5, linetype = 0.5)+
  #scale_fill_gradient2(low = '#BD4247', high =  '#42BDB8', na.value = '#fdfdfd', midpoint = 0)+
  #geom_boxplot()+
  geom_point(aes(color = season))+
  coord_flip()+
  #geom_boxplot()+
  labs( x = 'Alhaproteobacteria family', y = 'Response ratio', color = 'Season') +
  geom_hline(yintercept = 1)+
  facet_grid(vars(effects), scales = 'free', labeller =  effects_labels)+
  scale_color_manual(values = palette_seasons_4)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.text.x = element_text(size = 14), axis.text.x = element_text(angle = 90, size = 10), 
        legend.position = "right")

relative_effects_f_l %>% 
  filter(ratio > 10 | ratio < -10) %>%
  ggplot(aes(x = season, y = ratio, na.rm = TRUE))+
  #geom_tile(aes(fill = ratio), color = 'white', lwd = 1.5, linetype = 0.5)+
  #scale_fill_gradient2(low = '#BD4247', high =  '#42BDB8', na.value = '#fdfdfd', midpoint = 0)+
  geom_boxplot()+
  geom_point(aes(color = season))+
  #coord_flip()+
  #geom_boxplot()+
  labs( x = 'Season', y = 'Response ratio', color = 'Season') +
  #scale_y_continuous(breaks = c(-500, -200, -300,-50,  1, 50, 100,  300, 500))+
  geom_hline(yintercept = 1)+
  facet_wrap(vars(effects), scales = 'fixed', labeller =  effects_labels)+
  scale_color_manual(values = palette_seasons_4)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.text.x = element_text(size = 8), axis.text.x = element_text(angle = 0, size = 10), 
        legend.position = "right")

## Plotting relative effects at Family level-----
relative_effects_family_l %>%
  #filter(ratio > 1) %>%
  ggplot(aes(effects, ratio, color = phylum, na.rm = TRUE))+ #, na.rm = TRUE
  geom_boxplot(outlier.shape = NA)+
  scale_color_manual(values=palf_large(30))+
  scale_x_discrete(labels = effects_labels)+
  #scale_y_log10()+
  geom_hline(yintercept = 1)+
  scale_y_continuous(limits = quantile(na.rm = TRUE, relative_effects_family_l$ratio, c(0.1, 0.9)))+
  #facet_wrap(~effects, labeller =  effects, scales = 'free')+
  theme_bw()+
  theme(strip.text.x = element_text(size = 5),
        axis.text.x = element_text(angle = 0, size = 10), 
        legend.position = "right")

##Mean growth rates by families plots -----
##Alphaproteobacteria
relative_effects_family_l %>% 
  filter(class == 'Alphaproteobacteria') %>%
  #filter(ratio > 10 | ratio < -10) %>%
  ggplot(aes(x = order, y = ratio, na.rm = TRUE))+
  #geom_tile(aes(fill = ratio), color = 'white', lwd = 1.5, linetype = 0.5)+
  #scale_fill_gradient2(low = '#BD4247', high =  '#42BDB8', na.value = '#fdfdfd', midpoint = 0)+
  #geom_boxplot()+
  geom_point(aes(color = season))+
  coord_flip()+
  #geom_boxplot()+
  labs( x = 'Families', y = 'Response ratio', color = 'Season') +
  geom_hline(yintercept = 1)+
  facet_wrap(~effects, labeller =  effects_labels, scales = 'free')+
  #facet_grid(vars(effects), scales = 'free', labeller =  effects)+
  scale_color_manual(values = palette_seasons_4)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.text.x = element_text(size = 10), axis.text.x = element_text(angle = 90, size = 10), 
        legend.position = "right")

relative_effects_family_l %>% 
  filter(class == 'Gammaproteobacteria') %>%
  #filter(ratio > 10 | ratio < -10) %>%
  ggplot(aes(x = order, y = ratio, na.rm = TRUE))+
  #geom_tile(aes(fill = ratio), color = 'white', lwd = 1.5, linetype = 0.5)+
  #scale_fill_gradient2(low = '#BD4247', high =  '#42BDB8', na.value = '#fdfdfd', midpoint = 0)+
  #geom_boxplot()+
  geom_point(aes(color = season))+
  coord_flip()+
  #geom_boxplot()+
  labs( x = 'Families', y = 'Response ratio', color = 'Season') +
  geom_hline(yintercept = 1)+
  facet_wrap(~effects, labeller = effects_labels, scales = 'free')+
  #facet_wrap(~effects, labeller(multi_line = TRUE), scales = 'free')+
  #facet_grid(vars(effects), scales = 'free', labeller =  effects)+
  scale_color_manual(values = palette_seasons_4)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.text.x = element_text(size = 8), axis.text.x = element_text(angle = 90, size = 10), 
        legend.position = "right")

relative_effects_family_l %>% 
  filter(phylum == 'Actinobacteriota') %>%
  #filter(ratio > 10 | ratio < -10) %>%
  ggplot(aes(x = order, y = ratio, na.rm = TRUE))+
  #geom_tile(aes(fill = ratio), color = 'white', lwd = 1.5, linetype = 0.5)+
  #scale_fill_gradient2(low = '#BD4247', high =  '#42BDB8', na.value = '#fdfdfd', midpoint = 0)+
  #geom_boxplot()+
  geom_point(aes(color = season))+
  coord_flip()+
  #geom_boxplot()+
  labs( x = 'Families', y = 'Response ratio', color = 'Season') +
  geom_hline(yintercept = 1)+
  facet_wrap(~effects, labeller = effects_labels, scales = 'free')+
  #facet_wrap(~effects, labeller(multi_line = TRUE), scales = 'free')+
  #facet_grid(vars(effects), scales = 'free', labeller =  effects)+
  scale_color_manual(values = palette_seasons_4)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.text.x = element_text(size = 8), axis.text.x = element_text(angle = 90, size = 10), 
        legend.position = "right")

relative_effects_family_l %>% 
  filter(phylum == 'Bacteroidota') %>%
  #filter(ratio > 10 | ratio < -10) %>%
  ggplot(aes(x = order, y = ratio, na.rm = TRUE))+
  #geom_tile(aes(fill = ratio), color = 'white', lwd = 1.5, linetype = 0.5)+
  #scale_fill_gradient2(low = '#BD4247', high =  '#42BDB8', na.value = '#fdfdfd', midpoint = 0)+
  #geom_boxplot()+
  geom_point(aes(color = season))+
  coord_flip()+
  #geom_boxplot()+
  labs( x = 'Families', y = 'Response ratio', color = 'Season') +
  geom_hline(yintercept = 1)+
  facet_wrap(~effects, labeller = effects_labels, scales = 'free')+
  #facet_wrap(~effects, labeller(multi_line = TRUE), scales = 'free')+
  #facet_grid(vars(effects), scales = 'free', labeller =  effects)+
  scale_color_manual(values = palette_seasons_4)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.text.x = element_text(size = 8), axis.text.x = element_text(angle = 90, size = 10), 
        legend.position = "right")

relative_effects_family_l %>% 
  filter(phylum == 'Cyanobacteria') %>%
  #filter(ratio > 10 | ratio < -10) %>%
  ggplot(aes(x = order, y = ratio, na.rm = TRUE))+
  #geom_tile(aes(fill = ratio), color = 'white', lwd = 1.5, linetype = 0.5)+
  #scale_fill_gradient2(low = '#BD4247', high =  '#42BDB8', na.value = '#fdfdfd', midpoint = 0)+
  #geom_boxplot()+
  geom_point(aes(color = season))+
  coord_flip()+
  #geom_boxplot()+
  labs( x = 'Families', y = 'Response ratio', color = 'Season') +
  geom_hline(yintercept = 1)+
  facet_wrap(~effects, labeller = effects_labels, scales = 'free')+
  #facet_wrap(~effects, labeller(multi_line = TRUE), scales = 'free')+
  #facet_grid(vars(effects), scales = 'free', labeller =  effects)+
  scale_color_manual(values = palette_seasons_4)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.text.x = element_text(size = 8), axis.text.x = element_text(angle = 90, size = 10), 
        legend.position = "right")

relative_effects_family_l %>% 
  filter(phylum == 'Firmicutes') %>%
  #filter(ratio > 10 | ratio < -10) %>%
  ggplot(aes(x = order, y = ratio, na.rm = TRUE))+
  #geom_tile(aes(fill = ratio), color = 'white', lwd = 1.5, linetype = 0.5)+
  #scale_fill_gradient2(low = '#BD4247', high =  '#42BDB8', na.value = '#fdfdfd', midpoint = 0)+
  #geom_boxplot()+
  geom_point(aes(color = season))+
  coord_flip()+
  #geom_boxplot()+
  labs( x = 'Families', y = 'Response ratio', color = 'Season') +
  geom_hline(yintercept = 1)+
  facet_wrap(~effects, labeller = effects_labels, scales = 'free')+
  #facet_wrap(~effects, labeller(multi_line = TRUE), scales = 'free')+
  #facet_grid(vars(effects), scales = 'free', labeller =  effects)+
  scale_color_manual(values = palette_seasons_4)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.text.x = element_text(size = 8), axis.text.x = element_text(angle = 90, size = 10), 
        legend.position = "right")

relative_effects_family_l %>% 
  filter(phylum == 'Planctomycetota') %>%
  #filter(ratio > 10 | ratio < -10) %>%
  ggplot(aes(x = order, y = ratio, na.rm = TRUE))+
  #geom_tile(aes(fill = ratio), color = 'white', lwd = 1.5, linetype = 0.5)+
  #scale_fill_gradient2(low = '#BD4247', high =  '#42BDB8', na.value = '#fdfdfd', midpoint = 0)+
  #geom_boxplot()+
  geom_point(aes(color = season))+
  coord_flip()+
  geom_boxplot()+
  labs( x = 'Families', y = 'Response ratio', color = 'Season') +
  geom_hline(yintercept = 1)+
  facet_wrap(vars(effects), labeller = effects_labels, scales = 'free')+
  #facet_grid(vars(effects), scales = 'free', labeller =  effects)+
  scale_color_manual(values = palette_seasons_4)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.text.x = element_text(size = 8), axis.text.x = element_text(angle = 90, size = 10), 
        legend.position = "right")

relative_effects_family_l %>% 
  filter(phylum == 'Verrucomicrobiota') %>%
  #filter(ratio > 10 | ratio < -10) %>%
  ggplot(aes(x = order, y = ratio, na.rm = TRUE))+
  #geom_tile(aes(fill = ratio), color = 'white', lwd = 1.5, linetype = 0.5)+
  #scale_fill_gradient2(low = '#BD4247', high =  '#42BDB8', na.value = '#fdfdfd', midpoint = 0)+
  geom_boxplot()+
  geom_point(aes(color = season))+
  coord_flip()+
  #geom_boxplot()+
  labs( x = 'Families', y = 'Response ratio', color = 'Season') +
  geom_hline(yintercept = 1)+
  facet_wrap(~effects, labeller = effects_labels, scales = 'free')+
  #facet_wrap(~effects, labeller(multi_line = TRUE), scales = 'free')+
  #facet_grid(vars(effects), scales = 'free', labeller =  effects)+
  scale_color_manual(values = palette_seasons_4)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.text.x = element_text(size = 8), axis.text.x = element_text(angle = 90, size = 10), 
        legend.position = "right")

##Plotting effects at order taxonomic level-----
relative_effects_order %>%
  colnames()
relative_effects_order_l <- relative_effects_order %>%
  pivot_longer(cols = starts_with('effect'),
               names_to = 'effects',
               values_to = 'ratio') %>%
  as_tibble()

relative_effects_order_l %>%
  #filter(ratio > 1) %>%
  ggplot(aes(effects, ratio, color = phylum, na.rm = TRUE))+ #, na.rm = TRUE
  geom_boxplot(outlier.shape = NA)+
  scale_color_manual(values=palf_large(21))+
  scale_x_discrete(labels = effects_labels)+
  #scale_y_log10()+
  geom_hline(yintercept = 1)+
  scale_y_continuous(limits = quantile(na.rm = TRUE, relative_effects_family_l$ratio, c(0.1, 0.9)))+
  #facet_wrap(~effects, labeller =  effects, scales = 'free')+
  theme_bw()+
  theme(strip.text.x = element_text(size = 5),
        axis.text.x = element_text(angle = 0, size = 10), 
        legend.position = "right")

##Ploting top effects on families----
relative_effects_family_l %>%
  group_by(effects) %>%
  slice_max(ratio, n = 50) %>%
  ungroup() %>%
  ggplot(aes(x = family_f, y = ratio), na.rm = TRUE)+
  geom_point(aes(color = phylum), size = 3, alpha = 0.8)+
  scale_color_manual(values = palf_large(15))+
  coord_flip()+
  #geom_boxplot()+
  labs(x = 'Taxonomic family', y = 'Response ratio', shape = 'Effect', color = 'Phylum')+
  geom_hline(yintercept = 1)+
  #scale_y_continuous(limits = c(1,1000), breaks = c(1, 250, 500, 750, 1000))+
  #scale_y_discrete(labels = effects_labels)+
  facet_wrap(.~effects, labeller = effects_labels, scales = 'free')+
  #scale_shape_discrete(labels = effects_labels)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.text.x = element_text(size = 8), axis.text.x = element_text(angle = 90, size = 10), 
        legend.position = "right")

##color by class 
relative_effects_family_l %>%
  group_by(effects) %>%
  slice_max(ratio, n = 25) %>%
  ungroup() %>%
  ggplot(aes(x = family_f, y = ratio), na.rm = TRUE)+
  geom_point(aes(color = class), size = 2, alpha = 0.8)+
  scale_color_manual(values = palf_large(23))+
  coord_flip()+
  #geom_boxplot()+
  labs(x = 'Taxonomic family', y = 'Response ratio', shape = 'Effect', color = 'Class')+
  geom_hline(yintercept = 1)+
  #scale_y_continuous(limits = c(1,1000), breaks = c(1, 250, 500, 750, 1000))+
  #scale_y_discrete(labels = effects_labels)+
  facet_wrap(.~effects, labeller = effects_labels, scales = 'free_x')+
  #scale_shape_discrete(labels = effects_labels)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.text.x = element_text(size = 8), axis.text.x = element_text(angle = 90, size = 10), 
        axis.text.y = element_text(angle = 0, size = 8),
        legend.position = "right")

##Ploting top effects on orders-----
relative_effects_order_l %>%
  dim()

relative_effects_order_l %>%
  colnames()

relative_effects_order_l %>%
  glimpse()

relative_effects_order_l %>%
  group_by(effects) %>%
  slice_max(ratio, n = 25) %>%
  ungroup() %>%
  arrange(phylum, order) %>%
  ggplot(aes(x = order, y = ratio), na.rm = TRUE)+
  geom_point(aes(color = phylum), size = 3, alpha = 0.8)+
  scale_color_manual(values = palf_large(17))+
  coord_flip()+
  #geom_boxplot()+
  labs(x = 'Taxonomic order', y = 'Response ratio', shape = 'Effect', color = 'Phylum')+
  geom_hline(yintercept = 1)+
  #scale_y_continuous(limits = c(1,1000), breaks = c(1, 250, 500, 750, 1000))+
  #scale_y_discrete(labels = effects_labels)+
  facet_wrap(.~effects, labeller = effects_labels, scales = 'free')+
  #scale_shape_discrete(labels = effects_labels)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.text.x = element_text(size = 8), 
        axis.text.x = element_text(angle = 90, size = 10), 
        axis.text.y = element_text(angle = 0, size = 8),
        legend.position = "right")

##color by class 
relative_effects_order_l$order_f <-  factor(relative_effects_order_l$order_f, 
                                            levels=unique(relative_effects_order_l$order_f[order(relative_effects_order_l$class)]), 
                                            ordered=TRUE)

relative_effects_order_l %>%
  group_by(effects) %>%
  slice_max(ratio, n = 25) %>%
  ungroup() %>%
  arrange(phylum_f, order) %>%
  ggplot(aes(x = order_f, y = ratio), na.rm = TRUE)+
  geom_point(aes(color = class), size = 2, alpha = 0.8)+
  scale_color_manual(values = palf_large(26))+
  coord_flip()+
  #geom_boxplot()+
  labs(x = 'Taxonomic order', y = 'Response ratio', shape = 'Effect', color = 'Class')+
  geom_hline(yintercept = 1)+
  #scale_y_continuous(limits = c(1,1000), breaks = c(1, 250, 500, 750, 1000))+
  #scale_y_discrete(labels = effects_labels)+
  facet_wrap(.~effects, labeller = effects_labels, scales = 'free')+
  #guides(color = guide_legend(ncol=1))+
  #scale_shape_discrete(labels = effects_labels)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.text.x = element_text(size = 8), axis.text.x = element_text(angle = 90, size = 10), 
        axis.text.y = element_text(angle = 0, size = 8),
        legend.position = "right")

relative_effects_order_l %>%
  group_by(effects) %>%
  slice_max(ratio, n = 25) %>%
  ungroup() %>%
  ggplot(aes(x = order_f, y = ratio, shape = effects), na.rm = TRUE)+
  geom_point(aes(color = class), size = 2, alpha = 0.8)+
  scale_color_manual(values = palf_large(26))+
  coord_flip()+
  #geom_boxplot()+
  labs(x = 'Taxonomic order', y = 'Response ratio', shape = 'Effect', color = 'Class')+
  geom_hline(yintercept = 1)+
  guides(color = guide_legend(ncol=2))+
  #scale_y_continuous(limits = c(1,1000), breaks = c(1, 250, 500, 750, 1000))+
  #scale_y_discrete(labels = effects_labels)+
  #facet_wrap(.~effects, labeller = effects_labels, scales = 'free')+
  scale_shape_discrete(labels = effects_labels)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.text.x = element_text(size = 8), axis.text.x = element_text(angle = 90, size = 10), 
        axis.text.y = element_text(angle = 0, size = 8),
        legend.position = "right")

##For using class as axis but order as a group
#no s'ordenen els phylums per l'ordre correcte al de asv... descobrir perquè
effects_top25_order <- relative_effects_order_l %>%
  group_by(effects) %>%
  slice_max(ratio, n = 25) %>%
  ungroup() %>%
  ggplot(aes(x = class_f, y = ratio, color = phylum_f, group = order_f), na.rm = TRUE)+
  geom_point(aes(color = phylum_f), size = 2, alpha = 0.8)+
  scale_color_manual(values = palette_phylums)+
  coord_flip()+
  #geom_boxplot()+
  labs(x = 'Taxonomic class', y = 'Response ratio', shape = 'Effect', color = 'Phylum')+
  geom_hline(yintercept = 1)+
  #scale_y_continuous(limits = c(1,1000))+
  #scale_y_discrete(labels = effects_labels)+
  #facet_wrap(.~effects, labeller = effects_labels, scales = 'free')+
  facet_grid(.~effects, labeller = effects_labels, scales = 'free')+
  #scale_shape_discrete(labels = effects_labels)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.text.x = element_text(size = 8), axis.text.x = element_text(angle = 0, size = 10), 
        axis.text.y = element_text(angle = 0, size = 8),
        legend.text = element_text(size = 4),
        legend.position = "right",
        panel.spacing.x = unit(0.5,"line"))

effects_top25_asv <- relative_effects_l %>%
  group_by(effects) %>%
  slice_max(ratio, n = 25) %>%
  ungroup() %>%
  ggplot(aes(x = class_f, y = ratio, color = phylum_f, group = order_f), na.rm = TRUE)+
  geom_point(aes(color = phylum_f), size = 2, alpha = 0.8)+
  scale_color_manual(values = palette_phylums)+
  coord_flip()+
  #geom_boxplot()+
  labs(x = 'Taxonomic class', y = 'Response ratio', shape = 'Effect', color = 'Phylum')+
  geom_hline(yintercept = 1)+
  #scale_y_continuous(limits = c(1,1000))+
  #scale_y_discrete(labels = effects_labels)+
  #facet_wrap(.~effects, labeller = effects_labels, scales = 'free')+
  facet_grid(.~effects, labeller = effects_labels, scales = 'free')+
  #scale_shape_discrete(labels = effects_labels)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.text.x = element_text(size = 8), axis.text.x = element_text(angle = 0, size = 10), 
        axis.text.y = element_text(angle = 0, size = 8),
        legend.text = element_text(size = 4),
        legend.position = "right",
        panel.spacing.x = unit(0.5,"line"))

effects_top25_family <- relative_effects_family_l %>%
  group_by(effects) %>%
  slice_max(ratio, n = 25) %>%
  ungroup() %>%
  ggplot(aes(x = class_f, y = ratio, color = phylum_f, group = order_f), na.rm = TRUE)+
  geom_point(aes(color = phylum_f), size = 2, alpha = 0.8)+
  scale_color_manual(values = palette_phylums)+
  coord_flip()+
  #geom_boxplot()+
  labs(x = 'Taxonomic class', y = 'Response ratio', shape = 'Effect', color = 'Phylum')+
  geom_hline(yintercept = 1)+
  #scale_y_continuous(limits = c(1,1000))+
  #scale_y_discrete(labels = effects_labels)+
  #facet_wrap(.~effects, labeller = effects_labels, scales = 'free')+
  facet_grid(.~effects, labeller = effects_labels, scales = 'free')+
  #scale_shape_discrete(labels = effects_labels)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.text.x = element_text(size = 8), axis.text.x = element_text(angle = 0, size = 10), 
        axis.text.y = element_text(angle = 0, size = 8),
        legend.text = element_text(size = 4),
        legend.position = "right",
        panel.spacing.x = unit(0.5,"line"))

relative_effects_top25  <- multi_panel_figure(columns = 1, rows = 3, width = 250, height = 400, 
                                              row_spacing = 2, unit = 'mm',
                                              panel_label_type = 'none')

relative_effects_top25  %<>%
  fill_panel(effects_top25_asv, row = 1, col = 1) %<>%
  fill_panel(effects_top25_family, row = 2, col = 1) %<>%
  fill_panel(effects_top25_order, row = 3, col = 1)

ggsave('relative_effects_top25_2.pdf', relative_effects_top25, 
       path = "~/Documentos/Doctorat/REMEI/results/figures/raw/",
       width = 250,
       height = 400,
       units = 'mm')

##Ploting null efects (aprox 1)----
relative_effects_order_l$order_f <-  factor(relative_effects_order_l$order_f, 
                                            levels=unique(relative_effects_order_l$order_f[order(relative_effects_order_l$class)]), 
                                            ordered=TRUE)

relative_effects_order_l %>%
  group_by(effects) %>%
  #filter(between(ratio, 0.8, 1.2)) %>%
  ungroup() %>%
  ggplot(aes(x = order_f, y = ratio), na.rm = TRUE)+
  geom_point(aes(color = phylum), size = 2, alpha = 0.8)+
  scale_color_manual(values = palf_large(20))+
  coord_flip()+
  #geom_boxplot()+
  labs(x = 'Taxonomic order', y = 'Response ratio', shape = 'Effect', color = 'Class')+
  geom_hline(yintercept = 1)+
  #scale_y_continuous(limits = c(1,1000), breaks = c(1, 250, 500, 750, 1000))+
  #scale_y_discrete(labels = effects_labels)+
  facet_grid(.~effects, labeller = effects_labels, scales = 'free', switch = 'x', space = 'free')+
  #facet_wrap(.~effects, labeller = effects_labels, scales = 'free', strip.position = 'bottom')+
  #scale_shape_discrete(labels = effects_labels)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.text.x = element_text(size = 8), 
        axis.text.x = element_text(angle = 90, size = 10), 
        axis.text.y = element_text(angle = 0, size = 8),
        legend.position = "right")

###General heat map for all the effects on different taxonomic groups-----
library(scales)
relative_effects_phylum_l %>%
  ggplot(aes(effects, phylum_f, fill=ratio, group = phylum_f, na.rm = TRUE))+
  geom_tile(aes(fill = ratio), color = 'white', lwd = 1.5, 
            linetype = 0.5)+
  #scale_fill_gradient2(low = '#BD4247', high =  '#42BDB8', na.value = '#fdfdfd', midpoint = 1)+
  #scale_fill_gradientn(colors = c("#bd4247","#42bdb8","#fdfdfd","#84bd42","#7a42bd"))+
  scale_fill_gradientn(colors = c("#bd4247","#42bdb8", "#84bd42","#7a42bd"), 
                       na.value = '#fdfdfd',  values = rescale(c(-200,1,100)))+
  #geom_point(aes(color = phylum), size = 2, alpha = 0.8)+
  #scale_color_manual(values = palf_large(24))+
  #coord_flip()+
  #geom_boxplot()+
  labs(y = 'Taxonomic order', x = 'Response ratio', shape = 'Effect', color = 'Class')+
  #geom_hline(yintercept = 1)+
  #scale_y_continuous(limits = c(1,1000), breaks = c(1, 250, 500, 750, 1000))+
  scale_x_discrete(labels = effects_labels)+
  #facet_wrap(.~effects, labeller = effects_labels, scales = 'free')+
  #scale_shape_discrete(labels = effects_labels)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.text.x = element_text(size = 8), 
        axis.text.x = element_text(angle = 0, size = 8), 
        axis.text.y = element_text(angle = 0, size = 8),
        legend.position = "right")

relative_effects_order_l %>%
  filter(effects != is.na(effects)) %>% 
  ggplot(aes(effects, order_f, fill=ratio, group = phylum_f, na.rm = TRUE))+
  geom_tile(aes(fill = ratio), color = 'white', lwd = 1.5, 
            linetype = 0.5)+
  #scale_fill_gradient2(low = '#BD4247', high =  '#42BDB8', na.value = '#fdfdfd', midpoint = 1)+
  #scale_fill_gradientn(colors = c("#bd4247","#42bdb8","#fdfdfd","#84bd42","#7a42bd"))+
  scale_fill_gradientn(colors = c("#bd4247","#42bdb8"), 
                       na.value = '#fdfdfd',  values = rescale(c(-200,1,100)))+
  #geom_point(aes(color = phylum), size = 2, alpha = 0.8)+
  #scale_color_manual(values = palf_large(24))+
  #coord_flip()+
  #geom_boxplot()+
  labs(y = 'Taxonomic order', x = 'Response ratio', shape = 'Effect', color = 'Class')+
  #geom_hline(yintercept = 1)+
  #scale_y_continuous(limits = c(1,1000), breaks = c(1, 250, 500, 750, 1000))+
  scale_x_discrete(labels = effects_labels)+
  #facet_wrap(.~effects, labeller = effects_labels, scales = 'free')+
  #scale_shape_discrete(labels = effects_labels)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.text.x = element_text(size = 8), 
        axis.text.x = element_text(angle = 0, size = 8), 
        axis.text.y = element_text(angle = 0, size = 8),
        legend.position = "right")


relative_effects_family_l %>%
  ggplot(aes(effects, family_f, fill=ratio, na.rm = TRUE))+
  geom_tile(aes(fill = ratio), color = 'white', lwd = 1.5, linetype = 0.5)+
  #scale_fill_gradient2(low = '#BD4247', high =  '#42BDB8', na.value = '#fdfdfd', midpoint = 0)+
  scale_fill_gradientn(colors = c("#bd4247","#42bdb8"), 
                       na.value = '#fdfdfd',  values = rescale(c(-200,1,100)))+
  #geom_point(aes(color = phylum), size = 2, alpha = 0.8)+
  #scale_color_manual(values = palf_large(24))+
  #coord_flip()+
  #geom_boxplot()+
  labs(y = 'Taxonomic family', fill = 'Response ratio', x = 'Effect', color = 'Class')+
  #geom_hline(yintercept = 1)+
  #scale_y_continuous(limits = c(1,1000), breaks = c(1, 250, 500, 750, 1000))+
  scale_x_discrete(labels = effects_labels)+
  #facet_wrap(.~effects, labeller = effects_labels, scales = 'free')+
  #scale_shape_discrete(labels = effects_labels)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.text.x = element_text(size = 8), 
        axis.text.x = element_text(angle = 0, size = 8), 
        axis.text.y = element_text(angle = 0, size = 8),
        legend.position = "right")

asv_eff<- relative_effects_l %>%
  ggplot(aes(effects, asv_num, fill=ratio, na.rm = TRUE))+
  geom_tile(aes(fill = ratio), color = 'white', lwd = 1.5, linetype = 0.5)+
  scale_fill_gradient2(low = '#BD4247', high =  '#42BDB8', na.value = '#fdfdfd', midpoint = 0)+
  #geom_point(aes(color = phylum), size = 2, alpha = 0.8)+
  #scale_color_manual(values = palf_large(24))+
  #coord_flip()+
  #geom_boxplot()+
  labs(y = 'Taxonomic family', fill = 'Response ratio', x = 'Effect', color = 'Class')+
  #geom_hline(yintercept = 1)+
  #scale_y_continuous(limits = c(1,1000), breaks = c(1, 250, 500, 750, 1000))+
  scale_x_discrete(labels = effects_labels)+
  #facet_wrap(.~effects, labeller = effects_labels, scales = 'free')+
  #scale_shape_discrete(labels = effects_labels)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.text.x = element_text(size = 8), 
        axis.text.x = element_text(angle = 0, size = 8), 
        axis.text.y = element_text(angle = 0, size = 8),
        legend.position = "right")

#grid.arrange(order_eff, family_eff, asv_eff)

###By Seasons like Sanchez et al., 2020------
##comprobo que tingui el mateix nº de classes a tots els datasets per tenir els 
##mateixos colors
relative_effects_order_l %$% 
  class %>%
  unique()

relative_effects_family_l %$% 
  class %>%
  unique()

relative_effects_l %$% 
  class %>%
  unique()

effects_seasons_family<- relative_effects_family_l %>%
  ggplot(aes(season, ratio, na.rm = TRUE))+ #, na.rm = TRUE
  geom_point(aes(color = class), position = position_dodge(0.3))+
  geom_boxplot(alpha = 0.3)+
  scale_color_manual(values = palf_large(47))+
  facet_wrap(~effects, labeller =  effects_labels, scales = 'fixed', ncol = 1)+
  labs(x = 'Season', y = 'Response ratio', color = 'Taxonomic class')+
  geom_hline(yintercept = 1)+
  scale_y_continuous( limits = quantile(na.rm = TRUE, relative_effects_family_l$ratio, c(0.1, 0.9)))+
  theme_bw()+
  theme(strip.text.x = element_text(size = 5), 
        axis.text.x = element_text(angle = 0, size = 10), 
        legend.position = "none")

effects_seasons_asv <- relative_effects_l %>%
  ggplot(aes(season, ratio, na.rm = TRUE))+ #, na.rm = TRUE
  geom_point(aes(color = class), position = position_dodge(0.3))+
  geom_boxplot(alpha = 0.3)+
  scale_color_manual(values = palf_large(47))+
  facet_wrap(~effects, labeller =  effects_labels, scales = 'fixed', ncol = 1)+
  labs(x = 'Season', y = 'Response ratio', color = 'Taxonomic class')+
  geom_hline(yintercept = 1)+
  scale_y_continuous( limits = quantile(na.rm = TRUE, relative_effects_family_l$ratio, c(0.1, 0.9)))+
  theme_bw()+
  theme(strip.text.x = element_text(size = 5), 
        axis.text.x = element_text(angle = 0, size = 10), 
        legend.position = "none")

effects_seasons_order <- relative_effects_order_l %>%
  ggplot(aes(season, ratio, na.rm = TRUE))+ #, na.rm = TRUE
  geom_point(aes(color = class), position = position_dodge(0.3))+
  geom_boxplot(alpha = 0.3)+
  scale_color_manual(values = palf_large(47))+
  facet_wrap(~effects, labeller =  effects_labels, scales = 'fixed', ncol = 1)+
  labs(x = 'Season', y = 'Response ratio', color = 'Taxonomic\nclass')+
  geom_hline(yintercept = 1)+
  guides(color = guide_legend(ncol=2))+
  scale_y_continuous(limits = quantile(na.rm = TRUE, relative_effects_family_l$ratio, c(0.1, 0.9)))+
  theme_bw()+
  theme(strip.text.x = element_text(size = 5), 
        axis.text.x = element_text(angle = 0, size = 10), 
        legend.position = "right",
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12))

grid.arrange(effects_seasons_asv, 
             effects_seasons_family, 
             effects_seasons_order, ncol = 3)

relative_effects_seasons  <- multi_panel_figure(columns = 4, rows = 1, width = 410, height = 200, 
                                                col_spacing = 2, unit = 'mm',
                                                panel_label_type = 'none')

relative_effects_seasons  %<>%
  fill_panel(effects_seasons_asv, row = 1, col = 1) %<>%
  fill_panel(effects_seasons_family, row = 1, col = 2) %<>%
  fill_panel(effects_seasons_order, row = 1, col = 3:4)

ggsave('relative_effects_seasons_2.pdf', relative_effects_seasons, 
       path = "~/Documentos/Doctorat/REMEI/results/figures/raw/",
       width = 410,
       height = 200,
       units = 'mm')

### Violin plots for relative effects----
##For using class as axis but order as a group
#no s'ordenen els phylums per l'ordre correcte al de asv... descobrir perquè
relative_effects_order_l %>%
  group_by(effects) %>%
  #slice_max(ratio, n = 25) %>%
  ungroup() %>%
  ggplot(aes(x = class_f, y = ratio, color = phylum_f, group = order_f), na.rm = TRUE)+
  geom_point(aes(color = phylum_f), size = 2, alpha = 0.8)+
  scale_color_manual(values = palf_large(136))+
  coord_flip()+
  #geom_boxplot()+
  geom_violin()+
  labs(x = 'Taxonomic class', y = 'Response ratio', shape = 'Effect', color = 'Phylum')+
  geom_hline(yintercept = 1)+
  #scale_y_continuous(limits = c(1,1000))+
  #scale_y_discrete(labels = effects_labels)+
  #facet_wrap(.~effects, labeller = effects_labels, scales = 'free')+
  facet_grid(.~effects, labeller = effects_labels, scales = 'free')+
  #scale_shape_discrete(labels = effects_labels)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.text.x = element_text(size = 8), axis.text.x = element_text(angle = 0, size = 10), 
        axis.text.y = element_text(angle = 0, size = 8),
        legend.text = element_text(size = 4),
        legend.position = "right",
        panel.spacing.x = unit(0.5,"line"))

relative_effects_l %>%
  #group_by(effects) %>%
  #slice_max(ratio, n = 25) %>%
  ungroup() %>%
  ggplot(aes(x = class_f, y = ratio, color = phylum_f, group = order_f), na.rm = TRUE)+
  geom_point(aes(color = phylum_f), size = 2, alpha = 0.8)+
  scale_color_manual(values = palf_large(136))+
  coord_flip()+
  geom_violin()+
  #geom_boxplot()+
  labs(x = 'Taxonomic class', y = 'Response ratio', shape = 'Effect', color = 'Phylum')+
  geom_hline(yintercept = 1)+
  #scale_y_continuous(limits = c(1,1000))+
  #scale_y_discrete(labels = effects_labels)+
  #facet_wrap(.~effects, labeller = effects_labels, scales = 'free')+
  facet_grid(.~effects, labeller = effects_labels, scales = 'free')+
  #scale_shape_discrete(labels = effects_labels)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.text.x = element_text(size = 8), axis.text.x = element_text(angle = 0, size = 10), 
        axis.text.y = element_text(angle = 0, size = 8),
        legend.text = element_text(size = 4),
        legend.position = "right",
        panel.spacing.x = unit(0.5,"line"))

relative_effects_family_l %>%
  ggplot(aes(x = class_f, y = ratio, color = phylum_f), na.rm = TRUE)+
  geom_point(aes(color = phylum_f), size = 2, alpha = 0.8)+
  scale_color_manual(values = palf_large(136))+
  coord_flip()+
  #geom_boxplot()+
  geom_violin(aes(group = class_f))+
  #geom_boxplot(aes(group = class_f))+
  labs(x = 'Taxonomic class', y = 'Response ratio', shape = 'Effect', color = 'Phylum')+
  geom_hline(yintercept = 1)+
  #scale_y_continuous(limits = c(1,1000))+
  #scale_y_discrete(labels = effects_labels)+
  #facet_wrap(.~effects, labeller = effects_labels, scales = 'free')+
  facet_grid(.~effects, labeller = effects_labels, scales = 'free')+
  #scale_shape_discrete(labels = effects_labels)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.text.x = element_text(size = 8), axis.text.x = element_text(angle = 0, size = 10), 
        axis.text.y = element_text(angle = 0, size = 8),
        legend.text = element_text(size = 4),
        legend.position = "right",
        panel.spacing.x = unit(0.5,"line"))

relative_effects_top25  <- multi_panel_figure(columns = 1, rows = 3, width = 250, height = 400, 
                                              row_spacing = 2, unit = 'mm',
                                              panel_label_type = 'none')

relative_effects_top25  %<>%
  fill_panel(effects_top25_asv, row = 1, col = 1) %<>%
  fill_panel(effects_top25_family, row = 2, col = 1) %<>%
  fill_panel(effects_top25_order, row = 3, col = 1)

ggsave('relative_effects_top25.pdf', relative_effects_top25, 
       path = "~/Documentos/Doctorat/REMEI/results/figures/raw/",
       width = 250,
       height = 400,
       units = 'mm')

### Looking at control effects without calculating relative effects
reg_all_slopes_chosen_silva_tax$treatment <- reg_all_slopes_chosen_silva_tax$treatment %>% 
  as_factor() %>%
  factor(levels=(c("CL", "CD", "PL", "PD", "DL", "VL")))

reg_all_slopes_chosen_silva_tax$season <- reg_all_slopes_chosen_silva_tax$season %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))

reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05,
         slope_chosen_days != is.na(slope_chosen_days),
         slope_chosen_days > 0) %>%
  group_by(asv_num, treatment) %>%
  summarize(count = n()) %>%
  filter(count > 4)

reg_all_slopes_chosen_silva_tax %>%
  subset(asv_num == c('asv10',
         'asv55',
         'asv958')) %>%
  ggplot(aes(treatment, slope_chosen_days))+
  #geom_boxplot()+
  geom_point()+
  facet_wrap(vars(asv_num))+
  theme_bw()

reg_all_slopes_chosen_silva_tax %>%
  subset(asv_num == c('asv9', 'asv37', 'asv91')) %>%
  ggplot(aes(treatment, slope_chosen_days))+
  #geom_boxplot()+
  geom_point()+
  facet_wrap(vars(asv_num))+
  theme_bw()

reg_all_slopes_chosen_silva_tax %>%
  subset(asv_num == c('asv75', 'asv301')) %>%
  ggplot(aes(treatment, slope_chosen_days))+
  #geom_boxplot()+
  geom_point()+
  facet_grid(vars(asv_num))+
  theme_bw()

reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05,
         slope_chosen_days > 0) %>%
  group_by(family) %>%
  filter(n() > 10) %>%
  ggplot(aes(treatment, slope_chosen_days))+
  geom_point(aes(color = season))+
  geom_violin(alpha = 0.5)+
  facet_wrap(vars(family))+ #labeller = interaction(vars(class_f), vars(family_f), sep = ' ')
  scale_color_manual(values = palette_seasons_4)+
  theme_bw()+
  theme(legend.position = 'none')
  

##where is your maximal GR?
reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05,
         slope_chosen_days > 0) %>%
  select(treatment, asv_num, slope_chosen_days, phylum_f, class_f, order_f, family_f) %>%
  group_by(treatment, asv_num, phylum_f, class_f, order_f, family_f) %>%
  filter(n() > 4) %>%
  summarize(mean = mean(slope_chosen_days),
            sd = sd(slope_chosen_days)) %>%
  # pivot_wider(names_from = treatment, values_from = slope_chosen_days, values_fn = mean) %>%
  # pivot_longer(names_to = 'slope_chosen_days')
  ggplot(aes(interaction(fct_infreq(asv_num), family_f, sep = '\n'), mean, color = phylum_f))+
  geom_point()+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),
                size=.3,    # Thinner lines
                width=.2, color="grey")+
  facet_grid(vars(treatment))+
  scale_color_manual(values = palette_phylums)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))


##max gr table for treatments
reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05,
         slope_chosen_days > 0) %>%
  select(treatment, asv_num, slope_chosen_days, phylum_f, class_f, order_f, family_f, season) %>%
  group_by(treatment, asv_num, phylum_f, class_f, order_f, family_f, season) %>%
  summarize(mean = mean(slope_chosen_days),
            sd = sd(slope_chosen_days)) %>%
  group_by(treatment, season) %>%
  slice_max(mean, n = 10)  %>%
  # pivot_wider(names_from = treatment, values_from = slope_chosen_days, values_fn = mean) %>%
  # pivot_longer(names_to = 'slope_chosen_days')
  ggplot(aes(interaction(fct_infreq(family_f), order_f, sep = '\n'), mean, color = phylum_f, shape = treatment))+
  geom_point()+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),
                size=.3,    # Thinner lines
                width=.2, color="grey")+
  #facet_grid(vars(season))+
  scale_color_manual(values = palette_phylums)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))


##Bars from the center
relative_effects_family_l %>%
  colnames()
relative_effects_family_l %>%
  group_by( phylum_f, class_f, order_f, family_f, ratio) %>%
  filter(n() > 4) %>%
  ggplot(aes(ratio, family_f))+
  geom_point()+
  #facet_grid(. ~ vars(effects))+
  theme_bw()

##spider chart


##geom coordinate plot
### Reordering factors for ploting all datasets -----
relative_effects_phylum$season <-relative_effects_phylum$season %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))

relative_effects_class$season <-relative_effects_class$season %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))

relative_effects_order$season <-relative_effects_order$season %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))

relative_effects_family$season <-relative_effects_family$season %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))

relative_effects$season <-relative_effects$season %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))

relative_effects_phylum$effects <- relative_effects_phylum$effects %>% 
  factor(levels=(c("effect_bottom_up", "effect_top_down_grazers_dark", 
                   "effect_top_down_grazersight", "effect_top_down_virus",
                   "effectight_C", "effectight_P")))

relative_effects_class$effects <- relative_effects_class$effects %>% 
  factor(levels=(c("effect_bottom_up", "effect_top_down_grazers_dark", 
                   "effect_top_down_grazersight", "effect_top_down_virus",
                   "effectight_C", "effectight_P")))

relative_effects_order$effects <- relative_effects_order$effects %>% 
  factor(levels=(c("effect_bottom_up", "effect_top_down_grazers_dark", 
                   "effect_top_down_grazersight", "effect_top_down_virus",
                   "effectight_C", "effectight_P")))

relative_effects_family$effects <- relative_effects_family$effects %>% 
  factor(levels=(c("effect_bottom_up", "effect_top_down_grazers_dark", 
                   "effect_top_down_grazersight", "effect_top_down_virus",
                   "effectight_C", "effectight_P")))

relative_effects$effects <- relative_effects$effects %>% 
  factor(levels=(c("effect_bottom_up", "effect_top_down_grazers_dark", 
                   "effect_top_down_grazersight", "effect_top_down_virus",
                   "effectight_C", "effectight_P")))

library(GGally)
effects_difference_family %>%
  colnames()

ggparcoord(data = effects_difference_family, 
           columns = 19:24,
           groupColumn = 'family',
           showPoints = TRUE,
           boxplot = TRUE,
           alphaLines = 0.8,
           centerObsID = 1)+
  #facet_wrap(~phylum)+
  labs(x = 'Effects', y = 'Difference', color = 'Phylum')+
  scale_color_manual(values = palf_large_phylums(89))+
  #facet_grid(vars(phylum))+
  scale_x_discrete(labels = effects_labels)+
  theme_bw()+
  theme(legend.position = 'none')

relative_effects_family %>%
  colnames()

ggparcoord(data = relative_effects_family, 
           columns = 19:24,
           groupColumn = 'family',
           showPoints = TRUE,
           boxplot = TRUE,
           alphaLines = 0.8,
           centerObsID = 1)+
  facet_wrap(~phylum)+
  labs(x = 'Effects', y = 'Response ratio', color = 'Phylum')+
  scale_color_manual(values = palf_large_phylums(89))+
  #facet_grid(vars(phylum))+
  scale_x_discrete(labels = effects_labels)+
  theme_bw()+
  theme(legend.position = 'none')

relative_effects_class %>%
  colnames()

ggparcoord(data = relative_effects_class, 
           columns = 17:22,
           groupColumn = 'class',
           showPoints = TRUE,
           boxplot = FALSE,
           alphaLines = 0.8)+
  #facet_wrap(~phylum)+
  labs(x = 'Effects', y = 'Response ratio')+
  scale_color_manual(values = palette_large)+
  scale_x_discrete(labels = effects_labels)+
  theme_bw()

ggparcoord(data = relative_effects_order, 
           columns = 17:22,
           groupColumn = 'order',
           showPoints = TRUE,
           boxplot = TRUE,
           alphaLines = 0.8,
           centerObsID = 1)+
  #facet_wrap(~phylum)+
  labs(x = 'Effects', y = 'Response ratio')+
  #scale_color_manual(values = palette_phylums)+
  scale_x_discrete(labels = effects_labels)+
  theme_bw()+
  theme(legend.position = 'none')

relative_effects_phylum %>%
  colnames()

ggparcoord(data = relative_effects_phylum, 
           columns = 9:14,
           groupColumn = 'phylum',
           showPoints = TRUE,
           boxplot = TRUE,
           alphaLines = 0.8,
           centerObsID = 1)+
  facet_wrap(~phylum)+
  labs(x = 'Effects', y = 'Response ratio')+
  scale_color_manual(values = palette_phylums)+
  scale_x_discrete(labels = effects_labels)+
  theme_bw()

##
relative_effects_family %>%
  colnames()



#OPCIÓ 2: PROVEM DE FER RESTA ENTRE GR COMPTES DE DIVISIÓ (ACCEPTED) ------
effects_difference_family <- reg_all_slopes_chosen_silva_tax %>% 
  distinct(treatment, season, asv_num, domain, phylum, class, order, genus, species, tax_ed, asv_num, family, .keep_all = TRUE) %>%
  filter(pvalue_slope_chosen < 0.05 &
           slope_chosen_days > 0) %>% #només considero les growth rates positives
  group_by(treatment, season, domain, phylum, class, order, family) %>%
  summarise(slope_chosen_days_mean = mean(slope_chosen_days),
            slope_chosen_mean = mean(slope_chosen),
            slope_chosen_days_sd = sd(slope_chosen_days),
           # counts = add_count(),
            na.rm = TRUE) %>%
  pivot_wider(names_from = treatment, values_from = c(slope_chosen_days_mean, slope_chosen_days_sd), 
              id_cols = c(season, domain, phylum, class, order, family)) %>%
  as_tibble() %>%
  mutate(across(!c(domain, phylum, class, order, family, season), as.numeric)) %>%
  mutate(effect_top_down_virus = slope_chosen_days_mean_VL - slope_chosen_days_mean_DL,
         effect_top_down_grazers_dark = slope_chosen_days_mean_PD - slope_chosen_days_mean_CD,
         effect_top_down_grazers_light = slope_chosen_days_mean_PL - slope_chosen_days_mean_CL,
         effect_bottom_up = slope_chosen_days_mean_DL - slope_chosen_days_mean_PL,
         effect_light_C = slope_chosen_days_mean_CL - slope_chosen_days_mean_CD,
         effect_light_P = slope_chosen_days_mean_PL - slope_chosen_days_mean_PD)


# effects_difference_family %>%
#   colnames()
# 
# effects_difference_family_table <- effects_difference_family %>%
#   select(season, phylum, class, family, matches('mean'), matches('sd')) %>%
#   filter(slope_chosen_days_mean_CD != is.na(slope_chosen_days_mean_CD) &
#            slope_chosen_days_mean_CL != is.na(slope_chosen_days_mean_CL) &
#            slope_chosen_days_mean_PD != is.na(slope_chosen_days_mean_PD) &
#            slope_chosen_days_mean_PL != is.na(slope_chosen_days_mean_PL) &
#            slope_chosen_days_mean_DL != is.na(slope_chosen_days_mean_DL) &
#            slope_chosen_days_mean_VL != is.na(slope_chosen_days_mean_VL)) %>%
#   arrange(family)
# effects_difference_family_table %>%
#   colnames()
# 
# effects_difference_family_table <- effects_difference_family_table %>%
#   rename( 'Season' = season,
#           'Phylum' = phylum,
#           'Class' = class,
#           'Order' = order,
#           'Family' = family,
#     'Mean growth rate CD' = slope_chosen_days_mean_CD,
#           'Mean growth rate CL' = slope_chosen_days_mean_CL,
#           'Mean growth rate PL' = slope_chosen_days_mean_PL,
#           'Mean growth rate PD'= slope_chosen_days_mean_PD,
#           'Mean growth rate DL' = slope_chosen_days_mean_DL,
#           'Mean growth rate VL'= slope_chosen_days_mean_VL,
#           'sd growth rate CD' = slope_chosen_days_sd_CD,
#           'sd growth rate CL' = slope_chosen_days_sd_CL,
#           'sd growth rate PL' = slope_chosen_days_sd_PL,
#           'sd growth rate PD' = slope_chosen_days_sd_PD,
#           'sd growth rate DL' = slope_chosen_days_sd_DL,
#           'sd growth rate VL' = slope_chosen_days_sd_VL)
# 
# write.table(effects_difference_family_table, 'results/tables/mean_sd_gr_treatment_family.txt', sep = '\t')


effects_difference_family_l <- effects_difference_family %>%
  pivot_longer(cols = starts_with('effect'),
               names_to = 'effects',
               values_to = 'difference') %>%
  as_tibble()

effects_difference_family_l %>%
  colnames()

effects_difference_family_l$season <-effects_difference_family_l$season %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))

effects_difference_family_l$effects <- effects_difference_family_l$effects %>% 
  factor(levels=(c("effect_bottom_up", "effect_top_down_grazers_dark", 
                   "effect_top_down_grazers_light", "effect_top_down_virus",
                   "effect_light_C", "effect_light_P")))
effects_difference_family_l <- effects_difference_family_l %>%
  mutate(phylum_f = as_factor(phylum),
         family_f = as_factor(family),
         order_f = as_factor(order),
         class_f = as_factor(class))

effects_difference_family_l$class_f <-  factor(effects_difference_family_l$class_f, 
                                            levels=unique(effects_difference_family_l$class_f[order(effects_difference_family_l$phylum_f)]), 
                                            ordered=TRUE)

effects_difference_family_l$order_f <-  factor(effects_difference_family_l$order_f, 
                                            levels=unique(effects_difference_family_l$order_f[order(effects_difference_family_l$phylum_f,
                                                                                                    effects_difference_family_l$class_f)]), 
                                            ordered=TRUE)

effects_difference_family_l$family_f <-  factor(effects_difference_family_l$family_f, 
                                              levels=unique(effects_difference_family_l$family_f[order(effects_difference_family_l$phylum_f,
                                                                                                       effects_difference_family_l$class_f,
                                                                                                       effects_difference_family_l$order_f)]), 
                                              ordered=TRUE)
#interaction(class_f,phylum_f, sep = ' ')

# difference_gr_treatments_family <- effects_difference_family_l %>%
#   filter(difference != is.na(difference)) %>%
#   group_by(class) %>%
#   filter(n() > 10) %>%
#   group_by(family) %>%
#   filter(n() > 20) %>%
#   group_by(family, difference) %>%
#   mutate(counts = n()) %>%
#   ungroup() %>%
#   filter(difference > 0) %>%
#   ggplot(aes(season, effects, size = difference, color = phylum_f, alpha = counts))+
#   geom_point()+
#   scale_y_discrete(labels = effects_labels2)+
#   scale_size(range = c(0, 6), name="Growth rate difference\n between treatments")+
#   scale_color_manual(values = palette_phylums)+
#   scale_alpha(range = c(1,1), limits = c(0.1, 1))+
#   facet_wrap(vars(family_f))+ #, labeller = interaction(vars(class_f), vars(family_f), sep = ' ')
#   labs(y = 'Treatments growth rates differences', x = 'Season', color = 'Phylum')+
#   theme_bw()+
#   theme(legend.position = 'right', panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         axis.text.x = element_text(angle = 90), strip.text = element_text(size = 5), strip.background = element_blank())
# 
# ggsave('difference_gr_treatments_family3.pdf', difference_gr_treatments_family, 
#        path = "~/Documentos/Doctorat/REMEI/results/figures/raw/",
#        width = 240,
#        height = 180,
#        units = 'mm')

##color and size by difference in GR
palete_gradient <- c("#5e0000",
                     "#b24d5d",
                     '#FFFFFF',
                      "#4db2a2",
                      "#005a47") 
                     
difference_gr_treatments_family <- 
  effects_difference_family_l %>%
  filter(difference != is.na(difference)) %>%
  # group_by(class) %>%
  # filter(n() > 10) %>%
  group_by(family) %>%
  filter(n() > 10) %>%
  group_by(family, difference) %>%
  mutate(counts = n()) %>%
  # ungroup() %>%
  # filter(difference > 0) %>%
  ggplot(aes(season, effects, color = difference))+ #, alpha = counts
  geom_point(aes(color = difference, size = difference))+
  scale_y_discrete(labels = effects_labels2)+
  scale_size(range = c(0, 10), name ="Growth rate difference\n between treatments")+
  scale_colour_gradientn(colours = palete_gradient)+
   # scale_color_manual(values = palette_phylums)+
    #scale_color_gradient2(low = '#B24D5D', high =  '#4DB2A2',  midpoint = 0)+
  #scale_alpha(range = c(1,1), limits = c(1, 1))+
  facet_wrap(vars(family_f))+ #, labeller = interaction(vars(class_f), vars(family_f), sep = ' ')
  labs(y = 'Treatments growth rates differences', x = 'Season', color = '')+
  theme_bw()+
  theme(legend.position = 'right', panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 60, size = 5, hjust = 1), 
        axis.text.y = element_text(size = 5), legend.title = element_text(size = 7), axis.title = element_text(size = 7),
        strip.text = element_text(size = 7), strip.background = element_blank(), legend.text = element_text(size = 5))

ggsave('difference_gr_treatments_family4.pdf', difference_gr_treatments_family, 
       path = "~/Documentos/Doctorat/REMEI/results/figures/raw/",
       width = 260,
       height = 220,
       units = 'mm')

##order level bubble plot-----
effects_difference_order <- reg_all_slopes_chosen_silva_tax %>% 
  distinct(treatment, season, asv_num, domain, phylum, class, order, genus, species, tax_ed, asv_num, family, .keep_all = TRUE) %>%
  filter(pvalue_slope_chosen < 0.05 &
           slope_chosen_days > 0) %>%
  group_by(treatment, season, domain, phylum, class, order) %>%
  summarise(slope_chosen_days_mean = mean(slope_chosen_days),
            #slope_chosen_mean = mean(slope_chosen),
            slope_chosen_sd = sd(slope_chosen),
            na.rm = TRUE) %>%
  pivot_wider(names_from = treatment, values_from = c(slope_chosen_days_mean, slope_chosen_sd), #, slope_chosen_sd (per crear la taula)
              id_cols = c(season, domain, phylum, class, order)) %>%
  as_tibble() %>%
  mutate(across(!c(domain, phylum, class, order, season), as.numeric)) %>%
  mutate(effect_top_down_virus = slope_chosen_days_mean_VL - slope_chosen_days_mean_DL,
         effect_top_down_grazers_dark = slope_chosen_days_mean_PD - slope_chosen_days_mean_CD,
         effect_top_down_grazers_light = slope_chosen_days_mean_PL - slope_chosen_days_mean_CL,
         effect_bottom_up = slope_chosen_days_mean_DL - slope_chosen_days_mean_PL,
         effect_light_C = slope_chosen_days_mean_CL - slope_chosen_days_mean_CD,
         effect_light_P = slope_chosen_days_mean_PL - slope_chosen_days_mean_PD)

effects_difference_order %>%
  colnames()
# effects_difference_order_table <- effects_difference_order %>%
#   select(season, phylum, class, order, matches('mean'), matches('sd')) %>%
#   filter(slope_chosen_days_mean_CD != is.na(slope_chosen_days_mean_CD) &
#            slope_chosen_days_mean_CL != is.na(slope_chosen_days_mean_CL) &
#            slope_chosen_days_mean_PD != is.na(slope_chosen_days_mean_PD) &
#            slope_chosen_days_mean_PL != is.na(slope_chosen_days_mean_PL) &
#            slope_chosen_days_mean_DL != is.na(slope_chosen_days_mean_DL) &
#            slope_chosen_days_mean_VL != is.na(slope_chosen_days_mean_VL)) %>%
#   arrange(order)
# effects_difference_order_table %>%
#   colnames()
# 
# 
# effects_difference_order_table <- effects_difference_order_table %>%
#   rename( 'Mean growth rate CD' = slope_chosen_days_mean_CD,
#          'Mean growth rate CL' = slope_chosen_days_mean_CL,
#          'Mean growth rate PL' = slope_chosen_days_mean_PL,
#          'Mean growth rate PD'= slope_chosen_days_mean_PD,
#          'Mean growth rate DL' = slope_chosen_days_mean_DL,
#          'Mean growth rate VL'= slope_chosen_days_mean_VL,
#          'sd growth rate CD' = slope_chosen_sd_CD,
#          'sd growth rate CL' = slope_chosen_sd_CL,
#          'sd growth rate PL' = slope_chosen_sd_PL,
#          'sd growth rate PD' = slope_chosen_sd_PD,
#          'sd growth rate DL' = slope_chosen_sd_DL,
#          'sd growth rate VL' = slope_chosen_sd_VL)
# write.table(effects_difference_order_table, 'results/tables/mean_sd_gr_treatment_order.txt', sep = '\t')

effects_difference_order_l <- effects_difference_order %>%
  pivot_longer(cols = starts_with('effect'),
               names_to = 'effects',
               values_to = 'difference') %>%
  as_tibble()

effects_difference_order_l %>%
  colnames()

effects_difference_order %>%
  colnames()

effects_difference_order_l$season <-effects_difference_order_l$season %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))

effects_difference_order_l$effects <- effects_difference_order_l$effects %>% 
  factor(levels=(c("effect_bottom_up", "effect_top_down_grazers_dark", 
                   "effect_top_down_grazers_light", "effect_top_down_virus",
                   "effect_light_C", "effect_light_P")))

effects_difference_order_l <- effects_difference_order_l %>%
  mutate(phylum_f = as_factor(phylum),
         order_f = as_factor(order),
         class_f = as_factor(class))

effects_difference_order_l$class_f <-  factor(effects_difference_order_l$class_f, 
                                              levels=unique(effects_difference_order_l$class_f[order(effects_difference_order_l$phylum_f)]), 
                                              ordered=TRUE)

effects_difference_order_l$order_f <-  factor(effects_difference_order_l$order_f, 
                                              levels=unique(effects_difference_order_l$order_f[order(effects_difference_order_l$phylum_f,
                                                                                                     effects_difference_order_l$class_f)]), 
                                              ordered=TRUE)

difference_gr_treatments_order <- 
  effects_difference_order_l %>%
  filter(difference != is.na(difference)) %>%
  # group_by(class) %>%
  # filter(n() > 10) %>%
  group_by(order) %>%
  filter(n() > 15) %>%
  group_by(order, difference) %>%
  mutate(counts = n()) %>%
  ungroup() %>%
  # filter(difference > 0) %>%
  ggplot(aes(season, effects, color = difference))+ #, alpha = counts
  geom_point(aes(color = difference, size = difference))+
  scale_y_discrete(labels = effects_labels2)+
  scale_size(range = c(0, 12), name ="Growth rate difference\n between treatments")+
  scale_colour_gradientn(colours = palete_gradient)+
  # scale_color_manual(values = palette_phylums)+
  #scale_color_gradient2(low = '#B24D5D', high =  '#4DB2A2',  midpoint = 0)+
  #scale_alpha(range = c(1,1), limits = c(1, 1))+
  facet_wrap(vars(order_f))+ #, labeller = interaction(vars(class_f), vars(order_f), sep = ' ') #, labeller(order_f = ~ paste(order_f, class_f), .multi_line = TRUE)
  labs(y = 'Treatments growth rates differences', x = 'Season', color = '')+
  theme_bw()+
  theme(legend.position = 'right', panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 60, size = 5, hjust = 1), 
        axis.text.y = element_text(size = 5), legend.title = element_text(size = 7), axis.title = element_text(size = 7),
        strip.text = element_text(size = 7), strip.background = element_blank(), legend.text = element_text(size = 5))

ggsave('difference_gr_treatments_order5.pdf', difference_gr_treatments_order, 
       path = "~/Documentos/Doctorat/REMEI/results/figures/raw/",
       width = 280,
       height = 260,
       units = 'mm')

# relative_effects_family_l %>%
#   colnames()
# #relative_gr_treatments_family <- 
#   relative_effects_family_l %>%
#   filter(ratio != is.na(ratio)) %>%
#   group_by(class) %>%
#   filter(n() > 10) %>%
#   group_by(family) %>%
#   filter(n() > 10) %>%
#   # group_by(family, ratio) %>%
#   # mutate(counts = n()) %>%
#   # ungroup() %>%
#   filter(ratio > 1) %>%
#   ggplot(aes(season, effects, color = class_f, size = ratio))+#c(0, 40* sqrt(40)
#   geom_point()+
#   scale_y_discrete(labels = effects_labels2)+
#   #scale_radius(range = c(0, 40), name="Growth rate relative\n between treatments")+
#   scale_color_manual(values = palf_large_phylums(11))+
#   #scale_size_identity()+
#   #scale_alpha(range = c(0,40))+#, limits = c(0.1, 1)
#   facet_wrap(vars(family_f))+ #, labeller = interaction(vars(class_f), vars(family_f), sep = ' ')
#   labs(y = 'Treatments growth rates relatives', x = 'Season', color = 'Phylum')+
#   theme_bw()+
#   theme(legend.position = 'right', panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         axis.text.x = element_text(angle = 90), strip.text = element_text(size = 5), strip.background = element_blank())
# 
# ggsave('relative_gr_treatments_family2.pdf', relative_gr_treatments_family, 
#        path = "~/Documentos/Doctorat/REMEI/results/figures/raw/",
#        width = 240,
#        height = 180,
#        units = 'mm')

# effects_difference_family_l %>%
#   filter(difference != is.na(difference)) %>%
#   subset(phylum == 'Proteobacteria') %>%
#   group_by(class) %>%
#   filter(n() > 10) %>%
#   group_by(family) %>%
#   filter(n() > 21) %>%
#   ggplot(aes(season, effects, size = difference, color = class_f))+
#   geom_point()+
#   scale_y_discrete(labels = effects_labels2)+
#   scale_size(range = c(-3.5, 4.5), name=expression("Difference in \ngrowth rate day"^"-1"))+
#   #scale_color_manual(values = palette_large)+
#   facet_grid(~family_f)+
#   labs(y = 'Growth rate difference between treatments', x = 'Season', color = 'Phylum')+
#   theme_bw()+
#   theme(legend.position = 'right', panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         axis.text.x = element_text(angle = 90))
      
      # relative_effects_family_l %>%
      #   filter(ratio != is.na(ratio)) %>%
      #   subset(phylum == 'Proteobacteria') %>%
      #   group_by(class) %>%
      #   filter(n() > 10) %>%
      #   # group_by(family) %>%
      #   # filter(n() > 21) %>%
      #   ggplot(aes(season, effects, size = ratio, color = class_f))+
      #   geom_point()+
      #   scale_y_discrete(labels = effects_labels)+
      #   scale_size(range = c(-3.5, 4.5), name=expression("Difference in \ngrowth rate day"^"-1"))+
      #   scale_color_manual(values = palette_large)+
      #   facet_grid(~family_f)+
      #   labs(y = 'Growth rate difference between treatments', x = 'Season', color = 'Phylum')+
      #   theme_bw()+
      #   theme(legend.position = 'right', panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      #         axis.text.x = element_text(angle = 90))

#no funciona, volia afegir un histograma d'events al costat------
df2 <- effects_difference_family_l[rep(1:nrow(effects_difference_family_l), effects_difference_family_l$difference), c("x", "y") ]

xhist <- 
  axis_canvas(p, axis = "x") + 
  geom_histogram(data = df2, aes(x = x), color = 'lightgray')
yhist <-
  axis_canvas(p, axis = "y", coord_flip = TRUE) + 
  geom_histogram(data = df2, aes(x = y), color = 'lightgray') +
  coord_flip()


p %>%
  insert_xaxis_grob(xhist, grid::unit(1, "in"), position = "top") %>%
  insert_yaxis_grob(yhist, grid::unit(1, "in"), position = "right") %>%
  ggdraw()

library(ggExtra)
library(cowplot)
ggExtra::ggMarginal(p, type = 'densigram')
ggExtra::ggMarginal(p, data = effects_difference_family_l, type = 'histogram')

##testing different types of graphs----
effects_difference_family_l %>%
  filter(difference != is.na(difference)) %>%
  subset(phylum== 'Proteobacteria') %>%
  group_by(class_f) %>%
  filter(n() > 2) %>%
  group_by(order_f) %>%
  filter(n() > 10) %>%
  ggplot(aes(season, difference))+
  geom_violin(alpha = 0.5, draw_quantiles = c(0.25, 0.75))+
  geom_point(aes(color = effects), alpha = 0.8, position = position_jitter())+
  scale_color_manual(values = palette_treatments_remei)+
  facet_wrap(vars(order_f))+
  theme_bw()+
  theme(legend.position = 'right', panel.grid.major = element_blank(), panel.grid.minor = element_blank())

  
effects_difference_family_l %>%
  filter(difference != is.na(difference)) %>%
  subset(phylum == 'Proteobacteria') %>%
  group_by(class) %>%
  filter(n() > 10) %>%
  # group_by(family) %>%
  # filter(n() > 21) %>%
  ggplot(aes(season, effects))+
  geom_tile(aes(fill = difference))+
  scale_y_discrete(labels = effects_labels2)+
  #scale_size(range = c(-3.5, 4.5), name="Difference")+
  scale_fill_viridis_c()+
  facet_grid(~family_f)+
  labs(y = 'Growth rate difference between treatments', x = 'Season', color = 'Phylum')+
  theme_bw()+
  theme(legend.position = 'right', panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90))

effects_difference_family_l %>%
  filter(difference != is.na(difference)) %>%
  group_by(class) %>%
  filter(n() > 10) %>%
  group_by(family) %>%
  filter(n() > 20) %>%
  subset(phylum == 'Proteobacteria') %>%
  mutate(counts = n()) %>%
  ggplot(aes(difference, season, color = order_f))+
  geom_point(position = position_jitter_keep(0.25))+
  geom_boxplot(aes(group = season), alpha = 0.0)+
  #scale_y_discrete(labels = effects_labels2)+
  #scale_size(range = c(-3.5, 4.5), name="Growth rate difference\n between treatments")+
  scale_color_manual(values = palf_large(113))+
  geom_vline(xintercept =  0)+
  #facet_grid(effects ~class_f ~ family_f)+
  facet_wrap(vars(effects), nrow = 1, labeller = effects_labels2)+ #, labeller = interaction(vars(class_f), vars(family_f), sep = ' ')
  labs(y = 'Treatments growth rates differences', x = 'Season', color = 'Order')+
  theme_bw()+
  theme(legend.position = 'right', panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90), strip.text = element_text(size = 10), strip.background = element_blank())
  

##Just focusing on PREDATOR FREE TREATMENT-------
effects_difference_family_l %>%
  filter(difference != is.na(difference)) %>%
  group_by(class) %>%
  filter(n() > 10) %>%
  group_by(family) %>%
  filter(n() > 20) %>%
  subset(effects == 'effect_top_down_grazers_light') %>%
  mutate(counts = n()) %>%
  ggplot(aes(class, difference, color = season))+
  geom_point(position = position_jitter(0.25))+
  #geom_boxplot(aes(group = season), alpha = 0.0)+
  #scale_y_discrete(labels = effects_labels2)+
  #scale_size(range = c(-3.5, 4.5), name="Growth rate difference\n between treatments")+
  scale_color_manual(values = palette_seasons_4)+#values = palf_large(113)
  geom_vline(xintercept =  0)+
  #facet_grid(vars(season))+
  #facet_wrap(vars(class_f), nrow = 1)+ #, labeller = interaction(vars(class_f), vars(family_f), sep = ' ')
  labs(y = 'Treatments growth rates differences', x = 'Season', color = 'Order')+
  theme_bw()+
  theme(legend.position = 'right', panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90), strip.text = element_text(size = 10), strip.background = element_blank())

reg_all_slopes_chosen_silva_tax %>%
  #filter(difference != is.na(difference)) %>%
  # group_by(family) %>%
  # filter(n() > 20) %>%
  # filter(effects %in% c('effect_top_down_grazers_light',
  #                       'effect_top_down_grazers_dark')) %>%
#subset(treatment ==   c('CD',  'PD',  'PL',  'CL')
  subset(treatment %in%  c('PD', 'PL', 'CL', 'CD')) %>%
  filter(slope_chosen_days > 0) %>%
  group_by(order) %>%
  filter(n() > 10) %>%
  mutate(counts = n()) %>%
  ggplot(aes(fct_infreq(order), slope_chosen_days, color = phylum, shape = treatment))+
  geom_point(position = position_jitter(0.2), size = 2)+
  geom_violin(aes(group = order), alpha = 0.0)+
  #geom_boxplot(aes(group=order), alpha = 0)+
  #scale_y_discrete(labels = effects_labels2)+
  #scale_size(range = c(-3.5, 4.5), name="Growth rate difference\n between treatments")+
  scale_color_manual(values = palette_phylums)+#values = palf_large(113)
  geom_vline(xintercept =  0)+
  facet_grid(vars(treatment))+
  #facet_wrap(vars(class_f), nrow = 1)+ #, labeller = interaction(vars(class_f), vars(family_f), sep = ' ')
  labs(y = 'Treatments growth rates differences', x = 'Season', color = 'Order')+
  theme_bw()+
  theme(legend.position = 'right', panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90), strip.text = element_text(size = 10), strip.background = element_blank())

effects_difference_family_l %>%
  filter(difference != is.na(difference)) %>%
  # group_by(family) %>%
  # filter(n() > 20) %>%
  filter(effects %in% c('effect_top_down_grazers_light',
                      'effect_top_down_grazers_dark')) %>%
  filter(difference > 0) %>%
  group_by(order) %>%
  filter(n() > 4) %>%
  mutate(counts = n()) %>%
  ggplot(aes(fct_infreq(order), difference, color = phylum, shape = season))+
  geom_point(position = position_jitter(0.1), size = 3)+
  #geom_violin(aes(group = order), alpha = 0.0)+
  geom_boxplot(aes(group=order), alpha = 0)+
  #scale_y_discrete(labels = effects_labels2)+
  #scale_size(range = c(-3.5, 4.5), name="Growth rate difference\n between treatments")+
  scale_color_manual(values = palette_phylums)+#values = palf_large(113)
  geom_vline(xintercept =  0)+
  #facet_grid(vars(season))+
  #facet_wrap(vars(class_f), nrow = 1)+ #, labeller = interaction(vars(class_f), vars(family_f), sep = ' ')
  labs(y = 'Treatments growth rates differences', x = 'Season', color = 'Order')+
  theme_bw()+
  theme(legend.position = 'right', panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90), strip.text = element_text(size = 10), strip.background = element_blank())

##Plotting all treatments effects in general-----
summarize_effects_family_level_color_phylum <- effects_difference_family_l %>%
  filter(difference != is.na(difference)) %>%
  #subset(family == 'Pseudomonadaceae') %>%
  #group_by(class) %>%
  #filter(n() > 10) %>%
  #group_by(family) %>%
  #filter(n() > 20) %>%
  #subset(effects == 'effect_top_down_grazers_light') %>%
  #mutate(counts = n()) %>%
  ggplot(aes(effects, difference, color = phylum_f, group = effects))+ #, shape = effects
  geom_hline(yintercept =  0)+
  geom_point(position = position_jitter(0.15),  alpha = 0.8)+
  #geom_boxplot(aes(group = season), alpha = 0.0)+
  geom_violin(aes(group = effects), alpha = 0.0, draw_quantiles = c(0.25, 0.75))+
  stat_summary(fun = "mean", geom = "crossbar", 
               width = 0.6, colour = "black")+
  scale_x_discrete(labels = effects_labels2)+
  scale_y_continuous(limits = c(-5,5))+
  #scale_size(range = c(-3.5, 4.5), name="Growth rate difference\n between treatments")+
  scale_color_manual(values = palette_phylums_assigned)+#values = palf_large(113)
  #facet_wrap(vars(effects), labeller = effects_labels, nrow = 3)+
  #facet_wrap(vars(class_f), nrow = 1)+ #, labeller = interaction(vars(class_f), vars(family_f), sep = ' ')
  guides(color=guide_legend(ncol = 3))+
  labs(y = expression("Difference in growth rate day"^"-1"), x = 'Season', color = 'Phylum')+
  theme_bw()+
  theme(legend.position = 'right', panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1, size = 5), strip.text = element_text(size = 5), 
        axis.title = element_text(size = 7), axis.text.y = element_text(size = 5), legend.text = element_text(size = 5),
        legend.title = element_text(size = 7),
        strip.background = element_blank())

ggsave('summarize_effects_family_level_color_phylum.pdf', summarize_effects_family_level_color_phylum, 
       path = "~/Documentos/Doctorat/REMEI/results/figures/raw/",
       width = 180,
       height = 100,
       units = 'mm')

##scatter plot general (testing)-----
effects_difference_family %>%
  colnames()

effects_difference_family %>%
  ggplot(aes(slope_chosen_days_mean_CD, slope_chosen_days_mean_PD, color = season))+
  geom_point(position = position_jitter(0.25),  alpha = 0.7)+
  stat_ellipse(aes(group = season, color = season), 
               linetype = 2,
               lwd = 1.2)+
  scale_color_manual(values = palette_seasons_4)+#values = palf_large(113)
  theme_bw()+
  theme(legend.position = 'right', panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90), strip.text = element_text(size = 10), strip.background = element_blank())

effects_difference_family %>%
  # group_by(season, effects, family) %>%
  # summarize(count = n()) %>%
  # filter(count > 4) no puc filtrar perquè está la taula en format wide
  ggplot(aes(slope_chosen_days_mean_CD, slope_chosen_days_mean_PD, color = family))+
  geom_point(position = position_jitter(0.25),  alpha = 0.7)+
  stat_ellipse(aes(group = class, color = family), 
               linetype = 2)+ #               lwd = 1.2
  #scale_color_manual(values = palf_large(261))+#values = palf_large(113)
  facet_wrap(vars(season))+
  theme_bw()+
  theme(legend.position = 'none', panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90), strip.text = element_text(size = 10), strip.background = element_blank())


effects_difference_family %>%
  ggplot(aes(slope_chosen_days_mean_CL, slope_chosen_days_mean_PL, color = season))+
  geom_point(position = position_jitter(0.25),  alpha = 0.7)+
  scale_color_manual(values = palette_seasons_4)+#values = palf_large(113)
  theme_bw()+
  theme(legend.position = 'right', panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90), strip.text = element_text(size = 10), strip.background = element_blank())

effects_difference_family %>%
  ggplot(aes(slope_chosen_days_mean_CL, slope_chosen_days_mean_CD, color = season))+
  geom_point(position = position_jitter(0.25),  alpha = 0.7)+
  scale_color_manual(values = palette_seasons_4)+#values = palf_large(113)
  theme_bw()+
  theme(legend.position = 'right', panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90), strip.text = element_text(size = 10), strip.background = element_blank()) 

effects_difference_family %>%
  ggplot(aes(slope_chosen_days_mean_PL, slope_chosen_days_mean_PD, color = season))+
  geom_point(position = position_jitter(0.25),  alpha = 0.7)+
  scale_color_manual(values = palette_seasons_4)+#values = palf_large(113)
  theme_bw()+
  theme(legend.position = 'right', panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90), strip.text = element_text(size = 10), strip.background = element_blank()) 
  
effects_difference_family %>%
  ggplot(aes(slope_chosen_days_mean_VL, slope_chosen_days_mean_DL, color = season))+
  geom_point(position = position_jitter(0.25),  alpha = 0.7)+
  scale_color_manual(values = palette_seasons_4)+#values = palf_large(113)
  theme_bw()+
  theme(legend.position = 'right', panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90), strip.text = element_text(size = 10), strip.background = element_blank())  

effects_difference_family %>%
  subset(phylum == 'Proteobacteria')%>%
  ggplot(aes(slope_chosen_days_mean_DL, slope_chosen_days_mean_PL, color = season))+
  geom_point(position = position_jitter(0.25),  alpha = 0.7)+
  scale_color_manual(values = palette_seasons_4)+#values = palf_large(113)
  #geom_polygon(group = class, fill = class, stat = 'identity', position = 'dodge', width ='0.5', size = 0.5)+
  facet_wrap(vars(class))+
  #geom_smooth(method = 'loess', aes(group = order))+
  theme_bw()+
  theme(legend.position = 'right', panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90), strip.text = element_text(size = 10), strip.background = element_blank()) 


effects_difference_family %>%
  ggplot(aes(order, effect_top_down_grazers_dark, color = season))+
  geom_point(position = position_jitter(0.25),  alpha = 0.7)+
  scale_color_manual(values = palette_seasons_4)+#values = palf_large(113)
  theme_bw()+
  theme(legend.position = 'right', panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90), strip.text = element_text(size = 10), strip.background = element_blank()) 

##a nivell d'asv què passa?
reg_all_slopes_chosen_silva_tax %>%
  colnames() #
reg_all_slopes_chosen_silva_tax_filt %>%
 test %>%
  colnames()

reg_all_slopes_chosen_silva_tax_filt %>%
  distinct(treatment, season, asv_num, domain, phylum, class, order, family, genus, species, tax_ed, .keep_all = TRUE) %>%
  pivot_wider(names_from = treatment, values_from = c(slope_chosen_days, slope_chosen, pvalue_slope_chosen), 
              id_cols = c(asv_num, season, domain, phylum, class, order, family, genus, species, tax_ed)) %>%
  group_by(class, season) %>%
  filter(n() > 20) %>%
  #subset(phylum == 'Proteobacteria')%>%
  ggplot(aes(slope_chosen_days_DL, slope_chosen_days_PL, color = season))+
  geom_point(position = position_jitter(0.25),  alpha = 0.7)+
  stat_poly_line(color = 'black')+
  stat_poly_eq(color = 'black')+
  scale_color_manual(values = palette_seasons_4)+#values = palf_large(113)
  #geom_polygon(group = class, fill = class, stat = 'identity', position = 'dodge', width ='0.5', size = 0.5)+
  facet_wrap(vars(class))+
  stat_ellipse(aes(group = season),  #, color = family
               linetype = 2)+
  #geom_smooth(method = 'loess', aes(group = order))+
  theme_bw()+
  theme(legend.position = 'right', panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90), strip.text = element_text(size = 10), strip.background = element_blank()) 

##probo slope graphs
effects_difference_family_l %>%
  filter(difference != is.na(difference)) %>%
  #subset(class %in% c('Alphaproteobacteria', 'Gammaproteobacteria', 'Bacteroidia')) %>%
  group_by(class_f, season, effects) %>%
  #filter(n() > 25) %>%
  ggplot(aes(effects, difference, group = family_f))+
  geom_line(aes(color = family_f))+
  geom_point(aes(color = family_f))+
  #scale_color_manual(values = palf_large_phylums(116))+
  scale_x_discrete(labels = effects_labels)+
  #scale_y_continuous(limits = c(0,5))+
  facet_wrap(vars(class_f))+
  theme_bw()+
  theme(legend.position = 'none')

##posicionament de les families entre control top-down o bottom-up
relative_effects_family %>%
  colnames()

relative_effects_family %>%
  ggplot(aes(effect_top_down_grazers_light, effect_bottom_up, color = phylum))+
  geom_point()+
  scale_color_manual(values = palette_phylums)+
  scale_x_continuous(limits = c(0, 10))+
  scale_y_continuous(limits = c(0, 10))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  theme_bw()

  effects_difference_order %>%
  ggplot(aes(effect_top_down_grazers_light, effect_bottom_up, color = class))+
  geom_point(size = 4)+ #mida hauria de ser el número de OTUs dins de cada família
 # scale_color_manual(values = palf_large_phylums(60))+
  scale_x_continuous(limits = c(0, 3))+
  scale_y_continuous(limits = c(0, 3))+
  #facet_wrap(vars(season))+
  theme_bw()

#diferents combinacions i després les uneixo color per phylums
effects_labels2 <-  as_labeller(c(effect_top_down_virus = 'Top-down effect viruses (VL-DL)',
                                  effect_top_down_grazers_dark = 'Top-down effect grazers dark (PD-CD)',
                                  effect_top_down_grazers_light = 'Top-down effect grazers light (PL-CL)',
                                  effect_bottom_up = 'Bottom-up effect (DL-PL)',
                                  effect_light_C = 'Light effect (CL-CD)',
                                  effect_light_P = 'Light effect (PL-PD)'))

effects_difference_family$season <- effects_difference_family$season %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))


effects_difference_family_l$season <-effects_difference_family_l$season %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))

effects_difference_family_l$effects <- effects_difference_family_l$effects %>% 
  factor(levels=(c("effect_bottom_up", "effect_top_down_grazers_dark", 
                   "effect_top_down_grazers_light", "effect_top_down_virus",
                   "effect_light_C", "effect_light_P")))

effects_difference_family <- effects_difference_family %>%
  mutate(phylum_f = as_factor(phylum),
         family_f = as_factor(family),
         order_f = as_factor(order),
         class_f = as_factor(class))

effects_difference_family$class_f <-  factor(effects_difference_family$class_f, 
                                               levels=unique(effects_difference_family$class_f[order(effects_difference_family$phylum_f)]), 
                                               ordered=TRUE)

effects_difference_family$order_f <-  factor(effects_difference_family$order_f, 
                                               levels=unique(effects_difference_family$order_f[order(effects_difference_family$phylum_f,
                                                                                                       effects_difference_family$class_f)]), 
                                               ordered=TRUE)

effects_difference_family$family_f <-  factor(effects_difference_family$family_f, 
                                                levels=unique(effects_difference_family$family_f[order(effects_difference_family$phylum_f,
                                                                                                         effects_difference_family$class_f,
                                                                                                         effects_difference_family$order_f)]), 
                                                ordered=TRUE)

p1<- effects_difference_family %>%
  ggplot(aes(effect_top_down_grazers_light, effect_bottom_up, color = phylum))+
  scale_x_continuous(limits = c(-3, 3))+
  scale_y_continuous(limits = c(-3, 3))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  labs(x = 'Top-down effect grazers light (PL-CL)', y = 'Bottom-up effect (DL-PL)')+
  geom_point(size = 3, alpha =  0.8)+ #mida hauria de ser el número de OTUs dins de cada família
  scale_color_manual(values = palette_phylums_assigned)+
  facet_wrap(vars(season), nrow = 1)+
  theme_bw()+
  theme(legend.position = 'none', panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), strip.background = element_blank())

p2 <- effects_difference_family %>%
  ggplot(aes(effect_top_down_virus, effect_bottom_up, color = phylum))+
  scale_x_continuous(limits = c(-3, 3))+
  scale_y_continuous(limits = c(-3, 3))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  scale_color_manual(values = palette_phylums_assigned)+
  labs(x = 'Top-down effect viruses (VL-DL)', y = 'Bottom-up effect (DL-PL)')+
  geom_point(size = 3, alpha =  0.8)+ #mida hauria de ser el número de OTUs dins de cada família
  facet_wrap(vars(season), nrow = 1)+
  theme_bw()+
  theme(legend.position = 'none', panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_blank(),
        strip.text = element_blank())

p3 <- effects_difference_family %>%
  ggplot(aes(effect_top_down_virus, effect_top_down_grazers_light, color = phylum))+
  scale_x_continuous(limits = c(-3, 3))+
  scale_y_continuous(limits = c(-3, 3))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  geom_point(size = 3, alpha =  0.8)+ #mida hauria de ser el número de OTUs dins de cada família
  scale_color_manual(values = palette_phylums_assigned)+
  labs(x = 'Top-down effect viruses (VL-DL)', y = 'Top-down effect grazers light (PL-CL)')+
  facet_wrap(vars(season), nrow = 1)+
  theme_bw()+
  theme(legend.position = 'none', panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_blank(),
        strip.text = element_blank())

p4 <- effects_difference_family %>%
  ggplot(aes(effect_light_C, effect_light_P, color = phylum))+
  scale_x_continuous(limits = c(-3, 3))+
  scale_y_continuous(limits = c(-3, 3))+
  labs(x = 'Light effect (CL-CD)', y = 'Light effect (PL-PD)')+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  geom_point(size = 3, alpha =  0.8)+ #mida hauria de ser el número de OTUs dins de cada família
  scale_color_manual(values = palette_phylums_assigned)+
  facet_wrap(vars(season), nrow = 1)+
  theme_bw()+
  theme(legend.position = 'none', panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_blank(),
        strip.text = element_blank())

grid.arrange(p1, p2, p3, p4, ncol = 1)

##el mateix però color per classes----
effects_difference_family$season <- effects_difference_family$season %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))

p1 <- effects_difference_family %>%
  group_by(season, phylum_f) %>%
  filter(n() > 2) %>%
  group_by(season, class_f) %>%
  filter(n() > 2) %>%
  ggplot(aes(effect_top_down_grazers_light, effect_bottom_up, color = class_f))+
  scale_x_continuous(limits = c(-3, 3))+
  scale_y_continuous(limits = c(-3, 3))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  labs(x = 'Top-down effect grazers light (PL-CL)', y = 'Bottom-up effect (DL-PL)')+
  geom_point(size = 2, alpha =  0.8)+ #mida hauria de ser el número de OTUs dins de cada família
  scale_color_manual(values = palette_class_assigned)+
  facet_wrap(vars(season), nrow = 1)+
  theme_bw()+
  theme(legend.position = 'none',panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), strip.background = element_blank(), legend.text = element_text(size = 5),
        legend.title = element_text(size = 7), axis.title = element_text(size = 7), axis.text = element_text(size = 5))

p2 <- effects_difference_family %>%
  group_by(season, phylum_f) %>%
  filter(n() > 2) %>%
  group_by(season, class_f) %>%
  filter(n() > 2) %>%
  ggplot(aes(effect_top_down_virus, effect_bottom_up, color = class_f))+
  scale_x_continuous(limits = c(-3, 3))+
  scale_y_continuous(limits = c(-3, 3))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  scale_color_manual(values = palette_class_assigned)+
  labs(x = 'Top-down effect viruses (VL-DL)', y = 'Bottom-up effect (DL-PL)')+
  geom_point(size = 2, alpha =  0.8)+ #mida hauria de ser el número de OTUs dins de cada família
  facet_wrap(vars(season), nrow = 1)+
  theme_bw()+
  theme(legend.position = 'none', panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), strip.background = element_blank(), legend.text = element_text(size = 5),
        legend.title = element_text(size = 7), axis.title = element_text(size = 7), axis.text = element_text(size = 5))

p3 <- effects_difference_family %>%
  group_by(season, phylum_f) %>%
  filter(n() > 2) %>%
  group_by(season, class_f) %>%
  filter(n() > 2) %>%
  ggplot(aes(effect_top_down_virus, effect_top_down_grazers_light, color = class_f))+
  scale_x_continuous(limits = c(-3, 3))+
  scale_y_continuous(limits = c(-3, 3))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  geom_point(size = 2, alpha =  0.8)+ #mida hauria de ser el número de OTUs dins de cada família
  scale_color_manual(values = palette_class_assigned)+
  labs(x = 'Top-down effect viruses (VL-DL)', y = 'Top-down effect grazers light (PL-CL)')+
  facet_wrap(vars(season), nrow = 1)+
  theme_bw()+
  theme(legend.position = 'none',panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), strip.background = element_blank(), legend.text = element_text(size = 5),
        legend.title = element_text(size = 7), axis.title = element_text(size = 7), axis.text = element_text(size = 5))


p4 <- effects_difference_family %>%
  group_by(season, phylum_f) %>%
  filter(n() > 2) %>%
  group_by(season, class_f) %>%
  filter(n() > 2) %>%
  ggplot(aes(effect_light_C, effect_light_P, color = class_f))+
  scale_x_continuous(limits = c(-3, 3))+
  scale_y_continuous(limits = c(-3, 3))+
  labs(x = 'Light effect (CL-CD)', y = 'Light effect (PL-PD)', color = 'Taxonomic\n class')+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  geom_point(size = 2, alpha =  0.8)+ #mida hauria de ser el número de OTUs dins de cada família
  scale_color_manual(values = palette_class_assigned)+
  facet_wrap(vars(season), nrow = 1)+
  guides(color=guide_legend(ncol = 7, keyheight = 1))+
  theme_bw()+
  theme(legend.position = 'bottom', panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), strip.background = element_blank(), legend.text = element_text(size = 5),
      legend.title = element_text(size = 7), axis.title = element_text(size = 7), axis.text = element_text(size = 5))

grid.arrange(p1, p2, p3, p4, ncol = 1)

correlations_differences_family_treatments_class  <- multi_panel_figure(columns = 1, rows = 5, width = 200, height = 250, 
                                                                  row_spacing = 2, unit = 'mm',
                                                                  panel_label_type = 'none')

correlations_differences_family_treatments_class  %<>%
  fill_panel(p1, row = 1, col = 1) %<>%
  fill_panel(p2, row = 2, col = 1) %<>%
  fill_panel(p3, row = 3, col = 1) %<>%
  fill_panel(p4, row = 4:5, col = 1)

ggsave('correlations_differences_family_treatments_class_v2.pdf', correlations_differences_family_treatments_class, 
       path = "~/Documentos/Doctorat/REMEI/results/figures/raw/",
       width = 200,
       height = 250,
       units = 'mm')


    ##order level----
#diferents combinacions i després les uneixo color per phylums
effects_labels2 <-  as_labeller(c(effect_top_down_virus = 'Top-down effect viruses (VL-DL)',
                                  effect_top_down_grazers_dark = 'Top-down effect grazers dark (PD-CD)',
                                  effect_top_down_grazers_light = 'Top-down effect grazers light (PL-CL)',
                                  effect_bottom_up = 'Bottom-up effect (DL-PL)',
                                  effect_light_C = 'Light effect (CL-CD)',
                                  effect_light_P = 'Light effect (PL-PD)'))

effects_difference_order$season <- effects_difference_order$season %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))


effects_difference_order_l$season <-effects_difference_order_l$season %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))

effects_difference_order_l$effects <- effects_difference_order_l$effects %>% 
  factor(levels=(c("effect_bottom_up", "effect_top_down_grazers_dark", 
                   "effect_top_down_grazers_light", "effect_top_down_virus",
                   "effect_light_C", "effect_light_P")))

effects_difference_order <- effects_difference_order %>%
  mutate(phylum_f = as_factor(phylum),
         #family_f = as_factor(family),
         order_f = as_factor(order),
         class_f = as_factor(class))

effects_difference_order$class_f <-  factor(effects_difference_order$class_f, 
                                             levels=unique(effects_difference_order$class_f[order(effects_difference_order$phylum_f)]), 
                                             ordered=TRUE)

effects_difference_order$order_f <-  factor(effects_difference_order$order_f, 
                                             levels=unique(effects_difference_order$order_f[order(effects_difference_order$phylum_f,
                                                                                                   effects_difference_order$class_f)]), 
                                             ordered=TRUE)


p1 <- effects_difference_order %>%
  group_by(season, phylum_f) %>%
  filter(n() > 2) %>%
  group_by(season, class_f) %>%
  filter(n() > 2) %>%
  ggplot(aes(effect_top_down_grazers_light, effect_bottom_up, color = class))+
  scale_x_continuous(limits = c(-3, 3))+
  scale_y_continuous(limits = c(-3, 3))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  labs(x = 'Top-down effect grazers light (PL-CL)', y = 'Bottom-up effect (DL-PL)')+
  geom_point(size = 2, alpha =  0.8)+ #mida hauria de ser el número de OTUs dins de cada família
  scale_color_manual(values = palette_class_assigned)+
  facet_wrap(vars(season), nrow = 1)+
  theme_bw()+
  theme(legend.position = 'none',panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), strip.background = element_blank(), legend.text = element_text(size = 5),
        legend.title = element_text(size = 7), axis.title = element_text(size = 7), axis.text = element_text(size = 5))

p2 <- effects_difference_order %>%
  group_by(season, phylum_f) %>%
  filter(n() > 2) %>%
  group_by(season, class_f) %>%
  filter(n() > 2) %>%
  ggplot(aes(effect_top_down_virus, effect_bottom_up, color = class))+
  scale_x_continuous(limits = c(-3, 3))+
  scale_y_continuous(limits = c(-3, 3))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  scale_color_manual(values = palette_class_assigned)+
  labs(x = 'Top-down effect viruses (VL-DL)', y = 'Bottom-up effect (DL-PL)')+
  geom_point(size = 2, alpha =  0.8)+ #mida hauria de ser el número de OTUs dins de cada família
  facet_wrap(vars(season), nrow = 1)+
  theme_bw()+
  theme(legend.position = 'none',panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), strip.background = element_blank(), legend.text = element_text(size = 5),
        legend.title = element_text(size = 7), axis.title = element_text(size = 7), axis.text = element_text(size = 5))

p3 <- effects_difference_order %>%
  group_by(season, phylum_f) %>%
  filter(n() > 2) %>%
  group_by(season, class_f) %>%
  filter(n() > 2) %>%
  ggplot(aes(effect_top_down_virus, effect_top_down_grazers_light, color = class))+
  scale_x_continuous(limits = c(-3, 3))+
  scale_y_continuous(limits = c(-3, 3))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  geom_point(size = 2, alpha =  0.8)+ #mida hauria de ser el número de OTUs dins de cada família
  scale_color_manual(values = palette_class_assigned)+
  labs(x = 'Top-down effect viruses (VL-DL)', y = 'Top-down effect grazers light (PL-CL)')+
  facet_wrap(vars(season), nrow = 1)+
  theme_bw()+
  theme(legend.position = 'none',panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), strip.background = element_blank(), legend.text = element_text(size = 5),
        legend.title = element_text(size = 7), axis.title = element_text(size = 7), axis.text = element_text(size = 5))

p4 <- effects_difference_order %>%
  group_by(season, phylum_f) %>%
  filter(n() > 2) %>%
  group_by(season, class_f) %>%
  filter(n() > 2) %>%
  ggplot(aes(effect_light_C, effect_light_P, color =  class_f))+ #interaction(class_f, family_f,  sep = ' ')
  scale_x_continuous(limits = c(-3, 3))+
  scale_y_continuous(limits = c(-3, 3))+
  labs(x = 'Light effect (CL-CD)', y = 'Light effect (PL-PD)', color = 'Class')+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  guides(color=guide_legend(ncol = 3))+
  geom_point(size = 2, alpha =  0.8)+ #mida hauria de ser el número de OTUs dins de cada família
  scale_color_manual(values = palette_class_assigned)+
  facet_wrap(vars(season), nrow = 1)+
  guides(color=guide_legend(ncol = 5, keyheight = 0.9))+
  theme_bw()+
  theme(legend.position = 'bottom', panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), strip.background = element_blank(), legend.text = element_text(size = 5),
        legend.title = element_text(size = 7), axis.title = element_text(size = 7), axis.text = element_text(size = 5))

#grid.arrange(p1, p2, p3, p4, ncol = 1)

correlations_differences_order_treatments  <- multi_panel_figure(columns = 1, rows = 6, width = 180, height = 220, 
                                                                  row_spacing = 2, unit = 'mm',
                                                                  panel_label_type = 'none')

correlations_differences_order_treatments  %<>%
  fill_panel(p1, row = 1, col = 1) %<>%
  fill_panel(p2, row = 2, col = 1) %<>%
  fill_panel(p3, row = 3, col = 1) %<>%
  fill_panel(p4, row = 4:6, col = 1)

ggsave('correlations_differences_order_treatments_v2.pdf', correlations_differences_order_treatments, 
       path = "~/Documentos/Doctorat/REMEI/results/figures/raw/",
       width = 180,
       height = 220,
       units = 'mm')

##correlogram----
library(corrr)
effects_difference_family_l %>%
  colnames()

effects_difference_family %>%
  select(effect_top_down_virus, effect_top_down_grazers_dark, effect_top_down_grazers_light,
         effect_bottom_up, effect_light_C, effect_light_P) %>%
  correlate() %>%
  rplot()
  
##heatmap
effects_difference_family_l %>%
  filter(difference != is.na(difference)) %>%
  group_by(season, family_f) %>%
  filter(n() > 4) %>%
  ggplot(aes(effects, family_f, fill = difference))+
  geom_tile(color = 'white', lwd = 1.5, linetype = 0.5)+
  #geom_tile(aes(fill = ratio), color = 'white', lwd = 1.5, linetype = 0.5)+
  scale_fill_gradient2(low = '#BD4247', high =  '#42BDB8', na.value = '#fdfdfd', midpoint = 0)+
  #geom_point()+ #mida hauria de ser el número de OTUs dins de cada família
  #scale_color_manual(values = palf_large(100))+
  #scale_x_continuous(limits = c(0, 3))+
  #scale_y_continuous(limits = c(0, 3))+
  scale_x_discrete(labels = effects_labels)+
  #facet_wrap(vars(season), nrow = 1)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

##pels que están afectats per la llum
effects_difference_family_l %$%
  effects
effects_difference_family_l %>%
  filter(difference != is.na(difference)) %>%
  subset(effects == c('effect_light_P', 
         'effect_light_C')) %>%
  group_by(season, family_f) %>%
  filter(n() > 1) %>%
  filter(difference > 0) %>%
  ggplot(aes(effects, interaction(class_f, family_f,  sep = ' '), fill = difference))+
  geom_tile(color = 'white', lwd = 1.5, linetype = 0.5)+
  #geom_tile(aes(fill = ratio), color = 'white', lwd = 1.5, linetype = 0.5)+
  scale_fill_gradient2(low = '#BD4247', high =  '#42BDB8', na.value = '#fdfdfd', midpoint = 0)+
  #geom_point()+ #mida hauria de ser el número de OTUs dins de cada família
  #scale_color_manual(values = palf_large(100))+
  #scale_x_continuous(limits = c(0, 3))+
  #scale_y_continuous(limits = c(0, 3))+
  scale_x_discrete(labels = effects_labels)+
  facet_wrap(vars(season), nrow = 1)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


##filtrem només efectes positius
negatiu <- effects_difference_family_l %>%
  filter(difference != is.na(difference)) %>%
  filter(difference > 0) %>%
  group_by(season, family_f) %>%
  filter(n() > 3) %>%
  ggplot(aes(effects, interaction(family_f, class_f, sep = ' '), fill = difference))+
  geom_tile(color = 'white', lwd = 1.5, linetype = 0.5)+
  #geom_tile(aes(fill = ratio), color = 'white', lwd = 1.5, linetype = 0.5)+
  scale_fill_gradient2(low = '#BD4247', high =  '#42BDB8', na.value = '#fdfdfd', midpoint = 0)+
  #geom_point()+ #mida hauria de ser el número de OTUs dins de cada família
  #scale_color_manual(values = palf_large(100))+
  #scale_x_continuous(limits = c(0, 3))+
  #scale_y_continuous(limits = c(0, 3))+
  scale_x_discrete(labels = effects_labels)+
  facet_wrap(vars(season), nrow = 1)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

##filtrem només efectes negatius
positiu <- effects_difference_family_l %>%
  filter(difference != is.na(difference)) %>%
  filter(difference < 0) %>%
  group_by(season, family_f) %>%
  filter(n() > 3) %>%
  ggplot(aes(effects, interaction(family_f, class_f, sep = ' '), fill = difference))+
  geom_tile(color = 'white', lwd = 1.5, linetype = 0.5)+
  #geom_tile(aes(fill = ratio), color = 'white', lwd = 1.5, linetype = 0.5)+
  scale_fill_gradient2(low = '#BD4247', high =  '#42BDB8', na.value = '#fdfdfd', midpoint = 0)+
  #geom_point()+ #mida hauria de ser el número de OTUs dins de cada família
  #scale_color_manual(values = palf_large(100))+
  #scale_x_continuous(limits = c(0, 3))+
  #scale_y_continuous(limits = c(0, 3))+
  scale_x_discrete(labels = effects_labels)+
  facet_wrap(vars(season), nrow = 1)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

grid.arrange(positiu, negatiu, ncol = 1)

#REMIAU DATASET (WITH IN SITU DATA) (DISCARDED!) -------------------------
#OPCIÓ 2: PROVEM DE FER RESTA EN COMPTES DE DIVISIÓ ------
effects_difference_family <- reg_all_slopes_chosen_silva_tax %>% 
  distinct(treatment, season, asv_num, domain, phylum, class, order, genus, species, tax_ed, asv_num, family, .keep_all = TRUE) %>%
  filter(pvalue_slope_chosen < 0.05) %>%
  group_by(treatment, season, domain, phylum, class, order, family) %>%
  summarise(slope_chosen_days_mean = mean(slope_chosen_days),
            slope_chosen_mean = mean(slope_chosen),
            slope_chosen_days_sd = sd(slope_chosen_days),
            # counts = add_count(),
            na.rm = TRUE) %>%
  pivot_wider(names_from = treatment, values_from = c(slope_chosen_days_mean, slope_chosen_mean, slope_chosen_days_sd), 
              id_cols = c(season, domain, phylum, class, order, family)) %>%
  as_tibble() %>%
  mutate(across(!c(domain, phylum, class, order, family, season), as.numeric)) %>%
  mutate(effect_top_down_virus = slope_chosen_days_mean_VL - slope_chosen_days_mean_DL,
         effect_top_down_grazers_dark = slope_chosen_days_mean_PD - slope_chosen_days_mean_CD,
         effect_top_down_grazers_light = slope_chosen_days_mean_PL - slope_chosen_days_mean_CL,
         effect_bottom_up = slope_chosen_days_mean_DL - slope_chosen_days_mean_PL,
         effect_light_C = slope_chosen_days_mean_CL - slope_chosen_days_mean_CD,
         effect_light_P = slope_chosen_days_mean_PL - slope_chosen_days_mean_PD)

effects_difference_family %>%
  colnames()

effects_difference_family_l <- effects_difference_family %>%
  pivot_longer(cols = starts_with('effect'),
               names_to = 'effects',
               values_to = 'difference') %>%
  as_tibble()

effects_difference_family_l %>%
  colnames()


effects_difference_family_l$season <-effects_difference_family_l$season %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))

effects_difference_family_l$effects <- effects_difference_family_l$effects %>% 
  factor(levels=(c("effect_bottom_up", "effect_top_down_grazers_dark", 
                   "effect_top_down_grazers_light", "effect_top_down_virus",
                   "effect_light_C", "effect_light_P")))
effects_difference_family_l <- effects_difference_family_l %>%
  mutate(phylum_f = as_factor(phylum),
         family_f = as_factor(family),
         order_f = as_factor(order),
         class_f = as_factor(class))

effects_difference_family_l$class_f <-  factor(effects_difference_family_l$class_f, 
                                               levels=unique(effects_difference_family_l$class_f[order(effects_difference_family_l$phylum_f)]), 
                                               ordered=TRUE)

effects_difference_family_l$order_f <-  factor(effects_difference_family_l$order_f, 
                                               levels=unique(effects_difference_family_l$order_f[order(effects_difference_family_l$phylum_f,
                                                                                                       effects_difference_family_l$class_f)]), 
                                               ordered=TRUE)

effects_difference_family_l$family_f <-  factor(effects_difference_family_l$family_f, 
                                                levels=unique(effects_difference_family_l$family_f[order(effects_difference_family_l$phylum_f,
                                                                                                         effects_difference_family_l$class_f,
                                                                                                         effects_difference_family_l$order_f)]), 
                                                ordered=TRUE)
#interaction(class_f,phylum_f, sep = ' ')

difference_gr_treatments_family <- effects_difference_family_l %>%
  filter(difference != is.na(difference)) %>%
  group_by(class) %>%
  filter(n() > 10) %>%
  group_by(family) %>%
  filter(n() > 20) %>%
  group_by(family, difference) %>%
  mutate(counts = n()) %>%
  ungroup() %>%
  filter(difference > 0) %>%
  ggplot(aes(season, effects, size = difference, color = phylum_f, alpha = counts))+
  geom_point()+
  scale_y_discrete(labels = effects_labels2)+
  scale_size(range = c(0, 6), name="Growth rate difference\n between treatments")+
  scale_color_manual(values = palette_phylums)+
  scale_alpha(range = c(1,1), limits = c(0.1, 1))+
  facet_wrap(vars(family_f))+ #, labeller = interaction(vars(class_f), vars(family_f), sep = ' ')
  labs(y = 'Treatments growth rates differences', x = 'Season', color = 'Phylum')+
  theme_bw()+
  theme(legend.position = 'right', panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90), strip.text = element_text(size = 5), strip.background = element_blank())

# ggsave('remiau_difference_gr_treatments_family2.pdf', difference_gr_treatments_family, 
#        path = "~/Documentos/Doctorat/REMEI/results/figures/raw/",
#        width = 240,
#        height = 180,
#        units = 'mm')

##color and size by difference in GR
palete_gradient <- c("#5e0000",
                     "#b24d5d",
                     '#FFFFFF',
                     "#4db2a2",
                     "#005a47") 

difference_gr_treatments_family <- 
  effects_difference_family_l %>%
  filter(difference != is.na(difference)) %>%
  group_by(class) %>%
  filter(n() > 10) %>%
  group_by(family) %>%
  filter(n() > 15) %>%
  group_by(family, difference) %>%
  mutate(counts = n()) %>%
  # ungroup() %>%
  # filter(difference > 0) %>%
  ggplot(aes(season, effects, color = difference))+ #, alpha = counts
  geom_point(aes(color = difference, size = difference))+
  scale_y_discrete(labels = effects_labels2)+
  scale_size(range = c(0, 10), name ="Growth rate difference\n between treatments")+
  scale_colour_gradientn(colours = palete_gradient)+
  # scale_color_manual(values = palette_phylums)+
  #scale_color_gradient2(low = '#B24D5D', high =  '#4DB2A2',  midpoint = 0)+
  #scale_alpha(range = c(1,1), limits = c(1, 1))+
  facet_wrap(vars(family_f))+ #, labeller = interaction(vars(class_f), vars(family_f), sep = ' ')
  labs(y = 'Treatments growth rates differences', x = 'Season', color = '')+
  theme_bw()+
  theme(legend.position = 'right', panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90), strip.text = element_text(size = 11), strip.background = element_blank())

ggsave('remiau_difference_gr_treatments_family3.pdf', difference_gr_treatments_family, 
       path = "~/Documentos/Doctorat/REMEI/results/figures/raw/",
       width = 360,
       height = 290,
       units = 'mm')

relative_effects_family_l %>%
  colnames()
#relative_gr_treatments_family <- 
relative_effects_family_l %>%
  filter(ratio != is.na(ratio)) %>%
  group_by(class) %>%
  filter(n() > 10) %>%
  group_by(family) %>%
  filter(n() > 10) %>%
  # group_by(family, ratio) %>%
  # mutate(counts = n()) %>%
  # ungroup() %>%
  filter(ratio > 1) %>%
  ggplot(aes(season, effects, color = class_f, size = ratio))+#c(0, 40* sqrt(40)
  geom_point()+
  scale_y_discrete(labels = effects_labels2)+
  #scale_radius(range = c(0, 40), name="Growth rate relative\n between treatments")+
  scale_color_manual(values = palf_large_phylums(11))+
  #scale_size_identity()+
  #scale_alpha(range = c(0,40))+#, limits = c(0.1, 1)
  facet_wrap(vars(family_f))+ #, labeller = interaction(vars(class_f), vars(family_f), sep = ' ')
  labs(y = 'Treatments growth rates relatives', x = 'Season', color = 'Phylum')+
  theme_bw()+
  theme(legend.position = 'right', panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90), strip.text = element_text(size = 5), strip.background = element_blank())

ggsave('relative_gr_treatments_family2.pdf', relative_gr_treatments_family, 
       path = "~/Documentos/Doctorat/REMEI/results/figures/raw/",
       width = 240,
       height = 180,
       units = 'mm')


effects_difference_family_l %>%
  filter(difference != is.na(difference)) %>%
  subset(phylum == 'Proteobacteria') %>%
  group_by(class) %>%
  filter(n() > 10) %>%
  group_by(family) %>%
  filter(n() > 21) %>%
  ggplot(aes(season, effects, size = difference, color = class_f))+
  geom_point()+
  scale_y_discrete(labels = effects_labels2)+
  scale_size(range = c(-3.5, 4.5), name=expression("Difference in \ngrowth rate day"^"-1"))+
  scale_color_manual(values = palette_large)+
  facet_grid(~family_f)+
  labs(y = 'Growth rate difference between treatments', x = 'Season', color = 'Phylum')+
  theme_bw()+
  theme(legend.position = 'right', panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90))

relative_effects_family_l %>%
  filter(ratio != is.na(ratio)) %>%
  subset(phylum == 'Proteobacteria') %>%
  group_by(class) %>%
  filter(n() > 10) %>%
  # group_by(family) %>%
  # filter(n() > 21) %>%
  ggplot(aes(season, effects, size = ratio, color = class_f))+
  geom_point()+
  scale_y_discrete(labels = effects_labels)+
  scale_size(range = c(-3.5, 4.5), name=expression("Difference in \ngrowth rate day"^"-1"))+
  scale_color_manual(values = palette_large)+
  facet_grid(~family_f)+
  labs(y = 'Growth rate difference between treatments', x = 'Season', color = 'Phylum')+
  theme_bw()+
  theme(legend.position = 'right', panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90))


##order level bubble plot-----
effects_difference_order <- reg_all_slopes_chosen_silva_tax %>% 
  distinct(treatment, season, asv_num, domain, phylum, class, order, genus, species, tax_ed, asv_num, family, .keep_all = TRUE) %>%
  filter(pvalue_slope_chosen < 0.05 &
           slope_chosen_days > 0) %>%
  group_by(treatment, season, domain, phylum, class, order) %>%
  summarise(slope_chosen_days_mean = mean(slope_chosen_days),
            slope_chosen_mean = mean(slope_chosen),
            na.rm = TRUE) %>%
  pivot_wider(names_from = treatment, values_from = c(slope_chosen_days_mean, slope_chosen_mean), 
              id_cols = c(season, domain, phylum, class, order)) %>%
  as_tibble() %>%
  mutate(across(!c(domain, phylum, class, order, season), as.numeric)) %>%
  mutate(effect_top_down_virus = slope_chosen_days_mean_VL - slope_chosen_days_mean_DL,
         effect_top_down_grazers_dark = slope_chosen_days_mean_PD - slope_chosen_days_mean_CD,
         effect_top_down_grazers_light = slope_chosen_days_mean_PL - slope_chosen_days_mean_CL,
         effect_bottom_up = slope_chosen_days_mean_DL - slope_chosen_days_mean_PL,
         effect_light_C = slope_chosen_days_mean_CL - slope_chosen_days_mean_CD,
         effect_light_P = slope_chosen_days_mean_PL - slope_chosen_days_mean_PD)

effects_difference_order_l <- effects_difference_order %>%
  pivot_longer(cols = starts_with('effect'),
               names_to = 'effects',
               values_to = 'difference') %>%
  as_tibble()

effects_difference_order_l %>%
  colnames()

effects_difference_order %>%
  colnames()

effects_difference_order_l$season <-effects_difference_order_l$season %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))

effects_difference_order_l$effects <- effects_difference_order_l$effects %>% 
  factor(levels=(c("effect_bottom_up", "effect_top_down_grazers_dark", 
                   "effect_top_down_grazers_light", "effect_top_down_virus",
                   "effect_light_C", "effect_light_P")))

effects_difference_order_l <- effects_difference_order_l %>%
  mutate(phylum_f = as_factor(phylum),
         order_f = as_factor(order),
         class_f = as_factor(class))

effects_difference_order_l$class_f <-  factor(effects_difference_order_l$class_f, 
                                              levels=unique(effects_difference_order_l$class_f[order(effects_difference_order_l$phylum_f)]), 
                                              ordered=TRUE)

effects_difference_order_l$order_f <-  factor(effects_difference_order_l$order_f, 
                                              levels=unique(effects_difference_order_l$order_f[order(effects_difference_order_l$phylum_f,
                                                                                                     effects_difference_order_l$class_f)]), 
                                              ordered=TRUE)

difference_gr_treatments_order <- effects_difference_order_l %>%
  filter(difference != is.na(difference)) %>%
  group_by(class) %>%
  filter(n() > 10) %>%
  group_by(order) %>%
  filter(n() > 20) %>%
  group_by(order, difference) %>%
  mutate(counts = n()) %>%
  ggplot(aes(season, effects, size = difference, color = phylum_f, alpha = counts))+ #
  geom_point()+
  scale_y_discrete(labels = effects_labels2)+
  scale_size(range = c(-3.5, 4.5), name="Growth rate difference\n between treatments")+
  scale_alpha(range = c(1,1), limits = c(0.1, 1))+
  scale_color_manual(values = palette_phylums)+
  facet_wrap(vars(order_f))+ #, labeller = interaction(vars(class_f), vars(order_f), sep = ' ')
  labs(y = 'Treatments growth rates differences', x = 'Season', color = 'Phylum')+
  theme_bw()+
  theme(legend.position = 'right', panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90), strip.text = element_text(size = 5), strip.background = element_blank())

# ggsave('remiau_difference_gr_treatments_order2.pdf', difference_gr_treatments_order, 
#        path = "~/Documentos/Doctorat/REMEI/results/figures/raw/",
#        width = 240,
#        height = 180,
#        units = 'mm')

difference_gr_treatments_order <- 
  effects_difference_order_l %>%
  filter(difference != is.na(difference)) %>%
  group_by(class) %>%
  filter(n() > 10) %>%
  group_by(order) %>%
  filter(n() > 15) %>%
  group_by(order, difference) %>%
  mutate(counts = n()) %>%
  # ungroup() %>%
  # filter(difference > 0) %>%
  ggplot(aes(season, effects, color = difference))+ #, alpha = counts
  geom_point(aes(color = difference, size = difference))+
  scale_y_discrete(labels = effects_labels2)+
  scale_size(range = c(0, 12), name ="Growth rate difference\n between treatments")+
  scale_colour_gradientn(colours = palete_gradient)+
  # scale_color_manual(values = palette_phylums)+
  #scale_color_gradient2(low = '#B24D5D', high =  '#4DB2A2',  midpoint = 0)+
  #scale_alpha(range = c(1,1), limits = c(1, 1))+
  facet_wrap(vars(order_f))+ #, labeller = interaction(vars(class_f), vars(order_f), sep = ' ') #, labeller(order_f = ~ paste(order_f, class_f), .multi_line = TRUE)
  labs(y = 'Treatments growth rates differences', x = 'Season', color = '')+
  theme_bw()+
  theme(legend.position = 'right', panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 0), strip.text = element_text(size = 11), strip.background = element_blank(),
        axis.title = element_text(size = 18))

ggsave('remiau_difference_gr_treatments_order4.pdf', difference_gr_treatments_order, 
       path = "~/Documentos/Doctorat/REMEI/results/figures/raw/",
       width = 360,
       height = 290,
       units = 'mm')

##correlacions----
effects_difference_family$season <- effects_difference_family$season %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))


effects_difference_family <- effects_difference_family %>%
  mutate(phylum_f = as_factor(phylum),
         family_f = as_factor(family),
         order_f = as_factor(order),
         class_f = as_factor(class))

effects_difference_family$class_f <-  factor(effects_difference_family$class_f, 
                                             levels=unique(effects_difference_family$class_f[order(effects_difference_family$phylum_f)]), 
                                             ordered=TRUE)

effects_difference_family$order_f <-  factor(effects_difference_family$order_f, 
                                             levels=unique(effects_difference_family$order_f[order(effects_difference_family$phylum_f,
                                                                                                   effects_difference_family$class_f)]), 
                                             ordered=TRUE)

effects_difference_family$family_f <-  factor(effects_difference_family$family_f, 
                                              levels=unique(effects_difference_family$family_f[order(effects_difference_family$phylum_f,
                                                                                                     effects_difference_family$class_f,
                                                                                                     effects_difference_family$order_f)]), 
                                              ordered=TRUE)


p1<- effects_difference_family %>%
  ggplot(aes(effect_top_down_grazers_light, effect_bottom_up, color = phylum))+
  scale_x_continuous(limits = c(-3, 3))+
  scale_y_continuous(limits = c(-3, 3))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  labs(x = 'Top-down effect grazers light (PL-CL)', y = 'Bottom-up effect (DL-PL)')+
  geom_point(size = 3, alpha =  0.8)+ #mida hauria de ser el número de OTUs dins de cada família
  scale_color_manual(values = palette_phylums_assigned)+
  facet_wrap(vars(season), nrow = 1)+
  theme_bw()+
  theme(legend.position = 'none', panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_blank())

p2 <- effects_difference_family %>%
  ggplot(aes(effect_top_down_virus, effect_bottom_up, color = phylum))+
  scale_x_continuous(limits = c(-3, 3))+
  scale_y_continuous(limits = c(-3, 3))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  scale_color_manual(values = palette_phylums_assigned)+
  labs(x = 'Top-down effect viruses (VL-DL)', y = 'Bottom-up effect (DL-PL)')+
  geom_point(size = 3, alpha =  0.8)+ #mida hauria de ser el número de OTUs dins de cada família
  facet_wrap(vars(season), nrow = 1)+
  theme_bw()+
  theme(legend.position = 'none', panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_blank(),
        strip.text = element_blank())

p3 <- effects_difference_family %>%
  ggplot(aes(effect_top_down_virus, effect_top_down_grazers_light, color = phylum))+
  scale_x_continuous(limits = c(-3, 3))+
  scale_y_continuous(limits = c(-3, 3))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  geom_point(size = 3, alpha =  0.8)+ #mida hauria de ser el número de OTUs dins de cada família
  scale_color_manual(values = palette_phylums_assigned)+
  labs(x = 'Top-down effect viruses (VL-DL)', y = 'Top-down effect grazers light (PL-CL)')+
  facet_wrap(vars(season), nrow = 1)+
  theme_bw()+
  theme(legend.position = 'none', panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_blank(),
        strip.text = element_blank())

p4 <- effects_difference_family %>%
  ggplot(aes(effect_light_C, effect_light_P, color = phylum))+
  scale_x_continuous(limits = c(-3, 3))+
  scale_y_continuous(limits = c(-3, 3))+
  labs(x = 'Light effect (CL-CD)', y = 'Light effect (PL-PD)')+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  #stat_ellipse(aes(effect_light_C, effect_light_P,group = domain))+
  geom_point(size = 3, alpha =  0.8)+ #mida hauria de ser el número de OTUs dins de cada família
  scale_color_manual(values = palette_phylums_assigned)+
  facet_wrap(vars(season), nrow = 1)+
  theme_bw()+
  theme(legend.position = 'none', panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_blank(),
        strip.text = element_blank())

p4 %>%
  ggExtra::ggMarginal(aes(fill = phylum_f), groupFill = TRUE, type = "density")+#, color="grey") ##per afegir histograma al voltant del gràfic
  scale_fill_manual(values = palette_phylums_assigned)

grid.arrange(p1, p2, p3, p4, ncol = 1)

##el mateix però color per classes
p1<- effects_difference_family %>%
  group_by(season, phylum_f) %>%
  filter(n() > 2) %>%
  group_by(season, class_f) %>%
  filter(n() > 2) %>%
  ggplot(aes(effect_top_down_grazers_light, effect_bottom_up, color = class_f))+
  scale_x_continuous(limits = c(-3, 3))+
  scale_y_continuous(limits = c(-3, 3))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  labs(x = 'Top-down effect grazers light (PL-CL)', y = 'Bottom-up effect (DL-PL)')+
  geom_point(size = 4, alpha =  0.8)+ #mida hauria de ser el número de OTUs dins de cada família
  scale_color_manual(values = palf_large_phylums(10))+
  facet_wrap(vars(season), nrow = 1)+
  theme_bw()+
  theme(legend.position = 'none', panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_blank(), 
        axis.title = element_text(size = 16))

p2 <- effects_difference_family %>%
  group_by(season, phylum_f) %>%
  filter(n() > 2) %>%
  group_by(season, class_f) %>%
  filter(n() > 2) %>%
  ggplot(aes(effect_top_down_virus, effect_bottom_up, color = class_f))+
  scale_x_continuous(limits = c(-3, 3))+
  scale_y_continuous(limits = c(-3, 3))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  scale_color_manual(values = palf_large_phylums(10))+
  labs(x = 'Top-down effect viruses (VL-DL)', y = 'Bottom-up effect (DL-PL)')+
  geom_point(size = 4, alpha =  0.8)+ #mida hauria de ser el número de OTUs dins de cada família
  facet_wrap(vars(season), nrow = 1)+
  theme_bw()+
  theme(legend.position = 'none', panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_blank(),
        strip.text = element_blank(), axis.title = element_text(size = 16))

p3 <- effects_difference_family %>%
  group_by(season, phylum_f) %>%
  filter(n() > 2) %>%
  group_by(season, class_f) %>%
  filter(n() > 2) %>%
  ggplot(aes(effect_top_down_virus, effect_top_down_grazers_light, color = class_f))+
  scale_x_continuous(limits = c(-3, 3))+
  scale_y_continuous(limits = c(-3, 3))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  geom_point(size = 4, alpha =  0.8)+ #mida hauria de ser el número de OTUs dins de cada família
  scale_color_manual(values = palf_large_phylums(10))+
  labs(x = 'Top-down effect viruses (VL-DL)', y = 'Top-down effect grazers light (PL-CL)')+
  facet_wrap(vars(season), nrow = 1)+
  theme_bw()+
  theme(legend.position = 'none', panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_blank(),
        strip.text = element_blank(), axis.title = element_text(size = 16))

p4 <- effects_difference_family %>%
  group_by(season, phylum_f) %>%
  filter(n() > 2) %>%
  group_by(season, class_f) %>%
  filter(n() > 2) %>%
  ggplot(aes(effect_light_C, effect_light_P, color = class_f))+
  scale_x_continuous(limits = c(-3, 3))+
  scale_y_continuous(limits = c(-3, 3))+
  labs(x = 'Light effect (CL-CD)', y = 'Light effect (PL-PD)', color = 'Taxonomic\n class')+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  geom_point(size = 4, alpha =  0.8)+ #mida hauria de ser el número de OTUs dins de cada família
  scale_color_manual(values = palf_large_phylums(10))+
  facet_wrap(vars(season), nrow = 1)+
  guides(color=guide_legend(ncol = 1))+
  theme_bw()+
  theme(legend.position = 'bottom', panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), strip.background = element_blank(), legend.text = element_text(size = 20),
        legend.title = element_text(size = 24), axis.title = element_text(size = 16))

#grid.arrange(p1, p2, p3, p4, ncol = 1)

correlations_differences_family_treatments_class  <- multi_panel_figure(columns = 1, rows = 5, width = 400, height = 450, 
                                                                        row_spacing = 2, unit = 'mm',
                                                                        panel_label_type = 'none')

correlations_differences_family_treatments_class  %<>%
  fill_panel(p1, row = 1, col = 1) %<>%
  fill_panel(p2, row = 2, col = 1) %<>%
  fill_panel(p3, row = 3, col = 1) %<>%
  fill_panel(p4, row = 4:5, col = 1)

ggsave('remiau_correlations_differences_family_treatments_class.pdf', correlations_differences_family_treatments_class, 
       path = "~/Documentos/Doctorat/REMEI/results/figures/raw/",
       width = 400,
       height = 450,
       units = 'mm')


##order level
p1<- effects_difference_family %>%
  group_by(season, phylum_f) %>%
  filter(n() > 2) %>%
  group_by(season, class_f) %>%
  filter(n() > 2) %>%
  ggplot(aes(effect_top_down_grazers_light, effect_bottom_up, color = order))+
  scale_x_continuous(limits = c(-3, 3))+
  scale_y_continuous(limits = c(-3, 3))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  labs(x = 'Top-down effect grazers light (PL-CL)', y = 'Bottom-up effect (DL-PL)')+
  geom_point(size = 3, alpha =  0.8)+ #mida hauria de ser el número de OTUs dins de cada família
  scale_color_manual(values = palf_large_phylums(66))+
  facet_wrap(vars(season), nrow = 1)+
  theme_bw()+
  theme(legend.position = 'none', panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_blank())

p2 <- effects_difference_family %>%
  group_by(season, phylum_f) %>%
  filter(n() > 2) %>%
  group_by(season, class_f) %>%
  filter(n() > 2) %>%
  ggplot(aes(effect_top_down_virus, effect_bottom_up, color = order))+
  scale_x_continuous(limits = c(-3, 3))+
  scale_y_continuous(limits = c(-3, 3))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  scale_color_manual(values = palf_large_phylums(66))+
  labs(x = 'Top-down effect viruses (VL-DL)', y = 'Bottom-up effect (DL-PL)')+
  geom_point(size = 3, alpha =  0.8)+ #mida hauria de ser el número de OTUs dins de cada família
  facet_wrap(vars(season), nrow = 1)+
  theme_bw()+
  theme(legend.position = 'none', panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_blank(),
        strip.text = element_blank())

p3 <- effects_difference_family %>%
  group_by(season, phylum_f) %>%
  filter(n() > 2) %>%
  group_by(season, class_f) %>%
  filter(n() > 2) %>%
  ggplot(aes(effect_top_down_virus, effect_top_down_grazers_light, color = order))+
  scale_x_continuous(limits = c(-3, 3))+
  scale_y_continuous(limits = c(-3, 3))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  geom_point(size = 3, alpha =  0.8)+ #mida hauria de ser el número de OTUs dins de cada família
  scale_color_manual(values = palf_large_phylums(66))+
  labs(x = 'Top-down effect viruses (VL-DL)', y = 'Top-down effect grazers light (PL-CL)')+
  facet_wrap(vars(season), nrow = 1)+
  theme_bw()+
  theme(legend.position = 'none', panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_blank(),
        strip.text = element_blank())

p4 <- effects_difference_family %>%
  group_by(season, phylum_f) %>%
  filter(n() > 2) %>%
  group_by(season, class_f) %>%
  filter(n() > 2) %>%
  ggplot(aes(effect_light_C, effect_light_P, color = interaction(order_f, class_f, sep = ' ')))+ #interaction(class_f, family_f,  sep = ' ')
  scale_x_continuous(limits = c(-3, 3))+
  scale_y_continuous(limits = c(-3, 3))+
  labs(x = 'Light effect (CL-CD)', y = 'Light effect (PL-PD)', color = 'Order and Class')+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  guides(color=guide_legend(ncol = 3))+
  geom_point(size = 3, alpha =  0.8)+ #mida hauria de ser el número de OTUs dins de cada família
  scale_color_manual(values = palf_large_phylums(66))+
  facet_wrap(vars(season), nrow = 1)+
  theme_bw()+
  theme(legend.position = 'bottom', panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), strip.background = element_blank(), legend.text = element_text(size = 14))

#grid.arrange(p1, p2, p3, p4, ncol = 1)

correlations_differences_family_treatments  <- multi_panel_figure(columns = 1, rows = 6, width = 400, height = 450, 
                                                                  row_spacing = 2, unit = 'mm',
                                                                  panel_label_type = 'none')

correlations_differences_family_treatments  %<>%
  fill_panel(p1, row = 1, col = 1) %<>%
  fill_panel(p2, row = 2, col = 1) %<>%
  fill_panel(p3, row = 3, col = 1) %<>%
  fill_panel(p4, row = 4:6, col = 1)

ggsave('remiau_correlations_differences_family_treatments.pdf', correlations_differences_family_treatments, 
       path = "~/Documentos/Doctorat/REMEI/results/figures/raw/",
       width = 400,
       height = 450,
       units = 'mm')

##correlogram
library(corrr)
effects_difference_family_l %>%
  colnames()

effects_difference_family %>%
  select(effect_top_down_virus, effect_top_down_grazers_dark, effect_top_down_grazers_light,
         effect_bottom_up, effect_light_C, effect_light_P) %>%
  correlate() %>%
  rplot()

##heatmap
effects_difference_family_l %>%
  filter(difference != is.na(difference)) %>%
  group_by(season, family_f) %>%
  filter(n() > 4) %>%
  ggplot(aes(effects, family_f, fill = difference))+
  geom_tile(color = 'white', lwd = 1.5, linetype = 0.5)+
  #geom_tile(aes(fill = ratio), color = 'white', lwd = 1.5, linetype = 0.5)+
  scale_fill_gradient2(low = '#BD4247', high =  '#42BDB8', na.value = '#fdfdfd', midpoint = 0)+
  #geom_point()+ #mida hauria de ser el número de OTUs dins de cada família
  #scale_color_manual(values = palf_large(100))+
  #scale_x_continuous(limits = c(0, 3))+
  #scale_y_continuous(limits = c(0, 3))+
  scale_x_discrete(labels = effects_labels)+
  #facet_wrap(vars(season), nrow = 1)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

##pels que están afectats per la llum
effects_difference_family_l %$%
  effects
effects_difference_family_l %>%
  filter(difference != is.na(difference)) %>%
  subset(effects == c('effect_light_P', 
                      'effect_light_C')) %>%
  group_by(season, family_f) %>%
  filter(n() > 1) %>%
  filter(difference > 0) %>%
  ggplot(aes(effects, interaction(class_f, family_f,  sep = ' '), fill = difference))+
  geom_tile(color = 'white', lwd = 1.5, linetype = 0.5)+
  #geom_tile(aes(fill = ratio), color = 'white', lwd = 1.5, linetype = 0.5)+
  scale_fill_gradient2(low = '#BD4247', high =  '#42BDB8', na.value = '#fdfdfd', midpoint = 0)+
  #geom_point()+ #mida hauria de ser el número de OTUs dins de cada família
  #scale_color_manual(values = palf_large(100))+
  #scale_x_continuous(limits = c(0, 3))+
  #scale_y_continuous(limits = c(0, 3))+
  scale_x_discrete(labels = effects_labels)+
  facet_wrap(vars(season), nrow = 1)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


##filtrem només efectes positius
negatiu <- effects_difference_family_l %>%
  filter(difference != is.na(difference)) %>%
  filter(difference > 0) %>%
  group_by(season, family_f) %>%
  filter(n() > 3) %>%
  ggplot(aes(effects, interaction(family_f, class_f, sep = ' '), fill = difference))+
  geom_tile(color = 'white', lwd = 1.5, linetype = 0.5)+
  #geom_tile(aes(fill = ratio), color = 'white', lwd = 1.5, linetype = 0.5)+
  scale_fill_gradient2(low = '#BD4247', high =  '#42BDB8', na.value = '#fdfdfd', midpoint = 0)+
  #geom_point()+ #mida hauria de ser el número de OTUs dins de cada família
  #scale_color_manual(values = palf_large(100))+
  #scale_x_continuous(limits = c(0, 3))+
  #scale_y_continuous(limits = c(0, 3))+
  scale_x_discrete(labels = effects_labels)+
  facet_wrap(vars(season), nrow = 1)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

##filtrem només efectes negatius
positiu <- effects_difference_family_l %>%
  filter(difference != is.na(difference)) %>%
  filter(difference < 0) %>%
  group_by(season, family_f) %>%
  filter(n() > 3) %>%
  ggplot(aes(effects, interaction(family_f, class_f, sep = ' '), fill = difference))+
  geom_tile(color = 'white', lwd = 1.5, linetype = 0.5)+
  #geom_tile(aes(fill = ratio), color = 'white', lwd = 1.5, linetype = 0.5)+
  scale_fill_gradient2(low = '#BD4247', high =  '#42BDB8', na.value = '#fdfdfd', midpoint = 0)+
  #geom_point()+ #mida hauria de ser el número de OTUs dins de cada família
  #scale_color_manual(values = palf_large(100))+
  #scale_x_continuous(limits = c(0, 3))+
  #scale_y_continuous(limits = c(0, 3))+
  scale_x_discrete(labels = effects_labels)+
  facet_wrap(vars(season), nrow = 1)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

grid.arrange(positiu, negatiu, ncol = 1)

##comparació control-predator only
effects_difference_family_l %$%
  effects %>%
  unique()

effects_difference_family_l %>%
  filter(difference != is.na(difference)) %>%
  subset(effects == c('effect_top_down_grazers_dark', 
                      'effect_top_down_grazers_light')) %>%
  group_by(season, family_f) %>%
  filter(n() > 1) %>%
  #filter(difference > 0) %>%
  ggplot(aes(effects, interaction(class_f, family_f,  sep = ' '), fill = difference))+
  geom_tile(color = 'white', lwd = 1.5, linetype = 0.5)+
  #geom_tile(aes(fill = ratio), color = 'white', lwd = 1.5, linetype = 0.5)+
  scale_fill_gradient2(low = '#BD4247', high =  '#42BDB8', na.value = '#fdfdfd', midpoint = 0)+
  #geom_point()+ #mida hauria de ser el número de OTUs dins de cada família
  #scale_color_manual(values = palf_large(100))+
  #scale_x_continuous(limits = c(0, 3))+
  #scale_y_continuous(limits = c(0, 3))+
  scale_x_discrete(labels = effects_labels)+
  facet_wrap(vars(season), nrow = 1)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


