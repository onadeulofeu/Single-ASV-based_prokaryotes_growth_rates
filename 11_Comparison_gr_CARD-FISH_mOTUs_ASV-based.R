library(tidyverse)
library(magrittr)
library(ggplot2)
library(multipanelfigure)
library(ggpmisc) ##stat_poli_line and equation
library(ggpubr) #as_ggplot function, stat_cor
library(stringi) 
library(magrittr)
library(Hmisc) #for correlations

#Import data-----
## mOTUs-based gr
motus_env_l_w_fold_filt <- read.table('data/intermediate_files/motus_env_l_w_fold_filt.txt', sep = '\t')
## ASV-based gr
reg_all_slopes_chosen_silva_tax_filt <- read.csv("data/intermediate_files/reg_all_slopes_chosen_silva_tax.csv", sep=",") %>%
  filter(season != "Early_fall") %>%
  filter(slope_chosen_days > 0 &
           pvalue_slope_chosen < 0.05)
## abundance-based gr
abundance_gr <- read.delim2("data/envdata/GR_DAPIS_OS/GR_REMEI_DAPIS_OS_Ed.csv", sep = ";") %>%
  as_tibble()

abundance_gr <- GR_dapis_OS %>%
  filter(treatment %in% c("CL", "CD", "PL", "PD",  "DL", "VL")) %>%
  mutate(ALT = as.numeric(ALT),
         ROSE = as.numeric(ROSEO),
         GAMMA = as.numeric(GAMMA),
         NOR5= as.numeric(NOR5),
         CFB = as.numeric(CFB),
         SAR11 = as.numeric(SAR11)) 

##Correlations with groups from FISH-probes:ALT1413: Colwellia, Alteromonas --------
mOTUS_alteromonas <-  motus_env_l_w_fold_filt %>%
  # motus_env_l_w_h_growth_1perc %>% (només els que tenen més d'un 1%)
  dplyr::filter(gr_day > 0 &
                  consensus_taxonomy_mOTUs_num != 'unassigned') %>% 
  dplyr::filter(str_detect(consensus_taxonomy_mOTUs_num, "Alteromonas") | 
                      str_detect(consensus_taxonomy_mOTUs_num, "Colwellia")
  ) %>%
  group_by(season, treatment) %>%
  summarize(mean_motus = mean(gr_day))

asv_alteromonas <-  reg_all_slopes_chosen_silva_tax %>%
  dplyr::filter(season != "Early_fall") %>%
  dplyr::filter(slope_chosen_days > 0,
                pvalue_slope_chosen < 0.05) %>%
  dplyr::filter(str_detect(genus_f, "Alteromonas") | 
                str_detect(genus_f,'Colwellia')) %>%
  group_by(season, treatment) %>%
  summarize(mean_asv = mean(slope_chosen_days))

dapis_alt <-  GR_dapis_OS_filt %>%
  select(treatment, replicate, season, ALT) %>%
  group_by(treatment, season) %>%
  summarize(mean_dapis_alt = mean(ALT)) %>%
  filter(mean_dapis_alt >0)

gr_alt <- mOTUS_alteromonas %>%
  right_join(asv_alteromonas) %>%
  left_join(dapis_alt) 

gr_alt %>%
  ggplot(aes(mean_asv, mean_motus))+
  geom_point(aes(shape = treatment, color = season, size = 2))+
  scale_fill_manual(values = palette_seasons_4)+
  #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
  # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  #stat_poly_eq(color = 'black')+
  stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
  guides(size = 'none')+
  labs(x = 'Mean growth rates ASV based \nAlteromonas & Colwellia', 
       y = 'Mean growth rates mOTUs-based Alteromonas & Colwellia', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/3), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

gr_alt %>%
  #filter(season != 'Summer') %>%
  ggplot(aes(mean_dapis_alt, mean_asv))+
  geom_point(aes(shape = treatment, color = season, size = 2))+
  scale_fill_manual(values = palette_seasons_4)+
  #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
  # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  #stat_poly_eq(color = 'black')+
  stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
  guides(size = 'none')+
  labs(y = 'Mean growth rates ASV based\nAlteromonas & Colwellia', 
       x = 'Mean growth rates abundance-based\nAlteromonas & Colwellia', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/3), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

gr_alt %>%
  ggplot(aes(mean_dapis_alt, mean_motus))+
  geom_point(aes(shape = treatment, color = season, size = 2))+
  scale_fill_manual(values = palette_seasons_4)+
  #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
  # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  #stat_poly_eq(color = 'black')+
  stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
  guides(size = 'none')+
  labs(y = 'Mean growth rates mOTUs-based \n
       Alteromonas & Colwellia', x = 'Mean growth rates abundance-based Alteromonas & Colwellia', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/3), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

##gammaproteobacteria-------
mOTUS_gamma <-  motus_env_l_w_fold %>%
  # motus_env_l_w_h_growth_1perc %>% (només els que tenen més d'un 1%)
  dplyr::filter(gr_day > 0 &
                  consensus_taxonomy_mOTUs_num != 'unassigned') %>% 
  dplyr::filter(str_detect(consensus_taxonomy_mOTUs_num, "Gammaproteobacteria"),
                str_detect(consensus_taxonomy_mOTUs_num, "")) %>%
  group_by(season, treatment) %>%
  summarize(mean_motus = mean(gr_day))

asv_gamma <-  reg_all_slopes_chosen_silva_tax %>%
  dplyr::filter(season != "Early_fall") %>%
  dplyr::filter(slope_chosen_days > 0,
                pvalue_slope_chosen < 0.05) %>%
  dplyr::filter(str_detect(class_f, "Gammaproteobacteria")) %>%
  group_by(season, treatment) %>%
  summarize(mean_asv = mean(slope_chosen_days))

dapis_gamma <-  GR_dapis_OS_filt %>%
  dplyr::select(treatment, replicate, season, GAMMA) %>%
  group_by(treatment, season) %>%
  dplyr::summarize(mean_dapis_gamma = mean(as.numeric(GAMMA))) %>%
  dplyr::filter(mean_dapis_gamma >0)

gamma_gr <- mOTUS_gamma %>%
  right_join(asv_gamma) %>%
  left_join(dapis_gamma) 

gamma_gr %>%
  ggplot(aes(mean_asv, mean_motus))+
  geom_point(aes(shape = treatment, color = season, size = 2))+
  scale_fill_manual(values = palette_seasons_4)+
  #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
  # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  #stat_poly_eq(color = 'black')+
  stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
  guides(size = 'none')+
  labs(x = 'Mean growth rates ASV based \n
       Gammaproteobacteria', y = 'Mean growth rates mOTUs based Gammaproteobacteria', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/3), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

gamma42a_asv_gr_filt <- gamma_gr %>%
  filter(season != 'Summer') %>%
  ggplot(aes(mean_dapis_gamma, mean_asv))+
  geom_point(aes(shape = treatment, color = season, size = 2))+
  scale_fill_manual(values = palette_seasons_4)+
  #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
  # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  #stat_poly_eq(color = 'black')+
  stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
  guides(size = 'none')+
  labs(y = 'Mean growth rates ASV-based \nGammaproteobacteria', 
       x = 'Mean growth rates abundance-based Gammaproteobacteria', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/3), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

gamma_gr %>%
  ggplot(aes(mean_dapis_gamma, mean_motus))+
  geom_point(aes(shape = treatment, color = season, size = 2))+
  scale_fill_manual(values = palette_seasons_4)+
  #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
  # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  #stat_poly_eq(color = 'black')+
  stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
  guides(size = 'none')+
  labs(y = 'Mean growth rates mOTUs based Gammaproteobacteria',
       x = 'Mean growth rates abundance based Gammaproteobacteria', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/3), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

##Rhodobacteraceae-----
##Rhobobacter, Roseobacter
mOTUS_rhodo <-  motus_env_l_w_fold %>%
  # motus_env_l_w_h_growth_1perc %>% (només els que tenen més d'un 1%)
  dplyr::filter(gr_day > 0 &
                  consensus_taxonomy_mOTUs_num != 'unassigned') %>% 
  dplyr::filter(str_detect(consensus_taxonomy_mOTUs_num, "Rhodobacter") |
                  str_detect(consensus_taxonomy_mOTUs_num, "Roseobacter")) %>%
  group_by(season, treatment) %>%
  summarize(mean_motus = mean(gr_day))

asv_rhodo <-  reg_all_slopes_chosen_silva_tax %>%
  dplyr::filter(season != "Early_fall") %>%
  dplyr::filter(slope_chosen_days > 0,
                pvalue_slope_chosen < 0.05) %>%
  dplyr::filter(str_detect(genus_f, "Rhodobacter") |
                  str_detect(genus_f, 'Roseobacter')) %>% #clade NAC11-7 lineage
  group_by(season, treatment) %>%
  summarize(mean_asv = mean(slope_chosen_days))

dapis_rhodo <-  GR_dapis_OS_filt %>%
  dplyr::select(treatment, replicate, season, ROSEO) %>%
  group_by(treatment, season) %>%
  dplyr::summarize(mean_dapis_rhodo = mean(as.numeric(ROSEO))) %>%
  dplyr::filter(mean_dapis_rhodo >0)

rhodo_gr <- mOTUS_rhodo %>%
  left_join(asv_rhodo) %>%
  right_join(dapis_rhodo) 

rhodo_gr %>%
  ggplot(aes(mean_asv, mean_motus))+
  geom_point(aes(shape = treatment, color = season, size = 2))+
  scale_fill_manual(values = palette_seasons_4)+
  #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
  # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  #stat_poly_eq(color = 'black')+
  stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
  guides(size = 'none')+
  labs(x = 'Mean growth rates ASV based \n
       Roseobacter/Rhodobacter', y = 'Mean growth rates mOTUs based Roseobacter/Rhodobacter', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/3), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

rhodo_gr %>%
  ggplot(aes(mean_dapis_rhodo, mean_asv))+
  geom_point(aes(shape = treatment, color = season, size = 2))+
  scale_fill_manual(values = palette_seasons_4)+
  #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
  # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  #stat_poly_eq(color = 'black')+
  stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
  guides(size = 'none')+
  labs(y = 'Mean growth rates ASV based \n
       Roseobacter/Rhodobacter', x = 'Mean growth rates abundance-based Roseobacter/Rhodobacter', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/3), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

rhodo_gr %>%
  ggplot(aes(mean_dapis_rhodo, mean_motus))+
  geom_point(aes(shape = treatment, color = season, size = 2))+
  scale_fill_manual(values = palette_seasons_4)+
  #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
  # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  #stat_poly_eq(color = 'black')+
  stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
  guides(size = 'none')+
  labs(y = 'Mean growth rates mOTUs-based Roseobacter/Rhodobacter',
       x = 'Mean growth rates abundance-based Roseobacter/Rhodobacter', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/4), panel.grid.major = element_blank(), panel.grid.minor = element_blank())


###SAR11 clade----
mOTUS_sar11 <-  motus_env_l_w_fold %>%
  # motus_env_l_w_h_growth_1perc %>% (només els que tenen més d'un 1%)
  dplyr::filter(gr_day > 0 &
                  consensus_taxonomy_mOTUs_num != 'unassigned') %>% 
  dplyr::filter(str_detect(consensus_taxonomy_mOTUs_num, "SAR11 clade")) %>% ##no els troba!
  group_by(season, treatment) %>%
  summarize(mean_motus = mean(gr_day))

asv_sar11 <-  reg_all_slopes_chosen_silva_tax %>%
  dplyr::filter(season != "Early_fall") %>%
  dplyr::filter(slope_chosen_days > 0,
                pvalue_slope_chosen < 0.05) %>%
  dplyr::filter(str_detect(order_f, "SAR11 clade")) %>%
  group_by(season, treatment) %>%
  summarize(mean_asv = mean(slope_chosen_days))

dapis_sar11 <-  GR_dapis_OS_filt %>%
  dplyr::select(treatment, replicate, season, SAR11) %>%
  group_by(treatment, season) %>%
  dplyr::summarize(mean_dapis_sar11 = mean(as.numeric(SAR11))) %>%
  dplyr::filter(mean_dapis_sar11 >0)

sar11_gr <- #mOTUS_sar11 %>%
  asv_sar11 %>%
  left_join(dapis_sar11) 

# sar11_gr %>%
#   ggplot(aes(mean_asv, mean_motus))+
#   geom_point(aes(shape = treatment, color = season, size = 2))+
#   scale_fill_manual(values = palette_seasons_4)+
#   #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
#   # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
#   scale_color_manual(values = palette_seasons_4)+
#   stat_poly_line(color = 'black')+
#   #stat_poly_eq(color = 'black')+
#   stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
#   guides(size = 'none')+
#   labs(x = 'Mean growth rates ASV based \n
#        SAR11 clade', y = 'Mean growth rates mOTUs based SAR11 clade', color = 'Season', shape = 'Treatment')+
#   theme_bw()+
#   theme(aspect.ratio = (4/3), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

sar11_clade_asv_gr_filt <- sar11_gr %>%
  filter(season != 'Summer') %>%
  ggplot(aes(mean_dapis_sar11, mean_asv))+
  geom_point(aes(shape = treatment, color = season, size = 2))+
  scale_fill_manual(values = palette_seasons_4)+
  #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
  # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  #stat_poly_eq(color = 'black')+
  stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
  guides(size = 'none')+
  labs(y = 'Mean growth rates ASV-based \nSAR11 clade', 
       x = 'Mean growth rates abundance-based SAR11 clade', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/4), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# sar11_gr %>%
#   ggplot(aes(mean_dapis_sar11, mean_motus))+
#   geom_point(aes(shape = treatment, color = season, size = 2))+
#   scale_fill_manual(values = palette_seasons_4)+
#   #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
#   # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
#   scale_color_manual(values = palette_seasons_4)+
#   stat_poly_line(color = 'black')+
#   #stat_poly_eq(color = 'black')+
#   stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
#   guides(size = 'none')+
#   labs(y = 'Mean growth rates mOTUs based SAR11 clade',
#        x = 'Mean growth rates abundance based SAR11 clade', color = 'Season', shape = 'Treatment')+
#   theme_bw()+
#   theme(aspect.ratio = (4/3), panel.grid.major = element_blank(), panel.grid.minor = element_blank())


##CFB bacteroidetes----
mOTUS_cfb <-  motus_env_l_w_fold %>%
  # motus_env_l_w_h_growth_1perc %>% (només els que tenen més d'un 1%)
  dplyr::filter(gr_day > 0 &
                  consensus_taxonomy_mOTUs_num != 'unassigned') %>% 
  dplyr::filter(str_detect(consensus_taxonomy_mOTUs_num, "Bacteroidetes")) %>% ##no els troba!
  group_by(season, treatment) %>%
  summarize(mean_motus = mean(gr_day))

asv_cfb <-  reg_all_slopes_chosen_silva_tax %>%
  dplyr::filter(season != "Early_fall") %>%
  dplyr::filter(slope_chosen_days > 0,
                pvalue_slope_chosen < 0.05) %>%
  dplyr::filter(str_detect(phylum_f, "Bacteroidota")) %>%
  group_by(season, treatment) %>%
  summarize(mean_asv = mean(slope_chosen_days))

dapis_cfb <-  GR_dapis_OS_filt %>%
  dplyr::select(treatment, replicate, season, CFB) %>%
  group_by(treatment, season) %>%
  dplyr::summarize(mean_dapis_fcb = mean(as.numeric(CFB))) %>%
  dplyr::filter(mean_dapis_fcb >0)

cfb_gr <- mOTUS_cfb %>%
  right_join(asv_cfb) %>%
  left_join(dapis_cfb) 

# fcb_gr %>%
#   ggplot(aes(mean_asv, mean_motus))+
#   geom_point(aes(shape = treatment, color = season, size = 2))+
#   scale_fill_manual(values = palette_seasons_4)+
#   #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
#   # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
#   scale_color_manual(values = palette_seasons_4)+
#   stat_poly_line(color = 'black')+
#   #stat_poly_eq(color = 'black')+
#   stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
#   guides(size = 'none')+
#   labs(x = 'Mean growth rates ASV based \n
#        fcb clade', y = 'Mean growth rates mOTUs based fcb clade', color = 'Season', shape = 'Treatment')+
#   theme_bw()+
#   theme(aspect.ratio = (4/3), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

cf319a_asv_gr_filt <- cfb_gr %>%
  filter(season != 'Summer') %>%
  ggplot(aes(mean_dapis_fcb, mean_asv))+
  geom_point(aes(shape = treatment, color = season, size = 2))+
  scale_fill_manual(values = palette_seasons_4)+
  #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
  # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  #stat_poly_eq(color = 'black')+
  stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
  guides(size = 'none')+
  labs(y = 'Mean growth rates ASV based Bacteroidota', x = 'Mean growth rates abundance-based Bacteroidota', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/3), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# fcb_gr %>%
#   ggplot(aes(mean_dapis_fcb, mean_motus))+
#   geom_point(aes(shape = treatment, color = season, size = 2))+
#   scale_fill_manual(values = palette_seasons_4)+
#   #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
#   # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
#   scale_color_manual(values = palette_seasons_4)+
#   stat_poly_line(color = 'black')+
#   #stat_poly_eq(color = 'black')+
#   stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
#   guides(size = 'none')+
#   labs(y = 'Mean growth rates mOTUs based fcb clade',
#        x = 'Mean growth rates abundance based fcb clade', color = 'Season', shape = 'Treatment')+
#   theme_bw()+
#   theme(aspect.ratio = (4/3), panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# ##reorder treatments and seasons for plots
# GR_dapis_OS_filt$treatment <- GR_dapis_OS_filt$treatment %>% 
#   factor(levels=(c("CL", "CD", "PL", "PD", "DL", "VL")))
# GR_dapis_OS_filt$season <- GR_dapis_OS_filt$season %>% 
#   factor(levels=(c("Winter", "Spring", "Summer", "Fall")))

##NOR5-OM60------
mOTUS_nor5 <-  motus_env_l_w_fold %>%
  # motus_env_l_w_h_growth_1perc %>% (només els que tenen més d'un 1%)
  dplyr::filter(gr_day > 0 &
                  consensus_taxonomy_mOTUs_num != 'unassigned') %>% 
  dplyr::filter(str_detect(consensus_taxonomy_mOTUs_num, "nor5") |
                  str_detect(consensus_taxonomy_mOTUs_num, "nor5")) %>%
  group_by(season, treatment) %>%
  summarize(mean_motus = mean(gr_day))

asv_nor5 <-  reg_all_slopes_chosen_silva_tax %>%
  dplyr::filter(season != "Early_fall") %>%
   dplyr::filter(slope_chosen_days > 0,
                 pvalue_slope_chosen < 0.05) %>%
  dplyr::filter(str_detect(genus_f, "OM60")) %>%
  group_by(season, treatment) %>%
  summarize(mean_asv = mean(slope_chosen_days))

dapis_nor5 <-  GR_dapis_OS_filt %>%
  dplyr::select(treatment, replicate, season, NOR5) %>%
  group_by(treatment, season) %>%
  dplyr::summarize(mean_dapis_nor5 = mean(as.numeric(NOR5))) %>%
  dplyr::filter(mean_dapis_nor5 >0)

nor5_gr <- mOTUS_nor5 %>%
  right_join(asv_nor5) %>%
  right_join(dapis_nor5) 

# nor5_gr %>%
#   ggplot(aes(mean_asv, mean_motus))+
#   geom_point(aes(shape = treatment, color = season, size = 2))+
#   scale_fill_manual(values = palette_seasons_4)+
#   #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
#   # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
#   scale_color_manual(values = palette_seasons_4)+
#   stat_poly_line(color = 'black')+
#   #stat_poly_eq(color = 'black')+
#   stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
#   guides(size = 'none')+
#   labs(x = 'Mean growth rates ASV based \n
#        Roseobacter/nor5bacter', y = 'Mean growth rates mOTUs based Roseobacter/nor5bacter', color = 'Season', shape = 'Treatment')+
#   theme_bw()+
#   theme(aspect.ratio = (4/3), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

nor5_asv_gr_filt <- nor5_gr %>%
  filter(season != 'Summer') %>%
  ggplot(aes(mean_dapis_nor5, mean_asv))+
  geom_point(aes(shape = treatment, color = season, size = 2))+
  scale_fill_manual(values = palette_seasons_4)+
  #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
  # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  #stat_poly_eq(color = 'black')+
  stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
  guides(size = 'none')+
  labs(y = 'Mean growth rates ASV-based \nOM60(NOR5) clade', 
       x = 'Mean growth rates abundance-based NOR5/OM60', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/3), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# nor5_gr %>%
#   ggplot(aes(mean_dapis_nor5, mean_motus))+
#   geom_point(aes(shape = treatment, color = season, size = 2))+
#   scale_fill_manual(values = palette_seasons_4)+
#   #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
#   # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
#   scale_color_manual(values = palette_seasons_4)+
#   stat_poly_line(color = 'black')+
#   #stat_poly_eq(color = 'black')+
#   stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
#   guides(size = 'none')+
#   labs(y = 'Mean growth rates mOTUs-based Roseobacter/nor5bacter',
#        x = 'Mean growth rates abundance-based Roseobacter/nor5bacter', color = 'Season', shape = 'Treatment')+
#   theme_bw()+
#   theme(aspect.ratio = (4/4), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

##Bacteroidetes CF319a--------
mOTUS_Bacteroidetes <-  motus_env_l_w_fold %>%
  # motus_env_l_w_h_growth_1perc %>% (només els que tenen més d'un 1%)
  dplyr::filter(gr_day > 0 &
                  consensus_taxonomy_mOTUs_num != 'unassigned') %>% 
  dplyr::filter(str_detect(consensus_taxonomy_mOTUs_num, "Bacteroidetes")) %>%
  group_by(season, treatment) %>%
  summarize(mean_motus = mean(gr_day))

asv_Bacteroidetes <-  reg_all_slopes_chosen_silva_tax %>%
  dplyr::filter(season != "Early_fall") %>%
  dplyr::filter(slope_chosen_days > 0,
                pvalue_slope_chosen < 0.05) %>%
  dplyr::filter(str_detect(family_f, "Bacteroidetes")) %>%
  group_by(season, treatment) %>%
  summarize(mean_asv = mean(slope_chosen_days))

dapis_Bacteroidetes <-  GR_dapis_OS_filt %>%
  dplyr::select(treatment, replicate, season, CFB) %>%
  group_by(treatment, season) %>%
  dplyr::summarize(mean_dapis_Bacteroidetes = mean(as.numeric(CFB))) %>%
  dplyr::filter(mean_dapis_Bacteroidetes >0)

Bacteroidetes_gr <- mOTUS_Bacteroidetes %>%
  left_join(asv_Bacteroidetes) %>%
  right_join(dapis_Bacteroidetes) 

Bacteroidetes_gr %>%
  ggplot(aes(mean_asv, mean_motus))+
  geom_point(aes(shape = treatment, color = season, size = 2))+
  scale_fill_manual(values = palette_seasons_4)+
  #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
  # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  #stat_poly_eq(color = 'black')+
  stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
  guides(size = 'none')+
  labs(x = 'Mean growth rates ASV based \n
       Bacteroidetes', y = 'Mean growth rates mOTUs based Bacteroidetes', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/3), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

Bacteroidetes_gr %>%
  ggplot(aes(mean_dapis_Bacteroidetes, mean_asv))+
  geom_point(aes(shape = treatment, color = season, size = 2))+
  scale_fill_manual(values = palette_seasons_4)+
  #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
  # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  #stat_poly_eq(color = 'black')+
  stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
  guides(size = 'none')+
  labs(y = 'Mean growth rates ASV based \nOM60(Bacteroidetes) clade', 
       x = 'Mean growth rates abundance-based Bacteroidetes', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/4), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

Bacteroidetes_gr %>%
  ggplot(aes(mean_dapis_Bacteroidetes, mean_motus))+
  geom_point(aes(shape = treatment, color = season, size = 2))+
  scale_fill_manual(values = palette_seasons_4)+
  #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
  # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  #stat_poly_eq(color = 'black')+
  stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
  guides(size = 'none')+
  labs(y = 'Mean growth rates mOTUs-based Bacteroidetes',
       x = 'Mean growth rates abundance-based Bacteroidetes', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/4), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

##Eubacteria----
mOTUS_Eubacteria <-  motus_env_l_w_fold %>%
  # motus_env_l_w_h_growth_1perc %>% (només els que tenen més d'un 1%)
  dplyr::filter(gr_day > 0 &
                  consensus_taxonomy_mOTUs_num != 'unassigned') %>% 
  dplyr::filter(!str_detect(consensus_taxonomy_mOTUs_num, "archaea") |
                  !str_detect(consensus_taxonomy_mOTUs_num, "Archaea")) %>%
  group_by(season, treatment) %>%
  summarize(mean_motus = mean(gr_day))

asv_eubacateria <-  reg_all_slopes_chosen_silva_tax %>%
  dplyr::filter(season != "Early_fall") %>%
  dplyr::filter(slope_chosen_days > 0,
                pvalue_slope_chosen < 0.05) %>%
  dplyr::filter(!str_detect(domain_f, "Archaea")) %>%
  group_by(season, treatment) %>%
  summarize(mean_asv = mean(slope_chosen_days))

dapis_eubacteria <-  GR_dapis_OS_filt %>%
  dplyr::select(treatment, replicate, season, EUB) %>%
  group_by(treatment, season) %>%
  dplyr::summarize(mean_dapis_Bacteroidetes = mean(as.numeric(EUB))) %>%
  dplyr::filter(mean_dapis_Bacteroidetes >0)

eub_gr <- mOTUS_Eubacteria %>%
  left_join(asv_eubacateria) %>%
  right_join(dapis_eubacteria) 

eub_gr %>%
  ggplot(aes(mean_asv, mean_motus))+
  geom_point(aes(shape = treatment, color = season, size = 2))+
  scale_fill_manual(values = palette_seasons_4)+
  #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
  # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  #stat_poly_eq(color = 'black')+
  stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
  guides(size = 'none')+
  labs(x = 'Mean growth rates ASV based \n
       Bacteroidetes', y = 'Mean growth rates mOTUs based Bacteroidetes', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/3), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

eub_asv_gr_filt <- eub_gr %>%
  filter(season != 'Summer') %>%
  ggplot(aes(mean_dapis_Bacteroidetes, mean_asv))+
  geom_point(aes(shape = treatment, color = season, size = 2))+
  scale_fill_manual(values = palette_seasons_4)+
  #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
  # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  #stat_poly_eq(color = 'black')+
  stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
  guides(size = 'none')+
  labs(y = 'Mean growth rates ASV-based \nEubacteria', 
       x = 'Mean growth rates abundance-based Eubacteria', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/4), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

eub_gr %>%
  ggplot(aes(mean_dapis_Bacteroidetes, mean_motus))+
  geom_point(aes(shape = treatment, color = season, size = 2))+
  scale_fill_manual(values = palette_seasons_4)+
  #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
  # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  #stat_poly_eq(color = 'black')+
  stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
  guides(size = 'none')+
  labs(y = 'Mean growth rates mOTUs-based Bacteroidetes',
       x = 'Mean growth rates abundance-based Bacteroidetes', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/4), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

##COMPARING ONLY ASV-BASED GROWTH RATES VS ABUNDANCE-BASED GR
###ALt1413 filter ASVs that match with this probe-----
# mOTUS_alteromonas <-  motus_env_l_w_fold_filt %>%
#   # motus_env_l_w_h_growth_1perc %>% (només els que tenen més d'un 1%)
#   dplyr::filter(gr_day > 0 &
#                   consensus_taxonomy_mOTUs_num != 'unassigned') %>% 
#   dplyr::filter(str_detect(consensus_taxonomy_mOTUs_num, "Alteromonas") | 
#                   str_detect(consensus_taxonomy_mOTUs_num, "Colwellia")
#   ) %>%
#   group_by(season, treatment) %>%
#   summarize(mean_motus = mean(gr_day))

##check how many genus of Enteromonadaceae we have, to use only growth rates from the ASVs that match with ALT1413 FISH probe
reg_all_slopes_chosen_silva_tax %>%
  dplyr::filter(season != "Early_fall") %>%
  dplyr::filter(slope_chosen_days > 0,
                pvalue_slope_chosen < 0.05) %>%
  dplyr::filter(family_f == 'Alteromonadaceae') %$%
  genus_f %>%
  unique()

asv_alteromonadaceae <-  reg_all_slopes_chosen_silva_tax %>%
  dplyr::filter(season != "Early_fall") %>%
  dplyr::filter(slope_chosen_days > 0,
                pvalue_slope_chosen < 0.05) %>%
  dplyr::filter(str_detect(family_f, "Alteromonadaceae") #|
                  str_detect(family_f, "Colwelliaceae")
                  ) %>%
  group_by(season, treatment) %>%
  summarize(mean_asv = mean(slope_chosen_days))

dapis_alt <-  GR_dapis_OS_filt %>%
  select(treatment, replicate, season, ALT) %>%
  group_by(treatment, season) %>%
  summarize(mean_dapis_alt = mean(ALT)) %>%
  filter(mean_dapis_alt >0)

gr_alt <- asv_alteromonadaceae %>%
  left_join(dapis_alt) 

alt1413_asv_gr_filt <- gr_alt %>%
  filter(season != 'Summer') %>%
  ggplot(aes(mean_asv, mean_dapis_alt))+
  geom_point(aes(shape = treatment, color = season, size = 2))+
  scale_fill_manual(values = palette_seasons_4)+
  #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
  # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  #stat_poly_eq(color = 'black')+
  stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
  guides(size = 'none')+
  labs(x = 'Mean growth rates ASV-based \nAlteromonadaceae', 
       y = 'Mean growth rates abundance-based Alteromonadaceae', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/3), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# gr_alt %>%
#   ggplot(aes(mean_dapis_alt, mean_asv))+
#   geom_point(aes(shape = treatment, color = season, size = 2))+
#   scale_fill_manual(values = palette_seasons_4)+
#   #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
#   # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
#   scale_color_manual(values = palette_seasons_4)+
#   stat_poly_line(color = 'black')+
#   #stat_poly_eq(color = 'black')+
#   stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
#   guides(size = 'none')+
#   labs(y = 'Mean growth rates ASV based\nAlteromonas & Colwellia', 
#        x = 'Mean growth rates abundance-based\nAlteromonas & Colwellia', color = 'Season', shape = 'Treatment')+
#   theme_bw()+
#   theme(aspect.ratio = (4/3), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# gr_alt %>%
#   ggplot(aes(mean_dapis_alt, mean_motus))+
#   geom_point(aes(shape = treatment, color = season, size = 2))+
#   scale_fill_manual(values = palette_seasons_4)+
#   #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
#   # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
#   scale_color_manual(values = palette_seasons_4)+
#   stat_poly_line(color = 'black')+
#   #stat_poly_eq(color = 'black')+
#   stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
#   guides(size = 'none')+
#   labs(y = 'Mean growth rates mOTUs-based \n
#        Alteromonas & Colwellia', x = 'Mean growth rates abundance-based Alteromonas & Colwellia', color = 'Season', shape = 'Treatment')+
#   theme_bw()+
#   theme(aspect.ratio = (4/3), panel.grid.major = element_blank(), panel.grid.minor = element_blank())


##Rodobacteraceae (82.8%) match with probe ros537----
mOTUS_rhodo <-  motus_env_l_w_fold %>%
  # motus_env_l_w_h_growth_1perc %>% (només els que tenen més d'un 1%)
  dplyr::filter(gr_day > 0 &
                  consensus_taxonomy_mOTUs_num != 'unassigned') %>% 
  dplyr::filter(str_detect(consensus_taxonomy_mOTUs_num, "Rhodobacter") |
                  str_detect(consensus_taxonomy_mOTUs_num, "Roseobacteraceae")) %>%
  group_by(season, treatment) %>%
  summarize(mean_motus = mean(gr_day))

asv_rhodo <-  reg_all_slopes_chosen_silva_tax %>%
  dplyr::filter(season != "Early_fall") %>%
  dplyr::filter(slope_chosen_days > 0,
                pvalue_slope_chosen < 0.05) %>%
  dplyr::filter(str_detect(family_f, "Rhodobacteraceae")) %>% 
  group_by(season, treatment) %>%
  summarize(mean_asv = mean(slope_chosen_days))

dapis_rhodo <-  GR_dapis_OS_filt %>%
  dplyr::select(treatment, replicate, season, ROSEO) %>%
  group_by(treatment, season) %>%
  dplyr::summarize(mean_dapis_rhodo = mean(as.numeric(ROSEO))) %>%
  dplyr::filter(mean_dapis_rhodo >0)

rhodo_gr <- mOTUS_rhodo %>%
  right_join(asv_rhodo) %>%
  right_join(dapis_rhodo) 

rhodo_gr %>%
  ggplot(aes(mean_asv, mean_motus))+
  geom_point(aes(shape = treatment, color = season, size = 2))+
  scale_fill_manual(values = palette_seasons_4)+
  #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
  # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  #stat_poly_eq(color = 'black')+
  stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
  guides(size = 'none')+
  labs(x = 'Mean growth rates ASV based \n
       Rhodobacteraceae', y = 'Mean growth rates mOTUs-based Rhodobacteraceae', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/3), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ros537_asv_gr_filt <- rhodo_gr %>%
  filter(season != 'Summer') %>%
  ggplot(aes(mean_dapis_rhodo, mean_asv))+
  geom_point(aes(shape = treatment, color = season, size = 2))+
  scale_fill_manual(values = palette_seasons_4)+
  #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
  # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  #stat_poly_eq(color = 'black')+
  stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
  guides(size = 'none')+
  labs(y = 'Mean growth rates ASV-based\nRhodobacteraceae', 
       x = 'Mean growth rates abundance-based\nRhodobacteraceae', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/3), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

rhodo_gr %>%
  ggplot(aes(mean_dapis_rhodo, mean_motus))+
  geom_point(aes(shape = treatment, color = season, size = 2))+
  scale_fill_manual(values = palette_seasons_4)+
  #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
  # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  #stat_poly_eq(color = 'black')+
  stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
  guides(size = 'none')+
  labs(y = 'Mean growth rates mOTUs-based Rhodobacteraceae',
       x = 'Mean growth rates abundance-based Rhodobacteraceae', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/4), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

##Pannel construction----
card_fish_asv_gr  <- multi_panel_figure(columns = 2, rows = 4, width = 250, height = 400, 
                                              row_spacing = 2, unit = 'mm',
                                              panel_label_type = 'none')

card_fish_asv_gr  %<>%
  fill_panel(eub_asv_gr, row = 1, col = 1) %<>%
  fill_panel(gamma42a_asv_gr, row = 1, col = 2) %<>%
  fill_panel(cf319a_asv_gr, row = 2, col = 1) %<>%
  fill_panel(sar11_clade_asv_gr, row = 2, col = 2) %<>%
  fill_panel(nor5_asv_gr, row = 3, col = 1) %<>%
  fill_panel(alt1413_asv_gr, row = 3, col = 2) %<>%
  fill_panel(ros537_asv_gr, row = 4, col = 1) 

ggsave('card_fish_asv_gr.pdf', card_fish_asv_gr, 
       path = "~/Documentos/Doctorat/REMEI/results/figures/corrected_asv_num/",
       width = 250,
       height = 400,
       units = 'mm')

##seems like summer is affecting the correlations why?
card_fish_asv_gr_filt  <- multi_panel_figure(columns = 2, rows = 4, width = 250, height = 400, 
                                        row_spacing = 2, unit = 'mm',
                                        panel_label_type = 'none')

card_fish_asv_gr_filt  %<>%
  fill_panel(eub_asv_gr_filt, row = 1, col = 1) %<>%
  fill_panel(gamma42a_asv_gr_filt, row = 1, col = 2) %<>%
  fill_panel(cf319a_asv_gr_filt, row = 2, col = 1) %<>%
  fill_panel(sar11_clade_asv_gr_filt, row = 2, col = 2) %<>%
  fill_panel(nor5_asv_gr_filt, row = 3, col = 1) %<>%
  fill_panel(alt1413_asv_gr_filt, row = 3, col = 2) %<>%
  fill_panel(ros537_asv_gr_filt, row = 4, col = 1) 

ggsave('card_fish_asv_gr_filt.pdf', card_fish_asv_gr_filt, 
       path = "~/Documentos/Doctorat/REMEI/results/figures/corrected_asv_num/",
       width = 250,
       height = 400,
       units = 'mm')


###REPEAT WITH 16S rRNA copy number NORMALIZED--------------
#upload data -----
reg_all_slopes_chosen_silva_tax_cn_nor <- read.csv("data/intermediate_files/reg_all_slopes_chosen_silva_tax_cn_nor.csv", sep=",")

reg_all_slopes_chosen_silva_tax_cn_nor_filt <- reg_all_slopes_chosen_silva_tax_cn_nor %>%
  dplyr::filter(pvalue_slope_chosen < 0.05 &
                  slope_chosen > 0) %>%
  mutate(slope_chosen_days = slope_chosen*24)

##alteromonadaceae-----
asv_alteromonadaceae <-  reg_all_slopes_chosen_silva_tax_cn_nor_filt %>%
  dplyr::filter(str_detect(family, "Alteromonadaceae") |
                str_detect(family, "Colwelliaceae")
  ) %>%
  group_by(season, treatment) %>%
  summarize(mean_asv = mean(slope_chosen_days))

dapis_alt <-  GR_dapis_OS_filt %>%
  select(treatment, replicate, season, ALT) %>%
  group_by(treatment, season) %>%
  summarize(mean_dapis_alt = mean(ALT)) %>%
  filter(mean_dapis_alt >0)

gr_alt <- asv_alteromonadaceae %>%
  left_join(dapis_alt) 

alt1413_asv_gr_filt <- gr_alt %>%
  filter(season != 'Summer') %>%
  ggplot(aes(mean_asv, mean_dapis_alt))+
  geom_point(aes(shape = treatment, color = season, size = 2))+
  scale_fill_manual(values = palette_seasons_4)+
  #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
  # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  #stat_poly_eq(color = 'black')+
  stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
  guides(size = 'none')+
  labs(x = 'Mean growth rates ASV-based \nAlteromonadaceae', 
       y = 'Mean growth rates abundance-based Alteromonadaceae', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/3), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

##Rodobacteraceae (82.8%) match with probe ros537-----
mOTUS_rhodo <-  motus_env_l_w_fold %>%
  # motus_env_l_w_h_growth_1perc %>% (només els que tenen més d'un 1%)
  dplyr::filter(gr_day > 0 &
                  consensus_taxonomy_mOTUs_num != 'unassigned') %>% 
  dplyr::filter(str_detect(consensus_taxonomy_mOTUs_num, "Rhodobacter") |
                  str_detect(consensus_taxonomy_mOTUs_num, "Roseobacteraceae")) %>%
  group_by(season, treatment) %>%
  summarize(mean_motus = mean(gr_day))

asv_rhodo <-  reg_all_slopes_chosen_silva_tax_cn_nor_filt %>%
  dplyr::filter(str_detect(family, "Rhodobacteraceae")) %>% 
  group_by(season, treatment) %>%
  summarize(mean_asv = mean(slope_chosen_days))

dapis_rhodo <-  GR_dapis_OS_filt %>%
  dplyr::select(treatment, replicate, season, ROSEO) %>%
  group_by(treatment, season) %>%
  dplyr::summarize(mean_dapis_rhodo = mean(as.numeric(ROSEO))) %>%
  dplyr::filter(mean_dapis_rhodo >0)

rhodo_gr <- mOTUS_rhodo %>%
  right_join(asv_rhodo) %>%
  right_join(dapis_rhodo) 

# rhodo_gr %>%
#   ggplot(aes(mean_asv, mean_motus))+
#   geom_point(aes(shape = treatment, color = season, size = 2))+
#   scale_fill_manual(values = palette_seasons_4)+
#   #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
#   # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
#   scale_color_manual(values = palette_seasons_4)+
#   stat_poly_line(color = 'black')+
#   #stat_poly_eq(color = 'black')+
#   stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
#   guides(size = 'none')+
#   labs(x = 'Mean growth rates ASV based \n
#        Rhodobacteraceae', y = 'Mean growth rates mOTUs-based Rhodobacteraceae', color = 'Season', shape = 'Treatment')+
#   theme_bw()+
#   theme(aspect.ratio = (4/3), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ros537_asv_gr_filt <- rhodo_gr %>%
  filter(season != 'Summer') %>%
  ggplot(aes(mean_dapis_rhodo, mean_asv))+
  geom_point(aes(shape = treatment, color = season, size = 2))+
  scale_fill_manual(values = palette_seasons_4)+
  #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
  # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  #stat_poly_eq(color = 'black')+
  stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
  guides(size = 'none')+
  labs(y = 'Mean growth rates ASV-based\nRhodobacteraceae', 
       x = 'Mean growth rates abundance-based\nRhodobacteraceae', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/3), panel.grid.major = element_blank(), panel.grid.minor = element_blank())


##Eubacteria-----
asv_eubacateria <-  reg_all_slopes_chosen_silva_tax_cn_nor_filt %>%
  dplyr::filter(!str_detect(domain, "Archaea")) %>%
  group_by(season, treatment) %>%
  summarize(mean_asv = mean(slope_chosen_days))

dapis_eubacteria <-  GR_dapis_OS_filt %>%
  dplyr::select(treatment, replicate, season, EUB) %>%
  group_by(treatment, season) %>%
  dplyr::summarize(mean_dapis_Bacteroidetes = mean(as.numeric(EUB))) %>%
  dplyr::filter(mean_dapis_Bacteroidetes >0)

eub_gr <- asv_eubacateria%>%
  right_join(dapis_eubacteria) 

eub_asv_gr_filt <- eub_gr %>%
  filter(season != 'Summer') %>%
  ggplot(aes(mean_dapis_Bacteroidetes, mean_asv))+
  geom_point(aes(shape = treatment, color = season, size = 2))+
  scale_fill_manual(values = palette_seasons_4)+
  #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
  # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  #stat_poly_eq(color = 'black')+
  stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
  guides(size = 'none')+
  labs(y = 'Mean growth rates ASV-based \nEubacteria', 
       x = 'Mean growth rates abundance-based Eubacteria', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/4), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#Gammaproteobacteria------
asv_gamma <-  reg_all_slopes_chosen_silva_tax_cn_nor_filt %>%
  dplyr::filter(str_detect(class, "Gammaproteobacteria")) %>%
  group_by(season, treatment) %>%
  summarize(mean_asv = mean(slope_chosen_days))

dapis_gamma <-  GR_dapis_OS_filt %>%
  dplyr::select(treatment, replicate, season, GAMMA) %>%
  group_by(treatment, season) %>%
  dplyr::summarize(mean_dapis_gamma = mean(as.numeric(GAMMA))) %>%
  dplyr::filter(mean_dapis_gamma >0)

gamma_gr <- asv_gamma %>%
  left_join(dapis_gamma) 

gamma42a_asv_gr_filt <- gamma_gr %>%
  filter(season != 'Summer') %>%
  ggplot(aes(mean_dapis_gamma, mean_asv))+
  geom_point(aes(shape = treatment, color = season, size = 2))+
  scale_fill_manual(values = palette_seasons_4)+
  #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
  # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  #stat_poly_eq(color = 'black')+
  stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
  guides(size = 'none')+
  labs(y = 'Mean growth rates ASV-based \nGammaproteobacteria', 
       x = 'Mean growth rates abundance-based Gammaproteobacteria', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/3), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

##CFB bacteroidetes----
asv_cfb <-  reg_all_slopes_chosen_silva_tax_cn_nor_filt %>%
  dplyr::filter(str_detect(phylum, "Bacteroidota")) %>%
  group_by(season, treatment) %>%
  summarize(mean_asv = mean(slope_chosen_days))

dapis_cfb <-  GR_dapis_OS_filt %>%
  dplyr::select(treatment, replicate, season, CFB) %>%
  group_by(treatment, season) %>%
  dplyr::summarize(mean_dapis_fcb = mean(as.numeric(CFB))) %>%
  dplyr::filter(mean_dapis_fcb >0)

cfb_gr <- asv_cfb %>%
  left_join(dapis_cfb) 

# fcb_gr %>%
#   ggplot(aes(mean_asv, mean_motus))+
#   geom_point(aes(shape = treatment, color = season, size = 2))+
#   scale_fill_manual(values = palette_seasons_4)+
#   #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
#   # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
#   scale_color_manual(values = palette_seasons_4)+
#   stat_poly_line(color = 'black')+
#   #stat_poly_eq(color = 'black')+
#   stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
#   guides(size = 'none')+
#   labs(x = 'Mean growth rates ASV based \n
#        fcb clade', y = 'Mean growth rates mOTUs based fcb clade', color = 'Season', shape = 'Treatment')+
#   theme_bw()+
#   theme(aspect.ratio = (4/3), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

cf319a_asv_gr_filt <- cfb_gr %>%
  filter(season != 'Summer') %>%
  ggplot(aes(mean_dapis_fcb, mean_asv))+
  geom_point(aes(shape = treatment, color = season, size = 2))+
  scale_fill_manual(values = palette_seasons_4)+
  #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
  # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  #stat_poly_eq(color = 'black')+
  stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
  guides(size = 'none')+
  labs(y = 'Mean growth rates ASV based Bacteroidota', x = 'Mean growth rates abundance-based Bacteroidota', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/3), panel.grid.major = element_blank(), panel.grid.minor = element_blank())


##SAR11 clade----
asv_sar11 <-  reg_all_slopes_chosen_silva_tax_cn_nor_filt %>%
  dplyr::filter(str_detect(order, "SAR11 clade")) %>%
  group_by(season, treatment) %>%
  summarize(mean_asv = mean(slope_chosen_days))

dapis_sar11 <-  GR_dapis_OS_filt %>%
  dplyr::select(treatment, replicate, season, SAR11) %>%
  group_by(treatment, season) %>%
  dplyr::summarize(mean_dapis_sar11 = mean(as.numeric(SAR11))) %>%
  dplyr::filter(mean_dapis_sar11 >0)

sar11_gr <- #mOTUS_sar11 %>%
  asv_sar11 %>%
  left_join(dapis_sar11) 

# sar11_gr %>%
#   ggplot(aes(mean_asv, mean_motus))+
#   geom_point(aes(shape = treatment, color = season, size = 2))+
#   scale_fill_manual(values = palette_seasons_4)+
#   #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
#   # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
#   scale_color_manual(values = palette_seasons_4)+
#   stat_poly_line(color = 'black')+
#   #stat_poly_eq(color = 'black')+
#   stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
#   guides(size = 'none')+
#   labs(x = 'Mean growth rates ASV based \n
#        SAR11 clade', y = 'Mean growth rates mOTUs based SAR11 clade', color = 'Season', shape = 'Treatment')+
#   theme_bw()+
#   theme(aspect.ratio = (4/3), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

sar11_clade_asv_gr_filt <- sar11_gr %>%
  filter(season != 'Summer') %>%
  ggplot(aes(mean_dapis_sar11, mean_asv))+
  geom_point(aes(shape = treatment, color = season, size = 2))+
  scale_fill_manual(values = palette_seasons_4)+
  #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
  # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  #stat_poly_eq(color = 'black')+
  stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
  guides(size = 'none')+
  labs(y = 'Mean growth rates ASV-based \nSAR11 clade', 
       x = 'Mean growth rates abundance-based SAR11 clade', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/4), panel.grid.major = element_blank(), panel.grid.minor = element_blank())


##OM60----
asv_nor5 <-  reg_all_slopes_chosen_silva_tax_cn_nor_filt %>%
  dplyr::filter(str_detect(genus, "OM60")) %>%
  group_by(season, treatment) %>%
  summarize(mean_asv = mean(slope_chosen_days))

dapis_nor5 <-  GR_dapis_OS_filt %>%
  dplyr::select(treatment, replicate, season, NOR5) %>%
  group_by(treatment, season) %>%
  dplyr::summarize(mean_dapis_nor5 = mean(as.numeric(NOR5))) %>%
  dplyr::filter(mean_dapis_nor5 >0)

nor5_gr <-  asv_nor5%>%
  right_join(dapis_nor5) 

# nor5_gr %>%
#   ggplot(aes(mean_asv, mean_motus))+
#   geom_point(aes(shape = treatment, color = season, size = 2))+
#   scale_fill_manual(values = palette_seasons_4)+
#   #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
#   # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
#   scale_color_manual(values = palette_seasons_4)+
#   stat_poly_line(color = 'black')+
#   #stat_poly_eq(color = 'black')+
#   stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
#   guides(size = 'none')+
#   labs(x = 'Mean growth rates ASV based \n
#        Roseobacter/nor5bacter', y = 'Mean growth rates mOTUs based Roseobacter/nor5bacter', color = 'Season', shape = 'Treatment')+
#   theme_bw()+
#   theme(aspect.ratio = (4/3), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

nor5_asv_gr_filt <- nor5_gr %>%
  filter(season != 'Summer') %>%
  ggplot(aes(mean_dapis_nor5, mean_asv))+
  geom_point(aes(shape = treatment, color = season, size = 2))+
  scale_fill_manual(values = palette_seasons_4)+
  #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
  # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  #stat_poly_eq(color = 'black')+
  stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
  guides(size = 'none')+
  labs(y = 'Mean growth rates ASV-based \nOM60(NOR5) clade', 
       x = 'Mean growth rates abundance-based NOR5/OM60', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/3), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

##pannel construction-----
card_fish_asv_gr  <- multi_panel_figure(columns = 2, rows = 4, width = 250, height = 400, 
                                        row_spacing = 2, unit = 'mm',
                                        panel_label_type = 'none')

card_fish_asv_gr  %<>%
  fill_panel(eub_asv_gr_filt, row = 1, col = 1) %<>%
  fill_panel(gamma42a_asv_gr_filt, row = 1, col = 2) %<>%
  fill_panel(cf319a_asv_gr_filt, row = 2, col = 1) %<>%
  fill_panel(sar11_clade_asv_gr_filt, row = 2, col = 2) %<>%
  fill_panel(nor5_asv_gr_filt, row = 3, col = 1) %<>%
  fill_panel(alt1413_asv_gr_filt, row = 3, col = 2) %<>%
  fill_panel(ros537_asv_gr_filt, row = 4, col = 1) %<>%
  fill_panel(alt_colwelia, row = 4, col = 2) 

ggsave('card_fish_asv_gr_cn_nor.pdf', card_fish_asv_gr, 
       path = "~/Documentos/Doctorat/REMEI/results/figures/corrected_asv_num/",
       width = 250,
       height = 400,
       units = 'mm')

ggsave('card_fish_asv_gr_cn_nor_no_summer.pdf', card_fish_asv_gr, 
       path = "~/Documentos/Doctorat/REMEI/results/figures/corrected_asv_num/",
       width = 250,
       height = 400,
       units = 'mm')

##Al1413 just Alteromonas & Colwellia----
asv_alteromonas <-  reg_all_slopes_chosen_silva_tax_cn_nor_filt %>%
  dplyr::filter(str_detect(genus, "Alteromonas") | 
                  str_detect(genus,'Colwellia') |
                  str_detect(genus, 'Glaciecola')) %>%
  group_by(season, treatment) %>%
  summarize(mean_asv = mean(slope_chosen_days))

dapis_alt <-  GR_dapis_OS_filt %>%
  select(treatment, replicate, season, ALT) %>%
  group_by(treatment, season) %>%
  summarize(mean_dapis_alt = mean(ALT)) %>%
  filter(mean_dapis_alt >0)

gr_alt <-  asv_alteromonas %>%
  left_join(dapis_alt) 

alt_colwelia <- gr_alt %>%
  filter(season != 'Summer') %>%
  ggplot(aes(mean_dapis_alt, mean_asv))+
  geom_point(aes(shape = treatment, color = season, size = 2))+
  scale_fill_manual(values = palette_seasons_4)+
  #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
  # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  #stat_poly_eq(color = 'black')+
  stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
  guides(size = 'none')+
  labs(y = 'Mean growth rates ASV based\nAlteromonas & Colwellia', 
       x = 'Mean growth rates abundance-based\nAlteromonas & Colwellia', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/3), panel.grid.major = element_blank(), panel.grid.minor = element_blank())


##Another option: select only those ASVs that are dominating the community, since they are the one's that are dominating CARD_FISH (CORRECT)------
##probes

##first, calculate mean abundance in each SEASON and TREATMENT per ASV
###non normalized dataset
rem_relabun_melt <-  read.table("data/rem_relabun_melt.txt", sep="\t", header = TRUE) %>%
  filter(season != "Early_fall")

rem_relabun_melt %$%
  asv_num %>%
  unique() #4594 asv al meu dataset

rem_relabun_melt %>%
  colnames()

mean_abundance_asv <- rem_relabun_melt %>% 
  group_by(asv_num, treatment, replicate, season, domain, phylum, class, order, family, genus) %>% ##mitjana repliques
  dplyr::summarize(mean_rep = mean(Abundance))  %>%
  group_by(season, treatment, asv_num, domain, phylum, class, order, family, genus) %>%
  dplyr::summarize(mean_abund = mean(mean_rep)) 

###normalized dataset
asv_tab_cn_nor_rel_l_filt <- read.table('data/intermediate_files/asv_tab_cn_nor_rel_l_filt.txt', sep = '\t', header = TRUE) %>%
  filter(season != "Early_fall")
asv_tab_cn_nor_rel_l %>%
  colnames()

mean_abundance_asv_nor <- asv_tab_cn_nor_rel_l_filt %>% 
  group_by(asv_num, treatment, replicate, season, domain, phylum, class, order, family, genus) %>% ##mitjana repliques
  #filter(abundance != is.na(abundance)) %>%
  dplyr::summarize(mean_rep = mean(as.numeric(abundance)))  %>%
  #filter(mean_rep != is.na(mean_rep)) %>%
  group_by(season, treatment, asv_num, domain, phylum, class, order, family, genus) %>%
  dplyr::summarize(mean_abund = mean(as.numeric(mean_rep)))

##Most abundant alteromonas (16S rRNA normalized)----
x <- mean_abundance_asv_nor %>%
  dplyr::filter(genus %in% c('Alteromonas', 'Colwellia', 'Glaciecola')) %>%
  filter(mean_abund != is.na(mean_abund)) %>%
  #group_by(season, treatment) %>%
  group_by(family) %$%
  mean(mean_abund)

asv_alteromonas <- mean_abundance_asv_nor %>%
  dplyr::filter(genus %in% c('Alteromonas', 'Colwellia', 'Glaciecola')  & #
                  mean_abund > 0) %>%
  #filter(asv_num != 'asv53') %>%
  left_join(reg_all_slopes_chosen_silva_tax_cn_nor_filt, by = c('asv_num', 'season', 'treatment')) %>%
  filter(slope_chosen_days != is.na(slope_chosen_days)) %>%
  group_by(season, treatment) %>%
  summarize(mean_asv = mean(slope_chosen_days))

dapis_alt <-  GR_dapis_OS_filt %>%
  select(treatment, replicate, season, ALT) %>%
  group_by(treatment, season) %>%
  summarize(mean_dapis_alt = mean(ALT)) %>%
  filter(mean_dapis_alt >0)

gr_alt <-  asv_alteromonas %>%
  left_join(dapis_alt) 

alt1413_asv_gr_nor <- gr_alt %>%
  filter(season != 'Summer') %>%
  ggplot(aes(mean_dapis_alt, mean_asv))+
  geom_point(aes(shape = treatment, color = season, size = 0.5))+
  scale_fill_manual(values = palette_seasons_4)+
  #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
  # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  #stat_poly_eq(color = 'black')+
  stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
  guides(size = 'none')+
  labs(y = 'Mean growth rates ASV-based\nAlteromonas, Colwellia, Glaciecola', 
       x = 'Mean growth rates abundance-based\nALT1413', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/3), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = 'none', text= element_text(size = 5))

##look what is happening with AL1413 in summer, which ASVs are there?
x <- mean_abundance_asv_nor %>%
  dplyr::filter(genus %in% c('Alteromonas', 'Colwellia', 'Glaciecola')) %>%
  filter(mean_abund != is.na(mean_abund)) %>%
  #group_by(season, treatment) %>%
  group_by(family) %$%
  mean(mean_abund)

asv_alteromonas_summer <- mean_abundance_asv_nor %>%
  dplyr::filter(genus %in% c('Alteromonas', 'Colwellia', 'Glaciecola')  & #
                  mean_abund > x) %>%
  left_join(reg_all_slopes_chosen_silva_tax_cn_nor_filt, by = c('asv_num', 'season', 'treatment')) %>%
  filter(slope_chosen_days != is.na(slope_chosen_days)) %>%
  filter(season == 'Summer')

asv_alteromonas_seasons <- mean_abundance_asv_nor %>%
  dplyr::filter(genus %in% c('Alteromonas', 'Colwellia', 'Glaciecola')  & #
                  mean_abund > x) %>%
  left_join(reg_all_slopes_chosen_silva_tax_cn_nor_filt, by = c('asv_num', 'season', 'treatment')) %>%
  filter(slope_chosen_days != is.na(slope_chosen_days)) %>%
  filter(season %in% c('Winter', 'Spring', 'Fall')) %$%
  unique(asv_num)
# 
#   group_by(season, treatment) %>%
#   summarize(mean_asv = mean(slope_chosen_days))

##alteromonas, abundant ASVs without 16S rRNA normalization----
x <- mean_abundance_asv %>%
  dplyr::filter(genus %in% c('Alteromonas', 'Colwellia', 'Glaciecola')) %>%
  filter(mean_abund != is.na(mean_abund)) %>%
  #group_by(season, treatment) %>%
  group_by(family) %$%
  mean(mean_abund) 

reg_all_slopes_chosen_silva_tax_filt %>%
  colnames

asv_alteromonas_filt <- mean_abundance_asv %>%
  dplyr::filter(genus %in% c('Alteromonas', 'Colwellia', 'Glaciecola')  & #
                  mean_abund > x
                  ) %>%
  left_join(reg_all_slopes_chosen_silva_tax_filt, by = c('asv_num', 'season', 'treatment'), multiple = 'all') %>%
  filter(slope_chosen_days != is.na(slope_chosen_days)) %>%
  group_by(season, treatment) %>%
  dplyr::summarize(mean_asv_alt = mean(slope_chosen_days))

dapis_alt <-  GR_dapis_OS_filt %>%
  select(treatment, replicate, season, ALT) %>%
  group_by(treatment, season) %>%
  dplyr::summarize(mean_dapis_alt = mean(ALT)) %>%
  filter(mean_dapis_alt >0)

gr_alt_nn <-  asv_alteromonas_filt %>%
  left_join(dapis_alt) 

alt1413_asv_gr <- gr_alt %>%
  #filter(season != 'Summer') %>%
  ggplot(aes(mean_dapis_alt, mean_asv_alt))+
  geom_point(aes(shape = treatment, color = season, size = 0.5))+
  scale_fill_manual(values = palette_seasons_4)+
  #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
  # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  #stat_poly_eq(color = 'black')+
  stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
  guides(size = 'none')+
  labs(y = 'Mean growth rates ASV-based\nAlteromonas, Colwellia, Glaciecola', 
       x = 'Mean growth rates abundance-based\nALT1413', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/3), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = 'none', , text= element_text(size = 5))

##Most Gamma abundant 16S rRNA copy number normalized-----
x <- mean_abundance_asv_nor %>%
  dplyr::filter(class == 'Gammaproteobacteria') %>%
  filter(mean_abund != is.na(mean_abund)) %>%
  #group_by(season, treatment) %>%
  group_by(family) %$%
  mean(mean_abund)

asv_gamma <-  mean_abundance_asv_nor %>%
  dplyr::filter(class == 'Gammaproteobacteria' & #
                  mean_abund > x) %>%
  left_join(reg_all_slopes_chosen_silva_tax_cn_nor_filt, by = c('asv_num', 'season', 'treatment')) %>%
  filter(slope_chosen_days != is.na(slope_chosen_days)) %>%
  group_by(season, treatment) %>%
  summarize(mean_asv = mean(slope_chosen_days))

dapis_gamma <-  GR_dapis_OS_filt %>%
  dplyr::select(treatment, replicate, season, GAMMA) %>%
  group_by(treatment, season) %>%
  dplyr::summarize(mean_dapis_gamma = mean(as.numeric(GAMMA))) %>%
  dplyr::filter(mean_dapis_gamma >0)

gamma_gr <- asv_gamma %>%
  left_join(dapis_gamma) 

gamma42a_asv_gr_nor <- gamma_gr %>%
  #filter(season != 'Summer') %>%
  ggplot(aes(mean_dapis_gamma, mean_asv))+
  geom_point(aes(shape = treatment, color = season, size = 0.5))+
  scale_fill_manual(values = palette_seasons_4)+
  #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
  # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  #stat_poly_eq(color = 'black')+
  stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
  guides(size = 'none')+
  labs(y = 'Mean growth rates ASV-based \nGammaproteobacteria', 
       x = 'Mean growth rates abundance-based\nGAM42a', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/3), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = 'none', text= element_text(size = 5))

##Most Gamma abundant 16sRNA no normalized -----
x <- mean_abundance_asv %>%
  dplyr::filter(class == 'Gammaproteobacteria') %>%
  filter(mean_abund != is.na(mean_abund)) %>%
  #group_by(season, treatment) %>%
  group_by(family) %$%
  mean(mean_abund)

asv_gamma <-  mean_abundance_asv %>%
  dplyr::filter(class == 'Gammaproteobacteria' & #
                  mean_abund > x) %>%
  left_join(reg_all_slopes_chosen_silva_tax_filt, by = c('asv_num', 'season', 'treatment'), multiple = 'all') %>%
  filter(slope_chosen_days != is.na(slope_chosen_days)) %>%
  group_by(season, treatment) %>%
  summarize(mean_asv_gam = mean(slope_chosen_days))

dapis_gamma <-  GR_dapis_OS_filt %>%
  dplyr::select(treatment, replicate, season, GAMMA) %>%
  group_by(treatment, season) %>%
  dplyr::summarize(mean_dapis_gamma = mean(as.numeric(GAMMA))) %>%
  dplyr::filter(mean_dapis_gamma >0)

gamma_gr_nn <- asv_gamma %>%
  left_join(dapis_gamma) 

gamma42a_asv_gr <- gamma_gr %>%
  #filter(season != 'Summer') %>%
  ggplot(aes(mean_dapis_gamma, mean_asv_gam))+
  geom_point(aes(shape = treatment, color = season, size = 0.5))+
  scale_fill_manual(values = palette_seasons_4)+
  #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
  # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  #stat_poly_eq(color = 'black')+
  stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
  guides(size = 'none')+
  labs(y = 'Mean growth rates ASV-based \nGammaproteobacteria', 
       x = 'Mean growth rates abundance-based\nGAM42a', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/3), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = 'none', text= element_text(size = 5))

##Most Rhodobacteraceaae (16S rRNA normalized)----
##filter by those ASVs that are > than the mean realative abundance of the group we are analyzing
x <- mean_abundance_asv_nor %>%
  dplyr::filter(family %in% c('Rhodobacteraceae')) %>%
  filter(mean_abund != is.na(mean_abund)) %>%
  group_by(season, treatment) %>%
  group_by(family) %$%
  mean(mean_abund)

asv_rhodo <- mean_abundance_asv_nor %>%
  dplyr::filter(family %in% c('Rhodobacteraceae')  & #
                  mean_abund > x) %>%
  left_join(reg_all_slopes_chosen_silva_tax_cn_nor_filt, by = c('asv_num', 'season', 'treatment')) %>%
  filter(slope_chosen_days != is.na(slope_chosen_days)) %>%
  #filter(asv_num != 'asv114') %>%
  group_by(season, treatment) %>%
  summarize(mean_asv_rhodo = mean(slope_chosen_days))

dapis_rhodo <-  GR_dapis_OS_filt %>%
  dplyr::select(treatment, replicate, season, ROSEO) %>%
  group_by(treatment, season) %>%
  dplyr::summarize(mean_dapis_rhodo = mean(as.numeric(ROSEO))) %>%
  dplyr::filter(mean_dapis_rhodo >0)

rhodo_gr <-  asv_rhodo %>%
  right_join(dapis_rhodo) 

ros537_asv_gr_cn <- rhodo_gr %>%
  #filter(season != 'Summer') %>%
  ggplot(aes(mean_dapis_rhodo, mean_asv_rhodo))+
  geom_point(aes(shape = treatment, color = season, size = 0.5))+
  scale_fill_manual(values = palette_seasons_4)+
  #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
  # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  #stat_poly_eq(color = 'black')+
  stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
  guides(size = 'none')+
  labs(y = 'Mean growth rates ASV-based\nRhodobacteraceae', 
       x = 'Mean growth rates abundance-based\nROS537', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/3), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = 'none', text= element_text(size = 5))

##SUMMER, different tendency----
asv_rhodo_seasons <- mean_abundance_asv_nor %>%
  dplyr::filter(family %in% c('Rhodobacteraceae')  & #
                  mean_abund > x) %>%
  left_join(reg_all_slopes_chosen_silva_tax_cn_nor_filt, by = c('asv_num', 'season', 'treatment')) %>%
  filter(slope_chosen_days != is.na(slope_chosen_days)) %>%
  group_by(season) %>%
  summarize(asv_num = asv_num)


##Most Rhodobacteraceae no normalized----
x <- mean_abundance_asv %>%
  dplyr::filter(family %in% c('Rhodobacteraceae')) %>%
  filter(mean_abund != is.na(mean_abund)) %>%
  #group_by(season, treatment) %>%
  group_by(family) %$%
  mean(mean_abund)

asv_rhodo <- mean_abundance_asv %>%
  dplyr::filter(family %in% c('Rhodobacteraceae')  & #
                  mean_abund > x) %>%
  left_join(reg_all_slopes_chosen_silva_tax_filt, by = c('asv_num', 'season', 'treatment')) %>%
  filter(slope_chosen_days != is.na(slope_chosen_days)) %>%
  group_by(season, treatment) %>%
  summarize(mean_asv_rhodo = mean(slope_chosen_days))

dapis_rhodo <-  GR_dapis_OS_filt %>%
  dplyr::select(treatment, replicate, season, ROSEO) %>%
  group_by(treatment, season) %>%
  dplyr::summarize(mean_dapis_rhodo = mean(as.numeric(ROSEO))) %>%
  dplyr::filter(mean_dapis_rhodo >0)

rhodo_gr_nn <-  asv_rhodo %>%
  right_join(dapis_rhodo) 

ros537_asv_gr <- rhodo_gr %>%
  #filter(season != 'Summer') %>%
  ggplot(aes(mean_dapis_rhodo, mean_asv))+
  geom_point(aes(shape = treatment, color = season, size = 0.5))+
  scale_fill_manual(values = palette_seasons_4)+
  #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
  # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  #stat_poly_eq(color = 'black')+
  stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
  guides(size = 'none')+
  labs(y = 'Mean growth rates ASV-based\nRhodobacteraceae', 
       x = 'Mean growth rates abundance-based\nROS537', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/3), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = 'none', text= element_text(size = 5))

##Most Bacteroidetess (Cytophaga, Flavobacteria, Sphingobacteria, Empedobacter, Chryseobacterium, Bergeyella) norm-----
x <- mean_abundance_asv_nor %>%
  dplyr::filter(family %in% c('Flavobacteriaceae', 'Cytophagaceae', 'NS11-12 marine group')) %>%
  filter(mean_abund != is.na(mean_abund)) %>%
  #group_by(season, treatment) %>%
  group_by(family) %$%
  mean(mean_abund)

asv_cfb <- mean_abundance_asv_nor %>%
  dplyr::filter(family %in% c('Flavobacteriaceae', 'Cytophagaceae', 'NS11-12 marine group')  & #
                  mean_abund > x) %>%
  left_join(reg_all_slopes_chosen_silva_tax_cn_nor_filt, by = c('asv_num', 'season', 'treatment')) %>%
  filter(slope_chosen_days != is.na(slope_chosen_days)) %>%
  group_by(season, treatment) %>%
  summarize(mean_asv = mean(slope_chosen_days))

dapis_cfb <-  GR_dapis_OS_filt %>%
  dplyr::select(treatment, replicate, season, CFB) %>%
  group_by(treatment, season) %>%
  dplyr::summarize(mean_dapis_fcb = mean(as.numeric(CFB))) %>%
  dplyr::filter(mean_dapis_fcb >0)

cfb_gr <- asv_cfb %>%
  left_join(dapis_cfb) 


cf319a_asv_gr_cn <- cfb_gr %>%
  #filter(season != 'Summer') %>%
  ggplot(aes(mean_dapis_fcb, mean_asv))+
  geom_point(aes(shape = treatment, color = season, size = 0.5))+
  scale_fill_manual(values = palette_seasons_4)+
  #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
  # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  #stat_poly_eq(color = 'black')+
  stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
  guides(size = 'none')+
  labs(y = 'Mean growth rates ASV-based Bacteroidota', 
       x = 'Mean growth rates abundance-based\nCF319a', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/3), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = 'none', text= element_text(size = 5))

##Most Bacteroidetess (Cytophaga, Flavobacteria, Sphingobacteria, Empedobacter, Chryseobacterium, Bergeyella) no norm-----
x <- mean_abundance_asv %>%
  dplyr::filter(family %in% c('Flavobacteriaceae', 'Cytophagaceae', 'NS11-12 marine group')) %>%
  filter(mean_abund != is.na(mean_abund)) %>%
  #group_by(season, treatment) %>%
  group_by(family) %$%
  mean(mean_abund)

asv_cfb <- mean_abundance_asv %>%
  dplyr::filter(family %in% c('Flavobacteriaceae', 'Cytophagaceae', 'NS11-12 marine group')  & #
                  mean_abund > x) %>%
  left_join(reg_all_slopes_chosen_silva_tax_filt, by = c('asv_num', 'season', 'treatment'), multiple = 'all') %>%
  filter(slope_chosen_days != is.na(slope_chosen_days)) %>%
  group_by(season, treatment) %>%
  summarize(mean_asv_cfb = mean(slope_chosen_days))

dapis_cfb <-  GR_dapis_OS_filt %>%
  dplyr::select(treatment, replicate, season, CFB) %>%
  group_by(treatment, season) %>%
  dplyr::summarize(mean_dapis_cfb = mean(as.numeric(CFB))) %>%
  dplyr::filter(mean_dapis_cfb >0)

cfb_gr_nn <- asv_cfb %>%
  left_join(dapis_cfb) 

cf319a_asv_gr <- cfb_gr %>%
  #filter(season != 'Summer') %>%
  ggplot(aes(mean_dapis_fcb, mean_asv_cfb))+
  geom_point(aes(shape = treatment, color = season, size = 0.5))+
  scale_fill_manual(values = palette_seasons_4)+
  #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
  # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  #stat_poly_eq(color = 'black')+
  stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
  guides(size = 'none')+
  labs(y = 'Mean growth rates ASV-based Flavobacteriaceae,\nCytophagaceae, NS11-12 marine group', 
       x = 'Mean growth rates abundance-based\nCF319a', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/3), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = 'none', text= element_text(size = 5))

##MOST ABUNDANT SAR11 clade (normalized)-----
x <- mean_abundance_asv_nor %>%
  dplyr::filter(order %in% c('SAR11 clade') &
                  family %in% c('Clade I', 'Clade II', NA)) %>%
  filter(mean_abund != is.na(mean_abund)) %>%
  #group_by(season, treatment) %>%
  group_by(family) %$%
  mean(mean_abund) 

asv_sar11 <-  mean_abundance_asv_nor %>%
  dplyr::filter(str_detect(order, "SAR11 clade")&
                  family %in% c('Clade I', 'Clade II', NA) &
                  mean_abund > x) %>%
  left_join(reg_all_slopes_chosen_silva_tax_cn_nor_filt, by = c('asv_num', 'season', 'treatment')) %>%
  filter(slope_chosen_days != is.na(slope_chosen_days)) %>%
  group_by(season, treatment) %>%
  summarize(mean_asv = mean(slope_chosen_days))

dapis_sar11 <-  GR_dapis_OS_filt %>%
  dplyr::select(treatment, replicate, season, SAR11) %>%
  group_by(treatment, season) %>%
  dplyr::summarize(mean_dapis_sar11 = mean(as.numeric(SAR11))) %>%
  dplyr::filter(mean_dapis_sar11 >0)

sar11_gr <- #mOTUS_sar11 %>%
  asv_sar11 %>%
  left_join(dapis_sar11) 

sar11_411_clade_asv_gr_cn <- sar11_gr %>%
  # filter(season != 'Summer') %>%
  ggplot(aes(mean_dapis_sar11, mean_asv))+
  geom_point(aes(shape = treatment, color = season, size = 0.5))+
  scale_fill_manual(values = palette_seasons_4)+
  #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
  # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  #stat_poly_eq(color = 'black')+
  stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
  guides(size = 'none')+
  labs(y = 'Mean growth rates ASV-based \nSAR11 clade I & II', 
       x = 'Mean growth rates abundance-based SAR11-441R', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/4), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = 'none' , text= element_text(size = 5))

##MOST ABUNDANT SAR11 clade (no normalized)-----
x <- mean_abundance_asv %>%
  dplyr::filter(order %in% c('SAR11 clade') &
                  family %in% c('Clade I', 'Clade II', NA)) %>%
  filter(mean_abund != is.na(mean_abund)) %>%
  #group_by(season, treatment) %>%
  group_by(family) %$%
  mean(mean_abund) 

asv_sar11 <-  mean_abundance_asv %>%
  dplyr::filter(str_detect(order, "SAR11 clade")&
                  family %in% c('Clade I', 'Clade II', NA) &
                  mean_abund > x) %>%
  left_join(reg_all_slopes_chosen_silva_tax_filt, by = c('asv_num', 'season', 'treatment')) %>%
  filter(slope_chosen_days != is.na(slope_chosen_days)) %>%
  group_by(season, treatment) %>%
  summarize(mean_asv = mean(slope_chosen_days))

dapis_sar11 <-  GR_dapis_OS_filt %>%
  dplyr::select(treatment, replicate, season, SAR11) %>%
  group_by(treatment, season) %>%
  dplyr::summarize(mean_dapis_sar11 = mean(as.numeric(SAR11))) %>%
  dplyr::filter(mean_dapis_sar11 >0)

sar11_gr <- #mOTUS_sar11 %>%
  asv_sar11 %>%
  left_join(dapis_sar11) 

sar11_411_clade_asv_gr <- sar11_gr %>%
  #filter(season != 'Summer') %>%
  ggplot(aes(mean_dapis_sar11, mean_asv))+
  geom_point(aes(shape = treatment, color = season, size = 0.5))+
  scale_fill_manual(values = palette_seasons_4)+
  #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
  # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  #stat_poly_eq(color = 'black')+
  stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
  guides(size = 'none')+
  labs(y = 'Mean growth rates ASV-based \nSAR11 clade I & II', 
       x = 'Mean growth rates abundance-based SAR11-441R', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/4), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = 'none', text= element_text(size = 5))

##MOST ABUNDANT Eubacteria (normalized)-----
x <- mean_abundance_asv_nor %>%
  dplyr::filter(!str_detect(domain, "Archaea")) %>%
  filter(mean_abund != is.na(mean_abund)) %>%
  #group_by(season, treatment) %>%
  group_by(family) %$%
  mean(mean_abund) 

asv_eubacateria <- mean_abundance_asv_nor %>%
  dplyr::filter(!str_detect(domain, "Archaea")  & #
                  mean_abund > x) %>%
  left_join(reg_all_slopes_chosen_silva_tax_cn_nor_filt, by = c('asv_num', 'season', 'treatment')) %>%
  filter(slope_chosen_days != is.na(slope_chosen_days)) %>%
  group_by(season, treatment) %>%
  summarize(mean_asv = mean(slope_chosen_days))

dapis_eubacteria <-  GR_dapis_OS_filt %>%
  dplyr::select(treatment, replicate, season, EUB) %>%
  group_by(treatment, season) %>%
  dplyr::summarize(mean_dapis_Bacteroidetes = mean(as.numeric(EUB))) %>%
  dplyr::filter(mean_dapis_Bacteroidetes >0)


eub_gr <- asv_eubacateria %>%
  right_join(dapis_eubacteria) 

eub_asv_gr_cn <- eub_gr %>%
  #filter(season != 'Summer') %>%
  ggplot(aes(mean_dapis_Bacteroidetes, mean_asv))+
  geom_point(aes(shape = treatment, color = season, size = 0.5))+
  scale_fill_manual(values = palette_seasons_4)+
  #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
  # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  #stat_poly_eq(color = 'black')+
  stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
  guides(size = 'none')+
  labs(y = 'Mean growth rates ASV-based \nEubacteria', 
       x = 'Mean growth rates abundance-based\nEubacteria-I,-II,-III', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/4), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = 'none', text= element_text(size = 5))


##MOST ABUNDANT Eubacteria (NO normalized)-----
x <- mean_abundance_asv %>%
  dplyr::filter(!str_detect(domain, "Archaea")) %>%
  filter(mean_abund != is.na(mean_abund)) %>%
  #group_by(season, treatment) %>%
  group_by(family) %$%
  mean(mean_abund) 


asv_eubacateria <- mean_abundance_asv %>%
  dplyr::filter(!str_detect(domain, "Archaea")  & #
                  mean_abund > x) %>%
  left_join(reg_all_slopes_chosen_silva_tax_filt, by = c('asv_num', 'season', 'treatment'), multiple = 'all') %>%
  filter(slope_chosen_days != is.na(slope_chosen_days)) %>%
  group_by(season, treatment) %>%
  summarize(mean_asv_eub = mean(slope_chosen_days))


dapis_eubacteria <-  GR_dapis_OS_filt %>%
  dplyr::select(treatment, replicate, season, EUB) %>%
  group_by(treatment, season) %>%
  dplyr::summarize(mean_dapis_eub = mean(as.numeric(EUB))) %>%
  dplyr::filter(mean_dapis_eub >0)

eub_gr_nn <- asv_eubacateria %>%
  right_join(dapis_eubacteria) 

eub_asv_gr <- eub_gr %>%
  #filter(season != 'Summer') %>%
  ggplot(aes(mean_dapis_Bacteroidetes, mean_asv_eub))+
  geom_point(aes(shape = treatment, color = season))+
  scale_fill_manual(values = palette_seasons_4)+
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0))+
  scale_y_continuous(limits = c(0, 4.2), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  #stat_poly_eq(color = 'black')+
  stat_cor(aes(label = paste(..p.label..)),label.x = 0.25, p.digits = 1, r.digits = 2, p.accuracy = 0.001)+
  guides(size = 'none')+
  labs(y = 'Mean growth rates ASV-based \nEubacteria', 
       x = 'Mean growth rates abundance-based\nEubacteria-I,-II,-III probes', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/4), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = 'right', text= element_text(size = 5))

ggsave('eub_asv_gr.pdf', eub_asv_gr, 
       path = "~/Documentos/Doctorat/REMEI/results/figures/corrected_asv_num/",
       width = 88,
       height = 88,
       units = 'mm')

##MOST ABUNDANT NOR5-730 (NOR5/OM60) (normalized)----
x <- mean_abundance_asv_nor%>%
  dplyr::filter(family == 'Halieaceae') %>%
  filter(mean_abund != is.na(mean_abund)) %>%
  #group_by(season, treatment) %>%
  group_by(family) %$%
  mean(mean_abund)

asv_nor5 <-  mean_abundance_asv_nor %>%
  dplyr::filter(family == 'Halieaceae' & #
                  mean_abund > x) %>%
  left_join(reg_all_slopes_chosen_silva_tax_cn_nor_filt, by = c('asv_num', 'season', 'treatment')) %>%
  filter(slope_chosen_days != is.na(slope_chosen_days)) %>%
  group_by(season, treatment) %>%
  summarize(mean_asv = mean(slope_chosen_days))

dapis_nor5 <-  GR_dapis_OS_filt %>%
  dplyr::select(treatment, replicate, season, NOR5) %>%
  group_by(treatment, season) %>%
  dplyr::summarize(mean_dapis_nor5 = mean(as.numeric(NOR5))) %>%
  dplyr::filter(mean_dapis_nor5 >0)

nor5_gr <-  asv_nor5%>%
  right_join(dapis_nor5) 

nor5_asv_gr_cn <- nor5_gr  %>%
  #filter(season != 'Summer') %>%
  ggplot(aes(mean_dapis_nor5, mean_asv))+
  geom_point(aes(shape = treatment, color = season, size = 0.5))+
  scale_fill_manual(values = palette_seasons_4)+
  #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
  # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  #stat_poly_eq(color = 'black')+
  stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
  guides(size = 'none')+
  labs(y = 'Mean growth rates ASV-based \nHaliaceae', 
       x = 'Mean growth rates abundance-based NOR5-730', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/3), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = 'none', text= element_text(size = 5))

##MOST ABUNDANT NOR5-730 (NOR5/OM60) ( no normalized)----
x <- mean_abundance_asv %>%
  dplyr::filter(family == 'Halieaceae') %>%
  filter(mean_abund != is.na(mean_abund)) %>%
  #group_by(season, treatment) %>%
  group_by(family) %$%
  mean(mean_abund)

asv_nor5 <-  mean_abundance_asv %>%
  dplyr::filter(family == 'Halieaceae' & #
                  mean_abund > x) %>%
  left_join(reg_all_slopes_chosen_silva_tax_filt, by = c('asv_num', 'season', 'treatment'), multiple = 'all') %>%
  filter(slope_chosen_days != is.na(slope_chosen_days)) %>%
  group_by(season, treatment) %>%
  dplyr::summarize(mean_asv_nor5 = mean(slope_chosen_days))

dapis_nor5 <-  GR_dapis_OS_filt %>%
  dplyr::select(treatment, replicate, season, NOR5) %>%
  group_by(treatment, season) %>%
  dplyr::summarize(mean_dapis_nor5 = mean(as.numeric(NOR5))) %>%
  dplyr::filter(mean_dapis_nor5 >0)

nor5_gr_nn <-  asv_nor5%>%
  right_join(dapis_nor5) 

nor5_asv_gr <- nor5_gr %>%
  #filter(season != 'Summer') %>%
  ggplot(aes(mean_dapis_nor5, mean_asv))+
  geom_point(aes(shape = treatment, color = season, size = 0.5))+
  scale_fill_manual(values = palette_seasons_4)+
  #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
  # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  #stat_poly_eq(color = 'black')+
  stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
  guides(size = 'none')+
  labs(y = 'Mean growth rates ASV-based \nHaliaceae', 
       x = 'Mean growth rates abundance-based NOR5-730', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/3), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = 'none', text= element_text(size = 5))

##PANNEL CONSTRUCTION--------
#legend <- cowplot::get_legend(nor5_asv_gr) #extract legend
card_fish_asv_gr_abund  <- multi_panel_figure(columns = 2, rows = 7, width = 88, height = 325, 
                                        #row_spacing = 0, #unit = 'mm',
                                        #column_spacing = 0.1, 
                                     # text_grob(label = 'upper-alpha', x = 1, y = 0, just = label_just, face = "bold", size = 14),
                                        panel_label_type = 'upper-alpha'
                                        #font.label = list(size = 5)
                                        )

#panel_label <- text_grob(label = label, x = 1, y = 0, just = label_just, face = "bold", size = 14)

card_fish_asv_gr_abund  %<>%
  fill_panel(alt1413_asv_gr, row = 1, col = 1) %<>%
  fill_panel(alt1413_asv_gr_nor, row = 1, col = 2) %<>%
  fill_panel(ros537_asv_gr, row = 2, col = 1) %<>%
  fill_panel(ros537_asv_gr_cn, row = 2, col = 2) %<>%
  fill_panel(cf319a_asv_gr, row = 3, col = 1) %<>%
  fill_panel(cf319a_asv_gr_cn, row = 3, col = 2) %<>%
  fill_panel(gamma42a_asv_gr, row = 4, col = 1) %<>%
  fill_panel(gamma42a_asv_gr_nor, row = 4, col = 2) %<>%
  #fill_panel(sar11_411_clade_asv_gr, row = 5, col = 1) %<>%
  #fill_panel(sar11_411_clade_asv_gr_cn, row = 5, col = 2) %<>%
  fill_panel(nor5_asv_gr, row = 5, col = 1) %<>%
  fill_panel(nor5_asv_gr_cn, row = 5, col = 2) %<>%
  fill_panel(eub_asv_gr, row = 6, col = 1) %<>%
  fill_panel(eub_asv_gr_cn, row = 6, col = 2) %<>%
  fill_panel(legend, column = 1:2, row = 7)

ggsave('card_fish_asv_gr_cn_nor_non_abund.pdf', card_fish_asv_gr_abund, 
       path = "~/Documentos/Doctorat/REMEI/results/figures/corrected_asv_num/",
       width = 120,
       height = 325,
       units = 'mm')

##Correlation table for supplementary material (with no-normalized data)----
### https://rstudio-pubs-static.s3.amazonaws.com/240657_5157ff98e8204c358b2118fa69162e18.html good tutorial 
##edit colnames to identify which group of ASV are we referring to:
#colnames(cfb_gr) <- c("season", "treatment", "mean_asv_cfb","mean_dapis_cfb")

eub_gr_nn %>%
  dim() ##has the 24 conditions I use this dataframe as a template

corr_dapis_asv_gr <- eub_gr_nn %>%
  left_join(gamma_gr_nn) %>%
  left_join(gr_alt_nn) %>%
  left_join(rhodo_gr_nn) %>%
  left_join(cfb_gr_nn) %>%
  left_join(nor5_gr_nn)
  
corr_dapis_asv_gr %>%
  colnames()

##considerations before computing the correlation (ARE THEY FOLLOWING A NORMAL DISTRIBUTION?)
corr_dapis_asv_gr %>%
str()

ggscatter(corr_dapis_asv_gr, x = "mean_asv_eub", y = "mean_dapis_eub",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "mean gr ASV eub", ylab = "mean DAPIS eub")

# Shapiro-Wilk normality test for all variables pvalue > 0.05 normality can be assumed, if pvalue < 0.05 then no-normality
shapiro.test(as.numeric(corr_dapis_asv_gr$mean_asv_eub)) # => p-value = 0.28
ggqqplot(as.numeric(corr_dapis_asv_gr$mean_asv_eub))

shapiro.test(as.numeric(corr_dapis_asv_gr$mean_dapis_eub)) # =>  p-value = 0.1033
ggqqplot(as.numeric(corr_dapis_asv_gr$mean_asv_eub))

shapiro.test(as.numeric(corr_dapis_asv_gr$mean_asv_gam)) # => p-value = 0.07014
ggqqplot(as.numeric(corr_dapis_asv_gr$mean_asv_gam))

shapiro.test(as.numeric(corr_dapis_asv_gr$mean_dapis_gamma)) # => p-value = 0.0007678 (NO-NORMALITY)
ggqqplot(as.numeric(corr_dapis_asv_gr$mean_dapis_gamma))

shapiro.test(as.numeric(corr_dapis_asv_gr$mean_asv_alt)) # => p-value = 0.847
ggqqplot(as.numeric(corr_dapis_asv_gr$mean_asv_alt))

shapiro.test(as.numeric(corr_dapis_asv_gr$mean_dapis_alt)) # =>  p-value = 0.03223 (NO-NORMALITY)
ggqqplot(as.numeric(corr_dapis_asv_gr$mean_asv_alt))

shapiro.test(as.numeric(corr_dapis_asv_gr$mean_asv_alt)) # => p-value = 0.07014
ggqqplot(as.numeric(corr_dapis_asv_gr$mean_asv_alt))

shapiro.test(as.numeric(corr_dapis_asv_gr$mean_dapis_alt)) # => p-value = 0.0007678 (NO-NORMALITY)
ggqqplot(as.numeric(corr_dapis_asv_gr$mean_dapis_alt))

shapiro.test(as.numeric(corr_dapis_asv_gr$mean_asv_rhodo)) # =>  p-value = 0.223
ggqqplot(as.numeric(corr_dapis_asv_gr$mean_asv_rhodo))

shapiro.test(as.numeric(corr_dapis_asv_gr$mean_dapis_rhodo)) # => p-value = 0.4506
ggqqplot(as.numeric(corr_dapis_asv_gr$mean_dapis_rhodo))

shapiro.test(as.numeric(corr_dapis_asv_gr$mean_asv_cfb)) # => p-value = 0.3883
ggqqplot(as.numeric(corr_dapis_asv_gr$mean_asv_cfb))

shapiro.test(as.numeric(corr_dapis_asv_gr$mean_dapis_cfb)) # => p-value = 0.003669 (NO-NORMALITY)
ggqqplot(as.numeric(corr_dapis_asv_gr$mean_dapis_cfb))

shapiro.test(as.numeric(corr_dapis_asv_gr$mean_asv_nor5)) # => p-value = 0.08854
ggqqplot(as.numeric(corr_dapis_asv_gr$mean_asv_nor5))

shapiro.test(as.numeric(corr_dapis_asv_gr$mean_dapis_nor5)) # => p-value = 0.05104
ggqqplot(as.numeric(corr_dapis_asv_gr$mean_dapis_nor5))

##one by one----
cor(corr_dapis_asv_gr$mean_asv_eub, corr_dapis_asv_gr$mean_dapis_eub, method = 'kendall')
res <- cor.test(corr_dapis_asv_gr$mean_asv_eub, corr_dapis_asv_gr$mean_dapis_eub, method = 'spearman')
res$p.value


##all at the same time----
##for normal data
##function to extract a table
flat_cor_mat <- function(cor_r, cor_p){
  #This function provides a simple formatting of a correlation matrix
  #into a table with 4 columns containing :
  # Column 1 : row names (variable 1 for the correlation test)
  # Column 2 : column names (variable 2 for the correlation test)
  # Column 3 : the correlation coefficients
  # Column 4 : the p-values of the correlations
  library(tidyr)
  library(tibble)
  cor_r <- rownames_to_column(as.data.frame(cor_r), var = "row")
  cor_r <- gather(cor_r, column, cor, -1)
  cor_p <- rownames_to_column(as.data.frame(cor_p), var = "row")
  cor_p <- gather(cor_p, column, p, -1)
  cor_p_matrix <- left_join(cor_r, cor_p, by = c("row", "column"))
  cor_p_matrix
}

cor_pearson <- rcorr(as.matrix(corr_dapis_asv_gr[,3:12]), type = 'pearson')

my_cor_matrix <- flat_cor_mat(cor_pearson$r, cor_pearson$P)
head(my_cor_matrix)

my_cor_matrix_pear <- my_cor_matrix %>%
  as_tibble() %>%
  # filter(row %in% c('mean_dapis_gamma', 'mean_asv_gamma', 'mean_dapis_nor5', 'mean_asv_nor5')) %>%
  # filter(column != c('mean_dapis_gamma', 'mean_asv_gamma', 'mean_dapis_nor5', 'mean_asv_nor5')) %>%
  filter(row == 'mean_asv_eub' & column == 'mean_dapis_eub' |
         row == 'mean_asv_rhodo' & column == 'mean_dapis_rhodo' |
           row == 'mean_asv_cfb' & column == 'mean_dapis_cfb' )

##for non-normal data (pvalue < 0.05 in Shapiro test)
cor_spearman <- rcorr(as.matrix(corr_dapis_asv_gr[,3:14]), type = 'spearman')

my_cor_matrix_spear <- flat_cor_mat(cor_spearman$r, cor_spearman$P)
head(my_cor_matrix_spear)

my_cor_matrix_spear <- my_cor_matrix_spear %>%
  as_tibble() %>%
  # filter(row %in% c('mean_dapis_gamma', 'mean_asv_gamma', 'mean_dapis_nor5', 'mean_asv_nor5')) %>%
  # filter(column != c('mean_dapis_gamma', 'mean_asv_gamma', 'mean_dapis_nor5', 'mean_asv_nor5')) %>%
  filter(row == 'mean_asv_alt' & column == 'mean_dapis_alt' |
           row == 'mean_asv_gam' & column == 'mean_dapis_gamma' |
           row == 'mean_asv_nor5' & column == 'mean_dapis_nor5' )

##table with all the results for the correlations
my_correlation <- my_cor_matrix_pear %>%
  rbind(my_cor_matrix_spear)

colnames(my_correlation) = c('Mean ASV-based growth rates', 'Abundance-based growth rates', 'Correlation', 'p-value')

write.csv2(my_correlation, file = 'results/tables/correlations_card-fish_avs_gr.csv')
