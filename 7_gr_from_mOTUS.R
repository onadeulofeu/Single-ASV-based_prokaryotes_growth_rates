##Question: Are we overestimating GR for some groups with these method?
#Question: comparing GR from mOTUs with ASV-based growth rates

library(readxl)
library(tidyverse)
library(magrittr)
library(ggplot2)
library(phyloseq)
library(microbiome)
library(scales)
library(multipanelfigure)
library(seqinr)
library(Biostrings)
library(speedyseq) #mutate_tax_table
library(janitor)

#import data
motus_table <- read_tsv('data/mOTUs_Adrià/motus_table_ed.txt')
head(motus_table) 

motus_table %>%
  dim() ##33570 mOTUs and 66 samples

tax <- motus_table %>%
  select(consensus_taxonomy, mOTUs_num, consensus_taxonomy_mOTUs_num)

motus_table_ed <- motus_table %>%
  select(-consensus_taxonomy, -mOTUs_num)

##env data
env_data_mOTUs <- read_xlsx('data/envdata/metadata_remei_mOTUs.xlsx')

env_data_mOTUs %>%
  colnames()

env_data_mOTUs <- env_data_mOTUs %>% 
  rename(c("Treatment" = "treatment", "Light Regime" = "light_regime",
           "Replicate" = "replicate", "Time" = "time", "Hours" = "hours", "Season" = "season", 
           "Fraction" = "fraction", "Sample-Name" = "sample_name" ,"Sample Code" = "sample_code",  
           "amplicon_data_correspondence" = "amplicon_data_correspondence", "r1" = "r1", "r2" = "r2", 
           "Million\r\nread-pairs\r\n(or reads)" = 'reads', "Yield\r\n(Gb)" = 'yield'))

##edit t0 per duplicar t0 dark per light treatment també
change_t0s <- function(df){
  datch <- df %>% 
    filter( treatment %in% c('CT', 'PF'))
  
  
  datcl <- datch %>%  
    mutate( treatment = case_when( treatment == 'CT' ~ 'CL',
                                   treatment == 'PF' ~ 'PL',
                                   TRUE ~ treatment))
  
  
  datd <- datch %>%  
    mutate( treatment = case_when( treatment == 'CT' ~ 'CD',
                                   treatment == 'PF' ~ 'PD',
                                   TRUE ~ treatment))
  
  whole <- bind_rows(datcl, datd) %>% 
    mutate(selected_file_name_new = str_c('REMEIextra', 1:nrow(.)))
  
  df %>% 
    filter( !treatment %in% c('CT', 'PF')) %>% 
    bind_rows(whole) %>% 
    return()
}

env_data_mOTU_ed <- change_t0s(env_data_mOTUs)

env_data_mOTU_ed %>%
  colnames()

##add citometry to calculate pseudoabundances for mOTUs
rem_fc <- readRDS("./data/intermediate_files/remei_phyloseq_silva_fc.rds")
rem_fc_env <- rem_fc@sam_data %>%
  as_tibble()
rem_fc_env %>%
  colnames()

rem_fc_env_ed <- rem_fc_env %>%
  mutate(treatment = as.character(treatment),
         season = as.character(season)) %>%
  select(treatment, season, All, replicate, time, hours.x)

rem_fc_env %>%
  glimpse()

env_data_mOTU_ed %>%
  glimpse()

env_data_mOTU_ed_ct <- env_data_mOTU_ed %>%
  left_join(rem_fc_env_ed,
            by = c('treatment' = 'treatment', 
                   'season' = 'season', 
                   'replicate' = 'replicate', 
                   'time' = 'time'))
  

#dataset information----- 



##filter environmental dataset (failed)
env_data_mOTU_ed_ct <- 
  env_data_mOTU_ed_ct %>%
  filter(!is.na(amplicon_data_correspondence), amplicon_data_correspondence != "FAILED")

##filter mOTUs that are 0 in the whole dataset:
motus_table_ed %>%
  dim()

#filter mOTUs that sum 0
motus_table_ed %$%
  consensus_taxonomy_mOTUs_num %>%
  unique()

motus_table_ed %>%
  dim()

motus_table_ed %>%
  colnames()

motus_table_t <- motus_table_ed %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = 'sample_code') %>%
  row_to_names(row_number = 1, remove_row = T, remove_rows_above = F) 


motus_table_t <- motus_table_t %>%
  mutate(across(!consensus_taxonomy_mOTUs_num, as.numeric))

# motus_table_t %>%
#   select(!consensus_taxonomy_mOTUs_num) %>%
#   colSums()

motus_table_t_codes <- motus_table_t$consensus_taxonomy_mOTUs_num

motus_table_t_filt <- motus_table_t %>%
  select(-consensus_taxonomy_mOTUs_num) %>%
  select(where(~ sum(.) != 0)) %>%
  cbind(motus_table_t_codes)

motus_table_t_filt %>%
  colnames()

##add hours information
motus_env <- env_data_mOTU_ed_ct %>%
  select(hours, sample_code, time, season, treatment, replicate, yield, All) %>%
  left_join(motus_table_t_filt, by = c('sample_code' = 'motus_table_t_codes'))

env_data_mOTUs$sample_code == motus_table_t_filt$motus_table_t_codes

motus_table_t_filt %>%
  glimpse()

##pivot long to compare t0 and t4 and divide keep only positive values
motus_env %>%
  colnames() %>%
  head()

motus_env %$%
  treatment %>%
  unique()

motus_env %$%
  consensus_taxonomy_mOTUS_num %>%
  unique()

motus_env %$%
  time %>%
  unique()
motus_env %>%
  colnames()

motus_env_l <- motus_env %>%
  #select(-hours) %>%
  pivot_longer(!c(sample_code, time, season, treatment, replicate, hours, yield, All), 
               names_to = 'consensus_taxonomy_mOTUs_num', 
               values_to = 'rel_abund') %>%
  mutate(time_ed = str_replace_all(time, c('t4'= 'tf', 't3' = 'tf')))

###rel abund per citometria aquí!!!! per tenir pseudoabundances de les mOTUs!
motus_env_l_w <- motus_env_l %>%  
  select(treatment, season, consensus_taxonomy_mOTUs_num, All, time_ed, replicate, yield, rel_abund) %>%
  mutate(rel_abund_ed = ifelse(is.na(rel_abund), 0, rel_abund)) %>%
  mutate(pseudoabund = All*rel_abund_ed) %>%
  # subset(season == 'Winter' & 
  #          treatment == 'CD') %>%
  #select(-treatment, -season) %>%
  #group_by(season, treatment, consensus_taxonomy_mOTUs_num) %>%
  pivot_wider(id_cols = c(treatment, season, consensus_taxonomy_mOTUs_num, replicate), names_from = c(time_ed), values_from = rel_abund_reads)

##**Alternative linear regression using replicates to have pvalue (go to line 389)------
motus_env_l_w_hours <- env_data_mOTU_ed %>%
  select(treatment, season, replicate, hours, time) %>%
  mutate(time_ed = str_replace_all(time, c('t4'= 'tf', 't3' = 'tf'))) %>%
  subset(time_ed == 'tf') %>%
  left_join(motus_env_l_w)

motus_env_l_w_h_growth <- motus_env_l_w_hours %>%
  mutate(growth = (tf-t0)/hours)

#mean growth for each asv in replicates (no surt bé!, plotejar sense mean per rèpliques)
motus_env_l_w_h_growth_mean <- motus_env_l_w_h_growth %>%
  group_by(treatment, season, hours, time, consensus_taxonomy_mOTUs_num) %>%
  mutate(mean_growth_asv = mean(growth_days)) %>%
  distinct(treatment, season, hours, time, consensus_taxonomy_mOTUs_num, mean_growth_asv)

motus_env_l_w_h_growth_mean$consensus_taxonomy_mOTUs_num %>%
  unique() ##2715
  
##reorder treatments and seasons for plots
motus_env_l_w_h_growth_mean$treatment <- motus_env_l_w_h_growth_mean$treatment %>% 
  factor(levels=(c("CL", "CD", "PL", "PD", "DL", "VL")))
motus_env_l_w_h_growth_mean$season <- motus_env_l_w_h_growth_mean$season %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))

##transform slopes /hours to /days
motus_env_l_w_h_growth <- motus_env_l_w_h_growth %>%
  mutate(growth_days = growth*24)

##select only the most abundant >1% community at some point
remiau_relabun_melt <-  motus_table_ed %>%
  pivot_longer(cols = !consensus_taxonomy_mOTUs_num, values_to = 'rel_abund', names_to = 'sample_code')

remiau_relabun_melt %$%
  consensus_taxonomy_mOTUs_num %>%
  unique() #5111 asv al meu dataset

abundant_motus <- remiau_relabun_melt %>% 
  filter(rel_abund > 0.025) %>% #more than 1% of the community at some point
  dplyr::select(consensus_taxonomy_mOTUs_num) %>%
  unique() %>%
  as_tibble()

motus_env_l_w_h_growth_1perc <- motus_env_l_w_h_growth_mean %>% 
  right_join(abundant_motus, by = "consensus_taxonomy_mOTUs_num", copy = TRUE)

motus_gr_mean <- motus_env_l_w_h_growth_1perc %>%
  filter(growth > 0,
         consensus_taxonomy_mOTUs_num != 'unassigned') %>%
  group_by(treatment, season) %>% 
  mutate(mean_gr = mean(growth)) %>%
  distinct(treatment, season, mean_gr)

motus_table_ed %>%
  select(-consensus_taxonomy_mOTUs_num) %>%
  colSums()

##Single mOTUs based growth rates (FIGURE 1)----
motus_env_l_w_h_growth_1perc %>%
  colnames()

#statistics seasons----
# reg_all_slopes_chosen_silva_tax_1perc %>%
#   colnames()
# anova<-aov(slope_chosen_days ~season, data = reg_all_slopes_chosen_silva_tax_1perc)
# summary(anova) #p<0.05 ***
# TukeyHSD(anova)#para ver los grupos que son significativamente distintos
# aov_residuals <- residuals(object = anova)
# # plot(anova, 1)
# # plot(anova, 2)
# shapiro.test(x = aov_residuals )#From the output, the p-value > 0.05 implying that the distribution of the data are 
# #not significantly different from normal distribution. In other words, we can assume the normality.
# ##p-value < 2.2e-16 NO és normal 
# library(FSA)##test the d
# dunnTest(slope_chosen_days ~ season, data = reg_all_slopes_chosen_silva_tax_1perc,
#          method= 'bonferroni')
# results<-dunnTest(slope_chosen_days ~ season, data = reg_all_slopes_chosen_silva_tax_1perc,
#                   method='bonferroni')$res
# results<-results[results$P.adj<0.05,]
# #Differents: fall-spring, fall-summer, fall-winter, spring-winter, summer-winter
# X = results$P.adj <= 0.05
# names(X) = gsub(" ",  "",  results$Comparison)
# multcompLetters(X)

# Fall Spring Summer Winter 
# "a"    "b"    "b"    "c" 
#statistical groups for non-parammetric test
# letters_seas <- data.frame(multcompLetters(X)['Letters'])

#plot seasons boxplot----
#seas <- 
  motus_env_l_w_h_growth_1perc %>%
  filter(mean_growth_asv > 0 &
           consensus_taxonomy_mOTUs_num != 'unassigned') %>% 
  ggplot(aes(season, mean_growth_asv))+ #
  geom_point(aes(color = treatment), alpha = 0.7, position = position_jitter(0.25))+
  #geom_boxplot(alpha= 0.5)+
  geom_violin(alpha = 0.2, draw_quantiles = c(0.25, 0.75) )+
  labs(color = "Treatment")+
  stat_summary(fun = "mean", geom = "crossbar", 
               width = 0.6, colour = "black")+
  # geom_text(data = letters_seas, aes(y = 9, x = row.names(letters_seas), label = Letters),
  #           position = position_nudge(x = 0.2), hjust = 0.7, color = "black")+
  labs(x= "Season", y = expression("Growth rate day"^"-1"))+  #(μ) 
  #scale_x_discrete()+
  #scale_y_continuous(expand = c(0, 0), limits = c(0,10.4))+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_treatments_remei)+
  #facet_wrap(.~phylum, scales = "free")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#statistics treatments----
##Treatments
# anova<-aov(slope_chosen_days ~ treatment, data = reg_all_slopes_chosen_silva_tax_1perc)
# summary(anova) #p<0.05 ***
# TukeyHSD(anova)#para ver los grupos que son significativamente distintos
# aov_residuals <- residuals(object = anova)
# # plot(anova, 1)
# # plot(anova, 2)
# shapiro.test(x = aov_residuals )#From the output, the p-value > 0.05 implying that the distribution of the data are #not significantly different from normal distribution. In other words, we can assume the normality.
# ##p-value < 2.2e-16 NO és normal 
# #si no es normal:
# #hago test no parametrico
# kruskal.test(slope_chosen_days ~ treatment, data = reg_all_slopes_chosen_silva_tax_1perc) #
# #si sale p<0.05 hago dunn test para ver cuales son significativamente distintos
# library(FSA)##test the d
# dunnTest(slope_chosen_days ~ treatment, data = reg_all_slopes_chosen_silva_tax_1perc,
#          method= 'bonferroni')
# results<-dunnTest(slope_chosen_days ~ treatment, data = reg_all_slopes_chosen_silva_tax_1perc,
#                   method='bonferroni')$res
# #results<-results[results$P.adj<0.05,]
# ##Diferents: CD-PD, PD-VL
# X = results$P.adj <= 0.05
# 
# names(X) = gsub(" ",  "",  results$Comparison)
# multcompLetters(X)
# # CD   CL   DL   PD   PL   VL 
# # "a" "ab" "ab"  "b" "ab"  "a"
# letters_treat <- data.frame(multcompLetters(X)['Letters'])

#treat <- 
  motus_env_l_w_h_growth_1perc %>%
  filter(mean_growth_asv > 0.0 &
           consensus_taxonomy_mOTUs_num != 'unassigned') %>% 
  ggplot(aes(treatment, mean_growth_asv))+
  geom_point(aes(treatment, mean_growth_asv, color = season), alpha = 0.7, position = position_jitter(0.25))+ #position_jitter(0.1)  #, position = "identity"
  geom_violin(alpha = 0.2, draw_quantiles = c(0.25, 0.75), )+
  #geom_boxplot()+
  labs(color = "Season")+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.6,
               colour = "black")+
  # geom_text(data = letters_treat, aes(y = 9, x = row.names(letters_treat), label = Letters),
  #           position = position_nudge(x = 0.2), hjust = 0.7, color = "black")+
  labs(x= 'Season', y = expression("Growth rate day"^"-1"), color = "Season")+ #(μ)
  #scale_x_discrete()+
  #scale_y_continuous(expand = c(0, 0), limits = c(0,10.4))+ #forçar que comenci al 0 l'eix.
  #geom_smooth(method = "lm")
  scale_color_manual(values = palette_seasons_4)+
  #facet_wrap(.~phylum, scales = "free")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())


##Regression with DAPI-based growth rates & ASV-based growth rates-------
##ASV-based vs mOTUS-based gr
##asv gr mean filtered by negative values, NAs and unsignificant
general_mean_gr_asv_filt <- reg_all_slopes_chosen_silva_tax_1perc %>%
  filter(season != "Early_fall") %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  #group_by(treatment, season) %>%
  mutate(mean_asv = mean(as.numeric(slope_chosen_days))) %>%
  distinct(mean_asv)

mean_gr_asv_filt <- reg_all_slopes_chosen_silva_tax_1perc %>%
  filter(season != "Fall") %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  group_by(treatment, season) %>%
  mutate(mean_asv = mean(as.numeric(slope_chosen_days))) %>%
  distinct(season, treatment, mean_asv)

mean_gr_asv_filt_motus <-  mean_gr_asv_filt %>%
  left_join(motus_gr_mean)


library(ggpubr)

#rel_gr_asv_motus <-
  mean_gr_asv_filt_motus %>%
  as_tibble() %>%
  ggplot(aes(mean_asv, mean_gr))+
  #geom_smooth(method = 'lm', se = FALSE, color = 'grey')+
  geom_point(aes(shape = treatment, color = season, size = 2))+
  #scale_fill_manual(values = palette_seasons_4)+
  # scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
  # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  #stat_poly_eq(color = 'black')+
  stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
  guides(size = 'none')+
  labs(x = 'Mean growth rates ASV based \n
       (growing community)', y = 'Mean growth rates DAPI based', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/3), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

##comparison with dapis gr
mean_gr_dapis <- GR_dapis_OS_filt %>%
  group_by(treatment, season) %>%
  mutate(mean_dapis = mean(as.numeric(PRO))) %>%
  distinct(season, treatment, mean_dapis)

mean_gr_asv_motus_dapis <-  mean_gr_asv_filt_motus %>%
  left_join(mean_gr_dapis)

mean_gr_asv_motus_dapis %>%
  as_tibble() %>%
  ggplot(aes(mean_gr, mean_dapis))+
  #geom_smooth(method = 'lm', se = FALSE, color = 'grey')+
  geom_point(aes(shape = treatment, color = season, size = 2))+
  #scale_fill_manual(values = palette_seasons_4)+
  # scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
  # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  #stat_poly_eq(color = 'black')+
  stat_cor(aes(label = paste(..p.label..)), label.x = 0.0023, p.digits = 1, r.digits = 2)+
  guides(size = 'none')+
  labs(x = 'Mean growth rates mOTUs based \n
       (growing community)', y = 'Mean growth rates DAPI based', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/3), panel.grid.major = element_blank(), panel.grid.minor = element_blank())


###linear regression using replicates-------
motus_env_l_w %>%
  colnames()


##adapted from script 3:
#1. Prepare  datasets for calculating regressions------
motus_env_pseudo <- motus_env_l %>%  
  select(treatment, season, consensus_taxonomy_mOTUs_num, All, time_ed, replicate, yield, rel_abund) %>%
  mutate(rel_abund_ed = ifelse(is.na(rel_abund), 0, rel_abund)) %>%
  mutate(pseudoabund = All*rel_abund_ed)



motus_env_pseudo_hours <- env_data_mOTU_ed %>%
  select(treatment, season, replicate, hours, time) %>%
  mutate(time_ed = str_replace_all(time, c('t4'= 'tf', 't3' = 'tf'))) %>%
  # #subset(time_ed == 'tf') %>%
   left_join(motus_env_pseudo, by = c('season' = 'season', 'treatment' = 'treatment', 
                                   'replicate' = 'replicate', 'time_ed' = 'time_ed')) 

motus_env_pseudo_hours  %>%
  colnames()

###crear una nova variable que combina time and replicate per no tenir problemes amb els temps repetits de les rèpliques.
data <- motus_env_pseudo_hours %>% #pseudoabundances dataset form fc data in wide format with metadata.
  mutate(t_rep = paste(time_ed, replicate, sep="_")) %>%
  #select(-time, -time_ed) %>%
  #pivot_longer(cols = c('t0', 'tf'), values_to = 'pseudoabund', names_to = 'time') %>%
  mutate(across(pseudoabund, ~replace(., . == 0, "1"))) %>% ##cambiem 0 per 1 per poder fer la regressió
  mutate(pseudoabund = as.numeric(pseudoabund))

env <- env_data_mOTU_ed %>%
  as_tibble() %>%
  mutate(t_rep = paste(time, replicate, sep="_"))


###testing linear regressions for long data
library(dplyr)
library(purrr)
library(broom)

fit_model <- function(df) lm(hours ~ log(pseudoabund), data = df)
get_slope <- function(mod) tidy(mod)$estimate[2]
get_p_value <- function(mod) tidy(mod)$p.value[2]

fit_model(data)

reg_mOTUs <- data %>% 
  group_nest(treatment, season, consensus_taxonomy_mOTUs_num) %>% 
  mutate(model = map(data, fit_model),
         slope = map_dbl(model, get_slope),
         p_value = map_dbl(model, get_p_value))

# data %>%
#   nest_by(treatment, season, consensus_taxonomy_mOTUs_num) %>% 
#   mutate(model = list(lm(hours ~ as.numeric(pseudoabund), data))) %>% 
#   summarise(tidy(model)) %>% 
#   ungroup()

reg_mOTUs_filt <- reg_mOTUs %>%
  as_tibble() %>%
  select(treatment, season, slope, p_value, consensus_taxonomy_mOTUs_num) %>%
  mutate(pvalue = as.numeric(p_value)) %>%
  filter(slope > 0 &
           pvalue < 0.01 &
           consensus_taxonomy_mOTUs_num != 'unassigned')

##reorder treatments and seasons for plots
reg_mOTUs_filt$treatment <- reg_mOTUs_filt$treatment %>% 
  factor(levels=(c("CL", "CD", "PL", "PD", "DL", "VL")))
reg_mOTUs_filt$season <-reg_mOTUs_filt$season %>% 
  factor(levels=(c("Winter", "Spring", "Summer")))

#seas <- 
reg_mOTUs_filt %>%
  ggplot(aes(season, slope))+ #
  geom_point(aes(color = treatment), alpha = 0.7, position = position_jitter(0.25))+
  #geom_boxplot(alpha= 0.5)+
  geom_violin(alpha = 0.2, draw_quantiles = c(0.25, 0.75) )+
  labs(color = "Treatment")+
  stat_summary(fun = "mean", geom = "crossbar", 
               width = 0.6, colour = "black")+
  # geom_text(data = letters_seas, aes(y = 9, x = row.names(letters_seas), label = Letters),
  #           position = position_nudge(x = 0.2), hjust = 0.7, color = "black")+
  labs(x= "Season", y = expression("Growth rate day"^"-1"))+  #(μ) 
  #scale_x_discrete()+
  #scale_y_continuous(expand = c(0, 0), limits = c(0,10.4))+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_treatments_remei)+
  #facet_wrap(.~phylum, scales = "free")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#statistics treatments----
##Treatments
# anova<-aov(slope_chosen_days ~ treatment, data = reg_all_slopes_chosen_silva_tax_1perc)
# summary(anova) #p<0.05 ***
# TukeyHSD(anova)#para ver los grupos que son significativamente distintos
# aov_residuals <- residuals(object = anova)
# # plot(anova, 1)
# # plot(anova, 2)
# shapiro.test(x = aov_residuals )#From the output, the p-value > 0.05 implying that the distribution of the data are #not significantly different from normal distribution. In other words, we can assume the normality.
# ##p-value < 2.2e-16 NO és normal 
# #si no es normal:
# #hago test no parametrico
# kruskal.test(slope_chosen_days ~ treatment, data = reg_all_slopes_chosen_silva_tax_1perc) #
# #si sale p<0.05 hago dunn test para ver cuales son significativamente distintos
# library(FSA)##test the d
# dunnTest(slope_chosen_days ~ treatment, data = reg_all_slopes_chosen_silva_tax_1perc,
#          method= 'bonferroni')
# results<-dunnTest(slope_chosen_days ~ treatment, data = reg_all_slopes_chosen_silva_tax_1perc,
#                   method='bonferroni')$res
# #results<-results[results$P.adj<0.05,]
# ##Diferents: CD-PD, PD-VL
# X = results$P.adj <= 0.05
# 
# names(X) = gsub(" ",  "",  results$Comparison)
# multcompLetters(X)
# # CD   CL   DL   PD   PL   VL 
# # "a" "ab" "ab"  "b" "ab"  "a"
# letters_treat <- data.frame(multcompLetters(X)['Letters'])

#treat <- 
reg_mOTUs_filt %>%
  ggplot(aes(treatment, slope))+
  geom_point(aes(treatment, slope, color = season), alpha = 0.7, position = position_jitter(0.25))+ #position_jitter(0.1)  #, position = "identity"
  geom_violin(alpha = 0.2, draw_quantiles = c(0.25, 0.75), )+
  #geom_boxplot()+
  labs(color = "Season")+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.6,
               colour = "black")+
  # geom_text(data = letters_treat, aes(y = 9, x = row.names(letters_treat), label = Letters),
  #           position = position_nudge(x = 0.2), hjust = 0.7, color = "black")+
  labs(x= 'Season', y = expression("Growth rate day"^"-1"), color = "Season")+ #(μ)
  #scale_x_discrete()+
  #scale_y_continuous(expand = c(0, 0), limits = c(0,10.4))+ #forçar que comenci al 0 l'eix.
  #geom_smooth(method = "lm")
  scale_color_manual(values = palette_seasons_4)+
  #facet_wrap(.~phylum, scales = "free")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

##Regression with DAPI-based growth rates & ASV-based growth rates-------
##mean regressions mOTUS
mOTUs_relabun_melt <-  motus_table_ed %>%
  pivot_longer(cols = !consensus_taxonomy_mOTUs_num, values_to = 'rel_abund', names_to = 'sample_code')

mOTUs_relabun_melt %$%
  consensus_taxonomy_mOTUs_num %>%
  unique() #5111 asv al meu dataset

abundant_motus <- mOTUs_relabun_melt %>% 
  filter(rel_abund > 0.025) %>% #more than 1% of the community at some point
  dplyr::select(consensus_taxonomy_mOTUs_num) %>%
  unique() %>%
  as_tibble()

motus_env_l_w_h_growth_1perc <- reg_mOTUs_filt %>% 
  right_join(abundant_motus, by = "consensus_taxonomy_mOTUs_num", copy = TRUE)

motus_gr_mean <- motus_env_l_w_h_growth_1perc %>%
  # filter(slope > 0,
  #        pvalue < 0.05,
  #        consensus_taxonomy_mOTUs_num != 'unassigned') %>%
  group_by(treatment, season) %>% 
  mutate(mean_gr = mean(slope)) %>%
  distinct(treatment, season, mean_gr)
##ASV-based vs mOTUS-based gr
##asv gr mean filtered by negative values, NAs and unsignificant
general_mean_gr_asv_filt <- reg_all_slopes_chosen_silva_tax_1perc %>%
  filter(season != "Early_fall") %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  #group_by(treatment, season) %>%
  mutate(mean_asv = mean(as.numeric(slope_chosen_days))) %>%
  distinct(mean_asv)

mean_gr_asv_filt <- reg_all_slopes_chosen_silva_tax_1perc %>%
  filter(season != "Fall") %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  group_by(treatment, season) %>%
  mutate(mean_asv = mean(as.numeric(slope_chosen_days))) %>%
  distinct(season, treatment, mean_asv)

mean_gr_asv_filt_motus <-  mean_gr_asv_filt %>%
  left_join(motus_gr_mean)


library(ggpubr)

#rel_gr_asv_motus <-
mean_gr_asv_filt_motus %>%
  as_tibble() %>%
  ggplot(aes(mean_asv, mean_gr))+
  #geom_smooth(method = 'lm', se = FALSE, color = 'grey')+
  geom_point(aes(shape = treatment, color = season, size = 2))+
  #scale_fill_manual(values = palette_seasons_4)+
  # scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
  # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  #stat_poly_eq(color = 'black')+
  stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
  guides(size = 'none')+
  labs(x = 'Mean growth rates ASV based \n
       (growing community)', y = 'Mean growth rates mOTUs based', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/3), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

##comparison with dapis gr
mean_gr_dapis <- GR_dapis_OS_filt %>%
  group_by(treatment, season) %>%
  mutate(mean_dapis = mean(as.numeric(PRO))) %>%
  distinct(season, treatment, mean_dapis)

mean_gr_asv_motus_dapis <-  mean_gr_asv_filt_motus %>%
  left_join(mean_gr_dapis)

mean_gr_asv_motus_dapis %>%
  as_tibble() %>%
  ggplot(aes(mean_gr, mean_dapis))+
  #geom_smooth(method = 'lm', se = FALSE, color = 'grey')+
  geom_point(aes(shape = treatment, color = season, size = 2))+
  #scale_fill_manual(values = palette_seasons_4)+
  # scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
  # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  #stat_poly_eq(color = 'black')+
  stat_cor(aes(label = paste(..p.label..)), label.x = 0.0023, p.digits = 1, r.digits = 2)+
  guides(size = 'none')+
  labs(x = 'Mean growth rates mOTUs based \n
       (growing community)', y = 'Mean growth rates DAPI based', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/3), panel.grid.major = element_blank(), panel.grid.minor = element_blank())


