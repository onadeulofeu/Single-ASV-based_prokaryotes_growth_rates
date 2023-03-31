# This code is designed to calculate single-ASVs growth rates for the paper entitled 
# Growth rates of marine prokayotesare extremely diverse, even among closely related taxa
# Code developed by Ona Deulofeu-Capo.


sessionInfo()

##packages used
library(tidyverse)
library(magrittr)
library(reshape2)
library(plyr)
library(speedyseq)
library(ggplot2)
library(vegan)
library(scales) ##percent format for plot bars
library(microbiome)
library(forcats) #reorder categorical variables
library(gridExtra)
library(multipanelfigure) #merge plots with different sizes
library(ggforce) #ggplot2 extension
library(ggridges) #nice distribution plots 
library(ggpubr) #as_ggplot function, stat_cor
library(rstatix) #statistics
library(multcompView)
library(ggpmisc) ##stat_poli_line and equation
library(stringr)
library(FSA)##test the d
        

#define directory ----
setwd ("~/Documentos/Doctorat/REMEI/")

# importing functions ----
source("src/pseudoabundances_calculation_fc.R") ## for calculating pseudoabundances
source("src/filter.season.treatment.time.R") ## filter for calculating regressions
source("src/multiple.linear.regressions.R") ## multiple linear regressions
source("src/comparing.reg.R") ## comparisons between 3 times and 4 times regressions
source('src/growth.rates.distr.tax.ranks.ridges.divided.R') ##distribution plots
source('src/growth.rates.distr.tax.ranks.ridges.R') ##distribution plots divided by classes

#packages
library(lubridate)
library(scales)

# import data ----
rem_fc <- readRDS("./data/intermediate_files/remei_phyloseq_silva_fc.rds") #phyloseq object with otu_table, sample data and taxonomu
#sample data must have total abundance data to be able to calculate pseudoabundances

#palettes used ----
palette_seasons_4 <- c("Winter" = "#002562", 'Spring' = "#519741", 'Summer' = "#ffb900",'Fall' =  "#96220a")
palette_treatments_remei <- c("#cab2d6","#6a3d9a", "#a6cee3","#1f78b4", "#009E73", "#F2BB05")
palette_phylums_assigned <- c('Proteobacteria' = "#fcca46","Bacteroidota" = "#669bbc" , 'Actinobacteriota' = "#b0413e",
                              'Cyanobacteria' = "#009e73",'Crenarchaeota' = "#ffa737", 'Verrucomicrobiota' = "#cccccc",
                              'Planctomycetota' = "#69267e", 'Acidobacteriota' = "#1f78b4",'Bdellovibrionota' = "#8c789d",
                              'Firmicutes' = "#637b88", 'Myxococcota'= "#003029", 'Nitrospinota'= "#e3a6ce",
                              'Campilobacterota' = "#002960", 'Deinococcota'= "#ba9864",'Fusobacteriota' ="#fb9a99",
                              'Desulfobacterota' = "#005c69", 'Methylomirabilota' = "#000000" ,
                              'Gemmatimonadota' = "#c55e5c", 'Chloroflexi' = "#00d198", 'Spirochaetota' = "#5cb1d6",
                              'Calditrichota' = "#8a007a", 'Halobacterota' = "#b79c64", 'Nitrospirota' = "#41815f",
                              'Dependentiae' = "#5b95e5", 'Patescibacteria' = "#33af9c",'Cloacimonadota' = "#fbed5c",
                              'Synergistota' = "#ce7800", 'Abditibacteriota' = "#87878b", 'Deferribacterota' = "#4dbaa9")
palette_class_assigned <- c('Gammaproteobacteria' = '#fcca46', 'Alphaproteobacteria' = '#d95726', 'Zetaproteobacteria' = '#EBCD92',
                            'Bacteroidia' = '#669BBC','Rhodothermia' =  '#00425E',
                            'Ignavibacteria' = '#CAFFFF', 'Chlorobia' = '#5BD1FF', 'Kryptonia' = '#0071A8',
                            'Nitrososphaeria' = '#FFA737',
                            'Cyanobacteriia' = '#009E73', 'Vampirivibrionia' = '#00733C',
                            'Acidimicrobiia' = '#B0413E','Actinobacteria' = '#C55E5C',
                            'Coriobacteriia' = '#B44240', 'Thermoleophilia' = '#AB3F3D',
                            'Verrucomicrobiae' = '#CCCCCC', 'Kiritimatiellae' = '#6D6D6D',
                            'Lentisphaeria' = '#424242', 'Omnitrophia' = '#6D6D6D','Chlamydiae' = '#9B9B9B',
                            'Planctomycetes' = '#69267E', 'Phycisphaerae' = '#BE8DCB','Blastocatellia' = '#00ADFF',
                            'Holophagae' = '#86A4C6', 'Vicinamibacteria' = '#002671',
                            'Acidobacteriae' = '#1F78B4', 'Thermoanaerobaculia' = '#0051BF',
                            'Bdellovibrionia' = '#8C789D','Oligoflexia' = '#654584',
                            'Bacilli' = '#637B88',  'Clostridia' = '#384F5B',
                            'Negativicutes' = '#91AAB8', 'Syntrophomonadia' = '#C2DCEA',
                            'Desulfitobacteriia' =  '#0E2732', 'Thermoanaerobacteria' = '#005474',
                            'Myxococcia' = '#2E5A51','Polyangia' = '#003029',
                            'Nitrospinia' = '#e3a6ce', 'Campylobacteria' = '#002960',
                            'Deinococci' = '#ba9864','Fusobacteriia' = '#fb9a99',
                            'Desulfovibrionia' = '#005c69', 'Desulfobulbia' = '#3D518E',
                            'Desulfuromonadia' = '#6D7DBE', 'Desulfobacteria' = '#000036', 'Methylomirabilia' = '#000000',
                            'Gemmatimonadetes' = '#c55e5c',
                            'Anaerolineae' = '#00d198', 'Chloroflexia' = '#009F6A',
                            'Spirochaetia' = '#5CB1D6', 'Leptospirae' = '#005576', 
                            'Calditrichia'  = '#8a007a', 'Halobacteria' = '#b79c64', 'Nitrospiria' = '#41815f',
                            'Babeliae' = '#5b95e5', 'Saccharimonadia' = '#33af9c', 'Cloacimonadia' = '#fbed5c',
                            'Synergistia' = '#ce7800', 'Abditibacteria' = '#87878b', 'Deferribacteres' = '#4dbaa9')

##information from my dataset----
rem_fc %>% 
  summarize_phyloseq()

##check that all samples have their FC counts
rem_fc@sam_data %>% 
  head()

##extract some information from our dataset
rem_fc %>% 
  ntaxa() ##5670 asv silva // gtbd 4599
rem_fc %>% 
  nsamples() ##312 samples

##subseting by library size
nsamples(rem_fc) #312
rem_fc_filt <- prune_samples(sample_sums(rem_fc) > 5000, rem_fc)
nsamples(rem_fc_filt) #308 there were 4 samples with less than 5000 reads that were discarded

##filter all asv that sum 0
rem_fc_filt <- prune_taxa(taxa_sums(rem_fc_filt) >0, rem_fc_filt)
rem_fc_filt#4594 ASVs

#we add a column with asv number to plot them easily
rem_fc_filt@tax_table <- rem_fc_filt@tax_table %>% 
  mutate_tax_table(asv_num = str_c( 'asv' , 1:ncol(rem_fc_filt@otu_table)))

# ----- 1. RELATIVE ABUNDANCES FOR EACH ASV ################     ----------
rem_fc_relabun <- phyloseq::transform_sample_counts(rem_fc_filt, 
                                                    function(x)
                                                    {x / sum(x)})


# ----- 2. CALCULATE PSEUDOABUNDANCES  ################     ----------
##we melt the data to create a tidy dataset
rem_relabun_melt <- psmelt(rem_fc_relabun) %>%
  as_tibble

#write.table(rem_relabun_melt, "data/rem_relabun_melt.txt", sep= "\t")

## Pseudoabundances calculation
pseudoabund_df_fc <- rem_relabun_melt %>% 
  pseudoabun_fc(grp = c("asv_num", "treatment", "season", "time", "replicate", "hours.x"), 
                abundance = Abundance, fc = All)

pseudoabund_df_wide_fc <- pseudoabund_df_fc %>% ## data transformed to wide format
  pivot_wider(names_from = c("asv_num"), values_from = "pseudoabundance", 
              id_cols = c("treatment", "season", "time", "replicate", "hours.x"))

#write.table(pseudoabund_df_wide_fc, "data/intermediate_files/regressions_test_datasets/psueodabund_df_fc_wide_silva.txt", sep = "\t")


# ----- 3. MULTIPLE LINEAR REGRESSIONS TO CALCULATE THE SLOPE (I.E. GROWTH) FOR EACH ASV AT EACH CONDITION (SEASON-TREATMENT) ################     ----------

## 3.1 Prepare  datasets for calculating regressions
 
env <- rem_fc_filt@sam_data #metadata
env$treatment <- env$treatment %>% 
  factor(levels=(c("CL", "CD", "PL", "PD", "DL", "VL")))
env$season <- env$season %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Early_fall", "Fall")))

pseudoabund_df_wide_m <- pseudoabund_df_wide_fc %>% #asv_tab & metadata
  merge(env, by = c("treatment", "season", "time", "replicate", "hours.x"))

data <- pseudoabund_df_wide_m %>% #pseudoabundances dataset form fc data in wide format with metadata.
  mutate(t_rep = paste(time, replicate, sep="_")) #create a new variable combining replicate and time id.

env <- env %>% #create a new variable combining replicate and time id. for environmental data
  as_tibble() %>%
  mutate(t_rep = paste(time, replicate, sep="_"))

## 3.2 Filter for each treatment, season and time
asv_filt_t0234 <- filter.season.treatment.time(data = data,  
                                               treatment_col = treatment, treatment_keep = "VL", 
                                               season_col = season, season_keep = "Winter",
                                               time_col = time, time_keep = c("t0", "t2", "t3", "t4")) %>%
  select(starts_with("asv"), #we clean metadata columns that we don't need now and keep only the column t_rep
         matches(c("t_rep", "hours.x"))) %>%
  mutate(across(everything(), ~replace(., . == 0, "1"))) %>% ##0 values will be substituted by 1 because it is mandatory for transforming pseudoabundances to log
  mutate(across(!t_rep, as.numeric))

## 3.3 Calculating regressions (CHANGE MATRIX DIM WITH YOUR MATRIX DIM, IN THIS CASE: 4594)
regress <- apply(X = asv_filt_t0234[,c(1:4594)], MARGIN = 2, 
                 FUN = multiple.linear.regressions, env = asv_filt_t0234$hours.x) ##DIMENSIONS OF THE MATRIX SHOULD BE ADAPTED TO DATA!!!!!
regress <- as.data.frame(t(regress)) %>% 
  rownames_to_column(var = "asv_num")

asv_m_reg_Winter_VL_t0234 <- filter.season.treatment.time(data = env,  
                                                          treatment_col = treatment, treatment_keep = "VL", 
                                                          season_col = season, season_keep = "Winter",
                                                          time_col = time, time_keep = c("t0", "t2", "t3", "t4")) %>%
  select(starts_with("asv"), ".sample", "treatment", "replicate",  "time", "season", "sample_name", "sample_code", "selected_file_name",
         "reads", "light_regime.y", "hours.y", "LNA", "HNA", "All", "BM", "Leu.PB", "SGR", "TD", "t_rep")  %>%
  left_join(asv_filt_t023, by = "t_rep") %>%
  pivot_longer(cols = starts_with("asv"), names_to = "asv_num", values_to = "pseudoabundance") %>%
  left_join(regress, by = "asv_num") %>%
  add_column(regression_times = "t0234")

## GO TO THE STEP 3.2 AND CHANGE SEASON AND TIME UNTIL ALL REGRESSIONS ARE CALCULATED.
## calculate regressions with 3 points and 4 points.

##create objects for ech conditions
reg_winter <- bind_rows(asv_m_reg_Winter_CD_t0234, 
                        asv_m_reg_Winter_CL_t0234,
                        asv_m_reg_Winter_PL_t0234,
                        asv_m_reg_Winter_PD_t0234, 
                        asv_m_reg_Winter_DL_t0234,
                        asv_m_reg_Winter_VL_t0234)

reg_spring <- bind_rows(asv_m_reg_Spring_CD_t0234, 
                        asv_m_reg_Spring_CL_t0234,
                        asv_m_reg_Spring_PL_t0234,
                        asv_m_reg_Spring_PD_t0234, 
                        asv_m_reg_Spring_DL_t0234,
                        asv_m_reg_Spring_VL_t0234)

reg_summer <- bind_rows(asv_m_reg_Summer_CD_t0234, 
                        asv_m_reg_Summer_CL_t0234,
                        asv_m_reg_Summer_PL_t0234,
                        asv_m_reg_Summer_PD_t0234, 
                        asv_m_reg_Summer_DL_t0234,
                        asv_m_reg_Summer_VL_t0234)

reg_early_fall <- bind_rows(asv_m_reg_Early_fall_CD_t0234, 
                            asv_m_reg_Early_fall_CL_t0234,
                            asv_m_reg_Early_fall_PL_t0234,
                            asv_m_reg_Early_fall_PD_t0234, 
                            asv_m_reg_Early_fall_DL_t0234,
                            asv_m_reg_Early_fall_VL_t0234)

reg_fall <- bind_rows(asv_m_reg_Fall_CD_t0234, 
                      asv_m_reg_Fall_CL_t0234,
                      asv_m_reg_Fall_PL_t0234,
                      asv_m_reg_Fall_PD_t0234, 
                      asv_m_reg_Fall_DL_t0234,
                      asv_m_reg_Fall_VL_t0234)

## 4 times regressions
reg_all_t0234_silva <- bind_rows(reg_winter, 
                                 reg_spring,
                                 reg_summer,
                                 reg_early_fall,
                                 reg_fall)
## 3 times regressions
reg_all_t023_silva <- bind_rows(reg_winter, 
                                 reg_spring,
                                 reg_summer,
                                 reg_early_fall,
                                 reg_fall)

#save intermediate files
#write.csv(reg_all_t0234_silva, "data/intermediate_files/asv_reg_all_t0234_silva.csv")
#regressions 3 punts
#write.table(reg_all_t0234_silva, "data/intermediate_files/asv_reg_all_t0234_silva.txt", sep = "\t")


# ----- 4. COMPARISON BETWEEN 3 TIMES AND 4 TIMES SLOPES ################     ----------------

## import data 
# reg_all_t0234_silva <-  read.table("data/intermediate_files/asv_reg_all_t0234_silva.txt", sep="\t", header = T) %>%
#   as_tibble()
# reg_all_t023_silva <-  read.table("data/intermediate_files/asv_reg_all_t023_silva.txt", sep="\t", header = T) %>%
#   as_tibble()

reg_all_slopes_chosen_silva <- comparing.reg(df1 = reg_all_t023_silva, 
                                             df2 = reg_all_t0234_silva,
                                             treatment = treatment,
                                             season = season,
                                             asv_num = asv_num, 
                                             slope = slope, 
                                             slope.pval = slope.pval)

tax_table <- rem_fc_filt@tax_table %>%
  as_tibble()

reg_all_slopes_chosen_silva_tax <- reg_all_slopes_chosen_silva %>%
  left_join(tax_table, by = "asv_num", copy = TRUE) ##add taxonomy

#write.csv(reg_all_slopes_chosen_silva_tax, "data/intermediate_files/reg_all_slopes_chosen_silva_tax_corrected_asv_num.csv")

# ----- 5. VIOLIN PLOTS SINGLE-ASV BASED GROWTH RATES AND ABUNDANCE-BASED GROWTH RATES (FIGURE 1) ----------------

##filter only keep significant and positive slopes
reg_all_slopes_chosen_silva_tax <-  reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  filter(season != "Early_fall")

#filter by abundant_asv (at some point of the experiment >1% of relative abundance)
abundant_asv <- rem_relabun_melt %>% 
  filter(Abundance > 0.01) %>% #more than 1% of the community at some point
  dplyr::select(asv_num) %>%
  unique() %>%
  as_tibble()

reg_all_slopes_chosen_silva_tax_1perc <- reg_all_slopes_chosen_silva_tax %>% 
  right_join(abundant_asv, by = "asv_num", copy = TRUE) %>%
  dplyr::filter(slope_chosen != is.na(slope_chosen)) ##remove asv_num that are abundant but don't have significant slope

##transform slopes /hours to /days
reg_all_slopes_chosen_silva_tax_1perc <- reg_all_slopes_chosen_silva_tax_1perc %>%
  mutate(slope_chosen_days = slope_chosen*24)

##reorder treatments and seasons for plots
reg_all_slopes_chosen_silva_tax_1perc$treatment <- reg_all_slopes_chosen_silva_tax_1perc$treatment %>% 
  factor(levels=(c("CL", "CD", "PL", "PD", "DL", "VL")))
reg_all_slopes_chosen_silva_tax_1perc$season <- reg_all_slopes_chosen_silva_tax_1perc$season %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))

#statistics seasons----
anova<-aov(slope_chosen_days ~season, data = reg_all_slopes_chosen_silva_tax_1perc)
summary(anova) #p<0.05 ***
TukeyHSD(anova)#para ver los grupos que son significativamente distintos
aov_residuals <- residuals(object = anova)
# plot(anova, 1)
# plot(anova, 2)
shapiro.test(x = aov_residuals )#From the output, the p-value > 0.05 implying that the distribution of the data are 
#not significantly different from normal distribution. In other words, we can assume the normality.
##p-value < 2.2e-16 NO normal 
dunnTest(slope_chosen_days ~ season, data = reg_all_slopes_chosen_silva_tax_1perc,
         method= 'bonferroni')
results<-dunnTest(slope_chosen_days ~ season, data = reg_all_slopes_chosen_silva_tax_1perc,
                  method='bonferroni')$res
results<-results[results$P.adj<0.05,]
X = results$P.adj <= 0.05
names(X) = gsub(" ",  "",  results$Comparison)
multcompLetters(X)
# Fall Spring Summer Winter 
# "a"    "a"    "b"    "b" 

#statistical groups for non-parammetric test
letters_seas <- data.frame(multcompLetters(X)['Letters'])

## plot seasons violin plot----
seas <- reg_all_slopes_chosen_silva_tax_1perc %>%
  ggplot(aes(season, slope_chosen_days))+ 
  geom_point(aes(color = treatment), alpha = 0.7, position = position_jitter(0.25))+
  geom_violin(alpha = 0.2, draw_quantiles = c(0.25, 0.75) )+
  labs(color = "Treatment")+
  stat_summary(fun = "mean", geom = "crossbar", 
               width = 0.6, colour = "black")+
  geom_text(data = letters_seas, aes(y = 9, x = row.names(letters_seas), label = Letters),
            position = position_nudge(x = 0.2), hjust = 0.7, color = "black")+
  labs(x= "Season", y = expression("Growth rate day"^"-1"))+  
  scale_y_continuous(expand = c(0, 0), limits = c(0,10.4))+
  scale_color_manual(values = palette_treatments_remei)+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

## statistics treatments----
anova<-aov(slope_chosen_days ~ treatment, data = reg_all_slopes_chosen_silva_tax_1perc)
summary(anova) #p<0.05 ***
TukeyHSD(anova)
aov_residuals <- residuals(object = anova)
# plot(anova, 1)
# plot(anova, 2)
shapiro.test(x = aov_residuals )#From the output, the p-value > 0.05 implying that the distribution of the data are #not significantly different from normal distribution. In other words, we can assume the normality.
##p-value < 2.2e-16 NO és normal 
#if not normal: #non parametric test
kruskal.test(slope_chosen_days ~ treatment, data = reg_all_slopes_chosen_silva_tax_1perc) #
#si sale p<0.05 hago dunn test para ver cuales son significativamente distintos
dunnTest(slope_chosen_days ~ treatment, data = reg_all_slopes_chosen_silva_tax_1perc, ##test the d
         method= 'bonferroni')
results<-dunnTest(slope_chosen_days ~ treatment, data = reg_all_slopes_chosen_silva_tax_1perc,
                  method='bonferroni')$res
#results<-results[results$P.adj<0.05,]
X = results$P.adj <= 0.05
names(X) = gsub(" ",  "",  results$Comparison)
multcompLetters(X)
# CD   CL   DL   PD   PL   VL 
# "a"  "a" "bc" "bd"  "d"  "c" 
letters_treat <- data.frame(multcompLetters(X)['Letters'])

## plot treatments violin plot----
treat <- reg_all_slopes_chosen_silva_tax_1perc %>%
  ggplot(aes(treatment, slope_chosen_days))+ #
  geom_point(aes(treatment, slope_chosen_days, color = season), alpha = 0.7, position = position_jitter(0.25))+ 
  geom_violin(alpha = 0.2, draw_quantiles = c(0.25, 0.75), )+
  labs(color = "Season")+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.6,
               colour = "black")+
  geom_text(data = letters_treat, aes(y = 9, x = row.names(letters_treat), label = Letters),
            position = position_nudge(x = 0.2), hjust = 0.7, color = "black")+
  labs(x= 'Season', y = expression("Growth rate day"^"-1"), color = "Season")+
  scale_y_continuous(expand = c(0, 0), limits = c(0,10.4))+ 
  scale_color_manual(values = palette_seasons_4)+
  theme_bw()+
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

## Abundance-based growth rates ----

### Calculated by Sánchez O, Ferrera I, Mabrito I et al. Seasonal impact of grazing, viral mortality, resource availability and light on the group-specific growth rates of coastal Mediterranean bacterioplankton. Sci Rep 2020;10, DOI: 10.1038/s41598-020-76590-5.

GR_dapis_OS <- read.delim2("data/envdata/GR_DAPIS_OS/GR_REMEI_DAPIS_OS_Ed.csv", sep = ";") %>%
  as_tibble()

GR_dapis_OS_filt <- GR_dapis_OS %>%
  filter(treatment %in% c("CL", "CD", "PL", "PD",  "DL", "VL")) %>%
  mutate(PRO = as.numeric(PRO)) ##total abundance prokaryotes

#### reorder treatments and seasons for plots
GR_dapis_OS_filt$treatment <- GR_dapis_OS_filt$treatment %>% 
  factor(levels=(c("CL", "CD", "PL", "PD", "DL", "VL")))
GR_dapis_OS_filt$season <- GR_dapis_OS_filt$season %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))

### statistics seasons------
anova<-aov(PRO~season, data=GR_dapis_OS_filt)
summary(anova) #p<0.05 ***
tukey <- anova %>%
  TukeyHSD()
aov_residuals <- residuals(object = anova)
plot(anova, 1)
plot(anova, 2)
shapiro.test(x = aov_residuals )#From the output, the p-value > 0.05 implying that the distribution of the data are #not significantly different from normal distribution. In other words, we can assume the normality.
##p-value = 0.218 Normal 
#plot(TukeyHSD(anova, conf.level=.95), las = 2)

#Groups
cld <- multcompLetters4(anova, tukey)
# table with factors and 3rd quantile
dt <- GR_dapis_OS_filt %>%
  group_by(season) %>%
  summarise(w=mean(PRO), sd = sd(PRO)) %>%
  arrange(desc(w))

# extracting the compact letter display and adding to the Tk table
cld <- as.data.frame.list(cld$season)
dt$cld <- cld$Letters
print(dt)
# season     w     sd cld  
# <fct>  <dbl>  <dbl> <chr>
#   1 Summer 0.579 0.224  a    
# 2 Winter 0.438 0.296  a    
# 3 Fall   0.208 0.156  b    
# 4 Spring 0.172 0.0815 b 

gr_dapis_seas <- 
  GR_dapis_OS_filt %>%
  ggplot(aes(season, PRO))+ #
  geom_point(aes(color = treatment), alpha = 0.7, position = position_jitter(0.2))+
  geom_violin(alpha = 0.1, draw_quantiles = c(0.25, 0.75))+
  geom_text(data = dt, aes(y = 1.5, x = season, label = str_trim(cld)),
            position = position_nudge(x = 0), hjust = 0, color = "black")+
  stat_summary(fun = "mean", geom = "crossbar", 
               width = 0.3, colour = "black")+
  labs(x= "Season", y = expression("Growth rate day"^"-1"), color = "Treatment")+ #(μ)
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2.5))+
  scale_color_manual(values = palette_treatments_remei)+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0), legend.position = "right",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

##Treatments
anova<-aov(PRO~treatment, data = GR_dapis_OS_filt)
summary(anova) #p<0.05 ***
tuk<- TukeyHSD(anova)#para ver los grupos que son significativamente distintos Summer i Winter no son distintos y Fall i Spring tampoco
aov_residuals <- residuals(object = anova)
plot(anova, 1)
plot(anova, 2)
shapiro.test(x = aov_residuals )#From the output, the p-value > 0.05 implying that the distribution of the data are #not significantly different from normal distribution. In other words, we can assume the normality.
##p-value = 0.003229 NO és normal 
#si no es normal:
#test no parametrico
kruskal.test(PRO ~ season, data = GR_dapis_OS_filt) #
#si sale p<0.05 hago dunn test para ver cuales son significativamente distintos
#p-value = 3.598e-06
library(FSA)##test the d
duntest <- dunnTest(PRO~treatment, data = GR_dapis_OS_filt,
                    method= 'bonferroni')
results<-dunnTest(PRO~treatment, data = GR_dapis_OS_filt,
                  method='bonferroni')$res
#results<-results[results$P.adj<0.05,]
##Diferents: CD-DL, CL-DL, CL-PL, CD-VL, CL-VL
#statistical groups for non-parammetric test
X = results$P.adj <= 0.05

names(X) = gsub(" ",  "",  results$Comparison)
multcompLetters(X)
letters_dapis_treat <- data.frame(multcompLetters(X)['Letters'])

#    CD    CL    DL    PD    PL    VL 
# "ab"   "a"   "c" "abc"  "bc"   "c" 

gr_dapis_treat <- 
  GR_dapis_OS_filt %>%
  ggplot(aes(treatment, PRO))+ #
  geom_point(aes(color = season), alpha = 0.7, position = position_jitter(0.2))+
  geom_violin(alpha = 0.1, draw_quantiles = c(0.25, 0.75))+
  geom_text(data = letters_dapis_treat, aes(y = 1.5, x = row.names(letters_dapis_treat), label = Letters),
            position = position_nudge(x = 0), hjust = 0, color = "black")+
  stat_summary(fun = "mean", geom = "crossbar", 
               width = 0.3, colour = "black")+
  labs(x= "Treatment", y = expression("Growth rate day"^"-1"), color = "Season")+ #(μ)
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2.5))+
  scale_color_manual(values = palette_seasons_4)+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0), legend.position = "right",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())


## Pannel construction figure 1 ----
gr_asv_dapis <- grid.arrange(treat, gr_dapis_treat, seas, gr_dapis_seas)
# 
# ggsave('GR_asv_dapis_stat.pdf', gr_asv_dapis, path = "~/Documentos/Doctorat/REMEI/results/figures/corrected_asv_num/",
#        height = 175,
#        width = 250,
#        units = 'mm')


# 5. GROWTH RATES DISTRIBUTION PATTERNS AT DIFFERENT TAXONOMIC RANKS -----
reg_all_slopes_chosen_silva_tax$treatment <- reg_all_slopes_chosen_silva_tax$treatment %>% 
factor(levels=(c("CL", "CD", "PL", "PD", "DL", "VL")))
reg_all_slopes_chosen_silva_tax$season <- reg_all_slopes_chosen_silva_tax$season %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))

##transform slopes /hours to /days
reg_all_slopes_chosen_silva_tax <- reg_all_slopes_chosen_silva_tax %>%
  mutate(slope_chosen_days = slope_chosen*24)
