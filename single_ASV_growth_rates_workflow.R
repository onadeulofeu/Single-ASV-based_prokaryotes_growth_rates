# This code is designed to calculate single-ASVs growth rates for the paper entitled 
# Growth rates of marine prokayotes are extremely diverse, even among closely related taxa
# Code developed by Ona Deulofeu-Capo.

#Session information----
sessionInfo()

# R version 4.2.3 (2023-03-15)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Ventura 13.2.1
# 
# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] Biostrings_2.66.0      GenomeInfoDb_1.34.9    XVector_0.38.0         IRanges_2.32.0         S4Vectors_0.36.2      
# [6] BiocGenerics_0.44.0    seqinr_4.2-23          stringi_1.7.12         janitor_2.2.0          ggraph_2.1.0          
# [11] igraph_1.4.1           FSA_0.9.4              ggpmisc_0.5.2          ggpp_0.5.1             multcompView_0.1-8    
# [16] rstatix_0.7.2          ggpubr_0.6.0           ggridges_0.5.4         ggforce_0.4.1          multipanelfigure_2.1.2
# [21] gridExtra_2.3          microbiome_1.20.0      scales_1.2.1           vegan_2.6-4            lattice_0.20-45       
# [26] permute_0.9-7          speedyseq_0.5.3.9018   reshape2_1.4.4         magrittr_2.0.3         lubridate_1.9.2       
# [31] forcats_1.0.0          stringr_1.5.0          dplyr_1.1.1            purrr_1.0.1            readr_2.1.4           
# [36] tidyr_1.3.0            tibble_3.2.1           ggplot2_3.4.1          tidyverse_2.0.0        phyloseq_1.42.0       
# 
# loaded via a namespace (and not attached):
#   [1] Rtsne_0.16                 assertive.base_0.0-9       colorspace_2.1-0           ggsignif_0.6.4            
# [5] snakecase_0.11.0           rstudioapi_0.14            farver_2.1.1               graphlayouts_0.8.4        
# [9] MatrixModels_0.5-1         ggrepel_0.9.3              fansi_1.0.4                codetools_0.2-19          
# [13] splines_4.2.3              confintr_0.2.0             polyclip_1.10-4            dunn.test_1.3.5           
# [17] ade4_1.7-22                polynom_1.4-1              jsonlite_1.8.4             broom_1.0.4               
# [21] cluster_2.1.4              BiocManager_1.30.20        compiler_4.2.3             backports_1.4.1           
# [25] Matrix_1.5-3               cli_3.6.1                  tweenr_2.0.2               quantreg_5.94             
# [29] tools_4.2.3                gtable_0.3.3               glue_1.6.2                 GenomeInfoDbData_1.2.9    
# [33] Rcpp_1.0.10                carData_3.0-5              Biobase_2.58.0             cellranger_1.1.0          
# [37] vctrs_0.6.1                rhdf5filters_1.10.1        multtest_2.54.0            ape_5.7-1                 
# [41] nlme_3.1-162               iterators_1.0.14           assertive.files_0.0-2      timechange_0.2.0          
# [45] lifecycle_1.0.3            zlibbioc_1.44.0            MASS_7.3-58.3              tidygraph_1.2.3           
# [49] hms_1.1.3                  parallel_4.2.3             biomformat_1.26.0          rhdf5_2.42.0              
# [53] SparseM_1.81               foreach_1.5.2              boot_1.3-28.1              rlang_1.1.0               
# [57] pkgconfig_2.0.3            bitops_1.0-7               Rhdf5lib_1.20.0            labeling_0.4.2            
# [61] assertive.properties_0.0-5 cowplot_1.1.1              tidyselect_1.2.0           plyr_1.8.8                
# [65] R6_2.5.1                   magick_2.7.4               generics_0.1.3             pillar_1.9.0              
# [69] withr_2.5.0                mgcv_1.8-42                assertive.numbers_0.0-2    survival_3.5-5            
# [73] abind_1.4-5                RCurl_1.98-1.12            crayon_1.5.2               car_3.1-2                 
# [77] assertive.types_0.0-3      utf8_1.2.3                 tzdb_0.3.0                 viridis_0.6.2             
# [81] grid_4.2.3                 readxl_1.4.2               data.table_1.14.8          digest_0.6.31             
# [85] gridGraphics_0.5-1         munsell_0.5.0              viridisLite_0.4.1   


##packages used ----
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
library(lubridate)
library(igraph)
library(ggraph)
library(janitor)
library(stringi)
library(seqinr)
library(Biostrings)
library(Hmisc)
library(readxl)

#define directory ----
setwd ("~/Documentos/Doctorat/REMEI/")

# import functions ----
source("src/pseudoabundances_calculation_fc.R") ## for calculating pseudoabundances
source("src/filter.season.treatment.time.R") ## filter for calculating regressions
source("src/multiple.linear.regressions.R") ## multiple linear regressions
source("src/comparing.reg.R") ## comparisons between 3 times and 4 times regressions
source('src/growth.rates.distr.tax.ranks.ridges.divided.R') ##distribution plots
source('src/growth.rates.distr.tax.ranks.ridges.R') ##distribution plots divided by classes
source("src/calculate_common_asv_conditions.R") ##calculate common responsive ASVs between treatments and seasons
source("src/calculate_exclusive_asv_condition.R") ##calculate unique responsive ASVs at each treatment and season


# import data ----
rem_fc <- readRDS("./data/intermediate_files/remei_phyloseq_silva_fc.rds") #phyloseq object with otu_table, sample data and taxonomu
#sample data must have total abundance data to be able to calculate pseudoabundances

#palettes used (colors and shapes) ----
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

shapes_treatments <- c('CL' = 16, 'CD' = 1, 'PL' = 2, 'PD' = 17, 'DL' = 15, 'VL' = 18)

## Dataset information ----
rem_fc |> 
  summarize_phyloseq()

##check that all samples have their FC counts
rem_fc@sam_data |> 
  head()

##extract some information from our dataset
rem_fc |> 
  ntaxa() ##5670 asv silva // gtbd 4599
rem_fc |> 
  nsamples() ##312 samples

##subseting by library size
nsamples(rem_fc) #312
rem_fc_filt <- prune_samples(sample_sums(rem_fc) > 5000, rem_fc)
nsamples(rem_fc_filt) #308 there were 4 samples with less than 5000 reads that were discarded

##filter all asv that sum 0
rem_fc_filt <- prune_taxa(taxa_sums(rem_fc_filt) >0, rem_fc_filt)
rem_fc_filt#4594 ASVs

#we add a column with asv number to plot them easily
rem_fc_filt@tax_table <- rem_fc_filt@tax_table |> 
  mutate_tax_table(asv_num = str_c( 'asv' , 1:ncol(rem_fc_filt@otu_table)))


# ----- 1. RELATIVE ABUNDANCES FOR EACH ASV ################     ----------
rem_fc_relabun <- phyloseq::transform_sample_counts(rem_fc_filt, 
                                                    function(x)
                                                    {x / sum(x)})



# ----- 2. CALCULATE PSEUDOABUNDANCES  ################     ----------
##we melt the data to create a tidy dataset
rem_relabun_melt <- psmelt(rem_fc_relabun) |>
  as_tibble () |>
  dplyr::filter(season != 'Early_fall')

#write.table(rem_relabun_melt, "data/rem_relabun_melt.txt", sep= "\t")

## Pseudoabundances calculation
pseudoabund_df_fc <- rem_relabun_melt |> 
  pseudoabun_fc(grp = c("asv_num", "treatment", "season", "time", "replicate", "hours.x"), 
                abundance = Abundance, fc = All)

pseudoabund_df_wide_fc <- pseudoabund_df_fc |> ## data transformed to wide format
  pivot_wider(names_from = c("asv_num"), values_from = "pseudoabundance", 
              id_cols = c("treatment", "season", "time", "replicate", "hours.x"))

#write.table(pseudoabund_df_wide_fc, "data/intermediate_files/regressions_test_datasets/psueodabund_df_fc_wide_silva.txt", sep = "\t")



# ----- 3. MULTIPLE LINEAR REGRESSIONS TO CALCULATE THE SLOPE (I.E. GROWTH) FOR EACH ASV AT EACH CONDITION (SEASON-TREATMENT)################     ----------

## 3.1 Prepare  datasets for calculating regressions
 
env <- rem_fc_filt@sam_data #metadata
env$treatment <- env$treatment |> 
  factor(levels=(c("CL", "CD", "PL", "PD", "DL", "VL")))
env$season <- env$season |> 
  factor(levels=(c("Winter", "Spring", "Summer", "Early_fall", "Fall")))

pseudoabund_df_wide_m <- pseudoabund_df_wide_fc |> #asv_tab & metadata
  merge(env, by = c("treatment", "season", "time", "replicate", "hours.x"))

data <- pseudoabund_df_wide_m |> #pseudoabundances dataset form fc data in wide format with metadata.
  mutate(t_rep = paste(time, replicate, sep="_")) #create a new variable combining replicate and time id.

env <- env |> #create a new variable combining replicate and time id. for environmental data
  as_tibble() |>
  mutate(t_rep = paste(time, replicate, sep="_"))

## 3.2 Filter for each treatment, season and time
asv_filt_t0234 <- filter.season.treatment.time(data = data,  
                                               treatment_col = treatment, treatment_keep = "VL", 
                                               season_col = season, season_keep = "Winter",
                                               time_col = time, time_keep = c("t0", "t2", "t3", "t4")) |>
  select(starts_with("asv"), #we clean metadata columns that we don't need now and keep only the column t_rep
         matches(c("t_rep", "hours.x"))) |>
  mutate(across(everything(), ~replace(., . == 0, "1"))) |> ##0 values will be substituted by 1 because it is mandatory for transforming pseudoabundances to log
  mutate(across(!t_rep, as.numeric))

## 3.3 Calculating regressions (CHANGE MATRIX DIM WITH YOUR MATRIX DIM, IN THIS CASE: 4594)
regress <- apply(X = asv_filt_t0234[,c(1:4594)], MARGIN = 2, 
                 FUN = multiple.linear.regressions, env = asv_filt_t0234$hours.x) ##DIMENSIONS OF THE MATRIX SHOULD BE ADAPTED TO DATA!!!!!
regress <- as.data.frame(t(regress)) |> 
  rownames_to_column(var = "asv_num")

asv_m_reg_Winter_VL_t0234 <- filter.season.treatment.time(data = env,  
                                                          treatment_col = treatment, treatment_keep = "VL", 
                                                          season_col = season, season_keep = "Winter",
                                                          time_col = time, time_keep = c("t0", "t2", "t3", "t4")) |>
  select(starts_with("asv"), ".sample", "treatment", "replicate",  "time", "season", "sample_name", "sample_code", "selected_file_name",
         "reads", "light_regime.y", "hours.y", "LNA", "HNA", "All", "BM", "Leu.PB", "SGR", "TD", "t_rep")  |>
  left_join(asv_filt_t023, by = "t_rep") |>
  pivot_longer(cols = starts_with("asv"), names_to = "asv_num", values_to = "pseudoabundance") |>
  left_join(regress, by = "asv_num") |>
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
#write.table(reg_all_t0234_silva_ed, "data/intermediate_files/asv_reg_all_t0234_silva_v1.txt", sep = "\t")


# ----- 4. COMPARISON BETWEEN 3 TIMES AND 4 TIMES SLOPES ################     ----------------

## import data 
reg_all_t0234_silva <-  read.table("data/intermediate_files/asv_reg_all_t0234_silva_v1.txt", sep="\t", header = T) |>
  as_tibble()
reg_all_t023_silva <-  read.table("data/intermediate_files/asv_reg_all_t023_silva.txt", sep="\t", header = T) |>
  as_tibble()

reg_all_slopes_chosen_silva <- comparing.reg(df1 = reg_all_t023_silva, 
                                             df2 = reg_all_t0234_silva,
                                             treatment = treatment,
                                             season = season,
                                             asv_num = asv_num, 
                                             slope = slope, 
                                             slope.pval = slope.pval)

reg_all_slopes_chosen_silva_tax <- reg_all_slopes_chosen_silva |>
  left_join(tax_table, by = "asv_num", copy = TRUE) ##add taxonomy

#write.csv(reg_all_slopes_chosen_silva_tax, "data/intermediate_files/reg_all_slopes_chosen_silva_tax_corrected_asv_num.csv")
#write.csv(reg_all_slopes_chosen_silva_tax, "data/intermediate_files/reg_all_slopes_chosen_silva_tax_corrected_asv_num_v1.csv")
reg_all_slopes_chosen_silva_tax <- read.delim2('data/intermediate_files/reg_all_slopes_chosen_silva_tax_corrected_asv_num_v1.csv', sep = ',')|>
  as_tibble() |>
  dplyr::mutate(slope_chosen_days = as.numeric(slope_chosen_days),
                pvalue_slope_chosen = as.numeric(pvalue_slope_chosen))
  

## growth rates example calculation (figure 1) ----
pseudoabund_df_wide_fc <- read.table("data/intermediate_files/regressions_test_datasets/psueodabund_df_fc_wide_silva.txt", sep = "\t")

pseudoabund_df_wide_fc$season <-  pseudoabund_df_wide_fc$season |> factor( 
  levels = (c('Winter', 'Spring', 'Summer', 'Fall', 'Early_fall')))
pseudoabund_df_wide_fc$treatment <- pseudoabund_df_wide_fc$treatment |>
  factor(levels=(c("CL", "CD", "PL", "PD", "DL", "VL")))

## 4 times
t4_exemple <- 
  pseudoabund_df_wide_fc |>
  dplyr::select(1:5, 'asv2') |>
  subset(season != "Early_fall" & 
           treatment == 'DL') |>
  pivot_longer(cols = starts_with('asv'), names_to = 'asv_num', values_to = 'pseudoabund') |>
  left_join(reg_all_slopes_chosen_silva_tax, by = 'asv_num', relationship = "many-to-many") |>
  ggplot(aes(hours.x, log(pseudoabund), color = season.x, shape = treatment.x))+
  geom_point(size = 1)+
  stat_poly_line(aes(group = season.x), fm.values = TRUE)+
  stat_poly_eq(aes(
    label =  paste(after_stat(p.value.label))), 
    formula = x ~ log(y),
    method = stats::lm,
    p.digits = 2, 
    coef.keep.zeros = T, #npcx = 1, npcy =1,
    na.rm = FALSE)+
  scale_y_continuous(labels = scales::scientific, breaks = c(5e+00, 1e+01), limits = c(2e+00, 1.25e+01))+
  labs(y = 'ln(ASV2 pseudoabundance)', x = 'Time (h)')+
  scale_color_manual(values = palette_seasons_4)+
  scale_shape_manual(values = shapes_treatments)+
  theme_bw()+
  theme(strip.text = element_text(size = 6), axis.text.x = element_text(angle = 0, size= 6), legend.position = "none",
        strip.background = element_blank(), strip.placement = 'outside',
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin=unit(c(0.25,0.5,0.25,0.15),"cm"), axis.title = element_text(size=7), axis.text.y = element_text(size = 6), text = element_text(size = 6)) 

## 3 times
t3_exemple <- 
  pseudoabund_df_wide_fc |>
  dplyr::select(1:5, 'asv2') |>
  subset(season != "Early_fall" & 
           treatment == 'DL' &
           time != 't4')|>
  pivot_longer(cols = starts_with('asv'), names_to = 'asv_num', values_to = 'pseudoabund') |>
  left_join(reg_all_slopes_chosen_silva_tax, by = 'asv_num', relationship = "many-to-many") |>
  ggplot(aes(hours.x, log(pseudoabund), color = season.x, shape = treatment.x))+
  geom_point(size = 1)+
  scale_y_continuous(labels = scales::scientific, breaks = c(5e+00, 1e+01), limits = c(2e+00, 1.25e+01))+
  stat_poly_line(aes(group = season.x), method = 'lm')+
  labs(y = 'ln(ASV2 pseudoabundance)', x = 'Time (h)', color = 'Season', shape = 'Treatment')+
  scale_color_manual(values = palette_seasons_4)+
  scale_x_continuous(limits = c(0,50))+
  stat_poly_eq(aes(
    label =  paste(after_stat(p.value.label))), 
    formula = x ~ log(y),
    method = stats::lm,
    p.digits = 2, 
    coef.keep.zeros = T,
    na.rm = FALSE)+
  #scale_shape_manual(name = 'Treatment', labels = 'DL', values = 12)+
  scale_shape_manual(values = shapes_treatments)+
  guides(shape = 'none', color=guide_legend(override.aes=list(fill=NA)))+
  theme_bw()+
  theme(strip.text = element_blank(), axis.text.x = element_text(angle = 0, size= 7), legend.position = "right",
        axis.title.y = element_blank(), strip.background = element_blank(), axis.text.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        plot.margin=unit(c(0.25,0,0.25,0),"cm"), legend.text = element_text(size = 6), axis.title = element_text(size=7),
        legend.title = element_text(size=7), text = element_text(size = 6)) 


example_gr_calculation <- multi_panel_figure(columns = 2, rows = 1, width = 178, height = 110, 
                                             column_spacing = 0.2, unit = 'mm',
                                             panel_label_type = 'upper-alpha')

example_gr_calculation  %<>%
  fill_panel(t4_exemple, column = 1, row = 1) %<>%
  fill_panel(t3_exemple, column = 2, row = 1)

ggsave('example_gr_calculation3.pdf', example_gr_calculation ,
       path = "~/Documentos/Doctorat/REMEI/results/figures/corrected_asv_num/",
       width = 188,
       height = 120,
       units = 'mm')



# ----- 5. VIOLIN PLOTS SINGLE-ASV BASED GROWTH RATES AND ABUNDANCE-BASED GROWTH RATES ################ ----------------
## import data
reg_all_slopes_chosen_silva_tax <- read.delim2('data/intermediate_files/reg_all_slopes_chosen_silva_tax_corrected_asv_num_v1.csv', sep = ',')|>
  as_tibble() |>
  dplyr::mutate(slope_chosen_days = as.numeric(slope_chosen_days),
                pvalue_slope_chosen = as.numeric(pvalue_slope_chosen),
                slope_chosen = as.numeric(slope_chosen))

##filter only keep significant and positive slopes
reg_all_slopes_chosen_silva_tax <-  reg_all_slopes_chosen_silva_tax |>
  dplyr::filter(pvalue_slope_chosen < 0.05 ) |> #& intercept.pval < value_intercept_pval
  dplyr::filter(slope_chosen > 0) |>
  dplyr::filter(season != "Early_fall") #|>
  # dplyr::distinct(slope_chosen_days, slope_chosen, pvalue_slope_chosen, asv_num, treatment, season, domain, phylum, class, order, family, genus) |>
  # dplyr::mutate(slope_chosen = round(x = slope_chosen, digits = 5),
  #               slope_chosen_days = round(x = slope_chosen_days, digits = 5),
  #               pvalue_slope_chosen = round(x = pvalue_slope_chosen, digits = 5)) |>
  # dplyr::distinct(slope_chosen, pvalue_slope_chosen, slope_chosen_days, asv_num, treatment, season, domain, phylum, class, order, family, genus)

reg_all_slopes_chosen_silva_tax |>
  #group_by(treatment, season) |>
  dplyr::summarize(mean_gr = mean(slope_chosen_days),
                   sd = sd(slope_chosen_days),
                   median = median(slope_chosen_days))

# mean_gr sd median
#  2.56  1.53   2.23

treat_seas_stats <- reg_all_slopes_chosen_silva_tax |>
  group_by(treatment, season) |>
  dplyr::summarize(mean = mean(slope_chosen_days),
                   sd = sd(slope_chosen_days),
                   median = median(slope_chosen_days),
                   cv = sd/mean )

treat_stats <- reg_all_slopes_chosen_silva_tax |>
  group_by(treatment) |>
  dplyr::summarize(mean = mean(slope_chosen_days),
                   sd = sd(slope_chosen_days),
                   median = median(slope_chosen_days),
                   cv = sd/mean )
treat_seas_stats <- reg_all_slopes_chosen_silva_tax |>
  group_by(treatment, season) |>
  dplyr::summarize(mean = mean(slope_chosen_days),
                   sd = sd(slope_chosen_days),
                   median = median(slope_chosen_days),
                   cv = sd/mean )
seas_stats <- reg_all_slopes_chosen_silva_tax |>
  group_by(season) |>
  dplyr::summarize(mean = mean(slope_chosen_days),
                   sd = sd(slope_chosen_days),
                   median = median(slope_chosen_days),
                   cv = sd/mean )

##1 perc
treat_seas_stats <- reg_all_slopes_chosen_silva_tax_1perc |>
  group_by(treatment, season) |>
  dplyr::summarize(mean = mean(slope_chosen_days),
                   sd = sd(slope_chosen_days),
                   median = median(slope_chosen_days),
                   cv = sd/mean )

treat_stats <- reg_all_slopes_chosen_silva_tax_1perc |>
  group_by(treatment) |>
  dplyr::summarize(mean = mean(slope_chosen_days),
                   sd = sd(slope_chosen_days),
                   median = median(slope_chosen_days),
                   cv = sd/mean )

treat_seas_stats <- reg_all_slopes_chosen_silva_tax_1perc |>
  group_by(treatment, season) |>
  dplyr::summarize(mean = mean(slope_chosen_days),
                   sd = sd(slope_chosen_days),
                   median = median(slope_chosen_days),
                   cv = sd/mean )

seas_stats <- reg_all_slopes_chosen_silva_tax_1perc |>
  group_by(season) |>
  dplyr::summarize(mean = mean(slope_chosen_days),
                   sd = sd(slope_chosen_days),
                   median = median(slope_chosen_days),
                   cv = sd/mean )

# growth_rates_table <- reg_all_slopes_chosen_silva_tax_filt |>
#   dplyr::select(treatment, season, asv_num, slope_chosen_days, pvalue_slope_chosen, domain,
#                 phylum, class, order, family, genus)
# 
# write.csv2(growth_rates_table, 'results/tables/growth_rates_table.csv')

#filter by abundant_asv (at some point of the experiment >1% of relative abundance)
abundant_asv <- rem_relabun_melt |> 
  filter(Abundance > 0.01) |> #more than 1% of the community at some point
  dplyr::select(asv_num) |>
  unique() |>
  as_tibble()

reg_all_slopes_chosen_silva_tax_1perc <- reg_all_slopes_chosen_silva_tax |> 
  right_join(abundant_asv, by = "asv_num", copy = TRUE) |>
  dplyr::filter(slope_chosen_days != is.na(slope_chosen_days)) |> ##remove asv_num that are abundant but don't have significant slope
  dplyr::distinct(slope_chosen_days, pvalue_slope_chosen, asv_num, treatment, season, domain, phylum, class, order, family, genus)
  
##transform slopes /hours to /days
# reg_all_slopes_chosen_silva_tax_1perc <- reg_all_slopes_chosen_silva_tax_1perc |>
#   mutate(slope_chosen = as.numeric(slope_chosen)) |>
#   dplyr::mutate(slope_chosen_days = as.numeric(slope_chosen*24))

##reorder treatments and seasons for plots
reg_all_slopes_chosen_silva_tax_1perc$treatment <- reg_all_slopes_chosen_silva_tax_1perc$treatment |> 
  factor(levels=(c("CL", "CD", "PL", "PD", "DL", "VL")))
reg_all_slopes_chosen_silva_tax_1perc$season <- reg_all_slopes_chosen_silva_tax_1perc$season |> 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))

### statistics seasons----
anova<-aov(slope_chosen_days ~season, data = reg_all_slopes_chosen_silva_tax_1perc)
summary(anova) #p<0.05 ***
TukeyHSD(anova)#para ver los grupos que son significativamente distintos
aov_residuals <- residuals(object = anova)
# plot(anova, 1)
# plot(anova, 2)
shapiro.test(x = aov_residuals )#From the output, if the p-value > 0.05 implying that the distribution of the data are 
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
seas <- reg_all_slopes_chosen_silva_tax_1perc |>
  ggplot(aes(season, slope_chosen_days))+ 
  geom_point(aes(color = treatment), alpha = 0.7, position = position_jitter(0.25), size = 1.5)+
  geom_violin(alpha = 0.2, draw_quantiles = c(0.25, 0.75) )+
  labs(color = "Treatment")+
  stat_summary(fun = "mean", geom = "crossbar", 
               width = 0.6, colour = "black")+
  geom_text(data = letters_seas, aes(y = 9, x = row.names(letters_seas), label = Letters),
            position = position_nudge(x = 0.2), hjust = 0.7, color = "black")+
  labs(x= "Season", y = expression("Growth rate (d"^"-1" *')'))+  
  scale_y_continuous(expand = c(0, 0), limits = c(0,10.4))+
  scale_color_manual(values = palette_treatments_remei)+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 6))

### statistics treatments----
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
#if p<0.05 use dunn test to see which are significatively different
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
treat <- reg_all_slopes_chosen_silva_tax_1perc |>
  ggplot(aes(treatment, slope_chosen_days))+ #
  geom_point(aes(treatment, slope_chosen_days, color = season), 
             alpha = 0.7, position = position_jitter(0.25), size = 1.5)+ 
  geom_violin(alpha = 0.2, draw_quantiles = c(0.25, 0.75), )+
  labs(color = "Season")+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.6,
               colour = "black")+
  geom_text(data = letters_treat, aes(y = 9, x = row.names(letters_treat), label = Letters),
            position = position_nudge(x = 0.25), hjust = 0.7, color = "black")+
  labs(x= 'Treatment', y = expression("Growth rate (d"^"-1" *')'), color = "Season")+
  scale_y_continuous(expand = c(0, 0), limits = c(0,10.4))+ 
  scale_color_manual(values = palette_seasons_4)+
  theme_bw()+
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0), 
        legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text = element_text(size = 6))

## Abundance-based growth rates ----

### Calculated by Sánchez O, Ferrera I, Mabrito I et al. Seasonal impact of grazing, viral mortality, resource availability and light on the group-specific growth rates of coastal Mediterranean bacterioplankton. Sci Rep 2020;10, DOI: 10.1038/s41598-020-76590-5.

GR_dapis_OS <- read.delim2("data/envdata/GR_DAPIS_OS/GR_REMEI_DAPIS_OS_Ed.csv", sep = ";") |>
  as_tibble()

GR_dapis_OS_filt <- GR_dapis_OS |>
  filter(treatment %in% c("CL", "CD", "PL", "PD",  "DL", "VL")) |>
  mutate(PRO = as.numeric(PRO)) ##total abundance prokaryotes

#### reorder treatments and seasons for plots
GR_dapis_OS_filt$treatment <- GR_dapis_OS_filt$treatment |> 
  factor(levels=(c("CL", "CD", "PL", "PD", "DL", "VL")))
GR_dapis_OS_filt$season <- GR_dapis_OS_filt$season |> 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))

### statistics seasons------
anova<-aov(PRO~season, data=GR_dapis_OS_filt)
summary(anova) #p<0.05 ***
tukey <- anova |>
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
dt <- GR_dapis_OS_filt |>
  group_by(season) |>
  dplyr::summarise(w=mean(PRO), sd = sd(PRO)) |>
  dplyr::arrange(desc(w))

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
  GR_dapis_OS_filt |>
  ggplot(aes(season, PRO))+ #
  geom_point(aes(color = treatment), alpha = 0.7, position = position_jitter(0.2), size = 1.5)+
  geom_violin(alpha = 0.1, draw_quantiles = c(0.25, 0.75))+
  geom_text(data = dt, aes(y = 1.5, x = season, label = str_trim(cld), size = 1),
            position = position_nudge(x = 0), hjust = 0, color = "black")+
  stat_summary(fun = "mean", geom = "crossbar", 
               width = 0.2, colour = "black")+
  labs(x= "Season", y = expression("Growth rate (d"^"-1" *')'), color = "Treatment")+ #(μ)
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2.5))+
  scale_color_manual(values = palette_treatments_remei)+
  guides(size = 'none')+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), 
        axis.text.x = element_text(angle = 0), legend.position = "right",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text = element_text(size = 6))

### statistics treatments -----
anova<-aov(PRO~treatment, data = GR_dapis_OS_filt)
summary(anova) #p<0.05 ***
tuk<- TukeyHSD(anova)#
aov_residuals <- residuals(object = anova)
plot(anova, 1)
plot(anova, 2)
shapiro.test(x = aov_residuals )#From the output, the p-value > 0.05 implying that the distribution of the data are #not significantly different from normal distribution. In other words, we can assume the normality.
##p-value = 0.003229 NO normal 
#test no parametric
kruskal.test(PRO ~ treatment, data = GR_dapis_OS_filt) #
#p<0.05 hago dunn test 
#p-value = 1.375e-05
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

# CD    CL    DL    PD    PL    VL 
# "ab"   "a"   "c" "abc"  "bc"   "c" 

gr_dapis_treat <- 
  GR_dapis_OS_filt |>
  ggplot(aes(treatment, PRO))+ #
  geom_point(aes(color = season), alpha = 0.7, position = position_jitter(0.2), size = 1.5)+
  geom_violin(alpha = 0.1, draw_quantiles = c(0.25, 0.75))+
  geom_text(data = letters_dapis_treat, aes(y = 1.5, x = row.names(letters_dapis_treat), label = Letters),
            position = position_nudge(x = 0), hjust = 0, color = "black")+
  stat_summary(fun = "mean", geom = "crossbar", 
               width = 0.2, colour = "black")+
  labs(x= "Treatment", y = expression("Growth rate (d"^"-1" *')'), color = "Season")+ #(μ)
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2.5))+
  scale_color_manual(values = palette_seasons_4)+
  guides(size = 'none')+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0), 
        legend.position = "right",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        text = element_text(size = 6))


## mOTUs growth rates calculation ----
#import data
motus_table <- read_tsv('data/mOTUs_Adrià/motus_table_ed.txt')

# motus_table |>
#   dim() ##33570 mOTUs and 66 samples

tax <- motus_table |>
  dplyr::select(consensus_taxonomy, mOTUs_num, consensus_taxonomy_mOTUs_num)

motus_table_ed <- motus_table |>
  select(-consensus_taxonomy, -mOTUs_num)

##env data
env_data_mOTUs <- read_xlsx('data/envdata/metadata_remei_mOTUs.xlsx')

env_data_mOTUs <- env_data_mOTUs |> 
  rename(c("Treatment" = "treatment", "Light Regime" = "light_regime",
           "Replicate" = "replicate", "Time" = "time", "Hours" = "hours", "Season" = "season", 
           "Fraction" = "fraction", "Sample-Name" = "sample_name" ,"Sample Code" = "sample_code",  
           "amplicon_data_correspondence" = "amplicon_data_correspondence", "r1" = "r1", "r2" = "r2", 
           "Million\r\nread-pairs\r\n(or reads)" = 'reads', "Yield\r\n(Gb)" = 'yield'))

##edit t0 per duplicar t0 dark per light treatment també
change_t0s <- function(df){
  datch <- df %>% 
    dplyr::filter( treatment %in% c('CT', 'PF'))
  datcl <- datch %>%  
    dplyr::mutate( treatment = case_when( treatment == 'CT' ~ 'CL',
                                   treatment == 'PF' ~ 'PL',
                                   TRUE ~ treatment))
  datd <- datch %>%  
    dplyr::mutate( treatment = case_when( treatment == 'CT' ~ 'CD',
                                   treatment == 'PF' ~ 'PD',
                                   TRUE ~ treatment))
  whole <- bind_rows(datcl, datd) %>% 
    dplyr::mutate(selected_file_name_new = stringr::str_c('REMEIextra', 1:nrow(.)))
  df %>% 
    dplyr::filter( !treatment %in% c('CT', 'PF')) %>% 
    bind_rows(whole) %>% 
    return()
}

env_data_mOTU_ed <- change_t0s(env_data_mOTUs)

##add citometry to calculate pseudoabundances for mOTUs
rem_fc <- readRDS("./data/intermediate_files/remei_phyloseq_silva_fc.rds")
rem_fc_env <- rem_fc@sam_data |>
  as_tibble()

rem_fc_env_ed <- rem_fc_env |>
  dplyr::mutate(treatment = as.character(treatment),
         season = as.character(season)) |>
  dplyr::select(treatment, season, All, replicate, time, hours.x)

env_data_mOTU_ed_ct <- env_data_mOTU_ed |>
  left_join(rem_fc_env_ed,
            by = c('treatment' = 'treatment', 
                   'season' = 'season', 
                   'replicate' = 'replicate', 
                   'time' = 'time'))

##filter environmental dataset (failed)
env_data_mOTU_ed_ct <- 
  env_data_mOTU_ed_ct |>
  dplyr::filter(!is.na(amplicon_data_correspondence), amplicon_data_correspondence != "FAILED")

##filter mOTUs that were 0 in the whole dataset:
motus_table_t <- motus_table_ed |>
  t() |>
  as.data.frame() |>
  tibble::rownames_to_column(var = 'sample_code') |>
  janitor::row_to_names(row_number = 1, remove_row = T, remove_rows_above = F) 

motus_table_t <- motus_table_t |>
  dplyr::mutate(across(!consensus_taxonomy_mOTUs_num, as.numeric))

motus_table_t_codes <- motus_table_t$consensus_taxonomy_mOTUs_num

motus_table_t_filt <- motus_table_t |>
  select(-consensus_taxonomy_mOTUs_num) |>
  select(where(~ sum(.) != 0)) |>
  cbind(motus_table_t_codes)

##add hours information
motus_env <- env_data_mOTU_ed_ct |>
  select(hours, sample_code, time, season, treatment, replicate, yield, All) |>
  left_join(motus_table_t_filt, by = c('sample_code' = 'motus_table_t_codes'))

env_data_mOTUs$sample_code == motus_table_t_filt$motus_table_t_codes #check

##pivot long to compare t0 and t4 and divide keep only positive values
motus_env_l <- motus_env |>
  #select(-hours) |>
  pivot_longer(!c(sample_code, time, season, treatment, replicate, hours, yield, All), 
               names_to = 'consensus_taxonomy_mOTUs_num', 
               values_to = 'rel_abund') |>
  dplyr::mutate(time_ed = str_replace_all(time, c('t4'= 'tf', 't3' = 'tf')))

#calculate pseudoabundance for mOTUs
motus_env_l_w <- motus_env_l |>  
  dplyr::select(treatment, season, consensus_taxonomy_mOTUs_num, All, time_ed, replicate, yield, rel_abund, hours) |>
  #dplyr::mutate(rel_abund_ed = ifelse(is.na(rel_abund), 0, rel_abund)) |>
  dplyr::mutate(pseudoabund = All*rel_abund) |>
  group_by(treatment, season, consensus_taxonomy_mOTUs_num, time_ed, hours) |>
  dplyr::summarize(mean_pseudoabund = mean(pseudoabund)) |> ##replicates mean to calculate fold change
  distinct(treatment, season, consensus_taxonomy_mOTUs_num, time_ed, mean_pseudoabund) |>
  pivot_wider(id_cols = c(treatment, season, consensus_taxonomy_mOTUs_num), names_from = c(time_ed), values_from = mean_pseudoabund)

##calculate fold change and ln to calculate growth rate
motus_env_l_w_fold <- motus_env_l_w |>
  dplyr::filter(t0 != is.na(t0) &
           tf != is.na(tf)) |>
  dplyr::mutate(fold = tf/t0,
         gw = case_when(season == 'Winter' ~ log(tf/t0)/36,
                        season == 'Spring' ~ log(tf/t0)/48,
                        season == 'Summer' ~ log(tf/t0)/36)) |>
  dplyr::mutate(gr_day = gw*24) |>
  dplyr::filter(gw > 0 )

#write.csv(motus_env_l_w_fold, 'motus_env_l_w_fold_ln.csv')

motus_env_l_w_fold_mean <- 
  motus_env_l_w_fold |>
  dplyr::filter(gw > 0 &
                  consensus_taxonomy_mOTUs_num != 'unassigned') |>
  group_by(treatment, season) |>
  dplyr::summarize(mean_motus = mean(gr_day))

##select only the most abundant >1% community at some point------
motus_table_ed_melt <-  motus_table_ed |>
  pivot_longer(cols = !consensus_taxonomy_mOTUs_num, values_to = 'rel_abund', names_to = 'sample_code') %>%
  group_by(consensus_taxonomy_mOTUs_num) |>
  dplyr::filter(sum(rel_abund) > 0 &
                  consensus_taxonomy_mOTUs_num != 'unassigned') 

motus_table_ed_melt %$%
  consensus_taxonomy_mOTUs_num |>
  unique() #2715 mOTUs al meu dataset que no són 0

abundant_motus <- motus_table_ed_melt |> 
  ungroup() %>%
  dplyr::filter(rel_abund > 0.01) |> #more than 1% of the community at some point
  dplyr::select(consensus_taxonomy_mOTUs_num) |>
  unique() |>
  as_tibble()

motus_env_l_w_h_growth_1perc <- motus_env_l_w_fold |> 
  right_join(abundant_motus, by = "consensus_taxonomy_mOTUs_num", copy = TRUE) |>
  dplyr::filter(consensus_taxonomy_mOTUs_num != 'unassigned') 

mean_motus_1perc <- motus_env_l_w_h_growth_1perc |>
  dplyr::filter(gr_day > 0 &
         consensus_taxonomy_mOTUs_num != 'unassigned') |>
  group_by(treatment, season) |> 
  dplyr::mutate(mean_gr_motus_1perc = mean(gr_day)) |>
  dplyr::distinct(treatment, season, mean_gr_motus_1perc)

motus_gr_mean_1perc %>%
  left_join(motus_env_l_w_fold_mean)

##violin plot for mOTUs gr 1perc ----
motus_env_l_w_h_growth_1perc$season <- factor(motus_env_l_w_h_growth_1perc$season, c('Winter', 'Spring', 'Summer'))
motus_env_l_w_h_growth_1perc$treatment <- factor(motus_env_l_w_h_growth_1perc$treatment, c('CL', 'CD', 'PL', 'PD', 'DL', 'VL'))

seas_motus <- motus_env_l_w_h_growth_1perc |>
  dplyr::filter(gr_day > 0 &
                  consensus_taxonomy_mOTUs_num != 'unassigned') |>
  ggplot(aes(season, gr_day))+ 
  geom_point(aes(color = treatment), alpha = 0.7, position = position_jitter(0.25), size = 1.5)+
  geom_violin(alpha = 0.2, draw_quantiles = c(0.25, 0.75) )+
  labs(color = "Treatment")+
  labs(x= "Season", y = expression("Growth rate (d"^"-1" *')'))+  
  scale_y_continuous(expand = c(0, 0), limits = c(0,10.4))+
  scale_color_manual(values = palette_treatments_remei)+
  stat_summary(fun = "mean", geom = "crossbar", 
               width = 0.4, colour = "black")+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 16), axis.text.x = element_text(angle = 0), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 16))

treat_motus <- motus_env_l_w_h_growth_1perc |>
  ggplot(aes(treatment, gr_day))+ 
  geom_point(aes(color = season), alpha = 0.7, position = position_jitter(0.25), size = 1.5)+
  geom_violin(alpha = 0.2, draw_quantiles = c(0.25, 0.75) )+
  labs(color = "Season")+
  labs(x= "Treatment", y = expression("Growth rate (d"^"-1" *')'))+  
  scale_y_continuous(expand = c(0, 0), limits = c(0,10.4))+
  stat_summary(fun = "mean", geom = "crossbar", 
               width = 0.4, colour = "black")+
  scale_color_manual(values = palette_seasons_4)+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 16), axis.text.x = element_text(angle = 0), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 16))

#library(gridExtra)
motus_gr_1perc <- grid.arrange(treat_motus, seas_motus)
ggsave('motus_gr_1perc.pdf', motus_gr_1perc, path = "~/Documentos/Doctorat/REMEI/results/figures/corrected_asv_num/",
       width = 130,
       height = 180,
       units = 'mm' )

##violin plot for mOTUs gr----
motus_env_l_w_fold$season <- factor(motus_env_l_w_fold$season, c('Winter', 'Spring', 'Summer'))
motus_env_l_w_fold$treatment <- factor(motus_env_l_w_fold$treatment, c('CL', 'CD', 'PL', 'PD', 'DL', 'VL'))

seas_motus <- motus_env_l_w_fold |>
  ggplot(aes(season, gr_day))+ 
  geom_point(aes(color = treatment), alpha = 0.7, position = position_jitter(0.25), size = 1.5)+
  geom_violin(alpha = 0.2, draw_quantiles = c(0.25, 0.75) )+
  labs(color = "Treatment")+
  labs(x= "Season", y = expression("Growth rate (d"^"-1" *')'))+  
  scale_y_continuous(expand = c(0, 0), limits = c(0,10.4))+
  scale_color_manual(values = palette_treatments_remei)+
  stat_summary(fun = "mean", geom = "crossbar", 
               width = 0.4, colour = "black")+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 16), axis.text.x = element_text(angle = 0), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 16))

treat_motus <- motus_env_l_w_fold |>
  ggplot(aes(treatment, gr_day))+ 
  geom_point(aes(color = season), alpha = 0.7, position = position_jitter(0.25), size = 1.5)+
  geom_violin(alpha = 0.2, draw_quantiles = c(0.25, 0.75) )+
  labs(color = "Season")+
  labs(x= "Treatment", y = expression("Growth rate (d"^"-1" *')'))+  
  scale_y_continuous(expand = c(0, 0), limits = c(0,10.4))+
  stat_summary(fun = "mean", geom = "crossbar", 
               width = 0.4, colour = "black")+
  scale_color_manual(values = palette_seasons_4)+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 16), axis.text.x = element_text(angle = 0), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 16))

motus_gr <- grid.arrange(treat_motus, seas_motus)
ggsave('motus_gr_all.pdf', motus_gr, path = "~/Documentos/Doctorat/REMEI/results/figures/corrected_asv_num/",
       width = 130,
       height = 180,
       units = 'mm' )

##comparison with dapis gr and mOTUS gr CORRELATION PLOTS-----
mean_gr_dapis <- GR_dapis_OS_filt |>
  group_by(treatment, season) |>
  mutate(mean_dapis = mean(as.numeric(PRO))) |>
  distinct(season, treatment, mean_dapis)

mean_gr_asv <- reg_all_slopes_chosen_silva_tax |>
  group_by(treatment, season) |>
  mutate(mean_asv = mean(as.numeric(slope_chosen_days))) |>
  distinct(season, treatment, mean_asv)

mean_gr_asv_1perc <- reg_all_slopes_chosen_silva_tax_1perc |>
  dplyr::filter(season != "Early_fall") |>
  mutate(slope_chosen_days_all = case_when(is.na(slope_chosen_days) ~ '0',
                                           slope_chosen_days < 0 ~ '0',
                                           pvalue_slope_chosen > 0.05 ~ '0',
                                           TRUE ~ as.character(slope_chosen_days))) |>
  group_by(treatment, season) |>
  mutate(mean_asv_1perc = mean(as.numeric(slope_chosen_days_all))) |>
  distinct(season, treatment, mean_asv_1perc)

mean_gr_asv_motus_dapis <-  motus_gr_mean_1perc |>
  right_join(mean_gr_dapis) |>
  left_join(mean_gr_asv) |> #mean_gr_asv_filt (filtrat per >1%)
  left_join(motus_env_l_w_fold_mean) |> #comunitat general no amb l'1%.
  left_join(mean_gr_asv_1perc)

#write.table(mean_gr_asv_motus_dapis, 'results/tables/mean_treat_seas_motus_asv_dapi.txt', sep = '\t')

##check normality before correlation
shapiro.test(as.numeric(mean_gr_asv_motus_dapis$mean_gr_motus_1perc)) # => p-value = 0.1389 (NORMALITY)
ggqqplot(as.numeric(mean_gr_asv_motus_dapis$mean_gr_motus_1perc))

shapiro.test(as.numeric(mean_gr_asv_motus_dapis$mean_asv)) # => p-value = 0.1757 (NORMALITY)
ggqqplot(as.numeric(mean_gr_asv_motus_dapis$mean_asv))

shapiro.test(as.numeric(mean_gr_asv_motus_dapis$mean_asv_1perc)) # => p-value = 0.1818 (NORMALITY)
ggqqplot(as.numeric(mean_gr_asv_motus_dapis$mean_asv_1perc))

shapiro.test(as.numeric(mean_gr_asv_motus_dapis$mean_dapis)) # => p-value = 0.0349 (No NORMALITY)
ggqqplot(as.numeric(mean_gr_asv_motus_dapis$mean_dapis))

mean_gr_asv_motus_dapis$treatment <- factor(mean_gr_asv_motus_dapis$treatment, 
                                               levels = c('CL', 'CD', 'PL', 'PD', 'DL', 'VL'))

gr_dapis_asv <- mean_gr_asv_motus_dapis |>
  as_tibble() |>
  ggplot(aes(y = mean_asv, x = mean_dapis))+
  geom_point(aes(shape = treatment, color = season), size = 1.5)+
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0))+
  scale_y_continuous(limits = c(0, 4.5), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  scale_shape_manual(values= shapes_treatments)+
  stat_cor(aes(label = paste(..p.label..)), label.x = 0.1, label.y = 4, 
           p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman')+
  guides(size = 'none',
         color = 'none')+
  labs(y = 'Mean growth rates ASV-based', x = 'Mean growth rates abundance-based', 
       color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/4), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        text = element_text(size = 6), legend.position = 'none')

gr_asv_mOTUs <- mean_gr_asv_motus_dapis |>
  as_tibble() |>
  ggplot(aes(mean_gr_motus_1perc,  mean_asv_1perc))+ #mean_motus,
  geom_point(aes(shape = treatment, color = season), size = 1.5)+
  scale_y_continuous(limits = c(0, 4.5), expand = c(0, 0))+
  scale_x_continuous(limits = c(0, 2.1), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  stat_cor(aes(label = paste(..p.label..)),label.x = 0.2, label.y = 4, 
           p.digits = 0.01, digits = 2, 
           p.accuracy = 0.01, method = 'pearson')+
  guides(size = 'none',
         color= 'none')+
  scale_shape_manual(values= shapes_treatments)+
  labs(x = 'Mean growth rates mOTUs-based', y = 'Mean growth rates ASV-based', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/4), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        text = element_text(size = 6))

##calculate correlations ----
#for non-normal data (pvalue < 0.05 in Shapiro test)
mean_gr_asv_motus_dapis |>
  colnames()

cor_spearman <- rcorr(as.matrix(mean_gr_asv_motus_dapis[,c(4,5)]), type = 'spearman')

#for normal data (pvalue > 0.05 in Shapiro test)
cor_pear <- rcorr(as.matrix(mean_gr_asv_motus_dapis[,c(3,7)]), type = 'pearson')


# gr_asv_motus_all <- mean_gr_asv_motus_dapis |>
#   as_tibble() |>
#   ggplot(aes(mean_motus,  mean_asv))+ #mean_motus,
#   geom_point(aes(shape = treatment, color = season), size = 1.5)+
#   scale_y_continuous(limits = c(0, 4.5), expand = c(0, 0))+
#   scale_x_continuous(limits = c(0, 2.1), expand = c(0, 0))+
#   scale_color_manual(values = palette_seasons_4)+
#   stat_poly_line(color = 'black')+
#   stat_cor(aes(label = paste(..p.label..)),label.x = 0.2, label.y = 4, 
#            p.digits = 0.01, digits = 2, 
#            p.accuracy = 0.01, method = 'pearson')+
#   guides(size = 'none',
#          color= 'none')+
#   scale_shape_manual(values= shapes_treatments)+
#   labs(x = 'Mean growth rates mOTUs-based', y = 'Mean growth rates ASV-based', color = 'Season', shape = 'Treatment')+
#   theme_bw()+
#   theme(aspect.ratio = (4/4), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         text = element_text(size = 6))

# ggsave('GR_asv_mOTUS.pdf',  gr_asv_motus, path = "~/Documentos/Doctorat/REMEI/results/figures/corrected_asv_num/",
#        height = 88,
#        width = 88,
#        units = 'mm')
# # 
# ggsave('GR_dapis_mOTUS.pdf',  gr_dapis_mOTUs, path = "~/Documentos/Doctorat/REMEI/results/figures/corrected_asv_num/",
#        height = 88,
#        width = 88,
#        units = 'mm')

## Panel construction figure 2 ----
gr_asv_dapis <- grid.arrange(treat, gr_dapis_treat, seas, gr_dapis_seas, gr_dapis_asv, gr_asv_mOTUs)
# 
ggsave('GR_asv_dapis_stat_ed4.pdf', gr_asv_dapis, path = "~/Documentos/Doctorat/REMEI/results/figures/corrected_asv_num/manuscript_figures/",
       height = 190,
       width = 188,
       units = 'mm')
## Correlation between growth rates from 16S rRNA copy number normalization and 16S rRNA copy number no normalized ----
reg_all_slopes_chosen_silva_tax_cn_nor <- read.csv(
  "data/intermediate_files/reg_all_slopes_chosen_silva_tax_cn_nor.csv", 
  sep = ',')
#### prepare datasets
reg_all_slopes_chosen_silva_tax_cn_nor_filt <- reg_all_slopes_chosen_silva_tax_cn_nor |>
  dplyr::filter(pvalue_slope_chosen < 0.05 &
                  slope_chosen > 0) |>
  mutate(slope_chosen_days = slope_chosen*24) ##transform slopes /hours to /days

reg_all_slopes_chosen_silva_tax_cn_nor_filt %$%
  unique(asv_num) #1270 ASVs

reg_all_slopes_chosen_silva_tax_cn_nor_non <- reg_all_slopes_chosen_silva_tax_cn_nor_filt |>
  right_join(reg_all_slopes_chosen_silva_tax_filt, by = c('treatment', 'season', 'asv_num', 'phylum', 'class', 'order', 
                                                         'family'), suffix = c('.cn_nor', '.non')) 
##before correlation we check normality
shapiro.test(as.numeric(reg_all_slopes_chosen_silva_tax_cn_nor_non$slope_chosen_days.cn_nor)) # =>p-value < 2.2e-16 (NO NORMALITY)
ggqqplot(as.numeric(reg_all_slopes_chosen_silva_tax_cn_nor_non$slope_chosen_days.cn_nor))

shapiro.test(as.numeric(reg_all_slopes_chosen_silva_tax_cn_nor_non$slope_chosen_days.non)) # => p-value < 2.2e-16 (NO NORMALITY)
ggqqplot(as.numeric(reg_all_slopes_chosen_silva_tax_cn_nor_non$slope_chosen_days.non))

reg_all_slopes_chosen_silva_tax_cn_nor_non$season <- factor(reg_all_slopes_chosen_silva_tax_cn_nor_non$season,
                                                            levels = c('Winter', 'Spring', 'Summer', 'Fall'))

reg_all_slopes_chosen_silva_tax_cn_nor_non$treatment <- factor(reg_all_slopes_chosen_silva_tax_cn_nor_non$treatment,
                                                            levels = c('CL', 'CD', 'PL', 'PD', 'DL', 'VL'))

corr_cn_nor_non <- reg_all_slopes_chosen_silva_tax_cn_nor_non |>
  ggplot(aes(slope_chosen_days.cn_nor, slope_chosen_days.non))+
  geom_point(aes(shape = treatment, color = season), alpha = 0.6, size = 0.5)+
  scale_fill_manual(values = palette_seasons_4)+
  scale_y_continuous(limits = c(0, 10), expand = c(0, 0))+
  scale_x_continuous(limits = c(0, 10), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed')+
  scale_shape_manual(values = shapes_treatments)+
  stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, label.y = 9.5, p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman')+
  guides(size = 'none', shape = guide_legend(override.aes = list(size = 1.5)),
         color= guide_legend(override.aes = list(size = 1.5)))+
  labs(y = 'ASV-based growth rates', 
       x = 'ASV-based growth rates normalized\nby 16S rRNA copy number', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/4), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        text = element_text(size = 7))

# ggsave('correlation_gr_cn_nor_non_ed2.pdf', corr_cn_nor_non, 
#        path = "~/Documentos/Doctorat/REMEI/results/figures/corrected_asv_num/",
#        width = 88,
#        height = 88,
#        units = 'mm')


##for non-normal data (pvalue < 0.05 in Shapiro test)
reg_all_slopes_chosen_silva_tax_cn_nor_non |>
  colnames()
cor_spearman <- rcorr(as.matrix(reg_all_slopes_chosen_silva_tax_cn_nor_non[,c(18,23)]), type = 'spearman')

my_cor_matrix_spear <- flat_cor_mat(cor_spearman$r, cor_spearman$P)
head(my_cor_matrix_spear)

### by different treatments and seasons----
corr_cn_nor_non_seas <- reg_all_slopes_chosen_silva_tax_cn_nor_non |>
  ggplot(aes(slope_chosen_days.cn_nor, slope_chosen_days.non))+
  geom_point(aes(shape = treatment, color = season), alpha = 0.6, size = 0.5)+
  scale_fill_manual(values = palette_seasons_4)+
  scale_y_continuous(limits = c(0, 10), expand = c(0, 0))+
  scale_x_continuous(limits = c(0, 10), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(aes(group = season, color = season))+
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed')+
  scale_shape_manual(values = shapes_treatments)+
  stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, label.y = 9.5, p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman')+
  guides(size = 'none', shape = guide_legend(override.aes = list(size = 1.5)),
         color= guide_legend(override.aes = list(size = 1.5)))+
  labs(y = 'ASV-based growth rates', 
       x = 'ASV-based growth rates normalized\nby 16S rRNA copy number', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/4), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        text = element_text(size = 7))

# ggsave('correlation_gr_cn_nor_non_seasons_ed2.pdf', corr_cn_nor_non_seas, 
#        path = "~/Documentos/Doctorat/REMEI/results/figures/corrected_asv_num/",
#        width = 88,
#        height = 88,
#        units = 'mm')

corr_cn_nor_non_treat <- reg_all_slopes_chosen_silva_tax_cn_nor_non |>
  ggplot(aes(slope_chosen_days.cn_nor, slope_chosen_days.non))+
  geom_point(aes(shape = treatment, color = treatment), alpha = 0.6, size = 0.5)+
  scale_fill_manual(values = palette_seasons_4)+
  scale_y_continuous(limits = c(0, 10), expand = c(0, 0))+
  scale_x_continuous(limits = c(0, 10), expand = c(0, 0))+
  scale_color_manual(values = palette_treatments_remei)+
  stat_poly_line(aes(group = treatment, color = treatment))+
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed')+
  scale_shape_manual(values = shapes_treatments)+
  stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, label.y = 9.5, p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman')+
  guides(size = 'none', shape = guide_legend(override.aes = list(size = 1.5)),
         color= guide_legend(override.aes = list(size = 1.5)))+
  labs(y = 'ASV-based growth rates', 
       x = 'ASV-based growth rates normalized\nby 16S rRNA copy number', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/4), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        text = element_text(size = 7))

# ggsave('correlation_gr_cn_nor_non_treat_ed2.pdf', corr_cn_nor_non_treat, 
#        path = "~/Documentos/Doctorat/REMEI/results/figures/corrected_asv_num/",
#        width = 88,
#        height = 88,
#        units = 'mm')


# ----- 6. GROWTH RATES DISTRIBUTION PATTERNS AT DIFFERENT TAXONOMIC RANKS ################ ----------------
## prepare the dataset for plotting ----
## reorder factors for plots
reg_all_slopes_chosen_silva_tax$treatment <- reg_all_slopes_chosen_silva_tax$treatment |> 
factor(levels=(c("CL", "CD", "PL", "PD", "DL", "VL")))
reg_all_slopes_chosen_silva_tax$season <- reg_all_slopes_chosen_silva_tax$season |> 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))

# #transform slopes /hours to /days
reg_all_slopes_chosen_silva_tax <- reg_all_slopes_chosen_silva_tax |>
  dplyr::mutate(slope_chosen_days = slope_chosen*24)

## new factors for taxonomy to plot them correctly
reg_all_slopes_chosen_silva_tax <- reg_all_slopes_chosen_silva_tax |>
  mutate(domain_f = as_factor(domain),
         phylum_f = as_factor(phylum),
         class_f = as_factor(class),
         order_f = as_factor(order),
         family_f = as_factor(family),
         genus_f = as_factor(genus),
         asv_f = as_factor(asv_num))


reg_all_slopes_chosen_silva_tax |>
  colnames()

reg_all_slopes_chosen_silva_tax$phylum_f <-  factor(reg_all_slopes_chosen_silva_tax$phylum_f, 
                                                    levels=unique(
                                                      reg_all_slopes_chosen_silva_tax$phylum_f[
                                                        order(reg_all_slopes_chosen_silva_tax$domain_f)]), 
                                                    ordered=TRUE)
reg_all_slopes_chosen_silva_tax$phylum_f |> levels()

reg_all_slopes_chosen_silva_tax$class_f <-  factor(reg_all_slopes_chosen_silva_tax$class_f, 
                                                   levels=unique(
                                                     reg_all_slopes_chosen_silva_tax$class_f[
                                                       order(reg_all_slopes_chosen_silva_tax$domain_f,
                                                             reg_all_slopes_chosen_silva_tax$phylum_f)]), 
                                                   ordered=TRUE)

reg_all_slopes_chosen_silva_tax$order_f <-  factor(reg_all_slopes_chosen_silva_tax$order_f, 
                                                   levels=unique(
                                                     reg_all_slopes_chosen_silva_tax$order_f[
                                                       order(reg_all_slopes_chosen_silva_tax$domain_f,
                                                             reg_all_slopes_chosen_silva_tax$phylum_f,
                                                             reg_all_slopes_chosen_silva_tax$class_f)]), 
                                                   ordered=TRUE)

reg_all_slopes_chosen_silva_tax$family_f <-  factor(reg_all_slopes_chosen_silva_tax$family_f, 
                                                    levels=unique(
                                                      reg_all_slopes_chosen_silva_tax$family_f[
                                                        order(reg_all_slopes_chosen_silva_tax$domain_f,
                                                              reg_all_slopes_chosen_silva_tax$phylum_f,
                                                              reg_all_slopes_chosen_silva_tax$class_f,
                                                              reg_all_slopes_chosen_silva_tax$order_f)]), 
                                                    ordered=TRUE)

reg_all_slopes_chosen_silva_tax$genus_f <-  factor(reg_all_slopes_chosen_silva_tax$genus_f, 
                                                   levels=unique(
                                                     reg_all_slopes_chosen_silva_tax$genus_f[
                                                       order(reg_all_slopes_chosen_silva_tax$domain_f,
                                                             reg_all_slopes_chosen_silva_tax$phylum_f,
                                                             reg_all_slopes_chosen_silva_tax$class_f,
                                                             reg_all_slopes_chosen_silva_tax$order_f,
                                                             reg_all_slopes_chosen_silva_tax$family_f)]), 
                                                   ordered=TRUE)

reg_all_slopes_chosen_silva_tax$asv_f <-  factor(reg_all_slopes_chosen_silva_tax$asv_f, 
                                                 levels=unique(
                                                   reg_all_slopes_chosen_silva_tax$asv_f[
                                                     order(reg_all_slopes_chosen_silva_tax$domain_f,
                                                           reg_all_slopes_chosen_silva_tax$phylum_f,
                                                           reg_all_slopes_chosen_silva_tax$class_f,
                                                           reg_all_slopes_chosen_silva_tax$order_f,
                                                           reg_all_slopes_chosen_silva_tax$family_f,
                                                           reg_all_slopes_chosen_silva_tax$genus_f)]), 
                                                 ordered=TRUE)

## Distribution plots at different taxonomic ranks ----
reg_all_slopes_chosen_silva_tax <- reg_all_slopes_chosen_silva_tax  |>
  ungroup() |>
  mutate(phylum_f = fct_rev(fct_infreq(phylum_f)))

counts <-   reg_all_slopes_chosen_silva_tax |> #number of growth rates calculated by phylum
  group_by(phylum_f) |>
  dplyr::filter(n() >= 2) |>
  mutate(counts = n()) |>
  as_tibble() |>
  group_by(phylum_f, counts) |>
  dplyr::summarize() |>
  as_tibble()

ridge_ph <- # distribution plot at phylum taxonomic rank
  reg_all_slopes_chosen_silva_tax |>
  group_by(phylum_f) |>
  dplyr::filter(n() >= 2) |>
  mutate(counts = paste('n = ', n())) |> #
  ggplot(aes(y = fct_rev(fct_infreq(phylum_f)), x = slope_chosen_days, fill = phylum_f,  group = phylum_f, label = counts))+
  geom_density_ridges(alpha = 0.8, panel_scaling = TRUE, scale = 1,
                      jittered_points = TRUE,
                      point_shape = 21, point_size = 0.2, point_alpha = 0.0,
                      quantile_lines = TRUE,
                      quantile_fun = mean)+
  scale_fill_manual(values = palette_phylums_assigned)+
  scale_x_continuous(limits = c(0,11), expand = c(0,0))+
  coord_cartesian(clip = "off") +
  geom_text(nudge_x = 8.8, nudge_y = 0.3, check_overlap = TRUE, size = 3)+
  labs(y = '', fill= 'Phylum', x = expression("Growth rate (d"^"-1" *')'), title = 'Phylum')+
  theme_bw()+
  theme(legend.position = "none", strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0, size = 6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 6),
        plot.margin= unit(c(0.2, 0.5, 0.2, 0.2), "lines"), panel.border = element_blank(), 
        title = element_text(size = 8))

ridge_cl <- reg_all_slopes_chosen_silva_tax |> # distribution plot at class taxonomic rank
  group_by(phylum_f) |>
  dplyr::filter(n() >= 2) |>
  group_by(phylum_f) |>
  mutate(counts = paste('n = ', n())) |> #
  ggplot(aes(y = fct_rev(fct_infreq(phylum_f)), x = slope_chosen_days, fill = phylum_f, label = counts))+
  geom_density_ridges(aes(fill = phylum_f, group = class_f), scale = 1, alpha = 0.7,
                      jittered_points = TRUE,
                      point_shape = 21, point_size = 0.2, point_alpha = 0.0,
                      quantile_lines = TRUE,
                      quantile_fun = mean)+
  scale_fill_manual(values = palette_phylums_assigned)+
  labs(y = 'Class', fill= 'Phylum', x =expression("Growth rate (d"^"-1" *')'), title = 'Class')+
  scale_x_continuous(limits = c(0,11), expand = c(0,0))+
  geom_text(nudge_x = 7.8, nudge_y = 0.3, check_overlap = TRUE, size = 3)+
  theme_bw()+
  theme(legend.position = "none", strip.text.x = element_text(size = 6), 
        axis.text.x = element_text(angle = 0, size = 6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 0),
        axis.title.y = element_blank(), 
        plot.margin= unit(c(0.2, 0.5, 0.2, 0.2), "lines"), panel.border = element_blank(), 
        title = element_text(size = 8))

ridge_o <- reg_all_slopes_chosen_silva_tax |> # distribution plot at order taxonomic rank
  group_by(phylum_f) |>
  dplyr::filter(n() >= 2) |>
  ungroup() |>
  group_by(phylum_f) |>
  mutate(counts = paste('n = ', n())) |> #
  ggplot(aes(y = fct_rev(fct_infreq(phylum_f)), x = slope_chosen_days, fill = phylum_f, label = counts))+
  geom_density_ridges(aes(fill = phylum_f, group = order_f), scale = 0.8, alpha = 0.6,
                      jittered_points = TRUE,
                      point_shape = 21, point_size = 0.2, point_alpha = 0.0,
                      quantile_lines = TRUE,
                      quantile_fun = mean)+
  scale_fill_manual(values = palette_phylums_assigned)+
  labs(y = '', fill= 'Phylum', x = expression("Growth rate (d"^"-1" *')'), title = 'Order')+
  scale_x_continuous(limits = c(0,11), expand = c(0,0))+
  geom_text(nudge_x = 7.8, nudge_y = 0.3, check_overlap = TRUE, size = 3)+
  theme_bw()+
  theme(legend.position = "none", strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0, size = 6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 0),
        plot.margin= unit(c(0.2, 0.5, 0.2, 0.2), "lines"), panel.border = element_blank(), 
        title = element_text(size = 8))

ridge_f <- reg_all_slopes_chosen_silva_tax |> # distribution plot at order taxonomic rank
  group_by(phylum_f) |>
  dplyr::filter(n() >= 2) |>
  group_by(phylum_f) |>
  mutate(counts = paste('n = ', n())) |>
  ungroup() |>
  ggplot(aes(y = fct_rev(fct_infreq(phylum_f)), x = slope_chosen_days, fill = phylum_f, label = counts))+
  geom_density_ridges(aes(fill = phylum_f, group = family_f), scale = 0.8, alpha = 0.6,
                      jittered_points = TRUE,
                      point_shape = 21, point_size = 0.2, point_alpha = 0.0,
                      quantile_lines = TRUE,
                      quantile_fun = mean)+
  scale_fill_manual(values = palette_phylums_assigned)+
  labs(y = '', fill= 'Phylum', x = expression("Growth rate (d"^"-1" *')'), title = 'Family')+
  scale_x_continuous(limits = c(0,11), expand = c(0,0))+
  geom_text(nudge_x = 7.8, nudge_y = 0.3, check_overlap = TRUE, size = 3)+
  theme_bw()+
  theme(legend.position = "none", strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0, size = 6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 0),
        plot.margin= unit(c(0.2, 0, 0.2, 0.2), "lines"), panel.border = element_blank(), 
        title = element_text(size = 8))

## reorder ASVs by phylum
reg_all_slopes_chosen_silva_tax_asv_reorder <- reg_all_slopes_chosen_silva_tax |>
  group_by(phylum_f) |>
  mutate(counts = paste('n = ', n())) |> #
  ungroup()

reg_all_slopes_chosen_silva_tax_asv_reorder$phylum_f <- reg_all_slopes_chosen_silva_tax_asv_reorder$phylum_f |>
  factor( levels = c("Bacteroidota", "Proteobacteria", "Actinobacteriota",   "Cyanobacteria", "Planctomycetota",
                     "Verrucomicrobiota",  "Firmicutes", "Crenarchaeota", "Bdellovibrionota", "Campilobacterota", 
                     "Acidobacteriota",
                     "Nitrospinota",  "Nitrospirota",  "Myxococcota", "Desulfobacterota",
                     "Deinococcota" , "Chloroflexi" , "Fusobacteriota", 
                     "Spirochaetota",  "Abditibacteriota", "Latescibacterota",
                     "Methylomirabilota", "Halanaerobiaeota", "Sumerlaeota", "Calditrichota", "Gemmatimonadota"), 
          ordered = TRUE)

ridge_a <- reg_all_slopes_chosen_silva_tax_asv_reorder |> # distribution plot at ASV taxonomic rank
  group_by(phylum_f) |>
  dplyr::filter(n() >= 2) |>
  ggplot(aes(y = fct_rev(phylum_f), x = slope_chosen_days, label = counts))+
  geom_density_ridges(aes(fill = phylum_f, group = asv_f), scale = 0.6, alpha = 0.6,
                      jittered_points = TRUE,
                      point_shape = 21, point_size = 0.2, point_alpha = 0.0,
                      quantile_lines = TRUE,
                      quantile_fun = mean)+
  scale_fill_manual(values = palette_phylums_assigned)+
  labs(y = '', fill= 'Phylum', x = expression("Growth rate (d"^"-1" *')'), title = 'ASV')+
  scale_x_continuous(limits = c(0,11), expand = c(0,0))+
  geom_text(nudge_x = 7.2, nudge_y = 0.35, check_overlap = TRUE, size = 3)+
  coord_cartesian(clip = "off") +
  theme_bw()+
  theme(legend.position = "none", strip.text.x = element_text(size = 6), 
        axis.text.x = element_text(angle = 0, size = 6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 0),
        plot.margin= unit(c(0.2, 0.5, 0.2, 0.2), "lines"), panel.border = element_blank(), 
        title = element_text(size = 8))

## Percentage of significant growth rates calculated by Phylum ------ 
perc_gr_phylum <- reg_all_slopes_chosen_silva_tax |>
  group_by(phylum_f) |>
  dplyr::filter(n() > 2) |>
  group_by(phylum_f, slope_chosen_days) |> 
  dplyr::summarize(n = n()) |>
  group_by(phylum_f) |>
  dplyr::summarize(n = n()) |>
  mutate(sum = sum(n),
         perc = n/sum) |>
  ungroup() |>
  ggplot(aes(sum, perc, fill = phylum_f))+
  geom_col(aes(fill = fct_reorder(phylum_f, perc, .desc = TRUE)))+ 
  scale_y_continuous(labels = percent_format())+
  labs(y = 'Percentage of significant growth rates calculated per Phylum', fill = 'Phylum')+
  scale_fill_manual(values = palette_phylums_assigned)+
  coord_flip()+
  theme_bw()+
  theme(axis.text.y = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), aspect.ratio = 2/10, panel.border = element_blank(),
        axis.ticks = element_blank(), legend.position = 'none', axis.text.x = element_text(size = 6),
        axis.title.y = element_blank(), axis.text = element_text(size = 6), text = element_text(size = 6),
        plot.margin= unit(c(1, 0.5, 1, 0.5), "lines"),)

## Percentage of significant growth rates calculated at different taxonomic ranks ----
x <- reg_all_slopes_chosen_silva_tax |>
  group_by(phylum_f) |>
  dplyr::filter(n() > 2) |>
  dplyr::summarize(n = n()) |>
  dplyr::summarize(n = sum(n))
  
gr_asv_perc <- reg_all_slopes_chosen_silva_tax |>
  group_by(asv_f) |>
  dplyr::filter(n() > 2) |>
  dplyr::summarize(n = n()) |>
  mutate(sum = sum(n),
         asv_gr_perc = sum/x) |>
  dplyr::select(asv_gr_perc) |>
  unique() |>
  as_tibble()

gr_f_perc <- reg_all_slopes_chosen_silva_tax |>
  group_by(family) |>
  dplyr::filter(n() > 2) |>
  dplyr::summarize(n = n()) |>
  mutate(sum = sum(n),
         f_gr_perc = sum/x) |>
  dplyr::select(f_gr_perc) |>
  unique() |>
  as_tibble()

gr_o_perc <- reg_all_slopes_chosen_silva_tax |>
  group_by(order) |>
  dplyr::filter(n() > 2) |>
  dplyr::summarize(n = n()) |>
  mutate(sum = sum(n),
         order_gr_perc = sum/x) |>
  dplyr::select(order_gr_perc) |>
  unique() |>
  as_tibble()

gr_c_perc <- reg_all_slopes_chosen_silva_tax |>
  group_by(class) |>
  dplyr::filter(n() > 2) |>
  dplyr::summarize(n = n()) |>
  mutate(sum = sum(n),
         class_gr_perc = sum/x) |>
  dplyr::select(class_gr_perc) |>
  unique() |>
  as_tibble()

gr_p_perc <- reg_all_slopes_chosen_silva_tax |>
  group_by(phylum) |>
  dplyr::filter(n() > 2) |>
  dplyr::summarize(n = n()) |>
  mutate(sum = sum(n),
         phylum_gr_perc = sum/x) |>
  dplyr::select(phylum_gr_perc) |>
  unique() |>
  as_tibble() 

labels_perc_gr_rank <-  as_labeller(c(phylum_gr_perc = 'Phylum',
                                      class_gr_perc = 'Class',
                                      order_gr_perc = 'Order',
                                      f_gr_perc = 'Family',
                                      asv_gr_perc = 'ASV'))

gr_tax_ranks_perc <- cbind(gr_asv_perc, gr_f_perc, gr_o_perc, gr_c_perc, gr_p_perc) |>
  as_tibble() |>
  pivot_longer(cols = 1:5) |>
  ggplot(aes(x = fct_rev(as_factor(name)), y = value$n, group = 1))+
  geom_point()+
  #coord_flip()+
  geom_line()+
  labs(y='Percentage', x = 'Growth rates distribution plotted\nat different taxonomic ranks')+
  scale_x_discrete(labels = labels_perc_gr_rank)+
  scale_y_continuous(labels = percent_format())+
  theme_bw()+
  theme(aspect.ratio = 3/10, panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = 'none', axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 6), text = element_text(size = 6),
        plot.margin= unit(c(1, 0.5, 1, 0.5), "lines"),)

## mean growth rates at different taxonomic levels ----
labels_mean_rank <- as_labeller(c( mean_phylum = 'Phylum',
                                   mean_class = 'Class',
                                   mean_order = 'Order',
                                   mean_family = 'Family',
                                   mean_asv = 'ASV'))

mean_gr_ranks <- reg_all_slopes_chosen_silva_tax |>
  group_by(phylum_f) |>
  dplyr::filter(n() > 2) |>
  mutate(mean_phylum = mean(slope_chosen_days)) |>
  ungroup() |>
  group_by(class) |>
  mutate(mean_class = mean(slope_chosen_days)) |>
  ungroup() |>
  group_by(order) |>
  mutate(mean_order = mean(slope_chosen_days)) |>
  ungroup() |>
  group_by(family) |>
  mutate(mean_family = mean(slope_chosen_days)) |>
  ungroup() |>
  group_by(asv_num) |>
  mutate(mean_asv = mean(slope_chosen_days)) |>
  dplyr::select(starts_with('mean'), 'asv_num', 'family', 'order', 'class', 'phylum_f') |>
  pivot_longer(cols = starts_with('mean')) |>
  group_by(phylum_f) |>
  distinct(value, name, phylum_f) |>
  ggplot(aes(as_factor(name), value))+
  geom_point(aes(color = phylum_f), alpha = 0.6, 
             position = position_jitter(0.3))+
  geom_violin(alpha = 0.2, draw_quantiles = c(0.25, 0.75))+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.6,
               colour = "black")+
  scale_color_manual(values = palette_phylums_assigned)+
  labs(y = 'Mean growth rate\nat different\ntaxonomic ranks', x= 'Taxonomic rank', color = 'Phylum')+
  scale_x_discrete(labels = labels_mean_rank)+
  #coord_flip()+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = 'none', axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 6), text = element_text(size = 6),
        plot.margin= unit(c(1, 0.5, 1, 0.5), "lines"))
# 
# ggsave('mean_gr_ranks_ed.pdf', mean_gr_ranks,
# path = "~/Documentos/Doctorat/REMEI/results/figures/corrected_asv_num/",
# width = 150,
# height = 130,
# units = 'mm')


## Figure 3 creating panel design------
distribution_gr_plots_ridge <- multi_panel_figure(columns = 6, rows = 10, width = 188, height = 250, 
                                                  column_spacing = 0.2, unit = 'mm', row_spacing = 0.6,
                                                  panel_label_type = 'upper-alpha')

# distribution_gr_plots_ridge  %<>%
#   fill_panel(perc_gr_phylum, column = 1, row = 1) %<>%
#   fill_panel(gr_tax_ranks_perc, column = 1, row = 2) %<>%
#   fill_panel(mean_gr_ranks, column = 1, row = 3) %<>%
#   fill_panel(ridge_ph, column = 2:3, row = 1:3) %<>%
#   fill_panel(ridge_cl, column = 4, row = 1:3) %<>%
#   fill_panel(ridge_o, column = 5, row = 1:3) %<>%
#   fill_panel(ridge_f, column = 6, row = 1:3) %<>%
#   fill_panel(ridge_a, column = 7, row = 1:3)

distribution_gr_plots_ridge  %<>%
  fill_panel(perc_gr_phylum, column = 1:2, row = 1:2) %<>%
  fill_panel(gr_tax_ranks_perc, column = 3:4, row = 1:2) %<>%
  fill_panel(mean_gr_ranks, column = 5:6, row = 1:2) %<>%
  fill_panel(ridge_ph, column = 1:2, row = 3:10) %<>%
  fill_panel(ridge_cl, column = 3, row = 3:10) %<>%
  fill_panel(ridge_o, column = 4, row = 3:10) %<>%
  fill_panel(ridge_f, column = 5, row = 3:10) %<>%
  fill_panel(ridge_a, column = 6, row = 3:10)

# save the graph
ggsave('distribution_gr_plots_ridge_order_freq_perc_ed5.pdf', distribution_gr_plots_ridge,
       path = "~/Documentos/Doctorat/REMEI/results/figures/corrected_asv_num/",
       width = 188,
       height = 270,
       units = 'mm')

## Distribution plots separated by families (figure S3) ----

counts_divided <-   reg_all_slopes_chosen_silva_tax |>
  group_by(order_f) |>
  dplyr::filter(n() >= 10) |>
  group_by(phylum_f, treatment) |>
  mutate(counts = n()) |>
  as_tibble() |>
  group_by(phylum_f, counts, treatment) |>
  dplyr::summarize() |>
  as_tibble()

#write.table(counts_divided, 'results/tables/counts_divided_asv_num_corrected.txt', sep = '\t')

reg_all_slopes_chosen_silva_tax_asv_reorder$phylum_f <- reg_all_slopes_chosen_silva_tax_asv_reorder$phylum_f |>
  factor( levels = c("Bacteroidota", "Proteobacteria", "Actinobacteriota",   "Cyanobacteria", "Planctomycetota",
                     "Verrucomicrobiota",  "Firmicutes", "Crenarchaeota", "Campilobacterota", "Bdellovibrionota",
                     "Acidobacteriota",
                     "Nitrospinota",  "Nitrospirota",  "Myxococcota", "Desulfobacterota",
                     "Deinococcota" , "Chloroflexi" , "Fusobacteriota", 
                     "Spirochaetota",  "Abditibacteriota", "Latescibacterota",
                     "Methylomirabilota", "Halanaerobiaeota", "Sumerlaeota", "Calditrichota", "Gemmatimonadota"), 
          ordered = TRUE)

summary_gr_families_treatment <- reg_all_slopes_chosen_silva_tax_asv_reorder |>
  group_by(family_f) |>
  dplyr::filter(n() >= 8) |>
  group_by(treatment, phylum ) |>
  dplyr::summarize(n())

ridge_family_divided <- reg_all_slopes_chosen_silva_tax_asv_reorder |>
  group_by(family_f) |>
  dplyr::filter(n() >= 8) |>
  ggplot(aes(y = fct_rev(phylum_f), x = slope_chosen_days, label = counts))+
  geom_density_ridges(aes(fill = phylum_f, group = family_f), scale = 1.0, alpha = 0.6,
                      jittered_points = TRUE,
                      point_shape = 21, point_size = 0.2, point_alpha = 0.0,
                      quantile_lines = TRUE,
                      quantile_fun = mean)+
  scale_fill_manual(values = palette_phylums_assigned)+
  labs(y = '', fill= 'Phylum', x = expression("Growth rate (d"^"-1" *')'), title = 'Family')+
  scale_x_continuous(limits = c(0,11))+
  geom_text(nudge_x = 7.2, nudge_y = 0.35, check_overlap = TRUE, size = 3)+
  facet_grid(cols = vars(treatment), scales = 'free')+
  coord_cartesian(clip = "off") +
  theme_bw()+
  theme(legend.position = "none", strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 0),
        plot.margin= unit(c(1, 1, 0.2, 0.2), "lines"), panel.border = element_blank())

ggsave('distribution_gr_plots_ridge_family_treatment_divided.pdf', ridge_family_divided,
       path = "~/Documentos/Doctorat/REMEI/results/figures/corrected_asv_num/",
       width = 188,
       height = 150,
       units = 'mm')

## Distribution plots for Proteobacteria and Bacteroidota separated (figure 4 and S4) ----
ridge_proteo <- 
  growthrates.distribution.tax.rank.ridges.divided(
    data = reg_all_slopes_chosen_silva_tax_asv_reorder, 
    phylum_to_explore = 'Proteobacteria',
    title = 'Proteobacteria class',
    axis_y_title = '',
    x_c = 10,
    x_o = 10,
    x_f = 10,
    x_a = 8) |>
  as_ggplot()

ggsave('ridge_proteobacteria_filtered_v6.pdf', ridge_proteo,
       path = "~/Documentos/Doctorat/REMEI/results/figures/corrected_asv_num/",
       width = 230,
       height = 180,
       units = 'mm')

ridge_bacterio <- growthrates.distribution.tax.rank.ridges.divided(
  data = reg_all_slopes_chosen_silva_tax, 
  phylum_to_explore = 'Bacteroidota',
  title = 'Bacteroidota \n
  Class',
  axis_y_title = '',
  x_c = 4,
  x_o = 10,
  x_f = 10,
  x_a = 8) |>
  as_ggplot()

# ggsave('ridge_bacteroidota_filtered_v4.pdf', ridge_bacterio,
#        path = "~/Documentos/Doctorat/REMEI/results/figures/corrected_asv_num/",
#        width = 180,
#        height = 180,
#        units = 'mm')

## Distributions for the whole growth rates calculated (figure 8) ----

counts_treatments <-   reg_all_slopes_chosen_silva_tax |>
  group_by(domain) |>
  dplyr::filter(n() >= 8) |>
  group_by(treatment) |>
  mutate(counts = n()) |>
  as_tibble() |>
  group_by(treatment, counts) |>
  dplyr::summarize() |>
  as_tibble()

#write.csv2(counts_treatments, 'results/tables/counts_treatments.csv')

counts_season <-   reg_all_slopes_chosen_silva_tax |>
  group_by(domain) |>
  dplyr::filter(n() >= 8) |>
  group_by(season) |>
  mutate(counts = n()) |>
  as_tibble() |>
  group_by(season, counts) |>
  dplyr::summarize() |>
  as_tibble()

#write.csv2(counts_season, 'results/tables/counts_season.csv')

distribution_gr_domains <-
  reg_all_slopes_chosen_silva_tax |>
  group_by(domain) |>
  dplyr::filter(n() >= 2) |>
  mutate(counts = paste('n = ', n())) |> #
  ggplot(aes(y = domain, x = slope_chosen_days, label = counts))+
  geom_density_ridges(alpha = 0.8, panel_scaling = TRUE, scale = 3,
                      jittered_points = TRUE,
                      point_shape = 21, point_size = 0.2, point_alpha = 0.0,
                      quantile_lines = TRUE,
                      quantile_fun = mean)+
  scale_x_continuous(limits = c(0,11))+
  coord_cartesian(clip = "off") +
  guides(color=guide_legend(ncol = 2, size = 5))+
  geom_text(nudge_x = 8, nudge_y = 0.3, check_overlap = TRUE, size = 1)+
  facet_grid(rows = vars(domain), scales = 'free', margins = NULL)+
  labs(y = '', x = expression("Growth rate (d"^"-1" *')'), fill = 'Treatment')+
  theme_bw()+
  theme(legend.position = "bottom", strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0, size = 7),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 7),
        plot.margin= unit(c(1, 1, 0.2, 0.2), "lines"), panel.border = element_blank(), legend.text = element_text(size = 5),
        legend.title = element_text(size = 7), strip.background.y = element_blank(), strip.text.y = element_blank())

# ggsave('distribution_gr_domains_v2.pdf', distribution_gr_domains,
#        path = "~/Documentos/Doctorat/REMEI/results/figures/corrected_asv_num/",
#        width = 88,
#        height = 88,
#        units = 'mm')

distribution_gr_domains_seas <-
  reg_all_slopes_chosen_silva_tax |>
  group_by(domain) |>
  dplyr::filter(n() >= 8) |>
  mutate(counts = paste('n = ', n())) |> #
  ggplot(aes(y = season, x = slope_chosen_days, label = counts, fill = season))+
  geom_density_ridges(alpha = 0.8, panel_scaling = TRUE, scale = 3,
                      aes(group = season),
                      jittered_points = TRUE,
                      point_shape = 21, point_size = 0.2, point_alpha = 0.0,
                      quantile_lines = TRUE,
                      quantile_fun = mean)+
  scale_x_continuous(limits = c(0,11))+
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = palette_seasons_4)+
  guides(color=guide_legend(ncol = 2, size = 5))+
  geom_text(nudge_x = 8, nudge_y = 0.3, check_overlap = TRUE, size = 1)+
  facet_grid(rows = vars(season), scales = 'free', margins = NULL)+
  labs(y = '', x = expression("Growth rate (d"^"-1" *')'), fill = 'Season')+
  theme_bw()+
  theme(legend.position = "none", strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0, size = 7),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 7),
        plot.margin= unit(c(1, 1, 0.2, 0.2), "lines"), panel.border = element_blank(), legend.text = element_text(size = 5),
        legend.title = element_text(size = 7), strip.background.y = element_blank(), strip.text.y = element_blank())

ggsave('distribution_gr_domains_seasons_v2.pdf', distribution_gr_domains_seas,
       path = "~/Documentos/Doctorat/REMEI/results/figures/corrected_asv_num/",
       width = 88,
       height = 88,
       units = 'mm')

distribution_gr_domains_treat <-
  reg_all_slopes_chosen_silva_tax |>
  group_by(domain) |>
  dplyr::filter(n() >= 8) |>
  mutate(counts = paste('n = ', n())) |> #
  ggplot(aes(y = treatment, x = slope_chosen_days, label = counts, fill = treatment))+
  geom_density_ridges(alpha = 0.8, panel_scaling = TRUE, scale = 3,
                      aes(group = treatment),
                      jittered_points = TRUE,
                      point_shape = 21, point_size = 0.2, point_alpha = 0.0,
                      quantile_lines = TRUE,
                      quantile_fun = mean)+
  scale_x_continuous(limits = c(0,11))+
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = palette_treatments_remei)+
  guides(color=guide_legend(ncol = 2, size = 5))+
  geom_text(nudge_x = 8, nudge_y = 0.3, check_overlap = TRUE, size = 1)+
  facet_grid(rows = vars(treatment), scales = 'free', margins = NULL)+
  labs(y = '', x = expression("Growth rate (d"^"-1" *')'), fill = 'Treatment')+
  theme_bw()+
  theme(legend.position = "none", strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0, size = 7),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 7),
        plot.margin= unit(c(1, 1, 0.2, 0.2), "lines"), panel.border = element_blank(), legend.text = element_text(size = 5),
        legend.title = element_text(size = 7), strip.background.y = element_blank(), strip.text.y = element_blank())

ggsave('distribution_gr_domains_treatments_v2.pdf', distribution_gr_domains_treat,
       path = "~/Documentos/Doctorat/REMEI/results/figures/corrected_asv_num/",
       width = 88,
       height = 88,
       units = 'mm')

distribution_gr_domains_all <- multi_panel_figure(columns = 3, rows = 1, width = 188, height = 90, 
                                                  column_spacing = 0, unit = 'mm', row_spacing = 0.6,
                                                  panel_label_type = 'upper-alpha')

distribution_gr_domains_all  %<>%
  fill_panel(distribution_gr_domains, column = 1, row = 1) %<>%
  fill_panel(distribution_gr_domains_seas, column = 2, row = 1) %<>%
  fill_panel(distribution_gr_domains_treat, column = 3, row = 1) 

# ggsave('distribution_gr_domains_all.pdf', distribution_gr_domains_all,
#        path = "~/Documentos/Doctorat/REMEI/results/figures/corrected_asv_num/",
#        width = 188,
#        height = 100,
#        units = 'mm')

## Distribution for maximal growth rates only----
## Distribution plots at different taxonomic ranks ----
reg_all_slopes_chosen_silva_tax_max <- reg_all_slopes_chosen_silva_tax  |>
  group_by(asv_num) |>
  dplyr::filter(n() > 2) |>
  ungroup() |>
  slice_max(order_by = slope_chosen_days, n = 1, by = asv_num) |>
  mutate(phylum_f = fct_rev(fct_infreq(phylum_f)))

counts <-   reg_all_slopes_chosen_silva_tax |> #number of growth rates calculated by phylum
  group_by(asv_num) |>
  dplyr::filter(n() > 2) |>
  ungroup() |>
  slice_max(order_by = slope_chosen_days, n = 1, by = asv_num) |>
  group_by(phylum_f) |>
  dplyr::filter(n() >= 2) |>
  mutate(counts = n()) |>
  as_tibble() |>
  group_by(phylum_f, counts) |>
  dplyr::summarize() |>
  as_tibble()

ridge_ph <- # distribution plot at phylum taxonomic rank
  reg_all_slopes_chosen_silva_tax |>
  group_by(asv_num) |>
  dplyr::filter(n() >= 2) |>
  ungroup() |>
  slice_max(order_by = slope_chosen_days, n = 1, by = asv_num) |>
  group_by(phylum_f) |>
  dplyr::filter(n() >= 50) |>
  mutate(counts = paste('n = ', n())) |> #
  ggplot(aes(y = fct_rev(fct_infreq(phylum_f)), x = slope_chosen_days, fill = phylum_f,  group = phylum_f, label = counts))+
  geom_density_ridges(alpha = 0.8, panel_scaling = TRUE, scale = 1,
                      jittered_points = TRUE,
                      point_shape = 21, point_size = 0.2, point_alpha = 0.0,
                      quantile_lines = TRUE,
                      quantile_fun = mean)+
  scale_fill_manual(values = palette_phylums_assigned)+
  scale_x_continuous(limits = c(0,11), expand = c(0,0))+
  coord_cartesian(clip = "off") +
  geom_text(nudge_x = 8.8, nudge_y = 0.3, check_overlap = TRUE, size = 3)+
  labs(y = '', fill= 'Phylum', x = expression("Growth rate (d"^"-1" *')'), title = 'Phylum')+
  theme_bw()+
  theme(legend.position = "none", strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0, size = 6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 6),
        plot.margin= unit(c(0.2, 0.5, 0.2, 0.2), "lines"), panel.border = element_blank(), 
        title = element_text(size = 8))

ridge_cl <- reg_all_slopes_chosen_silva_tax |> # distribution plot at class taxonomic rank
  group_by(asv_num) |>
  dplyr::filter(n() >= 2) |>
  ungroup() |>
  slice_max(order_by = slope_chosen_days, n = 1, by = asv_num) |>
  group_by(phylum_f) |>
  dplyr::filter(n() >= 50) |>
  group_by(phylum_f) |>
  mutate(counts = paste('n = ', n())) |> #
  ggplot(aes(y = fct_rev(fct_infreq(phylum_f)), x = slope_chosen_days, fill = phylum_f, label = counts))+
  geom_density_ridges(aes(fill = phylum_f, group = class_f), scale = 1, alpha = 0.7,
                      jittered_points = TRUE,
                      point_shape = 21, point_size = 0.2, point_alpha = 0.0,
                      quantile_lines = TRUE,
                      quantile_fun = mean)+
  scale_fill_manual(values = palette_phylums_assigned)+
  labs(y = 'Class', fill= 'Phylum', x = expression("Growth rate (d"^"-1" *')'), title = 'Class')+
  scale_x_continuous(limits = c(0,11), expand = c(0,0))+
  geom_text(nudge_x = 7.8, nudge_y = 0.3, check_overlap = TRUE, size = 3)+
  theme_bw()+
  theme(legend.position = "none", strip.text.x = element_text(size = 6), 
        axis.text.x = element_text(angle = 0, size = 6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.y = element_text(size = 0),
        axis.title.y = element_blank(), 
        plot.margin= unit(c(0.2, 0.5, 0.2, 0.2), "lines"), panel.border = element_blank(), 
        title = element_text(size = 8))

ridge_o <- reg_all_slopes_chosen_silva_tax |> # distribution plot at order taxonomic rank
  group_by(asv_num) |>
  dplyr::filter(n() >= 2) |>
  ungroup() |>
  slice_max(order_by = slope_chosen_days, n = 1, by = asv_num) |>
  group_by(phylum_f) |>
  dplyr::filter(n() >= 50) |>
  ungroup() |>
  group_by(phylum_f) |>
  mutate(counts = paste('n = ', n())) |> #
  ggplot(aes(y = fct_rev(fct_infreq(phylum_f)), x = slope_chosen_days, fill = phylum_f, label = counts))+
  geom_density_ridges(aes(fill = phylum_f, group = order_f), scale = 0.8, alpha = 0.6,
                      jittered_points = TRUE,
                      point_shape = 21, point_size = 0.2, point_alpha = 0.0,
                      quantile_lines = TRUE,
                      quantile_fun = mean)+
  scale_fill_manual(values = palette_phylums_assigned)+
  labs(y = '', fill= 'Phylum', x = expression("Growth rate (d"^"-1" *')'), title = 'Order')+
  scale_x_continuous(limits = c(0,11), expand = c(0,0))+
  geom_text(nudge_x = 7.8, nudge_y = 0.3, check_overlap = TRUE, size = 3)+
  theme_bw()+
  theme(legend.position = "none", strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0, size = 6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 0),
        plot.margin= unit(c(0.2, 0.5, 0.2, 0.2), "lines"), panel.border = element_blank(), 
        title = element_text(size = 8))

ridge_f <- reg_all_slopes_chosen_silva_tax |> # distribution plot at order taxonomic rank
  group_by(asv_num) |>
  dplyr::filter(n() >= 2) |>
  ungroup() |>
  slice_max(order_by = slope_chosen_days, n = 1, by = asv_num) |>
  group_by(phylum_f) |>
  dplyr::filter(n() >= 50) |>
  group_by(phylum_f) |>
  mutate(counts = paste('n = ', n())) |>
  ungroup() |>
  ggplot(aes(y = fct_rev(fct_infreq(phylum_f)), x = slope_chosen_days, fill = phylum_f, label = counts))+
  geom_density_ridges(aes(fill = phylum_f, group = family_f), scale = 0.8, alpha = 0.6,
                      jittered_points = TRUE,
                      point_shape = 21, point_size = 0.2, point_alpha = 0.0,
                      quantile_lines = TRUE,
                      quantile_fun = mean)+
  scale_fill_manual(values = palette_phylums_assigned)+
  labs(y = '', fill= 'Phylum', x = expression("Growth rate (d"^"-1" *')'), title = 'Family')+
  scale_x_continuous(limits = c(0,11), expand = c(0,0))+
  geom_text(nudge_x = 7.8, nudge_y = 0.3, check_overlap = TRUE, size = 3)+
  theme_bw()+
  theme(legend.position = "none", strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0, size = 6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 0),
        plot.margin= unit(c(0.2, 0, 0.2, 0.2), "lines"), panel.border = element_blank(), 
        title = element_text(size = 8))

## reorder ASVs by phylum
reg_all_slopes_chosen_silva_tax_asv_reorder <- reg_all_slopes_chosen_silva_tax |>
  group_by(asv_num) |>
  dplyr::filter(n() >= 2) |>
  ungroup() |>
  slice_max(order_by = slope_chosen_days, n = 1, by = asv_num) |>
  group_by(phylum_f) |>
  mutate(counts = paste('n = ', n())) |> #
  ungroup()

reg_all_slopes_chosen_silva_tax_asv_reorder$phylum_f <- reg_all_slopes_chosen_silva_tax_asv_reorder$phylum_f |>
  factor( levels = c("Bacteroidota", "Proteobacteria", "Actinobacteriota",   "Cyanobacteria", "Planctomycetota",
                     "Verrucomicrobiota",  "Firmicutes", "Crenarchaeota", "Bdellovibrionota", "Campilobacterota", 
                     "Acidobacteriota",
                     "Nitrospinota",  "Nitrospirota",  "Myxococcota", "Desulfobacterota",
                     "Deinococcota" , "Chloroflexi" , "Fusobacteriota", 
                     "Spirochaetota",  "Abditibacteriota", "Latescibacterota",
                     "Methylomirabilota", "Halanaerobiaeota", "Sumerlaeota", "Calditrichota", "Gemmatimonadota"), 
          ordered = TRUE)

ridge_a <- reg_all_slopes_chosen_silva_tax_asv_reorder |> # distribution plot at ASV taxonomic rank
  group_by(asv_num) |>
  dplyr::filter(n() >= 2) |>
  ungroup() |>
  slice_max(order_by = slope_chosen_days, n = 1, by = asv_num) |>
  group_by(phylum_f) |>
  dplyr::filter(n() >= 50) |>
  ggplot(aes(y = fct_rev(phylum_f), x = slope_chosen_days, label = counts))+
  geom_density_ridges(aes(fill = phylum_f, group = asv_f), scale = 0.6, alpha = 0.6,
                      jittered_points = TRUE,
                      point_shape = 21, point_size = 0.2, point_alpha = 0.0,
                      quantile_lines = TRUE,
                      quantile_fun = mean)+
  scale_fill_manual(values = palette_phylums_assigned)+
  labs(y = '', fill= 'Phylum', x = expression("Growth rate (d"^"-1" *')'), title = 'ASV')+
  scale_x_continuous(limits = c(0,11), expand = c(0,0))+
  geom_text(nudge_x = 7.2, nudge_y = 0.35, check_overlap = TRUE, size = 3)+
  coord_cartesian(clip = "off") +
  theme_bw()+
  theme(legend.position = "none", strip.text.x = element_text(size = 6), 
        axis.text.x = element_text(angle = 0, size = 6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 0),
        plot.margin= unit(c(0.2, 0.5, 0.2, 0.2), "lines"), panel.border = element_blank(), 
        title = element_text(size = 8))

## Percentage of significant growth rates calculated by Phylum ------ 
perc_gr_phylum <- reg_all_slopes_chosen_silva_tax |>
  group_by(asv_num) |>
  dplyr::filter(n() >= 2) |>
  ungroup() |>
  slice_max(order_by = slope_chosen_days, n = 1, by = asv_num) |>
  group_by(phylum_f) |>
  dplyr::filter(n() > 2) |>
  group_by(phylum_f, slope_chosen_days) |> 
  dplyr::summarize(n = n()) |>
  group_by(phylum_f) |>
  dplyr::summarize(n = n()) |>
  mutate(sum = sum(n),
         perc = n/sum) |>
  ungroup() |>
  ggplot(aes(sum, perc, fill = phylum_f))+
  geom_col(aes(fill = fct_reorder(phylum_f, perc, .desc = TRUE)))+ 
  scale_y_continuous(labels = percent_format())+
  labs(y = 'Percentage of maximal significant growth rates calculated per Phylum', fill = 'Phylum')+
  scale_fill_manual(values = palette_phylums_assigned)+
  coord_flip()+
  theme_bw()+
  theme(axis.text.y = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), aspect.ratio = 2/10, panel.border = element_blank(),
        axis.ticks = element_blank(), legend.position = 'none', axis.text.x = element_text(size = 6),
        axis.title.y = element_blank(), axis.text = element_text(size = 6), text = element_text(size = 6),
        plot.margin= unit(c(1, 0.5, 1, 0.5), "lines"))

## Percentage of significant growth rates calculated at different taxonomic ranks ----
x <- reg_all_slopes_chosen_silva_tax |>
  group_by(asv_num) |>
  dplyr::filter(n() > 2) |>
  ungroup() |>
  slice_max(order_by = slope_chosen_days, n = 1, by = asv_num) |>
  group_by(phylum_f) |>
  dplyr::filter(n() > 2) |>
  dplyr::summarize(n = n()) |>
  dplyr::summarize(n = sum(n))

gr_asv_perc <- reg_all_slopes_chosen_silva_tax_filt |>
  group_by(asv_num) |>
  dplyr::filter(n() > 2) |>
  ungroup() |>
  slice_max(order_by = slope_chosen_days, n = 1, by = asv_num) |>
  group_by(asv_f) |>
  dplyr::filter(n() > 2) |>
  dplyr::summarize(n = n()) |>
  mutate(sum = sum(n),
         asv_gr_perc = sum/x) |>
  dplyr::select(asv_gr_perc) |>
  unique() |>
  as_tibble()

gr_f_perc <- reg_all_slopes_chosen_silva_tax |>
  group_by(asv_num) |>
  dplyr::filter(n() > 2) |>
  ungroup() |>
  slice_max(order_by = slope_chosen_days, n = 1, by = asv_num) |>
  group_by(family) |>
  dplyr::filter(n() > 2) |>
  dplyr::summarize(n = n()) |>
  mutate(sum = sum(n),
         f_gr_perc = sum/x) |>
  dplyr::select(f_gr_perc) |>
  unique() |>
  as_tibble()

gr_o_perc <- reg_all_slopes_chosen_silva_tax |>
  group_by(asv_num) |>
  dplyr::filter(n() > 2) |>
  ungroup() |>
  slice_max(order_by = slope_chosen_days, n = 1, by = asv_num) |>
  group_by(order) |>
  dplyr::filter(n() > 2) |>
  dplyr::summarize(n = n()) |>
  mutate(sum = sum(n),
         order_gr_perc = sum/x) |>
  dplyr::select(order_gr_perc) |>
  unique() |>
  as_tibble()

gr_c_perc <- reg_all_slopes_chosen_silva_tax |>
  group_by(asv_num) |>
  dplyr::filter(n() > 2) |>
  ungroup() |>
  slice_max(order_by = slope_chosen_days, n = 1, by = asv_num) |>
  group_by(class) |>
  dplyr::filter(n() > 2) |>
  dplyr::summarize(n = n()) |>
  mutate(sum = sum(n),
         class_gr_perc = sum/x) |>
  dplyr::select(class_gr_perc) |>
  unique() |>
  as_tibble()

gr_p_perc <- reg_all_slopes_chosen_silva_tax |>
  group_by(asv_num) |>
  dplyr::filter(n() > 2) |>
  ungroup() |>
  slice_max(order_by = slope_chosen_days, n = 1, by = asv_num) |>
  group_by(phylum) |>
  dplyr::filter(n() > 2) |>
  dplyr::summarize(n = n()) |>
  mutate(sum = sum(n),
         phylum_gr_perc = sum/x) |>
  dplyr::select(phylum_gr_perc) |>
  unique() |>
  as_tibble() 

labels_perc_gr_rank <-  as_labeller(c(phylum_gr_perc = 'Phylum',
                                      class_gr_perc = 'Class',
                                      order_gr_perc = 'Order',
                                      f_gr_perc = 'Family',
                                      asv_gr_perc = 'ASV'))

# gr_tax_ranks_perc <- cbind( gr_f_perc, gr_o_perc, gr_c_perc) |>
#   as_tibble() |>
#   pivot_longer(cols = 1:5) |>
#   ggplot(aes(x = fct_rev(as_factor(name)), y = value$n, group = 1))+
#   geom_point()+
#   #coord_flip()+
#   geom_line()+
#   labs(y='Percentage', x = 'Growth rates distribution plotted\nat different taxonomic ranks')+
#   scale_x_discrete(labels = labels_perc_gr_rank)+
#   scale_y_continuous(labels = percent_format())+
#   theme_bw()+
#   theme(aspect.ratio = 3/10, panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(),
#         legend.position = 'none', axis.text.y = element_text(size = 6),
#         axis.title.y = element_text(size = 6), text = element_text(size = 6),
#         plot.margin= unit(c(1, 0.5, 1, 0.5), "lines"))

## mean growth rates at different taxonomic levels ----
labels_mean_rank <- as_labeller(c( mean_phylum = 'Phylum',
                                   mean_class = 'Class',
                                   mean_order = 'Order',
                                   mean_family = 'Family',
                                   mean_asv = 'ASV'))

mean_gr_ranks_max <- reg_all_slopes_chosen_silva_tax |>
  group_by(asv_num) |>
  dplyr::filter(n() > 2) |>
  ungroup() |>
  slice_max(order_by = slope_chosen_days, n = 1, by = asv_num) |>
  group_by(phylum_f) |>
  dplyr::filter(n() > 2) |>
  mutate(mean_phylum = mean(slope_chosen_days)) |>
  ungroup() |>
  group_by(class) |>
  mutate(mean_class = mean(slope_chosen_days)) |>
  ungroup() |>
  group_by(order) |>
  mutate(mean_order = mean(slope_chosen_days)) |>
  ungroup() |>
  group_by(family) |>
  mutate(mean_family = mean(slope_chosen_days)) |>
  ungroup() |>
  group_by(asv_num) |>
  mutate(mean_asv = mean(slope_chosen_days)) |>
  dplyr::select(starts_with('mean'), 'asv_num', 'family', 'order', 'class', 'phylum_f') |>
  pivot_longer(cols = starts_with('mean')) |>
  group_by(phylum_f) |>
  distinct(value, name, phylum_f) |>
  ggplot(aes(as_factor(name), value))+
  geom_point(aes(color = phylum_f), alpha = 0.6, 
             position = position_jitter(0.3))+
  geom_violin(alpha = 0.2, draw_quantiles = c(0.25, 0.75))+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.6,
               colour = "black")+
  scale_color_manual(values = palette_phylums_assigned)+
  labs(y = 'Mean growth rate\nat different\ntaxonomic ranks', x= 'Taxonomic rank', color = 'Phylum')+
  scale_x_discrete(labels = labels_mean_rank)+
  #coord_flip()+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = 'none', axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 6), text = element_text(size = 6),
        plot.margin= unit(c(1, 0.5, 1, 0.5), "lines"))
# 
# ggsave('mean_gr_ranks_ed.pdf', mean_gr_ranks,
# path = "~/Documentos/Doctorat/REMEI/results/figures/corrected_asv_num/",
# width = 150,
# height = 130,
# units = 'mm')


## Figure 3 creating panel design------
distribution_gr_plots_ridge_max <- multi_panel_figure(columns = 5, rows = 1, width = 188, height = 150, 
                                                  column_spacing = 0.2, unit = 'mm', row_spacing = 0.6,
                                                  panel_label_type = 'upper-alpha')

# distribution_gr_plots_ridge  %<>%
#   fill_panel(perc_gr_phylum, column = 1, row = 1) %<>%
#   fill_panel(gr_tax_ranks_perc, column = 1, row = 2) %<>%
#   fill_panel(mean_gr_ranks, column = 1, row = 3) %<>%
#   fill_panel(ridge_ph, column = 2:3, row = 1:3) %<>%
#   fill_panel(ridge_cl, column = 4, row = 1:3) %<>%
#   fill_panel(ridge_o, column = 5, row = 1:3) %<>%
#   fill_panel(ridge_f, column = 6, row = 1:3) %<>%
#   fill_panel(ridge_a, column = 7, row = 1:3)

distribution_gr_plots_ridge_max  %<>%
  # fill_panel(perc_gr_phylum, column = 1:2, row = 1) %<>%
  # fill_panel(gr_tax_ranks_perc, column = 3:4, row = 1) %<>%
  # fill_panel(mean_gr_ranks_max, column = 5:6, row = 1) %<>%
  fill_panel(ridge_ph, column = 1:2, row = 1) %<>%
  fill_panel(ridge_cl, column = 3, row = 1) %<>%
  fill_panel(ridge_o, column = 4, row = 1) %<>%
  fill_panel(ridge_f, column = 5, row = 1) #%<>%
  # fill_panel(perc_gr_phylum, column = 1:2, row = 1) |>
  # fill_panel(mean_gr_ranks_max, column = 4:5, row = 1)

# save the graph
ggsave('distribution_gr_plots_ridge_order_freq_perc_max_ed2.pdf', distribution_gr_plots_ridge_max,
       path = "~/Documentos/Doctorat/REMEI/results/figures/corrected_asv_num/",
       width = 188,
       height = 150,
       units = 'mm')

## exploration at family level maximal growth rates
ridge_f <- reg_all_slopes_chosen_silva_tax |> # distribution plot at order taxonomic rank
  group_by(asv_num) |>
  dplyr::filter(n() >= 2) |>
  ungroup() |>
  slice_max(order_by = slope_chosen_days, n = 1, by = asv_num) |>
  group_by(family_f) |>
  dplyr::filter(n() >= 4) |>
  group_by(family_f) |>
  mutate(counts = paste('n = ', n())) |>
  ungroup() |>
  ggplot(aes(y = fct_rev(family_f), x = slope_chosen_days, fill = phylum_f, label = counts))+
  geom_density_ridges(aes(fill = class_f, group = family_f), scale = 1, alpha = 0.6,
                      jittered_points = TRUE,
                      point_shape = 21, point_size = 0.2, point_alpha = 0.0,
                      quantile_lines = TRUE,
                      quantile_fun = mean)+
  scale_fill_manual(values = palette_class_assigned)+
  labs(y = '', fill= 'Phylum', x = expression("Growth rate (d"^"-1" *')'), title = 'Family')+
  scale_x_continuous(limits = c(0,11), expand = c(0,0))+
  geom_text(nudge_x = 7.8, nudge_y = 0.3, check_overlap = TRUE, size = 3)+
  theme_bw()+
  theme(legend.position = "none", strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0, size = 6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 6),
        plot.margin= unit(c(0.2, 0, 0.2, 0.2), "lines"), panel.border = element_blank(), 
        title = element_text(size = 8))

ggsave('distribution_max_gr_family.pdf', ridge_f,
       path = "~/Documentos/Doctorat/REMEI/results/figures/corrected_asv_num/",
       width = 88,
       height = 180,
       units = 'mm')


ridge_f_winter <- reg_all_slopes_chosen_silva_tax |> # distribution plot at order taxonomic rank
  group_by(asv_num) |>
  dplyr::filter(n() >= 2) |>
  ungroup() |>
  slice_max(order_by = slope_chosen_days, n = 1, by = asv_num) |>
  dplyr::filter(season == 'Winter') |>
  group_by(family_f) |>
  dplyr::filter(n() >= 4) |>
  group_by(family_f) |>
  mutate(counts = paste('n = ', n())) |>
  ungroup() |>
  ggplot(aes(y = fct_rev(family_f), x = slope_chosen_days, fill = phylum_f, label = counts))+
  geom_density_ridges(aes(fill = class_f, group = family_f), scale = 0.8, alpha = 0.6,
                      jittered_points = TRUE,
                      point_shape = 21, point_size = 0.2, point_alpha = 0.0,
                      quantile_lines = TRUE,
                      quantile_fun = mean)+
  scale_fill_manual(values = palette_class_assigned)+
  labs(y = '', fill= 'Phylum', x = expression("Growth rate (d"^"-1" *')'), title = 'Family Winter')+
  scale_x_continuous(limits = c(0,11), expand = c(0,0))+
  geom_text(nudge_x = 7.8, nudge_y = 0.3, check_overlap = TRUE, size = 3)+
  theme_bw()+
  theme(legend.position = "none", strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0, size = 6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 6),
        plot.margin= unit(c(0.2, 0, 0.2, 0.2), "lines"), panel.border = element_blank(), 
        title = element_text(size = 8))

ridge_f_spring <- reg_all_slopes_chosen_silva_tax |> # distribution plot at order taxonomic rank
  group_by(asv_num) |>
  dplyr::filter(n() >= 2) |>
  ungroup() |>
  slice_max(order_by = slope_chosen_days, n = 1, by = asv_num) |>
  dplyr::filter(season == 'Spring') |>
  group_by(family_f) |>
  dplyr::filter(n() >= 4) |>
  group_by(family_f) |>
  mutate(counts = paste('n = ', n())) |>
  ungroup() |>
  ggplot(aes(y = fct_rev(family_f), x = slope_chosen_days, fill = phylum_f, label = counts))+
  geom_density_ridges(aes(fill = class_f, group = family_f), scale = 0.8, alpha = 0.6,
                      jittered_points = TRUE,
                      point_shape = 21, point_size = 0.2, point_alpha = 0.0,
                      quantile_lines = TRUE,
                      quantile_fun = mean)+
  scale_fill_manual(values = palette_class_assigned)+
  labs(y = '', fill= 'Phylum', x = expression("Growth rate (d"^"-1" *')'), title = 'Family Spring')+
  scale_x_continuous(limits = c(0,11), expand = c(0,0))+
  geom_text(nudge_x = 7.8, nudge_y = 0.3, check_overlap = TRUE, size = 3)+
  theme_bw()+
  theme(legend.position = "none", strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0, size = 6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 6),
        plot.margin= unit(c(0.2, 0, 0.2, 0.2), "lines"), panel.border = element_blank(), 
        title = element_text(size = 8))

ridge_f_summer <- reg_all_slopes_chosen_silva_tax |> # distribution plot at order taxonomic rank
  group_by(asv_num) |>
  dplyr::filter(n() >= 2) |>
  ungroup() |>
  slice_max(order_by = slope_chosen_days, n = 1, by = asv_num) |>
  dplyr::filter(season == 'Summer') |>
  group_by(family_f) |>
  dplyr::filter(n() >= 4) |>
  group_by(family_f) |>
  mutate(counts = paste('n = ', n())) |>
  ungroup() |>
  ggplot(aes(y = fct_rev(family_f), x = slope_chosen_days, fill = phylum_f, label = counts))+
  geom_density_ridges(aes(fill = class_f, group = family_f), scale = 0.8, alpha = 0.6,
                      jittered_points = TRUE,
                      point_shape = 21, point_size = 0.2, point_alpha = 0.0,
                      quantile_lines = TRUE,
                      quantile_fun = mean)+
  scale_fill_manual(values = palette_class_assigned)+
  labs(y = '', fill= 'Phylum', x = expression("Growth rate (d"^"-1" *')'), title = 'Family Summer')+
  scale_x_continuous(limits = c(0,11), expand = c(0,0))+
  geom_text(nudge_x = 7.8, nudge_y = 0.3, check_overlap = TRUE, size = 3)+
  theme_bw()+
  theme(legend.position = "none", strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0, size = 6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 6),
        plot.margin= unit(c(0.2, 0, 0.2, 0.2), "lines"), panel.border = element_blank(), 
        title = element_text(size = 8))

ridge_f_fall <- reg_all_slopes_chosen_silva_tax |> # distribution plot at order taxonomic rank
  group_by(asv_num) |>
  dplyr::filter(n() >= 2) |>
  ungroup() |>
  slice_max(order_by = slope_chosen_days, n = 1, by = asv_num) |>
  dplyr::filter(season == 'Fall') |>
  group_by(family_f) |>
  dplyr::filter(n() >= 4) |>
  group_by(family_f) |>
  mutate(counts = paste('n = ', n())) |>
  ungroup() |>
  ggplot(aes(y = fct_rev(family_f), x = slope_chosen_days, fill = class_f, label = counts))+
  geom_density_ridges(aes(fill = class_f, group = family_f), scale = 0.8, alpha = 0.6,
                      jittered_points = TRUE,
                      point_shape = 21, point_size = 0.2, point_alpha = 0.0,
                      quantile_lines = TRUE,
                      quantile_fun = mean)+
  scale_fill_manual(values = palette_class_assigned)+
  labs(y = '', fill= 'Phylum', x = expression("Growth rate (d"^"-1" *')'), title = 'Family Fall')+
  scale_x_continuous(limits = c(0,11), expand = c(0,0))+
  geom_text(nudge_x = 7.8, nudge_y = 0.3, check_overlap = TRUE, size = 3)+
  theme_bw()+
  theme(legend.position = "none", strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0, size = 6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 6),
        plot.margin= unit(c(0.2, 0, 0.2, 0.2), "lines"), panel.border = element_blank(), 
        title = element_text(size = 8))

ggarrange(ridge_f_winter, ridge_f_spring, ridge_f_summer, ridge_f_fall)

### treatments ---
ridge_f_treatments_cd <- reg_all_slopes_chosen_silva_tax |> # distribution plot at order taxonomic rank
  group_by(asv_num) |>
  dplyr::filter(n() >= 2) |>
  ungroup() |>
  slice_max(order_by = slope_chosen_days, n = 1, by = asv_num) |>
  dplyr::filter(treatment == 'CD') |>
  group_by(family_f) |>
  dplyr::filter(n() >= 4) |>
  group_by(family_f) |>
  mutate(counts = paste('n = ', n())) |>
  ungroup() |>
  ggplot(aes(y = fct_rev(family_f), x = slope_chosen_days, fill = phylum_f, label = counts))+
  geom_density_ridges(aes(fill = class_f, group = family_f), scale = 1, alpha = 0.6,
                      jittered_points = TRUE,
                      point_shape = 21, point_size = 0.2, point_alpha = 0.0,
                      quantile_lines = TRUE,
                      quantile_fun = mean)+
  scale_fill_manual(values = palette_class_assigned)+
  labs(y = '', fill= 'Phylum', x = expression("Growth rate (d"^"-1" *')'), title = 'Family')+
  scale_x_continuous(limits = c(0,11), expand = c(0,0))+
  geom_text(nudge_x = 7.8, nudge_y = 0.3, check_overlap = TRUE, size = 3)+
  theme_bw()+
  theme(legend.position = "none", strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0, size = 6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 6),
        plot.margin= unit(c(0.2, 0, 0.2, 0.2), "lines"), panel.border = element_blank(), 
        title = element_text(size = 8))

ridge_f_treatments_cl <- reg_all_slopes_chosen_silva_tax |> # distribution plot at order taxonomic rank
  group_by(asv_num) |>
  dplyr::filter(n() >= 2) |>
  ungroup() |>
  slice_max(order_by = slope_chosen_days, n = 1, by = asv_num) |>
  dplyr::filter(treatment == 'CL') |>
  group_by(family_f) |>
  dplyr::filter(n() >= 4) |>
  group_by(family_f) |>
  mutate(counts = paste('n = ', n())) |>
  ungroup() |>
  ggplot(aes(y = fct_rev(family_f), x = slope_chosen_days, fill = phylum_f, label = counts))+
  geom_density_ridges(aes(fill = class_f, group = family_f), scale = 1, alpha = 0.6,
                      jittered_points = TRUE,
                      point_shape = 21, point_size = 0.2, point_alpha = 0.0,
                      quantile_lines = TRUE,
                      quantile_fun = mean)+
  scale_fill_manual(values = palette_class_assigned)+
  labs(y = '', fill= 'Phylum', x = expression("Growth rate (d"^"-1" *')'), title = 'Family')+
  scale_x_continuous(limits = c(0,11), expand = c(0,0))+
  geom_text(nudge_x = 7.8, nudge_y = 0.3, check_overlap = TRUE, size = 3)+
  theme_bw()+
  theme(legend.position = "none", strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0, size = 6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 6),
        plot.margin= unit(c(0.2, 0, 0.2, 0.2), "lines"), panel.border = element_blank(), 
        title = element_text(size = 8))

ridge_f_treatments_pd <- reg_all_slopes_chosen_silva_tax |> # distribution plot at order taxonomic rank
  group_by(asv_num) |>
  dplyr::filter(n() >= 2) |>
  ungroup() |>
  slice_max(order_by = slope_chosen_days, n = 1, by = asv_num) |>
  dplyr::filter(treatment == 'PD') |>
  group_by(family_f) |>
  dplyr::filter(n() >= 4) |>
  group_by(family_f) |>
  mutate(counts = paste('n = ', n())) |>
  ungroup() |>
  ggplot(aes(y = fct_rev(family_f), x = slope_chosen_days, fill = phylum_f, label = counts))+
  geom_density_ridges(aes(fill = class_f, group = family_f), scale = 1, alpha = 0.6,
                      jittered_points = TRUE,
                      point_shape = 21, point_size = 0.2, point_alpha = 0.0,
                      quantile_lines = TRUE,
                      quantile_fun = mean)+
  scale_fill_manual(values = palette_class_assigned)+
  labs(y = '', fill= 'Phylum', x = expression("Growth rate (d"^"-1" *')'), title = 'Family')+
  scale_x_continuous(limits = c(0,11), expand = c(0,0))+
  geom_text(nudge_x = 7.8, nudge_y = 0.3, check_overlap = TRUE, size = 3)+
  theme_bw()+
  theme(legend.position = "none", strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0, size = 6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 6),
        plot.margin= unit(c(0.2, 0, 0.2, 0.2), "lines"), panel.border = element_blank(), 
        title = element_text(size = 8))

ridge_f_treatments_pl <- reg_all_slopes_chosen_silva_tax |> # distribution plot at order taxonomic rank
  group_by(asv_num) |>
  dplyr::filter(n() >= 2) |>
  ungroup() |>
  slice_max(order_by = slope_chosen_days, n = 1, by = asv_num) |>
  dplyr::filter(treatment == 'PL') |>
  group_by(family_f) |>
  dplyr::filter(n() >= 4) |>
  group_by(family_f) |>
  mutate(counts = paste('n = ', n())) |>
  ungroup() |>
  ggplot(aes(y = fct_rev(family_f), x = slope_chosen_days, fill = phylum_f, label = counts))+
  geom_density_ridges(aes(fill = class_f, group = family_f), scale = 1, alpha = 0.6,
                      jittered_points = TRUE,
                      point_shape = 21, point_size = 0.2, point_alpha = 0.0,
                      quantile_lines = TRUE,
                      quantile_fun = mean)+
  scale_fill_manual(values = palette_class_assigned)+
  labs(y = '', fill= 'Phylum', x = expression("Growth rate (d"^"-1" *')'), title = 'Family')+
  scale_x_continuous(limits = c(0,11), expand = c(0,0))+
  geom_text(nudge_x = 7.8, nudge_y = 0.3, check_overlap = TRUE, size = 3)+
  theme_bw()+
  theme(legend.position = "none", strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0, size = 6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 6),
        plot.margin= unit(c(0.2, 0, 0.2, 0.2), "lines"), panel.border = element_blank(), 
        title = element_text(size = 8))

ridge_f_treatments_dl <- reg_all_slopes_chosen_silva_tax |> # distribution plot at order taxonomic rank
  group_by(asv_num) |>
  dplyr::filter(n() >= 2) |>
  ungroup() |>
  slice_max(order_by = slope_chosen_days, n = 1, by = asv_num) |>
  dplyr::filter(treatment == 'DL') |>
  group_by(family_f) |>
  dplyr::filter(n() >= 4) |>
  group_by(family_f) |>
  mutate(counts = paste('n = ', n())) |>
  ungroup() |>
  ggplot(aes(y = fct_rev(family_f), x = slope_chosen_days, fill = phylum_f, label = counts))+
  geom_density_ridges(aes(fill = class_f, group = family_f), scale = 1, alpha = 0.6,
                      jittered_points = TRUE,
                      point_shape = 21, point_size = 0.2, point_alpha = 0.0,
                      quantile_lines = TRUE,
                      quantile_fun = mean)+
  scale_fill_manual(values = palette_class_assigned)+
  labs(y = '', fill= 'Phylum', x = expression("Growth rate (d"^"-1" *')'), title = 'Family')+
  scale_x_continuous(limits = c(0,11), expand = c(0,0))+
  geom_text(nudge_x = 7.8, nudge_y = 0.3, check_overlap = TRUE, size = 3)+
  theme_bw()+
  theme(legend.position = "none", strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0, size = 6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 6),
        plot.margin= unit(c(0.2, 0, 0.2, 0.2), "lines"), panel.border = element_blank(), 
        title = element_text(size = 8))

ridge_f_treatments_vl <- reg_all_slopes_chosen_silva_tax |> # distribution plot at order taxonomic rank
  group_by(asv_num) |>
  dplyr::filter(n() >= 2) |>
  ungroup() |>
  slice_max(order_by = slope_chosen_days, n = 1, by = asv_num) |>
  dplyr::filter(treatment == 'VL') |>
  group_by(family_f) |>
  dplyr::filter(n() >= 4) |>
  group_by(family_f) |>
  mutate(counts = paste('n = ', n())) |>
  ungroup() |>
  ggplot(aes(y = fct_rev(family_f), x = slope_chosen_days, fill = phylum_f, label = counts))+
  geom_density_ridges(aes(fill = class_f, group = family_f), scale = 1, alpha = 0.6,
                      jittered_points = TRUE,
                      point_shape = 21, point_size = 0.2, point_alpha = 0.0,
                      quantile_lines = TRUE,
                      quantile_fun = mean)+
  scale_fill_manual(values = palette_class_assigned)+
  labs(y = '', fill= 'Phylum', x = expression("Growth rate (d"^"-1" *')'), title = 'Family')+
  scale_x_continuous(limits = c(0,11), expand = c(0,0))+
  geom_text(nudge_x = 7.8, nudge_y = 0.3, check_overlap = TRUE, size = 3)+
  theme_bw()+
  theme(legend.position = "none", strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0, size = 6),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 6),
        plot.margin= unit(c(0.2, 0, 0.2, 0.2), "lines"), panel.border = element_blank(), 
        title = element_text(size = 8))

grid.arrange(ridge_f_treatments_cd, ridge_f_treatments_cl,
             ridge_f_treatments_pd, ridge_f_treatments_pl,
             ridge_f_treatments_dl, ridge_f_treatments_vl)




# ----- 7. GROWTH RATE TREATMENT RESPONSE ----------------
## labels----
effects_labels2 <-  as_labeller(c(effect_top_down_virus = 'Top-down effect viruses (VL-DL)',
                                  effect_top_down_grazers_dark = 'Top-down effect grazers dark (PD-CD)',
                                  effect_top_down_grazers_light = 'Top-down effect grazers light (PL-CL)',
                                  effect_bottom_up = 'Bottom-up effect (DL-PL)',
                                  effect_light_C = 'Light effect (CL-CD)',
                                  effect_light_P = 'Light effect (PL-PD)'))


## palette color and size by difference in growth rates----
palete_gradient <- c("#5e0000",
                     "#b24d5d",
                     '#FFFFFF' = 0,
                     "#4db2a2",
                     "#005a47") 

#240023

palete_gradient_cb <- c("#240023",
                     "#7d758a",
                     '#FFFFFF' = 0,
                     "#4db2a2",
                     "#005a47") 

palete_gradient_cb2 <- c("#2e0030",
                         "#797979",
                         '#FFFFFF' = 0,
                         "#79c295",
                         "#009759")

## Difference between growth rates at different treatments at ASV taxonomic rank (figure 5) -------
effects_difference_asv <- reg_all_slopes_chosen_silva_tax |> 
  distinct(treatment, season, asv_num, domain, phylum, class, order, genus, asv_num, family, .keep_all = TRUE) |>
  group_by(treatment, season, domain, phylum, class, order, family, asv_num) |>
  dplyr::summarise(slope_chosen_days_mean = mean(slope_chosen_days),
                   slope_chosen_mean = mean(slope_chosen),
                   slope_chosen_days_sd = sd(slope_chosen_days),
                   na.rm = TRUE) |>
  pivot_wider(names_from = treatment, values_from = c(slope_chosen_days_mean, slope_chosen_days_sd), 
              id_cols = c(season, domain, phylum, class, order, family, asv_num)) |>
  distinct() |>
  as_tibble() |>
  mutate(across(!c(domain, phylum, class, order, family, asv_num, season), as.numeric)) |>
  mutate(effect_top_down_virus = slope_chosen_days_mean_VL - slope_chosen_days_mean_DL,
         effect_top_down_grazers_dark = slope_chosen_days_mean_PD - slope_chosen_days_mean_CD,
         effect_top_down_grazers_light = slope_chosen_days_mean_PL - slope_chosen_days_mean_CL,
         effect_bottom_up = slope_chosen_days_mean_DL - slope_chosen_days_mean_PL,
         effect_light_C = slope_chosen_days_mean_CL - slope_chosen_days_mean_CD,
         effect_light_P = slope_chosen_days_mean_PL - slope_chosen_days_mean_PD)


##table with mean and sd 
# asv_bubble_plot <- effects_difference_asv_l |>
#   dplyr::filter(difference != is.na(difference)) |>
#   group_by(asv_num) |>
#   dplyr::filter(n() > 10) %$%
#   asv_num |>
#   unique()
# 
# effects_difference_asv_table <- effects_difference_asv |>
#   dplyr::filter(asv_num %in% asv_bubble_plot) |>
#   select(season, phylum, class, order, asv_num, family, matches('mean')) |> 
#   rename( 'Mean growth rate CD' = slope_chosen_days_mean_CD,
#           'Mean growth rate CL' = slope_chosen_days_mean_CL,
#           'Mean growth rate PL' = slope_chosen_days_mean_PL,
#           'Mean growth rate PD'= slope_chosen_days_mean_PD,
#           'Mean growth rate DL' = slope_chosen_days_mean_DL,
#           'Mean growth rate VL'= slope_chosen_days_mean_VL)
# 
# write.csv(effects_difference_asv_table, 'results/tables/mean_sd_gr_treatment_asv.csv')

#Heatmap plot ASV level
effects_difference_asv_l <- effects_difference_asv |>
  pivot_longer(cols = starts_with('effect'),
               names_to = 'effects',
               values_to = 'difference') |>
  as_tibble()

effects_difference_asv_l |>
  colnames()

effects_difference_asv_l$season <-effects_difference_asv_l$season |> 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))

effects_difference_asv_l$effects <- effects_difference_asv_l$effects |> 
  factor(levels=(c("effect_bottom_up", "effect_top_down_grazers_dark", 
                   "effect_top_down_grazers_light", "effect_top_down_virus",
                   "effect_light_C", "effect_light_P")))

effects_difference_asv_l <- effects_difference_asv_l |>
  mutate(phylum_f = as_factor(phylum),
         family_f = as_factor(family),
         order_f = as_factor(order),
         class_f = as_factor(class),
         asv_num_f = as_factor(asv_num))

effects_difference_asv_l$class_f <-  factor(effects_difference_asv_l$class_f, 
                                            levels=unique(effects_difference_asv_l$class_f[order(effects_difference_asv_l$phylum_f)]), 
                                            ordered=TRUE)

effects_difference_asv_l$order_f <-  factor(effects_difference_asv_l$order_f, 
                                            levels=unique(effects_difference_asv_l$order_f[order(effects_difference_asv_l$phylum_f,
                                                                                                 effects_difference_asv_l$class_f)]), 
                                            ordered=TRUE)

effects_difference_asv_l$family_f <-  factor(effects_difference_asv_l$family_f, 
                                             levels=unique(effects_difference_asv_l$family_f[order(effects_difference_asv_l$phylum_f,
                                                                                                   effects_difference_asv_l$class_f,
                                                                                                   effects_difference_asv_l$order_f)]), 
                                             ordered=TRUE)


effects_difference_asv_l$asv_num_f <-  factor(effects_difference_asv_l$asv_num_f, 
                                              levels=unique(effects_difference_asv_l$asv_num_f[order(effects_difference_asv_l$phylum_f,
                                                                                                     effects_difference_asv_l$class_f,
                                                                                                     effects_difference_asv_l$order_f,
                                                                                                     effects_difference_asv_l$family_f)]), 
                                              ordered=TRUE)

difference_gr_treatments_asv <- 
  effects_difference_asv_l |>
  dplyr::filter(difference != is.na(difference)) |>
  group_by(asv_num) |>
  dplyr::filter(n() > 10) |>
  group_by(asv_num, difference) |>
  mutate(counts = n(),
         family_asv_num = paste(family,'',asv_num)) |>
  ggplot(aes(season, effects))+
  geom_tile(aes(fill = difference), alpha = 1)+
  scale_y_discrete(labels = effects_labels2)+
  scale_size(range = c(0, 10), name ="Growth rate difference\n between treatments")+
  scale_fill_gradientn(colours = palete_gradient_cb)+
  #facet_grid(vars(family_asv_num))+
  facet_wrap(vars(family_asv_num))+
  labs(y = 'Treatments growth rates differences', x = 'Season', fill = '')+
  theme_bw()+
  theme(legend.position = 'right', panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 60, size = 5, hjust = 1), 
        axis.text.y = element_text(size = 5), legend.title = element_text(size = 5), axis.title = element_text(size = 7),
        strip.text = element_text(size = 5), strip.background = element_blank(), legend.text = element_text(size = 5))

# ggsave('difference_gr_treatments_asv_num_v5.pdf', difference_gr_treatments_asv,
#        path = "~/Documentos/Doctorat/REMEI/results/figures/corrected_asv_num/",
#        width = 180,
#        height = 160,
#        units = 'mm')

## Difference between growth rates at different treatments at family taxonomic rank (figure S6) ----
effects_difference_family <- reg_all_slopes_chosen_silva_tax |> 
  distinct(treatment, season, asv_num, domain, phylum, class, order, genus, asv_num, family, .keep_all = TRUE) |>
  group_by(treatment, season, domain, phylum, class, order, family) |>
  dplyr::summarise(slope_chosen_days_mean = mean(slope_chosen_days),
                   slope_chosen_mean = mean(slope_chosen),
                   slope_chosen_days_sd = sd(slope_chosen_days),
                   na.rm = TRUE) |>
  pivot_wider(names_from = treatment, values_from = c(slope_chosen_days_mean, slope_chosen_days_sd), 
              id_cols = c(season, domain, phylum, class, order, family)) |>
  distinct() |>
  as_tibble() |>
  mutate(across(!c(domain, phylum, class, order, family, season), as.numeric)) |>
  mutate(effect_top_down_virus = slope_chosen_days_mean_VL - slope_chosen_days_mean_DL,
         effect_top_down_grazers_dark = slope_chosen_days_mean_PD - slope_chosen_days_mean_CD,
         effect_top_down_grazers_light = slope_chosen_days_mean_PL - slope_chosen_days_mean_CL,
         effect_bottom_up = slope_chosen_days_mean_DL - slope_chosen_days_mean_PL,
         effect_light_C = slope_chosen_days_mean_CL - slope_chosen_days_mean_CD,
         effect_light_P = slope_chosen_days_mean_PL - slope_chosen_days_mean_PD)

# difference_gr_treatments_family <- effects_difference_family_l |>
#   filter(difference != is.na(difference)) |>
#   group_by(class) |>
#   filter(n() > 10) |>
#   group_by(family) |>
#   filter(n() > 20) |>
#   group_by(family, difference) |>
#   mutate(counts = n()) |>
#   ungroup() |>
#   filter(difference > 0) |>
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
#        path = "~/Documentos/Doctorat/REMEI/results/figures/corrected_asv_num/",
#        width = 240,
#        height = 180,
#        units = 'mm')

effects_difference_family_l <- effects_difference_family |>
  pivot_longer(cols = starts_with('effect'),
               names_to = 'effects',
               values_to = 'difference') |>
  as_tibble()

effects_difference_family_l$season <-effects_difference_family_l$season |> 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))

effects_difference_family_l$effects <- effects_difference_family_l$effects |> 
  factor(levels=(c("effect_bottom_up", "effect_top_down_grazers_dark", 
                   "effect_top_down_grazers_light", "effect_top_down_virus",
                   "effect_light_C", "effect_light_P")))

effects_difference_family_l <- effects_difference_family_l |>
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

difference_gr_treatments_family <- 
  effects_difference_family_l |>
  dplyr::filter(difference != is.na(difference)) |>
  group_by(family) |>
  dplyr::filter(n() > 10) |>
  group_by(family, difference) |>
  mutate(counts = n()) |>
  ggplot(aes(season, effects))+
  geom_tile(aes(fill = difference), alpha = 1)+
  scale_y_discrete(labels = effects_labels2)+
  scale_size(range = c(0, 10), name ="Growth rate difference\n between treatments")+
  scale_fill_gradientn(colours = palete_gradient_cb)+
  #facet_grid(vars(family_asv_num))+
  facet_wrap(vars(family_f))+
  labs(y = 'Treatments growth rates differences', x = 'Season', fill = '')+
  theme_bw()+
  theme(legend.position = 'right', panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 60, size = 5, hjust = 1), 
        axis.text.y = element_text(size = 5), legend.title = element_text(size = 5), axis.title = element_text(size = 7),
        strip.text = element_text(size = 5), strip.background = element_blank(), legend.text = element_text(size = 5))

ggsave('difference_gr_treatments_family6.pdf', difference_gr_treatments_family,
       path = "~/Documentos/Doctorat/REMEI/results/figures/corrected_asv_num/",
       width = 180,
       height = 160,
       units = 'mm')


## Difference between growth rates at different treatments at order taxonomic rank (figure S7) -----
effects_difference_order <- reg_all_slopes_chosen_silva_tax |> 
  distinct(treatment, season, asv_num, domain, phylum, class, order, genus, asv_num, family, .keep_all = TRUE) |>
  group_by(treatment, season, domain, phylum, class, order) |>
  dplyr::summarise(slope_chosen_days_mean = mean(slope_chosen_days),
                   slope_chosen_sd = sd(slope_chosen),
                   na.rm = TRUE) |>
  pivot_wider(names_from = treatment, values_from = c(slope_chosen_days_mean, slope_chosen_sd),
              id_cols = c(season, domain, phylum, class, order)) |>
  distinct() |>
  as_tibble() |>
  mutate(across(!c(domain, phylum, class, order, season), as.numeric)) |>
  mutate(effect_top_down_virus = slope_chosen_days_mean_VL - slope_chosen_days_mean_DL,
         effect_top_down_grazers_dark = slope_chosen_days_mean_PD - slope_chosen_days_mean_CD,
         effect_top_down_grazers_light = slope_chosen_days_mean_PL - slope_chosen_days_mean_CL,
         effect_bottom_up = slope_chosen_days_mean_DL - slope_chosen_days_mean_PL,
         effect_light_C = slope_chosen_days_mean_CL - slope_chosen_days_mean_CD,
         effect_light_P = slope_chosen_days_mean_PL - slope_chosen_days_mean_PD)

# effects_difference_order_table <- effects_difference_order |>
#   select(season, phylum, class, order, matches('mean'), matches('sd')) |>
#   filter(slope_chosen_days_mean_CD != is.na(slope_chosen_days_mean_CD) &
#            slope_chosen_days_mean_CL != is.na(slope_chosen_days_mean_CL) &
#            slope_chosen_days_mean_PD != is.na(slope_chosen_days_mean_PD) &
#            slope_chosen_days_mean_PL != is.na(slope_chosen_days_mean_PL) &
#            slope_chosen_days_mean_DL != is.na(slope_chosen_days_mean_DL) &
#            slope_chosen_days_mean_VL != is.na(slope_chosen_days_mean_VL)) |>
#   arrange(order)
# effects_difference_order_table |>
#   colnames()
# 
# 
# effects_difference_order_table <- effects_difference_order_table |>
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

effects_difference_order_l <- effects_difference_order |>
  pivot_longer(cols = starts_with('effect'),
               names_to = 'effects',
               values_to = 'difference') |>
  as_tibble()

effects_difference_order_l$season <-effects_difference_order_l$season |> 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))

effects_difference_order_l$effects <- effects_difference_order_l$effects |> 
  factor(levels=(c("effect_bottom_up", "effect_top_down_grazers_dark", 
                   "effect_top_down_grazers_light", "effect_top_down_virus",
                   "effect_light_C", "effect_light_P")))

effects_difference_order_l <- effects_difference_order_l |>
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
  effects_difference_order_l |>
  dplyr::filter(difference != is.na(difference)) |>
  group_by(order) |>
  dplyr::filter(n() > 12) |>
  group_by(order, difference) |>
  mutate(counts = n()) |>
  ungroup() |>
  ggplot(aes(season, effects))+ 
  geom_tile(aes(fill = difference), alpha = 1)+
  scale_y_discrete(labels = effects_labels2)+
  scale_size(range = c(0, 10), name ="Growth rate difference\n between treatments")+
  scale_fill_gradientn(colours = palete_gradient_cb)+
  #facet_grid(vars(family_asv_num))+
  facet_wrap(vars(order_f))+
  labs(y = 'Treatments growth rates differences', x = 'Season', fill = '')+
  theme_bw()+
  theme(legend.position = 'right', panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 60, size = 5, hjust = 1), 
        axis.text.y = element_text(size = 5), legend.title = element_text(size = 5), axis.title = element_text(size = 7),
        strip.text = element_text(size = 5), strip.background = element_blank(), legend.text = element_text(size = 5))

ggsave('difference_gr_treatments_order7.pdf', difference_gr_treatments_order,
       path = "~/Documentos/Doctorat/REMEI/results/figures/corrected_asv_num/",
       width = 180,
       height = 160,
       units = 'mm')


# ----- 8. ANALYSIS OF RESPONSIVE ASVs IN THE DIFFERENT TREATMENTS AND SEASONS ----------------
##palette ----
palete_gradient <- c(
  "#e6e6e6",
  "#003318") 


## Richness at the natural insitu community (Figure 6A) ----
miau_asv_tab <- read_rds('../MIAU_seqs/MIAU_runs_seqtab_final.rds')  |>
  as_tibble(rownames = 'sample_code')

tax_silva_miau <- readRDS("../MIAU_seqs/MIAUruns_assignTax_tax_assignation.rds") |> 
  as_tibble(rownames = 'sequence')

miau_asv_tab_insitu <- miau_asv_tab |>
  filter(sample_code %in% c("x254", "x255", "x256", "x257")) |>
  mutate(season = case_when(sample_code == 'x254' ~ 'Winter',
                            sample_code == 'x255' ~ 'Spring',
                            sample_code == 'x256' ~ 'Summer',
                            sample_code == 'x257' ~ 'Fall')) |>
  dplyr::select(-sample_code)

##reads per sample 
# miau_asv_tab_insitu_t <- miau_asv_tab_insitu |>
#   t() |>
#   as_tibble()
# 
# miau_asv_tab_insitu_t[-1,] |>
#   dplyr::mutate(across(c(V1, V2, V3, V4), as.numeric))|>
#   colSums()

miau_asv_tab_insitu_season <- miau_asv_tab_insitu$season

miau_asv_tab_insitu_filt <- miau_asv_tab_insitu |>
  dplyr::select(-season) |>
  dplyr::select(where(~ sum(.) != 0)) |>
  cbind(miau_asv_tab_insitu_season)

miau_seqs <- miau_asv_tab_insitu_filt |>
  t() |>
  as.data.frame() |>
  rownames_to_column(var = 'ASV_seq') |>
  row_to_names(row_number = 3874, remove_row = T, remove_rows_above = F)

miau_seqs |>
  colnames()

tax_silva_miau_filt <- tax_silva_miau |> 
  dplyr::filter(!is.na('Kingdom'), !is.na('Phylum')) |> #filter all results without domain assignation 
  dplyr::filter(Order !=  'Chloroplast') |>  #And the Chloros/Mitochondria seqs
  dplyr::filter(Family !=  'Mitochondria') 

#filter miau ASV tab by only the one's that are not unclassified or chloroplasts and mitochondria
miau_seqs_filt  <- miau_seqs |>
  filter(miau_seqs$miau_asv_tab_insitu_season %in% tax_silva_miau_filt$sequence) #2548 ASV

miau_seqs_filt_ed <- miau_seqs_filt |>
  mutate(across(!miau_asv_tab_insitu_season, as.numeric))

miau_seqs_filt_ed_fasta <- miau_seqs_filt_ed |>
  mutate(names = str_c( 'asv' , 1:nrow(miau_seqs_filt_ed))) |>
  select(-miau_asv_tab_insitu_season)

##calculate richness
miau_asv_tab_insitu_filt <- miau_seqs_filt_ed_fasta |>
  transmute(across(!'names', as.numeric)) 

miau_asv_tab_insitu_filt_ed <- otu_table(miau_asv_tab_insitu_filt, taxa_are_rows = T)

richness <- miau_asv_tab_insitu_filt_ed  |>
  estimate_richness() |>
  cbind(miau_asv_tab_insitu_season) |>
  rownames_to_column(var = 'season')

## Insitu richness plot----
richness$miau_asv_tab_insitu_season <- factor(richness$miau_asv_tab_insitu_season, levels = c('Winter', 'Spring', 'Summer',  'Fall'))
insitu_richness <- 
  richness |>
  ggplot(aes(miau_asv_tab_insitu_season, Observed))+
  geom_col(aes(fill = miau_asv_tab_insitu_season), alpha = 0.9)+#
  scale_fill_manual(values = palette_seasons_4)+
  #geom_text(aes(label = (treatment), y = (1000)))+
  coord_polar()+
  labs(x = '', y = '')+#, title = 'T0'
  annotate('text', x = 0, y = c(750 ,1000, 1300), label = c('750', '1000', '1300')) +
  #annotate('segment', x = 0, xend = 0, y = 0, yend =5.5)+
  scale_y_continuous(limits = c(0, 1350), breaks = c(750, 1000, 1300))+
  theme_bw()+
  theme(legend.position = 'none',
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_blank(),
        panel.grid.major.x = element_blank(),
        text = element_text(size = 5),
        plot.margin = margin(0,2,0,0, unit = 'cm' ))



## Total responsive ASVs GR>1(day-1) (figure 6B) ----
total_responsive_asvs <- reg_all_slopes_chosen_silva_tax |>
  filter(slope_chosen_days > 1) |>
  dplyr::select(treatment, season, asv_num) |> 
  mutate(seas_treat = paste(treatment, season)) |>
  group_by(seas_treat, asv_num) |> 
  distinct(seas_treat, asv_num) |>
  group_by(seas_treat) |>
  dplyr::summarize(num_responsive_asv = n()) |>
  as_tibble()

## Relative responsive ASVs (gr > 1) from total growing ASVs (gr > 0):
total_growing_asvs <- reg_all_slopes_chosen_silva_tax |>
  filter(slope_chosen_days > 0 &
           pvalue_slope_chosen < 0.05) |>
  dplyr::select(treatment, season, asv_num) |> #, slope_chosen_days, family_f
  #filter(treatment %in% c('CD', 'CL', 'PD', 'PL')) |>
  #dplyr::select(-slope_chosen_days, -family_f) |>
  mutate(seas_treat = paste(treatment, season)) |>
  group_by(seas_treat, asv_num) |> ##sembla que PD Winter está repetit al dataset original (trobar a on es duplica, té efecte aals gràfics?)
  distinct(seas_treat, asv_num) |>
  group_by(seas_treat) |>
  dplyr::summarize(num_growing_asv = n()) |>
  #mutate(num_responding_asv = as.numeric(num_responding_asv)) |>
  as_tibble()

total_responsive_asvs |>
  left_join(total_growing_asvs) |>
  mutate(relative_responsive_asvs = num_responsive_asv/num_growing_asv)

## plot relative responsive ASVs per treatment-season ----
total_responsive_asvs_plot_relative <- 
  total_responsive_asvs |>
  left_join(total_growing_asvs) |>
  mutate(relative_responsive_asvs = num_responsive_asv/num_growing_asv) |>
  separate('seas_treat', sep = ' ', into = c('Treatment', 'Season'), remove = F) |>
  as_tibble() |>
  mutate(seas_treat = fct_relevel(seas_treat, c("CD Winter", "CL Winter",  "PD Winter" , "PL Winter", "DL Winter", "VL Winter",
                                                "CD Spring", "CL Spring" , "PD Spring", "PL Spring", "DL Spring", "VL Spring", 
                                                "CD Summer", "CL Summer", "PD Summer", "PL Summer" , "DL Summer"  , "VL Summer",
                                                "CD Fall", "CL Fall", "PD Fall", "PL Fall", "DL Fall", "VL Fall"))) |>
  ggplot(aes(seas_treat, relative_responsive_asvs, fill = Season))+
  geom_col(alpha = 0.9)+
  geom_text(aes(label = (Treatment), y = (1.1)))+
  coord_polar()+
  labs(x = '', y = '')+
  scale_fill_manual(values = palette_seasons_4)+
  annotate('text', x = 0, y = c(0.2, 0.4, 0.6, 0.8, 1), label = c('0.2', '0.4', '0.6', '0.8', '1')) +
  scale_y_continuous(limits = c(0, 1.1))+
  theme_bw()+
  theme(legend.position = 'none',
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_blank(),
        panel.grid.major.x = element_blank(),
        text = element_text(size = 5),
        plot.margin = margin(0,2,0,0, unit = 'cm'))

# ggsave(filename = 'total_responsive_asvs_plot_relative.pdf', plot = total_responsive_asvs_plot_relative, 
#        path = 'results/figures/corrected_asv_num/',
#        width = 88, height = 88, units = 'mm')

total_responsive_asvs <-  total_responsive_asvs |>
  separate('seas_treat', sep = ' ', into = c('Treatment', 'Season'), remove = F) |>
  as_tibble() |>
  mutate(seas_treat = fct_relevel(seas_treat, c("CD Winter", "CL Winter",  "PD Winter" , "PL Winter", "DL Winter", "VL Winter",
                                                "CD Spring", "CL Spring" , "PD Spring", "PL Spring", "DL Spring", "VL Spring", 
                                                "CD Summer", "CL Summer", "PD Summer", "PL Summer" , "DL Summer"  , "VL Summer",
                                                "CD Fall", "CL Fall", "PD Fall", "PL Fall", "DL Fall", "VL Fall")))
total_responsive_asvs %$%
  num_responsive_asv |>
  range()

## plot absolute number of responsive ASVs ----
total_responsive_asvs_plot_absolut <- 
  total_responsive_asvs |>
  ggplot(aes(seas_treat, num_responsive_asv, fill = Season))+
  geom_col()+
  geom_text(aes(label = (Treatment), y = (300)))+
  coord_polar()+
  labs(x = '', y = '')+
  scale_fill_manual(values = palette_seasons_4)+
  annotate('text', x = 0, y = c(100, 200, 300), label = c('100', '200', '300')) +
  scale_y_continuous(limits = c(0, 300))+
  theme_bw()+
  theme(legend.position = 'none',
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        text = element_text(size = 5),
        plot.margin = margin(0,0,0,0))

# ggsave(filename = 'total_responsive_asvs_plot_absolut.pdf', plot = total_responsive_asvs_plot_absolut, 
#        path = 'results/figures/corrected_asv_num/',
#        width = 88, height = 88, units = 'mm')


## Number of exclusive ASVs growing GR>1(day-1) per sample----------------
data_for_exclusive_asvs <- reg_all_slopes_chosen_silva_tax |>
  filter(slope_chosen_days > 1) |>
  dplyr::select(treatment, season, asv_num) |> 
  mutate(seas_treat = paste(treatment, season)) |>
  group_by(seas_treat, asv_num) |> 
  distinct(seas_treat, asv_num) 

total_responding_asvs_sample  <-  data_for_exclusive_asvs |>
  group_by(seas_treat) |>
  dplyr::summarize(num_responding_asv = n()) |>
  mutate(num_responding_asv = as.numeric(num_responding_asv)) |>
  as_tibble()

exclusive_asvs_dataset <- list()
for(i in unique(data_for_exclusive_asvs$seas_treat)) {
  exclusive_asvs_dataset[[i]] <- exclusive.asvs(data_for_exclusive_asvs, sample = i)
}

exclusive_asvs_dataset <-  bind_rows(exclusive_asvs_dataset) |>
  group_by(seas_treat) |>
  dplyr::summarize(number_exclusive_asvs = n()) |>
  separate(seas_treat, ' ', into = c('Treatment', 'Season'), remove = F) |>
  mutate(seas_treat = fct_relevel(seas_treat, c("CD Winter", "CL Winter",  "PD Winter" , "PL Winter", "DL Winter", "VL Winter",
                                                "CD Spring", "CL Spring" , "PD Spring", "PL Spring", "DL Spring", "VL Spring", 
                                                "CD Summer", "CL Summer", "PD Summer", "PL Summer" , "DL Summer"  , "VL Summer",
                                                "CD Fall", "CL Fall", "PD Fall", "PL Fall", "DL Fall", "VL Fall")))
exclusive_asvs_dataset %$%
  number_exclusive_asvs |>
  range()

## Relative exclusive ASVs plot (Fig. 6C) ----
exclusive_asvs_relative <-  exclusive_asvs_dataset |>
  left_join(total_responding_asvs_sample, by = 'seas_treat') |>
  mutate(number_exclusive_asvs = as.numeric(number_exclusive_asvs)) |>
  as_tibble() |>
  mutate(relative_exclusive_asvs = number_exclusive_asvs/num_responding_asv) |>
  mutate(seas_treat = fct_relevel(seas_treat, c("CD Winter", "CL Winter",  "PD Winter" , "PL Winter", "DL Winter", "VL Winter",
                                                "CD Spring", "CL Spring" , "PD Spring", "PL Spring", "DL Spring", "VL Spring", 
                                                "CD Summer", "CL Summer", "PD Summer", "PL Summer" , "DL Summer"  , "VL Summer",
                                                "CD Fall", "CL Fall", "PD Fall", "PL Fall", "DL Fall", "VL Fall")))

exclusive_asvs_relative_plot <- exclusive_asvs_relative |>
  ggplot(aes(seas_treat, relative_exclusive_asvs, fill = Season))+
  geom_col()+
  geom_text(aes(label = (Treatment), y = (1)))+
  coord_polar()+
  labs(x = '', y = '')+
  scale_fill_manual(values = palette_seasons_4)+
 # annotate('text', x = 0, y = c(0.2, 0.4), label = c('0.2', '0.4')) +
  annotate('text', x = 0, y = c(0.2, 0.4, 0.6, 0.8, 1), label = c('0.2', '0.4', '0.6', '0.8', '1')) +
  scale_y_continuous(limits = c(0, 1.0))+
  #scale_y_continuous(limits = c(0, 0.55))+
  theme_bw()+
  theme(legend.position = 'none',
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_blank(),
        panel.grid.major.x = element_blank(),
        text = element_text(size = 5),
        plot.margin = margin(0,2,0,0, unit = 'cm'))

# ggsave(filename = 'exclusive_asvs_relative_1gr_ed2.pdf', plot = exclusive_asvs_relative_plot,
#        path = 'results/figures/corrected_asv_num/',
#        width = 188, height = 188, units = 'mm')

## Relative exclusive ASVs plot (Fig S8C) ----
exclusive_asvs_absolut <- exclusive_asvs_dataset |>
  ggplot(aes(seas_treat, number_exclusive_asvs, fill = Season))+
  geom_col()+
  geom_text(aes(label = (Treatment), y = (110)))+
  coord_polar()+
  labs(x = '', y = '')+
  scale_fill_manual(values = palette_seasons_4)+
  annotate('text', x = 0, y = c(20, 40, 60, 80, 100), label = c('20', '40', '60', '80', '100')) +
  scale_y_continuous(limits = c(0, 110))+
  theme_bw()+
  theme(legend.position = 'none',
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_blank(),
        panel.grid.major.x = element_blank())

# ggsave(filename = 'exclusive_asvs_absolut_1gr.pdf', plot = exclusive_asvs_absolut, 
#        path = 'results/figures/corrected_asv_num/',
#        width = 88, height = 88, units = 'mm')



####Total responsive ASVs (GR > 1)

total_responsive_asvs <- reg_all_slopes_chosen_silva_tax |>
  filter(slope_chosen_days > 1 &
           pvalue_slope_chosen < 0.05) |>
  dplyr::select(treatment, season, asv_num) |> #, slope_chosen_days, family_f
  #filter(treatment %in% c('CD', 'CL', 'PD', 'PL')) |>
  #dplyr::select(-slope_chosen_days, -family_f) |>
  mutate(seas_treat = paste(treatment, season)) |>
  group_by(seas_treat, asv_num) |> ##sembla que PD Winter está repetit al dataset original (trobar a on es duplica, té efecte aals gràfics?)
  distinct(seas_treat, asv_num) |>
  group_by(seas_treat) |>
  dplyr::summarize(num_responsive_asv = n()) |>
  #mutate(num_responding_asv = as.numeric(num_responding_asv)) |>
  as_tibble()

total_responsive_asvs <-  total_responsive_asvs |>
  #left_join(, by = 'seas_treat') |>
  separate('seas_treat', sep = ' ', into = c('Treatment', 'Season'), remove = F) |>
  #mutate(number_exclusive_asvs = as.numeric(number_exclusive_asvs)) |>
  as_tibble() |>
  #mutate(relative_exclusive_asvs = number_exclusive_asvs/num_responding_asv) |>
  mutate(seas_treat = fct_relevel(seas_treat, c("CD Winter", "CL Winter",  "PD Winter" , "PL Winter", "DL Winter", "VL Winter",
                                                "CD Spring", "CL Spring" , "PD Spring", "PL Spring", "DL Spring", "VL Spring", 
                                                "CD Summer", "CL Summer", "PD Summer", "PL Summer" , "DL Summer"  , "VL Summer",
                                                "CD Fall", "CL Fall", "PD Fall", "PL Fall", "DL Fall", "VL Fall")))
total_responsive_asvs %$%
  num_responsive_asv |>
  range()

total_responsive_asvs_plot_absolut <- 
  total_responsive_asvs |>
  ggplot(aes(seas_treat, num_responsive_asv, fill = Season))+
  geom_col()+
  geom_text(aes(label = (Treatment), y = (300)))+
  coord_polar()+
  labs(x = '', y = '')+
  scale_fill_manual(values = palette_seasons_4)+
  annotate('text', x = 0, y = c(100, 200, 300), label = c('100', '200', '300')) +
  scale_y_continuous(limits = c(0, 300))+
  theme_bw()+
  theme(legend.position = 'none',
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        text = element_text(size = 5),
        plot.margin = margin(0,0,0,0))

# ggsave(filename = 'total_responsive_asvs_plot_absolut.pdf', plot = total_responsive_asvs_plot_absolut, 
#        path = 'results/figures/corrected_asv_num/',
#        width = 88, height = 88, units = 'mm')


## Responsive ASVs: GR>1(day-1) -------------
data_for_common_asvs <- reg_all_slopes_chosen_silva_tax |>
  dplyr:: filter(slope_chosen_days > 1) |>
  dplyr::select(treatment, season, asv_num) |> 
  mutate(seas_treat = paste(treatment, season)) |>
  group_by(seas_treat, asv_num) |> 
  distinct(seas_treat, asv_num) 

num_asv_seas_treat <- data_for_common_asvs |>
  group_by(seas_treat) |>
  dplyr::summarize(n = n()) |>
  as_tibble()
## Calculating nº of common ASVs between treatments and seasons------
##between treatments from the same season
common_asvs_dataset <- list()
for(i in unique(data_for_common_asvs$seas_treat)) {
  common_asvs_dataset[[i]] <- common.asvs(data_for_common_asvs, sample1 = i, sample2 = i)
}

common_asv_sum <- bind_rows(common_asvs_dataset) |>
  separate(seas_treat.x, '_', into = c('from', 'to')) |>
  left_join(num_asv_seas_treat, by = c('from' = 'seas_treat')) |>
  left_join(num_asv_seas_treat, by = c('to' = 'seas_treat')) |>
  mutate(division = case_when( n.x < n.y ~ n.x,
                               n.y < n.x ~ n.y,
                               n.x == n.y ~ n.x)) |>
  mutate(relative_connectivity = num_common_asvs/division) |>
  mutate(sample = paste0(from, to))  

##prepare graph
edges <- common_asv_sum |>
  select(from, to) |>
  filter(from != is.na(from))
vertices <- common_asv_sum |>
  select(sample) |>
  mutate(size = 10) |>
  as_tibble()

vertices <- data.frame(rbind(as_tibble(common_asv_sum$from), 
                             as_tibble(common_asv_sum$to))) |>
  group_by(value) |>
  unique() |>
  mutate(size = 10) |>
  separate(value, ' ', into = c('Treatment', 'Season'), remove = F) |>
  mutate(Treatment = as.factor(Treatment),
         Season = as.factor(Season),
         value = as.factor(value)) |>
  select(value, Season, Treatment) |>
  filter(value != is.na(value)) |>
  as_tibble()

vertices$Season <- vertices$Season |> 
  factor(levels = c('Winter', 'Spring', 'Summer', 'Fall'))

vertices$Treatment <- vertices$Treatment |> 
  factor(levels = c('CD', 'CL', 'PD', 'PL', 'DL', 'VL'))

vertices <- vertices |>
  mutate(value = fct_relevel(value, c("CD Winter", "CL Winter",  "PD Winter" , "PL Winter", "DL Winter", "VL Winter",
                                      "CD Spring", "CL Spring" , "PD Spring", "PL Spring", "DL Spring", "VL Spring", 
                                      "CD Summer", "CL Summer", "PD Summer", "PL Summer" , "DL Summer"  , "VL Summer",
                                      "CD Fall", "CL Fall", "PD Fall", "PL Fall", "DL Fall", "VL Fall"))) |>
  
  arrange(value) |>
  as_tibble()

connections <- common_asv_sum |>
  separate(from, ' ', into = c('Treatment', 'Season'), remove = F) |>
  filter(relative_connectivity != is.na(relative_connectivity)) |>
  select(from, to, relative_connectivity, Treatment, Season) 

num_connections1 <- common_asv_sum |>
  group_by(from) |>
  dplyr::summarize(n = n()) |>
  right_join(vertices, by = c('from' = 'value'))

num_connections2 <- common_asv_sum |>
  group_by(to) |>
  dplyr::summarize(n = n()) |>
  right_join(vertices, by = c('to' = 'value'))

num_connections12 <- num_connections1 |>
  left_join(num_connections2, by=c('Season' = 'Season', 'Treatment' = 'Treatment')) |>
  mutate(n.x = case_when(#!is.na(n.x) ~ n.x,
    is.na(n.x) ~ '0',
    TRUE ~ as.character(n.x)),
    n.y = case_when(#!is.na(n.x) ~ n.x,
      is.na(n.y) ~ '0',
      TRUE ~ as.character(n.y))
  ) |>
  mutate(num_connections = (as.numeric(n.x)+as.numeric(n.y)))

common_asv_sum$from
common_asv_sum$to

# Build a network object from this dataset:
mygraph_rel <- graph_from_data_frame(edges, vertices = vertices)

# The connection object must refer to the ids of the leaves:
from = match( connections$from, vertices$value)
to = match( connections$to, vertices$value)

##CIRCULAR CONNECTIONS 
specific_range <- common_asv_sum %$%
  relative_connectivity |>
  range()

##Plot connections > 30%  (relative connections, figure 6D) ----
common_asvs_1gr_relative <- ggraph(mygraph_rel, layout="linear", circular = TRUE) + 
  geom_edge_arc(aes(width = as.numeric(connections$relative_connectivity), 
                    color = as.numeric(connections$relative_connectivity)),
  fold=FALSE,
  edge_alpha=0.8,
  lineend = 'square',
  linejoin = 'mitre', check_overlap = TRUE) +
  scale_edge_colour_gradientn(colours = palete_gradient, name = 'Relative\nconnectivity')+
  scale_edge_width_continuous(limits = c(0.3, 0.89), guide = 'none')+
  coord_fixed()+
  geom_node_point(aes(color = Season, size = num_connections12$num_connections/2), alpha = 0.8)+
  scale_size(range=c(11/2,19/2), limits = c(11/2,19/2)) +
  geom_node_text(aes(label=Treatment), vjust = -0.1, show.legend = TRUE) + 
  scale_color_manual(values = palette_seasons_4)+
  guides(color = guide_legend(ncol = 1, size = 10,
                              override.aes = aes(label = '')),
         shape = guide_legend(ncol = 1, size = 1),
         size = 'none')+
  theme_void() +
  theme(legend.text = element_text(size = 5), 
        legend.title = element_text(size =7),
        legend.box.spacing = unit(1, 'mm'))

# ggsave(filename = 'common_asvs_all_1gr_relative_02.pdf', plot = common_asvs_1gr_relative, 
#        path = 'results/figures/corrected_asv_num/',
#        width = 188, height = 188, units = 'mm')

## Plot connections > 30%  (absolut numbers, figure S8D) ----
edges <- common_asv_sum |>
  select(from, to) |>
  filter(from != is.na(from))
vertices <- common_asv_sum |>
  select(sample) |>
  mutate(size = 10) |>
  as_tibble()

vertices <- data.frame(rbind(as_tibble(common_asv_sum$from), 
                             as_tibble(common_asv_sum$to))) |>
  group_by(value) |>
  unique() |>
  mutate(size = 10) |>
  separate(value, ' ', into = c('Treatment', 'Season'), remove = F) |>
  mutate(Treatment = as.factor(Treatment),
         Season = as.factor(Season),
         value = as.factor(value)) |>
  select(value, Season, Treatment) |>
  filter(value != is.na(value)) |>
  as_tibble()

vertices$Season <- vertices$Season |> 
  factor(levels = c('Winter', 'Spring', 'Summer', 'Fall'))

vertices$Treatment <- vertices$Treatment |> 
  factor(levels = c('CD', 'CL', 'PD', 'PL', 'DL', 'VL'))

vertices <- vertices |>
  mutate(value = fct_relevel(value, c("CD Winter", "CL Winter",  "PD Winter" , "PL Winter", "DL Winter", "VL Winter",
                                      "CD Spring", "CL Spring" , "PD Spring", "PL Spring", "DL Spring", "VL Spring", 
                                      "CD Summer", "CL Summer", "PD Summer", "PL Summer" , "DL Summer"  , "VL Summer",
                                      "CD Fall", "CL Fall", "PD Fall", "PL Fall", "DL Fall", "VL Fall"))) |>
  
  arrange(value) |>
  as_tibble()

connections <- common_asv_sum |>
  separate(from, ' ', into = c('Treatment', 'Season'), remove = F) |>
  filter(num_common_asvs != is.na(num_common_asvs)) |>
  select(from, to, num_common_asvs, Treatment, Season) 

num_connections1 <- common_asv_sum |>
  group_by(from) |>
  dplyr::summarize(n = n()) |>
  right_join(vertices, by = c('from' = 'value'))

num_connections2 <- common_asv_sum |>
  group_by(to) |>
  dplyr::summarize(n = n()) |>
  right_join(vertices, by = c('to' = 'value'))

num_connections12 <- num_connections1 |>
  left_join(num_connections2, by=c('Season' = 'Season', 'Treatment' = 'Treatment')) |>
  mutate(n.x = case_when(#!is.na(n.x) ~ n.x,
    is.na(n.x) ~ '0',
    TRUE ~ as.character(n.x)),
    n.y = case_when(#!is.na(n.x) ~ n.x,
      is.na(n.y) ~ '0',
      TRUE ~ as.character(n.y))
  ) |>
  mutate(num_connections = (as.numeric(n.x)+as.numeric(n.y)))

# Build a network object from this dataset:
mygraph <- graph_from_data_frame(edges, vertices = vertices)

# The connection object must refer to the ids of the leaves:
from = match( connections$from, vertices$value)
to = match( connections$to, vertices$value)

##CIRCULAR CONNECTIONS 
specific_range <- common_asv_sum %$%
  num_common_asvs |>
  range() ##the 30% percent of 132 is 39,6 plot > 39.6 common ASVs

common_asvs_1gr_absolut <- ggraph(mygraph, layout="linear", circular = TRUE) +
  geom_edge_arc(aes(width = as.numeric(connections$num_common_asvs), 
                    color = as.numeric(connections$num_common_asvs)
  ),
  fold=FALSE,
  edge_alpha=0.8,
  lineend = 'square',
  linejoin = 'mitre', check_overlap = TRUE) + 
  scale_edge_colour_gradientn(colours = palete_gradient, name = 'Nº of\ncommon\nASVs')+
  scale_edge_width_continuous(limits = c(39.6, 132), guide = 'none')+
  coord_fixed()+
  geom_node_point(aes(color = Season, size = num_connections12$num_connections/2), alpha = 0.8)+#, size = 0.75
  scale_size(range=c(11/2,19/2), limits = c(11/2,19/2)) +
  geom_node_text(aes(label=Treatment), vjust = 0.1, show.legend = TRUE) + #, angle=65, hjust=1, nudge_y = -1.1, size=2.3
  scale_color_manual(values = palette_seasons_4)+
  guides(color = guide_legend(ncol = 1, size = 10,
                              override.aes = aes(label = '')),
         shape = guide_legend(ncol = 1, size = 1),
         size = 'none')+
  theme_void() +
  theme(legend.text = element_text(size = 5), 
        legend.title = element_text(size =7),
        legend.box.spacing = unit(1, 'mm'))
# 
# ggsave(filename = 'common_asvs_all_1gr_absolut_30perc.pdf', plot = common_asvs_1gr_absolut, 
#        path = 'results/figures/corrected_asv_num/',
#        width = 188, height = 188, units = 'mm')






# ----- 9. RELATIONSHIP BETWEEN ABUNDANCE AND GROWTH ----------------
# PREVIOUSLY WE PERFORMED A CLUSTERING AT 100% BETWEEN ASV SEQUENCES FROM IN SITU DATA & GROWTH RATES FROM THE EXPERIMENTS ---------
##for these purpuse we used vsearch function cluster_fast, which allows us to clusterize the fasta sequences in filename 

###ONLY ASV PRESENT IN SITU AND GROWING ASV (SIGNIFICATIVELY)
#Load results data from the clustering at 100% 
cluster_miau_remei <- read.table('data/cluster_remei_miau_asv/clustering_results.txt', header = F)
cluster_miau_remei |>
  dim() #5196 asv relacionades entre un dataset i l'altre 

colnames(cluster_miau_remei) <- c('v1', 'v2', 'v3', 'clustering', 'v5', 'v6', 'v7', 'v8', 'col1_asv', 'col2_asv')
str_count(cluster_miau_remei$col1_asv, "remei") |>
  sum() #3642
str_count(cluster_miau_remei$col1_asv, "miau") |>
  sum() #1554
str_count(cluster_miau_remei$col2_asv, "remei") |>
  sum() #1554
str_count(cluster_miau_remei$col2_asv, "miau") |>
  sum() #3642

cluster_miau_remei <- cluster_miau_remei |>
  separate(col1_asv, c('asv_num_col1', 'size_col1'), sep = ';') |>
  separate(col2_asv, c('asv_num_col2', 'size_col2'), sep = ';')

cluster_miau_remei <- cluster_miau_remei |>
  mutate(relation_miau_remei1 = paste(asv_num_col1,'-',asv_num_col2)
  )

cluster_miau_remei$relation_miau_remei1 |>
  unique() #5196 
cluster_miau_remei  |>
  dim()

## Add ASV sequences to relate both datasets
miau_seqs <- read.fasta('data/cluster_remei_miau_asv/MIAU_runs_seqtab_final.fasta')
remei_seqs <- read.fasta('data/cluster_remei_miau_asv/remei_1_2_pool_seqtab_final.fasta')

dna_REMEI<- readDNAStringSet('data/cluster_remei_miau_asv/remei_1_2_pool_seqtab_final.fasta')
seq_name = names(dna_REMEI)
sequence = paste(dna_REMEI)
remei_seqs <- data.frame(seq_name, sequence) |>
  as_tibble() |>
  mutate(asv_num_ed = str_replace(seq_name, ';', 'remei;')) |>
  separate(asv_num_ed, c('asv_num', 'size'), sep = ';')

dna_miau <- readDNAStringSet('data/cluster_remei_miau_asv/MIAU_runs_seqtab_final.fasta')
seq_name = names(dna_miau)
sequence = paste(dna_miau)

miau_seqs <- data.frame(seq_name, sequence) |>
  as_tibble() |>
  mutate(asv_num_ed = str_replace(seq_name, ';', 'miau;')) |>
  separate(asv_num_ed, c('asv_num', 'size'), sep = ';')

cluster_miau_remei <- cluster_miau_remei |>
  mutate(miau_asv_num = as.character(str_extract_all(relation_miau_remei1, pattern = ('asv[0-9]*(miau)+'))),
         remei_asv_num = as.character(str_extract_all(relation_miau_remei1, pattern = ('asv[0-9]*(remei)+')))) 

cluster_remei_seqs <- cluster_miau_remei |>
  left_join(remei_seqs, by = c('remei_asv_num' = 'asv_num'))

##unload plyr and dplyr if not working
cluster_remei_seqs <- cluster_remei_seqs %>%
  rename(seq_name = 'seq_name_remei',
         sequence = 'sequence_remei',
         size = 'size_remei')

cluster_miau_seqs <- cluster_miau_remei |>
  left_join(miau_seqs, by = c('miau_asv_num' = 'asv_num'))

cluster_miau_seqs <- cluster_miau_seqs |>
  rename(seq_name = 'seq_name_miau',
         sequence = 'sequence_miau',
         size = 'size_remei')

cluster_miau_remei %$%
  asv_num_col1 |>
  unique() |> #5196 
  length()

cluster_miau_remei %$%
  asv_num_col2 |>
  unique()|> #5196
  length()

#prepare growth rates experiments dataset
##taxonomy from experiments dataset
remei_tax <- rem_fc_filt@tax_table |>
  as_tibble()

##join ASVs sequences and growth rates
growth_rates_remei_seqs <- reg_all_slopes_chosen_silva_tax |>
  inner_join(remei_tax,  by = 'asv_num', copy = TRUE)|>
  as_tibble()

##Prepare insitu dataset
###Import data
# miau_seqs <- read.fasta('../MIAU_seqs/MIAU_runs_seqtab_final.fasta')

miau_asv_tab <- read_rds('../MIAU_seqs/MIAU_runs_seqtab_final.rds')  |>
  as_tibble(rownames = 'sample_code')

tax_silva_miau <- readRDS("../MIAU_seqs/MIAUruns_assignTax_tax_assignation.rds") |> 
  as_tibble(rownames = 'sequence')

miau_asv_tab_insitu <- miau_asv_tab |>
  filter(sample_code %in% c("x254", "x255", "x256", "x257"))

miau_asv_tab_insitu_codes <- miau_asv_tab_insitu$sample_code

miau_asv_tab_insitu_filt <- miau_asv_tab_insitu |>
  dplyr::select(-sample_code) |>
  dplyr::select(where(~ sum(.) != 0)) |>
  cbind(miau_asv_tab_insitu_codes)

miau_seqs <- miau_asv_tab_insitu_filt |>
  t() |>
  as.data.frame() |>
  rownames_to_column(var = 'ASV_seq') |>
  row_to_names(row_number = 3874, remove_row = T, remove_rows_above = F)

tax_silva_miau_filt <- tax_silva_miau |>  # filter all results without domain assignation 
  dplyr::filter(!is.na('Kingdom'), !is.na('Phylum')) |> #And the Chloros/Mitochondria seqs
  dplyr::filter(Order !=  'Chloroplast') |> 
  dplyr::filter(Family !=  'Mitochondria') 

miau_seqs_filt  <- miau_seqs |> #filter miau seqs by only the one's that are not unclassified or chloroplasts and mitochondria
  dplyr::filter(miau_asv_tab_insitu_codes %in% tax_silva_miau_filt$sequence) #2548 ASV

miau_seqs_filt_ed <- miau_seqs_filt |>
  mutate(across(!miau_asv_tab_insitu_codes, as.numeric))

miau_seqs_filt_ed_fasta <- miau_seqs_filt_ed |>
  mutate(names = str_c( 'asv' , 1:nrow(miau_seqs_filt_ed)))

miau_seqs_filt_rel_abund <- apply(miau_seqs_filt_ed[,c(2:5)], 2,   function(x){x / sum(x)}) |>
  cbind(miau_seqs_filt_ed$miau_asv_tab_insitu_codes) |>
  as_tibble()

colnames(miau_seqs_filt_rel_abund) <- c("x254", "x255", "x256", "x257", "asv_seqs")

miau_seqs_filt_rel_abund_long <- miau_seqs_filt_rel_abund |>
  pivot_longer(cols = starts_with('x'), names_to = 'sample', values_to = 'rel_abund') |>
  mutate(season = case_when(sample == 'x254' ~ 'Winter',
                            sample == 'x255' ~ 'Spring',
                            sample == 'x256' ~ 'Summer',
                            sample == 'x257' ~ 'Fall'),
         rel_abund = as.numeric(rel_abund))

miau_seqs_filt_rel_abund_long_tax <- miau_seqs_filt_rel_abund_long |>
  left_join(tax_silva_miau_filt, by = c('asv_seqs' = 'sequence')) ##add taxonomy to check if its the same for both datasets

## Join GR data with in situ rel abund and sequences 
miau_seqs_filt_rel_abund_long <- miau_seqs_filt_rel_abund_long |> ### Add a 0 to the growing ASV that don't have a match in situ
  mutate(miau_seqs_filt_rel_abund_long = ifelse(is.na(rel_abund), 0, rel_abund))

## Join GR data with in situ relative abundance and sequences 
sequences_rel_abund <- cluster_miau_seqs |>
  dplyr::filter(relation_miau_remei1 != is.na(relation_miau_remei1)) |>
  right_join(miau_seqs_filt_rel_abund_long, by = c('sequence_miau' = 'asv_seqs'))

sequences_rel_abund %$%
  relation_miau_remei1 |>
  unique() |>
  length() #1515 asv 

## Join gr data form the clustering and ASV sequences
sequences_gr <- cluster_remei_seqs |>
  dplyr::filter(relation_miau_remei1 != is.na(relation_miau_remei1)) |>
  right_join(growth_rates_remei_seqs, by = c('sequence_remei' = '.otu')) |>
  filter(slope_chosen_days != is.na(slope_chosen_days) &
           asv_num_col1 != is.na(asv_num_col1))

sequences_gr %$%
  asv_num_col1 |>
  unique() |>
  length() #833 asv

# Joing gr data with relative abundance 
seqs_gr_rel_abund <- sequences_gr |>
  inner_join(sequences_rel_abund, by = c('relation_miau_remei1', 'season')) |>
  dplyr::filter(slope_chosen_days != is.na(slope_chosen_days))

seqs_gr_rel_abund <- seqs_gr_rel_abund |>  
  mutate(rel_abund = ifelse(is.na(rel_abund), 0, rel_abund))

seqs_gr_rel_abund <- seqs_gr_rel_abund |>
  ungroup() |>
  mutate(rel_abund_ed = as.numeric(rel_abund))  |>
  mutate(rel_abund_ed2 = round(rel_abund, 3)) |>
  arrange(rel_abund_ed) |>
  mutate(rank_abund_2 = rank(-rel_abund_ed),
    behaviour = case_when(rel_abund_ed < 0.001 ~ 'Rare x < 0.1% ',
                          rel_abund_ed > 0.01 ~ 'Abundant x > 1%',
                          rel_abund_ed <0.01  & rel_abund >0.001 ~ 'Mid  0.1% < x < 1%')) 

seqs_gr_rel_abund$treatment <- seqs_gr_rel_abund$treatment |> 
  factor(levels=(c("CL", "CD", "PL", "PD", "DL", "VL", 'NA')))
seqs_gr_rel_abund$season <- seqs_gr_rel_abund$season |> 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))

seqs_gr_rel_abund |>
  dplyr::filter(slope_chosen_days != (is.na(slope_chosen_days))) %$% 
  relation_miau_remei1 |>
  unique() |>
  length()

## table with % of rare reactive ASVs (add to plots as geom_text) 
##from the total potential bloomers GR >2  & <0.1% community

total_asv_match_remei_miau <- seqs_gr_rel_abund |>
  dplyr::filter(rel_abund_ed < 0.01 & slope_chosen_days > 2) |>
  group_by(season, treatment) |>
  dplyr::summarize(total_common_asv_potential_bloomers = n())

total_asv_match_remei_miau_potential_bloomers <- seqs_gr_rel_abund |>
  group_by(season, treatment) |>
  dplyr::summarize(total_common_asv = n()) |>
  dplyr::select(-season) |>
  cbind(total_asv_match_remei_miau) |>
  dplyr::select(season...1, treatment...2, total_common_asv, total_common_asv_potential_bloomers)

colnames(total_asv_match_remei_miau_potential_bloomers) <- c('season', 'treatment', 'total_common_asv',
                                                             'total_common_asv_potential_bloomers')

rare_reactive_asv <- seqs_gr_rel_abund |>
  filter(rel_abund_ed < 0.001 & slope_chosen_days > 2) |>
  group_by(season, treatment) |>
  dplyr::summarize(rare_reactive_asv_sum = n())

summarize_match_remei_miau <- total_asv_match_remei_miau_potential_bloomers |>
  full_join(rare_reactive_asv, by = c('season' = 'season', 'treatment')) |>
  mutate(rare_reactive_perc = round(rare_reactive_asv_sum/total_common_asv_potential_bloomers*100, 2))
summarize_match_remei_miau <- summarize_match_remei_miau |>
  mutate(rare_reactive_perc_ed = ifelse(is.na(rare_reactive_perc), 0, rare_reactive_perc),
         rare_reactive_asv_sum = ifelse(is.na(rare_reactive_asv_sum), 0, rare_reactive_asv_sum))

growing_asv <- reg_all_slopes_chosen_silva_tax |>
  group_by(season, treatment) |>
  dplyr::summarize(total_growing_asv = n())

insitu_asv <- miau_seqs_filt_rel_abund_long |>
  filter(rel_abund > 0) |>
  group_by(season) |>
  dplyr::summarize(insitu_present_asv = n())

summarize_insitu_vs_growing <- summarize_match_remei_miau |>
  left_join(growing_asv, by = c('season' = 'season', 'treatment' = 'treatment')) |>
  left_join(insitu_asv, by = c('season' = 'season')) |>
  dplyr::select(-rare_reactive_perc)

#write.table(summarize_insitu_vs_growing, 'results/tables/dplyr::summarize_insitu_vs_growing_ed2.txt', sep='\t')


## plot relative abundance vs. growth rate (figure S9) -----
rel_abund_gr_relation_rank <- 
  seqs_gr_rel_abund |>
  mutate(bloomers_true = case_when(slope_chosen_days >2  ~ 'bloomers',
                                   slope_chosen_days <2 & rel_abund >= 0.01 ~ 'no')) |>
  ggplot(aes(log10(rank_abund_2), slope_chosen_days, color = ifelse(bloomers_true == 'bloomers', class.x, NA)))+ 
  geom_point(aes(shape =  behaviour), size = 2.5, alpha = 0.8)+
  geom_rect(aes(xmin = 0, xmax = Inf, ymin = 0, ymax = 2), 
            fill = '#D0CFC8', show.legend = FALSE, linejoin = NULL, color = NA, alpha = 0.02)+
  geom_rect(aes(xmin = 0, xmax = log10(rank_abund_2[match(0.01, round(rel_abund_ed2, digits = 2))]), ymin = 2, ymax = Inf), 
            fill = '#D0CFC8', show.legend = FALSE, linejoin = NULL, color = NA, alpha = 0.02)+
  geom_vline(aes(xintercept = log10(rank_abund_2[match(0.01, round(rel_abund_ed2, digits = 2))])), linetype = 'dashed')+
  geom_vline(aes(xintercept = log10(rank_abund_2[match(0.001, round(rel_abund_ed2, digits = 3))])), linetype = 'dashed')+
  geom_hline(yintercept = 2, linetype = 'dashed')+
  labs(y = expression("Growth rate (d"^"-1" *')'), 
       x = 'Log10(rank abundance)', color = 'Class',
       shape = 'In situ\nrelative\nabundance (%)')+
  guides(color=guide_legend(ncol = 4, size = 3),
         shape = guide_legend(ncol = 1, size = 3))+
  scale_color_manual(values = palette_class_assigned)+
  scale_x_continuous(expand = c(0, 0), limits = c(0, 3.6))+
  facet_grid(~season~treatment)+
  geom_text(data = summarize_match_remei_miau, mapping = aes(x = 3.25, y = 2.8, label = paste(round(rare_reactive_perc_ed, 0),'\n%')), 
            check_overlap = TRUE, na.rm = TRUE, show.legend = FALSE, nudge_x = 0, nudge_y = 6, 
            color = 'black', size = 2.5)+
  theme_bw()+
  theme(strip.text = element_text(size = 10), axis.text.x = element_text(angle = 0), legend.position = "bottom",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.text = element_text(size = 7),
        legend.title = element_text(size = 10), strip.background = element_blank(), #strip.text.x = element_blank(),
        panel.border = element_blank()) 

# ggsave(filename = 'rel_abund_gr_cluster100_rank_0rel_abund_ed4.pdf', 
#        plot = rel_abund_gr_relation_rank, device = NULL, 
#        path = 'results/figures/corrected_asv_num/',
#        width = 220, height = 188, units = 'mm')


## plot relative abundance vs. growth rate (figure 7) all_together figure-------
seqs_gr_rel_abund <- seqs_gr_rel_abund |>
  ungroup() |>
  mutate(rel_abund_ed = as.numeric(rel_abund))  |>
  mutate(rel_abund_ed2 = round(rel_abund, 3)) |>
  arrange(rel_abund_ed) |>
  mutate(rank_abund_2 = rev(row_number(rel_abund_ed)),
    behaviour = case_when(rel_abund_ed < 0.001 ~ 'Rare x < 0.1% ',
                          rel_abund_ed > 0.01 ~ 'Abundant x > 1%',
                          rel_abund_ed <0.01  & rel_abund >0.001 ~ 'Mid  0.1% < x < 1%')) 

seqs_gr_rel_abund$treatment <- seqs_gr_rel_abund$treatment |> 
  factor(levels=(c("CL", "CD", "PL", "PD", "DL", "VL", 'NA')))
seqs_gr_rel_abund$season <- seqs_gr_rel_abund$season |> 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))

##table with % of rare reactive ASVs (add to plots as geom_text) 
##from the total potential most-responsive <2 GR & <0.1% community
total_asv_match_remei_miau_potential_bloomers <- seqs_gr_rel_abund |>
  filter(rel_abund_ed < 0.01 & slope_chosen_days > 2) %$%
  asv_num |>
  unique() |>
  length() |>
  as_tibble_col(column_name = 'total_asv_match_remei_miau_potential_bloomers')

total_asv_match_remei_miau <- seqs_gr_rel_abund %$%
  asv_num |>
  unique() |>
  length() |>
  as_tibble_col(column_name = 'total_asv_match_remei_miau') |>
  cbind(total_asv_match_remei_miau_potential_bloomers) 

rare_reactive_asv <- seqs_gr_rel_abund |>
  filter(rel_abund_ed < 0.001 & slope_chosen_days > 2) %$%
  asv_num |>
  unique() |>
  length() |>
  as_tibble_col(column_name = 'rare_reactive_asv')

summarize_match_remei_miau <- total_asv_match_remei_miau |>
  cbind(rare_reactive_asv) |>
  dplyr::mutate(rare_reactive_perc = round(rare_reactive_asv/total_asv_match_remei_miau_potential_bloomers*100, 2))

rel_abund_gr_relation_rank <- 
  seqs_gr_rel_abund |>
  ggplot(aes(log10(rank_abund_2), slope_chosen_days))+ 
  geom_point(aes(color = '929292'), size = 1, alpha = 0.7)+
  geom_rect(aes(xmin = 0, xmax = Inf, ymin = 0, ymax = 2), 
            fill = '#D0CFC8', show.legend = FALSE, linejoin = NULL, color = NA, alpha = 0.02)+
  geom_rect(aes(xmin = 0, xmax = log10(rank_abund_2[match(0.01, round(rel_abund_ed2, digits = 2))]), ymin = 2, ymax = Inf), 
            fill = '#D0CFC8', show.legend = FALSE, linejoin = NULL, color = NA, alpha = 0.02)+
  geom_vline(aes(xintercept = log10(rank_abund_2[match(0.01, round(rel_abund_ed2, digits = 2))])), linetype = 'dashed')+
  geom_vline(aes(xintercept = log10(rank_abund_2[match(0.001, round(rel_abund_ed2, digits = 3))])), linetype = 'dashed')+ 
  geom_hline(yintercept = 2, linetype = 'dashed')+
  scale_x_continuous(expand = c(0, 0), limits = c(0, 3.4))+
  scale_y_continuous(expand = c(0,0), limits = c(0, 10.2))+
  labs(y = expression("Growth rate (d"^"-1" *')'), 
       x = 'Log10(rank abundance)', color = 'Class')+
  guides(color=guide_legend(ncol = 2, size = 7),
         shape = guide_legend(ncol = 2, size = 3))+
  scale_color_manual(values = palette_class_assigned)+
  geom_text(data = summarize_match_remei_miau, mapping = aes(x = 2.8, y = 3, label = paste(round(rare_reactive_perc, 0),'\n%')), 
            check_overlap = TRUE, na.rm = TRUE, show.legend = FALSE, nudge_x = 0.35, nudge_y = 6, 
            color = 'black', size = 5)+
  theme_bw()+
  theme(strip.text = element_text(size = 7), axis.text.x = element_text(angle = 0), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.text = element_text(size = 6),
        legend.title = element_text(size = 7), strip.background = element_blank(), #strip.text.x = element_blank(),
        panel.border = element_blank(), legend.key.size = unit(4, 'mm'), aspect.ratio = 38/39) 

# ggsave(filename = 'rel_abund_gr_cluster100_rank_0rel_abund_all_together_v2.pdf',
#        plot = rel_abund_gr_relation_rank, device = NULL, path = 'results/figures/corrected_asv_num/',
#        width = 188, height = 160, units = 'mm')

## plot most responsive ASVs GR>2, relative abundance t0 <1% and tf >1%  (figure S10)---------

abundant_asv_tf <- rem_relabun_melt |> 
  dplyr::filter(time %in% c('t3', 't4') &
           Abundance > 0.01) |> #more than 1% of the community at time 3 or t4
  dplyr::select(asv_num, Abundance, treatment, season, time) |>
  group_by(asv_num, treatment, season, time) |>
  dplyr::summarize(mean_abundance = mean(Abundance)) |>
  group_by(asv_num, treatment, season, time, mean_abundance) |>
  distinct() |>
  as_tibble() |>
  pivot_wider(names_from = time, values_from = mean_abundance, id_cols = c(asv_num, treatment, season)) |>
  dplyr::mutate(high_mean = case_when (t3 > t4 ~ t3,
                                t4 > t3 ~ t4, 
                                is.na(t3) ~ t4,
                                is.na(t4) ~ t3))

abundant_asv_t0 <- 
  rem_relabun_melt |> 
  dplyr::filter(time == 't0') |>
  dplyr::select(asv_num, Abundance, treatment, season, time) |>
  group_by(asv_num, treatment, season, time) |>
  dplyr::summarize(mean_abundance = mean(Abundance)) |>
  dplyr::filter(mean_abundance  < 0.01)

abundant_asv_tf_t0 <- abundant_asv_tf |>
  left_join(abundant_asv_t0, by = c('season', 'treatment', 'asv_num'))

reg_all_slopes_chosen_silva_tax_tf <- reg_all_slopes_chosen_silva_tax |> 
  right_join(abundant_asv_tf_t0, by = c('treatment' = 'treatment', 'season' = 'season', 'asv_num' = 'asv_num'), copy = TRUE) |> 
  mutate(rank_abund = rank(-high_mean)) |>
  distinct()

seqs_gr_rel_abund_filt <- seqs_gr_rel_abund |> 
  dplyr::select(asv_num, season, treatment, slope_chosen_days, rel_abund) |>
  distinct()

reg_all_slopes_chosen_silva_tax_tf |>
  colnames() 

real_bloomers_natural_community <- reg_all_slopes_chosen_silva_tax_tf |>
  dplyr::select(treatment, season, asv_num, slope_chosen_days, class, order, family, high_mean) |>
  left_join(seqs_gr_rel_abund_filt, by = c('asv_num' = 'asv_num', 
                                           'season' = 'season', 
                                           'treatment' = 'treatment',
                                           'slope_chosen_days' = 'slope_chosen_days')) |>
  distinct() |>
  mutate(rank_natural_community = rank(-rel_abund))

num_real_bloomers_natural_community <- real_bloomers_natural_community |>
  group_by(asv_num) |>
  dplyr::summarize(num_real_bloomers = n()) |>
  dplyr::summarize(num_real_bloomers = n())

summarize_match_remei_miau <- summarize_match_remei_miau |>
  cbind(num_real_bloomers_natural_community) |>
  dplyr::mutate(real_bloomers_percentage = num_real_bloomers/total_asv_match_remei_miau_potential_bloomers)

real_bloomers_experiments_remei_nat_com <- real_bloomers_natural_community |>
  filter(treatment != is.na(treatment)) |>
  mutate(bloomers_true = case_when(slope_chosen_days >2   ~ 'bloomers',
                                   slope_chosen_days <2 & rel_abund >= 0.01 ~ 'no')) |>
  ggplot(aes(log10(rank_natural_community), slope_chosen_days, color = ifelse(bloomers_true == 'bloomers', class, NA)))+ 
  geom_text(data = summarize_match_remei_miau, mapping = aes(x = 2.6, y = 2.8, label = paste(round(real_bloomers_percentage*100, 0), '%')), 
            check_overlap = TRUE, na.rm = TRUE, show.legend = FALSE, nudge_x = 0, nudge_y = 6.4, 
            color = 'black', size = 2.5)+
  geom_point(aes(), size = 1.5, alpha = 0.8)+#shape =  behaviour
  geom_rect(aes(xmin = 0, xmax = Inf, ymin = 0, ymax = 2), 
            fill = '#D0CFC8', show.legend = FALSE, linejoin = NULL, color = NA, alpha = 0.02)+
  geom_rect(aes(xmin = 0, xmax = log10(rank_natural_community[match(0.01, round(rel_abund, digits = 2))]), ymin = 2, ymax = Inf), 
            fill = '#D0CFC8', show.legend = FALSE, linejoin = NULL, color = NA, alpha = 0.02)+
  geom_vline(aes(xintercept = log10(rank_natural_community[match(0.01, round(rel_abund, digits = 2))])), linetype = 'dashed')+
  geom_vline(aes(xintercept = log10(rank_natural_community[match(0.001, round(rel_abund, digits = 3))])), linetype = 'dashed')+
  geom_hline(yintercept = 2, linetype = 'dashed')+
  labs(y = expression("Growth rate (d"^"-1" *')'), 
       x = 'Log10(rank abundance)', color = 'Class',
       shape = 'In situ\nrelative\nlog(rank_natural_community_abund) (%)')+ 
  guides(color=guide_legend(ncol = 3, size = 5),
         shape = guide_legend(ncol = 1, size = 2))+
  scale_color_manual(values = palette_class_assigned)+
  scale_x_continuous(expand = c(0, 0), limits = c(0, 2.8))+
theme_bw()+
  theme(strip.text = element_text(size = 7), axis.text.x = element_text(angle = 0, size = 0), legend.position = "bottom",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.text = element_text(size = 6),
        legend.title = element_text(size = 6), strip.background = element_blank(), #strip.text.x = element_blank(),
        panel.border = element_blank(), legend.key.size = unit(1, 'mm'), aspect.ratio = 38/39)

# ggsave(filename = 'real_bloomers_experiments_remei_natural_community_alltogether_ed3.pdf', plot = real_bloomers_experiments_remei_nat_com, device = NULL, path = 'results/figures/corrected_asv_num/',
#        width = 88, height = 88, units = 'mm')

##summarize class of real bloomers----
real_bloomers_natural_community <- real_bloomers_natural_community |>
  dplyr::filter(slope_chosen_days > 2 &
                  rel_abund < 0.001) |>
  dplyr::select(-rank_natural_community)

write.csv(real_bloomers_natural_community, 'results/tables/real_bloomers_natural_community.csv')

real_bloomers_natural_community |>
  colnames()

real_bloomers_family <- real_bloomers_natural_community |>
  dplyr::filter(slope_chosen_days > 2 &
                  rel_abund < 0.001) |>
  group_by(family) |>
  dplyr::filter(n() >= 5) |>
  ggplot(aes(rel_abund, slope_chosen_days, size = high_mean))+ 
  # geom_text(data = summarize_match_remei_miau, mapping = aes(x = 2.6, y = 2.8, label = paste(round(real_bloomers_percentage*100, 0), '%')),
  #           check_overlap = TRUE, na.rm = TRUE, show.legend = FALSE, nudge_x = 0, nudge_y = 6.4,
  #           color = 'black', size = 2.5)+
  geom_point(aes(shape = treatment, color = season), alpha = 0.8)+#shape =  behaviour
  # geom_rect(aes(xmin = 0, xmax = Inf, ymin = 0, ymax = 2), 
  #           fill = '#D0CFC8', show.legend = FALSE, linejoin = NULL, color = NA, alpha = 0.02)+
  # geom_rect(aes(xmin = 0, xmax = log10(rank_natural_community[match(0.01, round(rel_abund, digits = 2))]), ymin = 2, ymax = Inf), 
  #           fill = '#D0CFC8', show.legend = FALSE, linejoin = NULL, color = NA, alpha = 0.02)+
  # geom_vline(aes(xintercept = log10(rank_natural_community[match(0.01, round(rel_abund, digits = 2))])), linetype = 'dashed')+
  # geom_vline(aes(xintercept = log10(rank_natural_community[match(0.001, round(rel_abund, digits = 3))])), linetype = 'dashed')+
  # geom_hline(yintercept = 2, linetype = 'dashed')+
  labs(y = expression("Growth rate (d"^"-1" *')'), 
       x = 'Natural community relative abundance', color = 'Class',
       shape = 'Treatment',
       size = 'Relative abundance tf')+ 
  scale_shape_manual(values = shapes_treatments)+
  guides(color=guide_legend(ncol = 3, size = 5),
         shape = guide_legend(ncol = 1, size = 2))+
  scale_color_manual(values = palette_seasons_4)+
  facet_wrap(vars(family))+
  scale_x_continuous(expand = c(0, 0), labels = percent_format())+
  theme_bw()+
  theme(strip.text = element_text(size = 7), axis.text.x = element_text(angle = 0, size = 6), legend.position = "bottom",
        #panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), legend.text = element_text(size = 6),
        legend.title = element_text(size = 6), strip.background = element_blank(), #strip.text.x = element_blank(),
        panel.border = element_blank(), legend.key.size = unit(1, 'mm'), aspect.ratio = 38/39)


ggsave('real_bloomers_family.pdf', real_bloomers_family, path = 'results/figures/corrected_asv_num/',
        width = 150,
       height = 150,
       units = 'mm')



# ----- 10. CORRELATIONS WITH CARD-FISH BASED GROWTH RATES ----------------
## abundance-based gr dataset
abundance_gr <- read.delim2("data/envdata/GR_DAPIS_OS/GR_REMEI_DAPIS_OS_Ed.csv", sep = ";") |>
  as_tibble()

abundance_gr <- abundance_gr |>
  dplyr::filter(treatment %in% c("CL", "CD", "PL", "PD",  "DL", "VL")) |>
  dplyr::select(-AAPS) |>
  dplyr::mutate_at(c('PRO', 'EUB', 'ROSEO', 'ALT', 'GAMMA', 'NOR5', 'CFB', 'SAR11'), as.numeric) |>
  dplyr::mutate(experiment = 'REMEI') |>
  pivot_longer(cols = c('PRO', 'EUB', 'ROSEO', 'ALT', 'GAMMA', 'NOR5', 'CFB', 'SAR11'), names_to = 'taxonomy', values_to = 'gr_day') |>
  dplyr::filter(gr_day > 0)

## Alteromonas ASV-based growth rates vs. probe: ALT1413----
asv_alteromonas <- reg_all_slopes_chosen_silva_tax |>
  #dplyr::filter(family_f %in% c('Alteromonadaceae', 'Colwelliaceae'))|> 
  dplyr::filter(genus %in% c('Alteromonas', 'Colwellia', 'Glaciecola')) |>
  group_by(season, treatment) |>
  dplyr::summarize(mean_asv_alt = mean(slope_chosen_days))

dapis_alt <-  abundance_gr |>
  dplyr::filter(taxonomy == 'ALT') |>
  group_by(treatment, season) |>
  dplyr::summarize(mean_dapis_alt = mean(as.numeric(gr_day)))

gr_alt <-  asv_alteromonas |>
  left_join(dapis_alt) 

# gr_alt |>
#   ggplot(aes(mean_asv_alt, mean_dapis_alt))+
#   geom_point(aes(color = season))+
#   geom_smooth(method = 'lm')+
#   scale_color_manual(values = palette_seasons_4)+
#   theme_bw()

## Gammaproteobacteria ASV-based growth rates vs. probe: GAM42a-----
asv_gamma <-reg_all_slopes_chosen_silva_tax |>
  dplyr::filter(class == 'Gammaproteobacteria') |>
  group_by(season, treatment) |>
  dplyr::summarize(mean_asv_gam = mean(slope_chosen_days))

dapis_gamma <- abundance_gr |>
  dplyr::filter(taxonomy == 'GAMMA') |>
  group_by(treatment, season) |>
  dplyr::summarize(mean_dapis_gam = mean(as.numeric(gr_day)))

gamma_gr <- asv_gamma |>
  left_join(dapis_gamma) 

# gamma_gr |>
#   ggplot(aes(mean_asv_gam, mean_dapis_gam, color = season))+
#   geom_point()+
#   scale_color_manual(values = palette_seasons_4)+
#   theme_bw()

##Rhodobacteraceae ASV-based growth rates vs. probe: ROS537 ----
asv_rhodo <- reg_all_slopes_chosen_silva_tax  |>
  dplyr::filter(family %in% c('Rhodobacteraceae')) |># & #
  group_by(season, treatment) |>
  dplyr::summarize(mean_asv_rhodo = mean(slope_chosen_days))

dapis_rhodo <-  abundance_gr |>
  dplyr::filter(taxonomy == 'ROSEO') |>
  group_by(treatment, season) |>
  dplyr::summarize(mean_dapis_rhodo = mean(as.numeric(gr_day)))

rhodo_gr <-  asv_rhodo |>
  right_join(dapis_rhodo) 

# rhodo_gr |>
#   ggplot(aes(mean_asv_rhodo, mean_dapis_rhodo, color = season))+
#   geom_point()+
#   scale_color_manual(values = palette_seasons_4)+
#   theme_bw()

##Eubacteria ASV-based growth rates vs. probe: Eubacteria-I,-II,-III -----
asv_eubacateria <-  reg_all_slopes_chosen_silva_tax  |>
  dplyr::filter(!str_detect(domain, "Archaea")) |>
  group_by(season, treatment) |>
  dplyr::summarize(mean_asv_eub = mean(slope_chosen_days))

dapis_eubacteria <-  abundance_gr |>
  dplyr::filter(taxonomy == 'EUB') |>
  group_by(treatment, season) |>
  dplyr::summarize(mean_dapis_eub = mean(as.numeric(gr_day)))

eub_gr <- asv_eubacateria |>
  right_join(dapis_eubacteria) 

eub_gr$treatment <- factor(eub_gr$treatment,
                           levels = c('CL', 'CD', 'PD', 'PL', 'DL', 'VL'))

eub_gr$season <- factor(eub_gr$season,
                        levels = c('Winter', 'Spring', 'Summer', 'Fall')) 

eub_asv_gr <- eub_gr |>
  ggplot(aes(mean_dapis_eub, mean_asv_eub))+
  geom_point(aes(shape = treatment, color = season))+
  scale_fill_manual(values = palette_seasons_4)+
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0))+
  scale_y_continuous(limits = c(0, 4.2), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  scale_shape_manual(values = shapes_treatments)+
  stat_cor(aes(label = paste(..p.label..)),label.x = 0.25, p.digits = 1, r.digits = 2, p.accuracy = 0.001)+
  guides(size = 'none')+
  labs(y = 'Mean growth rates ASV-based \nEubacteria', 
       x = 'Mean growth rates abundance-based\nEubacteria-I,-II,-III probes', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/4), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = 'right', text= element_text(size = 5))

# ggsave('eub_asv_gr_ed5.pdf', eub_asv_gr,
#        path = "~/Documentos/Doctorat/REMEI/results/figures/corrected_asv_num/",
#        width = 88,
#        height = 88,
#        units = 'mm')

##NOR5-730 (NOR5/OM60) ASV-based growth rates vs. probe NOR5-730 ----
asv_nor5 <-  reg_all_slopes_chosen_silva_tax  |>
  dplyr::filter(family == 'Halieaceae') |>
  group_by(season, treatment) |>
  dplyr::summarize(mean_asv_nor5 = mean(slope_chosen_days))

dapis_nor5 <-  abundance_gr |>
  dplyr::filter(taxonomy == 'NOR5') |>
  group_by(treatment, season) |>
  dplyr::summarize(mean_dapis_nor5 = mean(as.numeric(gr_day)))

nor5_gr <-  asv_nor5 |>
  right_join(dapis_nor5) 

# nor5_gr|>
#   ggplot(aes(mean_asv_nor5, mean_dapis_nor5, color = season))+
#   geom_point()+
#   scale_color_manual(values = palette_seasons_4)+
#   theme_bw()

##Bacteroidetess (Cytophaga, Flavobacteria, Sphingobacteria, Empedobacter, Chryseobacterium, Bergeyella) vs. CF319a-----
asv_cfb <- reg_all_slopes_chosen_silva_tax  |>
  dplyr::filter(family %in% c('Flavobacteriaceae', 'Cytophagaceae', 'NS11-12 marine group') ) |>
  group_by(season, treatment) |>
  dplyr::summarize(mean_asv_cfb = mean(slope_chosen_days))

dapis_cfb <-  abundance_gr %>%
  dplyr::filter(taxonomy == 'CFB') |>
  group_by(treatment, season) |>
  dplyr::summarize(mean_dapis_cfb = mean(as.numeric(gr_day))) |>
  dplyr::filter(mean_dapis_cfb >0)

cfb_gr <- asv_cfb %>%
  left_join(dapis_cfb) 

# cfb_gr |> 
#   ggplot(aes(mean_asv_cfb, mean_dapis_cfb, color = season))+
#   geom_point()+
#   scale_color_manual(values = palette_seasons_4)+
#   theme_bw()

### SAR11 clade vs SAR11 probe
asv_sar11 <-  reg_all_slopes_chosen_silva_tax |>
  dplyr::filter(str_detect(order, "SAR11 clade") &
                  family %in% c('Clade I', 'Clade II') #', NA
  ) %>%
  group_by(season, treatment) %>%
  dplyr::summarize(mean_asv_sar11 = mean(slope_chosen_days))

dapis_sar11 <-  abundance_gr %>%
  dplyr::select(treatment, replicate, season, taxonomy, gr_day) %>%
  dplyr::filter(gr_day >0) |>
  dplyr::filter(taxonomy == 'SAR11') |>
  group_by(treatment, season) %>%
  dplyr::summarize(mean_dapis_sar11 = mean(as.numeric(gr_day)))

sar11_gr <- asv_sar11 %>%
  left_join(dapis_sar11) 

# sar11_411_clade_asv_gr <- sar11_gr %>%
#   #filter(season != 'Summer') %>%
#   ggplot(aes(mean_dapis_sar11, mean_asv))+
#   geom_point(aes(shape = treatment, color = season, size = 0.5))+
#   scale_fill_manual(values = palette_seasons_4)+
#   #scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
#   # scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
#   scale_color_manual(values = palette_seasons_4)+
#   stat_poly_line(color = 'black')+
#   #stat_poly_eq(color = 'black')+
#   stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
#   guides(size = 'none')+
#   labs(y = 'Mean growth rates ASV-based \nSAR11 clade I & II', 
#        x = 'Mean growth rates abundance-based SAR11-441R', color = 'Season', shape = 'Treatment')+
#   theme_bw()+
#   theme(aspect.ratio = (4/4), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         legend.position = 'none', text= element_text(size = 5))

## Correlation table for supplementary material (with no-normalized data)----
## edit colnames to identify which group of ASV are we referring to:
#colnames(cfb_gr) <- c("season", "treatment", "mean_asv_cfb","mean_dapis_cfb")
eub_gr |>
  dim() ##has the 24 conditions I use this dataframe as a template

corr_dapis_asv_gr <- eub_gr |>
  left_join(gamma_gr) |>
  left_join(gr_alt) |>
  left_join(rhodo_gr) |>
  left_join(nor5_gr) |>
  left_join(cfb_gr) |>
  left_join(sar11_gr)

##considerations before computing the correlation (ARE THEY FOLLOWING A NORMAL DISTRIBUTION?)
corr_dapis_asv_gr |>
  str()

# Shapiro-Wilk normality test for all variables pvalue > 0.05 normality can be assumed, if pvalue < 0.05 then no-normality
shapiro.test(as.numeric(corr_dapis_asv_gr$mean_asv_eub)) # => p-value = 0.1724
ggqqplot(as.numeric(corr_dapis_asv_gr$mean_asv_eub))

shapiro.test(as.numeric(corr_dapis_asv_gr$mean_dapis_eub)) # =>  p-value = 0.06861
ggqqplot(as.numeric(corr_dapis_asv_gr$mean_asv_eub))

shapiro.test(as.numeric(corr_dapis_asv_gr$mean_asv_gam)) # =>  p-value = 0.06861
ggqqplot(as.numeric(corr_dapis_asv_gr$mean_asv_gam))

shapiro.test(as.numeric(corr_dapis_asv_gr$mean_dapis_gam)) # => p-value = 0.0005865 (NO-NORMALITY)
ggqqplot(as.numeric(corr_dapis_asv_gr$mean_dapis_gam))

shapiro.test(as.numeric(corr_dapis_asv_gr$mean_asv_alt)) # => p-value = 0.0866
ggqqplot(as.numeric(corr_dapis_asv_gr$mean_asv_alt))

shapiro.test(as.numeric(corr_dapis_asv_gr$mean_dapis_alt)) # =>  p-value = 0.002147(NO-NORMALITY)
ggqqplot(as.numeric(corr_dapis_asv_gr$mean_asv_alt))

shapiro.test(as.numeric(corr_dapis_asv_gr$mean_asv_rhodo)) # =>  p-value = 0.4196
ggqqplot(as.numeric(corr_dapis_asv_gr$mean_asv_rhodo))

shapiro.test(as.numeric(corr_dapis_asv_gr$mean_dapis_rhodo)) # => p-value = 0.4506
ggqqplot(as.numeric(corr_dapis_asv_gr$mean_dapis_rhodo))

shapiro.test(as.numeric(corr_dapis_asv_gr$mean_asv_cfb)) # => p-value = 0.116
ggqqplot(as.numeric(corr_dapis_asv_gr$mean_asv_cfb))

shapiro.test(as.numeric(corr_dapis_asv_gr$mean_dapis_cfb)) # => p-value = 0.003669 (NO-NORMALITY)
ggqqplot(as.numeric(corr_dapis_asv_gr$mean_dapis_cfb))

shapiro.test(as.numeric(corr_dapis_asv_gr$mean_asv_nor5)) # => p-value = 0.02917
ggqqplot(as.numeric(corr_dapis_asv_gr$mean_asv_nor5))

shapiro.test(as.numeric(corr_dapis_asv_gr$mean_dapis_nor5)) # => p-value = 0.07689
ggqqplot(as.numeric(corr_dapis_asv_gr$mean_dapis_nor5))

shapiro.test(as.numeric(corr_dapis_asv_gr$mean_asv_sar11)) # => p-value = 0.7076
ggqqplot(as.numeric(corr_dapis_asv_gr$mean_asv_sar11))

shapiro.test(as.numeric(corr_dapis_asv_gr$mean_dapis_sar11)) # => p-value = 0.06008
ggqqplot(as.numeric(corr_dapis_asv_gr$mean_dapis_sar11))

### All at the same time----
## function to extract a table
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

cor_pearson <- rcorr(as.matrix(corr_dapis_asv_gr[,3:16]), type = 'pearson')

View(cor_pearson$n)

my_cor_matrix <- flat_cor_mat(cor_pearson$r, cor_pearson$P)
head(my_cor_matrix)

my_cor_matrix_pear <- my_cor_matrix |>
  as_tibble() |>
  # filter(row %in% c('mean_dapis_gamma', 'mean_asv_gamma', 'mean_dapis_nor5', 'mean_asv_nor5')) |>
  # filter(column != c('mean_dapis_gamma', 'mean_asv_gamma', 'mean_dapis_nor5', 'mean_asv_nor5')) |>
  dplyr::filter(row == 'mean_asv_eub' & column == 'mean_dapis_eub' |
                  row == 'mean_asv_rhodo' & column == 'mean_dapis_rhodo' |
                  row == 'mean_asv_nor5' & column == 'mean_dapis_nor5' |
                  row == 'mean_asv_sar11' & column == 'mean_dapis_sar11')

##for non-normal data (pvalue < 0.05 in Shapiro test)
cor_spearman <- rcorr(as.matrix(corr_dapis_asv_gr[,3:16]), type = 'spearman')

my_cor_matrix_spear <- flat_cor_mat(cor_spearman$r, cor_spearman$P)
head(my_cor_matrix_spear)

my_cor_matrix_spear <- my_cor_matrix_spear |>
  as_tibble() |>
  # filter(row %in% c('mean_dapis_gamma', 'mean_asv_gamma', 'mean_dapis_nor5', 'mean_asv_nor5')) |>
  # filter(column != c('mean_dapis_gamma', 'mean_asv_gamma', 'mean_dapis_nor5', 'mean_asv_nor5')) |>
  filter(row == 'mean_asv_alt' & column == 'mean_dapis_alt' |
           row == 'mean_asv_gam' & column == 'mean_dapis_gam' |
           row == 'mean_asv_cfb' & column == 'mean_dapis_cfb' )

##table with all the results for the correlations
my_correlation <- my_cor_matrix_pear |>
  rbind(my_cor_matrix_spear)

colnames(my_correlation) = c('Mean ASV-based growth rates', 'Abundance-based growth rates', 'Correlation', 'p-value')

#write.csv2(my_correlation, file = 'results/tables/correlations_card-fish_avs__all_gr.csv')
