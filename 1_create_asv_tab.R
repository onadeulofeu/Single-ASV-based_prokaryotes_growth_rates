##sessionInfo() versions of the loaded packages when this script was ran
#[1] Biobase_2.50.0      bit64_4.0.5         jsonlite_1.7.2      splines_4.0.5       foreach_1.5.1      
#[6] assertthat_0.2.1    stats4_4.0.5        blob_1.2.1          remotes_2.3.0       progress_1.2.2     
#[11] pillar_1.6.4        RSQLite_2.2.7       lattice_0.20-41     glue_1.6.0          XVector_0.30.0     
#[16] colorspace_2.0-2    Matrix_1.3-2        plyr_1.8.6          pkgconfig_2.0.3     zlibbioc_1.36.0    
#[21] purrr_0.3.4         scales_1.1.1        processx_3.5.2      tibble_3.1.6        mgcv_1.8-34        
#[26] generics_0.1.1      IRanges_2.24.1      ggplot2_3.3.5       ellipsis_0.3.2      cachem_1.0.4       
#[31] withr_2.4.3         BiocGenerics_0.36.1 cli_3.1.0           survival_3.2-10     magrittr_2.0.1     
#[36] crayon_1.4.2        memoise_2.0.0       ps_1.6.0            fansi_1.0.0         nlme_3.1-152       
#[41] MASS_7.3-53.1       pkgbuild_1.2.0      vegan_2.5-7         tools_4.0.5         data.table_1.14.2  
#[46] prettyunits_1.1.1   hms_1.0.0           lifecycle_1.0.1     stringr_1.4.0       Rhdf5lib_1.12.1    
#[51] S4Vectors_0.28.1    munsell_0.5.0       cluster_2.1.1       callr_3.7.0         Biostrings_2.58.0  
#[56] ade4_1.7-16         compiler_4.0.5      rlang_0.4.12        rhdf5_2.34.0        grid_4.0.5         
#[61] iterators_1.0.13    rhdf5filters_1.2.0  biomformat_1.18.0   rstudioapi_0.13     igraph_1.2.6       
#[66] gtable_0.3.0        codetools_0.2-18    multtest_2.46.0     DBI_1.1.1           curl_4.3.1         
#[71] reshape2_1.4.4      R6_2.5.1            dplyr_1.0.7         fastmap_1.1.0       bit_4.0.4          
#[76] utf8_1.2.2          rprojroot_2.0.2     permute_0.9-5       ape_5.6-1           stringi_1.7.6      
#[81] parallel_4.0.5      Rcpp_1.0.7          DECIPHER_2.18.1     vctrs_0.3.8         tidyselect_1.1.1 

sessionInfo()
##packages
library(tidyverse)
library(magrittr)
library(lubridate)
library(reshape2)
library(plyr)
library(speedyseq)
library(dplyr)
library(janitor)  #row_to_names

setwd ("~/Documentos/Doctorat/REMEI/")

##If you want to add abundance counts to the environmental data to be able to calculate pseudoabundances
##go to line 187 (DAPI) or 323 (FLOW CITOMETRY)

##load functions
#source("src/ggplot_theme_definition.R") ###same aspect in all graphs (theme definition)

#load input data
##consensus asv table
seqtab <- readRDS("data/dada2/02_nochimera_mergeruns/remei_1_2/remei_1_2_seqtab_final.rds")
tax_gtbd <- readRDS("data/dada2/03_taxonomy/remei_1_2_GTBD/remei_1_2_GTBD_tax_assignation.rds") %>% 
  as_tibble(rownames = 'sequence')
tax_silva <- readRDS("data/dada2/03_taxonomy/remei_1_2_SILVA/remei_1_2_SILVA_tax_assignation.rds") %>% 
  as_tibble(rownames = 'sequence')

##pool asv table (treballem amb aquesta taula, millor per a les asv rares: 
#pooling allows information to be shared across samples, which makes it easier to resolve rare variants that 
#were present as singletons or doubletone in one sample but were present many times across samples
seqtab_pool <- readRDS("data/dada2/02_nochimera_mergeruns/remei_1_2_pool//remei_1_2_pool_seqtab_final.rds")
tax_gtbd_pool <- readRDS("data/dada2/03_taxonomy/remei_1_2_GTBD_pool//remei_1_2_GTBD_pool_tax_assignation.rds") %>% 
  as_tibble(rownames = 'sequence')
tax_silva_pool <- readRDS("data/dada2/03_taxonomy/remei_1_2_silva_pool/remei_1_2_silva_pool_tax_assignation.rds") %>% 
  as_tibble(rownames = 'sequence')

##first work with concensus asv table
dim(seqtab_pool) 
##consensus : 14412 asv i 449 mostres.
##pool: 449 12198 asv

tax.filt.gtbd <- tax_gtbd_pool %>% 
  #Filter all results without domain assignation 
  filter(!is.na(domain), !is.na(phylum)) %>% 
  #And the Chloros/Mitochondria seqs
  filter( order !=  'Chloroplast') %>% 
  filter( family !=  'Mitochondria') 

tax.filt.silva <- tax_silva_pool %>% 
  #Filter all results without domain assignation 
  filter(!is.na(domain), !is.na(phylum)) %>% 
  #And the Chloros/Mitochondria seqs
  filter( order !=  'Chloroplast') %>% 
  filter( family !=  'Mitochondria')

# Environmental data ------------------------------------------------------
exp.data <- readxl::read_xlsx('data/envdata/REMEI_amplicon_sample_codes.xlsx',
                              skip = 1, col_names = T )  %>% 
  # We select only the ones in Selected filenames and without FAILED
  filter(!is.na(SelectedFileName), SelectedFileName != "FAILED") %>%
  filter(!is.na(Type) | !is.na(Treatment) | !is.na(Time)) %>% ##delete samples from another data set
  filter(Type != "Profile") %>%
  filter(Treatment !="ShortMock") %>%
  filter(Treatment !="LongMock") %>%
  filter(!duplicated(`Sample Code`)) ##Hem d'eliminar les mostres repetides que van sortir malament el primer cop 
##i que tenen 2 files a la metadata info
##Mirar que en env.data estiguin primer les  que  ens  volem quedar  (les noves) (potser buscar una manera més neta de fer-ho)
#total de mostres repetides 29, comprovar que les tinguem totes.
exp.data %>% 
  filter(str_detect(SelectedFileName, "6642")) ##29 rows!

##without filtering 369 samples 
dim(exp.data) ##292 samples
#edit column names to have a tidy one's
exp.data %>% 
  colnames()

library(plyr)
exp.data <- exp.data %>% 
  rename(c("Type" = "type", "Treatment" = "treatment", "Light Regime" = "light_regime",
                      "Replicate" = "replicate", "Time" = "time", "Hours" = "hours", "Season" = "season", 
                      "Fraction" = "fraction", "Sample-Name" = "sample_name" ,"Sample Code" = "sample_code",  
                      "SelectedFileName" = "selected_file_name", "Reads"= "reads", "FileName1" = "file_name1",
                      "Reads1" = "reads1", "Filename2" = "file_name2", "Reads2" = "reads2"))

# Filtering ---------------------------------------------------------------
# Keeping only samples and taxonomy of interest
ab.filt_gtbd <- seqtab_pool[ exp.data$selected_file_name, tax.filt.gtbd$sequence]
ab.filt_silva <- seqtab_pool[ exp.data$selected_file_name, tax.filt.silva$sequence]

# How many samples do we have? 
dim(ab.filt_gtbd) ##consensus:292 and 5154 asv, pool: 4599
dim(ab.filt_silva) ##consensus: 292 and 6610 asv,pool: 5670

# Putting control and predator-free t0 in the right place -----------------

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

newexp.data <- change_t0s(exp.data) 

samnames <- filter(newexp.data, !is.na(selected_file_name_new)) %>%
  pull(selected_file_name)

newab <- ab.filt_silva[samnames,] 
rownames(newab) <- filter(newexp.data, !is.na(selected_file_name_new)) %>%
  pull(selected_file_name_new)

old <- ab.filt_silva[!rownames(ab.filt_silva) %in% samnames,]

# The new names 
ab.filt_silva <- rbind(old, newab)

exp.data <- newexp.data %>% 
  mutate( selected_file_name = ifelse(is.na(selected_file_name_new), selected_file_name, 
                                      selected_file_name_new))
dim(newexp.data)
dim(exp.data) ##312 samples
dim(ab.filt_silva)

exp.data %>% 
  sample_names()
ab.filt_silva %>% 
  sample_names()

##If you want to add abundance counts to the environmental data to be able to calculate pseudoabundances-----
##go to line 191 (DAPI) or 326 (FLOW CITOMETRY)

# Phyloseq creation -------------------------------------------------------
# Putting cooler names to the project 
samnames <- exp.data$selected_file_name %>%
  str_split(pattern = "_",simplify = T) %>%
  # Keep only the first column
  .[,1]

length(samnames)
length(unique(samnames))
# Golden is the numbers are equal 
rownames(ab.filt_silva) <- samnames
rownames(exp.data) <- samnames

# Phyloseq
ASV <- otu_table(ab.filt_silva, taxa_are_rows = F)
TAX <- tax_table(as.matrix(tax.filt.silva %>% column_to_rownames(var = 'sequence')))
DAT <- sample_data(exp.data %>% as.data.frame()) 

rem_phy <- phyloseq(ASV,TAX,DAT)
saveRDS(rem_phy, 'data/intermediate_files/remei_phyloseq_pool_silva.rds')



##WE ADD DAPI TO THE EXP.DATA CREATE A PHYLOSEQ OBJECT WITHOUT CRUISE EXPERIMENT----------------------
# Environmental data 
exp.data <- readxl::read_xlsx('data/envdata/REMEI_amplicon_sample_codes.xlsx',
                              skip = 1, col_names = T )  %>% 
  # We select only the ones in Selected filenames and without FAILED
  filter(!is.na(SelectedFileName), SelectedFileName != "FAILED") %>%
  filter(!is.na(Type) | !is.na(Treatment) | !is.na(Time)) %>% ##delete samples from another data set
  filter(Type != "Profile") %>%
  filter(Treatment !="ShortMock") %>%
  filter(Treatment !="LongMock") %>%
  filter(!duplicated(`Sample Code`)) %>% 
  ##Hem d'eliminar les mostres repetides que van sortir malament el primer cop i que tenen 2 files a la metadata info
  filter(Season != "Cruise")

##Mirar que en env.data estiguin primer les  que ens  volem quedar  (les noves) 
##(potser buscar una manera més neta de fer-ho)
#total de mostres repetides 29, comprovar que les tinguem totes.
exp.data %>% filter(str_detect(SelectedFileName, "6642")) ##25 rows! (no hi són les del Cruise)

#edit column names to have a tidy one's
exp.data %>% 
  colnames()

exp.data <- exp.data %>% 
  rename(c("Type" = "type", "Treatment" = "treatment", "Light Regime" = "light_regime",
           "Replicate" = "replicate", "Time" = "time", "Hours" = "hours", "Season" = "season", 
           "Fraction" = "fraction", "Sample-Name" = "sample_name" ,"Sample Code" = "sample_code",  
           "SelectedFileName" = "selected_file_name", "Reads"= "reads", "FileName1" = "file_name1",
           "Reads1" = "reads1", "Filename2" = "file_name2", "Reads2" = "reads2"))

##without filtering 369 samples 
dim(exp.data) ##230 samples (sense el Cruise ara)

# Filtering ---------------------------------------------------------------
# Keeping only samples and taxonomy of interest
ab.filt_silva <- seqtab_pool[ exp.data$selected_file_name, tax.filt.silva$sequence]
ab.filt_silva <- seqtab_pool[ exp.data$selected_file_name, tax.filt.silva$sequence]

# How many samples do we have? 
dim(ab.filt_silva) ##230 and 5154 asv
dim(ab.filt_silva) ##230 and 6610 asv

# Putting control and predator-free t0 in the right place -----------------

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
    mutate( selected_file_name_new = str_c('REMEIextra', 1:nrow(.)))
  
  df %>% 
    filter( !treatment %in% c('CT', 'PF')) %>% 
    bind_rows(whole) %>% 
    return()
}

newexp.data <- change_t0s(exp.data) 

samnames <- filter(newexp.data, !is.na(selected_file_name_new)) %>%
  pull(selected_file_name)

newab <- ab.filt_silva[samnames,] 
rownames(newab) <- filter(newexp.data, !is.na(selected_file_name_new)) %>%
  pull(selected_file_name_new)

old <- ab.filt_silva[!rownames(ab.filt_silva) %in% samnames,]

# The new names 
ab.filt_silva <- rbind(old, newab)

exp.data <- newexp.data %>% 
  mutate( selected_file_name = ifelse(is.na(selected_file_name_new), selected_file_name, 
                                      selected_file_name_new))
dim(newexp.data)
dim(exp.data) ##246 samples
dim(ab.filt_silva) ##246

exp.data %>% 
  sample_names()
ab.filt_silva %>% 
  sample_names()

##add dapi counts to env.data
dapis <- readxl::read_xlsx('data/envdata/DAPI_remei_ed2.xlsx',
                           col_names = T )
exp.data %>% 
  colnames()
dapis %T>% 
  colnames() %>%
  head()
dapis <- dapis %>% 
  rename(c("DAPI" = "dapi", "Treatment" = "treatment", "Light Regime" = "light_regime",
           "Replicate" = "replicate", "Hours_dec" = "hours_dec", "Time" = "time",  "Season" = "season"))

exp.data <- exp.data %>% 
  left_join(dapis, by = c("treatment", "replicate", "time", "season"))
dim(exp.data)
exp.data %>% 
  View()

# Phyloseq creation with prokaryotic abundance data -------------------------------------------------------
# Putting cooler names to the project 
samnames <- exp.data$selected_file_name %>%
  str_split(pattern = "_",simplify = T) %>%
  # Keep only the first column
  .[,1]

length(samnames)
length(unique(samnames))
# Golden is the numbers are equal 

rownames(ab.filt_silva) <- samnames
rownames(exp.data) <- samnames

# Phyloseq
ASV <- otu_table(ab.filt_silva, taxa_are_rows = F)
TAX <- tax_table(as.matrix(tax.filt.silva %>% column_to_rownames(var = 'sequence')))
DAT <- sample_data(exp.data %>% as.data.frame()) 

rem.phy.dapis <- phyloseq(ASV,TAX,DAT)
saveRDS(rem.phy.dapis, 'data/intermediate_files/remei_phyloseq_silva_DAPI.rds')

##Add flow citometry counts to env.data-----
fc <- readxl::read_xlsx('data/envdata/Final_Data_FC_PB_Remei_experiments_ed4.xlsx',
                           col_names = T )

exp.data <- readxl::read_xlsx('data/envdata/REMEI_amplicon_sample_codes.xlsx',
                              skip = 1, col_names = T )  %>% 
  # We select only the ones in Selected filenames and without FAILED
  filter(!is.na(SelectedFileName), SelectedFileName != "FAILED") %>%
  filter(!is.na(Type) | !is.na(Treatment) | !is.na(Time)) %>% ##delete samples from another data set
  filter(Type != "Profile") %>%
  filter(Treatment !="ShortMock") %>%
  filter(Treatment !="LongMock") %>%
  filter(!duplicated(`Sample Code`))
  ##Hem d'eliminar les mostres repetides que van sortir malament el primer cop i que tenen 2 files a la metadata info


##Mirar que en env.data estiguin primer les  que ens  volem quedar  (les noves) 
##(potser buscar una manera més neta de fer-ho)
#total de mostres repetides 29, comprovar que les tinguem totes.
exp.data %>% 
  filter(str_detect(SelectedFileName, "6642")) ##25 rows! (no hi són les del Cruise)

#edit column names to have a tidy one's
exp.data %>% 
  colnames()

exp.data <- exp.data %>% 
  rename(c("Type" = "type", "Treatment" = "treatment", "Light Regime" = "light_regime",
           "Replicate" = "replicate", "Time" = "time", "Hours" = "hours", "Season" = "season", 
           "Fraction" = "fraction", "Sample-Name" = "sample_name" ,"Sample Code" = "sample_code",  
           "SelectedFileName" = "selected_file_name", "Reads"= "reads", "FileName1" = "file_name1",
           "Reads1" = "reads1", "Filename2" = "file_name2", "Reads2" = "reads2"))

##without filtering 369 samples 
dim(exp.data) ##292 samples

# Filtering ---------------------------------------------------------------
# Keeping only samples and taxonomy of interest
ab.filt_gtbd <- seqtab_pool[exp.data$selected_file_name, tax.filt.gtbd$sequence]

# How many samples do we have? 
dim(ab.filt_gtbd) ##292 and 5670 asv


# Putting control and predator-free t0 in the right place -----------------

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
    mutate( selected_file_name_new = str_c('REMEIextra', 1:nrow(.)))
  
  df %>% 
    filter( !treatment %in% c('CT', 'PF')) %>% 
    bind_rows(whole) %>% 
    return()
}

newexp.data <- change_t0s(exp.data) 

samnames <- filter(newexp.data, !is.na(selected_file_name_new)) %>%
  pull(selected_file_name)

newab <- ab.filt_gtbd[samnames,] 
rownames(newab) <- filter(newexp.data, !is.na(selected_file_name_new)) %>%
  pull(selected_file_name_new)

old <- ab.filt_gtbd[!rownames(ab.filt_gtbd) %in% samnames,]

# The new names 
ab.filt_gtbd <- rbind(old, newab)

exp.data <- newexp.data %>% 
  mutate( selected_file_name = ifelse(is.na(selected_file_name_new), selected_file_name, 
                                      selected_file_name_new))
dim(newexp.data)
dim(exp.data) ##312 samples
dim(ab.filt_gtbd) ##312

exp.data %>% 
  sample_names()
ab.filt_gtbd %>% 
  sample_names()

exp.data %>% 
  colnames()
fc %T>% 
  colnames() %>%
  head()
# fc <- fc %>% 
#   rename(c("All" = "fc", "Treatment" = "treatment", "Light Regime" = "light_regime",
#            "Replicate" = "replicate", "Hours_dec" = "hours_dec", "Time" = "time",  "Season" = "season"))

exp.data <- exp.data %>% 
  left_join(fc, by = c("treatment", "replicate", "time", "season"))
dim(exp.data)
# exp.data %>% 
#   View()

# Phyloseq creation -------------------------------------------------------
# Putting cooler names to the project 
samnames <- exp.data$selected_file_name %>%
  str_split(pattern = "_",simplify = T) %>%
  # Keep only the first column
  .[,1]

length(samnames)
length(unique(samnames))
# Golden is the numbers are equal 

rownames(ab.filt_gtbd) <- samnames
rownames(exp.data) <- samnames

# Phyloseq
ASV <- otu_table(ab.filt_gtbd, taxa_are_rows = F)
TAX <- tax_table(as.matrix(tax.filt.gtbd %>% column_to_rownames(var = 'sequence')))
DAT <- sample_data(exp.data %>% as.data.frame()) 

rem.phy.fcs <- phyloseq(ASV,TAX,DAT)
saveRDS(rem.phy.fcs, 'data/intermediate_files/remei_phyloseq_gtbd_fc.rds')

##add t0 experiments from BBMO
#keep samples start with BL

#Create ASV clustered TAB with experiments in situ samples from BBMO---------
##miro primer si tenim el mateixn nº de ASV amb el clustering que sense (tene el mateix per tant utilitzo la unclust)
# seqtab_clust <- readRDS("data/dada2/remei_1_2_miau_1_seqtab_clust100.rds")
# seqtab_clust %>%
#   dim()

seqtab_unclust <- readRDS("data/dada2/02_nochimera_mergeruns/remei_1_2_miau_1_pool/remei_1_2_miau_1_pool_seqtab_final.rds")
seqtab_unclust %>%
  dim() ##711 samples with 33099 ASV

tax_silva_remiau <- readRDS("data/dada2/03_taxonomy/remei_1_2_miau_1_silva_pool/remei_1_2_miau_1_silva_pool_tax_assignation.rds") %>% 
  as_tibble(rownames = 'sequence')

##treballo amb SILVA perquè va millor pels amplicons però si vull comparar-ho amb MAGs tinc també la classifació GTBD
tax_silva_remiau <- tax_silva_remiau %>% 
  #Filter all results without domain assignation 
  filter(!is.na(domain), !is.na(phylum)) %>% 
  #And the Chloros/Mitochondria seqs
  filter( order !=  'Chloroplast') %>% 
  filter( family !=  'Mitochondria')

tax_silva_remiau %>%
dim() ##14325 ASV de 33099

# Environmental data MIAU sequencing but for in situ REMEI experiments ------------------------------------------------------
exp.data_miau <- read.table('/Users/onadeulofeucapo/Documentos/Doctorat/MIAU_seqs/sample_labels.txt', sep = '\t')

##keep only the samples from BBMO in situ
exp.data_miau_filt <- 
  exp.data_miau %>%
  filter(str_detect(V1, 'BL_')) %>%
  mutate(
    season = case_when(.$V1 == 'BL_170220Exp. REMEI-W' ~ 'Winter',
                            .$V1 == 'BL_170425Exp REMEI-SP' ~ 'Spring',
                            .$V1 == 'BL_170704Exp REMEI-SU' ~ 'Summer',
                            .$V1 == 'BL_171106EXP REMEI-Fall' ~ 'Fall'),
         time ='t-2',
         hours = as.double('-2'),
         light_regime = 'L',
         selected_file_name = V3)

exp.data_miau_filt %>%
  head()

##REMEI
exp.data <- readxl::read_xlsx('data/envdata/REMEI_amplicon_sample_codes.xlsx',
                              skip = 1, col_names = T )  %>% 
  # We select only the ones in Selected filenames and without FAILED
  filter(!is.na(SelectedFileName), SelectedFileName != "FAILED") %>%
  filter(!is.na(Type) | !is.na(Treatment) | !is.na(Time)) %>% ##delete samples from another data set
  filter(Type != "Profile") %>%
  filter(Treatment !="ShortMock") %>%
  filter(Treatment !="LongMock") %>%
  filter(!duplicated(`Sample Code`)) ##Hem d'eliminar les mostres repetides que van sortir malament el primer cop 
##i que tenen 2 files a la metadata info
##Mirar que en env.data estiguin primer les  que  ens  volem quedar  (les noves) (potser buscar una manera més neta de fer-ho)
#total de mostres repetides 29, comprovar que les tinguem totes.
exp.data %>% 
  filter(str_detect(SelectedFileName, "6642")) ##29 rows!

##without filtering 369 samples 
dim(exp.data) ##292 samples
#edit column names to have a tidy one's
exp.data %>% 
  colnames()

#activar plyr, desactivar dplyr 
library(plyr)
exp.data <- exp.data %>% 
  rename(c("Type" = "type", "Treatment" = "treatment", "Light Regime" = "light_regime",
           "Replicate" = "replicate", "Time" = "time", "Hours" = "hours", "Season" = "season", 
           "Fraction" = "fraction", "Sample-Name" = "sample_name" ,"Sample Code" = "sample_code",  
           "SelectedFileName" = "selected_file_name", "Reads"= "reads", "FileName1" = "file_name1",
           "Reads1" = "reads1", "Filename2" = "file_name2", "Reads2" = "reads2"))
detach('package:plyr', unload= T)

exp.data_remiau <- exp.data %>% 
  bind_rows(exp.data_miau_filt) ##296 samples REMEI + t0 de MIAU 4
 
exp.data_remiau %>%
  dim()

# Filtering ---------------------------------------------------------------
# Keeping only samples and taxonomy of interest
dim(seqtab_unclust) #711 samples 33099 asv
dim(exp.data_remiau) # 296 samples
dim(tax_silva_remiau) #14325 asv (sense cloroplasts and mitocondria)
#ab.filt_silva <- seqtab_pool[ exp.data$selected_file_name, tax.filt.silva$sequence]
##error subíndice fuera de los límites
#ab.filt_silva <- seqtab_clust[exp.data_remiau$selected_file_name, tax_silva_remiau$sequence]

#filter by samples of interest
seqtab_filt <- seqtab_unclust %>%
  as_tibble(rownames = 'selected_file_name') %>%
  right_join(exp.data_remiau, by = 'selected_file_name') ##296 samples

# seqtab_unclust$selected_file_names
# exp.data_remiau$selected_file_name
# seqtab_unclust %>%
#   rownames()

# exp.data_remiau %$%
#   selected_file_name %>%
#   unique()
# 
# exp.data_remiau %>%
#   dim()
# seqtab_clust_filt %>%
#   colnames()

##filter by taxonomy of interest
ab.filt_silva <- seqtab_filt %>%
  t() %>%
  row_to_names(row_number = 1) %>%
  as_tibble(rownames = 'sequence') %>%
  right_join(tax_silva_remiau, by = 'sequence') %>%
  select(-c(298:304)) %>% #remove cols with taxonomy information
  t() %>%
  row_to_names(row_number = 1) %>%
  as_tibble(rownames = NA)

# How many samples do we have? 
#dim(ab.filt_gtbd) ##consensus:292 and 5154 asv, pool: 4599
dim(ab.filt_silva) ##pool:  296 and 6610 asv,pool: 14325

# Putting control and predator-free t0 in the right place -----------------

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

newexp.data <- change_t0s(exp.data_remiau) 

samnames <- filter(newexp.data, !is.na(selected_file_name_new)) %>%
  pull(selected_file_name)

newab <- ab.filt_silva[samnames,] 
rownames(newab) <- filter(newexp.data, !is.na(selected_file_name_new)) %>%
  pull(selected_file_name_new)

old <- ab.filt_silva[!rownames(ab.filt_silva) %in% samnames,]

# The new names 
ab.filt_silva <- rbind(old, newab)

exp.data <- newexp.data %>% 
  mutate( selected_file_name = ifelse(is.na(selected_file_name_new), selected_file_name, 
                                      selected_file_name_new))
dim(newexp.data)
dim(exp.data) ##316 samples
dim(ab.filt_silva)

exp.data %>% 
  sample_names()
ab.filt_silva %>% 
  sample_names()

# fc %T>% 
#   colnames() %>%
#   head()
# fc <- fc %>% 
#   rename(c("All" = "fc", "Treatment" = "treatment", "Light Regime" = "light_regime",
#            "Replicate" = "replicate", "Hours_dec" = "hours_dec", "Time" = "time",  "Season" = "season"))
fc <- readxl::read_xlsx('data/envdata/Final_Data_FC_PB_Remei_experiments_ed5.xlsx',
                        col_names = T )

exp.data <- exp.data %>% 
  left_join(fc, by = c("treatment", "replicate", "time", "season"))
dim(exp.data)
# exp.data %>% 
#   View()

# Phyloseq creation -------------------------------------------------------
# Putting cooler names to the project 
samnames <- exp.data$selected_file_name %>%
  str_split(pattern = "_",simplify = T) %>%
  # Keep only the first column
  .[,1]

length(samnames)
length(unique(samnames))
# Golden is the numbers are equal 

rownames(ab.filt_silva) <- samnames
rownames(exp.data) <- samnames

ab.filt_silva_ed <- ab.filt_silva %>%
  rownames_to_column(var = 'selected_file_name') %>%
  mutate(across(c(where(is.character), -selected_file_name), as.numeric)) %>%
  as_tibble()

# Phyloseq
ab.filt_silva_ed %>%
  glimpse()
ASV <- otu_table(ab.filt_silva_ed, taxa_are_rows = F)
TAX <- tax_table(as.matrix(tax_silva_remiau %>% column_to_rownames(var = 'sequence')))
DAT <- sample_data(exp.data %>% as.data.frame()) 

remiau.phy.fcs <- phyloseq(ASV,TAX,DAT)
saveRDS(remiau.phy.fcs, 'data/intermediate_files/remiau_phyloseq_silva_fc.rds')

#MIAU DATA subset ALL MIAU experiments (common ASV names with REMEI experiment)------------
exp.data_miau <- read.table('/Users/onadeulofeucapo/Documentos/Doctorat/MIAU_seqs/sample_labels.txt', sep = '\t')

exp.data_miau %>%
  head()

##abundance table
seqtab_unclust <- readRDS("data/dada2/02_nochimera_mergeruns/remei_1_2_miau_1_pool/remei_1_2_miau_1_pool_seqtab_final.rds")
seqtab_unclust %>%
  dim() ##711 samples with 33099 ASV

tax_silva_remiau <- readRDS("data/dada2/03_taxonomy/remei_1_2_miau_1_silva_pool/remei_1_2_miau_1_silva_pool_tax_assignation.rds") %>% 
  as_tibble(rownames = 'sequence')

##treballo amb SILVA perquè va millor pels amplicons però si vull comparar-ho amb MAGs tinc també la classifació GTBD
tax_silva_miau <- tax_silva_remiau %>% 
  #Filter all results without domain assignation 
  filter(!is.na(domain), !is.na(phylum)) %>% 
  #And the Chloros/Mitochondria seqs
  filter( order !=  'Chloroplast') %>% 
  filter( family !=  'Mitochondria')

tax_silva_miau %>%
  dim() ##14325 ASV de 33099

# Environmental data MIAU sequencing but for in situ REMEI experiments ------------------------------------------------------
exp.data_miau <- read.table('/Users/onadeulofeucapo/Documentos/Doctorat/MIAU_seqs/sample_labels.txt', sep = '\t') %>%
  as_tibble()

colnames(exp.data_miau)  = exp.data_miau[1,]
exp.data_miau <- exp.data_miau[-1,]
exp.data_miau %>%
  head()

# Filtering ---------------------------------------------------------------
# Keeping only samples and taxonomy of interest
dim(seqtab_unclust) #711 samples 33099 asv
dim(exp.data_miau) # 263 samples
dim(tax_silva_remiau) #14325 asv (sense cloroplasts and mitocondria)
#ab.filt_silva <- seqtab_pool[ exp.data$selected_file_name, tax.filt.silva$sequence]
##error subíndice fuera de los límites
#ab.filt_silva <- seqtab_clust[exp.data_remiau$selected_file_name, tax_silva_remiau$sequence]

#filter by samples of interest
seqtab_filt <- seqtab_unclust %>%
  as_tibble(rownames = 'selected_file_name') %>%
  right_join(exp.data_miau, by =c('selected_file_name'=  'seq_code')) ##263 samples

seqtab_filt %>%
  dim()
# seqtab_unclust$selected_file_names
# exp.data_remiau$selected_file_name
# seqtab_unclust %>%
#   rownames()

# exp.data_remiau %$%
#   selected_file_name %>%
#   unique()
# 
# exp.data_remiau %>%
#   dim()
# seqtab_clust_filt %>%
#   colnames()

##filter by taxonomy of interest
ab.filt_silva <- seqtab_filt %>%
  t() %>%
  row_to_names(row_number = 1) %>%
  as_tibble(rownames = 'sequence') %>%
  right_join(tax_silva_miau, by = 'sequence') %>%
  select(-c(263:269)) %>% #remove cols with taxonomy information
  t() %>%
  row_to_names(row_number = 1) %>%
  as_tibble(rownames = NA)

# How many samples do we have? 
#dim(ab.filt_gtbd) ##consensus:292 and 5154 asv, pool: 4599
dim(ab.filt_silva) ##pool:  262 and 14325 asv


# fc %T>% 
#   colnames() %>%
#   head()
# fc <- fc %>% 
#   rename(c("All" = "fc", "Treatment" = "treatment", "Light Regime" = "light_regime",
#            "Replicate" = "replicate", "Hours_dec" = "hours_dec", "Time" = "time",  "Season" = "season"))

#PENDENT NO TINC LES CITOMETRIES ENDREÇADES DE MIAU
fc <- readxl::read_xlsx('data/envdata/Final_Data_FC_PB_Remei_experiments_ed4.xlsx',
                        col_names = T )

exp.data <- exp.data_miau_filt %>% 
  left_join(fc, by = c("treatment", "replicate", "time", "season"))
dim(exp.data)
# exp.data %>% 
#   View()

# Phyloseq creation -------------------------------------------------------
# Putting cooler names to the project 
exp.data_miau %>%
  head()
samnames <- exp.data_miau$Sample_ID 
#%>%
  #str_split(pattern = "_",simplify = T) %>%
  # Keep only the first column
  #.[,1]

length(samnames)
length(unique(samnames))
# Golden is the numbers are equal 

rownames(ab.filt_silva) <- samnames
rownames(exp.data_miau) <- samnames

ab.filt_silva_ed <- ab.filt_silva %>%
  rownames_to_column(var = 'selected_file_name') %>%
  mutate(across(c(where(is.character), -selected_file_name), as.numeric)) %>%
  as_tibble()

##SAMPLE NAMES DO NOT MATCH
# Phyloseq
ASV <- otu_table(ab.filt_silva_ed, taxa_are_rows = F)
TAX <- tax_table(as.matrix(tax_silva_miau %>% column_to_rownames(var = 'sequence')))
DAT <- sample_data(exp.data_miau %>% as.data.frame(), sample_names(seq_code)) 

miau.phy.fcs <- phyloseq(ASV,TAX,DAT)

saveRDS(miau.phy.fcs, '../MIAU_seqs/intermediate_files/miau_phyloseq_silva.rds')

