

library(tidyverse)
library(speedyseq)
library(janitor) #row_to_names
library(ggpmisc) ##stat_poli_line and equation
library(ggpubr) #as_ggplot function, stat_cor

## web https://rrndb.umms.med.umich.edu/estimate/run_classifier

source("src/pseudoabundances_calculation_fc.R") 

##SCRIPT D'EXEMPLE MARKEL:

# Load tables:
# fpath <- "C:/Users/Usuario/Documents/marta/tax.RDS"
# tax <- readRDS(fpath)

#fpath <- "C:/Users/Usuario/Documents/marta/tax.RDS"
tax <- remei_tax #mirar d'on surt aquesta taula

fpath <- "data/intermediate_files/16S_rna_copy_num_normalization/cnadjusted_remei_asvs.hier.tsv"
cn <- read.table(file = fpath, header = TRUE, sep = "\t", stringsAsFactors = FALSE)


# We need to add $cnumber to tax, using the data that gets as close as possible to the taxonomy of each ASV in tax.

# We need to somehow pair the data from both table considering the taxonomic resolution available for each ASV.
# If phylum, phylum; if genus, genus; etc. So if we have domainA_phylumB, we need to match the cnumber 
# corresponding to domainA_phylumB, not only to domainA.

# Create a modified lineage variable where only the actual taxonomy is left (remove "root", "domain", etc.)
cn$id.lin <- gsub(pattern = "Root;rootrank;|domain;|phylum;|class;|order;|family;|genus;", replacement = "", x = cn$lineage)
# check:
# cn$id.lin[4]

# Remove the last ";" and replace the others with "_" so that the id is "NAME_NAME_NAME" etc.
# Well no there are some taxa with "_" already in their name so better use mmm "." maybe?
cn$id.lin <- gsub(";", ".", gsub(";$", "", cn$id.lin))
# check:
# cn$id.lin[4]

# Similarly, merge the taxonomic assignations in tax the same way:
tax %>%
  colnames()
tax$id.tax <- apply(tax, MARGIN = 1, FUN = function(x){paste(x[c("domain", "phylum", "class", "order", "family", "genus")], collapse = ".")})
# check:
# tax$id.tax[1]

# Remove the ".NA"s introduced by unidentified taxonomy levels:
# check:
# tax$id.tax[4]
# gsub(".NA", "", tax$id.tax[4])
tax$id.tax <- gsub(".NA", "", tax$id.tax)


# Now they have shared IDs.
#
# Add cnumber to taxa:
tax$cnumber <- cn$cnumber[match(tax$id.tax, cn$id.lin)]

#no funciona el match faig aquesta funció
tax <- tax %>%
  left_join(cn, by = c('id.tax' = 'id.lin')) %>%
  mutate(cnumber = remei_asvs.fasta)

# How many ASVs have cnumber?
table(!is.na(tax$cnumber))
# From the top 200:
table(!is.na(tax$cnumber[1:200]))

# And in terms of unique taxa? (to avoid repetitions as ASVs may share the same taxonomy)
table(!is.na(tax$cnumber[!duplicated(tax$id.tax)]))

# So many are missing... An option would be that if there is no match, seek values of previous taxonomic
# levels. E.g., if a genus has no cnumber, see if its family has.

# But first, write down which ASVs have exact matches:
tax$cnumber.exact <- !is.na(tax$cnumber)

cn <- 
 cn %>%
  mutate(cnumber = remei_asvs.fasta)

# Find close results. Now this is more tricky.
# We can maybe do a loop and go entry by entry. If no exact match, crop the id one level, which equals going one up.
# That could work.
# can check with i = 3 (exact match) 6 (genus to family), 46 (genus to order)
for(i in 1:nrow(tax)){
  
  # Is there a cnumber entry for ASV i?
  if(tax$id.tax[i] %in% cn$id.lin){
    
    # If there is, get the cnumber:
    tax$cnumber[i] <- cn$cnumber[cn$id.lin==tax$id.tax[i]]
    
    # Write down match type. In this case, it is an exact match:
    tax$cnumber.match[i] <- "exact"
    
  } else {
    
    # If not, find the closest value.
    # We can do that by cropping the ID one level back and checking again.
    # But never going over "order".
    # To know if we exceed order or not, we can check the id length.
    # Ids look like:
    # domain.phylum.class.order.family.genus
    # So when split by ".", it has to have at least four elements (to be at least order):
    # domain.phylum.class.order
    
    # How long is the ID we are checking?    
    l <- length(strsplit(tax$id.tax[i], split = "\\.")[[1]])
    
    # We can only seek for alternatives if its last level is family or genus, i.e., l>4
    if(l==5){
      
      # If it is family (l=5), we can only try going one level up until reaching order:
      alt <- paste(strsplit(tax$id.tax[i], split = "\\.")[[1]][1:4], collapse = ".")
      
      # Does its order have an entry?
      if(alt %in% cn$id.lin){ 
        
        # if yes:
        tax$cnumber[i] <- cn$cnumber[cn$id.lin==alt]
        tax$cnumber.match[i] <- "from.order"
        
      } else {
        
        # if not, set default value to 1.
        tax$cnumber[i] <- 1
        tax$cnumber.match[i] <- "default.to.1"
        
      }
      
    } else if(l==6){
      
      # If it is genus (l=6), we can try going up twice.
      # First, try with family:
      alt <- paste(strsplit(tax$id.tax[i], split = "\\.")[[1]][1:5], collapse = ".")
      
      # Does its family have an entry?
      if(alt %in% cn$id.lin){ 
        
        # if yes, 
        tax$cnumber[i] <- cn$cnumber[cn$id.lin==alt]
        tax$cnumber.match[i] <- "from.family"
        
      } else {
        
        # if not, try its order:
        alt2 <- paste(strsplit(tax$id.tax[i], split = "\\.")[[1]][1:4], collapse = ".")
        
        if(alt2 %in% cn$id.lin){ 
          
          # if yes:
          tax$cnumber[i] <- cn$cnumber[cn$id.lin==alt2]
          tax$cnumber.match[i] <- "from.order"
          
        } else {
          
          # if not, set default value to 1.
          tax$cnumber[i] <- 1
          tax$cnumber.match[i] <- "default.to.1"
          
        }
      }
    } else {
      
      # If it already was order (or class/phylum) (l<=4), do not seek alternatives.
      # Set instead to 1 as default.
      tax$cnumber[i] <- 1
      tax$cnumber.match[i] <- "default.to.1"
      
    }
  }
}


##add this to find information for genus that have a different taxonomic rute-------
tax$cnumber_match_genus<-cn[match(tax$genus, cn$name),7]
cn_curated<-NULL
for (i in 1:dim(tax)[1]){
  if(is.na(tax[i,12]))
    cn_curated[i]<-tax[i, 9]
  else  cn_curated[i]<-tax[i,12]
}
cn_curated
tax$cn_curated<-cn_curated
cnumber_match_curated<-NULL
for (i in 1:dim(tax)[1]){
  if(is.na(tax[i,12]))
    cnumber_match_curated[i]<-tax[i, 11]
  else  cnumber_match_curated[i]<- "exact"
}
cnumber_match_curated
tax$cnumber_match_curated<-cnumber_match_curated


##in case cnumber match curated !na use this value, otherwise use cnumber
tax %>%
  colnames()

tax <- tax %>%
  mutate(cn_final = case_when(!is.na(cnumber_match_genus) ~ cnumber_match_genus,
                              is.na(cnumber_match_genus) ~ cnumber))

##count nº of ASVs that have 1 as copy number (default)
tax %>%
  filter(cn_final == 1) %$%
  asv_num %>%
  unique() #2099

tax %>%
  filter(cn_final != 1) %$%
  asv_num %>%
  unique() #2495

##Normalization of ASV_tab by 16S rRNA copy number (asv_tab in reads not in relative abundances)-----
##This phyloseq object has in sample data a column with FC counts to calculate pseudoabundances
rem_fc <- readRDS("./data/intermediate_files/remei_phyloseq_silva_fc.rds")
#rem_fc <- readRDS("./data/intermediate_files/remei_phyloseq_gtbd_fc.rds")

##subseting by library size----
nsamples(rem_fc) #312
rem_fc_filt <- prune_samples(sample_sums(rem_fc) > 5000, rem_fc)
nsamples(rem_fc_filt) #308 hi ha 4 samples per sota dels 5,000 silva // amb gtbd n'hi ha 7
#write_rds(rem_fc_filt, "data/intermediate_files/rem_fc_filt.rds")

##filter all asv that sum 0-----
# rem_fc_filt@otu_table %>% 
#   colSums()
rem_fc_filt <- prune_taxa(taxa_sums(rem_fc_filt) >0, rem_fc_filt)
rem_fc_filt#4594 asv silva // gtbd 3810

#we add a column with asv number to plot them easily---
rem_fc_filt@tax_table <- rem_fc_filt@tax_table %>% 
  mutate_tax_table(asv_num = str_c( 'asv' , 1:ncol(rem_fc_filt@otu_table)))

c <- rem_fc_filt@otu_table %>%
  as_tibble(pivot = FALSE) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = 'ASV_seq') %>%
  row_to_names(row_number = 1, remove_row = T, remove_rows_above = F)

##add copy number normalized for each ASV
tax %>%
  colnames()

asv_tab_remei %>%
  colnames()
tax_filt <- tax %>%
  dplyr::select(.otu, cn_final)

asv_tab_cn <- tax %>%
  dplyr::select(.otu, cn_final) %>%
  right_join(asv_tab_remei, by = c('.otu' = '.sample'))

asv_tab_cn %>%
  dim()
asv_tab_cn %>%
  glimpse()

asv_tab_cn_nor <- asv_tab_cn  %>%
  # dplyr::select(-'.otu') %>%
  #mutate(across(where(is.character), as.numeric))
  mutate(across(!'.otu', as.numeric)) %>%
  #test2 <- test %>%
  #select(-cn_final) %>%
  transmute(across(!'.otu' & !'cn_final', funs (./asv_tab_cn$cn_final))) %>% #:310
  cbind(tax_filt)
  #select(-cn_final)
 
##relative abundance for reads normalized by 16S rRNA copy number
asv_tab_cn_nor %>%
  dim()

asv_tab_cn_nor %>%
  apply(1:308, function(x)
{x / sum(x)})

asv_tab_cn_nor_rel <- apply(X = asv_tab_cn_nor[,c(1:308)], MARGIN = 2, function(x) {x / sum(x)}) 

asv_tab_cn_nor_rel %>%
  colnames()

asv_tab_cn_nor_rel %>%
  # as_tibble() %>%
  # dplyr::select(-V309) %>%
  colSums()
asv_tab_cn_nor_rel <- asv_tab_cn_nor_rel %>%
  cbind(asv_tab_cn_nor$.otu)

##Calculate pseudoabundace from flow cytometry data
#add abundance data from flow cytometry to the dataset
sam_data_remei <- rem_fc_filt@sam_data %>%
  as_tibble()
##we melt the data to create a tidy dataset
asv_tab_cn_nor_rel_l <- asv_tab_cn_nor_rel %>%
  as_tibble() %>%
  pivot_longer(cols = starts_with(c('R', '66')), names_to = 'sample_id', values_to = 'abundance') %>%
  mutate(sample_id_ed = str_replace(sample_id, '_/', '')) %>%
  left_join(sam_data_remei, by = c('sample_id_ed' = '.sample')) %>%
  mutate(pseudoabundance = as.numeric(abundance)*as.numeric(All))

#  asv_tab_cn_nor_rel_l_filt <- asv_tab_cn_nor_rel_l %>% ##for comparison with CARD-FISH data
#    left_join(tax, by = c('V309' = '.otu')) %>%
#    select(asv_num, treatment, season, replicate, abundance, time, domain, phylum, class, order, family, genus)
# # 
#   write.table(asv_tab_cn_nor_rel_l_filt, 'data/intermediate_files/asv_tab_cn_nor_rel_l_filt.txt', sep = '\t')

asv_tab_cn_nor_rel_l %>%
  colnames()
tax %>%
  colnames()

#pivot wide and add asv_num
pseudoabund_cn_nor_w <- asv_tab_cn_nor_rel_l %>%
  select(V309, sample_id_ed, treatment, replicate, time, hours.x, season, pseudoabundance) %>%
  left_join(tax, by = c('V309' = '.otu')) %>%
  select(asv_num, treatment, replicate, time, hours.x, season, pseudoabundance) %>%
  pivot_wider(id_cols = c('treatment', 'replicate', 'time', 'hours.x', 'season'), 
              names_from = 'asv_num', 
              values_from = 'pseudoabundance')

#Recalculate growth rates with 16S rRNA copy nº normalization---------
#functions needed
source("src/filter.season.treatment.time.R") ###filtrar les pseudoabund per treatment season i time per calcular regressions
source("src/multiple.linear.regressions.R") ##per fer multiples regressions d'un dataframe
source("src/comparing.reg.R") ##per comparar les taxes de creiement entre 4 i 3 punts (les regressions)

#objecte de phyloseq del que traiem la metadata i la tax
#rem_fc_filt_GTBD <- readRDS("data/intermediate_files/remei_phyloseq_gtbd_fc.rds") 
#rem_fc_silva <- readRDS("data/intermediate_files/remei_phyloseq_silva_fc.rds") 
env <- rem_fc_filt@sam_data
#nsamples(rem_fc_filt_SILVA) #312 (sense filtrar, pot ser un problema amb la taxonomia i el número d'ASV)
#substituint silva per GTBD o al revés fas els càlculs d'un o l'altre.
# #asv's taxonomy this new column is important for ploting asv with their taxonomic information
# nsamples(rem_fc_silva) #312
# rem_fc_filt <- prune_samples(sample_sums(rem_fc_silva) > 5000, rem_fc_silva)
# nsamples(rem_fc_filt) #308 hi ha 4 samples per sota dels 5,000 silva // amb gtbd n'hi ha 7
# #write_rds(rem_fc_filt, "data/intermediate_files/rem_fc_filt.rds")
# 
# ##filter all asv that sum 0-----
# # rem_fc_filt@otu_table %>% 
# #   colSums()
# rem_fc_filt <- prune_taxa(taxa_sums(rem_fc_filt) >0, rem_fc_filt)
# rem_fc_filt#4594 asv silva // gtbd 3810
# 
# rem_fc_filt@tax_table <- rem_fc_filt@tax_table %>% 
#   mutate_tax_table(asv_num = str_c( 'asv' , 1:ncol(rem_fc_filt@otu_table)))
# tax_table <- rem_fc_filt@tax_table %>%
#   mutate_tax_table(tax_ed = (paste(asv_num, class, family, genus, species, sep = "_")))
# 
# #metadata 
# env <- rem_fc_filt@sam_data
# env$treatment <- env$treatment %>% 
#   factor(levels=(c("CL", "CD", "PL", "PD", "DL", "VL")))
# env$season <- env$season %>% 
#   factor(levels=(c("Winter", "Spring", "Summer", "Early_fall", "Fall")))
# 
# #asv_tab & metadata
# pseudoabund_df_wide_m <-pseudoabund_cn_nor_w %>%
#   merge(env, by = c("treatment", "season", "time", "replicate", "hours.x"))

#1. Prepare  datasets for calculating regressions------
###crear una nova variable que combina time and replicate per no tenir problemes amb els temps repetits de les rèpliques.
data <- pseudoabund_cn_nor_w %>% #pseudoabundances dataset form fc data in wide format with metadata.
  mutate(t_rep = paste(time, replicate, sep="_"))

env <- env %>%
  as_tibble() %>%
  mutate(t_rep = paste(time, replicate, sep="_"))

data %>%
  glimpse()

## 2. Filtering for each treatment, season and time
##intento crear una funció

asv_filt_t0234 <- filter.season.treatment.time(data = data,  
                                              treatment_col = treatment, treatment_keep = "PD", 
                                              season_col = season, season_keep = "Winter",
                                              time_col = time, time_keep = c("t0", "t2", "t3", 't4')) %>%
  select(starts_with("asv"), #we clean metadata columns that we don't need now and keep only the column t_rep
         matches(c("t_rep", "hours.x"))) %>%
  mutate(across(everything(), ~replace(., . == 0, "1"))) %>% ##0 values will be substituted by 1 because it is mandatory for transforming pseudoabundances to log
  mutate(across(!t_rep, as.numeric))
# 
#  asv_filt_t023 %>%
#     glimpse() #las dos ultimas columnas són t_rep i hours.x 

## 3. Calculating regressions (CHANGE MATRIX DIM FOR gtbd: 3810, SILVA: 4594)
regress <- apply(X = asv_filt_t0234[,c(1:4594)], MARGIN = 2, FUN = multiple.linear.regressions, env = asv_filt_t0234$hours.x) ##DIMENSIONS OF THE MATRIX SHOULD BE ADAPTED TO DATA!!!!!
regress <- as.data.frame(t(regress)) %>% # le damos la vuelta y convertimos en data frame
  rownames_to_column(var = "asv_num")

# resultado %T>% #si surten NAs en el resultat és perquè les dimensions a funció d'apply están malament.
#   dim() %>%
#   head()

#cambiar nombre  del objeto en función del tratamiento y la estació
##vull un dataset amb tota la informació important el problema és que em quedo sense memòria de l'ordinador haig de simCLificar el procés
asv_m_reg_Winter_PD_t0234 <- filter.season.treatment.time(data = env,  
                                                         treatment_col = treatment, treatment_keep = "PD", 
                                                         season_col = season, season_keep = "Winter",
                                                         time_col = time, time_keep = c("t0", "t2", "t3", "t4")) %>%
  select(starts_with("asv"), ".sample", "treatment", "replicate",  "time", "season", "sample_name", 
         "sample_code", "selected_file_name",
         "reads", "light_regime.y", "hours.y", "LNA", "HNA", "All", "BM", "Leu.PB", "SGR", "TD", "t_rep")  %>%
  left_join(asv_filt_t023, by = "t_rep") %>%
  pivot_longer(cols = starts_with("asv"), names_to = "asv_num", values_to = "pseudoabundance") %>%
  left_join(regress, by = "asv_num") %>%
  add_column(regression_times = "t0234")


##tornar al punt 2 substituir per les noves condicions amb les que volem calcular les regressions.

##dataset proba per provar de comparar regressions
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
# 
# reg_early_fall <- bind_rows(asv_m_reg_Early_fall_CD_t0234, 
#                             asv_m_reg_Early_fall_CL_t0234,
#                             asv_m_reg_Early_fall_PL_t0234,
#                             asv_m_reg_Early_fall_PD_t0234, 
#                             asv_m_reg_Early_fall_DL_t0234,
#                             asv_m_reg_Early_fall_VL_t0234)

reg_fall <- bind_rows(asv_m_reg_Fall_CD_t0234, 
                      asv_m_reg_Fall_CL_t0234,
                      asv_m_reg_Fall_PL_t0234,
                      asv_m_reg_Fall_PD_t0234, 
                      asv_m_reg_Fall_DL_t0234,
                      asv_m_reg_Fall_VL_t0234)

##si vull fer totes les regressions amb diferents temps el que hauria de fer és afegir una nova columna que digui 
## de quin temps a quin temps es correspon cada regressió.

##unir tots els datasets
reg_all_t0234_silva <- bind_rows(reg_winter, 
                                 reg_spring,
                                 reg_summer,
                                 #reg_early_fall,
                                 reg_fall)

#regressions 4 punts
write.csv(reg_all_t0234_silva, "data/intermediate_files/asv_reg_cn_nor_t0234_silva.csv")
#regressions 3 punts
write.table(reg_all_t023_silva, "data/intermediate_files/asv_reg_cn_nor_t023_silva.txt", sep = "\t")

###Regressions t3 and t4
reg_all_t0234_silva <-  read.table("data/intermediate_files/asv_reg_cn_nor_t0234_silva.csv", sep="\t", header = T) %>%
  as_tibble()
reg_all_t023_silva <-  read.table("data/intermediate_files/asv_reg_cn_nor_t023_silva.txt", sep="\t", header = T) %>%
  as_tibble()

#function needed
source("src/comparing.reg.R") ##per comparar les taxes de creiement entre 4 i 3 punts

#primer eliminem duplicats perquè tinguin les mateixes files els dos datasets
# df1 <- reg_all_t023_silva %>% 
#   distinct(treatment, season, asv_num, slope, slope.pval)
# df2 <- reg_all_t0234_silva %>% 
#   distinct(treatment, season, asv_num, slope, slope.pval)
# comp_df <- df1 %>% ##funció per comparar les taxes de creixement (slopes) entre 3 i 4 temps.
#   left_join(df2, by = c("treatment", "season", "asv_num")) %>%
#   mutate(slopes_compared = case_when(slope.x > slope.y & slope.pval.x < 0.05 ~ slope.x,
#                                      slope.x < slope.y ~ slope.y),
#          pvalue_slope_kept =case_when(slope.x > slope.y & slope.pval.x < 0.05 ~ slope.pval.x,
#                                       slope.x < slope.y ~ slope.pval.y))
# comp_df_tax <- comp_df %>% 
#   left_join(tax_table, by= "asv_num", copy = TRUE) #afegim la tax pels plots

reg_all_t023_silva %>%
  dim()
reg_all_t0234_silva %>%
  dim()

#df1 ha de ser 3 temps i df2 4 temps pel problema del PL spring (sol·lucionar-ho)
reg_all_slopes_chosen_silva_cn_nor <- comparing.reg(df1 = reg_all_t023_silva, 
                                             df2 = reg_all_t0234_silva,
                                             treatment = treatment,
                                             season = season,
                                             asv_num = asv_num, 
                                             slope = slope, 
                                             slope.pval = slope.pval)

#Afegim taxonomia i fem boxplots generals per tenir un overview de la coherència que hi ha de les ----
##taxes de creiement per grups taxonomics
tax_red <- tax %>%
  select(domain, phylum, class, order, family, genus, species, asv_num)

reg_all_slopes_chosen_silva_tax_cn_nor <- reg_all_slopes_chosen_silva_cn_nor %>%
  left_join(tax_red, by = "asv_num", copy = TRUE) 

write.csv(reg_all_slopes_chosen_silva_tax_cn_nor, "data/intermediate_files/reg_all_slopes_chosen_silva_tax_cn_nor.csv")

###Correlation with growth rates no-normalized by 16S rRNA copy genes
reg_all_slopes_chosen_silva_tax <- read.csv("data/intermediate_files/reg_all_slopes_chosen_silva_tax_corrected_asv_num.csv", sep=",") %>%
  filter(season != "Early_fall") ##la taxonomia és incrrecte perque l'asv_num no estava filtrat anteriorment i no es corresponia.

reg_all_slopes_chosen_silva_tax %>%
  colnames()
##prepare datasets
reg_all_slopes_chosen_silva_tax_filt <- reg_all_slopes_chosen_silva_tax %>%
  dplyr::filter(pvalue_slope_chosen < 0.05 &
                  slope_chosen > 0) %>%
  mutate(slope_chosen_days = slope_chosen*24) ##transform slopes /hours to /days

reg_all_slopes_chosen_silva_tax_cn_nor_filt <- reg_all_slopes_chosen_silva_tax_cn_nor %>%
  dplyr::filter(pvalue_slope_chosen < 0.05 &
                   slope_chosen > 0) %>%
  mutate(slope_chosen_days = slope_chosen*24) ##transform slopes /hours to /days

reg_all_slopes_chosen_silva_tax_filt %>%
  dim()

reg_all_slopes_chosen_silva_tax_cn_nor_filt %>%
  dim()

reg_all_slopes_chosen_silva_tax_cn_nor_non <- reg_all_slopes_chosen_silva_tax_cn_nor_filt %>%
  left_join(reg_all_slopes_chosen_silva_tax_filt, by = c('treatment', 'season', 'asv_num', 'phylum', 'class', 'order', 
                                                        'family', 'species', 'genus'), suffix = c('.cn_nor', '.non')) 
##before correlation we check normality
shapiro.test(as.numeric(reg_all_slopes_chosen_silva_tax_cn_nor_non$slope_chosen_days.cn_nor)) # =>p-value < 2.2e-16 (NO NORMALITY)
ggqqplot(as.numeric(reg_all_slopes_chosen_silva_tax_cn_nor_non$slope_chosen_days.cn_nor))

shapiro.test(as.numeric(reg_all_slopes_chosen_silva_tax_cn_nor_non$slope_chosen_days.non)) # => p-value < 2.2e-16 (NO NORMALITY)
ggqqplot(as.numeric(reg_all_slopes_chosen_silva_tax_cn_nor_non$slope_chosen_days.non))

corr_cn_nor_non <- reg_all_slopes_chosen_silva_tax_cn_nor_non %>%
  ggplot(aes(slope_chosen_days.cn_nor, slope_chosen_days.non))+
  geom_point(aes(shape = treatment, color = season), alpha = 0.6, size = 0.5)+
  scale_fill_manual(values = palette_seasons_4)+
  scale_y_continuous(limits = c(0, 10), expand = c(0, 0))+
  scale_x_continuous(limits = c(0, 10), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  # stat_poly_eq(aes(
  #   label =  paste(after_stat(p.value.label))), 
  #   #formula = x ~ y,
  #   method = stats::lm,
  #   p.digits = 2, 
  #   coef.keep.zeros = T)+
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed')+
  #stat_poly_eq(color = 'black')+
  stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, label.y = 9.5, p.digits = 0.01, digits = 2, p.accuracy = 0.01, method = 'spearman')+
  guides(size = 'none', shape = guide_legend(override.aes = list(size = 1.5)),
         color= guide_legend(override.aes = list(size = 1.5)))+
  labs(y = 'ASV-based growth rates', 
       x = 'ASV-based growth rates normalized\nby 16S rRNA copy number', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/3), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        text = element_text(size = 7))
  
ggsave('correlation_gr_cn_nor_non_ed.pdf', corr_cn_nor_non, 
       path = "~/Documentos/Doctorat/REMEI/results/figures/corrected_asv_num/",
       width = 88,
       height = 88,
       units = 'mm')

##Violin plots----
reg_all_slopes_chosen_silva_tax_cn_nor_filt$treatment <- reg_all_slopes_chosen_silva_tax_cn_nor_filt$treatment %>% 
  factor(levels=(c("CL", "CD", "PL", "PD", "DL", "VL")))
reg_all_slopes_chosen_silva_tax_cn_nor_filt$season <- reg_all_slopes_chosen_silva_tax_cn_nor_filt$season %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))

seas <- reg_all_slopes_chosen_silva_tax_cn_nor_filt %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  filter(season != "Early_fall") %>%
  #filter(pseudoabundance > 1000) %>%
  ggplot(aes(season, slope_chosen_days))+ #
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
  scale_y_continuous(expand = c(0, 0), limits = c(0,10.4))+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_treatments_remei)+
  #facet_wrap(.~phylum, scales = "free")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

treat <- reg_all_slopes_chosen_silva_tax_cn_nor_filt %>%
  # filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  # filter(slope_chosen_days > 0) %>%
  # filter(season != "Early_fall") %>%
  #filter(pseudoabundance > 1000) %>%
  ggplot(aes(treatment, slope_chosen_days))+ #
  geom_point(aes(treatment, slope_chosen_days, color = season), alpha = 0.7, position = position_jitter(0.25))+ #position_jitter(0.1)  #, position = "identity"
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
  scale_y_continuous(expand = c(0, 0), limits = c(0,10.4))+ #forçar que comenci al 0 l'eix.
  #geom_smooth(method = "lm")
  scale_color_manual(values = palette_seasons_4)+
  #facet_wrap(.~phylum, scales = "free")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

##number of GR calculated
reg_all_slopes_chosen_silva_tax_cn_nor_filt %>%
  group_by(slope_chosen_days, asv_num) %>%
  summarize(n = n()) %>%
  ungroup()%>%
  summarize(n = sum(n)) #3847 growth rates 

##PLOT for only those ASVs that represent >1% at some point of the community-----
abundant_asv_cn_nor <- asv_tab_cn_nor_rel_l_filt %>% 
  filter(abundance > 0.01) %>% #more than 1% of the community at some point
  dplyr::select(asv_num) %>%
  unique() %>%
  as_tibble()

#filter by abundant_asv----
reg_all_slopes_chosen_silva_tax_filt_cn_1perc <- reg_all_slopes_chosen_silva_tax_cn_nor_filt %>% 
  filter(asv_num %in% abundant_asv_cn_nor$asv_num)

seas_1perc_cn_nor <- reg_all_slopes_chosen_silva_tax_filt_cn_1perc %>%
  #filter(pseudoabundance > 1000) %>%
  ggplot(aes(season, slope_chosen_days))+ #
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
  scale_y_continuous(expand = c(0, 0), limits = c(0,10.4))+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_treatments_remei)+
  #facet_wrap(.~phylum, scales = "free")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

treat_1perc_cn_nor <- reg_all_slopes_chosen_silva_tax_filt_cn_1perc %>%
# filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
# filter(slope_chosen_days > 0) %>%
# filter(season != "Early_fall") %>%
#filter(pseudoabundance > 1000) %>%
ggplot(aes(treatment, slope_chosen_days))+ #
  geom_point(aes(treatment, slope_chosen_days, color = season), alpha = 0.7, position = position_jitter(0.25))+ #position_jitter(0.1)  #, position = "identity"
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
  scale_y_continuous(expand = c(0, 0), limits = c(0,10.4))+ #forçar que comenci al 0 l'eix.
  #geom_smooth(method = "lm")
  scale_color_manual(values = palette_seasons_4)+
  #facet_wrap(.~phylum, scales = "free")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

gr_1perc_cn_nor <- ggarrange(seas_1perc_cn_nor, treat_1perc_cn_nor)

ggsave('gr_1perc_cn_nor.pdf', gr_1perc_cn_nor, 
       path = "~/Documentos/Doctorat/REMEI/results/figures/corrected_asv_num/",
       width = 188,
       height = 88,
       units = 'mm')
