sessionInfo()
##packages
library(tidyverse)
library(magrittr)
library(lubridate)
library(reshape2)
library(plyr)
library(speedyseq)
library(ggplot2)
library(vegan)
library(scales) ##percent format for plot bars
library(microbiome)

setwd ("~/Documentos/Doctorat/REMEI/")

##importing functions
source("src/pseudoabundances_calculation.R") ###for calculating pseudoabundances
source("src/pseudoabundances_calculation_fc.R") 
##we work without Cruise experiment data since we don't have DAPI data for this experiment but for FC we have them.

##DATA
##Línia 796 per fer els càlculs amb les dades de FLOW CITOMETRY
##This phyloseq object has in sample data a column with DAPI counts to calculate pseudoabundances
rem_dapis <- readRDS("./data/intermediate_files/remei_phyloseq_gtbd_DAPI.rds")

#palettes------
palette_large <- c("#a6cee3","#1f78b4", "#009E73","#943132", "#fb9a99", "#FFE867", "#fdbf6f","#ff7f00","#cab2d6",
                   "#6a3d9a", "#d33f6a", "#5562bd",
                   "#c6dec7", "#ef9708", "#c398d9", "#CC6666", "#CABFAB", "#637b88", "#e3a6ce", 
                   "#cee3a6", "#a6b0e3", "#a6e3da", "#e3bba6", "#e0c643", "#84c680", "#e49273", "#004643", 
                   "#80ada0", "#5f5566", 
                   "#773344", "#F2BB05", "#cbf3f0", "#ffbc42", "#cd5334", "#005c69", 
                   "#2c1a1d", "#d496a7", "#a6a15e", "#b0413e", "#1d7874", "#fffffc", 
                   "#247ba0", "#fcca46", "#183059", "#ffa737", "#669bbc", "#bcc4db", "#000000")

palette_treatments_remei <- c("#cab2d6","#6a3d9a", "#a6cee3","#1f78b4", "#009E73", "#F2BB05")

palette_seasons <- c( "#005c69",  "#80ada0", "#ffbc42", "#cd5334", "#2c1a1d")

#this function allows us to create many colors from our palette
palf_large <- colorRampPalette(palette_large) ##to use it scale_color_manual(values=palf_large(x))+ x=number of colours you want

##information from my dataset----
rem_dapis %>% 
  summarize_phyloseq()

##reorder the environmental data
rem_dapis@sam_data$treatment <- rem_dapis@sam_data$treatment %>% 
  factor(levels=(c("CL", "CD", "PL", "PD", "DL", "VL")))
rem_dapis@sam_data$season <- rem_dapis@sam_data$season %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))

##check that all samples have their DAPI counts
rem_dapis@sam_data %>% 
  head()

##extract some information from our dataset----
rem_dapis %>% 
  ntaxa() ##4599 asv
rem_dapis %>% 
  nsamples() ##246 samples

# rem_dapis %>% 
#   sample_names() ##how do our sample names look like?
# rem_dapis %>% 
#   taxa_sums()
# rem_dapis %>% 
#   rank_names()

##explore variables----
# rem_dapis %>% 
#   get_variable(varName = "season")

##information about concrete asv or samples
# rem_dapis %>%
#   get_sample(i="TACGAAGGGGGCTAGCGTTGTTCGGATTTACTGGGCGTAAAGCGCACGTAGGCGGATTGTTAAGTCAGGGGTGAAATCCTGGAGCTCAACT
#              CCAGAACTGCCTTTGATACTGGCAATCTCGAGTCCGGAAGAGGTAAGTGGAACTCCTAGTGTAGAGGTGGAATTCGTAGATATTAGGAAGAACA
#              CCAGTGGCGAAGGCGGCTTACTGGTCCGGAACTGACGCTGAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGT
#              AAACTATGAGAGCTAGCCGTTGGAGGGTTTACCCTTCAGTGGCGCAGCTAACGCATTAAGCTCTCCGCCTGGGGAGTACGGTCGCAAGATTA") 

#abundance values of asv i in all samples
# rem_dapis %>%
#   get_taxa(i="REMEI024")  ##abundance values of all asv in sample i


##subseting by library size----
nsamples(rem_dapis) #246
rem_dapis_filt <- prune_samples(sample_sums(rem_dapis) > 5000, rem_dapis)
nsamples(rem_dapis_filt) #239 hi ha 7 samples per sota dels 5,000

##filter all asv that sum 0-----
# rem_dapis_filt@otu_table %>% 
#   colSums()
rem_dapis_filt <- prune_taxa(taxa_sums(rem_dapis_filt) >0, rem_dapis_filt)
rem_dapis_filt#3736 asv

#we add a column with asv number to plot them easily-----
rem_dapis_filt@tax_table <- rem_dapis_filt@tax_table %>% 
  mutate_tax_table(asv_num = str_c( 'asv' , 1:ncol(rem_dapis_filt@otu_table)))

###FIRST STEP: RELATIVE ABUNDANCES FOR EACH ASV-------
rem_dapis_relabun <- phyloseq::transform_sample_counts(rem_dapis_filt, 
                                                         function(x)
                                                           {x / sum(x)})
# rem_dapis_filt@otu_table %>%
#   dim()
# rem_dapis_relabun@otu_table %>% 
#   dim()
rem_dapis_relabun@tax_table %>% 
  head()
# rem_dapis_relabun@otu_table %>%
#   dim()

##SECOND STEP: CALCULATING PSEUDOABUNDANCES-----
##we melt the data to create a tidy dataset
rem_relabun_melt <- psmelt(rem_dapis_relabun) %>%
  as_tibble
rem_relabun_melt %>% 
  head()

#take a look at the data we have
rem_relabun_melt %T>% 
  colnames() %>% 
  str()

#Pseudoabundances calculation
pseudoabund_df <- rem_relabun_melt %>% 
  pseudoabun(grp = c("asv_num", "treatment", "season", "time", "replicate", "hours_dec"), 
            abundance = Abundance, dapi = dapi)
#he afegit unnest perquè es crea un dataframe dins de les pseudoabund i no em deixa calcular la següen funció

str(pseudoabund_df)

##plot massa gran no es dibuixa amb el meu ordinador es queda penjat------
# pseudoabund_df %>% 
#   #subset(season == "Winter") %>%
#   #subset(treatment == "CL") %>%
#   group_by(season, treatment, time, hours_dec, asv_num) %>% 
#   summarize(mean=mean(pseudoabundance)) %>%
#   ggplot(aes(hours_dec, mean, color = asv_num))+
#   labs(y = "Pseudoabundances mean", x=  "Time (h)")+
#   #scale_y_continuous(labels = percent_format())+
#   geom_point()+
#   geom_smooth(method = 'loess', se = FALSE)+
#   scale_color_manual(values = palf_large(3736))+
#   facet_wrap(~season~treatment)+
#   theme(legend.position = 'none') %>%
#   theme_bw()

##we subset by pseudoabund changes----
##data need to be in wide form
pseudoabund_df_wide <- pseudoabund_df %>%
  pivot_wider(names_from = c("asv_num"), values_from = "pseudoabundance", 
              id_cols = c("treatment", "season", "time", "replicate", "hours_dec"))

#########Now, we were interested in those asvs that showed the largest changes in pseudoabund
#between all times so we identify this particular behavior as follows:
shif<-as.data.frame(pseudoabund_df_wide) #when this function doesn't work is because the  num of samples in env.data != shif samples
shif %>%
  str()

##necessitem afegir el num de mostra per poder transposar i correr el loop 
exp.data.dapis <- rem_dapis_relabun@sam_data
exp.data.dapis$treatment <- exp.data.dapis$treatment %>% 
  factor(levels=(c("CL", "CD", "PL", "PD", "DL", "VL")))
exp.data.dapis$season <- exp.data.dapis$season %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))

exp.data.dapis %>% 
  str()

shif$season <- shif$season %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))

exp.data.dapis$sample_code %>% 
  unique()

shif_m <- shif %>%
  merge(exp.data.dapis, by = c("treatment", "season", "time", "replicate", "hours_dec"))
shif_m %>%
  dim()

shif_filt <- shif_m %>%
  select(starts_with(c('asv', "selected_file_name"))) %>%
  as.data.frame()
str(shif_filt)
shif_filt %>%
  dim()
# colnames(shif_filt)
# dim(shif_filt)
# View(shif_filt)
row.names(shif_filt) <- unique(shif_filt[,3737])
shif_filt <- shif_filt[,-(3737:3738)] #eliminem columnes selected_file_name

#shif_filt_col_names <- row.names(shif_filt)
# shif_filt_ed <- shif_filt_col_names %>% 
#   cbind(shif_filt)
# shif_filt_ed %T>% 
#   dim() %>%
#   View()
# shif_filt_ed <- shif_filt_ed %>%
#   rename(c("." = "selected_file_name"))
# #crec que no fa falta transposar
# shif_filt_t <-  t(shif_filt)
# dim(shif_filt_t)
# View(shif_filt_t)
# dim(shif)
# dim(exp.data.dapis)
shif %>% 
  str()

class <- NULL
for (i in 1:length(names(shif_filt))){
  if (mean(subset(shif_filt[,i],subset= exp.data.dapis$time == "t0")) >mean(subset(shif_filt[,i],subset= exp.data.dapis$time == "t4")))
    class[i]<-"down" #these are the exp.data.dapisting asvs that are more abundant in t0 than in t4
  else if (mean(subset(shif_filt[,i],subset=exp.data.dapis$time == "t0")) >mean(subset(shif_filt[,i],subset=exp.data.dapis$time == "t3")))
    class[i]<-"down" #these are the exp.data.dapisting asvs that are more abundant in t0 than in t3
  else if (mean(subset(shif_filt[,i],subset=exp.data.dapis$time == "t0")) >mean(subset(shif_filt[,i],subset=exp.data.dapis$time == "t2")))
    class[i]<-"down" #these are the exp.data.dapisting asvs that are more abundant in t0 than in t2
  else if (mean(subset(shif_filt[,i],subset=exp.data.dapis$time == "t2")) >mean(subset(shif_filt[,i],subset=exp.data.dapis$time == "t3")))
    class[i]<-"down" #these are the exp.data.dapisting asvs that are more abundant in t2 than in t3
  else if (mean(subset(shif_filt[,i],subset=exp.data.dapis$time == "t3")) >mean(subset(shif_filt[,i],subset=exp.data.dapis$time == "t4")))
    class[i]<-"down" #these are the exp.data.dapisting asvs that are more abundant in t3 than in t4
  else if (mean(subset(shif_filt[,i],subset=exp.data.dapis$time == "t2")) >mean(subset(shif_filt[,i],subset=exp.data.dapis$time == "t4")))
    class[i]<-"down" #these are the exp.data.dapisting asvs that are more abundant in t3 than in t4
  else if (mean(subset(shif_filt[,i],subset=exp.data.dapis$time == "t4")) >mean(subset(shif_filt[,i],subset=exp.data.dapis$time == "t0")))
    class[i]<-"up" #these are the exp.data.dapisting asvs that are more abundant in t4 than in t0
  else if (mean(subset(shif_filt[,i],subset=exp.data.dapis$time == "t4")) >mean(subset(shif_filt[,i],subset=exp.data.dapis$time == "t3")))
    class[i]<-"up" #these are the exp.data.dapisting asvs that are more abundant in t4 than in t3
  else if (mean(subset(shif_filt[,i],subset=exp.data.dapis$time == "t4")) >mean(subset(shif_filt[,i],subset=exp.data.dapis$time == "t2")))
    class[i]<-"up" #these are the exp.data.dapisting asvs that are more abundant in t4 than in t2
  else if (mean(subset(shif_filt[,i],subset=exp.data.dapis$time == "t3")) >mean(subset(shif_filt[,i],subset=exp.data.dapis$time == "t2")))
    class[i]<-"up" #these are the exp.data.dapisting asvs that are more abundant in t3 than in t2
  else if (mean(subset(shif_filt[,i],subset=exp.data.dapis$time == "t2")) >mean(subset(shif_filt[,i],subset=exp.data.dapis$time == "t0")))
    class[i]<-"up" #these are the exp.data.dapisting asvs that are more abundant in t2 than in t0
  else if (mean(subset(shif_filt[,i],subset=exp.data.dapis$time == "t3")) >mean(subset(shif_filt[,i],subset=exp.data.dapis$time == "t0")))
    class[i]<-"up" #these are the exp.data.dapisting asvs that are more abundant in t3 than in t0
  else class[i]<-"Other"}	#these are the rest of the exp.data.dapisting asvs that do not fulfill any of the previous conditions

shifgroup<-data.frame(cbind(names(shif_filt),class))
rownames(shifgroup)<-as.character(shifgroup$V1)
shifgroup$V1<-NULL
shif_category <- shif_filt
colnames(shif_category) <-shifgroup$class

#check if there are some Other asvs in the dataset
shifgroup$class %>%
  subset(shifgroup$class == "other")
##FUNCIONA FINS AQUÍ!
#Separate these two groups of asvs with different behaviors:
down <-shif_filt[,colnames(shif_category) == "down"]##he canviat la coma de lloc crec que potser s'ha de transposar la taula
up <-shif_filt[,colnames(shif_category) == "up"]

dim(down) #3341 asvs  shifters que disminuyen mucho su abundancia al final
dim(up) #asvs que aumentan al final 395

#we need to change row.names so they are the same as in up
sample_names(exp.data.dapis) <-  sample_data(exp.data.dapis)$selected_file_name

up_m <-  merge(up, exp.data.dapis, by.x = 0, by.y= 0) #nomes fan merge els t0
# row.names(up)
# row.names(exp.data.dapis)
# dim(up_m)
# View(exp.data.dapis)
# View(up)
# up_m %>% 
#   colnames()
# up %>% 
#   colnames()
# View(up_m)
# up_m %>%
#   dim()

up_m_melt <- up_m %>% 
  melt(id.vars = c("Row.names", "type", "treatment", "light_regime.x", "replicate", "time", "hours", "season", "fraction", "sample_name",
                   "sample_code", "selected_file_name", "reads", "file_name1", "reads1", "file_name2", "reads2", 
                   "selected_file_name_new", "dapi", "light_regime.y", "hours_dec"))

up_m_melt %T>% 
  head() %>%
  str()
#up_m_melt %>% View()

#plot
##add tax to plots
tax_table <- rem_dapis_filt@tax_table %>%
  mutate_tax_table(tax_ed = (paste(asv_num, class, family, genus, species, sep = "_")))

up_m_melt_tax <- up_m_melt %>% 
  merge(tax_table, by.x = "variable",  by.y = "asv_num" )

up_m_melt$treatment <- up_m_melt$treatment %>% 
  factor(levels=(c("CL", "CD", "PL", "PD", "DL", "VL")))
up_m_melt$season <- up_m_melt$season %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Cruise", "Fall")))

up_m_melt %>%  
  ggplot(aes(hours, value, color = variable))+
  geom_point(alpha=0.8)+
  labs(y="Up asv pseudoabundances")+
  scale_color_manual(values = palf_large(395))+
  scale_y_sqrt()+
  geom_smooth(method = "loess", se=FALSE)+
  facet_grid(~season~treatment)+
  guides(color = guide_legend(ncol = 4))+
  theme_bw()

##subset for testing regressions
# up_m_melt %>%
#   class()
# asv_128_pseudoabund <- up_m_melt %>%
#   subset(season == "Winter") %>%
#   subset(treatment == "PL") %>%
#   subset(variable == "asv128")
# 
# write.csv2(asv_128_pseudoabund,"data/intermediate_files/remei_asv_128_pseudoabund_W_PL.csv", row.names = TRUE)
# asv_128_pseudoabund %>%
#   ggplot(aes(hours, value))+
#            geom_point()+
#            geom_smooth(method = 'loess', se = FALSE)+
#            labs(y = "Up asv pseudoabundances",x=  "Time (h)")+
#            #scale_y_continuous(labels = percent_format())+
#            scale_color_manual(values = palf_large(15))+
#            #facet_grid(~season~treatment)+
#            theme_bw()

##filter by  max increase >1000----
more_than_1000<- up_m_melt %>%
  group_by(variable) %>%
  filter(value >1000, .preserve = TRUE)  %>%
  distinct(variable)
more_than_1000 %>%
  dim()
more_than_1000_all <- up_m_melt_tax %>% 
  inner_join(more_than_1000, by = 'variable')

more_than_1000_all %>%
  ggplot(aes(hours, value, color=tax_ed))+
  geom_point()+
  geom_smooth(method = 'loess', se = FALSE)+
  labs(y = "Up asv pseudoabundances",x=  "Time (h)")+
  #scale_y_continuous(labels = percent_format())+
  scale_color_manual(values = palf_large(15))+
  facet_grid(~season~treatment)+
  theme_bw()

##filter by  max increase >999---- <500----
more_than_500 <-  up_m_melt_tax %>%
  group_by(variable) %>%
  filter(between(value, 500, 999), .preserve = TRUE ) %>%
  distinct(variable)

up_m_melt_tax$value %>% 
  range()
dim(more_than_500)

more_than_500_fil <- more_than_500 %>% 
  anti_join(more_than_1000, by = 'variable')

more_than_500_all <- up_m_melt_tax %>% 
  inner_join(more_than_500_fil, by = 'variable')

more_than_500_fil$variable %>% 
  unique() #we know the number  of colors needed

more_than_500_all %>%
  ggplot(aes(hours, value, color=tax_ed))+
  geom_point()+
  geom_smooth(method = 'loess', se = FALSE)+
  labs(y = "Up asv pseudoabundances",x=  "Time (h)")+
  scale_color_manual(values = palf_large(18))+
  facet_grid(~season~treatment)+
  guides(color = guide_legend(ncol = 2))+
  theme_bw()

##filter by  max increase >499---- <200----
more_than_200 <-  up_m_melt_tax %>%
  group_by(variable) %>%
  filter(between(value, 200, 499), .preserve = TRUE ) %>%
  distinct(variable)

up_m_melt_tax$value %>% 
  range()
dim(more_than_200)

more_than_200_fil <- more_than_200 %>% 
  anti_join(more_than_500, by = 'variable')

more_than_200_all <- up_m_melt_tax %>% 
  inner_join(more_than_200_fil, by = 'variable')

more_than_200_fil$variable %>% 
  unique() #we know the number  of colors needed

more_than_200_all %>%
  ggplot(aes(hours, value, color=tax_ed))+
  geom_point()+
  geom_smooth(method = 'loess', se = FALSE)+
  labs(y = "Up asv pseudoabundances",x=  "Time (h)")+
  scale_color_manual(values = palf_large(48))+
  facet_grid(~season~treatment)+
  guides(color = guide_legend(ncol = 2))+
  theme_bw()

#--------
pseudoabund_df %>% colnames()
pseudoabund_df %>% str()
rem_relabun_melt %>% str()
pseudoabund_df %>% View()
##CALCULATING MEAN AND SD BY REPLICATES
mean.sd <-function(data, grp, Value){
  #if("package:plyr" %in% search()) detach("package:plyr", unload=TRUE) 
  ##unload library plyr (it doesn't work with this library uploaded) because 
  ##plyr::summari[sz]e on a grouped tbl_df seems to destroy the grouping.
  pseudoabund_mean_sd <- data %>%
    group_by_at(grp) %>% 
    summarise_if(is.numeric, mean = mean(Value, na.rm = TRUE), sd = sd(Value, na.rm = TRUE), .groups = "keep") %>%
    ungroup()
  print(head(pseudoabund_mean_sd))
  return(pseudoabund_mean_sd)
}

pseudoabund_df %T>% class() %>% str()
rem_relabun_melt %>% class()

pseudoabund_df <- pseudoabund_df %>% as_tibble


Replicate <- pseudoabund_df %$% 
  Replicate %>% 
  as.character()
pseudoabund_df_test <- pseudoabund_df %>% cbind(Replicate) %>% [,-(5)] 

View(pseudoabund_df)
pseudoabund_df_test <-  pseudoabund_df %>% 
  mutate(Season = as.character(Season), 
         Treatment = as.character(Treatment),
         Replicate = as.character(Replicate),
         Pseudoabundance = as.numeric(Pseudoabundance)) 

str(pseudoabund_df)
#pseudoabund_mean_sd <- 
pseudoabund_mean_sd <- pseudoabund_df %>% 
  mean.sd(grp = c("OTU", "Treatment", "Season", "Time", "Hours_dec"), 
          Value = "Pseudoabundance")

View(pseudoabund_mean_sd)
pseudoabund_df %T>% colnames() %>% groups()
pseudoabund_df %>% groups()

##sense  funció sí que funciona
##amb el summarize sí que m'ho fa per grups llegir article de perquè no hem de calcular amb summarize
pseudoabund_df_test <- pseudoabund_df %>% 
  group_by(OTU, Treatment, Season, Time, Hours_dec) %>%
  summarize(Mean = mean(Pseudoabundance), SD = sd(Pseudoabundance), .groups = "keep")


##provo amb mutate perquè diuen que summarize no calcula bé a veure si funciona
##crec no funciona perquè manté el replicate  per tant no está calculant la  mean
test <- pseudoabund_df %>% 
  group_by(OTU, Treatment, Season, Time, Hours_dec) %>%
  mutate(Mean = mean(Pseudoabundance), SD = sd(Pseudoabundance), .groups = "keep")

test <- pseudoabund_df %>% 
  group_by(OTU, Treatment, Season, Time, Hours_dec) %>%
  aggregate(Mean = mean(Pseudoabundance), SD = sd(Pseudoabundance), .groups = "keep")

View(test)

test %>% tally(sort = TRUE) ##sí que está unint per repliques perque té 3 per cada group
test %>% count(sort = TRUE)
test %>% dim()
test_ed %>% str()

test_ed <- test %>% mutate(Pseudoabundance = as.numeric(Pseudoabundance))
test_ed <- test %>% ungroup()
test_ed <- test %>% as_tibble()

##FUNCIONA
pseudoabund_df_test_mean <- pseudoabund_df_test %>% 
  group_by(OTU, Treatment, Season, Time, Hours_dec) %>%
  summarize(Mean = mean(Pseudoabundance), SD = sd(Pseudoabundance), .groups = "keep")

##intento passar-ho a funció
mean.sd <-function(data, grp, a){
  #if("package:plyr" %in% search()) detach("package:plyr", unload=TRUE) 
  ##unload library plyr (it doesn't work with this library uploaded) because 
  ##plyr::summari[sz]e on a grouped tbl_df seems to destroy the grouping.
  pseudoabund_mean_sd <- data %>%
    group_by_at(grp) %>%  
    summarize(mean = mean(a), .groups = "keep") #, sd = sd(a), 
  print(head(pseudoabund_mean_sd))
  return(pseudoabund_mean_sd)
}

pseudoabund_df_test <- pseudoabund_df_test %>% unnest(cols = "Pseudoabundance")

str(pseudoabund_df_test)
pseudoabund_df_test %>% mean.sd(grp = c("OTU", "Season"), a = "Pseudoabundance")

pseudoabund_df$Pseudoabundance

warnings()

head(pseudoabund_df_test)
View(pseudoabund_df_test)
test_mean %>% dim()
test %>% dim()
View(test_mean)
test %>% colnames()

pseudoabund_df_test_function <- pseudoabund_df_test %>% 
  as_tibble() %>%
  mean.sd(grp = c("OTU", "Treatment", "Season", "Time", "Hours_dec"))

pseudoabund_mean_sd <- pseudoabund_df_test %>% 
  #as_tibble %>%
  #unnest() %>%
  mean.sd(grp = c("OTU", "Treatment", "Season", "Time", "Hours_dec"), 
          Value = Pseudoabundance)

##treballo amb el subset  d'hivern per fer més ràpids els càlculs
pseudoabund_df_w <- pseudoabund_df %>%
  subset(Season == "Winter")

pseudoabund_df_w %>% dim()

group_by(OTU, Treatment, Season, Time, Hours_dec) %>%
  summarize(Mean = mean(Pseudoabundance), SD = sd(Pseudoabundance), .groups = "keep")
  #summarize(Mean=mean(Abundance), SD = sd(Abundance), .groups = "keep") %>%
  #group_by(OTU, Treatment, Season, Time) %>%


##FILTERING BY DEGREE OF PSEUDOABUNDANCES INCREASE 
##TRY TO PIVOT_WIDER
pseudoabund_df_test %>% 
  dim() ##494784
pseudoabund_df %>% 
  dim() ##1267884
pseudoabund_df_test %>% 
  head()
pseudoabund_df_winter <- pseudoabund_df_test %>% subset(Season == "Winter")
View(pseudoabund_df_winter)

pseudoabund_df_winter %>% colnames()
##es queden duplicated OTUs i  per tant no puc fer la resta entre temps...
pseudoabund_wide <- pseudoabund_df_winter %>% 
  pivot_wider(id_cols = c("OTU"), names_from = c("Time", "Treatment"), values_from = "Mean")
pseudoabund_wide %>% 
  head()
View(pseudoabund_wide)
##finde if we have duplicated OTUs as  rownames
pseudoabund_wide %$% OTU %>% duplicated()
dim(pseudoabund_wide)

pseudoabund_df_winter %$% OTU %>% duplicated()

pseudoabund_increase <- pseudoabund_wide %>% 
  mutate(increase_tf = (t4_CL-t0_CL))
View(pseudoabund_increase)

##trobar la manera de filtrar per un gran incrase i recuperar els valors dels altres temps
##calcular els increases per totes les combinacions de temps i mirar les que tenen major increase
##aconseguir  poder treballar amb tot el dataset conjunt no haber de fer tants subsets 
##per poder veure si es comparteixen les asv que creixen més
##aconseguir fer un pivot_wider sense haber de tenir tantes columnes
##tenir una taxonomia més fàcil d'entendre en comptes de les seqüencies de les asv

pseudoabund_increase %>% 
  filter(abs(increase_tf) > 15000) %>%
  melt(values="increase_tf")%>%
  filter(str_detect(variable, "CL")) %>%
  ggplot(aes(variable, value, colour= OTU))+
  geom_point(aes(color=OTU))+
  scale_colour_manual(values= c(palette_large))+
  theme_bw()+
  theme(legend.position = "none")


##growth rate
pseudoabund_gr <- pseudoabund_wide %>% mutate(growth_rate = ((t4-t0)/Hours_dec))

##treballo amb el subset  d'hivern per fer més ràpids els càlculs
rem_dapis_rel_abund_winter@otu_table %>%head()

increase <- rem_dapis_rel_abund_winter %>%
  group_by(OTU, Treatment, Season) %>%
  summarize(Pseudoabundance = (Abundance*DAPI), .groups = "keep") %>%
  summarize(Increase = ())

##testing how to filter the  data set for  the highest increases
rem_dapis_rel_abund_winter %>% filter_taxa()
#filtest<-rowSums(abs(otu_table(comp))>3)>0

rem_dapis_rel_abund_winter %>% 
  sample_sums()

filtest<- rowSums((otu_table(rem_dapis_rel_abund_winter)>0))

dim(filtest)

head(filtest)
View(filtest)
View(tb0)
pseudabund_df %T>% colnames() %>% dim()

##afegir errorbars 
pseudabund_df %>%
  subset(Season == "Winter") %>%
  ggplot(aes(x = Hours_dec, y = Pseudoabundance)) + #
  # geom_bar(stat = "identity", aes(fill="OTU"))+
  facet_wrap(~Treatment~Replicate~Season)+
  geom_point(aes(fill="OTU"))+
  #scale_color_manual(values=paired_genus)+
  #scale_y_continuous(labels = percent_format())+
  #geom_bar(aes(color= "OTU")) +
  labs(x = "Hours", y = "Pseudoabundance") +
  theme_bw()

test %>% colnames()
test %>%
  subset(Season == "Winter") %>%
  ggplot(aes(x = Hours_dec, y = Mean)) + #
  # geom_bar(stat = "identity", aes(fill="OTU"))+
  facet_wrap(~Treatment)+
  geom_point(aes(fill="OTU"))+
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD),
              size=.3,    # Thinner lines
              width=.2, color="grey")+
  #scale_color_manual(values=paired_genus)+
  #scale_y_continuous(labels = percent_format())+
  #geom_bar(aes(color= "OTU")) +
  labs(x = "Hours", y = "Pseudoabundance") +
  theme_bw()

##plot de les means de les pseudoabundances i les sd
tb0 %>%
  #subset(Season == "Winter") %>%
  ggplot(aes(x = Hours_dec, y = Mean)) + #
  # geom_bar(stat = "identity", aes(fill="OTU"))+
  facet_wrap(~Treatment~Season)+
  geom_point(aes(fill="OTU"))+
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD),
                size=.3,    # Thinner lines
                width=.2, color="grey")+
  #scale_color_manual(values=paired_genus)+
  #scale_y_continuous(labels = percent_format())+
  #geom_bar(aes(color= "OTU")) +
  labs(x = "Hours", y = "Pseudoabundance") +
  theme_bw()

##FILTER BY HIGHEST INCREASES TRY ALL COMBINATIONS
##t1-t2

##t2-t3

##t3-t4

##t2-t4

##t1-t3

##t0-t4



###FIFTH PLOTTING PSEUDOABUNDANCES OR GROWTH RATES

ps_winter_phylum_melt %>%
  ggplot(aes(x = Hours, y = Abundance, fill="OTU")) +
  # geom_bar(stat = "identity", aes(fill="OTU"))+
  facet_wrap(~Treatment~Replicate)+
  geom_point(aes(fill="OTU"))+
  scale_color_manual(values=paired_genus)+
  scale_y_continuous(labels = percent_format())+
  #geom_bar(aes(color= "OTU")) +
  labs(x = "Hours", y = "Abundance\n") +
  theme_bw()


##HEAT MAPS
#EXAMPLE
gpt <- subset_taxa(GlobalPatterns, Kingdom=="Bacteria")
gpt <- prune_taxa(names(sort(taxa_sums(gpt),TRUE)[1:300]), gpt)
plot_heatmap(gpt, sample.label="SampleType")


rem_dapis_rel_abund@otu_table %>% head()

##we can change the season and see it more clearly (but not much)
rem_dapis_rel_abund %>% 
  subset_samples( Season=="Winter") %>%
  plot_heatmap()
  
##plot only the most abundant ASV (300)
rem_dapis_rel_abund_winter <- rem_dapis_rel_abund %>% 
  subset_samples( Season=="Winter")
rem_dapis_rel_abund_winter_most <-  
  prune_taxa(names(sort(taxa_sums(rem_dapis_rel_abund_winter),TRUE)[1:300], rem_dapis_rel_abund_winter)) ##NOT WORKING
rem_dapis_rel_abund_winter %>% prune_taxa()
rem_dapis_rel_abund_winter_most  %>% plot_heatmap()

gpt <- prune_taxa(names(sort(taxa_sums(gpt),TRUE)[1:300]), gpt)


##PROVO DE CALCULAR EL LOG2 FOLD CHANGE https://rdrr.io/github/yeguanhuav/visual16S/src/R/log2fc.R 
##no em funciona
log2fc <- function(
  phyloseq,
  feature,
  level = NA,
  p_value = 0.05,
  save_res = FALSE,
  reference = NA,
  treatment = NA
) {
  set.seed(99)
  ## Step 1: Construt table for DESeq2
  # Create a string to parse feature argument to DESeq
  feature_formula <- paste0("~ ", feature)
  # Create a DESeq Data Set
  if (is.na(level)) {
    dds <- phyloseq_to_deseq2(phyloseq, design = as.formula(feature_formula))
    cts <- counts(dds)
    geoMeans <- apply(
      cts, 1, function(row) if (all(row == 0)) 0 else exp(sum(log(row[row != 0]))/length(row))
    )
    dds <- estimateSizeFactors(dds, geoMeans = geoMeans)
  } else {
    countData <- construct_otu_table(phyloseq, level) %>% t()
    colData <- extract_metadata_phyloseq(phyloseq = phyloseq, feature = feature)
    dds <- DESeqDataSetFromMatrix(
      countData = countData, colData = colData, design = as.formula(feature_formula)
    )
  }
  dds <- DESeq(dds)
  ## Step 2: Calculate log2 fold change
  if (!is.na(reference)) {
    res <- results(dds, contrast = c(feature, treatment, reference))
  } else {
    res <- results(dds)
  }
  res <- res[order(res$padj),]
  ## Step 3: Plot fold change
  # Create table for plotting
  ifelse(is.na(level), var <- "OTU", var <- level)
  log2fc <- res %>% as.data.frame() %>% rownames_to_column(var = var) %>%
    filter(!is.na(padj))
  # Identify differential expression
  log2fc$significant = ifelse(
    log2fc$padj > p_value, paste0('p-value > ', p_value), paste0('p-value < ', p_value)
  )
  if (save_res) {
    # save res object to current working directory
    saveRDS(log2fc, "DESeq2_result.rds")
    print('DESeq2_result.rds has been saved to current working directory.')
  }
  # Extract table that meet the p-value threshold to print to the screen
  log2fc_print <- log2fc %>%
    .[.$padj < p_value,] %>%
    select(!!var, log2FoldChange, padj) %>%
    .[order(.[,2], decreasing = TRUE),]
  # EnhancedVolcano plot
  print(res@elementMetadata$description[2])
  print(log2fc_print)
  EnhancedVolcano(
    log2fc, lab = log2fc[[var]], x = "log2FoldChange", y = "padj", FCcutoff = 1, pCutoff = p_value,
    gridlines.major = FALSE, gridlines.minor = FALSE, transcriptPointSize = 3,
    col = c("darkgrey", "#00AFBB", "#FC4E07", "red2"),
    legend = c("Not Significant",
               paste0("p > ", p_value, ", log2 fold change > 1"),
               paste0("p < ", p_value, ", log2 fold change < 1"),
               paste0("p < ", p_value, ", log2 fold change > 1"))
  )
}

rem_dapis_rel_abund_winter %>% log2fc(feature = "Time")
install.packages("DESeq2")

##probo la llibreria Deseq per  calcular els increases
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

library(DESeq2)
pseudoabund_df_test %>% colnames()

rem.phy_w <- rem.phy %>% subset_samples(Season = "Winter")
##
rem_d <- rem.phy %>% phyloseq_to_deseq2(~ Treatment + Time + Replicate ) ##deal with 0 tornar-ho a provar a veure si funciona amb el Replicate
rem_d  <- rem_d %>% estimateSizeFactors(type = "iterate") ##Error in match.arg(method) : 'arg' must be of length 1

rem_dds <- rem_d %>% DESeq()

###PSEUDOABUNDANCES WITH FLOW CITOMETRY ABUNDANCE-----------
##DATA
##This phyloseq object has in sample data a column with FC counts to calculate pseudoabundances
rem_fc <- readRDS("./data/intermediate_files/remei_phyloseq_silva_fc.rds")
rem_fc <- readRDS("./data/intermediate_files/remei_phyloseq_gtbd_fc.rds")

##information from my dataset----
rem_fc %>% 
  summarize_phyloseq()

##reorder the environmental data
rem_fc@sam_data$treatment <- rem_fc@sam_data$treatment %>% 
  factor(levels=(c("CL", "CD", "PL", "PD", "DL", "VL")))
rem_fc@sam_data$season <- rem_fc@sam_data$season %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Early_fall", "Fall")))

##check that all samples have their FC counts
rem_fc@sam_data %>% 
  head()

##extract some information from our dataset----
rem_fc %>% 
  ntaxa() ##5670 asv silva // gtbd 4599
rem_fc %>% 
  nsamples() ##312 samples


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

#we add a column with asv number to plot them easily-----
rem_fc_filt@tax_table <- rem_fc_filt@tax_table %>% 
  mutate_tax_table(asv_num = str_c( 'asv' , 1:ncol(rem_fc_filt@otu_table)))

###FIRST STEP: RELATIVE ABUNDANCES FOR EACH ASV-------
rem_fc_relabun <- phyloseq::transform_sample_counts(rem_fc_filt, 
                                                       function(x)
                                                       {x / sum(x)})

# rem_fc_filt@otu_table %>%
#   dim()
# rem_fc_relabun@otu_table %>% 
#   dim()
rem_fc_relabun@tax_table %>% 
  head()
# rem_fc_relabun@otu_table %>%
#   dim()

##SECOND STEP: CALCULATING PSEUDOABUNDANCES-----
##we melt the data to create a tidy dataset
rem_relabun_melt <- psmelt(rem_fc_relabun) %>%
  as_tibble
#write.table(rem_relabun_melt, "data/rem_relabun_melt.txt", sep= "\t")
rem_relabun_melt %>% 
  head()

#take a look at the data we have
rem_relabun_melt %T>% 
  colnames() %>% 
  str()

#Pseudoabundances calculation
##he hagut de crear una funció nova perquè no entenia que dapi = fc veure com arreglar-ho
pseudoabund_df_fc <- rem_relabun_melt %>% 
  pseudoabun_fc(grp = c("asv_num", "treatment", "season", "time", "replicate", "hours.x"), 
             abundance = Abundance, fc = All)

#he afegit unnest perquè es crea un dataframe dins de les pseudoabund i no em deixa calcular la següen funció

str(pseudoabund_df_fc)

##plot massa gran no es dibuixa amb el meu ordinador es queda penjat------
# pseudoabund_df %>% 
#   #subset(season == "Winter") %>%
#   #subset(treatment == "CL") %>%
#   group_by(season, treatment, time, hours_dec, asv_num) %>% 
#   summarize(mean=mean(pseudoabundance)) %>%
#   ggplot(aes(hours_dec, mean, color = asv_num))+
#   labs(y = "Pseudoabundances mean", x=  "Time (h)")+
#   #scale_y_continuous(labels = percent_format())+
#   geom_point()+
#   geom_smooth(method = 'loess', se = FALSE)+
#   scale_color_manual(values = palf_large(3736))+
#   facet_wrap(~season~treatment)+
#   theme(legend.position = 'none') %>%
#   theme_bw()

##we subset by pseudoabund changes----
##data need to be in wide form
pseudoabund_df_wide_fc <- pseudoabund_df_fc %>%
  pivot_wider(names_from = c("asv_num"), values_from = "pseudoabundance", 
              id_cols = c("treatment", "season", "time", "replicate", "hours.x"))

write.table(pseudoabund_df_wide_fc, "data/intermediate_files/regressions_test_datasets/psueodabund_df_fc_wide_silva.txt", sep = "\t")
write.table(pseudoabund_df_wide_fc, "data/intermediate_files/regressions_test_datasets/psueodabund_df_fc_wide_gtbd.txt", sep = "\t")

##PSEUDOABUNDANCES WITH FLOW CITOMETRY DATA WITH COMPLETE DATASET--------
### aquí tinc les dades de les mostres in situ
##DATA
##This phyloseq object has in sample data a column with FC counts to calculate pseudoabundances
remiau_fc <- readRDS("./data/intermediate_files/remiau_phyloseq_silva_fc.rds")


##information from my dataset----
remiau_fc %>% 
  summarize_phyloseq()

##reorder the environmental data
remiau_fc@sam_data$treatment <- remiau_fc@sam_data$treatment %>% 
  factor(levels=(c("CL", "CD", "PL", "PD", "DL", "VL", "NA")))
remiau_fc@sam_data$season <- remiau_fc@sam_data$season %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Early_fall", "Fall")))

##check that all samples have their FC counts
remiau_fc@sam_data %>% 
  head()

##extract some information from our dataset----
remiau_fc %>% 
  ntaxa() ##14325 asv silva 
remiau_fc %>% 
  nsamples() ##316 samples (té les 4 in situ)

##subseting by library size----
nsamples(remiau_fc) #316
remiau_fc_filt <- prune_samples(sample_sums(remiau_fc) > 5000, remiau_fc)
nsamples(remiau_fc_filt) #312 hi ha 4 samples per sota dels 5,000 silva
#write_rds(remiau_fc_filt, "data/intermediate_files/remiau_fc_filt.rds")
#remiau_fc_filt@sam_data %$% treatment

##filter all asv that sum 0-----
# remiau_fc_filt@otu_table %>% 
#   colSums()
remiau_fc_filt <- prune_taxa(taxa_sums(remiau_fc_filt) >0, remiau_fc_filt)
remiau_fc_filt#5111 asv silva // gtbd (no hi estic treballant)

#we add a column with asv number to plot them easily-----
remiau_fc_filt@tax_table <- remiau_fc_filt@tax_table %>% 
  mutate_tax_table(asv_num = str_c( 'asv' , 1:ncol(remiau_fc_filt@otu_table)))

###FIRST STEP: RELATIVE ABUNDANCES FOR EACH ASV-------
remiau_fc_relabun <- phyloseq::transform_sample_counts(remiau_fc_filt, 
                                                    function(x)
                                                    {x / sum(x)})

#remiau_fc_relabun@sam_data %$% treatment
# remiau_fc_filt@otu_table %>%
#   dim()
# remiau_fc_relabun@otu_table %>% 
#   dim()
remiau_fc_relabun@tax_table %>% 
  head()
# remiau_fc_relabun@otu_table %>%
#   dim()

##SECOND STEP: CALCULATING PSEUDOABUNDANCES-----
##we melt the data to create a tidy dataset
remiau_relabun_melt <- psmelt(remiau_fc_relabun) %>%
  as_tibble()
#write.table(remiau_relabun_melt, "data/remiau_relabun_melt.txt", sep= "\t")
remiau_relabun_melt %>% 
  head()
#remiau_relabun_melt$treatment
#take a look at the data we have
remiau_relabun_melt %T>% 
  colnames() %>% 
  str()

#Pseudoabundances calculation
##he hagut de crear una funció nova perquè no entenia que dapi = fc veure com arreglar-ho
pseudoabund_df_fc <- remiau_relabun_melt %>% 
  pseudoabun_fc(grp = c("asv_num", "treatment", "season", "time", "replicate", "hours.x"), 
                abundance = Abundance, fc = All)

#pseudoabund_df_fc %$% treatment
#he afegit unnest perquè es crea un dataframe dins de les pseudoabund i no em deixa calcular la següen funció

str(pseudoabund_df_fc)

##plot massa gran no es dibuixa amb el meu ordinador es queda penjat------
# pseudoabund_df %>% 
#   #subset(season == "Winter") %>%
#   #subset(treatment == "CL") %>%
#   group_by(season, treatment, time, hours_dec, asv_num) %>% 
#   summarize(mean=mean(pseudoabundance)) %>%
#   ggplot(aes(hours_dec, mean, color = asv_num))+
#   labs(y = "Pseudoabundances mean", x=  "Time (h)")+
#   #scale_y_continuous(labels = percent_format())+
#   geom_point()+
#   geom_smooth(method = 'loess', se = FALSE)+
#   scale_color_manual(values = palf_large(3736))+
#   facet_wrap(~season~treatment)+
#   theme(legend.position = 'none') %>%
#   theme_bw()

##we subset by pseudoabund changes----
##data need to be in wide form
pseudoabund_df_wide_fc <- pseudoabund_df_fc %>%
  pivot_wider(names_from = c("asv_num"), values_from = "pseudoabundance", 
              id_cols = c("treatment", "season", "time", "replicate", "hours.x"))

write.table(pseudoabund_df_wide_fc, "data/intermediate_files/regressions_test_datasets/remiau_psueodabund_df_fc_wide_silva.txt", sep = "\t")
#write.table(pseudoabund_df_wide_fc, "data/intermediate_files/regressions_test_datasets/psueodabund_df_fc_wide_gtbd.txt", sep = "\t")

