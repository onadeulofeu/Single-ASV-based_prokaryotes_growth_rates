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
library(janitor) #row_to_names
library(speedyseq)

setwd ("~/Documentos/Doctorat/REMEI/")
##remmiau dataset---------
remiau <- readRDS("./data/intermediate_files/remiau_phyloseq_silva_fc.rds")

##Question: Are ASV growing at highest growth rates the ones that were dominating 
##the natural community (in situ samples)?
###GO TO LINE 991 FOR THE FIGURES
##Which % of the community is growing actively? (FALTA CONTESTAR)

##information from my dataset----
remiau %>% 
  summarize_phyloseq()

##reorder the environmental data
remiau@sam_data$treatment <- remiau@sam_data$treatment %>% 
  factor(levels=(c("CL", "CD", "PL", "PD", "DL", "VL", 'NA')))
remiau@sam_data$season <- remiau@sam_data$season %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))

##check that all samples have their DAPI counts
remiau@sam_data %>% 
  head()

##extract some information from our dataset----
remiau %>% 
  ntaxa() ##14325 asv

remiau %>% 
  nsamples() ##316 samples

# remiau %>% 
#   sample_names() ##how do our sample names look like?
# remiau %>% 
#   taxa_sums()
# remiau %>% 
#   rank_names()

##explore variables----
# remiau %>% 
#   get_variable(varName = "season")

##information about concrete asv or samples
# remiau %>%
#   get_sample(i="TACGAAGGGGGCTAGCGTTGTTCGGATTTACTGGGCGTAAAGCGCACGTAGGCGGATTGTTAAGTCAGGGGTGAAATCCTGGAGCTCAACT
#              CCAGAACTGCCTTTGATACTGGCAATCTCGAGTCCGGAAGAGGTAAGTGGAACTCCTAGTGTAGAGGTGGAATTCGTAGATATTAGGAAGAACA
#              CCAGTGGCGAAGGCGGCTTACTGGTCCGGAACTGACGCTGAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGT
#              AAACTATGAGAGCTAGCCGTTGGAGGGTTTACCCTTCAGTGGCGCAGCTAACGCATTAAGCTCTCCGCCTGGGGAGTACGGTCGCAAGATTA") 

#abundance values of asv i in all samples
# remiau %>%
#   get_taxa(i="REMEI024")  ##abundance values of all asv in sample i


##subseting by library size----
nsamples(remiau) #316
remiau_filt <- prune_samples(sample_sums(remiau) > 5000, remiau)
nsamples(remiau_filt) #312 hi ha 4 samples per sota dels 5,000

##filter all asv that sum 0-----
# remiau_filt@otu_table %>% 
#   colSums()
remiau_filt <- prune_taxa(taxa_sums(remiau_filt) >0, remiau_filt)
remiau_filt#5111 asv

#we add a column with asv number to plot them easily-----
remiau_filt@tax_table <- remiau_filt@tax_table %>% 
  mutate_tax_table(asv_num = str_c( 'asv' , 1:ncol(remiau_filt@otu_table)))

###FIRST STEP: RELATIVE ABUNDANCES FOR EACH ASV-------
remiau_relabun <- phyloseq::transform_sample_counts(remiau_filt, 
                                                       function(x)
                                                       {x / sum(x)})
# remiau_filt@otu_table %>%
#   dim()
# remiau_relabun@otu_table %>% 
#   dim()
remiau_relabun@tax_table %>% 
  head()
# remiau_relabun@otu_table %>%
#   dim()

remiau_relabun@sam_data %>%
  head()
metadata <- remiau_relabun@sam_data
### No puc calcular pseudoabundances d'aquestes dades perquè no tinc les abundàncies dels insitu data
insitu_phy_obj  <- remiau_relabun %>%
  subset_samples(selected_file_name %in% c("x254", "x255", "x256", "x257"))

insitu_phy_obj_melt <- psmelt(insitu_phy_obj) %>%
  as_tibble

plot_bar(insitu_phy_obj, fill="class")

remiau_relabun %>%
  subset_samples(season %in% c("Winter")) %>%
  plot_bar(x= 'time', fill = 'phylum', facet_grid = 'treatment')+
  scale_color_manual(values = new_palette_large)

remiau_relabun %>%
  class()
metadata%>%
  colnames()

remiau_relabun %>%
  group_by() %>% 
  summarise_if(is.numeric, mean = mean(Value, na.rm = TRUE), sd = sd(Value, na.rm = TRUE), .groups = "keep") %>%
  ungroup()

##from the phyloseq object create the three datasets (metadata, asv_tab & taxonomy)
instiu_metadata <- insitu_phy_obj@sam_data
insitu_asv_tab <- insitu_phy_obj@otu_table
insitu_tax <- insitu_phy_obj@tax_table

insitu_phy_obj_melt %>%
  colnames()


##ordenem els factors
insitu_phy_obj_melt <- insitu_phy_obj_melt %>%
  mutate(domain_f = as_factor(domain),
         phylum_f = as_factor(phylum),
         class_f = as_factor(class),
         order_f = as_factor(order),
         family_f = as_factor(family),
         genus_f = as_factor(genus),
         asv_f = as_factor(asv_num))

insitu_phy_obj_melt %>%
  colnames()

insitu_phy_obj_melt$phylum_f <-  factor(insitu_phy_obj_melt$phylum_f, 
                                                          levels=unique(
                                                            insitu_phy_obj_melt$phylum_f[
                                                              order(insitu_phy_obj_melt$domain_f)]), 
                                                          ordered=TRUE)

insitu_phy_obj_melt$class_f <-  factor(insitu_phy_obj_melt$class_f, 
                                                         levels=unique(
                                                           insitu_phy_obj_melt$class_f[
                                                             order(insitu_phy_obj_melt$domain_f,
                                                                   insitu_phy_obj_melt$phylum_f)]), 
                                                         ordered=TRUE)
insitu_phy_obj_melt$order_f <-  factor(insitu_phy_obj_melt$order_f, 
                                                         levels=unique(
                                                           insitu_phy_obj_melt$order_f[
                                                             order(insitu_phy_obj_melt$domain_f,
                                                                   insitu_phy_obj_melt$phylum_f,
                                                                   insitu_phy_obj_melt$class_f)]), 
                                                         ordered=TRUE)

insitu_phy_obj_melt$family_f <-  factor(insitu_phy_obj_melt$family_f, 
                                                          levels=unique(
                                                            insitu_phy_obj_melt$family_f[
                                                              order(insitu_phy_obj_melt$domain_f,
                                                                    insitu_phy_obj_melt$phylum_f,
                                                                    insitu_phy_obj_melt$class_f,
                                                                    insitu_phy_obj_melt$order_f)]), 
                                                          ordered=TRUE)

insitu_phy_obj_melt$genus_f <-  factor(insitu_phy_obj_melt$genus_f, 
                                                         levels=unique(
                                                           insitu_phy_obj_melt$genus_f[
                                                             order(insitu_phy_obj_melt$domain_f,
                                                                   insitu_phy_obj_melt$phylum_f,
                                                                   insitu_phy_obj_melt$class_f,
                                                                   insitu_phy_obj_melt$order_f,
                                                                   insitu_phy_obj_melt$family_f)]), 
                                                         ordered=TRUE)
insitu_phy_obj_melt$asv_f <-  factor(insitu_phy_obj_melt$asv_f, 
                                                       levels=unique(
                                                         insitu_phy_obj_melt$asv_f[
                                                           order(insitu_phy_obj_melt$domain_f,
                                                                 insitu_phy_obj_melt$phylum_f,
                                                                 insitu_phy_obj_melt$class_f,
                                                                 insitu_phy_obj_melt$order_f,
                                                                 insitu_phy_obj_melt$family_f,
                                                                 insitu_phy_obj_melt$genus_f)]), 
                                                       ordered=TRUE)

##gràfic de barres dels temps in situ
insitu_phy_obj_melt %>%
  ggplot(aes(season, Abundance, fill = fct_infreq(phylum)))+
  geom_col(aes())+
  scale_fill_manual(values = palette_phylums)+
  scale_y_continuous(labels = percent_format())+
  #scale_y_continuous(percent_format())+
  theme_bw()

insitu_phy_obj_melt %>%
  colnames()

insitu_phy_obj_melt %$%
  phylum %>%
  unique()

unique_cl_phy <- insitu_phy_obj_melt %>%
  distinct(class, phylum)

insitu_phy_obj_melt$class_f <- insitu_phy_obj_melt %>%
  interaction(factor(class_f), phylum, sep = '\n')

insitu_phy_obj_melt$class_f <- insitu_phy_obj_melt$class_f %>%
  factor(levels = c('Gammaproteobacteria', 'Alphaproteobacteria', 'Zetaproteobacteria',  'Bacteroidia', 'Rhodothermia',  
                    'Ignavibacteria', 'Chlorobia', 'Kryptonia', 'Nitrososphaeria',   'Cyanobacteriia', 
                    'Vampirivibrionia',   'Acidimicrobiia',   'Actinobacteria', 'Coriobacteriia',  'Thermoleophilia',
                    'Verrucomicrobiae',  'Kiritimatiellae', 'Lentisphaeria', 'Omnitrophia', 'Chlamydiae', 
                    'Planctomycetes',  'Phycisphaerae',   'Blastocatellia',   'Holophagae', 'Vicinamibacteria', 
                    'Acidobacteriae',  'Thermoanaerobaculia',  'Bdellovibrionia',  'Oligoflexia',  'Bacilli', 
                    'Clostridia', 'Negativicutes', 'Syntrophomonadia', 'Desulfitobacteriia', 'Thermoanaerobacteria', 
                    'Myxococcia', 'Polyangia',  'Nitrospinia',  'Campylobacteria',   'Deinococci', 
                    'Fusobacteriia', 'Desulfovibrionia', 'Desulfobulbia', 'Desulfuromonadia',  'Desulfobacteria',  
                    'Methylomirabilia', 'Gemmatimonadetes', 'Anaerolineae', 'Chloroflexia', 'Spirochaetia', 
                    'Leptospirae', 'Calditrichia', 'Halobacteria',  'Nitrospiria', 'Babeliae', 
                    'Saccharimonadia', 'Cloacimonadia', 'Synergistia', 'Abditibacteria',  'Deferribacteres'), 
         ordered = TRUE)

insitu_phy_obj_melt %>%
  

insitu_community <- insitu_phy_obj_melt %>%
  group_by(season, class_f, time) %>%
  summarize(sum_abundance = sum(Abundance)) %>%
  ggplot(aes(time, sum_abundance, fill = class_f))+
  geom_col(aes(fill = fct_rev(class_f)))+
  #scale_fill_manual(values = palette_phylums_assigned)+
  scale_fill_manual(values = palette_class_assigned)+
  # scale_fill_manual(values = palette_class_assigned, 
  # labels = interaction( reg_all_slopes_chosen_silva_tax_filt_rel_gr$phylum_f,
  #                       reg_all_slopes_chosen_silva_tax_filt_rel_gr$class_f, sep = '\n'))+
  scale_y_continuous(labels = percent_format())+
  facet_grid(season~.)+
  labs(x = 'In situ community\n composition')+
  #scale_y_continuous(percent_format())+
  theme_minimal()+
  theme(axis.text.x  = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.text = element_blank(),
        legend.position = 'none')

##ploting max growht rates in relative abundances
reg_all_slopes_chosen_silva_tax %>%
  colnames()

reg_all_slopes_chosen_silva_tax$season <-  reg_all_slopes_chosen_silva_tax$season %>% factor( 
  levels = (c('Winter', 'Spring', 'Summer', 'Fall', 'Early_fall')))
reg_all_slopes_chosen_silva_tax$treatment <- reg_all_slopes_chosen_silva_tax$treatment %>%
  factor(levels=(c("CL", "CD", "PL", "PD", "DL", "VL")))

reg_all_slopes_chosen_silva_tax %$%
  class %>%
  unique()

# growth_rates_treatment <- reg_all_slopes_chosen_silva_tax %>%
#   filter(pvalue_slope_chosen < 0.05 &
#            slope_chosen_days > 0) %>%
#   ggplot(aes(treatment, slope_chosen_days, fill = interaction(fct_infreq(class_f),phylum_f, sep = '\n')))+
#   geom_col(aes(fill = fct_infreq(class_f)))+
#   scale_fill_manual(values = palf_large_phylums(40))+
#   #scale_y_continuous(labels = percent_format())+
#   facet_grid(season~.)+
#   #scale_y_continuous(labels = percent_format())+
#   labs(x = 'Treatments', fill = 'Phylum')+
#   #scale_y_continuous(percent_format())+
#   theme_minimal()+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = 'right')

##intento dibuixar la gr relativa (la contribució de cada ASV al creixement de la comunitat)
sum_gr <- reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05 &
           slope_chosen_days > 0) %>%
  group_by(treatment, season) %>%
  summarize(sum_gr = sum(slope_chosen_days))

reg_all_slopes_chosen_silva_tax_filt_rel_gr <- reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05 &
           slope_chosen_days > 0) %>%
  left_join(sum_gr) %>%
  mutate(relative_growth_rate = slope_chosen_days/sum_gr)

reg_all_slopes_chosen_silva_tax_filt_rel_gr %>%
  dim()

reg_all_slopes_chosen_silva_tax_filt_rel_gr %>%
  colnames()
reg_all_slopes_chosen_silva_tax_filt_rel_gr$class_f %>%
  levels()

insitu_phy_obj_melt$class_f %>%
  levels()

reg_all_slopes_chosen_silva_tax_filt_rel_gr %>%
  left_join(insitu_phy_obj_melt, keep = TRUE)

order_df <- insitu_phy_obj_melt %>%
  arrange(class_f)%>%
  pull(class_f)

reg_all_slopes_chosen_silva_tax_filt_rel_gr %$%
  class_f %>%
  levels()

##reorder classes manually perquè estiguin igual les de les growth rates i les dels in situ i a totes les figures
reg_all_slopes_chosen_silva_tax_filt_rel_gr$phylum_f <- reg_all_slopes_chosen_silva_tax_filt_rel_gr$phylum_f %>%
  factor(levels = c("Proteobacteria",  "Bacteroidota", "Planctomycetota", "Actinobacteriota", 
                     "Crenarchaeota", "Verrucomicrobiota", "Cyanobacteria", "Firmicutes", "Acidobacteriota",
                     "Nitrospinota", "Bdellovibrionota", "Nitrospirota",  "Myxococcota", "Desulfobacterota",
                     "Deinococcota" , "Chloroflexi" , "Fusobacteriota", 
                     "Spirochaetota", "Campilobacterota", "Abditibacteriota", "Latescibacterota",
                     "Methylomirabilota", "Halanaerobiaeota", "Sumerlaeota", "Calditrichota", "Gemmatimonadota"), 
          ordered = TRUE)
reg_all_slopes_chosen_silva_tax_filt_rel_gr %$%
  class_f %>%
  unique()

reg_all_slopes_chosen_silva_tax_filt_rel_gr$class_f <- reg_all_slopes_chosen_silva_tax_filt_rel_gr$class_f %>%
  factor(levels = c('Gammaproteobacteria', 'Alphaproteobacteria', 'Zetaproteobacteria',  'Bacteroidia', 'Rhodothermia',  
                    'Ignavibacteria', 'Chlorobia', 'Kryptonia', 'Nitrososphaeria',   'Cyanobacteriia', 
                    'Vampirivibrionia',   'Acidimicrobiia',   'Actinobacteria', 'Coriobacteriia',  'Thermoleophilia',
                    'Verrucomicrobiae',  'Kiritimatiellae', 'Lentisphaeria', 'Omnitrophia', 'Chlamydiae', 
                    'Planctomycetes',  'Phycisphaerae',   'Blastocatellia',   'Holophagae', 'Vicinamibacteria', 
                    'Acidobacteriae',  'Thermoanaerobaculia',  'Bdellovibrionia',  'Oligoflexia',  'Bacilli', 
                    'Clostridia', 'Negativicutes', 'Syntrophomonadia', 'Desulfitobacteriia', 'Thermoanaerobacteria', 
                    'Myxococcia', 'Polyangia',  'Nitrospinia',  'Campylobacteria',   'Deinococci', 
                    'Fusobacteriia', 'Desulfovibrionia', 'Desulfobulbia', 'Desulfuromonadia',  'Desulfobacteria',  
                    'Methylomirabilia', 'Gemmatimonadetes', 'Anaerolineae', 'Chloroflexia', 'Spirochaetia', 
                    'Leptospirae', 'Calditrichia', 'Halobacteria',  'Nitrospiria', 'Babeliae', 
                    'Saccharimonadia', 'Cloacimonadia', 'Synergistia', 'Abditibacteria',  'Deferribacteres'), 
         ordered = TRUE)

##VALORS REPETITS A LA LLEGENDA
reg_all_slopes_chosen_silva_tax_filt_rel_gr %>%
  colnames()
relative_gr <- 
  reg_all_slopes_chosen_silva_tax_filt_rel_gr %>%
  ungroup() %>%
  # group_by(season, class_f, treatment, relative_growth_rate) %>%
  # summarize(sum_rel_gr = sum(relative_growth_rate)) %>%
  left_join(calculated_uncalculated_gr) %>%
  #distinct(treatment, season, asv_num, relative_growth_rate, class_f, percentages, phylum_f) %>% #fct_reorder(class_f, phylum))
  ggplot(aes(treatment, relative_growth_rate))+##fct_infreq  interaction(class_f,  phylum_f, sep = '\n')) interaction(fct_relevel(class_f, as.character(order_df))), phylum_f, drop=TRUE)
  geom_col(aes(fill = fct_rev(class_f)))+
  scale_fill_manual(values = palette_class_assigned)+
                    # labels = interaction( reg_all_slopes_chosen_silva_tax_filt_rel_gr$phylum_f, 
                    #                       reg_all_slopes_chosen_silva_tax_filt_rel_gr$class_f, sep = '\n'))+
  scale_y_continuous(labels = percent_format())+
  facet_grid(season~.)+
  labs(x = 'Treatments', fill = 'Phylum - Class', y='Relative growth rate')+
  #guides(fill=guide_legend(ncol = 3))+
  #legend(legend = levels(factor(class_f)))+
  #scale_y_continuous(percent_format())+
  geom_point(aes(treatment, percentages), color = 'black')+ #intento afegir els percentatges a sobre del gràfic
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = 'right',
        legend.spacing = unit(10, 'cm'), legend.text = element_text(size = 8),
        legend.spacing.x = unit(0.2, 'cm'))

reg_all_slopes_chosen_silva_tax %$%
  phylum %>%
  unique()

##growth community vs uncalculable gr-------

reg_all_slopes_chosen_silva_tax

insitu_phy_obj_melt %$%
  phylum %>%
  unique()

x<- reg_all_slopes_chosen_silva_tax %>%
  subset(slope_chosen_days == is.na(slope_chosen_days) | 
           pvalue_slope_chosen > 0.05) %>%
  group_by(treatment, season, phylum_f) %>%
  summarize(nas = n())

y <- reg_all_slopes_chosen_silva_tax %>%
  subset(slope_chosen_days != is.na(slope_chosen_days) |
           pvalue_slope_chosen < 0.05 ) %>%
  group_by(treatment, season, phylum_f) %>%
  summarize(non_nas = n())

dim(x)
dim(y)

calculated_uncalculated_gr <- x %>% left_join(y) %>%
  mutate(percentage_nas = nas/(nas + non_nas),
         percentage_non_nas = non_nas/(nas + non_nas)) %>%
  pivot_longer(names_to = 'gr_calculated_percentage', 
               values_to = 'percentages',
               cols = starts_with('percentage')) %>%
  filter(gr_calculated_percentage == 'percentage_non_nas') #em quedo amb el percentage que he pogut calcular de cada comunitat 

##en general (no per phylums)----
x<- reg_all_slopes_chosen_silva_tax %>%
  subset(slope_chosen_days == is.na(slope_chosen_days) | 
           pvalue_slope_chosen > 0.05) %>%
  group_by(treatment, season) %>%
  summarize(nas = n())

y <- reg_all_slopes_chosen_silva_tax %>%
  subset(slope_chosen_days != is.na(slope_chosen_days) |
           pvalue_slope_chosen < 0.05 ) %>%
  group_by(treatment, season) %>%
  summarize(non_nas = n())

dim(x)
dim(y)

calculated_uncalculated_gr <- x %>% left_join(y) %>%
  mutate(percentage_nas = nas/(nas + non_nas),
         percentage_non_nas = non_nas/(nas + non_nas)) %>%
  pivot_longer(names_to = 'gr_calculated_percentage', 
               values_to = 'percentages',
               cols = starts_with('percentage'))

##figure construction -----
gr_vs_insitu <- multi_panel_figure(columns = 4, rows = 1, width = 250, height = 150, 
                                             column_spacing = 0.0, unit = 'mm',
                                             panel_label_type = 'none')

# gr_vs_insitu  %<>%
#   fill_panel(insitu_community, column = 1, row = 1) %<>%
#   fill_panel(growth_rates_treatment, column = 2:3, row = 1)

#en relatiu
gr_vs_insitu  %<>%
  fill_panel(insitu_community, column = 1, row = 1) %<>%
  fill_panel(relative_gr, column = 2:4, row = 1)

ggsave('rel_gr_vs_insitu_assigned_colors.pdf', gr_vs_insitu  , 
       path = "~/Documentos/Doctorat/REMEI/results/figures/raw/",
       width = 250,
       height = 150,
       units = 'mm')


 #### A nivell d'ordres/classes -------
reg_all_slopes_chosen_silva_tax_filt %$%
  class %>%
  unique() #25

reg_all_slopes_chosen_silva_tax_filt %$%
  order %>%
  unique() #75

##intento dibuixar la gr relativa (la contribució de cada ASV al creixement de la comunitat)
sum_gr <- reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05 &
           slope_chosen_days > 0) %>%
  group_by(treatment, season) %>%
  summarize(sum_gr = sum(slope_chosen_days))

reg_all_slopes_chosen_silva_tax_filt_rel_gr <- reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05 &
           slope_chosen_days > 0) %>%
  left_join(sum_gr) %>%
  mutate(relative_growth_rate = slope_chosen_days/sum_gr)

reg_all_slopes_chosen_silva_tax_filt_rel_gr %>%
  dim()
relative_gr <- reg_all_slopes_chosen_silva_tax_filt_rel_gr %>%
  ungroup() %>%
  left_join(calculated_uncalculated_gr) %>%
  #distinct(treatment, season, asv_num, relative_growth_rate, class_f, percentages) %>%
  ggplot(aes(treatment, relative_growth_rate, fill = fct_infreq(class_f)))+
  geom_col(aes())+
  scale_fill_manual(values = palf_large_phylums(40))+
  scale_y_continuous(labels = percent_format())+
  facet_grid(season~.)+
  labs(x = 'Treatments', fill = 'Phylum')+
  #scale_y_continuous(percent_format())+
  geom_point(aes(treatment, percentages), color = 'black')+ #intento afegir els percentatges a sobre del gràfic
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = 'right')

reg_all_slopes_chosen_silva_tax %$%
  phylum %>%
  unique()

##growth community vs uncalculable gr-------

reg_all_slopes_chosen_silva_tax

insitu_phy_obj_melt %$%
  phylum %>%
  unique()

x<- reg_all_slopes_chosen_silva_tax %>%
  subset(slope_chosen_days == is.na(slope_chosen_days) | 
           pvalue_slope_chosen > 0.05) %>%
  group_by(treatment, season, phylum_f) %>%
  summarize(nas = n())

y <- reg_all_slopes_chosen_silva_tax %>%
  subset(slope_chosen_days != is.na(slope_chosen_days) |
           pvalue_slope_chosen < 0.05 ) %>%
  group_by(treatment, season, phylum_f) %>%
  summarize(non_nas = n())

dim(x)
dim(y)

calculated_uncalculated_gr <- x %>% left_join(y) %>%
  mutate(percentage_nas = nas/(nas + non_nas),
         percentage_non_nas = non_nas/(nas + non_nas)) %>%
  pivot_longer(names_to = 'gr_calculated_percentage', 
               values_to = 'percentages',
               cols = starts_with('percentage')) %>%
  filter(gr_calculated_percentage == 'percentage_non_nas') #em quedo amb el percentage que he pogut calcular de cada comunitat 

##en general (no per phylums)
x<- reg_all_slopes_chosen_silva_tax %>%
  subset(slope_chosen_days == is.na(slope_chosen_days) | 
           pvalue_slope_chosen > 0.05) %>%
  group_by(treatment, season) %>%
  summarize(nas = n())

y <- reg_all_slopes_chosen_silva_tax %>%
  subset(slope_chosen_days != is.na(slope_chosen_days) |
           pvalue_slope_chosen < 0.05 ) %>%
  group_by(treatment, season) %>%
  summarize(non_nas = n())

dim(x)
dim(y)

calculated_uncalculated_gr <- x %>% left_join(y) %>%
  mutate(percentage_nas = nas/(nas + non_nas),
         percentage_non_nas = non_nas/(nas + non_nas)) %>%
  pivot_longer(names_to = 'gr_calculated_percentage', 
               values_to = 'percentages',
               cols = starts_with('percentage'))

#growth_rates_treatment <- 
calculated_uncalculated_gr %>%
  group_by(phylum_f) %>%
  filter(n() > 10) %>%
  ggplot(aes(treatment, percentages, fill = gr_calculated_percentage))+ #fct_infreq(phylum_f)
  geom_col(aes())+
  scale_fill_manual(values = palette_seasons)+
  #scale_y_continuous(labels = percent_format())+
  facet_grid(phylum_f~season)+
  scale_y_continuous(labels = percent_format())+
  labs(x = 'Treatments', fill = 'Calculated uncalculated')+
  #scale_y_continuous(percent_format())+
  theme_minimal()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = 'right')

##figure 
gr_vs_insitu <- multi_panel_figure(columns = 3, rows = 1, width = 250, height = 150, 
                                   column_spacing = 0.0, unit = 'mm',
                                   panel_label_type = 'none')

gr_vs_insitu  %<>%
  fill_panel(insitu_community, column = 1, row = 1) %<>%
  fill_panel(growth_rates_treatment, column = 2:3, row = 1)

#en relatiu
gr_vs_insitu  %<>%
  fill_panel(insitu_community, column = 1, row = 1) %<>%
  fill_panel(relative_gr, column = 2:3, row = 1)

ggsave('rel_gr_vs_insitu.pdf', gr_vs_insitu  , 
       path = "~/Documentos/Doctorat/REMEI/results/figures/raw/",
       width = 250,
       height = 150,
       units = 'mm')

##Sankey  diagram-----
library(networkD3)

sankeyNetwork(Links = Energy$links, Nodes = Energy$nodes, Source = "source",
              Target = "target", Value = "value", NodeID = "name",
              units = "TWh", fontSize = 12, nodeWidth = 30)

##entre seasons i treatments
data_long <- reg_all_slopes_chosen_silva_tax_filt %>%
  #mutate(season_treatment = paste(season,'_',treatment)) %>%
  select(season, treatment, slope_chosen) #phylum_f, class_f, order_f, family_f

data_long %>%
  head()

colnames(data_long) <- c("source", "target", "value")
data_long$target <- paste(data_long$target, " ", sep="")
# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(name=c(as.character(data_long$source), as.character(data_long$target)) %>% unique())

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
data_long$IDsource=match(data_long$source, nodes$name)-1 
data_long$IDtarget=match(data_long$target, nodes$name)-1

# prepare colour scale
ColourScal ='d3.scaleOrdinal() .range(["#FDE725FF","#B4DE2CFF","#6DCD59FF","#35B779FF","#1F9E89FF","#26828EFF","#31688EFF","#3E4A89FF","#482878FF","#440154FF"])'

# Make the Network
sankeyNetwork(Links = data_long, Nodes = nodes,
              Source = "IDsource", Target = "IDtarget",
              Value = "value", NodeID = "name", 
              sinksRight=FALSE, colourScale=ColourScal, nodeWidth=40, fontSize=13, nodePadding=20)

##Chord diagram
##vull veure quin nº d'ASV es comparteixen entre treatments-seasons-insitu
library(circlize)

numbers <- sample(c(1:1000), 100, replace = T)
data <- matrix(numbers, ncol=5)
rownames(data) <- paste0("orig-", seq(1,20))
colnames(data) <- paste0("dest-", seq(1,5))

# Make the circular plot
chordDiagram(data, transparency = 0.5)

# Create an edge list: a list of connections between 10 origin nodes, and 10 destination nodes:
origin <- paste0("orig ", sample(c(1:10), 20, replace = T))
destination <- paste0("dest ", sample(c(1:10), 20, replace = T))
data <- data.frame(origin, destination)

# Transform input data in a adjacency matrix
adjacencyData <- with(data, table(origin, destination))

# Charge the circlize library
library(circlize)

# Make the circular plot
chordDiagram(adjacencyData, transparency = 0.5)

#Another try
library(tidyverse)
library(viridis)
library(patchwork)
library(hrbrthemes)
library(circlize)
library(chorddiag) 
reg_all_slopes_chosen_silva_tax_filt %>%
  colnames()
data_long <- reg_all_slopes_chosen_silva_tax_filt %>%
  mutate(season_treatment = paste(season,'_',treatment)) %>%
  select(season_treatment, slope_chosen) #phylum_f, class_f, order_f, family_f

data_long %>%
  head()

data_long %$%
  season_treatment %>%
  unique()

# parameters
circos.clear()
circos.par(start.degree = 90, gap.degree = 0, track.margin = c(-0.1, 0.1), points.overflow.warning = FALSE)
par(mar = rep(0, 4))

# color palette
mycolor <- viridis(24, alpha = 1, begin = 0, end = 1, option = "D")
mycolor <- mycolor[sample(1:24)]

chordDiagram(
  x = data_long, 
  grid.col = 24,
  transparency = 0.25,
  directional = 1,
  direction.type = c("arrows", "diffHeight"), 
  diffHeight  = -0.04,
  annotationTrack = "grid", 
  annotationTrackHeight = c(0.05, 0.1),
  link.arr.type = "big.arrow", 
  link.sort = TRUE, 
  link.largest.ontop = TRUE)

# Add text and axis
circos.trackPlotRegion(
  track.index = 1, 
  bg.border = NA, 
  panel.fun = function(x, y) {
    
    xlim = get.cell.meta.data("xlim")
    sector.index = get.cell.meta.data("sector.index")
    
    # Add names to the sector. 
    circos.text(
      x = mean(xlim), 
      y = 3.2, 
      labels = sector.index, 
      facing = "bending", 
      cex = 0.8
    )
    
    # Add graduation on axis
    circos.axis(
      h = "top", 
      major.at = seq(from = 0, to = xlim[2], by = ifelse(test = xlim[2]>10, yes = 2, no = 1)), 
      minor.ticks = 1, 
      major.tick.percentage = 0.5,
      labels.niceFacing = FALSE)
  }
)

##chord diagram https://jokergoo.github.io/circlize_book/book/advanced-usage-of-chorddiagram.html -------
##computacionalment lent de crear el gràfic-----
circos.clear()
circos.info() # chordDiagram() creates two tracks, one track for labels and one track for grids with axes

#These two tracks can be controlled by annotationTrack argument. Available values for this argument are grid, 
#name and axis. The height of annotation tracks can be set by annotationTrackHeight which is the percentage 
#to the radius of unit circle and can be set by mm_h() function with an absolute unit. Axes are only added 
#if grid is set in annotationTrack
#From verion 0.4.10 of the circlize package, there is a new group argument in chordDiagram() function which 
#is very convenient for making multiple-group Chord diagrams.

#dcast converts a dataframe to a matrix (we need a matrix for a chordplot)
data_long <- reg_all_slopes_chosen_silva_tax_filt %>%
  mutate(season_treatment = paste(season,'_',treatment)) %>%
  select(season_treatment, slope_chosen, asv_num, season, treatment) #phylum_f, class_f, order_f, family_f

data_long %>%
  head()
library(reshape2)
data_long %>%
  pivot_wider(id_cols = c('asv_num', 'season', 'treatment'), values_from = slope_chosen, names_from = c('asv_num', 'season', 'treatment'))

nm = unique(unlist(dimnames(data_long)))
group = structure(gsub("\\d", "", nm), names = nm)
group



#una altra opció------
### Define ranges of circos sectors and their colors (both of the sectors and the links)-----
data_long$xmin <- 0
data_long$xmax <- rowSums(data_long) + colSums(data_long)
n <- nrow(df1)
df1$rcol<-rgb(df1$r, df1$g, df1$b, max = 255)
df1$lcol<-rgb(df1$r, df1$g, df1$b, alpha=200, max = 255)

### Plot sectors (outer part)
par(mar=rep(0,4))
circos.clear()
### Basic circos graphic parameters
circos.par(cell.padding=c(0,0,0,0), track.margin=c(0,0.15), start.degree = 90, gap.degree =4)

### Sector details
circos.initialize(factors = data_long$season, x = data_long$slope_chosen) #, xlim = cbind(data_long$slope_chosen, data_long$slope_chosen)

### Plot sectors
circos.trackPlotRegion(ylim = c(0, 1), factors = df1$country, track.height=0.1,
                       #panel.fun for each sector
                       panel.fun = function(x, y) {
                         #select details of current sector
                         name = get.cell.meta.data("sector.index")
                         i = get.cell.meta.data("sector.numeric.index")
                         xlim = get.cell.meta.data("xlim")
                         ylim = get.cell.meta.data("ylim")
                         
                         #text direction (dd) and adjusmtents (aa)
                         theta = circlize(mean(xlim), 1.3)[1, 1] %% 360
                         dd <- ifelse(theta < 90 || theta > 270, "clockwise", "reverse.clockwise")
                         aa = c(1, 0.5)
                         if(theta < 90 || theta > 270)  aa = c(0, 0.5)
                         
                         #plot country labels
                         circos.text(x=mean(xlim), y=1.7, labels=name, facing = dd, cex=0.6,  adj = aa)
                         
                         #plot main sector
                         circos.rect(xleft=xlim[1], ybottom=ylim[1], xright=xlim[2], ytop=ylim[2], 
                                     col = df1$rcol[i], border=df1$rcol[i])
                         
                         #blank in part of main sector
                         circos.rect(xleft=xlim[1], ybottom=ylim[1], xright=xlim[2]-rowSums(m)[i], ytop=ylim[1]+0.3, 
                                     col = "white", border = "white")
                         
                         #white line all the way around
                         circos.rect(xleft=xlim[1], ybottom=0.3, xright=xlim[2], ytop=0.32, col = "white", border = "white")
                         
                         #plot axis
                         circos.axis(labels.cex=0.6, direction = "outside", major.at=seq(from=0,to=floor(df1$xmax)[i],by=5), 
                                     minor.ticks=1, labels.away.percentage = 0.15)
                       })

### Plot links (inner part)
### Add sum values to df1, marking the x-position of the first links
### out (sum1) and in (sum2). Updated for further links in loop below.
df1$sum1 <- colSums(m)
df1$sum2 <- numeric(n)

### Create a data.frame of the flow matrix sorted by flow size, to allow largest flow plotted first
df2 <- cbind(as.data.frame(m),orig=rownames(m),  stringsAsFactors=FALSE)
df2 <- reshape(df2, idvar="orig", varying=list(1:n), direction="long",
               timevar="dest", time=rownames(m),  v.names = "m")
df2 <- arrange(df2,desc(m))

### Keep only the largest flows to avoid clutter
df2 <- subset(df2, m > quantile(m,0.6))

### Plot links
for(k in 1:nrow(df2)){
  #i,j reference of flow matrix
  i<-match(df2$orig[k],df1$country)
  j<-match(df2$dest[k],df1$country)
  
  #plot link
  circos.link(sector.index1=df1$country[i], point1=c(df1$sum1[i], df1$sum1[i] + abs(m[i, j])),
              sector.index2=df1$country[j], point2=c(df1$sum2[j], df1$sum2[j] + abs(m[i, j])),
              col = df1$lcol[i])
  
  #update sum1 and sum2 for use when plotting the next link
  df1$sum1[i] = df1$sum1[i] + abs(m[i, j])
  df1$sum2[j] = df1$sum2[j] + abs(m[i, j])
}

### nested proportions parallel sets plot


####Relation ASV in situ - ASV growing in REMEI--------
reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05 &
           slope_chosen_days > 0)

insitu_phy_obj %>%
  dim()

insitu_phy_obj_melt %>%
  head()
reg_all_slopes_chosen_silva_tax_filt %>%
  dim()

reg_all_slopes_chosen_silva_tax_filt %>%
  select(treatment, season, asv_num, slope_chosen_days, phylum) %>%
  left_join(insitu_phy_obj_melt, by = c('season', 'asv_num', 'phylum')) %>%
  ggplot(aes(Abundance, slope_chosen_days, color = phylum))+
  geom_point(aes())+
  scale_color_manual(values = palette_phylums_assigned)+
  scale_x_continuous(labels = percent_format())+
  facet_grid(~season)+
  theme_bw()

abundance_gr_dataset <- reg_all_slopes_chosen_silva_tax_1perc %>%
  filter(pvalue_slope_chosen < 0.05 &
           slope_chosen_days > 0) %>%
  select(treatment, season, asv_num, slope_chosen_days, phylum, class) %>%
  right_join(insitu_phy_obj_melt, by = c('season', 'asv_num', 'phylum', 'class')) %>%
  filter(treatment.x != is.na(treatment.x)) %>%
  mutate(behaviour = case_when(Abundance < 0.01 & slope_chosen_days > 1 ~ 'Rare responsive',
                               Abundance > 0.01 & slope_chosen_days > 1 ~ 'Abundant responsive',
                               Abundance < 0.01 & slope_chosen_days < 1 ~ 'Rare unresponsive',
                               Abundance > 0.01 & slope_chosen_days < 1 ~ 'Abundant unresponsive'))
         
  #          pvalue_slope_kept =case_when(slope.x > slope.y & slope.pval.x < 0.05 ~ slope.pval.x,
  #                                       slope.x < slope.y ~ slope.pval.y))
  ###vull saber quin % d'ASV rares creixen als tractaments 
  # mutate(number_non_rare_growing = count(Abundance < 0.01 & slope_chosen_days > 1)),
  #        total = n(Abundance != is.na(Abundance)))
  # filter(Abundance > 0.01 & slope_chosen_days > 1)
abundance_gr_dataset %>%
  filter(pvalue_slope_chosen < 0.05 & slope_chosen_days > 0) %>%
  ggplot(aes(Abundance, slope_chosen_days, color = class, group = behaviour))+
  geom_point(aes(), size = 2, alpha = 0.8)+
  labs(y = expression("Growth rate day"^"-1"), x = 'ASV relative abundance in situ community', fill = 'Class')+
  guides(color=guide_legend(ncol = 5))+
  scale_color_manual(values = palette_class_assigned)+
  scale_x_continuous(labels = percent_format(), expand = c(0,0))+
  geom_hline(yintercept = 1)+
  geom_vline(xintercept = 0.01)+
  scale_y_continuous(expand = c(0, 0), limits = c(0,14))+ #forçar que comenci al 0 l'eix
  #geom_polygon(aes(group = behaviour))+
  #geom_polygon(aes(x=0.01, y = 1), fill = 'black')+
  facet_grid(~season~treatment.x)+
  #annotate(geom = 'text', x = 4, y = 5, label = 'Abundant responsive')+#, element_text(size = 2)
  #annotate(geom = 'text', x = 4, y = 0.7, label = 'Abundant unresponsive')+
  theme_bw()+
  theme(strip.text.x = element_text(size = 12), axis.text.x = element_text(angle = 0), legend.position = "bottom",
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.text = element_text(size = 12),
        legend.title = element_text(size = 16))


##amb els non growing or not significant (canviats a 0)-------
reg_all_slopes_chosen_silva_tax_1perc %>%
  colnames()

abundance_gr_dataset <- reg_all_slopes_chosen_silva_tax_1perc %>%
  mutate(all_gr = case_when(pvalue_slope_chosen > 0.05 ~ 0,
           TRUE ~ slope_chosen_days)) %>%
  select(treatment, season, asv_num, all_gr, phylum, class, pvalue_slope_chosen, slope_chosen_days) %>%
  right_join(insitu_phy_obj_melt, by = c('season', 'asv_num', 'phylum', 'class')) %>%
  filter(treatment.x != is.na(treatment.x)) %>%
  mutate(behaviour = case_when(Abundance < 0.01 & slope_chosen_days > 1 ~ 'Rare responsive',
                               Abundance > 0.01 & slope_chosen_days > 1 ~ 'Abundant responsive',
                               Abundance < 0.01 & slope_chosen_days < 1 ~ 'Rare unresponsive',
                               Abundance > 0.01 & slope_chosen_days < 1 ~ 'Abundant unresponsive'))

#          pvalue_slope_kept =case_when(slope.x > slope.y & slope.pval.x < 0.05 ~ slope.pval.x,
#                                       slope.x < slope.y ~ slope.pval.y))
###vull saber quin % d'ASV rares creixen als tractaments 
# mutate(number_non_rare_growing = count(Abundance < 0.01 & slope_chosen_days > 1)),
#        total = n(Abundance != is.na(Abundance)))
# filter(Abundance > 0.01 & slope_chosen_days > 1)
abundance_gr_dataset %>%
  filter(pvalue_slope_chosen < 0.05) %>%
  ggplot(aes(Abundance, slope_chosen_days, color = class, group = behaviour))+
  geom_point(aes(), size = 2, alpha = 0.8)+
  labs(y = expression("Growth rate day"^"-1"), x = expression(paste('ASV relative abundance', italic(' in situ '), 'community')), fill = 'Class')+
  guides(color=guide_legend(ncol = 5))+
  scale_color_manual(values = palette_class_assigned)+
  scale_x_continuous(labels = percent_format())+
  geom_hline(yintercept = 0.05)+
  geom_vline(xintercept = 0.01)+
  #geom_area()+
  #geom_polygon(aes(x=0.01, y = 1), fill = 'black')+
  facet_grid(~season~treatment.x)+
  #annotate(geom = 'text', x = 4, y = 5, label = 'Abundant responsive')+#, element_text(size = 2)
  #annotate(geom = 'text', x = 4, y = 0.7, label = 'Abundant unresponsive')+
  theme_bw()+
  theme(strip.text.x = element_text(size = 12), axis.text.x = element_text(angle = 0), legend.position = "bottom",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.text = element_text(size = 12),
        legend.title = element_text(size = 16))

test <- abundance_gr_dataset %>%
  filter(pvalue_slope_chosen < 0.05) %>%
  group_by(season, slope_chosen_days, Abundance, asv_num, treatment.x) %>%
  filter(slope_chosen_days != is.na(slope_chosen_days) &
    
           Abundance > 0) %>%
  group_by(asv_num, season) %T>%
  mutate(n = n()) %>%
  summarize(n = n())

#blast ASV insitu amb ASV REMEI dataset align profiles (mirar si hi ha molta diferència entre  els 2 gràfics)--------
library(DECIPHER)

##import data (treballo amb les asv de miau dataset sol (el que va processar la Marta), no amb el conjunt meu REMIAU).
remei_seqs <- read.fasta('./data/dada2/02_nochimera_mergeruns/remei_1_2_pool/remei_1_2_pool_seqtab_final.fasta')
miau_seqs <- read.fasta('../MIAU_seqs/MIAU_runs_seqtab_final.fasta')

##?entrar dades en .fasta

dna_REMEI<- readDNAStringSet(remei_seqs)
dna_MIAU <- readDNAStringSet(miau_seqs)

dna_REMEI_MIAU <- paste(dna_REMEI, dna_MIAU)
dna_REMEI_MIAU <- readDNAStringSet(dna_REMEI_MIAU)

alignment_1 <- 
  DECIPHER::AlignSeqs(dna_REMEI_MIAU, 
                      processors = NULL)

library()

head(dna_REMEI)

test <- match(dna_MIAU, dna_REMEI, nomatch = 0) %>%
  as_tibble() ##no funciona  bé

test <- match(remei_seqs, miau_seqs, nomatch = 0) %>%
  as_tibble()

left_join(remei_tax, remiau_tax)

remei_tax <- remei_tax %>%
  as_tibble()

remiau_tax <- remiau_tax %>%
  as_tibble()

remiau_tax %>%
  head()

test <- match(remei_tax$.otu, remiau_tax$.otu, nomatch = 0) %>%
  as_tibble(9)

test <- remei_tax %>%
  filter(.otu %in% remiau_tax$.otu)

remei_tax %>%
  dim()

test %>%
  dim()


#FIGURE WITH INNER_JOIN BETWEEN ASV SEQUENCES FROM IN SITU DATA (MIAU DATASET) & GROWTH RATES FROM REMEI---------
###ONLY ASV PRESENT IN SITU AND GROWING ASV (SIGNIFICATIVELY)
##Prepare MIAU dataset-----
###Import data
miau_seqs <- read.fasta('../MIAU_seqs/MIAU_runs_seqtab_final.fasta')

miau_asv_tab <- read_rds('../MIAU_seqs/MIAU_runs_seqtab_final.rds')  %>%
  as_tibble(rownames = 'sample_code')

tax_silva_miau <- readRDS("../MIAU_seqs/MIAUruns_assignTax_tax_assignation.rds") %>% 
  as_tibble(rownames = 'sequence')

##filtro primer les asv de MIAU que no són 0 a cap dels meus temps in situ per fer més fàcil la relació entre datasets
miau_asv_tab_insitu <- miau_asv_tab %>%
  filter(sample_code %in% c("x254", "x255", "x256", "x257"))

miau_asv_tab_insitu_codes <- miau_asv_tab_insitu$sample_code

miau_asv_tab_insitu_filt <- miau_asv_tab_insitu %>%
  select(-sample_code) %>%
  select(where(~ sum(.) != 0)) %>%
  cbind(miau_asv_tab_insitu_codes)

miau_seqs <- miau_asv_tab_insitu_filt %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = 'ASV_seq') %>%
  row_to_names(row_number = 3874, remove_row = T, remove_rows_above = F)

#filtrem seqs per taxonomia (sense phylum o que son chloroplasts o mitocondria)
tax_silva_miau_filt <- tax_silva_miau %>% 
  #Filter all results without domain assignation 
  dplyr::filter(!is.na('Kingdom'), !is.na('Phylum')) %>% 
  #And the Chloros/Mitochondria seqs
  dplyr::filter(Order !=  'Chloroplast') %>% 
  dplyr::filter(Family !=  'Mitochondria') 

#filter miau seqs by only the one's that are not unclassified or chloroplasts and mitochondria
miau_seqs %>%
  colnames()
tax_silva_miau_filt %>%
  colnames()
miau_seqs_filt  <- miau_seqs %>%
  filter(miau_asv_tab_insitu_codes %in% tax_silva_miau_filt$sequence) #2548 ASV

miau_seqs_filt_ed <- miau_seqs_filt %>%
  mutate(across(!miau_asv_tab_insitu_codes, as.numeric))

miau_seqs_filt_ed_tax <- miau_seqs_filt_ed %>%
  left_join(tax_silva_miau_filt, by = c('miau_asv_tab_insitu_codes' = 'sequence')) ##add taxonomy to compare with exp tax

##fem relative abundances at in situ time in MIAU dataset after filtering ASV
#miau_seqs_filt_relabund <- 
 # miau_seqs_filt_ed %>%
 #   colnames()
# miau_seqs_filt %>%
#   glimpse()
miau_seqs_filt_rel_abund <- apply(miau_seqs_filt_ed[,c(2:5)], 2,   function(x){x / sum(x)}) %>%
  cbind(miau_seqs_filt_ed$miau_asv_tab_insitu_codes) %>%
  as_tibble()

colnames(miau_seqs_filt_rel_abund) <- c("x254", "x255", "x256", "x257", "asv_seqs")
# miau_seqs_filt_rel_abund %>% 
#   colnames()

miau_seqs_filt_rel_abund_long <- miau_seqs_filt_rel_abund %>%
  pivot_longer(cols = starts_with('x'), names_to = 'sample', values_to = 'rel_abund') %>%
  mutate(season = case_when(sample == 'x254' ~ 'Winter',
                            sample == 'x255' ~ 'Spring',
                            sample == 'x256' ~ 'Summer',
                            sample == 'x257' ~ 'Fall'),
         rel_abund = as.numeric(rel_abund))

#prepare REMEI dataset -----
rem_fc_filt_silva <- readRDS("data/intermediate_files/remei_phyloseq_silva_fc.rds") 

rem_fc_filt_silva@tax_table <- rem_fc_filt_silva@tax_table %>% 
  mutate_tax_table(asv_num = str_c( 'asv' , 1:ncol(rem_fc_filt_silva@otu_table)))

remei_tax <- rem_fc_filt_silva@tax_table %>%
  as_tibble()
# head(remei_tax)

##afegim les dades de growth rates
reg_all_slopes_chosen_silva_tax <- read.csv("data/intermediate_files/reg_all_slopes_chosen_silva_tax.csv", sep=",") %>%
  filter(season != "Early_fall")
reg_all_slopes_chosen_silva_tax_filt <- reg_all_slopes_chosen_silva_tax %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05)

##unim les seqüències de les asv a les dades de growth rates
growth_rates_remei_sequences <- reg_all_slopes_chosen_silva_tax_filt %>%
  inner_join(remei_tax,  by = 'asv_num', copy = TRUE)%>%
  as_tibble()

# growth_rates_remei_sequences %>%
#   class()
# 
# reg_all_slopes_chosen_silva_tax_filt %$%
#   asv_num %>%
#   unique() #1863

##filter MIAU ASV sequences by REMEI -----
# remei_tax %>%
#   colnames()
# remei_tax %>%
#   head()
# miau_seqs_filt_rel_abund %>%
#   class()
# remei_tax %>%
#   dim()
# miau_seqs_filt_rel_abund %>%
#   colnames()
# miau_seqs_filt_rel_abund_long %>%
#   colnames()
# growth_rates_remei_sequences%>%
#   colnames()
# miau_seqs_filt_rel_abund_long %>%
#   dim()
# 
# miau_seqs_filt_rel_abund_long %$% 
#   asv_seqs %>%
#   unique() #2548
# 
# growth_rates_remei_sequences %>%
#   dim() ##5543
# 
# growth_rates_remei_sequences %$%
#   asv_num %>%
#   unique() #1863

# growth_rates_remei_sequences %>%
#    colnames()
#  miau_seqs_filt_rel_abund_long %>%
#    colnames()
# miau_seqs_filt_rel_abund_long %>%
#   dim()
# growth_rates_remei_sequences %>%
#   dim()
##volem que totes les gr tinguin la seva abundancia a in situ time (si no la tenen posarem un 0)
remei_gr_in_situ_data  <- miau_seqs_filt_rel_abund_long %>%
  right_join(growth_rates_remei_sequences, by = c('asv_seqs' = '.otu', 'season' = 'season')) 

remei_gr_in_situ_data %>%
  filter(treatment != is.na(treatment))%$% 
  asv_num %>%
  unique() #1288
  
# remei_gr_in_situ_data %$%
#   asv_num %>%
#   unique() #1815

#filter(miau_seqs_filt_rel_abund$asv_seqs %in%  growth_rates_remei_sequences$.otu)

# remei_gr_in_situ_data %>%
#   colnames()
# remei_gr_in_situ_data %>%
#   glimpse()

remei_gr_in_situ_data <- remei_gr_in_situ_data  %>%
    filter(treatment != is.na(treatment),
           rel_abund_ed != is.na(rel_abund_ed)) %>%
  ungroup() %>%
  mutate(rel_abund_ed = as.numeric(rel_abund_ed))  %>%
  mutate(rel_abund_ed2 = round(rel_abund_ed, 3)) %>%
  #arrange(rel_abund_ed, .by_group = FALSE) %>%#season.x, treatment, #ordeno per la rel abundance a in situ
  #group_by(treatment, season.x) %>% ##season només 1 mostra per tant no cal que agrupem per estació
  mutate(#rank_abund = as.numeric(order(rel_abund_ed, decreasing = TRUE)), ##general per tot el dataset #no está ben ordenat
         rank_abund_2 = rank(-rel_abund_ed),
         behaviour = case_when(rel_abund_ed < 0.001 ~ 'Rare x < 0.1% ',
                               rel_abund_ed > 0.01 ~ 'Abundant x > 1%',
                               rel_abund_ed < 0.01  & rel_abund_ed >0.001 ~ 'Mid  0.1% < x < 1%')) 


remei_gr_in_situ_data$treatment <- remei_gr_in_situ_data$treatment %>% 
  factor(levels=(c("CL", "CD", "PL", "PD", "DL", "VL", 'NA')))
remei_gr_in_situ_data$season <- remei_gr_in_situ_data$season %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))

# remei_gr_in_situ_data %>%
#   colnames()
# remei_gr_in_situ_data %>%
#   filter(rel_abund != is.na(rel_abund) &
#            slope_chosen_days != is.na(slope_chosen_days)) %$%
#   asv_num %>%
#   unique()
# remei_gr_in_situ_data %>%
#   colnames()
#  
# sequences_rel_abund_gr %>%
#   colnames()
# 
# sequences_rel_abund_gr %$%
#   asv_num_col2 %>%
#   unique()

# seqs_gr_rel_abund %>%
#   filter(rel_abund != (is.na(rel_abund)),
#          slope_chosen_days != (is.na(slope_chosen_days))) %$% ##RECUPERO 391 ASV GR and IN SITU 
#   relation_miau_remei1 %>%
#   unique()
# test <- seqs_gr_rel_abund %>%
#   subset(treatment == 'CD' & season.x == 'Spring')
# seqs_gr_rel_abund$rel_abund_ed2 %>%
#   class()
# seqs_gr_rel_abund$rel_abund_ed == 0.01
# seqs_gr_rel_abund %$% 
#   rank_abund %>%
#   class()
# seqs_gr_rel_abund %>%
#   glimpse()

##table with % of rare reactive ASVs (add to plots as geom_text)-------
##in relation to all potential bloomers
remei_gr_in_situ_data %>%
  colnames()

total_asv_match_remei_miau <- remei_gr_in_situ_data %>%
  filter(rel_abund_ed < 0.01 & slope_chosen_days > 2) %>%
  group_by(season, treatment) %>%
  summarize(total_common_asv_potential_bloomers = n())

total_asv_match_remei_miau_potential_bloomers <- remei_gr_in_situ_data %>%
  group_by(season, treatment) %>%
  summarize(total_common_asv = n()) %>%
  select(-season) %>%
  cbind(total_asv_match_remei_miau) %>%
  select(season...1, treatment...2, total_common_asv, total_common_asv_potential_bloomers)

colnames(total_asv_match_remei_miau_potential_bloomers) <- c('season', 'treatment', 'total_common_asv',
                                                             'total_common_asv_potential_bloomers')

rare_reactive_asv <- remei_gr_in_situ_data %>%
  filter(rel_abund_ed < 0.001 & slope_chosen_days > 2) %>%
  group_by(season, treatment) %>%
  summarize(rare_reactive_asv_sum = n())

summarize_match_remei_miau <- total_asv_match_remei_miau_potential_bloomers %>%
  full_join(rare_reactive_asv) %>%
  mutate(rare_reactive_perc = round(rare_reactive_asv_sum/total_common_asv_potential_bloomers*100, 2))
# as_tibble() %>%
# mutate(rare_reactive_asv_ed = as.numeric(rare_reactive_asv))  %>%
# mutate( ~replace_na(rare_reactive_asv_sum, 0))
summarize_match_remei_miau <- summarize_match_remei_miau %>%
  mutate(rare_reactive_perc_ed = ifelse(is.na(rare_reactive_perc), 0, rare_reactive_perc),
         rare_reactive_asv_sum = ifelse(is.na(rare_reactive_asv_sum), 0, rare_reactive_asv_sum))

growing_asv <- reg_all_slopes_chosen_silva_tax_filt %>%
  group_by(season, treatment) %>%
  summarize(total_growing_asv = n())

miau_seqs_filt_rel_abund_long %>%
  colnames()

insitu_asv <- miau_seqs_filt_rel_abund_long %>%
  filter(rel_abund > 0) %>%
  group_by(season) %>%
  summarize(insitu_present_asv = n())

summarize_insitu_vs_growing <- summarize_match_remei_miau %>%
  left_join(growing_asv, by = c('season' = 'season', 'treatment' = 'treatment')) %>%
  left_join(insitu_asv, by = c('season' = 'season')) %>%
  select(-rare_reactive_perc)

#write.table(summarize_insitu_vs_growing, 'results/tables/summarize_insitu_vs_growing_ed2.txt', sep='\t')

summarize_insitu_vs_growing %>%
  colnames()

summarize_insitu_vs_growing %>%
  mutate(perc_common_asv = total_common_asv/insitu_present_asv) %>%
  select(-rare_reactive_perc_ed, - perc_common_asv) %>%
  pivot_longer(cols = c(total_common_asv, rare_reactive_asv_sum, total_growing_asv, insitu_present_asv)) %>%
  ggplot(aes(name, value))+
  geom_col()+
  facet_grid(season~treatment)+
  theme_bw()+
  theme(strip.text = element_text(size = 16), axis.text.x = element_text(angle = 0), legend.position = "bottom",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.text = element_text(size = 10),
        legend.title = element_text(size = 12), strip.background = element_blank()) 


summarize_insitu_vs_growing$treatment <- summarize_insitu_vs_growing$treatment %>% 
  factor(levels=(c("CL", "CD", "PL", "PD", "DL", "VL")))

summarize_insitu_vs_growing$season <- summarize_insitu_vs_growing$season %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))

percentage_shared_asv_insitu_experiments <- summarize_insitu_vs_growing %>%
  mutate(perc_common_asv = total_common_asv/insitu_present_asv,
         perc_rare_reactive_asv = rare_reactive_asv_sum/total_common_asv_potential_bloomers)%>%
  #select(starts_with('perc'), season.x...1, treatment...2) %>%
  select(perc_common_asv, season, treatment)%>%
  pivot_longer(cols = starts_with('perc')) %>%
  ggplot(aes(x= season, y = value, fill = season))+
  geom_bar(stat='identity', width = 1)+
  scale_fill_manual(values = palette_seasons_4)+
  #coord_polar('y', start = 0)+
  scale_y_continuous(labels = percent_format())+
  facet_grid(.~ treatment)+
  labs(fill = 'Season', y = 'shared ASVs between in situ and experiments', x= 'Season')+
  theme_bw()+
  theme(strip.text = element_text(size = 16), legend.position = "top",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.text = element_text(size = 10),
        legend.title = element_text(size = 12), strip.background = element_blank(), axis.title.y = element_text(size = 6),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title = element_blank(), 
        axis.text.y = element_text(size = 6), plot.margin= unit(c(1, 2, 0, 0.2), "lines"), panel.border = element_blank()) 

#PLOT rank abund and gr -----
# summarize_match_remei_miau$rare_reactive_perc_ed
# summarize_match_remei_miau %>%
#   glimpse()
# remei_gr_in_situ_data %>%
#   colnames()
# 
# remei_gr_in_situ_data %>%
#   dim()

rel_abund_gr_relation_rank <- remei_gr_in_situ_data %>%
  # left_join(summarize_match_remei_miau) %>%
  #filter(pvalue_slope_chosen < 0.05) %>%
  #filter(season.x != is.na(season.x)) %>%
  ggplot(aes(log10(rank_abund_2), slope_chosen_days, color = class.x))+ #, group = behaviour
  geom_text(data = summarize_match_remei_miau, mapping = aes(x = 3, y = 2.8, 
                                                             label = as.character(rare_reactive_perc_ed)), 
            check_overlap = TRUE, na.rm = TRUE, show.legend = FALSE, nudge_x = 0, nudge_y = 6, 
            color = 'black', size = 4)+
  geom_point(aes(shape =  behaviour), size = 2.5, alpha = 0.8)+
  geom_rect(aes(xmin = 0, xmax = Inf, ymin = 0, ymax = 2), 
            fill = '#D0CFC8', show.legend = FALSE, linejoin = NULL,  alpha = 0.02, color = NA)+
  geom_rect(aes(xmin = 0, xmax = log10(rank_abund_2[match(0.01, rel_abund_ed2)]), ymin = 2, ymax = Inf), color = NA, 
            fill = '#D0CFC8', show.legend = FALSE, linejoin = NULL,  alpha = 0.02)+
  geom_vline(aes(xintercept = log10(rank_abund_2[match(0.01, rel_abund_ed2)])), linetype = 'dashed')+
  geom_vline(aes(xintercept = log10(rank_abund_2[match(0.001, rel_abund_ed2)])), linetype = 'dashed')+
  # geom_rect(aes(xmin = 0, xmax = log10(rank_abund[match(0.001, rel_abund_ed2)]), ymin = 2.5, ymax = Inf), 
  #           fill = '#D0CFC8', show.legend = FALSE, linejoin = NULL)+
  # geom_rect(aes(xmin = 0, xmax = log10(rank_abund[match(0.001, rel_abund_ed2)]), ymin = 0, ymax = 2.5), 
  #           fill = '#D0CFC8', show.legend = FALSE, linejoin = NULL, color = '#D0CFC8')+
  geom_hline(yintercept = 2, linetype = 'dashed')+
  labs(y = expression("Growth rate day"^"-1"), 
       x = 'Log10(rank abundance)', color = 'Class',
       shape = 'In situ\nrelative\nabundance (%)')+ #expression(paste(italic('In situ'),'\n relative abundance (%)' em queda mal col·locat
  guides(color=guide_legend(ncol = 4, size = 5),
         shape = guide_legend(ncol = 1, size = 3))+
  scale_color_manual(values = palette_class_assigned)+
  #scale_x_continuous(expand = c(0, 0), limits = c(0, 4))+
  # geom_ribbon(aes(xmin = aes(log10(rank_abund[match(0.001, rel_abund_ed2)])), ymin = 2.5, 
  #             xmax = rank_abund, ymax = slope_chosen_days))+
  #scale_x_continuous(labels = percent_format())+
  #geom_hline(yintercept = 0.05)+
  # geom_vline(xintercept = log(663))+
  # geom_vline(xintercept = log(185))+
  # geom_vline(xintercept = seqs_gr_rel_abund$rel_abund_ed2 = 0.010)+
  #geom_vline(xintercept = seqs_gr_rel_abund$rank_abund[match(0.001, seqs_gr_rel_abund$rel_abund_ed2)])+
  #geom_area(aes(as.numeric(slope_chosen_days) > 2.5))+
  #geom_polygon(aes(x=0.01, y = 1), fill = 'black')+
  facet_grid(~season~treatment)+
  #annotate(geom = 'text', x = 4, y = 5, label = 'Abundant responsive')+#, element_text(size = 2)
  #annotate(geom = 'text', x = 4, y = 0.7, label = 'Abundant unresponsive')+
  theme_bw()+
  theme(strip.text = element_text(size = 16), axis.text.x = element_text(angle = 0), legend.position = "bottom",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.text = element_text(size = 10),
        legend.title = element_text(size = 14), strip.background = element_blank(), #strip.text.x = element_blank(),
        panel.border = element_blank()) 
# 
# ggsave(filename = 'shared_asv_insitu_experiments_innerjoin_ed2.pdf', plot = rel_abund_gr_relation_rank, device = NULL, path = 'results/figures/raw/',
#        width = 33, height = 30, units = 'cm')  
# 

##taxonomy check-------
remei_gr_in_situ_data %>%
  colnames()
phylum_check <- remei_gr_in_situ_data$phylum.x == remei_gr_in_situ_data$phylum.y
summary(phylum_check)
class <- remei_gr_in_situ_data$class.x == remei_gr_in_situ_data$class.y
summary(class)
order <- remei_gr_in_situ_data$order.x == remei_gr_in_situ_data$order.y
summary(order)
family <- remei_gr_in_situ_data$family.x == remei_gr_in_situ_data$family.y
summary(family)
remei_gr_in_situ_data$asv_seqs %>%
  unique()
remei_gr_in_situ_data$sequence_miau %>%
  unique()

taxonomy_check <- remei_gr_in_situ_data$sequence_remei == remei_gr_in_situ_data$sequence_miau
summary(taxonomy_check)
#Same graph but adding a 0 to the growing ASV that don't have a match in situ-----
##considering that this could be 
##due to sequencing depth
##volem que totes les gr tinguin la seva abundancia a in situ time (si no la tenen posarem un 0)
remei_gr_in_situ_data  <- miau_seqs_filt_rel_abund_long %>%
  right_join(growth_rates_remei_sequences, by = c('asv_seqs' = '.otu', 'season' = 'season')) 

remei_gr_in_situ_data  <- remei_gr_in_situ_data %>%
  mutate(rel_abund_ed = if_else(is.na(rel_abund), 0, rel_abund)) #si no tens rel_abun la posem a 0 considerant que has de ser-hi per créixer

remei_gr_in_situ_data <- remei_gr_in_situ_data  %>%
  # filter(#treatment != is.na(treatment),
  #        rel_abund_ed != is.na(rel_abund_ed)) %>%
  ungroup() %>%
  mutate(rel_abund_ed = as.numeric(rel_abund_ed))  %>%
  mutate(rel_abund_ed2 = round(rel_abund_ed, 3)) %>%
  #arrange(rel_abund_ed, .by_group = FALSE) %>%#season.x, treatment, #ordeno per la rel abundance a in situ
  #group_by(treatment, season.x) %>% ##season només 1 mostra per tant no cal que agrupem per estació
  mutate(#rank_abund = as.numeric(order(rel_abund_ed, decreasing = TRUE)), ##general per tot el dataset (no está ben ordenat)
         rank_abund_2 = rank(-rel_abund_ed),
         behaviour = case_when(rel_abund_ed < 0.001 ~ 'Rare x < 0.1% ',
                               rel_abund_ed > 0.01 ~ 'Abundant x > 1%',
                               rel_abund_ed < 0.01  & rel_abund_ed >0.001 ~ 'Mid  0.1% < x < 1%')) 

##table with % of rare reactive ASVs (add to plots as geom_text)-------
##in relation to all potential bloomers
##in relation to all potential bloomers
remei_gr_in_situ_data %>%
  colnames()

total_asv_match_remei_miau <- remei_gr_in_situ_data %>%
  filter(rel_abund_ed < 0.01 & slope_chosen_days > 2) %>%
  group_by(season, treatment) %>%
  summarize(total_common_asv_potential_bloomers = n())

total_asv_match_remei_miau_potential_bloomers <- remei_gr_in_situ_data %>%
  group_by(season, treatment) %>%
  summarize(total_common_asv = n()) %>%
  select(-season) %>%
  cbind(total_asv_match_remei_miau) %>%
  select(season...1, treatment...2, total_common_asv, total_common_asv_potential_bloomers)

colnames(total_asv_match_remei_miau_potential_bloomers) <- c('season', 'treatment', 'total_common_asv',
                                                             'total_common_asv_potential_bloomers')

rare_reactive_asv <- remei_gr_in_situ_data %>%
  filter(rel_abund_ed < 0.001 & slope_chosen_days > 2) %>%
  group_by(season, treatment) %>%
  summarize(rare_reactive_asv_sum = n())

summarize_match_remei_miau <- total_asv_match_remei_miau_potential_bloomers %>%
  full_join(rare_reactive_asv) %>%
  mutate(rare_reactive_perc = round(rare_reactive_asv_sum/total_common_asv_potential_bloomers*100, 2))
# as_tibble() %>%
# mutate(rare_reactive_asv_ed = as.numeric(rare_reactive_asv))  %>%
# mutate( ~replace_na(rare_reactive_asv_sum, 0))
summarize_match_remei_miau <- summarize_match_remei_miau %>%
  mutate(rare_reactive_perc_ed = ifelse(is.na(rare_reactive_perc), 0, rare_reactive_perc),
         rare_reactive_asv_sum = ifelse(is.na(rare_reactive_asv_sum), 0, rare_reactive_asv_sum))

growing_asv <- reg_all_slopes_chosen_silva_tax_filt %>%
  group_by(season, treatment) %>%
  summarize(total_growing_asv = n())

miau_seqs_filt_rel_abund_long %>%
  colnames()

insitu_asv <- miau_seqs_filt_rel_abund_long %>%
  filter(rel_abund > 0) %>%
  group_by(season) %>%
  summarize(insitu_present_asv = n())

summarize_insitu_vs_growing <- summarize_match_remei_miau %>%
  left_join(growing_asv, by = c('season' = 'season', 'treatment' = 'treatment')) %>%
  left_join(insitu_asv, by = c('season' = 'season')) %>%
  select(-rare_reactive_perc)

#write.table(summarize_insitu_vs_growing, 'results/tables/summarize_insitu_vs_growing_ed2.txt', sep='\t')

summarize_insitu_vs_growing %>%
  colnames()

summarize_insitu_vs_growing %>%
  mutate(perc_common_asv = total_common_asv/insitu_present_asv) %>%
  select(-rare_reactive_perc_ed, - perc_common_asv) %>%
  pivot_longer(cols = c(total_common_asv, rare_reactive_asv_sum, total_growing_asv, insitu_present_asv)) %>%
  ggplot(aes(name, value))+
  geom_col()+
  facet_grid(season~treatment)+
  theme_bw()+
  theme(strip.text = element_text(size = 16), axis.text.x = element_text(angle = 0), legend.position = "bottom",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.text = element_text(size = 10),
        legend.title = element_text(size = 12), strip.background = element_blank()) 

summarize_insitu_vs_growing$treatment <- summarize_insitu_vs_growing$treatment %>% 
  factor(levels=(c("CL", "CD", "PL", "PD", "DL", "VL")))

summarize_insitu_vs_growing$season <- summarize_insitu_vs_growing$season %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))

percentage_shared_asv_insitu_experiments <- summarize_insitu_vs_growing %>%
  mutate(perc_common_asv = total_common_asv/insitu_present_asv,
         perc_rare_reactive_asv = rare_reactive_asv_sum/total_common_asv_potential_bloomers)%>%
  #select(starts_with('perc'), season.x...1, treatment...2) %>%
  select(perc_common_asv, season, treatment)%>%
  pivot_longer(cols = starts_with('perc')) %>%
  ggplot(aes(x= season, y = value, fill = season))+
  geom_bar(stat='identity', width = 1)+
  scale_fill_manual(values = palette_seasons_4)+
  #coord_polar('y', start = 0)+
  scale_y_continuous(labels = percent_format())+
  facet_grid(.~ treatment)+
  labs(fill = 'Season', y = 'shared ASVs between in situ and experiments', x= 'Season')+
  theme_bw()+
  theme(strip.text = element_text(size = 16), legend.position = "top",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.text = element_text(size = 10),
        legend.title = element_text(size = 12), strip.background = element_blank(), axis.title.y = element_text(size = 6),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title = element_blank(), 
        axis.text.y = element_text(size = 6), plot.margin= unit(c(1, 2, 0, 0.2), "lines"), panel.border = element_blank()) 

summarize_insitu_vs_growing$treatment <- summarize_insitu_vs_growing$treatment %>% 
  factor(levels=(c("CL", "CD", "PL", "PD", "DL", "VL")))

summarize_insitu_vs_growing$season <- summarize_insitu_vs_growing$season %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))

percentage_shared_asv_insitu_experiments <- summarize_insitu_vs_growing %>%
  mutate(perc_common_asv = total_common_asv/insitu_present_asv,
         perc_rare_reactive_asv = rare_reactive_asv_sum/total_common_asv_potential_bloomers)%>%
  #select(starts_with('perc'), season.x...1, treatment...2) %>%
  select(perc_common_asv, season, treatment)%>%
  pivot_longer(cols = starts_with('perc')) %>%
  ggplot(aes(x= season, y = value, fill = season))+
  geom_bar(stat='identity', width = 1)+
  scale_fill_manual(values = palette_seasons_4)+
  #coord_polar('y', start = 0)+
  scale_y_continuous(labels = percent_format())+
  facet_grid(.~ treatment)+
  labs(fill = 'Season', y = 'shared ASVs between in situ and experiments', x= 'Season')+
  theme_bw()+
  theme(strip.text = element_text(size = 16), legend.position = "top",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.text = element_text(size = 10),
        legend.title = element_text(size = 12), strip.background = element_blank(), axis.title.y = element_text(size = 6),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title = element_blank(), 
        axis.text.y = element_text(size = 6), plot.margin= unit(c(1, 2, 0, 0.2), "lines"), panel.border = element_blank()) 

#rel_abund_gr_relation_rank <-
# summarize_match_remei_miau$rare_reactive_perc_ed
# summarize_match_remei_miau %>%
#   glimpse()
# remei_gr_in_situ_data %>%
#   colnames()
# remei_gr_in_situ_data %>%
#   dim()

##plot rank abund vs gr with 0 rel_abund in situ-------
rel_abund_gr_relation_rank <- remei_gr_in_situ_data %>%
  # left_join(summarize_match_remei_miau) %>%
  #filter(pvalue_slope_chosen < 0.05) %>%
  #filter(season.x != is.na(season.x)) %>%
  ggplot(aes(log10(rank_abund_2), slope_chosen_days, color = class.x))+ #, group = behaviour
  geom_text(data = summarize_match_remei_miau, mapping = aes(x = 3, y = 2.8, 
                                                             label = as.character(rare_reactive_perc_ed)), 
            check_overlap = TRUE, na.rm = TRUE, show.legend = FALSE, nudge_x = 0, nudge_y = 6, 
            color = 'black', size = 4)+
  geom_point(aes(shape =  behaviour), size = 2.5, alpha = 0.8)+
  geom_rect(aes(xmin = 0, xmax = Inf, ymin = 0, ymax = 2), 
            fill = '#D0CFC8', show.legend = FALSE, linejoin = NULL,  alpha = 0.02, color = NA)+
  geom_rect(aes(xmin = 0, xmax = log10(rank_abund_2[match(0.01, rel_abund_ed2)]), ymin = 2, ymax = Inf), color = NA, 
            fill = '#D0CFC8', show.legend = FALSE, linejoin = NULL,  alpha = 0.02)+
  geom_vline(aes(xintercept = log10(rank_abund_2[match(0.01, rel_abund_ed2)])), linetype = 'dashed')+
  geom_vline(aes(xintercept = log10(rank_abund_2[match(0.001, rel_abund_ed2)])), linetype = 'dashed')+
  # geom_rect(aes(xmin = 0, xmax = log10(rank_abund[match(0.001, rel_abund_ed2)]), ymin = 2.5, ymax = Inf), 
  #           fill = '#D0CFC8', show.legend = FALSE, linejoin = NULL)+
  # geom_rect(aes(xmin = 0, xmax = log10(rank_abund[match(0.001, rel_abund_ed2)]), ymin = 0, ymax = 2.5), 
  #           fill = '#D0CFC8', show.legend = FALSE, linejoin = NULL, color = '#D0CFC8')+
  geom_hline(yintercept = 2, linetype = 'dashed')+
  labs(y = expression("Growth rate day"^"-1"), 
       x = 'Log10(rank abundance)', color = 'Class',
       shape = 'In situ\nrelative\nabundance (%)')+ #expression(paste(italic('In situ'),'\n relative abundance (%)' em queda mal col·locat
  guides(color=guide_legend(ncol = 4, size = 5),
         shape = guide_legend(ncol = 1, size = 3))+
  scale_color_manual(values = palette_class_assigned)+
  #scale_x_continuous(expand = c(0, 0), limits = c(0, 4))+
  # geom_ribbon(aes(xmin = aes(log10(rank_abund[match(0.001, rel_abund_ed2)])), ymin = 2.5, 
  #             xmax = rank_abund, ymax = slope_chosen_days))+
  #scale_x_continuous(labels = percent_format())+
  #geom_hline(yintercept = 0.05)+
  # geom_vline(xintercept = log(663))+
  # geom_vline(xintercept = log(185))+
  # geom_vline(xintercept = seqs_gr_rel_abund$rel_abund_ed2 = 0.010)+
  #geom_vline(xintercept = seqs_gr_rel_abund$rank_abund[match(0.001, seqs_gr_rel_abund$rel_abund_ed2)])+
  #geom_area(aes(as.numeric(slope_chosen_days) > 2.5))+
  #geom_polygon(aes(x=0.01, y = 1), fill = 'black')+
facet_grid(~season~treatment)+
  #annotate(geom = 'text', x = 4, y = 5, label = 'Abundant responsive')+#, element_text(size = 2)
  #annotate(geom = 'text', x = 4, y = 0.7, label = 'Abundant unresponsive')+
  theme_bw()+
  theme(strip.text = element_text(size = 16), axis.text.x = element_text(angle = 0), legend.position = "bottom",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.text = element_text(size = 10),
        legend.title = element_text(size = 14), strip.background = element_blank(), #strip.text.x = element_blank(),
        panel.border = element_blank()) 
# 
ggsave(filename = 'shared_asv_insitu_experiments_fulljoin_0rel_abund_ed2.pdf', plot = rel_abund_gr_relation_rank, device = NULL, path = 'results/figures/raw/',
       width = 33, height = 30, units = 'cm')




#FIGURE WITH CLUSTERING AT 100% BETWEEN ASV SEQUENCES FROM IN SITU DATA (MIAU DATASET) & GROWTH RATES FROM REMEI---------
##for these purpuse we used vsearch function cluster_fast, which allows us to clusterize the fasta sequences in filename, 
##automatically sort by decreasing sequence length beforehand.
###ONLY ASV PRESENT IN SITU AND GROWING ASV (SIGNIFICATIVELY)
#Load results data from the clustering at 100% 
cluster_miau_remei <- read.table('data/cluster_remei_miau_asv/clustering_results.txt', header = F)
cluster_miau_remei %>%
  dim() #5196 asv relacionades entre un dataset i l'altre 
#(moltes més que en l'altre cas, perquè filtro abans del cluster, podria ser)
colnames(cluster_miau_remei) <- c('v1', 'v2', 'v3', 'clustering', 'v5', 'v6', 'v7', 'v8', 'col1_asv', 'col2_asv')
str_count(cluster_miau_remei$col1_asv, "remei") %>%
  sum() #3642
str_count(cluster_miau_remei$col1_asv, "miau") %>%
  sum() #1554
str_count(cluster_miau_remei$col2_asv, "remei") %>%
  sum() #1554
str_count(cluster_miau_remei$col2_asv, "miau") %>%
  sum() #3642

# cluster_miau_remei %>%
#   head()
# cluster_miau_remei %$%
#   col1_asv %>%
#   unique()
# cluster_miau_remei %$%
# asv_num_col2 %>%
#   unique()
cluster_miau_remei <- cluster_miau_remei %>%
  separate(col1_asv, c('asv_num_col1', 'size_col1'), sep = ';') %>%
  separate(col2_asv, c('asv_num_col2', 'size_col2'), sep = ';')   #hi han asv de remei a la columna 1 i 2 i el mateix per miau


cluster_miau_remei <- cluster_miau_remei %>%
  mutate(relation_miau_remei1 = paste(asv_num_col1,'-',asv_num_col2)
         #relation_miau_remei2 = paste(asv_num_col2,'-',asv_num_col1)
         )

#cluster_miau_remei$relation_miau_remei1 == cluster_miau_remei$relation_miau_remei2
cluster_miau_remei$relation_miau_remei1 %>%
  unique() #5196
# cluster_miau_remei$relation_miau_remei2 %>%
#   unique()
cluster_miau_remei  %>%
  dim()
##col1 and col2 tenen els mateixos match però estan repetits

##haig d'afegir la seqüencia al resultat del clustering per poder relacionar-ho amb les rel_abun i les gr 
##(ja que el número de l'asv podria no ser el mateix amb el que estic treballant en el meu dataset)
miau_seqs <- read.fasta('data/cluster_remei_miau_asv/MIAU_runs_seqtab_final.fasta')
remei_seqs <- read.table('data/cluster_remei_miau_asv/remei_1_2_seqtab_final.fasta')

dna_REMEI<- readDNAStringSet('data/cluster_remei_miau_asv/remei_1_2_seqtab_final.fasta')
seq_name = names(dna_REMEI)
sequence = paste(dna_REMEI)
remei_seqs <- data.frame(seq_name, sequence) %>%
  as_tibble() %>%
  mutate(asv_num_ed = str_replace(seq_name, ';', 'remei;')) %>%
  separate(asv_num_ed, c('asv_num', 'size'), sep = ';')

dna_miau <- readDNAStringSet('data/cluster_remei_miau_asv/MIAU_runs_seqtab_final.fasta')
seq_name = names(dna_miau)
sequence = paste(dna_miau)
miau_seqs <- data.frame(seq_name, sequence) %>%
  as_tibble() %>%
  mutate(asv_num_ed = str_replace(seq_name, ';', 'miau;')) %>%
  separate(asv_num_ed, c('asv_num', 'size'), sep = ';')
# miau_seqs %>%
#   head()
# 
# remei_seqs %>%
#   head()
# 
# remei_seqs %>%
#   colnames()
# 
# cluster_miau_remei %>%
#   head()
# 
# cluster_miau_remei %>%
#   colnames()
# 
# cluster_miau_remei %$%
#   remei_asv
# 
# remei_seqs %$%
#   asv_num_ed

#faig una columna al cluster miau_remei només amb els asv num de remei i una amb només asv num de miau
#així puc unir les sequencies sense que se'm dupliquin les columnes.
library(stringi)
cluster_miau_remei <- cluster_miau_remei %>%
  mutate(miau_asv_num = as.character(str_extract_all(relation_miau_remei1, pattern = ('asv[0-9]*(miau)+'))),
         remei_asv_num = as.character(str_extract_all(relation_miau_remei1, pattern = ('asv[0-9]*(remei)+')))) 
           
cluster_remei_seqs <- cluster_miau_remei %>%
  #full_join(c, by = c('asv_num_col1' = 'asv_num')) %>%
  left_join(remei_seqs, by = c('remei_asv_num' = 'asv_num'))

cluster_remei_seqs %>%
  colnames()

cluster_remei_seqs <- cluster_remei_seqs %>%
  rename(seq_name = 'seq_name_remei',
         sequence = 'sequence_remei',
         size = 'size_remei')
# cluster_remei_seqs  %$%
#   asv_num_col1 %$%
#   str_detect(str_detect('remei'))
#   # full_join(remei_seqs, by = c('asv_num_col2' = 'asv_num')) %>%
#   # filter(clustering != is.na(clustering),
#   #        sequence.x !=  is.na(sequence.x))

cluster_miau_seqs <- cluster_miau_remei %>%
  #full_join(c, by = c('asv_num_col1' = 'asv_num')) %>%
  left_join(miau_seqs, by = c('miau_asv_num' = 'asv_num'))

cluster_miau_seqs <- cluster_miau_seqs %>%
  rename(seq_name = 'seq_name_miau',
         sequence = 'sequence_miau',
         size = 'size_remei')

cluster_miau_remei %$%
  asv_num_col1 %>%
  unique() #5196

cluster_miau_remei %$%
  asv_num_col2 %>%
  unique() #5196

#prepare growth rates REMEI dataset -----
rem_fc_filt_silva <- readRDS("data/intermediate_files/remei_phyloseq_silva_fc.rds") 

rem_fc_filt_silva@tax_table <- rem_fc_filt_silva@tax_table %>% 
  mutate_tax_table(asv_num = str_c( 'asv' , 1:ncol(rem_fc_filt_silva@otu_table)))

remei_tax <- rem_fc_filt_silva@tax_table %>%
  as_tibble()
# head(remei_tax)

##afegim les dades de growth rates
reg_all_slopes_chosen_silva_tax <- read.csv("data/intermediate_files/reg_all_slopes_chosen_silva_tax.csv", sep=",") %>%
  filter(season != "Early_fall")
reg_all_slopes_chosen_silva_tax_filt <- reg_all_slopes_chosen_silva_tax %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05)

##unim les seqüències de les asv a les dades de growth rates
growth_rates_remei_seqs <- reg_all_slopes_chosen_silva_tax_filt %>%
  inner_join(remei_tax,  by = 'asv_num', copy = TRUE)%>%
  as_tibble()

growth_rates_remei_seqs %>%
  colnames()
##Prepare MIAU dataset-----
###Import data
# miau_seqs <- read.fasta('../MIAU_seqs/MIAU_runs_seqtab_final.fasta')

miau_asv_tab <- read_rds('../MIAU_seqs/MIAU_runs_seqtab_final.rds')  %>%
  as_tibble(rownames = 'sample_code')

tax_silva_miau <- readRDS("../MIAU_seqs/MIAUruns_assignTax_tax_assignation.rds") %>% 
  as_tibble(rownames = 'sequence')

##filtro primer les asv de MIAU que no són 0 a cap dels meus temps in situ per fer més fàcil la relació entre datasets
miau_asv_tab_insitu <- miau_asv_tab %>%
  filter(sample_code %in% c("x254", "x255", "x256", "x257"))

miau_asv_tab_insitu_codes <- miau_asv_tab_insitu$sample_code

miau_asv_tab_insitu_filt <- miau_asv_tab_insitu %>%
  select(-sample_code) %>%
  select(where(~ sum(.) != 0)) %>%
  cbind(miau_asv_tab_insitu_codes)

miau_seqs <- miau_asv_tab_insitu_filt %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = 'ASV_seq') %>%
  row_to_names(row_number = 3874, remove_row = T, remove_rows_above = F)

#filtrem seqs per taxonomia (sense phylum o son chloroplasts o mitocondria)
tax_silva_miau_filt <- tax_silva_miau %>% 
  #Filter all results without domain assignation 
  dplyr::filter(!is.na('Kingdom'), !is.na('Phylum')) %>% 
  #And the Chloros/Mitochondria seqs
  dplyr::filter(Order !=  'Chloroplast') %>% 
  dplyr::filter(Family !=  'Mitochondria') 

#filter miau seqs by only the one's that are not unclassified or chloroplasts and mitochondria
miau_seqs_filt  <- miau_seqs %>%
  filter(miau_asv_tab_insitu_codes %in% tax_silva_miau_filt$sequence) #2548 ASV

miau_seqs_filt_ed <- miau_seqs_filt %>%
  mutate(across(!miau_asv_tab_insitu_codes, as.numeric))

##fem relative abundances at in situ time in MIAU dataset after filtering ASV
#miau_seqs_filt_relabund <- 
# miau_seqs_filt_ed %>%
#   colnames()
# miau_seqs_filt %>%
#   glimpse()
miau_seqs_filt_rel_abund <- apply(miau_seqs_filt_ed[,c(2:5)], 2,   function(x){x / sum(x)}) %>%
  cbind(miau_seqs_filt_ed$miau_asv_tab_insitu_codes) %>%
  as_tibble()

colnames(miau_seqs_filt_rel_abund) <- c("x254", "x255", "x256", "x257", "asv_seqs")
# miau_seqs_filt_rel_abund %>% 
#   colnames()

miau_seqs_filt_rel_abund_long <- miau_seqs_filt_rel_abund %>%
  pivot_longer(cols = starts_with('x'), names_to = 'sample', values_to = 'rel_abund') %>%
  mutate(season = case_when(sample == 'x254' ~ 'Winter',
                            sample == 'x255' ~ 'Spring',
                            sample == 'x256' ~ 'Summer',
                            sample == 'x257' ~ 'Fall'),
         rel_abund = as.numeric(rel_abund))


miau_seqs_filt_rel_abund_long_tax <- miau_seqs_filt_rel_abund_long %>%
  left_join(tax_silva_miau_filt, by = c('asv_seqs' = 'sequence')) ##add taxonomy to check if its the same for both datasets

##Unim GR data amb in situ rel abund amb les seqüències-------
##unim les  relative abundances a les dades del cluster amb les seqs de MIAU
sequences_rel_abund <- cluster_miau_seqs %>%
  right_join(miau_seqs_filt_rel_abund_long_tax, by = c('sequence_miau' = 'asv_seqs')) %>%
  filter(relation_miau_remei1 != is.na(relation_miau_remei1)) #eliminem les que no tenen relació amb l'altre dataset

miau_seqs_filt_rel_abund_long_tax %>%
  dim()
sequences_rel_abund %>%
  dim() 

cluster_miau_seqs %>%
  colnames()
miau_seqs_filt_rel_abund_long %>%
  colnames()
sequences_rel_abund %>%
  colnames()
sequences_rel_abund %>%
  dim()
sequences_rel_abund %$%
  relation_miau_remei1 %>%
  unique() ###1515 relations amb GR

##unim les gr a les dades del cluster amb les seqs de REMEI
sequences_gr <- cluster_remei_seqs %>%
  right_join(growth_rates_remei_seqs, by = c('sequence_remei' = '.otu')) %>%
  filter(relation_miau_remei1 != is.na(relation_miau_remei1)) #eliminem les que no tenen relació amb l'altre dataset

sequences_gr %>%
  colnames()

sequences_gr %$%
  asv_num_col1 %>%
  unique() #759 asv

# sequences_rel_abund %$% 
#   sequence
# 
# sequences_gr %$% 
#   sequence

#unim les growth rates amb les dades de rel_abundance 
seqs_gr_rel_abund <- sequences_gr %>%
  full_join(sequences_rel_abund, by = 'relation_miau_remei1') %>%
  filter(rel_abund != is.na(rel_abund) &
           slope_chosen_days != is.na(slope_chosen_days)) ##es podria subsitutir el is.na per 0 (mirar al següent script)

##taxonomy check-------
phylum_check <- seqs_gr_rel_abund$Phylum == seqs_gr_rel_abund$phylum.x
summary(phylum_check)
class <- seqs_gr_rel_abund$Class == seqs_gr_rel_abund$class.x
summary(class)
order <- seqs_gr_rel_abund$Order == seqs_gr_rel_abund$order.x
summary(order)
seqs_gr_rel_abund$Family == seqs_gr_rel_abund$family.x
seqs_gr_rel_abund$sequence_remei %>%
  unique()
seqs_gr_rel_abund$sequence_miau %>%
  unique()

taxonomy_check <- seqs_gr_rel_abund$sequence_remei == seqs_gr_rel_abund$sequence_miau
summary(taxonomy_check)
#---------
#seqs_gr_rel_abund1$relation_miau_remei1

# seqs_gr_rel_abund %>%
#   colnames()

# seqs_gr_rel_abund2 <- seqs_gr_rel_abund1 %>%
#   full_join(sequences_rel_abund, by = c('relation_miau_remei1' = 'relation_miau_remei1')) %>%
#   filter(rel_abund.x != is.na(rel_abund.x) &
#            slope_chosen_days != is.na(slope_chosen_days))

seqs_gr_rel_abund %$%
  asv_num_col1.x %>%
  unique() ##391 ASV comunes amb gr i abundancia a in situ

# sequences_gr %>%
#   colnames()
# # seqs_gr_rel_abund1 %>%
#   colnames()

# seqs_rel_abund_gr %>%
#   dim()
miau_seqs_filt_rel_abund_long %>%
  colnames()
reg_all_slopes_chosen_silva_tax_filt %>%
  colnames()

# sequences_rel_abund_gr %>%
#   colnames()

# seqs_gr_rel_abund %>%
#   colnames()
# 
# seqs_gr_rel_abund %>%
#   glimpse()

seqs_gr_rel_abund <- seqs_gr_rel_abund %>%
  ungroup() %>%
  mutate(rel_abund_ed = as.numeric(rel_abund))  %>%
  mutate(rel_abund_ed2 = round(rel_abund, 3)) %>%
  arrange(rel_abund_ed) %>%#season.x, treatment, #ordeno per la rel abundance a in situ
  #group_by(treatment, season.x) %>% ##season només 1 mostra per tant no cal que agrupem per estació
  mutate(#rank_abund = as.numeric(order(rel_abund_ed, decreasing = TRUE)), ##general per tot el dataset
         rank_abund_2 = rank(-rel_abund_ed),
         behaviour = case_when(rel_abund_ed < 0.001 ~ 'Rare x < 0.1% ',
                               rel_abund_ed > 0.01 ~ 'Abundant x > 1%',
                               rel_abund_ed <0.01  & rel_abund >0.001 ~ 'Mid  0.1% < x < 1%')) 

seqs_gr_rel_abund$treatment <- seqs_gr_rel_abund$treatment %>% 
  factor(levels=(c("CL", "CD", "PL", "PD", "DL", "VL", 'NA')))
seqs_gr_rel_abund$season.x <- seqs_gr_rel_abund$season.x %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))

# sequences_rel_abund_gr %>%
#   colnames()
# 
# sequences_rel_abund_gr %$%
#   asv_num_col2 %>%
#   unique()

seqs_gr_rel_abund %>%
  filter(rel_abund != (is.na(rel_abund)),
         slope_chosen_days != (is.na(slope_chosen_days))) %$% ##RECUPERO 391 ASV GR and IN SITU 
  relation_miau_remei1 %>%
  unique()
# test <- seqs_gr_rel_abund %>%
#   subset(treatment == 'CD' & season.x == 'Spring')
# seqs_gr_rel_abund$rel_abund_ed2 %>%
#   class()
# seqs_gr_rel_abund$rel_abund_ed == 0.01
# seqs_gr_rel_abund %$% 
#   rank_abund %>%
#   class()
# seqs_gr_rel_abund %>%
#   glimpse()

##table with % of rare reactive ASVs (add to plots as geom_text)-------
##from the total potential bloomers >2 GR & <1% in situ community
seqs_gr_rel_abund %>%
  colnames()

total_asv_match_remei_miau <- seqs_gr_rel_abund %>%
  filter(rel_abund_ed < 0.01 & slope_chosen_days > 2) %>%
  group_by(season.x, treatment) %>%
  summarize(total_common_asv_potential_bloomers = n())

total_asv_match_remei_miau_potential_bloomers <- seqs_gr_rel_abund %>%
  group_by(season.x, treatment) %>%
  summarize(total_common_asv = n()) %>%
  select(-season.x) %>%
  cbind(total_asv_match_remei_miau) %>%
  select(season.x...1, treatment...2, total_common_asv, total_common_asv_potential_bloomers)

colnames(total_asv_match_remei_miau_potential_bloomers) <- c('season.x', 'treatment', 'total_common_asv',
                                                             'total_common_asv_potential_bloomers')

rare_reactive_asv <- seqs_gr_rel_abund %>%
  filter(rel_abund_ed < 0.001 & slope_chosen_days > 2) %>%
  group_by(season.x, treatment) %>%
  summarize(rare_reactive_asv_sum = n())

summarize_match_remei_miau <- total_asv_match_remei_miau_potential_bloomers %>%
  full_join(rare_reactive_asv, by = c('season.x' = 'season.x', 'treatment')) %>%
  mutate(rare_reactive_perc = round(rare_reactive_asv_sum/total_common_asv_potential_bloomers*100, 2))
# as_tibble() %>%
# mutate(rare_reactive_asv_ed = as.numeric(rare_reactive_asv))  %>%
# mutate( ~replace_na(rare_reactive_asv_sum, 0))
summarize_match_remei_miau <- summarize_match_remei_miau %>%
  mutate(rare_reactive_perc_ed = ifelse(is.na(rare_reactive_perc), 0, rare_reactive_perc),
         rare_reactive_asv_sum = ifelse(is.na(rare_reactive_asv_sum), 0, rare_reactive_asv_sum))

growing_asv <- reg_all_slopes_chosen_silva_tax_filt %>%
  group_by(season, treatment) %>%
  summarize(total_growing_asv = n())

miau_seqs_filt_rel_abund_long %>%
  colnames()

insitu_asv <- miau_seqs_filt_rel_abund_long %>%
  filter(rel_abund > 0) %>%
  group_by(season) %>%
  summarize(insitu_present_asv = n())

summarize_insitu_vs_growing <- summarize_match_remei_miau %>%
  left_join(growing_asv, by = c('season.x' = 'season', 'treatment' = 'treatment')) %>%
  left_join(insitu_asv, by = c('season.x' = 'season')) %>%
  select(-rare_reactive_perc)

#write.table(summarize_insitu_vs_growing, 'results/tables/summarize_insitu_vs_growing_ed2.txt', sep='\t')

summarize_insitu_vs_growing %>%
  mutate(perc_common_asv = total_common_asv/insitu_present_asv) %>%
  select(-rare_reactive_perc_ed, - perc_common_asv) %>%
  pivot_longer(cols = c(total_common_asv, rare_reactive_asv_sum, total_growing_asv, insitu_present_asv)) %>%
  ggplot(aes(name, value))+
  geom_col()+
  facet_grid(season.x~treatment)+
  theme_bw()+
  theme(strip.text = element_text(size = 16), axis.text.x = element_text(angle = 0), legend.position = "bottom",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.text = element_text(size = 10),
        legend.title = element_text(size = 12), strip.background = element_blank()) 


summarize_insitu_vs_growing$treatment <- summarize_insitu_vs_growing$treatment %>% 
  factor(levels=(c("CL", "CD", "PL", "PD", "DL", "VL")))

summarize_insitu_vs_growing$season.x <- summarize_insitu_vs_growing$season.x %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))

percentage_shared_asv_insitu_experiments <- summarize_insitu_vs_growing %>%
  mutate(perc_common_asv = total_common_asv/insitu_present_asv,
         perc_rare_reactive_asv = rare_reactive_asv_sum/total_common_asv_potential_bloomers)%>%
  #select(starts_with('perc'), season.x...1, treatment...2) %>%
  select(perc_common_asv, season.x, treatment)%>%
  pivot_longer(cols = starts_with('perc')) %>%
  ggplot(aes(x= season.x, y = value, fill = season.x))+
  geom_bar(stat='identity', width = 1)+
  scale_fill_manual(values = palette_seasons_4)+
  #coord_polar('y', start = 0)+
  scale_y_continuous(labels = percent_format())+
  facet_grid(.~ treatment)+
  labs(fill = 'Season', y = 'shared ASVs between in situ and experiments', x= 'Season')+
  theme_bw()+
  theme(strip.text = element_text(size = 16), legend.position = "top",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.text = element_text(size = 10),
        legend.title = element_text(size = 12), strip.background = element_blank(), axis.title.y = element_text(size = 6),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title = element_blank(), 
        axis.text.y = element_text(size = 6), plot.margin= unit(c(1, 2, 0, 0.2), "lines"), panel.border = element_blank()) 

##detect bloomers: generalist ASV, specific ASV per treatment and seasonal bloomers
seqs_gr_rel_abund %>%
  filter(slope_chosen_days > 2) %>%
  group_by(season.x, asv_num) %>%
  summarize(n = n()) %>%
  filter(n > 18)

seqs_gr_rel_abund %>%
  filter(slope_chosen_days > 2) %>%
  group_by(treatment, asv_num) %>%
  summarize(n = n()) %>%
  filter(n > 15)


##PLOT GR vs rel abund in situ (clustering)-----
#library(ggrepel)
rel_abund_gr_relation_rank <- 
  seqs_gr_rel_abund %>%
  # left_join(summarize_match_remei_miau) %>%
  #filter(pvalue_slope_chosen < 0.05) %>%
  #filter(season.x != is.na(season.x)) %>%
  ggplot(aes(log10(rank_abund_2), slope_chosen_days, color = class.x))+ #, group = behaviour
  geom_text(data = summarize_match_remei_miau, mapping = aes(x = 3.4, y = 2.8, label = rare_reactive_perc_ed), 
            check_overlap = TRUE, na.rm = TRUE, show.legend = FALSE, nudge_x = 0.35, nudge_y = 6, 
            color = 'black', size = 2.5)+
  geom_point(aes(shape =  behaviour), size = 2.5, alpha = 0.8)+
  geom_rect(aes(xmin = 0, xmax = Inf, ymin = 0, ymax = 2), 
             fill = '#D0CFC8', show.legend = FALSE, linejoin = NULL, color = NA, alpha = 0.02)+
  geom_rect(aes(xmin = 0, xmax = log10(rank_abund_2[match(0.01, rel_abund_ed2)]), ymin = 2, ymax = Inf), 
            fill = '#D0CFC8', show.legend = FALSE, linejoin = NULL, color = NA, alpha = 0.02)+
  geom_vline(aes(xintercept = log10(rank_abund_2[match(0.01, rel_abund_ed2)])), linetype = 'dashed')+
  geom_vline(aes(xintercept = log10(rank_abund_2[match(0.001, rel_abund_ed2)])), linetype = 'dashed')+
  # geom_rect(aes(xmin = 0, xmax = log10(rank_abund[match(0.001, rel_abund_ed2)]), ymin = 2.5, ymax = Inf), 
  #           fill = '#D0CFC8', show.legend = FALSE, linejoin = NULL)+
  # geom_rect(aes(xmin = 0, xmax = log10(rank_abund[match(0.001, rel_abund_ed2)]), ymin = 0, ymax = 2.5), 
  #           fill = '#D0CFC8', show.legend = FALSE, linejoin = NULL, color = '#D0CFC8')+
  geom_hline(yintercept = 2, linetype = 'dashed')+
  labs(y = expression("Growth rate day"^"-1"), 
       x = 'Log10(rank abundance)', color = 'Class',
       shape = 'In situ\nrelative\nabundance (%)')+ #expression(paste(italic('In situ'),'\n relative abundance (%)' em queda mal col·locat
  guides(color=guide_legend(ncol = 4, size = 5),
         shape = guide_legend(ncol = 1, size = 3))+
  scale_color_manual(values = palette_class_assigned)+
  scale_x_continuous(expand = c(0, 0), limits = c(0, 4))+
  # geom_ribbon(aes(xmin = aes(log10(rank_abund[match(0.001, rel_abund_ed2)])), ymin = 2.5, 
  #             xmax = rank_abund, ymax = slope_chosen_days))+
  #scale_x_continuous(labels = percent_format())+
  #geom_hline(yintercept = 0.05)+
  # geom_vline(xintercept = log(663))+
  # geom_vline(xintercept = log(185))+
 # geom_vline(xintercept = seqs_gr_rel_abund$rel_abund_ed2 = 0.010)+
  #geom_vline(xintercept = seqs_gr_rel_abund$rank_abund[match(0.001, seqs_gr_rel_abund$rel_abund_ed2)])+
  #geom_area(aes(as.numeric(slope_chosen_days) > 2.5))+
  #geom_polygon(aes(x=0.01, y = 1), fill = 'black')+
  facet_grid(~season.x~treatment)+
  #annotate(geom = 'text', x = 4, y = 5, label = 'Abundant responsive')+#, element_text(size = 2)
  #annotate(geom = 'text', x = 4, y = 0.7, label = 'Abundant unresponsive')+
  theme_bw()+
  theme(strip.text = element_text(size = 16), axis.text.x = element_text(angle = 0), legend.position = "bottom",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.text = element_text(size = 10),
        legend.title = element_text(size = 14), strip.background = element_blank(), #strip.text.x = element_blank(),
        panel.border = element_blank()) 
# 
ggsave(filename = 'rel_abund_gr_cluster100_rank_ed5.pdf', plot = rel_abund_gr_relation_rank, device = NULL, path = 'results/figures/raw/',
        width = 33, height = 30, units = 'cm')

# rel_abund_gr_relation_perc <- seqs_gr_rel_abund %>%
#   #filter(pvalue_slope_chosen < 0.05) %>%
#   #filter(season.x != is.na(season.x)) %>%
#   ggplot(aes(rel_abund, slope_chosen_days, color = class.x))+ #, group = behaviour
#   geom_text(data = summarize_match_remei_miau, mapping = aes(x = 0.006, y = 3, label = rare_reactive_perc_ed),
#             check_overlap = TRUE, na.rm = TRUE, show.legend = FALSE, nudge_x = 0, nudge_y = 3,
#             color = 'black', size = 2)+
#   # geom_rect(aes(xmin = 0, xmax = Inf, ymin = 0, ymax = 2.5),
#   #           fill = '#D0CFC8', show.legend = FALSE, linejoin = NULL, color = '#D3D7DB')+
#   # geom_rect(aes(xmin = 0.01, xmax = Inf, ymin = 2.5, ymax = Inf),
#   #           fill = '#D0CFC8', show.legend = FALSE, linejoin = NULL, color = '#D3D7DB')+
#   geom_vline(xintercept = 0.01)+
#   geom_vline(xintercept = 0.001)+
#   geom_hline(yintercept = 2.5)+
#   geom_point(aes(shape =  behaviour), size = 2, alpha = 0.8)+
#   labs(y = expression("Growth rate day"^"-1"),
#        x = expression(paste('ASV relative abundance',
#                             italic(' in situ '), 'community')), color = 'Class',
#        shape = 'Abundance')+
#   guides(color=guide_legend(ncol = 5),
#          shape = guide_legend(ncol = 1))+
#   scale_color_manual(values = palette_class_assigned)+
#   scale_x_continuous(labels = percent_format(), expand = c(0, 0))+
#   #geom_area()+
#   #geom_polygon(aes(x=0.01, y = 1), fill = 'black')+
#   facet_grid(~season.x~treatment)+
#   #annotate(geom = 'text', x = 4, y = 5, label = 'Abundant responsive')+#, element_text(size = 2)
#   #annotate(geom = 'text', x = 4, y = 0.7, label = 'Abundant unresponsive')+
#   theme_bw()+
#   theme(strip.text = element_text(size = 16), axis.text.x = element_text(angle = 0), legend.position = "bottom",
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.text = element_text(size = 10),
#         legend.title = element_text(size = 12), strip.background = element_blank(), strip.text.x = element_blank(), panel.border = element_blank())
# 
# ggsave(filename = 'rel_abund_gr_cluster100_perc_ed2.pdf', plot = rel_abund_gr_relation_perc, device = NULL, path = 'results/figures/raw/',
#        width = 33, height = 30, units = 'cm')

# ggsave(filename = 'percentage_shared_asv_insitu_experiments.pdf', plot = percentage_shared_asv_insitu_experiments, device = NULL, path = 'results/figures/raw/',
#        width = 33, height = 5, units = 'cm')

# shared_asv_insitu_experiments <- multi_panel_figure(columns = 1, rows = 7,  width = 33, height = 35, 
#                                                   column_spacing = 0.1, unit = 'cm',
#                                                   panel_label_type = 'none')
# 
# shared_asv_insitu_experiments  %<>%
#   fill_panel(percentage_shared_asv_insitu_experiments, column = 1, row = 1) %<>%
#   fill_panel(rel_abund_gr_relation_rank, column = 1, row = 2:7)
# 
# ggsave(filename = 'shared_asv_insitu_experiments_pannel_ed2.pdf', plot = shared_asv_insitu_experiments, device = NULL, path = 'results/figures/raw/',
#        width = 33, height = 35, units = 'cm')

#Same graph but adding a 0 to the growing ASV that don't have a match in situ-------
##Graph without filtering by relative abundance = 0 at in situ time
##SUBSTITUEIXO ELS IS.NA PER ABUNDÀNCIA 0
miau_seqs_filt_rel_abund_long <- miau_seqs_filt_rel_abund_long %>%
  mutate(miau_seqs_filt_rel_abund_long = ifelse(is.na(rel_abund), 0, rel_abund))

##Unim GR data amb in situ rel abund amb les seqüències-------
##unim les  relative abundances a les dades del cluster amb les seqs de MIAU
sequences_rel_abund <- cluster_miau_seqs %>%
  filter(relation_miau_remei1 != is.na(relation_miau_remei1)) %>%
  right_join(miau_seqs_filt_rel_abund_long, by = c('sequence_miau' = 'asv_seqs'))

cluster_miau_seqs %>%
  colnames()
miau_seqs_filt_rel_abund_long %>%
  colnames()
sequences_rel_abund %>%
  colnames()
sequences_rel_abund %>%
  dim()
sequences_rel_abund %$%
  relation_miau_remei1 %>%
  unique() ###1514 asv 

##unim les gr a les dades del cluster amb les seqs de REMEI
sequences_gr <- cluster_remei_seqs %>%
  filter(relation_miau_remei1 != is.na(relation_miau_remei1)) %>%
  right_join(growth_rates_remei_seqs, by = c('sequence_remei' = '.otu')) %>%
  filter(slope_chosen_days != is.na(slope_chosen_days) &
           asv_num_col1 != is.na(asv_num_col1))

sequences_gr %>%
  colnames()

sequences_gr %$%
  asv_num_col1 %>%
  unique() #759 asv

# sequences_rel_abund %$% 
#   sequence
# 
# sequences_gr %$% 
#   sequence
sequences_rel_abund %>%
  head()
#unim les growth rates amb les dades de rel_abundance 
seqs_gr_rel_abund <- sequences_gr %>%
  inner_join(sequences_rel_abund, by = 'relation_miau_remei1') %>%
  filter(#rel_abund != is.na(rel_abund) &
           slope_chosen_days != is.na(slope_chosen_days))

seqs_gr_rel_abund <- seqs_gr_rel_abund %>%  
  mutate(rel_abund = ifelse(is.na(rel_abund), 0, rel_abund))
#seqs_gr_rel_abund1$relation_miau_remei1 
# sequences_rel_abund$relation_miau_remei1 #8208
# 
# seqs_gr_rel_abund %>%
#   colnames()

# seqs_gr_rel_abund2 <- seqs_gr_rel_abund1 %>%
#   full_join(sequences_rel_abund, by = c('relation_miau_remei1' = 'relation_miau_remei1')) %>%
#   filter(rel_abund.x != is.na(rel_abund.x) &
#            slope_chosen_days != is.na(slope_chosen_days))
# 
# seqs_gr_rel_abund %$%
#   asv_num_col1.x %>%
#   unique() ##400 ASV comunes amb gr i abundancia a in situ
# 
# sequences_gr %>%
#   colnames()
# seqs_gr_rel_abund1 %>%
#   colnames()

# seqs_rel_abund_gr %>%
#   dim()
# miau_seqs_filt_rel_abund_long %>%
#   colnames()
# reg_all_slopes_chosen_silva_tax_filt %>%
#   colnames()

# sequences_rel_abund_gr %>%
#   colnames()
# seqs_gr_rel_abund %>%
#   colnames()
# 
# seqs_gr_rel_abund %>%
#   glimpse()

seqs_gr_rel_abund <- seqs_gr_rel_abund %>%
  ungroup() %>%
  mutate(rel_abund_ed = as.numeric(rel_abund))  %>%
  mutate(rel_abund_ed2 = round(rel_abund, 3)) %>%
  arrange(rel_abund_ed) %>%#season.x, treatment, #ordeno per la rel abundance a in situ
  #group_by(treatment, season.x) %>% ##season només 1 mostra per tant no cal que agrupem per estació
  mutate(#rank_abund = as.numeric(order(rel_abund_ed, decreasing = FALSE)), ##general per tot el dataset
    rank_abund_2 = rank(-rel_abund_ed),
    behaviour = case_when(rel_abund_ed < 0.001 ~ 'Rare x < 0.1% ',
                          rel_abund_ed > 0.01 ~ 'Abundant x > 1%',
                          rel_abund_ed <0.01  & rel_abund >0.001 ~ 'Mid  0.1% < x < 1%')) 

seqs_gr_rel_abund$treatment <- seqs_gr_rel_abund$treatment %>% 
  factor(levels=(c("CL", "CD", "PL", "PD", "DL", "VL", 'NA')))
seqs_gr_rel_abund$season.x <- seqs_gr_rel_abund$season.x %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))

# sequences_rel_abund_gr %>%
#   colnames()
# 
# sequences_rel_abund_gr %$%
#   asv_num_col2 %>%
#   unique()

seqs_gr_rel_abund %>%
  filter(#rel_abund != (is.na(rel_abund)),
    slope_chosen_days != (is.na(slope_chosen_days))) %$% ##RECUPERO 392 ASV GR and IN SITU 
  relation_miau_remei1 %>%
  unique()
# test <- seqs_gr_rel_abund %>%
#   subset(treatment == 'CD' & season.x == 'Spring')
# seqs_gr_rel_abund$rel_abund_ed2 %>%
#   class()
# seqs_gr_rel_abund$rel_abund_ed == 0.01
# seqs_gr_rel_abund %$% 
#   rank_abund %>%
#   class()
# seqs_gr_rel_abund %>%
#   glimpse()

##table with % of rare reactive ASVs (add to plots as geom_text)-------
##from the total potential bloomers <2 GR & <0.1% community
seqs_gr_rel_abund %>%
  colnames()

total_asv_match_remei_miau <- seqs_gr_rel_abund %>%
  filter(rel_abund_ed < 0.01 & slope_chosen_days > 2) %>%
  group_by(season.x, treatment) %>%
  summarize(total_common_asv_potential_bloomers = n())

total_asv_match_remei_miau_potential_bloomers <- seqs_gr_rel_abund %>%
  group_by(season.x, treatment) %>%
  summarize(total_common_asv = n()) %>%
  select(-season.x) %>%
  cbind(total_asv_match_remei_miau) %>%
  select(season.x...1, treatment...2, total_common_asv, total_common_asv_potential_bloomers)

colnames(total_asv_match_remei_miau_potential_bloomers) <- c('season.x', 'treatment', 'total_common_asv',
                                                             'total_common_asv_potential_bloomers')

rare_reactive_asv <- seqs_gr_rel_abund %>%
  filter(rel_abund_ed < 0.001 & slope_chosen_days > 2) %>%
  group_by(season.x, treatment) %>%
  summarize(rare_reactive_asv_sum = n())

summarize_match_remei_miau <- total_asv_match_remei_miau_potential_bloomers %>%
  full_join(rare_reactive_asv, by = c('season.x' = 'season.x', 'treatment')) %>%
  mutate(rare_reactive_perc = round(rare_reactive_asv_sum/total_common_asv_potential_bloomers*100, 2))
# as_tibble() %>%
# mutate(rare_reactive_asv_ed = as.numeric(rare_reactive_asv))  %>%
# mutate( ~replace_na(rare_reactive_asv_sum, 0))
summarize_match_remei_miau <- summarize_match_remei_miau %>%
  mutate(rare_reactive_perc_ed = ifelse(is.na(rare_reactive_perc), 0, rare_reactive_perc),
         rare_reactive_asv_sum = ifelse(is.na(rare_reactive_asv_sum), 0, rare_reactive_asv_sum))

growing_asv <- reg_all_slopes_chosen_silva_tax_filt %>%
  group_by(season, treatment) %>%
  summarize(total_growing_asv = n())

miau_seqs_filt_rel_abund_long %>%
  colnames()

insitu_asv <- miau_seqs_filt_rel_abund_long %>%
  filter(rel_abund > 0) %>%
  group_by(season) %>%
  summarize(insitu_present_asv = n())

summarize_insitu_vs_growing <- summarize_match_remei_miau %>%
  left_join(growing_asv, by = c('season.x' = 'season', 'treatment' = 'treatment')) %>%
  left_join(insitu_asv, by = c('season.x' = 'season')) %>%
  select(-rare_reactive_perc)

#write.table(summarize_insitu_vs_growing, 'results/tables/summarize_insitu_vs_growing_ed2.txt', sep='\t')

summarize_insitu_vs_growing %>%
  mutate(perc_common_asv = total_common_asv/insitu_present_asv) %>%
  select(-rare_reactive_perc_ed, - perc_common_asv) %>%
  pivot_longer(cols = c(total_common_asv, rare_reactive_asv_sum, total_growing_asv, insitu_present_asv)) %>%
  ggplot(aes(name, value))+
  geom_col()+
  facet_grid(season.x~treatment)+
  theme_bw()+
  theme(strip.text = element_text(size = 16), axis.text.x = element_text(angle = 0), legend.position = "bottom",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.text = element_text(size = 10),
        legend.title = element_text(size = 12), strip.background = element_blank()) 


summarize_insitu_vs_growing$treatment <- summarize_insitu_vs_growing$treatment %>% 
  factor(levels=(c("CL", "CD", "PL", "PD", "DL", "VL")))

summarize_insitu_vs_growing$season.x <- summarize_insitu_vs_growing$season.x %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))

percentage_shared_asv_insitu_experiments <- summarize_insitu_vs_growing %>%
  mutate(perc_common_asv = total_common_asv/insitu_present_asv,
         perc_rare_reactive_asv = rare_reactive_asv_sum/total_common_asv_potential_bloomers)%>%
  #select(starts_with('perc'), season.x...1, treatment...2) %>%
  select(perc_common_asv, season.x, treatment)%>%
  pivot_longer(cols = starts_with('perc')) %>%
  ggplot(aes(x= season.x, y = value, fill = season.x))+
  geom_bar(stat='identity', width = 1)+
  scale_fill_manual(values = palette_seasons_4)+
  #coord_polar('y', start = 0)+
  scale_y_continuous(labels = percent_format())+
  facet_grid(.~ treatment)+
  labs(fill = 'Season', y = 'shared ASVs between in situ and experiments', x= 'Season')+
  theme_bw()+
  theme(strip.text = element_text(size = 16), legend.position = "top",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.text = element_text(size = 10),
        legend.title = element_text(size = 12), strip.background = element_blank(), axis.title.y = element_text(size = 6),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title = element_blank(), 
        axis.text.y = element_text(size = 6), plot.margin= unit(c(1, 2, 0, 0.2), "lines"), panel.border = element_blank()) 

##detect bloomers: generalist ASV, specific ASV per treatment and seasonal bloomers
seqs_gr_rel_abund %>%
  filter(slope_chosen_days > 2) %>%
  group_by(season.x, asv_num) %>%
  summarize(n = n()) %>%
  filter(n > 18)

seqs_gr_rel_abund %>%
  filter(slope_chosen_days > 2) %>%
  group_by(treatment, asv_num) %>%
  summarize(n = n()) %>%
  filter(n > 15)

rel_abund_gr_relation_rank <- 
  seqs_gr_rel_abund %>%
  #left_join(summarize_match_remei_miau) %>%
  #filter(pvalue_slope_chosen < 0.05) %>%
  #filter(season.x != is.na(season.x)) %>%
  ggplot(aes(log10(rank_abund_2), slope_chosen_days, color = class.x))+ #, group = behaviour
  geom_point(aes(shape =  behaviour), size = 2.5, alpha = 0.8)+
  geom_rect(aes(xmin = 0, xmax = Inf, ymin = 0, ymax = 2), 
            fill = '#D0CFC8', show.legend = FALSE, linejoin = NULL, color = NA, alpha = 0.02)+
  geom_rect(aes(xmin = 0, xmax = log10(rank_abund_2[match(0.01, rel_abund_ed2)]), ymin = 2, ymax = Inf), 
            fill = '#D0CFC8', show.legend = FALSE, linejoin = NULL, color = NA, alpha = 0.02)+
  geom_vline(aes(xintercept = log10(rank_abund_2[match(0.01, rel_abund_ed2)])), linetype = 'dashed')+
  geom_vline(aes(xintercept = log10(rank_abund_2[match(0.001, rel_abund_ed2)])), linetype = 'dashed')+
  # geom_rect(aes(xmin = 0, xmax = log10(rank_abund[match(0.001, rel_abund_ed2)]), ymin = 2.5, ymax = Inf), 
  #           fill = '#D0CFC8', show.legend = FALSE, linejoin = NULL)+
  # geom_rect(aes(xmin = 0, xmax = log10(rank_abund[match(0.001, rel_abund_ed2)]), ymin = 0, ymax = 2.5), 
  #           fill = '#D0CFC8', show.legend = FALSE, linejoin = NULL, color = '#D0CFC8')+
  geom_hline(yintercept = 2, linetype = 'dashed')+
  labs(y = expression("Growth rate day"^"-1"), 
       x = 'Log10(rank abundance)', color = 'Class',
       shape = 'In situ\nrelative\nabundance (%)')+ #expression(paste(italic('In situ'),'\n relative abundance (%)' em queda mal col·locat
  guides(color=guide_legend(ncol = 4, size = 5),
         shape = guide_legend(ncol = 1, size = 3))+
  scale_color_manual(values = palette_class_assigned)+
  scale_x_continuous(expand = c(0, 0), limits = c(0, 4))+
  # geom_ribbon(aes(xmin = aes(log10(rank_abund[match(0.001, rel_abund_ed2)])), ymin = 2.5, 
  #             xmax = rank_abund, ymax = slope_chosen_days))+
  #scale_x_continuous(labels = percent_format())+
  #geom_hline(yintercept = 0.05)+
  # geom_vline(xintercept = log(663))+
  # geom_vline(xintercept = log(185))+
  # geom_vline(xintercept = seqs_gr_rel_abund$rel_abund_ed2 = 0.010)+
  #geom_vline(xintercept = seqs_gr_rel_abund$rank_abund[match(0.001, seqs_gr_rel_abund$rel_abund_ed2)])+
  #geom_area(aes(as.numeric(slope_chosen_days) > 2.5))+
  #geom_polygon(aes(x=0.01, y = 1), fill = 'black')+
  facet_grid(~season.x~treatment)+
  geom_text(data = summarize_match_remei_miau, mapping = aes(x = 3.4, y = 2.8, label = rare_reactive_perc_ed), 
            check_overlap = TRUE, na.rm = TRUE, show.legend = FALSE, nudge_x = 0.35, nudge_y = 6, 
            color = 'black', size = 2.5)+
  #annotate(geom = 'text', x = 4, y = 5, label = 'Abundant responsive')+#, element_text(size = 2)
  #annotate(geom = 'text', x = 4, y = 0.7, label = 'Abundant unresponsive')+
  theme_bw()+
  theme(strip.text = element_text(size = 16), axis.text.x = element_text(angle = 0), legend.position = "bottom",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.text = element_text(size = 10),
        legend.title = element_text(size = 14), strip.background = element_blank(), #strip.text.x = element_blank(),
        panel.border = element_blank()) 

ggsave(filename = 'rel_abund_gr_cluster100_rank_0rel_abund_ed5.pdf', plot = rel_abund_gr_relation_rank, device = NULL, path = 'results/figures/raw/',
       width = 33, height = 30, units = 'cm')

# rel_abund_gr_relation_perc <- seqs_gr_rel_abund %>%
#   #filter(pvalue_slope_chosen < 0.05) %>%
#   #filter(season.x != is.na(season.x)) %>%
#   ggplot(aes(rel_abund, slope_chosen_days, color = class.x))+ #, group = behaviour
#   geom_text(data = summarize_match_remei_miau, mapping = aes(x = 0.006, y = 3, label = rare_reactive_perc_ed), 
#             check_overlap = TRUE, na.rm = TRUE, show.legend = FALSE, nudge_x = 0, nudge_y = 3, 
#             color = 'black', size = 2)+
#   # geom_rect(aes(xmin = 0, xmax = Inf, ymin = 0, ymax = 2.5), 
#   #           fill = '#D0CFC8', show.legend = FALSE, linejoin = NULL, color = '#D3D7DB')+
#   # geom_rect(aes(xmin = 0.01, xmax = Inf, ymin = 2.5, ymax = Inf), 
#   #           fill = '#D0CFC8', show.legend = FALSE, linejoin = NULL, color = '#D3D7DB')+
#   geom_vline(xintercept = 0.01)+
#   geom_vline(xintercept = 0.001)+
#   geom_hline(yintercept = 2.5)+
#   geom_point(aes(shape =  behaviour), size = 2, alpha = 0.8)+
#   labs(y = expression("Growth rate day"^"-1"), 
#        x = expression(paste('ASV relative abundance',
#                             italic(' in situ '), 'community')), color = 'Class',
#        shape = 'Abundance')+
#   guides(color=guide_legend(ncol = 5),
#          shape = guide_legend(ncol = 1))+
#   scale_color_manual(values = palette_class_assigned)+
#   scale_x_continuous(labels = percent_format(), expand = c(0, 0))+
#   #geom_area()+
#   #geom_polygon(aes(x=0.01, y = 1), fill = 'black')+
#   facet_grid(~season.x~treatment)+
#   #annotate(geom = 'text', x = 4, y = 5, label = 'Abundant responsive')+#, element_text(size = 2)
#   #annotate(geom = 'text', x = 4, y = 0.7, label = 'Abundant unresponsive')+
#   theme_bw()+
#   theme(strip.text = element_text(size = 16), axis.text.x = element_text(angle = 0), legend.position = "bottom",
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.text = element_text(size = 10),
#         legend.title = element_text(size = 12), strip.background = element_blank(), strip.text.x = element_blank(),
#         panel.border = element_blank())  
# 
# ggsave(filename = 'rel_abund_gr_cluster100_perc_0rel_abund_ed2.pdf', plot = rel_abund_gr_relation_perc, device = NULL, path = 'results/figures/raw/',
#        width = 33, height = 30, units = 'cm')
# 
# rel_abund_gr_relation_log10 <- 
# seqs_gr_rel_abund %>%
#   # left_join(summarize_match_remei_miau) %>%
#   #filter(pvalue_slope_chosen < 0.05) %>%
#   #filter(season.x != is.na(season.x)) %>%
#   ggplot(aes(log10(rel_abund), slope_chosen_days, color = class.x))+ #, group = behaviour
#   geom_text(data = summarize_match_remei_miau, mapping = aes(x = -4, y = 4, label = rare_reactive_perc_ed), 
#             check_overlap = TRUE, na.rm = TRUE, show.legend = FALSE, nudge_x = 0, nudge_y = 5, 
#             color = 'black', size = 4)+
#   # geom_rect(aes(xmin = -Inf, xmax = 0, ymin = 0, ymax = 2.5), 
#   #           fill = '#D0CFC8', show.legend = FALSE, linejoin = NULL, color = '#D3D7DB')+
#   # geom_rect(aes(xmax = log10(rel_abund[match(0.001, rel_abund_ed2)]), xmin = -Inf, ymin = 2.5, ymax = Inf), 
#   #           fill = '#D0CFC8', show.legend = FALSE, linejoin = NULL, color = '#D3D7DB')+
#   geom_vline(aes(xintercept = log10(0.01)))+
#   geom_vline(aes(xintercept = log10(0.001)))+
#   # geom_rect(aes(xmin = 0, xmax = log10(rank_abund[match(0.001, rel_abund_ed2)]), ymin = 2.5, ymax = Inf), 
#   #           fill = '#D0CFC8', show.legend = FALSE, linejoin = NULL)+
#   # geom_rect(aes(xmin = 0, xmax = log10(rank_abund[match(0.001, rel_abund_ed2)]), ymin = 0, ymax = 2.5), 
#   #           fill = '#D0CFC8', show.legend = FALSE, linejoin = NULL, color = '#D0CFC8')+
#   geom_hline(yintercept = 2.5)+
#   geom_point(aes(shape =  behaviour), size = 2.5, alpha = 0.7)+
#   labs(y = expression("Growth rate day"^"-1"), 
#        x = 'Log10(relative abundance)', color = 'Class',
#        shape = 'In situ\nrelative\nabundance (%)')+ #expression(paste(italic('In situ'),'\n relative abundance (%)' em queda mal col·locat
#   guides(color=guide_legend(ncol = 4, size = 5),
#          shape = guide_legend(ncol = 1, size = 3))+
#   scale_color_manual(values = palette_class_assigned)+
#   #scale_x_continuous(expand = c(0, 0), limits = c(0, 3.4))+
#   # geom_ribbon(aes(xmin = aes(log10(rank_abund[match(0.001, rel_abund_ed2)])), ymin = 2.5, 
#   #             xmax = rank_abund, ymax = slope_chosen_days))+
#   #scale_x_continuous(labels = percent_format())+
#   #geom_hline(yintercept = 0.05)+
#   # geom_vline(xintercept = log(663))+
#   # geom_vline(xintercept = log(185))+
#   # geom_vline(xintercept = seqs_gr_rel_abund$rel_abund_ed2 = 0.010)+
#   #geom_vline(xintercept = seqs_gr_rel_abund$rank_abund[match(0.001, seqs_gr_rel_abund$rel_abund_ed2)])+
#   #geom_area(aes(as.numeric(slope_chosen_days) > 2.5))+
#   #geom_polygon(aes(x=0.01, y = 1), fill = 'black')+
#   facet_grid(~season.x~treatment)+
#   #annotate(geom = 'text', x = 4, y = 5, label = 'Abundant responsive')+#, element_text(size = 2)
#   #annotate(geom = 'text', x = 4, y = 0.7, label = 'Abundant unresponsive')+
#   theme_bw()+
#   theme(strip.text = element_text(size = 16), axis.text.x = element_text(angle = 0), legend.position = "bottom",
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.text = element_text(size = 10),
#         legend.title = element_text(size = 14), strip.background = element_blank(), strip.text.x = element_blank(),
#         panel.border = element_blank()) 
# 
# ggsave(filename = 'rel_abund_gr_cluster100_log10relabund_0rel_abund_ed2.pdf', plot = rel_abund_gr_relation_log10, device = NULL, path = 'results/figures/raw/',
#        width = 33, height = 30, units = 'cm')
# 
# # ggsave(filename = 'percentage_shared_asv_insitu_experiments.pdf', plot = percentage_shared_asv_insitu_experiments, device = NULL, path = 'results/figures/raw/',
# #        width = 33, height = 5, units = 'cm')
# 
# shared_asv_insitu_experiments <- multi_panel_figure(columns = 1, rows = 7,  width = 33, height = 35, 
#                                                     column_spacing = 0.1, unit = 'cm',
#                                                     panel_label_type = 'none')
# 
# shared_asv_insitu_experiments  %<>%
#   fill_panel(percentage_shared_asv_insitu_experiments, column = 1, row = 1) %<>%
#   fill_panel(rel_abund_gr_relation_rank, column = 1, row = 2:7)
# 
# ggsave(filename = 'shared_asv_insitu_experiments_pannel_0rel_abund_ed2.pdf', plot = shared_asv_insitu_experiments, device = NULL, path = 'results/figures/raw/',
#        width = 33, height = 35, units = 'cm')
# 
# 
#SAME GRAPH BUT AT 99% OF IDENTITY??------
