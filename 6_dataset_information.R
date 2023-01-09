
library(gt) #per fer les taules
library(webshot2) #per guardar les taules
library(markdown)

####LINE 722 DATASET REMIAU-----
####REMEI DATASET (WITHOUT IN SITU DATA)
##Data import------
rem_relabun_melt <-  read.table("data/rem_relabun_melt.txt", sep="\t", header = TRUE) %>%
  filter(season != "Early_fall")
reg_all_slopes_chosen_silva_tax <- read.csv("data/intermediate_files/reg_all_slopes_chosen_silva_tax.csv", sep=",") %>%
  filter(season != "Early_fall")

rem_relabun_melt 

%>%
  colnames()

rem_relabun_melt %$%
  asv_num %>%
  unique() #4594 asv al meu dataset

reg_all_slopes_chosen_silva_tax$treatment <- reg_all_slopes_chosen_silva_tax$treatment %>% 
  factor(levels=(c("CL", "CD", "PL", "PD", "DL", "VL")))
reg_all_slopes_chosen_silva_tax$season <- reg_all_slopes_chosen_silva_tax$season %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))

##Information of my dataset ------
rem_relabun_melt_1perc <- rem_relabun_melt_1perc %>%
  mutate(replicate = as.factor(replicate))

abundance_1perc_asv <- rem_relabun_melt_1perc %>%
  group_by(treatment, season, replicate, hours.x) %>%
  summarize(abundance_per_sample = sum(Abundance))

abundance_1perc_asv %>%
  ungroup() %>%
  summarize(mean(abundance_per_sample), #0.832
            min(abundance_per_sample), #0.614
            max(abundance_per_sample)) #0.944

abundance_1perc_asv %>%
  #group_by(treatment, season, replicate, hours.x) %>%
  #summarize(Abundance) %>%
  ggplot(aes(x=hours.x, y=abundance_per_sample))+
  geom_col()+
  geom_hline(yintercept =0.5)+
  facet_grid(.~treatment~season~replicate)+
  theme_bw()

reg_all_slopes_chosen_silva_tax_1perc %$%
  asv_num %>%
  unique() #167

reg_all_slopes_chosen_silva_tax %$%
  asv_num %>%
  unique() #4594

167/4594 #3.64% of total ASV

#Quantes GR hem pogut calcular de tot el dataset?
no_nas_sig <- reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05) %>%
  summarise(total_non_na = sum(!is.na(slope_chosen_days)))

nas <- reg_all_slopes_chosen_silva_tax %>%
  summarise(total_na = sum(is.na(slope_chosen_days)))

total <- reg_all_slopes_chosen_silva_tax %>%
  dim()
no_nas_sig/total #5.018% de GR calculades

nas/total#82.95% no calculades

#D'aquestes que arriben almenys a l'1% de la comunitat en algun experiment?
no_nas_sig <- reg_all_slopes_chosen_silva_tax_1perc %>%
  filter(pvalue_slope_chosen < 0.05) %>%
  summarise(total_non_na = sum(!is.na(slope_chosen_days)))

nas <- reg_all_slopes_chosen_silva_tax_1perc %>%
  summarise(total_na = sum(is.na(slope_chosen_days)))

total <- reg_all_slopes_chosen_silva_tax_1perc %>%
  dim()
no_nas_sig/total #52.44% de GR calculades però 28.87% significatives

nas/total#47.56% no calculades

#Amb quantes ASV treballo quan selecciono només les ASV que arriben a l'1% de la comunitat en algun moment?
reg_all_slopes_chosen_silva_tax %>%
  group_by(asv_num) %>%
  count(asv_num) #4594

reg_all_slopes_chosen_silva_tax_1perc %>%
  group_by(asv_num) %>%
  count(asv_num) #167

167/4594 # 0.0379 = 3.635%

##Quin % de reads representen aquestes ASV de l'1%?
reg_all_slopes_chosen_silva_tax %>%
  colnames()

rem_relabun_melt %>%
  colnames()

reg_all_slopes_chosen_silva_tax_1perc %T>% 
  colnames() %>%
  dim()

percentage_rel_abund_1perc_dataset <- rem_relabun_melt %>%
  group_by(season, treatment, asv_num) %>%
  mutate(mean_rel_abund = mean(Abundance)) %>%
  distinct(mean_rel_abund) %>%
  right_join(reg_all_slopes_chosen_silva_tax_1perc, keep = TRUE, by = c('treatment', 'asv_num', 'season')) %>%
  distinct(mean_rel_abund, treatment.x, season.x, asv_num.x) %>%
  group_by(treatment.x, season.x) %>%
  mutate(abundance_1perc = sum(mean_rel_abund)) %>%
  distinct(abundance_1perc) %>%
  as_tibble()

#create a table
library(gt)
library(gtExtras)
library(webshot2)

percentage_rel_abund_1perc_dataset$treatment.x <- percentage_rel_abund_1perc_dataset$treatment.x %>% 
  factor(levels=(c("CL", "CD", "PL", "PD", "DL", "VL")))
percentage_rel_abund_1perc_dataset$season.x <- percentage_rel_abund_1perc_dataset$season.x %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))

table_percentage_rel_abund_1perc <- percentage_rel_abund_1perc_dataset %>%
  arrange(season.x, treatment.x, abundance_1perc) %>%
  mutate(abundance_1perc = abundance_1perc*100) %>%
  gt() %>%
  fmt_number(columns = abundance_1perc, decimals = 1) %>%
  cols_label(season.x = 'Season', treatment.x = 'Treatment', abundance_1perc = 'Abundance represented by all the ASV \n that represent \n >1%') %>%
  cols_width('treatment.x' ~ px(80),
             'season.x' ~ px(80),
    everything() ~ px(200))
  
table_percentage_rel_abund_1perc %>%
  class()

table_percentage_rel_abund_1perc %>%
  gtsave(filename = 'results/tables/table_percentage_rel_abund_1perc.png', expand = 3)


#How many asv are present in different treatments and we could calculate different growth rates?
reg_all_slopes_chosen_silva_tax_1perc  %>%
  colnames()
reg_all_slopes_chosen_silva_tax_1perc %>%
  filter(!is.na(slope_chosen_days)) %>%
  filter(pvalue_slope_chosen > 0.05) %>%
  group_by(asv_num) %>%
  count(asv_num) %>%
  as_tibble() %>%
  filter(n > 2) #139

reg_all_slopes_chosen_silva_tax %>%
  filter(!is.na(slope_chosen_days)) %>%
  filter(pvalue_slope_chosen > 0.05) %>%
  group_by(asv_num) %>%
  count(asv_num) %>%
  as_tibble() %>%
  filter(n > 2) ##2234


# pivot_wider(names_from = 'asv_num', values_from = 'slope_chosen_days', 
#             id_cols = c('treatment','season','domain' , 'phylum', 'class', 
#                         'order', 'family', 'genus', 'species')) %>%
# summarize(!is.na(20:187))

##Deviation of GR between treatments for each ASV
reg_all_slopes_chosen_silva_tax_1perc %>%
  group_by(treatment, season, asv_num) %>%
  summarize(deviation = sd(slope_chosen_days)) %>%
  filter(!is.na(deviation)) %>%
  ggplot(aes(asv_num, deviation))+
  geom_point()+
  theme_bw()

#Percentage of unclassified at different levels and by different taxonomical groups
##whole dataset
rem_relabun_melt %>%
  group_by(asv_num) %>%
  summarise(is_na = sum(is.na(genus))) %>%
  filter(is_na > 0) %>%
  nrow()
rem_relabun_melt %>%
  colnames()
rem_relabun_melt %$%
  asv_num %>%
  unique()
1819/4594
#species 4594 100% unclassified
#genus 1819 39.59% unclassified
#family 0 0% unclassified
#order 0
#class 0

#only 1perc group
rem_relabun_melt_1perc %>%
  group_by(asv_num) %>%
  summarise(is_na = sum(is.na(species))) %>%
  filter(is_na > 0) %>%
  nrow()

59/167
#species 167  (all) 100% unclassified
#genus 56 35.33% unclassified
#family 0
#order 0
#class 0

reg_all_slopes_chosen_silva_tax_1perc %>% group_by()

#Shared asv between treatments and seasons?
#General barplots distributions

#Summarizing general dataset
reg_all_slopes_chosen_silva_tax %>%
  ungroup() %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  group_by(phylum) %>%
  count(slope_chosen_days)
dplyr::summarize(max_growht = max(slope_chosen_days), 
                 num_asv = n()) 

reg_all_slopes_chosen_silva_tax %>%
  ungroup() %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  group_by(phylum) %>%
  tally(slope_chosen_days)

test<- reg_all_slopes_chosen_silva_tax %>%
  ungroup() %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  group_by(phylum) %>%
  add_tally(slope_chosen_days)

reg_all_slopes_chosen_silva_tax %>%
  group_by(treatment, phylum) %>%
  summarize(num_asv = n())

reg_all_slopes_chosen_silva_tax %>%
  colnames()

reg_all_slopes_chosen_silva_tax %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  group_by(treatment, season, phylum, slope_chosen_days) %>%
  summarize(n = n())

#max gr of the dataset for each phylum
max_gr_table_phylum <- reg_all_slopes_chosen_silva_tax %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  group_by(treatment, phylum) %>%
  summarize(n = max(slope_chosen_days)) %>%
  pivot_wider(values_from = n, names_from = treatment)

max_gr_table_class <- reg_all_slopes_chosen_silva_tax %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  group_by(treatment, phylum, class) %>%
  summarize(n = max(slope_chosen_days)) %>%
  pivot_wider(values_from = n, names_from = treatment)

max_gr_table_order <- reg_all_slopes_chosen_silva_tax %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  group_by(treatment, phylum, class, order) %>%
  summarize(n = max(slope_chosen_days)) %>%
  pivot_wider(values_from = n, names_from = treatment)

max_gr_table_family <- reg_all_slopes_chosen_silva_tax %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  group_by(treatment, phylum, class, order, family) %>%
  summarize(n = max(slope_chosen_days)) %>%
  pivot_wider(values_from = n, names_from = treatment)

max_gr_table_asv <- reg_all_slopes_chosen_silva_tax %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  group_by(treatment, phylum, class, order, family, asv_num) %>%
  summarize(n = max(slope_chosen_days)) %>%
  pivot_wider(values_from = n, names_from = treatment)

reg_all_slopes_chosen_silva_tax %>% 
  subset(asv_num == 'asv81') %>%
  ggplot(aes(treatment, slope_chosen_days))+
  geom_point()+
  theme_bw()

reg_all_slopes_chosen_silva_tax %>%
  colnames()

reg_all_slopes_chosen_silva_tax %>%
  subset(family == 'Alteromonadaceae') %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  ggplot(aes(treatment, slope_chosen_days))+
  geom_point(alpha = 0.7, position = position_jitter(0.25), aes(color = season))+
  geom_violin(alpha = 0.2, draw_quantiles = c(0.25, 0.75))+
  geom_line()+
  labs(x = "Treatment", y = expression("Growth rate day"^"-1"), color = 'Season')+
  scale_color_manual(values = palette_seasons_4)+
  #facet_wrap(.~genus, scales = "free")+
  theme_bw()+
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

reg_all_slopes_chosen_silva_tax %>%
  subset(family == 'Alteromonadaceae') %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  ggplot(aes(treatment, slope_chosen_days))+
  geom_point(alpha = 0.7, position = position_jitter(0.25), aes(color = season))+
  geom_violin(alpha = 0.2, draw_quantiles = c(0.25, 0.75))+
  geom_line()+
  scale_y_continuous(limits = c(0,11))+
  labs(x = "Treatment", y = expression("Growth rate day"^"-1"), color = 'Season')+
  scale_color_manual(values = palette_seasons_4)+
  #facet_wrap(.~genus, scales = "free")+
  theme_bw()+
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())


#distributions SAR11 clade i Alteromonadaceae comprovar que no siguin més altres a SAR i perquè.
reg_all_slopes_chosen_silva_tax %>%
  subset(order == 'SAR11 clade' |
           family == 'Alteromonadaceae') %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  ggplot(aes(family, slope_chosen_days))+
  geom_point(alpha = 0.7, position = position_jitter(0.25), aes(color = season))+
  geom_violin(alpha = 0.2, draw_quantiles = c(0.25, 0.75))+
  geom_line()+
  scale_y_continuous(limits = c(0,11))+
  labs(x = "Family", y = expression("Growth rate day"^"-1"), color = 'Season')+
  scale_color_manual(values = palette_seasons_4)+
  #facet_wrap(.~genus, scales = "free")+
  theme_bw()+
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())


sar_altero <- reg_all_slopes_chosen_silva_tax %>%
  subset(order == 'SAR11 clade' |
           family == 'Alteromonadaceae') %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  ggplot(aes(slope_chosen_days, family))+
  #geom_point(alpha = 0.7, position = position_jitter(0.25), aes(color = season))+
  #geom_violin(alpha = 0.2, draw_quantiles = c(0.25, 0.75))+
  #geom_line()+
  scale_x_continuous(limits = c(0,11))+
  geom_density_ridges()+
  labs(y = "Family", x = expression("Growth rate day"^"-1"), color = 'Season')+
  scale_color_manual(values = palette_seasons_4)+
  #facet_wrap(.~genus, scales = "free")+
  theme_bw()+
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

sar_alter2<- reg_all_slopes_chosen_silva_tax %>%
  subset(order == 'SAR11 clade' |
           family == 'Alteromonadaceae') %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  ggplot(aes(y = slope_chosen_days, x = family))+
  geom_point(alpha = 0.7, position = position_jitter(0.25), aes(color = season))+
  geom_violin(alpha = 0.2, draw_quantiles = c(0.25, 0.75))+
  geom_line()+
  scale_y_continuous(limits = c(0,11))+
  #geom_density_ridges()+
  labs(x = "Family", y = expression("Growth rate day"^"-1"), color = 'Season')+
  scale_color_manual(values = palette_seasons_4)+
  #facet_wrap(.~genus, scales = "free")+
  theme_bw()+
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

sar11_altero <- grid.arrange(sar_altero, sar_alter2)

##Supplementary tables-------
library(gtsummary)
##phylums
table_mean_gr_phylum <- reg_all_slopes_chosen_silva_tax %>%
  filter(slope_chosen_days != is.na(slope_chosen_days)) %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  group_by(phylum) %>%
  mutate(counts = n()) %>%
  group_by(phylum) %>%
  summarize(mean_growth_rate = mean(slope_chosen_days),
            sd_growth_rate = sd(slope_chosen_days),
            counts = counts) %>%
  distinct(phylum, mean_growth_rate, sd_growth_rate, counts) %>%
  arrange(-mean_growth_rate)  %>%
  filter(counts > 1) %>%
  ungroup() %>%
  gt() %>%
  fmt_number(columns = contains('growth_rate'), decimals = 1) %>%
  cols_label(mean_growth_rate = 'Mean growth rate/day', sd_growth_rate = 'sd growth rate', phylum = 'Phylum', 
             counts = 'Number of significant gr calculated') %>%
  cols_width(everything() ~ px(130))

table_mean_gr_phylum <- 
  tableGrob(table_mean_gr_phylum, 
    theme = theme_transparent()

library(glue)
save(file = table_mean_gr_phylum, 'table_mean_gr_phylum.pdf')
pdf(table_mean_gr_phylum, file = 'table_mean_gr_phylum.pdf', height = 10, width = 8)
table_mean_gr_phylum <- grid.table(table_mean_gr_phylum)
dev.off()

table_mean_gr_phylum %>%
  ggsave(filename = 'table_mean_gr_phylum.png', expand = 3, path = 'results/tables/', dpi = 300)

table_percentage_rel_abund_1perc %>%
  gtsave(filename = 'results/tables/table_percentage_rel_abund_1perc.png', expand = 3)
table_mean_gr_phylum %>%
  grid.table() %>%
save(table_mean_gr_phylum, file = 'table_mean_gr_phylum.pdf')
##class
reg_all_slopes_chosen_silva_tax %>%
  filter(slope_chosen_days != is.na(slope_chosen_days)) %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  group_by(class) %>%
  mutate(counts = n()) %>%
  group_by(phylum_f, class) %>%
  summarize(mean_growth_rate = mean(slope_chosen_days),
            sd_growth_rate = sd(slope_chosen_days),
            counts = counts) %>%
  distinct(class, mean_growth_rate, sd_growth_rate, counts) %>%
  filter(counts > 1) %>%
  ungroup() %>%
  arrange(-mean_growth_rate)  %>%
  gt() %>%
  fmt_number(columns = contains('growth_rate'), decimals = 1) %>%
  cols_label(mean_growth_rate = 'Mean growth rate/day', sd_growth_rate = 'sd growth rate', phylum_f = 'Phylum', 
             class = 'Class',
             counts = 'Number of significant gr calculated') %>%
  cols_width(everything() ~ px(130))


##order
reg_all_slopes_chosen_silva_tax %>%
  filter(slope_chosen_days != is.na(slope_chosen_days)) %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  group_by(order) %>%
  mutate(counts = n()) %>%
  group_by(phylum_f, class_f, order) %>%
  summarize(mean_growth_rate = mean(slope_chosen_days),
            sd_growth_rate = sd(slope_chosen_days),
            counts = counts) %>%
  distinct(phylum_f, class_f, order, mean_growth_rate, sd_growth_rate, counts) %>%
  filter(counts > 1) %>%
  ungroup() %>%
  arrange(-mean_growth_rate)  %>%
  gt() %>%
  fmt_number(columns = contains('growth_rate'), decimals = 1) %>%
  cols_label(mean_growth_rate = 'Mean growth rate/day', sd_growth_rate = 'sd growth rate', phylum_f = 'Phylum', 
             class_f = 'Class', order = 'Order',
             counts = 'Number of significant gr calculated') %>%
  cols_width(everything() ~ px(130))


#family
reg_all_slopes_chosen_silva_tax %>%
  filter(slope_chosen_days != is.na(slope_chosen_days)) %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  group_by(family) %>%
  mutate(counts = n()) %>%
  group_by(phylum_f, class_f, order_f, family) %>%
  summarize(mean_growth_rate = mean(slope_chosen_days),
            sd_growth_rate = sd(slope_chosen_days),
            counts = counts) %>%
  distinct(family, mean_growth_rate, sd_growth_rate, counts) %>%
  filter(counts > 1) %>%
  ungroup() %>%
  arrange(-mean_growth_rate)  %>%
  gt() %>%
  fmt_number(columns = contains('growth_rate'), decimals = 1) %>%
  cols_label(mean_growth_rate = 'Mean growth rate/day', sd_growth_rate = 'sd growth rate', phylum_f = 'Phylum', 
             class_f = 'Class', order_f = 'Order', family = 'Family',
             counts = 'Number of significant gr calculated') %>%
  cols_width(everything() ~ px(130))


#asv_num
table_mean_gr_asv <-  
  reg_all_slopes_chosen_silva_tax %>%
  filter(slope_chosen_days != is.na(slope_chosen_days)) %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  group_by(asv_num) %>%
  mutate(counts = n()) %>%
  group_by(phylum_f, class_f, order_f, family_f, asv_num) %>%
  summarize(mean_growth_rate = mean(slope_chosen_days),
            sd_growth_rate = sd(slope_chosen_days),
            counts_growth_rate = counts) %>%
  distinct(asv_num, mean_growth_rate, sd_growth_rate, counts_growth_rate) %>%
  filter(counts_growth_rate > 1) %>%
  ungroup() %>%
  arrange(-mean_growth_rate)  %>%
  gt() %>%
  fmt_number(columns = contains('growth_rate'), decimals = 1) %>%
  cols_label(mean_growth_rate = 'Mean growth rate/day', sd_growth_rate = 'sd growth rate', phylum_f = 'Phylum', 
             class_f = md('Class'), order_f = 'Order', family_f = 'Family', asv_num = 'ASV number',
             counts_growth_rate = 'Nº of significant growth rates calculated') %>%
  cols_align(align = "left", columns = contains('_f')) %>%
  cols_align(align = "center", columns = contains('growth_rate')) %>%
  cols_width(contains('growth_rate') ~ px(80),
             everything() ~ px(180))
  
# table_mean_gr_asv %>%
#     glimpse()
# setwd ("~/Documentos/Doctorat/REMEI/")

# write.table(table_mean_gr_asv, 'results/tables/table_mean_gr_asv.txt', sep = '\t')

table_mean_gr_asv %>%
  gt::gtsave(file = 'table_mean_gr_asv.png', path = 'results/tables/') ##error attempt to apply non-function. 
##Error due to not specifing a mathematiccal operator or a function betwen the two values.

library(kableExtra)

##mean per family
table_mean_gr_family <-  
  reg_all_slopes_chosen_silva_tax %>%
  filter(slope_chosen_days != is.na(slope_chosen_days)) %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  group_by(family) %>%
  mutate(counts = n()) %>%
  group_by(phylum_f, class_f, order_f, family_f) %>%
  summarize(mean_growth_rate = mean(slope_chosen_days),
            sd_growth_rate = sd(slope_chosen_days),
            counts_growth_rate = counts) %>%
  distinct(family_f, mean_growth_rate, sd_growth_rate, counts_growth_rate) %>%
  filter(counts_growth_rate > 1) %>%
  ungroup() %>%
  arrange(-mean_growth_rate)  %>%
  gt() %>%
  fmt_number(columns = contains('growth_rate'), decimals = 1) %>%
  cols_label(mean_growth_rate = 'Mean growth rate/day', sd_growth_rate = 'sd growth rate', phylum_f = 'Phylum', 
             class_f = md('Class'), order_f = 'Order', family_f = 'Family',
             counts_growth_rate = 'Nº of significant growth rates calculated') %>%
  cols_align(align = "left", columns = contains('_f')) %>%
  cols_align(align = "center", columns = contains('growth_rate')) %>%
  cols_width(contains('growth_rate') ~ px(80),
             everything() ~ px(180))

# table_mean_gr_asv %>%
#     glimpse()
# setwd ("~/Documentos/Doctorat/REMEI/")

# write.table(table_mean_gr_asv, 'results/tables/table_mean_gr_asv.txt', sep = '\t')

table_mean_gr_family %>%
  gt::gtsave(file = 'table_mean_gr_family.png', path = 'results/tables/') 


##mean per order
table_mean_gr_order <-  
  reg_all_slopes_chosen_silva_tax %>%
  filter(slope_chosen_days != is.na(slope_chosen_days)) %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  group_by(order) %>%
  mutate(counts = n()) %>%
  group_by(phylum_f, class_f, order_f) %>%
  summarize(mean_growth_rate = mean(slope_chosen_days),
            sd_growth_rate = sd(slope_chosen_days),
            counts_growth_rate = counts) %>%
  distinct(order_f, mean_growth_rate, sd_growth_rate, counts_growth_rate) %>%
  filter(counts_growth_rate > 1) %>%
  ungroup() %>%
  arrange(-mean_growth_rate)  %>%
  gt() %>%
  fmt_number(columns = contains('growth_rate'), decimals = 1) %>%
  cols_label(mean_growth_rate = 'Mean growth rate/day', sd_growth_rate = 'sd growth rate', phylum_f = 'Phylum', 
             class_f = md('Class'), order_f = 'Order',
             counts_growth_rate = 'Nº of significant growth rates calculated') %>%
  cols_align(align = "left", columns = contains('_f')) %>%
  cols_align(align = "center", columns = contains('growth_rate')) %>%
  cols_width(contains('growth_rate') ~ px(90),
             everything() ~ px(200))

# table_mean_gr_asv %>%
#     glimpse()
# setwd ("~/Documentos/Doctorat/REMEI/")

# write.table(table_mean_gr_asv, 'results/tables/table_mean_gr_asv.txt', sep = '\t')

table_mean_gr_order %>%
  gt::gtsave(file = 'table_mean_gr_order.png', path = 'results/tables/') 

##tables filtered by those groups that at least have 12 GR calculated-----------------
#asv_num
table_mean_gr_asv <-  
  reg_all_slopes_chosen_silva_tax %>%
  filter(slope_chosen_days != is.na(slope_chosen_days)) %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  group_by(asv_num) %>%
  mutate(counts = n()) %>%
  group_by(phylum_f, class_f, order_f, family_f, asv_num) %>%
  summarize(mean_growth_rate = mean(slope_chosen_days),
            sd_growth_rate = sd(slope_chosen_days),
            counts_growth_rate = counts) %>%
  distinct(asv_num, mean_growth_rate, sd_growth_rate, counts_growth_rate) %>%
  filter(counts_growth_rate > 12) %>%
  ungroup() %>%
  arrange(-mean_growth_rate)  %>%
  gt() %>%
  fmt_number(columns = contains('growth_rate'), decimals = 1) %>%
  cols_label(mean_growth_rate = 'Mean growth rate/day', sd_growth_rate = 'sd growth rate', phylum_f = 'Phylum', 
             class_f = 'Class', order_f = 'Order', family_f = 'Family', asv_num = 'ASV number',
             counts_growth_rate = 'Nº of significant growth rates calculated') %>%
  cols_align(align = "left", columns = contains('_f')) %>%
  cols_align(align = "center", columns = contains('growth_rate')) )
# %>% si corro aquest últim punt no em funciona la creació de la taula.
#   cols_width(contains('growth_rate') ~ px(80),

#              everything() ~ px(200)) 

# table_mean_gr_asv %>%
#     glimpse()
# setwd ("~/Documentos/Doctorat/REMEI/")

#write.table(table_mean_gr_asv, 'results/tables/table_mean_gr_asv.txt', sep = '\t')

table_mean_gr_asv %>%
  gt::gtsave(file = 'table_mean_gr_asv_filt.png', path = 'results/tables/') ##error attempt to apply non-function. 
##Error due to not specifing a mathematiccal operator or a function betwen the two values.

library(kableExtra)

##mean per family
table_mean_gr_family <-  
  reg_all_slopes_chosen_silva_tax %>%
  filter(slope_chosen_days != is.na(slope_chosen_days)) %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  group_by(family) %>%
  mutate(counts = n()) %>%
  group_by(phylum_f, class_f, order_f, family_f) %>%
  summarize(mean_growth_rate = mean(slope_chosen_days),
            sd_growth_rate = sd(slope_chosen_days),
            counts_growth_rate = counts) %>%
  distinct(family_f, mean_growth_rate, sd_growth_rate, counts_growth_rate) %>%
  filter(counts_growth_rate > 12) %>%
  ungroup() %>%
  arrange(-mean_growth_rate)  %>%
  gt() %>%
  fmt_number(columns = contains('growth_rate'), decimals = 1) %>%
  cols_label(mean_growth_rate = 'Mean growth rate/day', sd_growth_rate = 'sd growth rate', phylum_f = 'Phylum', 
             class_f = md('Class'), order_f = 'Order', family_f = 'Family',
             counts_growth_rate = 'Nº of significant growth rates calculated') %>%
  cols_align(align = "left", columns = contains('_f')) %>%
  cols_align(align = "center", columns = contains('growth_rate')) %>%
  cols_width(contains('growth_rate') ~ px(80),
             everything() ~ px(180))

# table_mean_gr_asv %>%
#     glimpse()
# setwd ("~/Documentos/Doctorat/REMEI/")

# write.table(table_mean_gr_asv, 'results/tables/table_mean_gr_asv.txt', sep = '\t')

table_mean_gr_family %>%
  gt::gtsave(file = 'table_mean_gr_family_filt.png', path = 'results/tables/') 


##mean per order
table_mean_gr_order <-  
  reg_all_slopes_chosen_silva_tax %>%
  filter(slope_chosen_days != is.na(slope_chosen_days)) %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  group_by(order) %>%
  mutate(counts = n()) %>%
  group_by(phylum_f, class_f, order_f) %>%
  summarize(mean_growth_rate = mean(slope_chosen_days),
            sd_growth_rate = sd(slope_chosen_days),
            counts_growth_rate = counts) %>%
  distinct(order_f, mean_growth_rate, sd_growth_rate, counts_growth_rate) %>%
  filter(counts_growth_rate > 12) %>%
  ungroup() %>%
  arrange(-mean_growth_rate)  %>%
  gt() %>%
  fmt_number(columns = contains('growth_rate'), decimals = 1) %>%
  cols_label(mean_growth_rate = 'Mean growth rate/day', sd_growth_rate = 'sd growth rate', phylum_f = 'Phylum', 
             class_f = md('Class'), order_f = 'Order',
             counts_growth_rate = 'Nº of significant growth rates calculated') %>%
  cols_align(align = "left", columns = contains('_f')) %>%
  cols_align(align = "center", columns = contains('growth_rate')) %>%
  cols_width(contains('growth_rate') ~ px(90),
             everything() ~ px(200))

# table_mean_gr_asv %>%
#     glimpse()
# setwd ("~/Documentos/Doctorat/REMEI/")

# write.table(table_mean_gr_asv, 'results/tables/table_mean_gr_asv.txt', sep = '\t')

table_mean_gr_order %>%
  gt::gtsave(file = 'table_mean_gr_order_filt.png', path = 'results/tables/')


######REMIAU DATA SET (WITH IN SITU DATA)------
##Data import-------
rem_relabun_melt <-  read.table("data/remiau_relabun_melt.txt", sep="\t", header = TRUE) %>%
  filter(season != "Early_fall")
reg_all_slopes_chosen_silva_tax <- read.csv("data/intermediate_files/remiau_reg_all_slopes_chosen_silva_tax.csv", sep=",") %>%
  filter(season != "Early_fall")

rem_relabun_melt %$%
  asv_num %>%
  unique() #5111 asv al meu dataset

reg_all_slopes_chosen_silva_tax$treatment <- reg_all_slopes_chosen_silva_tax$treatment %>% 
  factor(levels=(c("CL", "CD", "PL", "PD", "DL", "VL")))
reg_all_slopes_chosen_silva_tax$season <- reg_all_slopes_chosen_silva_tax$season %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))

##Information of my dataset ------
rem_relabun_melt_1perc <- rem_relabun_melt %>%
  right_join(abundant_asv, by = "asv_num", copy = TRUE) #per comprovar amb quin % de relabund estic treballant

rem_relabun_melt_1perc <- rem_relabun_melt_1perc %>%
  mutate(replicate = as.factor(replicate))

abundance_1perc_asv <- rem_relabun_melt_1perc %>%
  group_by(treatment, season, replicate, hours.x) %>%
  summarize(abundance_per_sample = sum(Abundance))

abundance_1perc_asv %>%
  ungroup() %>%
  summarize(mean(abundance_per_sample), #0.832
            min(abundance_per_sample), #0.617
            max(abundance_per_sample)) #0.944

abundance_1perc_asv %>%
  #group_by(treatment, season, replicate, hours.x) %>%
  #summarize(Abundance) %>%
  ggplot(aes(x=hours.x, y=abundance_per_sample))+
  geom_col()+
  geom_hline(yintercept =0.5)+
  facet_grid(.~treatment~season~replicate)+
  theme_bw()

reg_all_slopes_chosen_silva_tax_1perc %$%
  asv_num %>%
  unique() #173

reg_all_slopes_chosen_silva_tax %$%
  asv_num %>%
  unique() #5111 

173/5111 #3.38% of total ASV

#Quantes GR hem pogut calcular de tot el dataset?
no_nas_sig <- reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05) %>%
  summarise(total_non_na = sum(!is.na(slope_chosen_days)))

nas <- reg_all_slopes_chosen_silva_tax %>%
  summarise(total_na = sum(is.na(slope_chosen_days)))

total <- reg_all_slopes_chosen_silva_tax %>%
  dim()
no_nas_sig/total #4.66% de GR calculades i significatives

nas/total#75.51% no calculades

(total-nas-no_nas_sig)/total ##19.83% calculades però no significatives.

#D'aquestes que arriben a almenys a l'1% de la comunitat en algun experiment?
no_nas_sig <- reg_all_slopes_chosen_silva_tax_1perc %>%
  filter(pvalue_slope_chosen < 0.05) %>%
  summarise(total_non_na = sum(!is.na(slope_chosen_days)))

nas <- reg_all_slopes_chosen_silva_tax_1perc %>%
  summarise(total_na = sum(is.na(slope_chosen_days)))

total <- reg_all_slopes_chosen_silva_tax_1perc %>%
  dim()
no_nas_sig/total #52.47% de GR calculades però 27.63% significatives
(total-nas-no_nas_sig)/total ##24.84 calculades i no significatives.
nas/total#47.63% no calculades

#Amb quantes ASV treballo quan selecciono només les ASV que arriben a l'1% de la comunitat en algun moment?
reg_all_slopes_chosen_silva_tax %>%
  group_by(asv_num) %>%
  count(asv_num) #5111

reg_all_slopes_chosen_silva_tax_1perc %>%
  group_by(asv_num) %>%
  count(asv_num) #173

173/5111 #3.38% of total ASV

##Quin % de reads representen aquestes ASV de l'1%?
reg_all_slopes_chosen_silva_tax %>%
  colnames()

rem_relabun_melt %>%
  colnames()

reg_all_slopes_chosen_silva_tax_1perc %T>% 
  colnames() %>%
  dim()

dataset_mean_abundance <- rem_relabun_melt %>%
  group_by(season, treatment, asv_num) %>%
  mutate(mean_rel_abund = mean(Abundance),
         min = min(Abundance),
         max = max(Abundance))
  # summarize(mean_abundance = mean(Abundance), #0.832
  #           min(Abundance), #0.614
  #           max(Abundance)) %>%
dataset_mean_abundance %>%
  distinct(mean_rel_abund) ##no funciona, perquè vull la mean per cada treatment/season

percentage_rel_abund_1perc_dataset <- rem_relabun_melt %>%
  group_by(season, treatment, asv_num) %>%
  mutate(mean_rel_abund = mean(Abundance)) %>%
  distinct(mean_rel_abund) %>%
  right_join(reg_all_slopes_chosen_silva_tax_1perc, keep = TRUE, by = c('treatment', 'asv_num', 'season')) %>%
  distinct(mean_rel_abund, treatment.x, season.x, asv_num.x) %>%
  group_by(treatment.x, season.x) %>%
  mutate(abundance_1perc = sum(mean_rel_abund)) %>%
  distinct(abundance_1perc) %>%
  as_tibble()

#create a table
library(gt)
#library(gtExtras) #no tinc el package instalat
library(webshot2)

percentage_rel_abund_1perc_dataset$treatment.x <- percentage_rel_abund_1perc_dataset$treatment.x %>% 
  factor(levels=(c("CL", "CD", "PL", "PD", "DL", "VL")))
percentage_rel_abund_1perc_dataset$season.x <- percentage_rel_abund_1perc_dataset$season.x %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))

table_percentage_rel_abund_1perc <- percentage_rel_abund_1perc_dataset %>%
  arrange(season.x, treatment.x, abundance_1perc) %>%
  mutate(abundance_1perc = abundance_1perc*100) %>%
  gt() %>%
  fmt_number(columns = abundance_1perc, decimals = 1) %>%
  cols_label(season.x = 'Season', treatment.x = 'Treatment', abundance_1perc = 
               'Abundance represented by all the ASV \n that represent \n >1%') %>%
  cols_width('treatment.x' ~ px(80),
             'season.x' ~ px(80),
             everything() ~ px(200))

table_percentage_rel_abund_1perc %>%
  class()

table_percentage_rel_abund_1perc %>%
  gtsave(filename = 'results/tables/remiau_table_percentage_rel_abund_1perc.png', expand = 3)

#How many asv are present in different treatments and we could calculate different growth rates?
reg_all_slopes_chosen_silva_tax_1perc  %>%
  colnames()
reg_all_slopes_chosen_silva_tax_1perc %>%
  filter(!is.na(slope_chosen_days)) %>%
  filter(pvalue_slope_chosen > 0.05) %>%
  group_by(asv_num) %>%
  count(asv_num) %>%
  as_tibble() %>%
  filter(n > 2) #138

reg_all_slopes_chosen_silva_tax %>%
  filter(!is.na(slope_chosen_days)) %>%
  filter(pvalue_slope_chosen > 0.05) %>%
  group_by(asv_num) %>%
  count(asv_num) %>%
  as_tibble() %>%
  filter(n > 2) ##2268


# pivot_wider(names_from = 'asv_num', values_from = 'slope_chosen_days', 
#             id_cols = c('treatment','season','domain' , 'phylum', 'class', 
#                         'order', 'family', 'genus', 'species')) %>%
# summarize(!is.na(20:187))

##Deviation of GR between treatments for each ASV
reg_all_slopes_chosen_silva_tax_1perc %>%
  group_by(treatment, season, asv_num) %>%
  summarize(deviation = sd(slope_chosen_days)) %>%
  filter(!is.na(deviation)) %>%
  ggplot(aes(asv_num, deviation))+
  geom_point()+
  theme_bw()

#Percentage of unclassified at different levels and by different taxonomical groups
##whole dataset
rem_relabun_melt %>%
  group_by(asv_num) %>%
  summarise(is_na = sum(is.na(species))) %>%
  filter(is_na > 0) %>%
  nrow()
rem_relabun_melt %>%
  colnames()
rem_relabun_melt %$%
  asv_num %>%
  unique()
2018/5111
#species 5111 100% unclassified
#genus 1819 39.48% unclassified
#family 0 0% unclassified
#order 0
#class 0

#only 1perc group
rem_relabun_melt_1perc %>%
  group_by(asv_num) %>%
  summarise(is_na = sum(is.na(genus))) %>%
  filter(is_na > 0) %>%
  nrow()

57/173
#species 167  (all) 100% unclassified
#genus 56 32.94% unclassified
#family 0
#order 0
#class 0

reg_all_slopes_chosen_silva_tax_1perc %>% 
  group_by()

#Shared asv between treatments and seasons?
#General barplots distributions

#Summarizing general dataset
reg_all_slopes_chosen_silva_tax %>%
  ungroup() %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  group_by(phylum_f, class_f) %>%
  count(slope_chosen_days) %>%
  dplyr::summarize(max_growht = max(slope_chosen_days), 
                 num_asv = n(),
                 mean_gr = mean(slope_chosen_days),
                 sd_gr = sd(slope_chosen_days)) 

reg_all_slopes_chosen_silva_tax %>%
  ungroup() %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  group_by(phylum) %>%
  tally(slope_chosen_days)

test <- reg_all_slopes_chosen_silva_tax %>%
  ungroup() %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  group_by(phylum) %>%
  add_tally(slope_chosen_days)

reg_all_slopes_chosen_silva_tax %>%
  group_by(treatment, phylum) %>%
  summarize(num_asv = n())

reg_all_slopes_chosen_silva_tax %>%
  colnames()

reg_all_slopes_chosen_silva_tax %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  group_by(treatment, season, phylum, slope_chosen_days) %>%
  summarize(n = n())

#max gr of the dataset for each phylum
max_gr_table_phylum <- reg_all_slopes_chosen_silva_tax %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  group_by(treatment, phylum) %>%
  summarize(n = max(slope_chosen_days)) %>%
  pivot_wider(values_from = n, names_from = treatment)

max_gr_table_class <- reg_all_slopes_chosen_silva_tax %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  group_by(treatment, phylum, class) %>%
  summarize(n = max(slope_chosen_days)) %>%
  pivot_wider(values_from = n, names_from = treatment)

max_gr_table_order <- reg_all_slopes_chosen_silva_tax %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  group_by(treatment, phylum, class, order) %>%
  summarize(n = max(slope_chosen_days)) %>%
  pivot_wider(values_from = n, names_from = treatment)

max_gr_table_family <- reg_all_slopes_chosen_silva_tax %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  group_by(treatment, phylum, class, order, family) %>%
  summarize(n = max(slope_chosen_days)) %>%
  pivot_wider(values_from = n, names_from = treatment)

max_gr_table_asv <- reg_all_slopes_chosen_silva_tax %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  group_by(treatment, phylum, class, order, family, asv_num) %>%
  summarize(n = max(slope_chosen_days)) %>%
  pivot_wider(values_from = n, names_from = treatment)

reg_all_slopes_chosen_silva_tax %>% 
  subset(asv_num == 'asv81') %>%
  ggplot(aes(treatment, slope_chosen_days))+
  geom_point()+
  theme_bw()

reg_all_slopes_chosen_silva_tax %>%
  colnames()

reg_all_slopes_chosen_silva_tax %>%
  subset(family == 'Alteromonadaceae') %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  ggplot(aes(treatment, slope_chosen_days))+
  geom_point(alpha = 0.7, position = position_jitter(0.25), aes(color = season))+
  geom_violin(alpha = 0.2, draw_quantiles = c(0.25, 0.75))+
  geom_line()+
  labs(x = "Treatment", y = expression("Growth rate day"^"-1"), color = 'Season')+
  scale_color_manual(values = palette_seasons_4)+
  #facet_wrap(.~genus, scales = "free")+
  theme_bw()+
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

reg_all_slopes_chosen_silva_tax %>%
  subset(family == 'Alteromonadaceae') %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  ggplot(aes(treatment, slope_chosen_days))+
  geom_point(alpha = 0.7, position = position_jitter(0.25), aes(color = season))+
  geom_violin(alpha = 0.2, draw_quantiles = c(0.25, 0.75))+
  geom_line()+
  scale_y_continuous(limits = c(0,11))+
  labs(x = "Treatment", y = expression("Growth rate day"^"-1"), color = 'Season')+
  scale_color_manual(values = palette_seasons_4)+
  #facet_wrap(.~genus, scales = "free")+
  theme_bw()+
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())


#distributions SAR11 clade i Alteromonadaceae comprovar que no siguin més altres a SAR i perquè.
reg_all_slopes_chosen_silva_tax %>%
  subset(order == 'SAR11 clade' |
           family == 'Alteromonadaceae') %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  ggplot(aes(family, slope_chosen_days))+
  geom_point(alpha = 0.7, position = position_jitter(0.25), aes(color = season))+
  geom_violin(alpha = 0.2, draw_quantiles = c(0.25, 0.75))+
  geom_line()+
  scale_y_continuous(limits = c(0,11))+
  labs(x = "Family", y = expression("Growth rate day"^"-1"), color = 'Season')+
  scale_color_manual(values = palette_seasons_4)+
  #facet_wrap(.~genus, scales = "free")+
  theme_bw()+
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

library(ggridges)
library(gridExtra)
sar_altero <- reg_all_slopes_chosen_silva_tax %>%
  subset(order == 'SAR11 clade' |
           family == 'Alteromonadaceae') %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  ggplot(aes(slope_chosen_days, family))+
  #geom_point(alpha = 0.7, position = position_jitter(0.25), aes(color = season))+
  #geom_violin(alpha = 0.2, draw_quantiles = c(0.25, 0.75))+
  #geom_line()+
  scale_x_continuous(limits = c(0,11))+
  geom_density_ridges()+
  labs(y = "Family", x = expression("Growth rate day"^"-1"), color = 'Season')+
  scale_color_manual(values = palette_seasons_4)+
  #facet_wrap(.~genus, scales = "free")+
  theme_bw()+
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

sar_alter2<- reg_all_slopes_chosen_silva_tax %>%
  subset(order == 'SAR11 clade' |
           family == 'Alteromonadaceae') %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  ggplot(aes(y = slope_chosen_days, x = family))+
  geom_point(alpha = 0.7, position = position_jitter(0.25), aes(color = season))+
  geom_violin(alpha = 0.2, draw_quantiles = c(0.25, 0.75))+
  geom_line()+
  scale_y_continuous(limits = c(0,11))+
  #geom_density_ridges()+
  labs(x = "Family", y = expression("Growth rate day"^"-1"), color = 'Season')+
  scale_color_manual(values = palette_seasons_4)+
  #facet_wrap(.~genus, scales = "free")+
  theme_bw()+
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

sar11_altero <- grid.arrange(sar_altero, sar_alter2)

##Supplementary tables (no funciona,  revisar-ho)-------
library(gtsummary)
##phylums
table_mean_gr_phylum <- reg_all_slopes_chosen_silva_tax %>%
  filter(slope_chosen_days != is.na(slope_chosen_days)) %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  group_by(phylum) %>%
  mutate(counts = n()) %>%
  group_by(phylum) %>%
  summarize(mean_growth_rate = mean(slope_chosen_days),
            sd_growth_rate = sd(slope_chosen_days),
            counts = counts) %>%
  distinct(phylum, mean_growth_rate, sd_growth_rate, counts) %>%
  arrange(-mean_growth_rate)  %>%
  filter(counts > 1) %>%
  ungroup() %>%
  gt() %>%
  fmt_number(columns = contains('growth_rate'), decimals = 1) %>%
  cols_label(mean_growth_rate = 'Mean growth rate/day', sd_growth_rate = 'sd growth rate', phylum = 'Phylum', 
             counts = 'Number of significant gr calculated') %>%
  cols_width(everything() ~ px(130))

table_mean_gr_phylum <- 
  tableGrob(table_mean_gr_phylum, 
            theme = theme_transparent())
            
library(glue)
#save(file = table_mean_gr_phylum, 'remiau_table_mean_gr_phylum.pdf')
pdf(table_mean_gr_phylum, file = 'remiau_table_mean_gr_phylum.pdf', height = 10, width = 8)
table_mean_gr_phylum <- grid.table(table_mean_gr_phylum)
dev.off()

table_mean_gr_phylum %>%
  ggsave(filename = 'remiau_table_mean_gr_phylum.png', expand = 3, path = 'results/tables/', dpi = 300)

table_percentage_rel_abund_1perc %>%
  gtsave(filename = 'results/tables/table_percentage_rel_abund_1perc.png', expand = 3)
table_mean_gr_phylum %>%
  grid.table() %>%
  save(table_mean_gr_phylum, file = 'table_mean_gr_phylum.pdf')
##class
reg_all_slopes_chosen_silva_tax %>%
  filter(slope_chosen_days != is.na(slope_chosen_days)) %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  group_by(class) %>%
  mutate(counts = n()) %>%
  group_by(phylum_f, class) %>%
  summarize(mean_growth_rate = mean(slope_chosen_days),
            sd_growth_rate = sd(slope_chosen_days),
            counts = counts) %>%
  distinct(class, mean_growth_rate, sd_growth_rate, counts) %>%
  filter(counts > 1) %>%
  ungroup() %>%
  arrange(-mean_growth_rate)  %>%
  gt() %>%
  fmt_number(columns = contains('growth_rate'), decimals = 1) %>%
  cols_label(mean_growth_rate = 'Mean growth rate/day', sd_growth_rate = 'sd growth rate', phylum_f = 'Phylum', 
             class = 'Class',
             counts = 'Number of significant gr calculated') %>%
  cols_width(everything() ~ px(130))


##order
reg_all_slopes_chosen_silva_tax %>%
  filter(slope_chosen_days != is.na(slope_chosen_days)) %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  group_by(order) %>%
  mutate(counts = n()) %>%
  group_by(phylum_f, class_f, order) %>%
  summarize(mean_growth_rate = mean(slope_chosen_days),
            sd_growth_rate = sd(slope_chosen_days),
            counts = counts) %>%
  distinct(phylum_f, class_f, order, mean_growth_rate, sd_growth_rate, counts) %>%
  filter(counts > 1) %>%
  ungroup() %>%
  arrange(-mean_growth_rate)  %>%
  gt() %>%
  fmt_number(columns = contains('growth_rate'), decimals = 1) %>%
  cols_label(mean_growth_rate = 'Mean growth rate/day', sd_growth_rate = 'sd growth rate', phylum_f = 'Phylum', 
             class_f = 'Class', order = 'Order',
             counts = 'Number of significant gr calculated') %>%
  cols_width(everything() ~ px(130))


#family
reg_all_slopes_chosen_silva_tax %>%
  filter(slope_chosen_days != is.na(slope_chosen_days)) %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  group_by(family) %>%
  mutate(counts = n()) %>%
  group_by(phylum_f, class_f, order_f, family) %>%
  summarize(mean_growth_rate = mean(slope_chosen_days),
            sd_growth_rate = sd(slope_chosen_days),
            counts = counts) %>%
  distinct(family, mean_growth_rate, sd_growth_rate, counts) %>%
  filter(counts > 1) %>%
  ungroup() %>%
  arrange(-mean_growth_rate)  %>%
  gt() %>%
  fmt_number(columns = contains('growth_rate'), decimals = 1) %>%
  cols_label(mean_growth_rate = 'Mean growth rate/day', sd_growth_rate = 'sd growth rate', phylum_f = 'Phylum', 
             class_f = 'Class', order_f = 'Order', family = 'Family',
             counts = 'Number of significant gr calculated') %>%
  cols_width(everything() ~ px(130))


#asv_num
table_mean_gr_asv <-  
  reg_all_slopes_chosen_silva_tax %>%
  filter(slope_chosen_days != is.na(slope_chosen_days)) %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  group_by(asv_num) %>%
  mutate(counts = n()) %>%
  group_by(phylum_f, class_f, order_f, family_f, asv_num) %>%
  summarize(mean_growth_rate = mean(slope_chosen_days),
            sd_growth_rate = sd(slope_chosen_days),
            counts_growth_rate = counts) %>%
  distinct(asv_num, mean_growth_rate, sd_growth_rate, counts_growth_rate) %>%
  filter(counts_growth_rate > 1) %>%
  ungroup() %>%
  arrange(-mean_growth_rate)  %>%
  gt() %>%
  fmt_number(columns = contains('growth_rate'), decimals = 1) %>%
  cols_label(mean_growth_rate = 'Mean growth rate/day', sd_growth_rate = 'sd growth rate', phylum_f = 'Phylum', 
             class_f = md('Class'), order_f = 'Order', family_f = 'Family', asv_num = 'ASV number',
             counts_growth_rate = 'Nº of significant growth rates calculated') %>%
  cols_align(align = "left", columns = contains('_f')) %>%
  cols_align(align = "center", columns = contains('growth_rate')) %>%
  cols_width(contains('growth_rate') ~ px(80),
             everything() ~ px(180))

# table_mean_gr_asv %>%
#     glimpse()
# setwd ("~/Documentos/Doctorat/REMEI/")

# write.table(table_mean_gr_asv, 'results/tables/table_mean_gr_asv.txt', sep = '\t')

table_mean_gr_asv %>%
  gt::gtsave(file = 'table_mean_gr_asv.png', path = 'results/tables/') ##error attempt to apply non-function. 
##Error due to not specifing a mathematiccal operator or a function betwen the two values.

library(kableExtra)

##mean per family
table_mean_gr_family <-  
  reg_all_slopes_chosen_silva_tax %>%
  filter(slope_chosen_days != is.na(slope_chosen_days)) %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  group_by(family) %>%
  mutate(counts = n()) %>%
  group_by(phylum_f, class_f, order_f, family_f) %>%
  summarize(mean_growth_rate = mean(slope_chosen_days),
            sd_growth_rate = sd(slope_chosen_days),
            counts_growth_rate = counts) %>%
  distinct(family_f, mean_growth_rate, sd_growth_rate, counts_growth_rate) %>%
  filter(counts_growth_rate > 1) %>%
  ungroup() %>%
  arrange(-mean_growth_rate)  %>%
  gt() %>%
  fmt_number(columns = contains('growth_rate'), decimals = 1) %>%
  cols_label(mean_growth_rate = 'Mean growth rate/day', sd_growth_rate = 'sd growth rate', phylum_f = 'Phylum', 
             class_f = md('Class'), order_f = 'Order', family_f = 'Family',
             counts_growth_rate = 'Nº of significant growth rates calculated') %>%
  cols_align(align = "left", columns = contains('_f')) %>%
  cols_align(align = "center", columns = contains('growth_rate')) %>%
  cols_width(contains('growth_rate') ~ px(80),
             everything() ~ px(180))

# table_mean_gr_asv %>%
#     glimpse()
# setwd ("~/Documentos/Doctorat/REMEI/")

# write.table(table_mean_gr_asv, 'results/tables/table_mean_gr_asv.txt', sep = '\t')

table_mean_gr_family %>%
  gt::gtsave(file = 'table_mean_gr_family.png', path = 'results/tables/') 


##mean per order
table_mean_gr_order <-  
  reg_all_slopes_chosen_silva_tax %>%
  filter(slope_chosen_days != is.na(slope_chosen_days)) %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  group_by(order) %>%
  mutate(counts = n()) %>%
  group_by(phylum_f, class_f, order_f) %>%
  summarize(mean_growth_rate = mean(slope_chosen_days),
            sd_growth_rate = sd(slope_chosen_days),
            counts_growth_rate = counts) %>%
  distinct(order_f, mean_growth_rate, sd_growth_rate, counts_growth_rate) %>%
  filter(counts_growth_rate > 1) %>%
  ungroup() %>%
  arrange(-mean_growth_rate)  %>%
  gt() %>%
  fmt_number(columns = contains('growth_rate'), decimals = 1) %>%
  cols_label(mean_growth_rate = 'Mean growth rate/day', sd_growth_rate = 'sd growth rate', phylum_f = 'Phylum', 
             class_f = md('Class'), order_f = 'Order',
             counts_growth_rate = 'Nº of significant growth rates calculated') %>%
  cols_align(align = "left", columns = contains('_f')) %>%
  cols_align(align = "center", columns = contains('growth_rate')) %>%
  cols_width(contains('growth_rate') ~ px(90),
             everything() ~ px(200))

# table_mean_gr_asv %>%
#     glimpse()
# setwd ("~/Documentos/Doctorat/REMEI/")

# write.table(table_mean_gr_asv, 'results/tables/table_mean_gr_asv.txt', sep = '\t')

table_mean_gr_order %>%
  gt::gtsave(file = 'table_mean_gr_order.png', path = 'results/tables/') 

##tables filtered by those groups that at least have 12 GR calculated-----------------
#asv_num
table_mean_gr_asv <-  
  reg_all_slopes_chosen_silva_tax %>%
  filter(slope_chosen_days != is.na(slope_chosen_days)) %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  group_by(asv_num) %>%
  mutate(counts = n()) %>%
  group_by(phylum_f, class_f, order_f, family_f, asv_num) %>%
              summarize(mean_growth_rate = mean(slope_chosen_days),
                        sd_growth_rate = sd(slope_chosen_days),
                        counts_growth_rate = counts) %>%
              distinct(asv_num, mean_growth_rate, sd_growth_rate, counts_growth_rate) %>%
              filter(counts_growth_rate > 12) %>%
              ungroup() %>%
              arrange(-mean_growth_rate)  %>%
              gt() %>%
              fmt_number(columns = contains('growth_rate'), decimals = 1) %>%
              cols_label(mean_growth_rate = 'Mean growth rate/day', sd_growth_rate = 'sd growth rate', phylum_f = 'Phylum', 
                         class_f = 'Class', order_f = 'Order', family_f = 'Family', asv_num = 'ASV number',
                         counts_growth_rate = 'Nº of significant growth rates calculated') %>%
              cols_align(align = "left", columns = contains('_f')) %>%
              cols_align(align = "center", columns = contains('growth_rate')) 
            # %>% si corro aquest últim punt no em funciona la creació de la taula.
            #   cols_width(contains('growth_rate') ~ px(80),
            #              everything() ~ px(200)) 
            
            # table_mean_gr_asv %>%
            #     glimpse()
            # setwd ("~/Documentos/Doctorat/REMEI/")
            
            #write.table(table_mean_gr_asv, 'results/tables/table_mean_gr_asv.txt', sep = '\t')
            
            table_mean_gr_asv %>%
              gt::gtsave(file = 'table_mean_gr_asv_filt.png', path = 'results/tables/') ##error attempt to apply non-function. 
            ##Error due to not specifing a mathematiccal operator or a function betwen the two values.
            
            library(kableExtra)
            
            ##mean per family
            table_mean_gr_family <-  
              reg_all_slopes_chosen_silva_tax %>%
              filter(slope_chosen_days != is.na(slope_chosen_days)) %>%
              filter(slope_chosen_days > 0,
                     pvalue_slope_chosen < 0.05) %>%
              group_by(family) %>%
              mutate(counts = n()) %>%
              group_by(phylum_f, class_f, order_f, family_f) %>%
              summarize(mean_growth_rate = mean(slope_chosen_days),
                        sd_growth_rate = sd(slope_chosen_days),
                        counts_growth_rate = counts) %>%
              distinct(family_f, mean_growth_rate, sd_growth_rate, counts_growth_rate) %>%
              filter(counts_growth_rate > 12) %>%
              ungroup() %>%
              arrange(-mean_growth_rate)  %>%
              gt() %>%
              fmt_number(columns = contains('growth_rate'), decimals = 1) %>%
              cols_label(mean_growth_rate = 'Mean growth rate/day', sd_growth_rate = 'sd growth rate', phylum_f = 'Phylum', 
                         class_f = md('Class'), order_f = 'Order', family_f = 'Family',
                         counts_growth_rate = 'Nº of significant growth rates calculated') %>%
              cols_align(align = "left", columns = contains('_f')) %>%
              cols_align(align = "center", columns = contains('growth_rate')) %>%
              cols_width(contains('growth_rate') ~ px(80),
                         everything() ~ px(180))
            
            # table_mean_gr_asv %>%
            #     glimpse()
            # setwd ("~/Documentos/Doctorat/REMEI/")
            
            # write.table(table_mean_gr_asv, 'results/tables/table_mean_gr_asv.txt', sep = '\t')
            
            table_mean_gr_family %>%
              gt::gtsave(file = 'table_mean_gr_family_filt.png', path = 'results/tables/') 
            
            
            ##mean per order
            table_mean_gr_order <-  
              reg_all_slopes_chosen_silva_tax %>%
              filter(slope_chosen_days != is.na(slope_chosen_days)) %>%
              filter(slope_chosen_days > 0,
                     pvalue_slope_chosen < 0.05) %>%
              group_by(order) %>%
              mutate(counts = n()) %>%
              group_by(phylum_f, class_f, order_f) %>%
              summarize(mean_growth_rate = mean(slope_chosen_days),
                        sd_growth_rate = sd(slope_chosen_days),
                        counts_growth_rate = counts) %>%
              distinct(order_f, mean_growth_rate, sd_growth_rate, counts_growth_rate) %>%
              filter(counts_growth_rate > 12) %>%
              ungroup() %>%
              arrange(-mean_growth_rate)  %>%
              gt() %>%
              fmt_number(columns = contains('growth_rate'), decimals = 1) %>%
              cols_label(mean_growth_rate = 'Mean growth rate/day', sd_growth_rate = 'sd growth rate', phylum_f = 'Phylum', 
                         class_f = md('Class'), order_f = 'Order',
                         counts_growth_rate = 'Nº of significant growth rates calculated') %>%
              cols_align(align = "left", columns = contains('_f')) %>%
              cols_align(align = "center", columns = contains('growth_rate')) %>%
              cols_width(contains('growth_rate') ~ px(90),
                         everything() ~ px(200))
            
            # table_mean_gr_asv %>%
            #     glimpse()
            # setwd ("~/Documentos/Doctorat/REMEI/")
            
            # write.table(table_mean_gr_asv, 'results/tables/table_mean_gr_asv.txt', sep = '\t')
            
            table_mean_gr_order %>%
              gt::gtsave(file = 'table_mean_gr_order_filt.png', path = 'results/tables/')
            

