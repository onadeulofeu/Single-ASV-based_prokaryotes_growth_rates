library(tidyverse)
library(magrittr)
library(phyloseq)
library(speedyseq)
library(forcats) #reorder categorical variables
library(gridExtra)
library(multipanelfigure) #unir plots amb diferents mides i distribucions
library(ggforce) #ampliació de ggplot2
library(ggplot2)
library(ggridges) #per fer els plots de les distribucions density més bonics
library(ggpubr) #as_gglot function
library('rstatix') #estadística
library(multcompView)
library(ggpmisc) ##stat_poli_line and equation
library(stringr)

#create table packages
library(gt)
library(gtExtras)
library(webshot2)
library(scales) #percent format in an axis
###LINIA 1565 COMENÇA SCRIPT BO AMB LES DADES DE CITOMETRIA PSEUDOABUND SIMPLIFICAT I NET 1754 BOXPLOTS-----

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
#palette_seasons_4 <- c("#3377a3", "#519741", "#ffbc42", "#cd5334") ##nous colors prova 
palette_seasons_4 <- c("#002562", "#519741", "#ffb900", "#96220a") #verd-blau no es confonen ara en principi
palette_phylums <- c("#fcca46", "#b0413e", "#669bbc", "#009e73", "#ffa737", "#cccccc",
                     "#69267e", "#fb9a99", "#1f78b4", "#8c789d", "#637b88", "#003029",
                     "#e3a6ce", "#002960", "#ba9864", "#005c69")

palette_phylums_complete <- c("#fcca46", "#669bbc", "#b0413e", "#009e73", "#ffa737", "#cccccc", "#69267e", "#1f78b4", 
                              "#8c789d", "#637b88", "#003029", "#e3a6ce", "#002960", "#ba9864", "#fb9a99", "#005c69", 
                              "#000000", "#c55e5c", "#00d198", "#5cb1d6", "#8a007a")

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

palette_class_assigned <- c('Gammaproteobacteria' = '#fcca46', 'Alphaproteobacteria' = '#d95726', 'Zetaproteobacteria' = '#EBCD92',
                            'Bacteroidia' = '#669BBC','Rhodothermia' =  '#00425E',
                            'Ignavibacteria' = '#CAFFFF', 'Chlorobia' = '#5BD1FF', 'Kryptonia' = '#0071A8',
                            'Nitrososphaeria' = '#B0413E',
                            'Cyanobacteriia' = '#009E73', 'Vampirivibrionia' = '#00733C',
                            'Acidimicrobiia' = '#FFA737','Actinobacteria' = '#DA6D00',
                            'Coriobacteriia' = '#C77800', 'Thermoleophilia' = '#A63B00',
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

reg_all_slopes_chosen_silva_tax %$%
  class_f %>%
  unique()

# [1] Cyanobacteriia      Alphaproteobacteria Bacteroidia         Gammaproteobacteria Planctomycetes      Vicinamibacteria    Bdellovibrionia    
# [8] Nitrososphaeria     Acidimicrobiia      Lentisphaeria       Verrucomicrobiae    Actinobacteria      Phycisphaerae       Rhodothermia       
# [15] Bacilli             Nitrospiria         Campylobacteria     Nitrospinia         Myxococcia          Kiritimatiellae     Clostridia         
# [22] Armatimonadia       Desulfuromonadia    Negativicutes       Thermoanaerobaculia Polyangia           Deinococci          Coriobacteriia     
# [29] Fusobacteriia       Desulfobacteria     Desulfobulbia       Blastocatellia      Chloroflexia        Desulfovibrionia    Ignavibacteria     
# [36] Gemmatimonadetes    Chlamydiae          Spirochaetia        Acidobacteriae      Holophagae 

#fixing colors by phylum----
reg_all_slopes_chosen_silva_tax_filt <- reg_all_slopes_chosen_silva_tax %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  group_by(phylum_f) %>%
  #summarize(n = n())
  filter(n() > 2)
# reg_all_slopes_chosen_silva_tax_filt%>%
#   dim()
names(palette_phylums) = unique(reg_all_slopes_chosen_silva_tax_filt$phylum_f)

#fixing  colors by all phylums-----
names(palette_phylums_complete) = unique(reg_all_slopes_chosen_silva_tax$phylum_f)

#this function allows us to create many colors from our palette
palf_large <- colorRampPalette(new_palette_large) 
palf_large_phylums <- colorRampPalette(palette_phylums)
palf.large.test <- colorRampPalette("#002960")
##to use it scale_color_manual(values=palf_large(x))+ x=number of colors you want

#Prepare data------ 
##Pseudoabundances from DAPI counts------
##We need pseudoabundances 
pseudoabund_df_wide <- pseudoabund_df %>%
  pivot_wider(names_from = c("asv_num"), values_from = "pseudoabundance", 
              id_cols = c("treatment", "season", "time", "replicate", "hours_dec"))

pseudoabund_df_wide_m <- pseudoabund_df_wide %>%
  merge(exp.data.dapis, by = c("treatment", "season", "time", "replicate", "hours_dec"))

##subset per season, treatment 
# winter_cd <- pseudoabund_df_wide_m %>% 
#   filter(c(treatment == "DL" & season == "Winter"))
# 
# winter_cd <- winter_cd %>% 
#   mutate(t_rep = paste0(time, sep = "_", replicate))
# 
# winter_cd %$% t_rep
# 
# winter_cd_ed <- winter_cd %>%
#   select(starts_with("asv"),
#          matches("t_rep"))
 # select(-selected_file_name_new)  #quedar-me només amb selected file name

##
asv <- winter_cd_ed ##las abundancias seran pseudoabundancias en log. asv cols time rows

# row.names(asv)  <- asv[,4502]
#  str(asv)
# asv <- asv %>%
#   select(!"t_rep")

env <- exp.data.dapis %>%
  subset(c(treatment == "DL" & season == "Winter")) %>%
  as_tibble()

env %>% 
  colnames()

env <- env %>% 
  mutate(t_rep = paste0(time, sep = "_", replicate))

env %$%
  t_rep

#asv_t<-t(asv)
dim(asv)
dim(env)

#taxonomy
tax_table <- rem_dapis_filt@tax_table %>%
  mutate_tax_table(tax_ed = (paste(asv_num, class, family, genus, species, sep = "_")))

##Necesitamos un merge para que se correspondan tiempos con las muestras-----

##0 values will be substituted by 1 because it is mandatory for transforming pseudoabundances to log
# asv_ed <- asv %>% 
#   mutate(across(everything(), ~replace(., . == 0, "1"))) %>%
#   mutate(across(!t_rep, as.numeric))
# 
# asv_m <- left_join(asv_ed, env, by = "t_rep") #to make sure that time is correctly matched with samples
# str(asv_m)
# 
# lm(log(asv_m[,1])~asv_m$hours_dec)
# 
# fun<-function(abund, env.var){
#   test<-lm(log(abund)~env.var)
#   res<-c(summary(test)$coefficients[2,1],
#          summary(test)$coefficients[2,2],
#          summary(test)$coefficients[2,2]*qt(.975, df = summary(test)$df[2]),
#          summary(test)$coefficients[2,4],
#          summary(test)$coefficients[1,1],
#          summary(test)$coefficients[1,2],
#          summary(test)$coefficients[1,2]*qt(.975, df = summary(test)$df[2]),
#          summary(test)$coefficients[1,4],
#          summary(test)$adj.r.squared)
#   names(res)<-c("slope","slope.se","slope.ci","slope.pval","intercept","intercept.se","intercept.ci","intercept.pval","r2.adj")
#   res
#   #return(res)
# }
# 
# # La podemos probar con la primera especie, como hicimos arriba manualmente:
# # asv_m %>%
# #   as_tibble()
# # 
# # fun(abund = asv_m[,1], env.var = asv_m$hours_dec)
# 
# # Una vez esto funciona, podemos usar "apply()" para ejecutar esta funci??n para todas las columnas, i.e. todas las especies.
# # Apply toma una matriz y ejecuta para cada fila o columna la funci??n que le indiquemos. X es la matriz, FUN indica la funci??n y
# # MARGIN indica si la función se ejecuta para cada fila (1) o columna (2). Cualquier otro argumento va luego.
# # Este mismo paso se podr??a hacer con un loop, pero as?? es mejor.
# 
# resultado <- apply(X = asv_m[,c(1:3736)], MARGIN = 2, FUN = fun, env.var= asv_m$hours_dec)
# resultado <- as.data.frame(t(resultado)) # le damos la vuelta y convertimos en data frame
# 
# resultado #hay NAs averiguar porqué
# resultado %T>% 
#   dim() %>%
#   head()
# 
# ##classify by significant p-value
# resultado_filt <- resultado %>% 
#   subset(slope.pval > 0.01) %>%
#   subset(intercept.pval > 0.01) %>%
#   subset(slope > 0.08)
# 
# significant_asv <- resultado_filt %>% #1428
#   row.names() %>%
#   as_tibble()
# 
# asv_ed_filt <- asv_ed %>%
#   subset( select = significant_asv$value) %>%
#   cbind(t_rep = asv_ed$t_rep )
# 
# #filter all asv that were 0 before transforming 0 to 1. this value will be the number of rows for each experimental
# #subset
# asv_m_filt <- asv_ed_filt %>%
#   select_if(~any (. > 11.0)) %>% #no es filtra cap
#   left_join(env, by = "t_rep")
# 
# asv_m_filt %>%
#   colnames()
# 
# ##plots de nuestras regresiones-----
# plot(log(asv_m[,1]) ~ asv_m$hours_dec)
# env %>% colnames()
# asv_m_long <- asv_m_filt %>% 
#   pivot_longer(cols = starts_with("asv")) #to plot all asv at the same time to get an overview
# 
# asv_m_long %>% 
#   colnames()
# 
# asv_m_filt %>%
#   colnames()
# tax_table %>%
#   colnames()
# 
# asv_m_long <- asv_m_long %>% 
#   left_join(tax_table, by = c("name"= "asv_num"), copy = TRUE) 
# 
# asv_m_long %>% ggplot(aes(hours_dec, value))+
#   geom_point()+
#   geom_smooth(method = "lm")+
#   facet_wrap(.~tax_ed)+
#   theme_bw()

##regression using 3 points


##comparison 3 points regression vs 4 points (sd increases? ) make a decision


##I edit the script so that it becames more authomatic (creating functions)----------
##FUNCTIONS USED (AIXÍ ÉS MÉS FÀCIL D'EDITAR)-------
filter.season.treatment <- function(data, treatment_col, treatment_keep, season_col, season_keep){
  my_var1 <- rlang::enquo(treatment_col)
  my_var2 <- rlang::enquo(season_col)
  data_filt <- data %>%
    filter(!!my_var1 == treatment_keep & !!my_var2 == season_keep)
  return(data_filt)
}

#function to create a new column combining time and replicate 
new.col <- function(data, old_col1, old_col2, new_col){
  old_col1 <- enquo(old_col1)
  old_col2 <- enquo(old_col2)
  new_col <- quo_name(new_col)
  data_ed <- dplyr::mutate(data, !!new_col := (paste0(!!old_col1, sep = "_", !!old_col2))) 
  return(data_ed)
}

fun<-function(abund, env.var){
  test<-lm(log(abund)~env.var)
  res<-c(summary(test)$coefficients[2,1],
         summary(test)$coefficients[2,2],
         summary(test)$coefficients[2,2]*qt(.975, df = summary(test)$df[2]),
         summary(test)$coefficients[2,4],
         summary(test)$coefficients[1,1],
         summary(test)$coefficients[1,2],
         summary(test)$coefficients[1,2]*qt(.975, df = summary(test)$df[2]),
         summary(test)$coefficients[1,4],
         summary(test)$adj.r.squared)
  names(res)<-c("slope","slope.se","slope.ci","slope.pval","intercept","intercept.se","intercept.ci","intercept.pval","r2.adj")
  res
  #return(res)
}

filter.significant.reg <- function(reg_result, pval, slope_filt, asv_df){
  reg_filt <- reg_result %>% 
    subset(slope.pval > pval) %>%
    subset(intercept.pval > pval) %>%
    subset(slope > slope_filt)
  
  significant_asv_reg <- reg_filt %>% #1428
    row.names() %>%
    as_tibble()
  
  asv_ed_filt <- asv_df %>%
    subset(select = significant_asv_reg$value) %>%
    cbind(t_rep = asv_ed$t_rep )
  print(significant_asv_reg)
  return(asv_ed_filt)
}


plot.reg <- function(data){
  data %>%  
    ggplot(aes(hours_dec, value, color = treatment, shape = season))+
    geom_point()+
    labs(x= "Time (h)", y = "log(pseudoabundance)")+
    geom_smooth(method = "lm")+
    scale_color_manual(values = palette_treatments_remei)+
    facet_wrap(.~tax_ed, scales = "free_y")+
    theme_bw()+
    theme(strip.text.x = element_text(size = 3))
  
}

##We need an asv_tab with pseudoabundances & metadata------
#asv_tab pseudoabundances (recalculate with pseudoabundances from flow citometry data)
pseudoabund_df_wide <- pseudoabund_df %>%
  pivot_wider(names_from = c("asv_num"), values_from = "pseudoabundance", 
              id_cols = c("treatment", "season", "time", "replicate", "hours_dec"))

#metadata 
exp.data.dapis <- rem_dapis_relabun@sam_data
 

#asv_tab & metadata
pseudoabund_df_wide_m <- pseudoabund_df_wide %>%
  merge(exp.data.dapis, by = c("treatment", "season", "time", "replicate", "hours_dec"))

#asv's taxonomy this new column is important for ploting asv with their taxonomic information
tax_table <- rem_dapis_filt@tax_table %>%
  mutate_tax_table(tax_ed = (paste(asv_num, class, family, genus, species, sep = "_")))
tax_table %>% 
  View()

tax_table_tab <- tax_table %>%
  as_tibble()
tax_table_tab %>%
  head()

##subset per season & treatment. We need to calculate regressions for each treatment and season------
asv_filt <- filter.season.treatment(data = pseudoabund_df_wide_m,  
                                treatment_col = treatment, treatment_keep = "VL", 
                                season_col = season, season_keep = "Winter")

env_filt <- exp.data.dapis %>%
  as_tibble() %>%
  filter.season.treatment(treatment_col = treatment, 
                                        treatment_keep = "VL", season_col = season, 
                                        season_keep = "Winter")

env_filt_ed <- env_filt %>%
  new.col(old_col1 = time, old_col2 = replicate, new_col = "t_rep")

asv_filt_ed <- asv_filt %>%
  new.col(old_col1 = time, old_col2 = replicate, new_col = "t_rep") %>%
  select(starts_with("asv"), #we clean metadata columns that we don't need now and keep only the column t_rep
         matches("t_rep"))

###
asv <- asv_filt_ed ##simplyfy name. Las abundancias seran pseudoabundancias en log. ASV cols time_rep rows
env <- env_filt_ed
asv %>% 
  dim() == env_filt_ed %>% 
  dim() #true for rows, false for columns should be

##0 values will be substituted by 1 because it is mandatory for transforming pseudoabundances to log
asv_ed <- asv %>% 
  mutate(across(everything(), ~replace(., . == 0, "1"))) %>%
  mutate(across(!t_rep, as.numeric))

asv_m <- left_join(asv_ed, env, by = "t_rep") #to make sure that time is correctly matched with samples
str(asv_m)

# lm(log(asv_m[,1])~asv_m$hours_dec)

# La podemos probar con la primera especie, como hicimos arriba manualmente:
# asv_m %>%
#   as_tibble()
# 
# fun(abund = asv_m[,1], env.var = asv_m$hours_dec)

# Una vez esto funciona, podemos usar "apply()" para ejecutar esta funci??n para todas las columnas, i.e. todas las especies.
# Apply toma una matriz y ejecuta para cada fila o columna la funci??n que le indiquemos. X es la matriz, FUN indica la funci??n y
# MARGIN indica si la función se ejecuta para cada fila (1) o columna (2). Cualquier otro argumento va luego.
# Este mismo paso se podr??a hacer con un loop, pero as?? es mejor.

dim(asv_m)
asv_m %>%
  colnames()
env %>%
  dim()
  
asv_m$hours_dec
env %>%
  colnames()
resultado <- apply(X = asv_m[,c(1:4502)], MARGIN = 2, FUN = fun, env.var= asv_m$hours_dec)
resultado <- as.data.frame(t(resultado)) # le damos la vuelta y convertimos en data frame
#resultado #hay NAs averiguar porqué
resultado %T>% 
  dim() %>%
  head()

asv_ed_filt <- filter.significant.reg(reg_result = resultado, pval = 0.01, slope_filt = 0.08, asv_df = asv_ed)

#filter all asv that were 0 before transforming 0 to 1. this value will be the number of rows for each experimental
#subset
asv_m_filt <- asv_ed_filt %>%
  select_if(~any (. > 11.0)) %>% #estem filtrant per les files que sumen més que si totes tinguessin 1 (que seria 0) comprovar que sigui correcte per totes les situuacions
  left_join(env, by = "t_rep") %>% 
  pivot_longer(cols = starts_with("asv")) %>% #to plot all asv at the same time to get an overview
  left_join(tax_table, by = c("name"= "asv_num"), copy = TRUE) 
# dim(asv_ed_filt)

asv_m_filt %>% 
  plot.reg()

##plots de nuestras regresiones ALL treatments-----
# plot(log(asv_m[,1]) ~ asv_m$hours_dec)
# env %>% colnames()
#PER EVITAR AIXÒ EL QUE HAURIA DE FER ÉS ANOMENAR EL DATAFRAME EN FUNCIÓ DE LA VARIABLE TREATMENT DINS D'UNA FUNCIÓ------
asv_m_long_cd <- asv_m_filt
asv_m_long_cl <- asv_m_filt 
asv_m_long_PL <- asv_m_filt
asv_m_long_pl <- asv_m_filt 
asv_m_long_dl <- asv_m_filt 
asv_m_long_vl <- asv_m_filt

# asv_m_long %>% 
#   colnames()
# asv_m_filt %>%
#   colnames()
# tax_table %>%
#   colnames()
asv_reg_all_treat_w <- bind_rows(asv_m_long_cd, 
                                 asv_m_long_cl, 
                                 asv_m_long_PL, 
                                 asv_m_long_pl, 
                                 asv_m_long_dl, 
                                 asv_m_long_vl)

asv_reg_all_treat_seas <- bind_rows(asv_reg_all_treat_w, 
                                    asv_reg_all_treat_sp,
                                    asv_reg_all_treat_su, 
                                    asv_reg_all_treat_f)

winter_reg <- asv_reg_all_treat_w %>% 
  plot.reg()
winter_reg

#all seasons all treatments
asv_reg_all_treat_seas %>%
  plot.reg()

asv_reg_all_treat_seas  %>% 
  ggplot(aes(hours_dec, value, color = treatment, shape = season))+
  geom_point()+
  labs(x= "Time (h", y = "log(pseudoabundance)")+
  geom_smooth(method = "lm")+
  scale_color_manual(values = palette_treatments_remei)+
  facet_wrap(.~season, scales = "free_y")+
  theme_bw()+
  theme(strip.text.x = element_text(size = 8))

asv_reg_all_treat_seas  %>% 
  ggplot(aes(hours_dec, value, color = treatment, shape = season))+
  geom_point()+
  labs(x= "Time (h", y = "log(pseudoabundance)")+
  geom_smooth(method = "lm")+
  scale_color_manual(values = palette_treatments_remei)+
  facet_wrap(.~genus, scales = "free_y")+
  theme_bw()+
  theme(strip.text.x = element_text(size = 8))


##subset by highest values because we don't see it properly.
# 
# function(reg_result, pval, slope_filt, asv_df){
#   value_filt <- asv_reg_all_treat_seas %>% 
#     subset(value > 1000)
#   
#   high_slope <- value_filt %$% #1428
#     name %>%
#     as_tibble()
#   
#   asv_ed_filt <- asv_reg_all_treat_seas %>%
#     filter(name %in% high_slope) %>%
#     cbind(t_rep = asv_ed$t_rep )
#   print(significant_asv_reg)
#   return(asv_ed_filt)
# }
# 
# asv_reg_all_treat_seas %>%
#   filter(eval(parse(text = paste(high_slope, collapse = "&"))))# %in% high_slope)
#  
# asv_reg_all_treat_seas %>% 
#   colnames()
# high_slope <- high_slope %>% 
#   rename("name"  = "value")
# 
# asv_reg_all_treat_seas %>% 
#   filter(!!! rlang::parse_exprs(high_slope))

value_filt <- asv_reg_all_treat_seas %>% 
  subset(value > 1000)
  high_slope <- value_filt %$% #1428
  name %>%
  as_tibble()

test <- asv_reg_all_treat_seas %>% 
  right_join(high_slope, copy = TRUE, by= "name") #10878 es multipliquen x10 les files

test %>% 
  colnames()

test %>% 
  ggplot(aes(hours_dec, value, color = treatment, shape = season))+
  geom_point()+
  labs(x= "Time (h)", y = "log(pseudoabundance)")+
  geom_smooth(method = "lm")+
  scale_color_manual(values = palette_treatments_remei)+
  facet_wrap(.~genus, scales = "free_y")+
  theme_bw()+
  theme(strip.text.x = element_text(size = 6))
  
# bind.df <- function(...){
#   global_df <- bind_rows(...,)
# return(global_df)
# }
# 
# bind.df(asv_reg_all_treat, asv_reg_all_treat_sp)

all_reg <- asv_reg_all_treat_seas  %>% 
  ggplot(aes(hours_dec, value, color = treatment, shape = season))+
  geom_point()+
  labs(x= "Time (h", y = "log(pseudoabundance)")+
  geom_smooth(method = "lm")+
  scale_color_manual(values = palette_treatments_remei)+
  facet_wrap(.~tax_ed, scales = "free_y")+
  theme_bw()+
  theme(strip.text.x = element_text(size = 3))
all_reg

summer_reg <- asv_reg_all_treat %>% 
  ggplot(aes(hours_dec, value, color = treatment, shape = season))+
  geom_point()+
  labs(x= "Time (h", y = "log(pseudoabundance)")+
  geom_smooth(method = "lm")+
  scale_color_manual(values = palette_treatments_remei)+
  facet_wrap(.~tax_ed, scales = "free_y")+
  theme_bw()+
  theme(strip.text.x = element_text(size = 3))
summer_reg

fall_reg <- asv_reg_all_treat %>% 
  ggplot(aes(hours_dec, value, color = treatment, shape = season))+
  geom_point()+
  labs(x= "Time (h", y = "log(pseudoabundance)")+
  geom_smooth(method = "lm")+
  scale_color_manual(values = palette_treatments_remei)+
  facet_wrap(.~tax_ed, scales = "free_y")+
  theme_bw()+
fall_reg

###Provo de fer group_by per no separar els datasets i plotejar-ho tot de cop (no em surt de moment)-------

#function to create a new column combining time and replicate---
new.col <- function(data, old_col1, old_col2, new_col){
  old_col1 <- enquo(old_col1)
  old_col2 <- enquo(old_col2)
  new_col <- quo_name(new_col)
  data_ed <- dplyr::mutate(data, !!new_col := (paste0(!!old_col1, sep = "_", !!old_col2))) 
  return(data_ed)
}

env_ed <- exp.data.dapis %>%
  as_tibble() %>%
  new.col(old_col1 = time, old_col2 = replicate, new_col = "t_rep")
env_ed %>%
  colnames()

asv_ed <- pseudoabund_df_wide_m %>%
  new.col(old_col1 = time, old_col2 = replicate, new_col = "t_rep") %>%
  select(starts_with("asv"), #we clean metadata columns that we don't need now and keep only the column t_rep
         matches("t_rep"))

asv_filt_ed %$% t_rep 

###
asv <- asv_ed ##simplify name. Las abundancias seran pseudoabundancias en log. ASV cols time_rep rows
env <- env_ed
asv %>% 
  dim() == env_filt_ed %>% 
  dim() #true for rows, false for columns should be

##0 values will be substituted by 1 because it is mandatory for transforming pseudoabundances to log
asv_ed <- asv %>% 
  mutate(across(everything(), ~replace(., . == 0, "1"))) %>%
  mutate(across(!t_rep, as.numeric))

asv_m <- left_join(asv_ed, env, by = "t_rep") #to make sure that time is correctly matched with samples
str(asv_m)

# lm(log(asv_m[,1])~asv_m$hours_dec)

fun<-function(abund, env.var){
  test<-lm(log(abund)~env.var)
  res<-c(summary(test)$coefficients[2,1],
         summary(test)$coefficients[2,2],
         summary(test)$coefficients[2,2]*qt(.975, df = summary(test)$df[2]),
         summary(test)$coefficients[2,4],
         summary(test)$coefficients[1,1],
         summary(test)$coefficients[1,2],
         summary(test)$coefficients[1,2]*qt(.975, df = summary(test)$df[2]),
         summary(test)$coefficients[1,4],
         summary(test)$adj.r.squared)
  names(res)<-c("slope","slope.se","slope.ci","slope.pval","intercept","intercept.se","intercept.ci","intercept.pval","r2.adj")
  res
  #return(res)
}

# La podemos probar con la primera especie, como hicimos arriba manualmente:
# asv_m %>%
#   as_tibble()
# 
# fun(abund = asv_m[,1], env.var = asv_m$hours_dec)

# Una vez esto funciona, podemos usar "apply()" para ejecutar esta funci??n para todas las columnas, i.e. todas las especies.
# Apply toma una matriz y ejecuta para cada fila o columna la funci??n que le indiquemos. X es la matriz, FUN indica la funci??n y
# MARGIN indica si la función se ejecuta para cada fila (1) o columna (2). Cualquier otro argumento va luego.
# Este mismo paso se podr??a hacer con un loop, pero as?? es mejor.

resultado <- apply(X = asv_m[,c(1:3736)], MARGIN = 2, FUN = fun, env.var= asv_m$hours_dec)
resultado <- as.data.frame(t(resultado)) # le damos la vuelta y convertimos en data frame

resultado #hay NAs averiguar porqué
resultado %T>% 
  dim() %>%
  head()

##classify by significant p-value, and highest slopes
resultado_filt <- resultado %>% 
  subset(slope.pval > 0.01) %>%
  subset(intercept.pval > 0.01) %>%
  subset(slope > 0.08)

significant_asv <- resultado_filt %>% #1428
  row.names() %>%
  as_tibble()

asv_ed_filt <- asv_ed %>%
  subset( select = significant_asv$value) %>%
  cbind(t_rep = asv_ed$t_rep )

#filter all asv that were 0 before transforming 0 to 1. this value will be the number of rows for each experimental
#subset
asv_m_filt <- asv_ed_filt %>%
  select_if(~any (. > 11.0)) %>% #no es filtra cap
  left_join(env, by = "t_rep") %>% 
  pivot_longer(cols = starts_with("asv")) %>% #to plot all asv at the same time to get an overview
  left_join(tax_table, by = c("name"= "asv_num"), copy = TRUE) 

asv_m_filt %>%
  ggplot(aes(hours_dec, value))+
  geom_point()+
  geom_smooth(method = "lm")+
  facet_wrap(.~tax_ed)+
  theme_bw()

##IMPORTING DATASET FOR CALCULATING THE REGRESSIONS
##PSEUDOABUNDANCES FROM FC DATA------------
#asv_tab pseudoabundances (recalculated with pseudoabundances from flow citometry data)
pseudoabund_df_wide_fc <- read.table("data/intermediate_files/regressions_test_datasets/pseudoabund_df_wide_fc.txt", 
                                     header = TRUE, sep = "\t")

#metadata 
exp.data.fc <- rem_fc_relabun@sam_data
exp.data.fc$treatment <- exp.data.fc$treatment %>% 
  factor(levels=(c("CL", "CD", "PL", "PD", "DL", "VL")))
exp.data.fc$season <- exp.data.fc$season %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))

#asv_tab & metadata
pseudoabund_df_wide_m <- pseudoabund_df_wide_fc %>%
  merge(exp.data.fc, by = c("treatment", "season", "time", "replicate", "hours.x"))

#asv's taxonomy this new column is important for ploting asv with their taxonomic information
tax_table <- rem_fc_filt@tax_table %>%
  mutate_tax_table(tax_ed = (paste(asv_num, class, family, genus, species, sep = "_")))

tax_table_tab <- tax_table %>%
  as_tibble()
tax_table_tab %>%
  head()

##subset per season & treatment. We need to calculate regressions for each treatment and season------
asv_filt <- filter.season.treatment(data = pseudoabund_df_wide_m,  
                                    treatment_col = treatment, treatment_keep = "CD", 
                                    season_col = season, season_keep = "Winter")

env_filt <- exp.data.fc %>%
  as_tibble() %>%
  filter.season.treatment(treatment_col = treatment, 
                          treatment_keep = "CD", season_col = season, 
                          season_keep = "Winter")

env_filt_ed <- env_filt %>%
  new.col(old_col1 = time, old_col2 = replicate, new_col = "t_rep")

asv_filt_ed <- asv_filt %>%
  new.col(old_col1 = time, old_col2 = replicate, new_col = "t_rep") %>%
  select(starts_with("asv"), #we clean metadata columns that we don't need now and keep only the column t_rep
         matches("t_rep"))

###
asv <- asv_filt_ed ##simplyfy name. Las abundancias seran pseudoabundancias en log. ASV cols time_rep rows
env <- env_filt_ed
asv %>% 
  dim() == env_filt_ed %>% 
  dim() #true for rows, false for columns should be

##0 values will be substituted by 1 because it is mandatory for transforming pseudoabundances to log
asv_ed <- asv %>% 
  mutate(across(everything(), ~replace(., . == 0, "1"))) %>%
  mutate(across(!t_rep, as.numeric))

asv_m <- left_join(asv_ed, env, by = "t_rep") #to make sure that time is correctly matched with samples
str(asv_m)
dim(asv_m)
dim(asv_ed)
##regressions
##MATRIX DIMENSIONS SHOULD BE CHECKED BEFORE APPLYING FUNCTION
resultado <- apply(X = asv_m[,c(1:4501)], MARGIN = 2, FUN = fun, env.var= asv_m$hours.x)
resultado <- as.data.frame(t(resultado)) # le damos la vuelta y convertimos en data frame
#resultado #hay NAs averiguar porqué
resultado %T>% 
  dim() %>%
  head()

asv_ed_filt <- filter.significant.reg(reg_result = resultado, pval = 0.01, slope_filt = 0.08, asv_df = asv_ed)

#filter all asv that were 0 before transforming 0 to 1. this value will be the number of rows for each experiment
#de moment he eliminat aquesta fila del script pensar bé com s'hauria de fer, potser es pot obviar i que no afecti.
#subset
asv_m_filt <- asv_ed_filt %>%
  #select_if(~any (. > 11.0)) %>% #estem filtrant per les files que sumen més que si totes tinguessin 1 (que seria 0) comprovar que sigui correcte per totes les situacions
  left_join(env, by = "t_rep") %>% 
  pivot_longer(cols = starts_with("asv")) %>% #to plot all asv at the same time to get an overview
  left_join(tax_table, by = c("name"= "asv_num"), copy = TRUE) 
# dim(asv_ed_filt)

asv_m_filt %>% 
  plot.reg.fc
# plot.reg.fc  <- function(data){
#     data %>%
#       ggplot(aes(hours.x, value, color = treatment, shape = season))+
#       geom_point()+
#       labs(x= "Time (h)", y = "log(pseudoabundance)")+
#       geom_smooth(method = "lm")+
#       scale_color_manual(values = palette_treatments_remei)+
#       facet_wrap(.~tax_ed, scales = "free_y")+
#       theme_bw()+
#       theme(strip.text.x = element_text(size = 3))
# 
#   }

##plots de nuestras regresiones ALL treatments-----
##tornar a la línia 609 i fer els càlculs per un altre treatment-season.
#PER EVITAR AIXÒ EL QUE HAURIA DE FER ÉS ANOMENAR EL DATAFRAME EN FUNCIÓ DE LA VARIABLE TREATMENT DINS D'UNA FUNCIÓ-----

asv_m_long_cd <- asv_m_filt
asv_m_long_cl <- asv_m_filt 
asv_m_long_PL <- asv_m_filt
asv_m_long_pl <- asv_m_filt 
asv_m_long_dl <- asv_m_filt 
asv_m_long_vl <- asv_m_filt

asv_reg_all_treat_f <- bind_rows(asv_m_long_cd, 
                                 asv_m_long_cl, 
                                 asv_m_long_PL, 
                                 asv_m_long_pl, 
                                 asv_m_long_dl, 
                                 asv_m_long_vl)

asv_reg_all_treat_seas <- bind_rows(asv_reg_all_treat_w, 
                                    asv_reg_all_treat_sp,
                                    asv_reg_all_treat_su, 
                                    asv_reg_all_treat_f)

spring_reg <- asv_reg_all_treat_sp %>% 
  plot.reg.fc()
winter_reg

write.table(asv_reg_all_treat_seas, "data/intermediate_files/regressions_test_datasets/asv_reg_all_treat_seas.txt", sep = "\t")
write.table(tax_table,"data/intermediate_files/regressions_test_datasets/tax_table_silva.txt", sep = "\t")

#all seasons all treatments
asv_reg_all_treat_seas %>%
  plot.reg.fc()

asv_reg_all_treat_seas  %>% 
  ggplot(aes(hours_dec, value, color = treatment, shape = season))+
  geom_point()+
  labs(x= "Time (h)", y = "log(pseudoabundance)")+
  geom_smooth(method = "lm")+
  scale_color_manual(values = palette_treatments_remei)+
  facet_wrap(.~season, scales = "free_y")+
  theme_bw()+
  theme(strip.text.x = element_text(size = 8))

asv_reg_all_treat_seas  %>% 
  ggplot(aes(hours.x, value, color = tax_ed, shape = season))+
  geom_point()+
  labs(x= "Time (h", y = "log(pseudoabundance)")+
  geom_smooth(method = "lm")+
  scale_color_manual(values = palf_large(405))+
  facet_wrap(.~treatment~season, scales = "free_y")+
  theme_bw()+
  theme(strip.text.x = element_text(size = 14), legend.position = "none")

##Study regressions by Genus
asv_reg_all_treat_seas %$% 
  genus %>% 
  unique()

asv_reg_all_treat_seas %>%
  subset(genus == "Fluviicola") %>%
  plot.reg.fc()+
  theme(strip.text.x = element_text(size = 8))

##interessants: Pseudoalteromonas

#Study regressions by Family
asv_reg_all_treat_seas %$% 
  family %>% 
  unique()

asv_reg_all_treat_seas %>%
  subset(family == "Candidatus Hepatincola") %>%
  ggplot(aes(hours, value, color = treatment, shape = season))+
  geom_point()+
  labs(x= "Time (h)", y = "log(pseudoabundance)")+
  geom_smooth(method = "lm")+
  scale_color_manual(values = palette_treatments_remei)+
  facet_wrap(.~genus, scales = "free_y")+
  theme_bw()+
  theme(strip.text.x = element_text(size = 14))

##interessants: Halieaceae, Rhodobacteraceae, Pseudoalteromonadaceae, Colwelliaceae
#Cryomorphaceae, Cyclobacteriaceae, Sphingomonadaceae, Marinomonadaceae, Hyphomonadaceae
#Parvularculaceae, Saprospiraceae, Balneolaceae, Rhodothermaceae

## [1] "Flavobacteriaceae", "Vibrionaceae", "Halieaceae", "Microtrichaceae"          
# [5] "NS9 marine group", "Porticoccaceae","Clade I", "Actinomycetaceae", 
# [9] "Crocinitomicaceae", "Spongiibacteraceae","MWH-UniP1 aquatic group""Rhodobacteraceae", 
# [13] "Comamonadaceae","Pirellulaceae","Arenicellaceae","Shewanellaceae", 
# [17] "Pseudoalteromonadaceae","Alcanivoracaceae1", "Colwelliaceae","Caulobacteraceae", 
# [21] "Cryomorphaceae","Alteromonadaceae", "Cyclobacteriaceae", "Litoricolaceae", 
# [25] "Cellvibrionaceae", "UBA7415", "Sphingomonadaceae", "Beijerinckiaceae", 
# [29] "Marinomonadaceae", "SAR116 clade", "Thioglobaceae","Microbacteriaceae",
# [33] "NS11-12 marine group", "Hyphomonadaceae", "Parvularculaceae", "Saprospiraceae", 
# [37] "Halomonadaceae","Moraxellaceae","Burkholderiaceae", "Balneolaceae",
# [41] "Xanthobacteraceae", "Phycisphaeraceae", "Devosiaceae", "Enterobacteriaceae",
# [45] "Saccharospirillaceae", "S25-593", "PS1 clade","AEGEAN-169 marine group", 
# [49] "Pseudomonadaceae", "Marinobacteraceae", "Methylophilaceae", "Rubinisphaeraceae",
# [53] "Bacteriovoracaceae","Cyanobiaceae", "Sporolactobacillaceae", "Puniceicoccaceae", 
# [57] "Bacillaceae", "Oxalobacteraceae", "Rhodothermaceae", "Azospirillaceae", 
# [61] "Neisseriaceae","Clade IV", "Arcobacteraceae", "Aeromonadaceae", 
# [65] "Desulfocapsaceae", "Kordiimonadaceae", "Stappiaceae", "Thalassospiraceae",
# [69] "Phormidesmiaceae", "Rhizobiaceae", "Chitinophagaceae", "Xenococcaceae",
# [73] "Granulosicoccaceae","Pseudohongiellaceae","Clade III","Ruminococcaceae", 
# [77] "NS7 marine group", "Rubritaleaceae","Salinisphaeraceae", "Veillonellaceae", 
# [81] "Erysipelatoclostridiaceae" "Lachnospiraceae", "Unknown Family_3", "Xanthomonadaceae", 
# [85] "Diplorickettsiaceae","Micrococcaceae","Microscillaceae", "Candidatus Hepatincola" 

##hi ha asv que no es dibuixen, busco perquè. asv978 el problema són les hours.x pq?
asv_reg_all_treat_seas %>%
  subset(name == "asv978") %>%
  ggplot(aes(hours, value, color = treatment, shape = season))+
  geom_point()+
  labs(x= "Time (h)", y = "log(pseudoabundance)")+
  geom_smooth(method = "lm")+
  scale_color_manual(values = palette_treatments_remei)+
  facet_wrap(.~name, scales = "free_y")+
  theme_bw()+
  theme(strip.text.x = element_text(size = 14), legend.position = "none")

##boxplots by genus-----
asv_reg_all_treat_seas %$% 
  genus %>% 
  unique()

asv_reg_all_treat_seas %>%
  subset(genus == "Pseudoalteromonas") %>%
  ggplot(aes(hours, value, group=time, color =  name))+
  geom_boxplot(aes(hours, value))+
  labs(x= "Time (h)", y = "log(pseudoabundance)")+
  geom_smooth(method = "lm")+
  scale_color_manual(values = palf_large(10))+
  facet_grid(.~treatment~season, scales = "free_y")+
  theme_bw()+
  theme(strip.text.x = element_text(size = 14))

asv_reg_all_treat_seas %>%
  subset(genus == "Marinomonas") %>%
  ggplot(aes(hours, value, group=time, color =  name))+
  geom_boxplot(aes(hours, value))+
  labs(x= "Time (h)", y = "log(pseudoabundance)")+
  geom_smooth(method = "lm")+
  scale_color_manual(values = palf_large(10))+
  facet_grid(.~treatment~season, scales = "free_y")+
  theme_bw()+
  theme(strip.text.x = element_text(size = 14))

asv_reg_all_treat_seas %>%
  subset(genus == "Alteromonas") %>%
  ggplot(aes(hours, value, group=time, color =  name))+
  geom_boxplot(aes(hours, value))+
  labs(x= "Time (h)", y = "log(pseudoabundance)")+
  geom_smooth(method = "lm")+
  scale_color_manual(values = palf_large(10))+
  facet_grid(.~treatment~season, scales = "free_y")+
  theme_bw()+
  theme(strip.text.x = element_text(size = 14))

asv_reg_all_treat_seas %>%
  subset(genus == "Thalassotalea") %>%
  ggplot(aes(hours, value, group=time, color =  name))+
  geom_boxplot(aes(hours, value))+
  labs(x= "Time (h)", y = "log(pseudoabundance)")+
  geom_smooth(method = "lm")+
  scale_color_manual(values = palf_large(10))+
  facet_grid(.~treatment~season, scales = "free_y")+
  theme_bw()+
  theme(strip.text.x = element_text(size = 14))

asv_reg_all_treat_seas %>%
  subset(genus == "Colwellia") %>%
  ggplot(aes(hours, value, group=time, color =  name))+
  geom_boxplot(aes(hours, value))+
  labs(x= "Time (h)", y = "log(pseudoabundance)")+
  geom_smooth(method = "lm")+
  scale_color_manual(values = palf_large(10))+
  facet_grid(.~treatment~season, scales = "free_y")+
  theme_bw()+
  theme(strip.text.x = element_text(size = 14))

asv_reg_all_treat_seas %>%
  subset(genus == "Vicingus") %>%
  ggplot(aes(hours, value, group=time, color =  name))+
  geom_boxplot(aes(hours, value))+
  labs(x= "Time (h)", y = "log(pseudoabundance)")+
  geom_smooth(method = "lm")+
  scale_color_manual(values = palf_large(10))+
  facet_grid(.~treatment~season, scales = "free_y")+
  theme_bw()+
  theme(strip.text.x = element_text(size = 14))

asv_reg_all_treat_seas %>%
  subset(genus == "Litoricola") %>%
  ggplot(aes(hours, value, group=time, color =  name))+
  geom_boxplot(aes(hours, value))+
  labs(x= "Time (h)", y = "log(pseudoabundance)")+
  geom_smooth(method = "lm")+
  scale_color_manual(values = palf_large(10))+
  facet_grid(.~treatment~season, scales = "free_y")+
  theme_bw()+
  theme(strip.text.x = element_text(size = 14))

asv_reg_all_treat_seas %>%
  subset(genus == "Methylobacterium-Methylorubrum") %>%
  ggplot(aes(hours, value, group=time, color =  name))+
  geom_boxplot(aes(hours, value))+
  labs(x= "Time (h)", y = "log(pseudoabundance)")+
  geom_smooth(method = "lm")+
  scale_color_manual(values = palf_large(10))+
  facet_grid(.~treatment~season, scales = "free_y")+
  theme_bw()+
  theme(strip.text.x = element_text(size = 14))

asv_reg_all_treat_seas %>%
  subset(genus == "Candidatus Puniceispirillum") %>%
  ggplot(aes(hours, value, group=time, color =  name))+
  geom_boxplot(aes(hours, value))+
  labs(x= "Time (h)", y = "log(pseudoabundance)")+
  geom_smooth(method = "lm")+
  scale_color_manual(values = palf_large(10))+
  facet_grid(.~treatment~season, scales = "free_y")+
  theme_bw()+
  theme(strip.text.x = element_text(size = 14))

asv_reg_all_treat_seas %>%
  subset(genus == "Parvularcula") %>%
  ggplot(aes(hours, value, group=time, color =  name))+
  geom_boxplot(aes(hours, value))+
  labs(x= "Time (h)", y = "log(pseudoabundance)")+
  geom_smooth(method = "lm")+
  scale_color_manual(values = palf_large(10))+
  facet_grid(.~treatment~season, scales = "free_y")+
  theme_bw()+
  theme(strip.text.x = element_text(size = 14))

asv_reg_all_treat_seas %>%
  subset(genus == "Balneola") %>%
  ggplot(aes(hours, value, group=time, color =  name))+
  geom_boxplot(aes(hours, value))+
  labs(x= "Time (h)", y = "log(pseudoabundance)")+
  geom_smooth(method = "lm")+
  scale_color_manual(values = palf_large(10))+
  facet_grid(.~treatment~season, scales = "free_y")+
  theme_bw()+
  theme(strip.text.x = element_text(size = 14))

asv_reg_all_treat_seas %>%
  subset(genus == "Rubrivirga") %>%
  ggplot(aes(hours, value, group=time, color =  name))+
  geom_boxplot(aes(hours, value))+
  labs(x= "Time (h)", y = "log(pseudoabundance)")+
  geom_smooth(method = "lm")+
  scale_color_manual(values = palf_large(10))+
  facet_grid(.~treatment~season, scales = "free_y")+
  theme_bw()+
  theme(strip.text.x = element_text(size = 14))



##NEW IDEA OF SCRIPT: CALCULATE REGRESSIONS (NO FILTERING), HUGE DATATABLE WITH ALL REGRESSIONS-----
##load functions
source("src/filter.season.treatment.time.R") ###filtrar les pseudoabund per treatment season i time per calcular regressions
source("src/multiple.linear.regressions.R")

##subset per season & treatment & time. We need to calculate regressions for each treatment and season------
filter.season.treatment.time <- function(data, treatment_col, treatment_keep, season_col, season_keep, time_col, time_keep){
  my_var1 <- rlang::enquo(treatment_col)
  my_var2 <- rlang::enquo(season_col)
  my_var3 <- rlang::enquo(time_col)
  data_filt <- data %>%
    filter(!!my_var1 == treatment_keep & !!my_var2 == season_keep & !!my_var3 %in% time_keep)
  return(data_filt)
}

#asv_tab pseudoabundances (recalculated with pseudoabundances from flow citometry data)----
pseudoabund_df_wide_fc <- read.table("data/intermediate_files/regressions_test_datasets/pseudoabund_df_wide_fc.txt", 
                                     header = TRUE, sep = "\t")

#metadata 
exp.data.fc <- rem_fc_relabun@sam_data
exp.data.fc$treatment <- exp.data.fc$treatment %>% 
  factor(levels=(c("CL", "CD", "PL", "PL", "DL", "VL")))
exp.data.fc$season <- exp.data.fc$season %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))

#asv_tab & metadata
pseudoabund_df_wide_m <- pseudoabund_df_wide_fc %>%
  merge(exp.data.fc, by = c("treatment", "season", "time", "replicate", "hours.x"))

#asv's taxonomy this new column is important for ploting asv with their taxonomic information
tax_table <- rem_fc_filt@tax_table %>%
  mutate_tax_table(tax_ed = (paste(asv_num, class, family, genus, species, sep = "_")))

tax_table_tab <- tax_table %>%
  as_tibble()
tax_table_tab %>%
  head()

##subset per season & treatment. We need to calculate regressions for each treatment and season------
##com que són les regressions de 4 temps agafo els 4 a time_keep
asv_filt <- filter.season.treatment.time(data = pseudoabund_df_wide_m,  
                                              treatment_col = treatment, treatment_keep = "CD", 
                                              season_col = season, season_keep = "Winter",
                                              time_col = time, time_keep = c("t0", "t2", "t3", "t4"))

env_filt <- exp.data.fc %>%
  as_tibble() %>%
  filter.season.treatment(treatment_col = treatment, 
                          treatment_keep = "PL", season_col = season, 
                          season_keep = "Fall")

env_filt_ed <- env_filt %>%
  new.col(old_col1 = time, old_col2 = replicate, new_col = "t_rep")

asv_filt_ed <- asv_filt %>%
  new.col(old_col1 = time, old_col2 = replicate, new_col = "t_rep") %>%
  select(starts_with("asv"), #we clean metadata columns that we don't need now and keep only the column t_rep
         matches("t_rep"))

###
asv <- asv_filt_ed ##simplyfy name. Las abundancias seran pseudoabundancias en log. ASV cols time_rep rows
env <- env_filt_ed
asv %>% 
  dim() == env_filt_ed %>% 
  dim() #It should be true for rows, false for columns 

##0 values will be substituted by 1 because it is mandatory for transforming pseudoabundances to log
asv_ed <- asv %>% 
  mutate(across(everything(), ~replace(., . == 0, "1"))) %>%
  mutate(across(!t_rep, as.numeric))

asv_m <- left_join(asv_ed, env, by = "t_rep") #to make sure that time is correctly matched with samples

asv_ed %>%
  dim()
##regressions
##DIMENSIONS OF THE MATRIX SHOULD BE ADAPTED TO DATA!!!!!
resultado <- apply(X = asv_m[,c(1:4501)], MARGIN = 2, FUN = fun, env.var= asv_m$hours.x)
resultado <- as.data.frame(t(resultado)) # le damos la vuelta y convertimos en data frame
resultado  <- resultado %>%
  rownames_to_column(var = "asv_num")
#resultado #hay NAs averiguar porqué crec que és per les dimensions de la taula estaven malament.
resultado %T>% 
  dim() %>%
  head()

##merge resultado with asv_m 
#cambiar nombre  del objeto en función del tratamiento y la estación
asv_m_reg_fall_pl <- asv_m %>% 
  select(starts_with("asv"), ".sample", "treatment", "replicate",  "time", "season", "sample_name", "sample_code", "selected_file_name",
         "reads", "light_regime", "hours", "LNA.y", "HNA.y", "fc", "BM.y", "Leu.PB.y", "SGR.y", "TD.y", "t_rep")  %>%
  pivot_longer(cols = starts_with("asv"),names_to = "asv_num", values_to = "pseudoabundance") %>%
  left_join(resultado, by = "asv_num")

asv_reg_all_treat_f <- bind_rows(asv_m_reg_fall_cd, 
                                 asv_m_reg_fall_cl, 
                                 asv_m_reg_fall_PL, 
                                 asv_m_reg_fall_pl, 
                                 asv_m_reg_fall_dl, 
                                 asv_m_reg_fall_vl)

asv_reg_all_treat_seas_no_filter <- bind_rows(asv_reg_all_treat_w, 
                                    asv_reg_all_treat_sp,
                                    asv_reg_all_treat_su, 
                                    asv_reg_all_treat_f)

write.table(asv_reg_all_treat_seas_no_filter, "asv_reg_all_treat_seas_no_filter.txt", sep ="\t")

##UPLOAD DATA totes les regressions de tots els treatments amb 4 temps
asv_reg_all_treat_seas_no_filter <- read.table("data/intermediate_files/regressions_test_datasets/asv_reg_all_treat_seas_no_filter.txt", 
                                               header = TRUE, sep = "\t") %>% 
  as_tibble()

asv_reg_all_treat_seas_no_filter %>%
  dim()
##plot by families we add taxonomy first
asv_reg_all_treat_seas_no_filter <- asv_reg_all_treat_seas_no_filter %>%
  left_join(tax_table, by = "asv_num", copy = TRUE) 

asv_reg_all_treat_seas_no_filter %>%
  colnames()
asv_reg_all_treat_seas_no_filter %$% 
  family %>%
  unique()

##busco les families amb major taxa de creixement
asv_reg_all_treat_seas_no_filter %>%
  filter(slope > 0.1) %>%
  filter(slope.pval < 0.05 & intercept.pval < 0.05) %$%
  family %>%
  unique()

# [1] "Vibrionaceae"  "Alteromonadaceae"   "Flavobacteriaceae"  "Marinomonadaceae"      
# [5] "Rhodobacteraceae" "Pseudoalteromonadaceae" "Beijerinckiaceae"  "Colwelliaceae"         
# [9] "Phycisphaeraceae"    "Cryomorphaceae"   "Saprospiraceae" "NS9 marine group"      
# [13] "Rhizobiaceae"  "Sphingomonadaceae" 

asv_reg_all_treat_seas_no_filter %$% 
  genus %>%
  unique()

##necessita molta memòria per això necessitem subset
##families
plot.tax.filtered(data = asv_reg_all_treat_seas_no_filter,
                  value_slope_pval = 0.05, 
                  value_intercept_pval = 0.05, 
                  tax_level = family, 
                  select_tax_level  = "Rhodobacteraceae")

# plot.tax.filtered(data = asv_reg_all_treat_seas_no_filter,
#                   value_slope_pval = 0.05, 
#                   value_intercept_pval = 0.05, 
#                   tax_level = family, 
#                   select_tax_level  = "Thalassobaculaceae")

plot.tax.filtered(data = asv_reg_all_treat_seas_no_filter,
                  value_slope_pval = 0.05, 
                  value_intercept_pval = 0.05, 
                  tax_level = family, 
                  select_tax_level  = "Halieaceae")

plot.tax.filtered(data = asv_reg_all_treat_seas_no_filter,
                  value_slope_pval = 0.05, 
                  value_intercept_pval = 0.05, 
                  tax_level = family, 
                  select_tax_level  = "Pseudoalteromonadaceae")

plot.tax.filtered(data = asv_reg_all_treat_seas_no_filter,
                  value_slope_pval = 0.05, 
                  value_intercept_pval = 0.05, 
                  tax_level = family, 
                  select_tax_level  = "Colwelliaceae")

plot.tax.filtered(data = asv_reg_all_treat_seas_no_filter,
                  value_slope_pval = 0.05, 
                  value_intercept_pval = 0.05, 
                  tax_level = family, 
                  select_tax_level  = "Cryomorphaceae")

plot.tax.filtered(data = asv_reg_all_treat_seas_no_filter,
                  value_slope_pval = 0.05, 
                  value_intercept_pval = 0.05, 
                  tax_level = family, 
                  select_tax_level  = "Cyclobacteriaceae")

plot.tax.filtered(data = asv_reg_all_treat_seas_no_filter,
                  value_slope_pval = 0.05, 
                  value_intercept_pval = 0.05, 
                  tax_level = family, 
                  select_tax_level  = "Sphingomonadaceae")

plot.tax.filtered(data = asv_reg_all_treat_seas_no_filter,
                  value_slope_pval = 0.05, 
                  value_intercept_pval = 0.05, 
                  tax_level = family, 
                  select_tax_level  = "Marinomonadaceae")

plot.tax.filtered(data = asv_reg_all_treat_seas_no_filter,
                  value_slope_pval = 0.05, 
                  value_intercept_pval = 0.05, 
                  tax_level = family, 
                  select_tax_level  = "Parvularculaceae")

plot.tax.filtered(data = asv_reg_all_treat_seas_no_filter,
                  value_slope_pval = 0.05, 
                  value_intercept_pval = 0.05, 
                  tax_level = family, 
                  select_tax_level  = "Saprospiraceae")

plot.tax.filtered(data = asv_reg_all_treat_seas_no_filter,
                  value_slope_pval = 0.05, 
                  value_intercept_pval = 0.05, 
                  tax_level = family, 
                  select_tax_level  = "Balneolaceae")

plot.tax.filtered(data = asv_reg_all_treat_seas_no_filter,
                  value_slope_pval = 0.05, 
                  value_intercept_pval = 0.05, 
                  tax_level = family, 
                  select_tax_level  = "Rhodothermaceae")

plot.tax.filtered(data = asv_reg_all_treat_seas_no_filter,
                  value_slope_pval = 0.05, 
                  value_intercept_pval = 0.05, 
                  tax_level = family, 
                  select_tax_level  = "Vibrionaceae")

plot.tax.filtered(data = asv_reg_all_treat_seas_no_filter,
                  value_slope_pval = 0.05, 
                  value_intercept_pval = 0.05, 
                  tax_level = family, 
                  select_tax_level  = "Flavobacteriaceae")

plot.tax.filtered(data = asv_reg_all_treat_seas_no_filter,
                  value_slope_pval = 0.05, 
                  value_intercept_pval = 0.05, 
                  tax_level = family, 
                  select_tax_level  = "Rhodobacteraceae")

plot.tax.filtered(data = asv_reg_all_treat_seas_no_filter,
                  value_slope_pval = 0.05, 
                  value_intercept_pval = 0.05, 
                  tax_level = family, 
                  select_tax_level  = "Alteromonadaceae")

##genus
#busco gèneres amb major taxa de creixement
asv_reg_all_treat_seas_no_filter %>%
  filter(slope > 0.1) %>%
  filter(slope.pval < 0.05 & intercept.pval < 0.05) %$%
  genus %>%
  unique

# [1] NA                               "Pseudofulvibacter"              "Polaribacter"                  
# [4] "Glaciecola"                     "Dokdonia"                       "Marinomonas"                   
# [7] "Nonlabens"                      "Tenacibaculum"                  "Pseudoalteromonas"             
# [10] "Vibrio"                         "Methylobacterium-Methylorubrum" "NS3a marine group"             
# [13] "Octadecabacter"                 "Colwellia"                      "Palleronia-Pseudomaribius"     
# [16] "Formosa"                        "NS5 marine group"               "CL500-3"                       
# [19] "Owenweeksia"                    "NS2b marine group"              "Fulvimarina"                   
# [22] "Alteromonas"                    "Thalassotalea"   

plot.tax.filtered <- function(data, value_slope_pval, value_intercept_pval, tax_level, select_tax_level){
  my_var1 <- rlang::enquo(tax_level)
  data %>%
    filter(slope.pval < value_slope_pval & intercept.pval < value_intercept_pval) %>%
    filter(!!my_var1 == select_tax_level) %>%
    ggplot(aes(hours, pseudoabundance, color = tax_ed, shape = season))+
    geom_point()+
    labs(x= "Time (h)", y = "log(pseudoabundance)")+
    geom_smooth(method = "lm")+
    scale_color_manual(values = palf_large(70))+
    facet_grid(.~treatment~season, scales = "fixed")+
    guides(color=guide_legend(ncol = 1))+
    theme_bw()+
    theme(strip.text.x = element_text(size = 14))
}

library(rlang)
# plot.tax.filtered <- function(data, value_slope_pval, value_intercept_pval, tax_level, select_tax_level, x, y){
#   my_var1 <- rlang::enquo(tax_level)
#   data %>%
#     filter(slope.pval < value_slope_pval & intercept.pval < value_intercept_pval) %>%
#     filter(!!my_var1 == select_tax_level) %>%
#     ggplot(aes(hours, pseudoabundance, color = tax_ed, shape = season))+
#     geom_point()+
#     labs(x= "Time (h)", y = "log(pseudoabundance)")+
#     geom_smooth(method = "lm")+
#     scale_color_manual(values = palf_large(35))+
#     #facet_grid(vars, scales = "fixed")+
#     facet_grid(.~{{x}}~{{y}}, scales = "fixed")+
#     guides(color=guide_legend(ncol = 1))+
#     theme_bw()+
#     theme(strip.text.x = element_text(size = 14))
# }

# wrap_by <- function(...){
#   x <- enquo(facet1)
#   y <- enquo(facet2)
#   facet_wrap(.~x~y, scales = "fixed")
# }
#NO FUNCIONA LA PART DE FACET_GRID TROBAR LA MANERA 

plot.tax.filtered(data = asv_reg_all_treat_seas_no_filter,
                    value_slope_pval = 0.05, 
                    value_intercept_pval = 0.05, 
                    tax_level = genus, 
                    select_tax_level  = "Pseudofulvibacter")

plot.tax.filtered(data = asv_reg_all_treat_seas_no_filter,
                    value_slope_pval = 0.05, 
                    value_intercept_pval = 0.05, 
                    tax_level = genus, 
                    select_tax_level  = "Polaribacter")

plot.tax.filtered(data = asv_reg_all_treat_seas_no_filter,
                    value_slope_pval = 0.05, 
                    value_intercept_pval = 0.05, 
                    tax_level = genus, 
                    select_tax_level  = "Glaciecola")

plot.tax.filtered(data = asv_reg_all_treat_seas_no_filter,
                    value_slope_pval = 0.05, 
                    value_intercept_pval = 0.05, 
                    tax_level = genus, 
                    select_tax_level  = "Dokdonia")

plot.tax.filtered(data = asv_reg_all_treat_seas_no_filter,
                    value_slope_pval = 0.05, 
                    value_intercept_pval = 0.05, 
                    tax_level = genus, 
                    select_tax_level  = "Marinomonas")

plot.tax.filtered(data = asv_reg_all_treat_seas_no_filter,
                    value_slope_pval = 0.05, 
                    value_intercept_pval = 0.05, 
                    tax_level = genus, 
                    select_tax_level  = "Nonlabens")

plot.tax.filtered(data = asv_reg_all_treat_seas_no_filter,
                    value_slope_pval = 0.05, 
                    value_intercept_pval = 0.05, 
                    tax_level = genus, 
                    select_tax_level  = "Tenacibaculum")

plot.tax.filtered(data = asv_reg_all_treat_seas_no_filter,
                    value_slope_pval = 0.05, 
                    value_intercept_pval = 0.05, 
                    tax_level = genus, 
                    select_tax_level  = "Pseudoalteromonas")

plot.tax.filtered(data = asv_reg_all_treat_seas_no_filter,
                    value_slope_pval = 0.05, 
                    value_intercept_pval = 0.05, 
                    tax_level = genus, 
                    select_tax_level  = "Vibrio")

plot.tax.filtered(data = asv_reg_all_treat_seas_no_filter,
                    value_slope_pval = 0.05, 
                    value_intercept_pval = 0.05, 
                    tax_level = genus, 
                    select_tax_level  = "Methylobacterium-Methylorubrum")

plot.tax.filtered(data = asv_reg_all_treat_seas_no_filter,
                    value_slope_pval = 0.05, 
                    value_intercept_pval = 0.05, 
                    tax_level = genus, 
                    select_tax_level  = "NS3a marine group")

plot.tax.filtered(data = asv_reg_all_treat_seas_no_filter,
                    value_slope_pval = 0.05, 
                    value_intercept_pval = 0.05, 
                    tax_level = genus, 
                    select_tax_level  = "Octadecabacter")

plot.tax.filtered(data = asv_reg_all_treat_seas_no_filter,
                    value_slope_pval = 0.05, 
                    value_intercept_pval = 0.05, 
                    tax_level = genus, 
                    select_tax_level  = "Colwellia")

plot.tax.filtered(data = asv_reg_all_treat_seas_no_filter,
                    value_slope_pval = 0.05, 
                    value_intercept_pval = 0.05, 
                    tax_level = genus, 
                    select_tax_level  = "Palleronia-Pseudomaribius")

plot.tax.filtered(data = asv_reg_all_treat_seas_no_filter,
                    value_slope_pval = 0.05, 
                    value_intercept_pval = 0.05, 
                    tax_level = genus, 
                    select_tax_level  = "Formosa")

plot.tax.filtered(data = asv_reg_all_treat_seas_no_filter,
                    value_slope_pval = 0.05, 
                    value_intercept_pval = 0.05, 
                    tax_level = tax, 
                    select_tax_level  = "NS5 marine group")

plot.tax.filtered(data = asv_reg_all_treat_seas_no_filter,
                    value_slope_pval = 0.05, 
                    value_intercept_pval = 0.05, 
                    tax_level = tax, 
                    select_tax_level  = "CL500-3")

plot.tax.filtered(data = asv_reg_all_treat_seas_no_filter,
                    value_slope_pval = 0.05, 
                    value_intercept_pval = 0.05, 
                    tax_level = tax, 
                    select_tax_level  = "Owenweeksia")

plot.tax.filtered(data = asv_reg_all_treat_seas_no_filter,
                    value_slope_pval = 0.05, 
                    value_intercept_pval = 0.05, 
                    tax_level = tax, 
                    select_tax_level  = "NS2b marine group")

plot.tax.filtered(data = asv_reg_all_treat_seas_no_filter,
                    value_slope_pval = 0.05, 
                    value_intercept_pval = 0.05, 
                    tax_level = tax, 
                    select_tax_level  = "NS2b marine group")

plot.tax.filtered(data = asv_reg_all_treat_seas_no_filter,
                    value_slope_pval = 0.05, 
                    value_intercept_pval = 0.05, 
                    tax_level = tax, 
                    select_tax_level  = "Fulvimarina")

plot.tax.filtered(data = asv_reg_all_treat_seas_no_filter,
                    value_slope_pval = 0.05, 
                    value_intercept_pval = 0.05, 
                    tax_level = tax, 
                    select_tax_level  = "Alteromonas")

plot.tax.filtered(data = asv_reg_all_treat_seas_no_filter,
                    value_slope_pval = 0.05, 
                    value_intercept_pval = 0.05, 
                    tax_level = tax, 
                    select_tax_level  = "Thalassotalea")


###REGRESSIONS WITH 3 TIMES INSTEAD OF 4---------
##load functions----
source("src/filter.season.treatment.time.R") ###filtrar les pseudoabund per treatment season i time per calcular regressions

##Possibles combinacions: 
  ##t0-t2-t3
  ##t0-t3-t4
  ##t2-t3-t4
  ##t0-t2-t4  

filter.season.treatment.time <- function(data, treatment_col, treatment_keep, season_col, season_keep, time_col, time_keep){
  my_var1 <- rlang::enquo(treatment_col)
  my_var2 <- rlang::enquo(season_col)
  my_var3 <- rlang::enquo(time_col)
  data_filt <- data %>%
    filter(!!my_var1 == treatment_keep & !!my_var2 == season_keep & !!my_var3 %in% time_keep)
  return(data_filt)
}

##subset per season & treatment & time. We need to calculate regressions for each treatment and season and time------
##PENSAR SI HI HA UNA MANERA DE FER EL SCRIPT MÉS SENZILL O AUTOMÀTIC
##Prepare  datasets for calculating regressions------
asv_filt_t024 <- filter.season.treatment.time(data = pseudoabund_df_wide_m,  
                                    treatment_col = treatment, treatment_keep = "CD", 
                                    season_col = season, season_keep = "Winter",
                                    time_col = time, time_keep = c("t0", "t2", "t4")) %>%
  new.col(old_col1 = time, old_col2 = replicate, new_col = "t_rep") %>%
  select(starts_with("asv"), #we clean metadata columns that we don't need now and keep only the column t_rep
         matches("t_rep"))

env_filt_t024 <- exp.data.fc %>%
  as_tibble() %>%
  filter.season.treatment.time(treatment_col = treatment, 
                          treatment_keep = "CD", season_col = season, 
                          season_keep = "Winter",
                          time_col = time, time_keep = c("t0", "t2", "t4")) %>%
  new.col(old_col1 = time, old_col2 = replicate, new_col = "t_rep")

##simplyfy names
asv <- asv_filt_t024 ## Las abundancias seran pseudoabundancias en log. ASV cols time_rep rows
env <- env_filt_t024
asv %>% 
  dim() == env_filt_t024 %>% 
  dim() #It should be true for rows, false for columns 

##0 values will be substituted by 1 because it is mandatory for transforming pseudoabundances to log
asv_ed <- asv %>% 
  mutate(across(everything(), ~replace(., . == 0, "1"))) %>%
  mutate(across(!t_rep, as.numeric))

asv_m <- left_join(asv_ed, env, by = "t_rep") #to make sure that time is correctly matched with samples

asv_ed %>%
  dim()
##Calculating regressions
##DIMENSIONS OF THE MATRIX SHOULD BE ADAPTED TO DATA!!!!!
##PROVO de calcular les regressions utilitzant group_by per poder-ho fer tot més automàtic.
#script anterior
resultado <- apply(X = asv_m[,c(1:4501)], MARGIN = 2, FUN = fun, env.var= asv_m$hours.x)
resultado <- as.data.frame(t(resultado)) # le damos la vuelta y convertimos en data frame
resultado  <- resultado %>%
  rownames_to_column(var = "asv_num")  

pseudoabund_df_wide_m %>%
  dim()
pseudoabund_df_wide_m %>%
  colnames()
test <- pseudoabund_df_wide_fc %>% 
  pivot_longer(cols = starts_with("asv"), names_to = "asv", values_to = "pseudoabund") %>%
  as_tibble() %>%
  group_by(treatment, season, asv)
test %>%
  group_by(treatment, season, asv) %>%
  fun(abund = test$pseudoabund, env.var = test$hours.x)

  apply(X = pseudoabund_df_wide_m[,7], MARGIN = 2, FUN = fun, env.var = pseudoabund_df_wide_m$hours.x,)
test %>%
  str()
test %>%
  colnames()
test %>%
  glimpse()
test %>%
  as_tibble() %>%
  lm(log(test$pseudoabund)~test$hours.x, data = subset(test, treatment == "CD" & season == "Winter"))

test %$%
  asv %>%
  unique()
test %>%
  ggplot(aes(x = hours.x, y = pseudoabund, color = asv))+
  geom_boxplot()+
  scale_color_manual(values= palf_large(4501))

  res<-c(summary(test)$coefficients[2,1],
           summary(test)$coefficients[2,2],
           summary(test)$coefficients[2,2]*qt(.975, df = summary(test)$df[2]),
           summary(test)$coefficients[2,4],
           summary(test)$coefficients[1,1],
           summary(test)$coefficients[1,2],
           summary(test)$coefficients[1,2]*qt(.975, df = summary(test)$df[2]),
           summary(test)$coefficients[1,4],
           summary(test)$adj.r.squared)
    names(res)<-c("slope","slope.se","slope.ci","slope.pval","intercept","intercept.se","intercept.ci","intercept.pval","r2.adj")
    res
    #return(res)
 
pseudoabund_df_wide_m %>%
  str()
asv_m %>%
resultado <- pseudoabund_df_wide_m %>% 
  group_by(treatment, season) %>%
  group_modify(.f = fun, pseudoabund_df_wide_m[,c(6:4506)], pseudoabund_df_wide_m$hours.x)

test <- pseudoabund_df_wide_m %>% 
  group_by(treatment, season) %>%
  summarise( mod = list(lm(log(pseudoabund_df_wide_m[,c(6:4506)])) ~ pseudoabund_df_wide_m$hours.x)) 
test$mod[[1]]$coefficients
  
test$mod
apply(abund = pseudoabund_df_wide_m[,c(6:4506)], MARGIN = 2, FUN = fun, env.var= pseudoabund_df_wide_m$hours.x)

library(broom)
resultado <- pseudoabund_df_wide_m %>% 
  group_by(treatment, season) %>%
  do()
  group_modify(.f = fun, pseudoabund_df_wide_m[,c(6:4506)], pseudoabund_df_wide_m$hours.x)

   ( X = pseudoabund_df_wide_m[,c(6:4506)], MARGIN = 2, FUN = fun, env.var= pseudoabund_df_wide_m$hours.x)

resultado <- as.data.frame(t(resultado)) # le damos la vuelta y convertimos en data frame
resultado  <- resultado %>%
  rownames_to_column(var = "asv_num")
#resultado #hay NAs averiguar porqué crec que és per les dimensions de la taula estaven malament.
resultado %T>% 
  dim() %>%
  head()

##merge resultado with asv_m 
#cambiar nombre  del objeto en función del tratamiento y la estación
asv_m_reg_fall_pl <- asv_m %>% 
  select(starts_with("asv"), ".sample", "treatment", "replicate",  "time", "season", "sample_name", "sample_code", "selected_file_name",
         "reads", "light_regime", "hours", "LNA.y", "HNA.y", "fc", "BM.y", "Leu.PB.y", "SGR.y", "TD.y", "t_rep")  %>%
  pivot_longer(cols = starts_with("asv"),names_to = "asv_num", values_to = "pseudoabundance") %>%
  left_join(resultado, by = "asv_num")

asv_reg_all_treat_f <- bind_rows(asv_m_reg_fall_cd, 
                                 asv_m_reg_fall_cl, 
                                 asv_m_reg_fall_PL, 
                                 asv_m_reg_fall_pl, 
                                 asv_m_reg_fall_dl, 
                                 asv_m_reg_fall_vl)

asv_reg_all_treat_seas_no_filter <- bind_rows(asv_reg_all_treat_w, 
                                              asv_reg_all_treat_sp,
                                              asv_reg_all_treat_su, 
                                              asv_reg_all_treat_f)

write.table(asv_reg_all_treat_seas_no_filter, "asv_reg_all_treat_seas_no_filter.txt", sep ="\t")



##SIMPLYFYING THE SCRIPT AND OPTIMIZING FC DATA INCLUDING GdC EXPERIMENT---- 
#functions needed
source("src/filter.season.treatment.time.R") ###filtrar les pseudoabund per treatment season i time per calcular regressions
source("src/multiple.linear.regressions.R") ##per fer multiples regressions d'un dataframe
source("src/comparing.reg.R") ##per comparar les taxes de creiement entre 4 i 3 punts (les regressions)

#import data
#asv_tab pseudoabundances (recalculated with pseudoabundances from flow citometry data)
pseudoabund_df_wide_fc_gtbd <- read.table("data/intermediate_files/regressions_test_datasets/psueodabund_df_fc_wide_gtbd.txt", 
                                     header = TRUE, sep = "\t")
pseudoabund_df_wide_fc_silva <- read.table("data/intermediate_files/regressions_test_datasets/psueodabund_df_fc_wide_silva.txt", 
                                     header = TRUE, sep = "\t")
#objecte de phyloseq del que traiem la metadata i la tax
rem_fc_filt_GTBD <- readRDS("data/intermediate_files/remei_phyloseq_gtbd_fc.rds") 
rem_fc_filt_SILVA <- readRDS("data/intermediate_files/remei_phyloseq_silva_fc.rds") 

#substituint silva per GTBD o al revés fas els càlculs d'un o l'altre.
#metadata 
env <- rem_fc_filt_SILVA@sam_data
env$treatment <- env$treatment %>% 
  factor(levels=(c("CL", "CD", "PL", "PD", "DL", "VL")))
env$season <- env$season %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Early_fall", "Fall")))

#asv_tab & metadata
pseudoabund_df_wide_m <- pseudoabund_df_wide_fc_silva %>%
  merge(env, by = c("treatment", "season", "time", "replicate", "hours.x"))

#asv's taxonomy this new column is important for ploting asv with their taxonomic information
rem_fc_filt_SILVA@tax_table <- rem_fc_filt_SILVA@tax_table %>% 
  mutate_tax_table(asv_num = str_c( 'asv' , 1:ncol(rem_fc_filt_SILVA@otu_table)))
tax_table <- rem_fc_filt_SILVA@tax_table %>%
  mutate_tax_table(tax_ed = (paste(asv_num, class, family, genus, species, sep = "_")))

#1. Prepare  datasets for calculating regressions------
###crear una nova variable que combina time and replicate per no tenir problemes amb els temps repetits de les rèpliques.
data <- pseudoabund_df_wide_m %>% #pseudoabundances dataset form fc data in wide format with metadata.
  mutate(t_rep = paste(time, replicate, sep="_"))

env <- env %>%
  as_tibble() %>%
  mutate(t_rep = paste(time, replicate, sep="_"))

## 2. Filtering for each treatment, season and time
asv_filt_t023 <- filter.season.treatment.time(data = data,  
                                              treatment_col = treatment, treatment_keep = "VL", 
                                              season_col = season, season_keep = "Winter",
                                              time_col = time, time_keep = c("t0", "t2", "t3")) %>%
  select(starts_with("asv"), #we clean metadata columns that we don't need now and keep only the column t_rep
         matches(c("t_rep", "hours.x"))) %>%
  mutate(across(everything(), ~replace(., . == 0, "1"))) %>% ##0 values will be substituted by 1 because it is mandatory for transforming pseudoabundances to log
  mutate(across(!t_rep, as.numeric))
# 
#  asv_filt_t023 %>%
#     glimpse() #las dos ultimas columnas són t_rep i hours.x 

## 3. Calculating regressions (CHANGE MATRIX DIM FOR gtbd: 3810, SILVA: 4594)
regress <- apply(X = asv_filt_t023[,c(1:4594)], MARGIN = 2, FUN = multiple.linear.regressions, env = asv_filt_t023$hours.x) ##DIMENSIONS OF THE MATRIX SHOULD BE ADAPTED TO DATA!!!!!
regress <- as.data.frame(t(regress)) %>% # le damos la vuelta y convertimos en data frame
  rownames_to_column(var = "asv_num")

# resultado %T>% #si surten NAs en el resultat és perquè les dimensions a funció d'apply están malament.
#   dim() %>%
#   head()

#cambiar nombre  del objeto en función del tratamiento y la estación
#NO ENTENC QUÊ FAIG AQUÍ ara no tinc el dataset asv_m perquè l'he eliminat. Com necessito que sigui per fer els CLots?
##vull un dataset amb tota la informació important el problema és que em quedo sense memòria de l'ordinador haig de simCLificar el procés
asv_m_reg_Winter_VL_t023 <- filter.season.treatment.time(data = env,  
                               treatment_col = treatment, treatment_keep = "VL", 
                               season_col = season, season_keep = "Winter",
                               time_col = time, time_keep = c("t0", "t2", "t3")) %>%
  select(starts_with("asv"), ".sample", "treatment", "replicate",  "time", "season", "sample_name", "sample_code", "selected_file_name",
         "reads", "light_regime.y", "hours.y", "LNA", "HNA", "All", "BM", "Leu.PB", "SGR", "TD", "t_rep")  %>%
  left_join(asv_filt_t023, by = "t_rep") %>%
  pivot_longer(cols = starts_with("asv"), names_to = "asv_num", values_to = "pseudoabundance") %>%
  left_join(regress, by = "asv_num") %>%
  add_column(regression_times = "t023")

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

##si vull fer totes les regressions amb diferents temps el que hauria de fer és afegir una nova columna que digui 
## de quin temps a quin temps es correspon cada regressió.

##unir tots els datasets
reg_all_t0234_silva <- bind_rows(reg_winter, 
                     reg_spring,
                     reg_summer,
                     reg_early_fall,
                     reg_fall)
## TINC UN PROBLEMA AMB PL-SPRING PERQUÈ LES REG HAN ESTAT CALCULADES DUES VEGADES I ESTAN DUPLICADES EN EL DATASET
#per això hi ha un pas on les elimino

#regressions 4 punts
write.csv(reg_all_t0234_silva, "data/intermediate_files/asv_reg_all_t0234_silva.csv")
#regressions 3 punts
write.table(reg_all_t0234_silva, "data/intermediate_files/asv_reg_all_t0234_silva.txt", sep = "\t")

##Comparison columna slope 4 temps i 3 temps, decidim quina ens quedem per cada asv ------
##import complete regression datasets 
###GTBD
reg_all_t0234_gtbd <- read.csv("data/intermediate_files/regressions_test_datasets/asv_reg_all_t0234_gtbd.csv", sep="\t", header = T) %>%
  as_tibble()
reg_all_t023_gtbd <- read.csv("data/intermediate_files/regressions_test_datasets/asv_reg_all_t023_gtbd.csv", sep="\t", header = T) %>%
  as_tibble()
###SILVA
reg_all_t0234_silva <-  read.table("data/intermediate_files/asv_reg_all_t0234_silva.txt", sep="\t", header = T) %>%
  as_tibble()
reg_all_t023_silva <-  read.table("data/intermediate_files/asv_reg_all_t023_silva.txt", sep="\t", header = T) %>%
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

#df1 ha de ser 3 temps i df2 4 temps pel problema del PL spring (sol·lucionar-ho)
reg_all_slopes_chosen_silva <- comparing.reg(df1 = reg_all_t023_silva, 
                                                        df2 = reg_all_t0234_silva,
                                                        treatment = treatment,
                                                        season = season,
                                                        asv_num = asv_num, 
                                                        slope = slope, 
                                                        slope.pval = slope.pval)


#Afegim taxonomia i fem boxplots generals per tenir un overview de la coherència que hi ha de les ----
##taxes de creiement per grups taxonomics
reg_all_slopes_chosen_silva_tax <- reg_all_slopes_chosen_silva %>%
  left_join(tax_table, by = "asv_num", copy = TRUE) 

write.csv(reg_all_slopes_chosen_silva_tax, "data/intermediate_files/reg_all_slopes_chosen_silva_tax.csv")

#es podria afegir la rel_abund per filtrar per exemple pels més abundants en general dels experiments.

## Plotting boxplots for general exploration of the dataset------------
## Import data
### Cambiar GTBD por SILVA en funció del dataset amb el que estigui treballant.
reg_all_slopes_chosen_gtbd_tax <- read.csv("data/intermediate_files/reg_all_slopes_chosen_gtbd_tax.csv", sep=",") %>%
  filter(season != "Early_fall")
reg_all_slopes_chosen_silva_tax <- read.csv("data/intermediate_files/reg_all_slopes_chosen_silva_tax.csv", sep=",") %>%
  filter(season != "Early_fall")

##filter by asv present at least >1% relative abundance at some treatment or station or replicate
#rem_fc_relabun <- read_rds("data/rem_fc_relabun.rds")
rem_relabun_melt <-  read.table("data/rem_relabun_melt.txt", sep="\t", header = TRUE) %>%
  filter(season != "Early_fall")

rem_relabun_melt %$%
  asv_num %>%
  unique() #4594 asv al meu dataset

abundant_asv <- rem_relabun_melt %>% 
  filter(Abundance > 0.01) %>% #more than 1% of the community at some point
  dplyr::select(asv_num) %>%
  unique() %>%
  as_tibble()

# rem_relabun_melt_1perc <- rem_relabun_melt %>%
#   right_join(abundant_asv, by = "asv_num", copy = TRUE) #per comprovar amb quin % de relabund estic treballant

# reg_all_slopes_chosen_silva_tax_1perc %>%
#   colnames()

reg_all_slopes_chosen_silva_tax_1perc <- reg_all_slopes_chosen_silva_tax %>% 
  right_join(abundant_asv, by = "asv_num", copy = TRUE)

##reorder treatments and seasons for plots
reg_all_slopes_chosen_silva_tax_1perc$treatment <- reg_all_slopes_chosen_silva_tax_1perc$treatment %>% 
  factor(levels=(c("CL", "CD", "PL", "PD", "DL", "VL")))
reg_all_slopes_chosen_silva_tax_1perc$season <- reg_all_slopes_chosen_silva_tax_1perc$season %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))

##transform slopes /hours to /days
reg_all_slopes_chosen_silva_tax_1perc <- reg_all_slopes_chosen_silva_tax_1perc %>%
  mutate(slope_chosen_days = slope_chosen*24)

##Single ASV based growth rates (FIGURE 1)----
reg_all_slopes_chosen_silva_tax %>%
  colnames()

#statistics seasons----
reg_all_slopes_chosen_silva_tax_1perc %>%
  colnames()
anova<-aov(slope_chosen_days ~season, data = reg_all_slopes_chosen_silva_tax_1perc)
summary(anova) #p<0.05 ***
TukeyHSD(anova)#para ver los grupos que son significativamente distintos
aov_residuals <- residuals(object = anova)
# plot(anova, 1)
# plot(anova, 2)
shapiro.test(x = aov_residuals )#From the output, the p-value > 0.05 implying that the distribution of the data are 
#not significantly different from normal distribution. In other words, we can assume the normality.
##p-value < 2.2e-16 NO és normal 
library(FSA)##test the d
dunnTest(slope_chosen_days ~ season, data = reg_all_slopes_chosen_silva_tax_1perc,
         method= 'bonferroni')
results<-dunnTest(slope_chosen_days ~ season, data = reg_all_slopes_chosen_silva_tax_1perc,
                  method='bonferroni')$res
results<-results[results$P.adj<0.05,]
#Differents: fall-spring, fall-summer, fall-winter, spring-winter, summer-winter
X = results$P.adj <= 0.05
names(X) = gsub(" ",  "",  results$Comparison)
multcompLetters(X)

# Fall Spring Summer Winter 
# "a"    "b"    "b"    "c" 
#statistical groups for non-parammetric test
letters_seas <- data.frame(multcompLetters(X)['Letters'])

#plot seasons violin plot----
seas <- reg_all_slopes_chosen_silva_tax_1perc %>%
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
  geom_text(data = letters_seas, aes(y = 9, x = row.names(letters_seas), label = Letters),
            position = position_nudge(x = 0.2), hjust = 0.7, color = "black")+
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

#statistics treatments----
##Treatments
anova<-aov(slope_chosen_days ~ treatment, data = reg_all_slopes_chosen_silva_tax_1perc)
summary(anova) #p<0.05 ***
TukeyHSD(anova)#para ver los grupos que son significativamente distintos
aov_residuals <- residuals(object = anova)
# plot(anova, 1)
# plot(anova, 2)
shapiro.test(x = aov_residuals )#From the output, the p-value > 0.05 implying that the distribution of the data are #not significantly different from normal distribution. In other words, we can assume the normality.
##p-value < 2.2e-16 NO és normal 
#si no es normal:
#hago test no parametrico
kruskal.test(slope_chosen_days ~ treatment, data = reg_all_slopes_chosen_silva_tax_1perc) #
#si sale p<0.05 hago dunn test para ver cuales son significativamente distintos
library(FSA)##test the d
dunnTest(slope_chosen_days ~ treatment, data = reg_all_slopes_chosen_silva_tax_1perc,
         method= 'bonferroni')
results<-dunnTest(slope_chosen_days ~ treatment, data = reg_all_slopes_chosen_silva_tax_1perc,
                  method='bonferroni')$res
#results<-results[results$P.adj<0.05,]
##Diferents: CD-PD, PD-VL
X = results$P.adj <= 0.05

names(X) = gsub(" ",  "",  results$Comparison)
multcompLetters(X)
# CD   CL   DL   PD   PL   VL 
# "a" "ab" "ab"  "b" "ab"  "a"
letters_treat <- data.frame(multcompLetters(X)['Letters'])

treat <- reg_all_slopes_chosen_silva_tax_1perc %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  filter(season != "Early_fall") %>%
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
  geom_text(data = letters_treat, aes(y = 9, x = row.names(letters_treat), label = Letters),
            position = position_nudge(x = 0.2), hjust = 0.7, color = "black")+
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

 gr_asv_dapis <- grid.arrange(treat, gr_dapis_treat, seas, gr_dapis_seas)
# 
 ggsave('GR_asv_dapis_stat.pdf', gr_asv_dapis, path = "~/Documentos/Doctorat/REMEI/results/figures/raw/",
       height = 175,
       width = 250,
       units = 'mm')

##Mean growth rate by season and treatment considering non-growing prokaryotes (discarded figure)-----
seas_bulk_com <- reg_all_slopes_chosen_silva_tax_1perc %>%
  filter(season != "Early_fall") %>%
  mutate(slope_chosen_days_all = case_when(is.na(slope_chosen_days) ~ '0',
                                           slope_chosen_days < 0 ~ '0',
                                           pvalue_slope_chosen > 0.05 ~ '0',
                                           TRUE ~ as.character(slope_chosen_days))) %>%
  ggplot(aes(season, as.numeric(slope_chosen_days_all)))+ #
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

treat_bulk_com <- reg_all_slopes_chosen_silva_tax_1perc %>%
  filter(season != "Early_fall") %>%
  mutate(slope_chosen_days_all = case_when(is.na(slope_chosen_days) ~ '0',
                                           slope_chosen_days < 0 ~ '0',
                                           pvalue_slope_chosen > 0.05 ~ '0',
                                           TRUE ~ as.character(slope_chosen_days))) %>%
  ggplot(aes(treatment, as.numeric(slope_chosen_days_all)))+ #
  geom_point(aes(color = season), alpha = 0.7, position = position_jitter(0.25))+
  #geom_boxplot(alpha= 0.5)+
  geom_violin(alpha = 0.2, draw_quantiles = c(0.25, 0.75) )+
  labs(color = "Season")+
  stat_summary(fun = "mean", geom = "crossbar", 
               width = 0.6, colour = "black")+
  # geom_text(data = letters_seas, aes(y = 9, x = row.names(letters_seas), label = Letters),
  #           position = position_nudge(x = 0.2), hjust = 0.7, color = "black")+
  labs(x= "Treatment", y = expression("Growth rate day"^"-1"))+  #(μ) 
  #scale_x_discrete()+
  scale_y_continuous(expand = c(0, 0), limits = c(0,10.4))+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_seasons_4)+
  #facet_wrap(.~phylum, scales = "free")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

##mean DAPIS vs mean ASV (filtered & bulk community) ------
general_mean_gr_asv <- reg_all_slopes_chosen_silva_tax_1perc %>%
  filter(season != "Early_fall") %>%
  mutate(slope_chosen_days_all = case_when(is.na(slope_chosen_days) ~ '0',
                                           slope_chosen_days < 0 ~ '0',
                                           pvalue_slope_chosen > 0.05 ~ '0',
                                           TRUE ~ as.character(slope_chosen_days))) %>%
  #group_by(treatment, season) %>%
  mutate(mean_asv = mean(as.numeric(slope_chosen_days_all))) %>%
  distinct(mean_asv)

mean_gr_asv <- reg_all_slopes_chosen_silva_tax %>%
  filter(season != "Early_fall") %>%
  mutate(slope_chosen_days_all = case_when(is.na(slope_chosen_days) ~ '0',
                                           slope_chosen_days < 0 ~ '0',
                                           pvalue_slope_chosen > 0.05 ~ '0',
                                           TRUE ~ as.character(slope_chosen_days))) %>%
  group_by(treatment, season) %>%
  mutate(mean_asv = mean(as.numeric(slope_chosen_days_all))) %>%
  distinct(season, treatment, mean_asv)

##table mean gw general community with 0 (ho faig amb tota la comunitat general o només amb 1perc com en els gràfics)
mean_gr_asv <- reg_all_slopes_chosen_silva_tax_1perc %>%
  filter(season != "Early_fall") %>%
  mutate(slope_chosen_days_all = case_when(is.na(slope_chosen_days) ~ '0',
                                           slope_chosen_days < 0 ~ '0',
                                           pvalue_slope_chosen > 0.05 ~ '0',
                                           TRUE ~ as.character(slope_chosen_days))) %>%
  group_by(season) %>%
  mutate(mean_asv = mean(as.numeric(slope_chosen_days_all)),
         sd_asv = sd(as.numeric(slope_chosen_days_all))) %>%
  distinct(season, mean_asv, sd_asv)


general_mean_gr_dapis <- GR_dapis_OS_filt %>%
  #group_by(treatment, season) %>%
  mutate(mean_dapis = mean(as.numeric(PRO))) %>%
  distinct(mean_dapis)

mean_gr_dapis <- GR_dapis_OS_filt %>%
  group_by(treatment, season) %>%
  mutate(mean_dapis = mean(as.numeric(PRO))) %>%
  distinct(season, treatment, mean_dapis)

mean_gr_asv_dapis <-  mean_gr_asv %>%
  left_join(mean_gr_dapis)

rel_gr_dapis_asv_bulk <- mean_gr_asv_dapis %>%
  as_tibble() %>%
  ggplot(aes(mean_asv, mean_dapis))+
  #geom_smooth(method = 'lm', se = FALSE, color = 'grey')+
  geom_point(aes(shape = treatment, color = season, size = 2))+
  #scale_fill_manual(values = palette_seasons_4)+
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  stat_poly_eq(color = 'black')+
  guides(size = 'none')+
  labs(x = 'Mean growth rates ASV based \n
       (whole community)', y = 'Mean growth rates DAPI based', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (3/4), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

##asv gr mean filtered by negative values, NAs and unsignificant
general_mean_gr_asv_filt <- reg_all_slopes_chosen_silva_tax_1perc %>%
  filter(season != "Early_fall") %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  #group_by(treatment, season) %>%
  mutate(mean_asv = mean(as.numeric(slope_chosen_days))) %>%
  distinct(mean_asv)

mean_gr_asv_filt <- reg_all_slopes_chosen_silva_tax_1perc %>%
  filter(season != "Early_fall") %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  group_by(treatment, season) %>%
  mutate(mean_asv = mean(as.numeric(slope_chosen_days))) %>%
  distinct(season, treatment, mean_asv)

mean_gr_asv_filt_dapis <-  mean_gr_asv_filt %>%
  left_join(mean_gr_dapis)

rel_gr_dapis_asv <- mean_gr_asv_filt_dapis %>%
  as_tibble() %>%
  ggplot(aes(mean_asv, mean_dapis))+
  #geom_smooth(method = 'lm', se = FALSE, color = 'grey')+
  geom_point(aes(shape = treatment, color = season, size = 2))+
  #scale_fill_manual(values = palette_seasons_4)+
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
  scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  #stat_poly_eq(color = 'black')+
  stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
  guides(size = 'none')+
  labs(x = 'Mean growth rates ASV based \n
       (growing community)', y = 'Mean growth rates DAPI based', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/3), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

##general means graphic
cbind(general_mean_gr_asv, general_mean_gr_asv_filt, general_mean_gr_dapis) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = 'mean_gr') %>%
  as_tibble() %>%
  ggplot(aes(mean_gr, V1))+
  geom_point()+
  theme_bw()

##FIGURE 1 panel construction------
gr_asv_community <- grid.arrange(treat_bulk_com, seas_bulk_com)
gr_asv_dapis_relations <- grid.arrange(rel_gr_dapis_asv, rel_gr_dapis_asv_bulk, nrow = 1)
 
gr_asv_dapis_relations <- grid.arrange(treat, seas, gr_dapis_treat, gr_dapis_seas,  treat_bulk_com, seas_bulk_com,
                                        rel_gr_dapis_asv, rel_gr_dapis_asv_bulk, ncol = 2)
 
ggsave('gr_asv_community.pdf', gr_asv_community, path = "~/Documentos/Doctorat/REMEI/results/figures/raw/",
       height = 200,
       width = 200,
       units = 'mm')

ggsave('gr_asv_dapis_relations2.pdf', gr_asv_dapis_relations, path = "~/Documentos/Doctorat/REMEI/results/figures/raw/",
       height = 125,
       width = 200,
       units = 'mm')

ggsave('gr_asv_dapis.pdf', rel_gr_dapis_asv, path = "~/Documentos/Doctorat/REMEI/results/figures/raw/",
       height = 125,
       width = 125,
       units = 'mm')


#Box plots----
#Single ASV based mean growth rate plots by different taxonomic ranks (Suppelementary)-------
reg_all_slopes_chosen_silva_tax %>%
  filter(slope.pval.x >0.05 & slope_chosen_days > 0) %>%  ##hauria de ser < 0.05
  summarize(slope_chosen_days_mean = mean(slope_chosen_days), #0.89
            slope_chosen_days= sd(slope_chosen_days)) #0.78
#mean GR by class
reg_all_slopes_chosen_silva_tax %>%
  filter(slope.pval.x >0.05 & slope_chosen_days > 0) %>% 
  group_by(class) %>%
  summarize(slope_chosen_days_mean = mean(slope_chosen_days), 
            slope_chosen_days= sd(slope_chosen_days))

##Order axis x by higher taxonomy levels (of example: order Classes by Phylum) 1%-----
reg_all_slopes_chosen_silva_tax_1perc %>%
  glimpse()

reg_all_slopes_chosen_silva_tax_1perc <- reg_all_slopes_chosen_silva_tax_1perc %>%
  mutate(domain_f = as_factor(domain),
         phylum_f = as_factor(phylum),
         class_f = as_factor(class),
         order_f = as_factor(order),
         family_f = as_factor(family),
         genus_f = as_factor(genus),
         asv_f = as_factor(asv_num))

reg_all_slopes_chosen_silva_tax_1perc %>%
  colnames()

reg_all_slopes_chosen_silva_tax_1perc$phylum_f <-  factor(reg_all_slopes_chosen_silva_tax_1perc$phylum_f, 
                                                         levels=unique(
                                                           reg_all_slopes_chosen_silva_tax_1perc$phylum_f[
                                                             order(reg_all_slopes_chosen_silva_tax_1perc$domain_f)]), 
                                                         ordered=TRUE)

reg_all_slopes_chosen_silva_tax_1perc$class_f <-  factor(reg_all_slopes_chosen_silva_tax_1perc$class_f, 
                                                          levels=unique(
                                                            reg_all_slopes_chosen_silva_tax_1perc$class_f[
                                                              order(reg_all_slopes_chosen_silva_tax_1perc$domain_f,
                                                                    reg_all_slopes_chosen_silva_tax_1perc$phylum_f)]), 
                                                          ordered=TRUE)
reg_all_slopes_chosen_silva_tax_1perc$order_f <-  factor(reg_all_slopes_chosen_silva_tax_1perc$order_f, 
                                                         levels=unique(
                                                           reg_all_slopes_chosen_silva_tax_1perc$order_f[
                                                             order(reg_all_slopes_chosen_silva_tax_1perc$domain_f,
                                                                   reg_all_slopes_chosen_silva_tax_1perc$phylum_f,
                                                                   reg_all_slopes_chosen_silva_tax_1perc$class_f)]), 
                                                         ordered=TRUE)

reg_all_slopes_chosen_silva_tax_1perc$family_f <-  factor(reg_all_slopes_chosen_silva_tax_1perc$family_f, 
                                                         levels=unique(
                                                           reg_all_slopes_chosen_silva_tax_1perc$family_f[
                                                             order(reg_all_slopes_chosen_silva_tax_1perc$domain_f,
                                                                   reg_all_slopes_chosen_silva_tax_1perc$phylum_f,
                                                                   reg_all_slopes_chosen_silva_tax_1perc$class_f,
                                                                   reg_all_slopes_chosen_silva_tax_1perc$order_f)]), 
                                                         ordered=TRUE)

reg_all_slopes_chosen_silva_tax_1perc$genus_f <-  factor(reg_all_slopes_chosen_silva_tax_1perc$genus_f, 
                                                          levels=unique(
                                                            reg_all_slopes_chosen_silva_tax_1perc$genus_f[
                                                              order(reg_all_slopes_chosen_silva_tax_1perc$domain_f,
                                                                    reg_all_slopes_chosen_silva_tax_1perc$phylum_f,
                                                                    reg_all_slopes_chosen_silva_tax_1perc$class_f,
                                                                    reg_all_slopes_chosen_silva_tax_1perc$order_f,
                                                                    reg_all_slopes_chosen_silva_tax_1perc$family_f)]), 
                                                          ordered=TRUE)
reg_all_slopes_chosen_silva_tax_1perc$asv_f <-  factor(reg_all_slopes_chosen_silva_tax_1perc$asv_f, 
                                                         levels=unique(
                                                           reg_all_slopes_chosen_silva_tax_1perc$asv_f[
                                                             order(reg_all_slopes_chosen_silva_tax_1perc$domain_f,
                                                                   reg_all_slopes_chosen_silva_tax_1perc$phylum_f,
                                                                   reg_all_slopes_chosen_silva_tax_1perc$class_f,
                                                                   reg_all_slopes_chosen_silva_tax_1perc$order_f,
                                                                   reg_all_slopes_chosen_silva_tax_1perc$family_f,
                                                             reg_all_slopes_chosen_silva_tax_1perc$genus_f)]), 
                                                         ordered=TRUE)

##dataset without filtering by 1%------
reg_all_slopes_chosen_silva_tax <- reg_all_slopes_chosen_silva_tax %>%
  mutate(domain_f = as_factor(domain),
         phylum_f = as_factor(phylum),
         class_f = as_factor(class),
         order_f = as_factor(order),
         family_f = as_factor(family),
         genus_f = as_factor(genus),
         asv_f = as_factor(asv_num))


reg_all_slopes_chosen_silva_tax %>%
  colnames()

reg_all_slopes_chosen_silva_tax$phylum_f <-  factor(reg_all_slopes_chosen_silva_tax$phylum_f, 
                                                          levels=unique(
                                                            reg_all_slopes_chosen_silva_tax$phylum_f[
                                                              order(reg_all_slopes_chosen_silva_tax$domain_f)]), 
                                                          ordered=TRUE)
reg_all_slopes_chosen_silva_tax$phylum_f %>% levels()

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

reg_all_slopes_chosen_silva_tax$phylum_f %>%
  levels()
##exploration at class level---
### el vull al supplementary
#by seasons
plot_class_seasons <- reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  group_by(class_f) %>%
  filter(n() > 20) %>%
  #filter(pseudoabundance > 1000) %>%
  ggplot(aes(x=interaction(class_f,phylum_f, sep = '\n'), slope_chosen_days, group = class_f))+ #
  geom_point(aes(color = season), alpha = 0.8, position = position_jitter(0.4))+#
  geom_boxplot(alpha= 0.1)+
  labs(color = "Season")+
  #geom_violin(alpha = 0.5, draw_quantiles = c(0.25, 0.75))+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.3,
               colour = "black")+
  labs(x= "Class\n Pylum", y = expression("Growth rate day"^"-1"))+#(μ)
  # geom_label(position = position_stack(vjust = 0.99))+
  #geom_text(position = position_dodge(width = 1), aes(x=class, y=-50))+
  #scale_y_continuous(limits = c(0, 10))+
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_seasons_4)+#palf_large(94)+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 10.4))+
  #coord_flip()+
  #facet_wrap(.~treatment, scales = "fixed")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 4), 
        axis.text.x = element_text(angle = 0, size = 5), 
        axis.title = element_text(size = 7),
        legend.position = "right",
        axis.text.y = element_text(angle = 0, size = 5),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 7),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        legend.margin=margin(0,0,0,0.15))

ggsave('boxplot_class_seasons.pdf', plot_class_seasons, path = "~/Documentos/Doctorat/REMEI/results/figures/raw/",
       height = 150,
       width = 180,
       units = 'mm')

#by treatmens
plot_class_treatments <- reg_all_slopes_chosen_silva_tax_1perc %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  group_by(class_f) %>%
  filter(n() > 2) %>%
  #filter(pseudoabundance > 1000) %>%
  ggplot(aes(x=interaction(class_f,phylum_f, sep = '\n'), slope_chosen_days, group = class_f))+ #
  geom_point(aes(color = treatment), alpha = 0.8, position = position_jitter())+#
  geom_boxplot(alpha= 0.3)+
  labs(color = "Treatment")+
  #geom_violin(alpha = 0.5, draw_quantiles = c(0.25, 0.75))+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.3,
               colour = "black")+
  labs(x= "Class\n Pylum", y = expression("Growth rate day"^"-1"))+#(μ)
  # geom_label(position = position_stack(vjust = 0.99))+
  #geom_text(position = position_dodge(width = 1), aes(x=class, y=-50))+
  #scale_y_continuous(limits = c(0, 10))+
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_treatments_remei)+#palf_large(94)+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 10.4))+
  #coord_flip()+
  #facet_wrap(.~treatment, scales = "fixed")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 4), 
        axis.text.x = element_text(angle = 0, size = 5), 
        axis.title = element_text(size = 7),
        legend.position = "right",
        axis.text.y = element_text(angle = 0, size = 5),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 7),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        legend.margin=margin(0,0,0,0.15))

ggsave('boxplot_class_treatments_1perc.pdf', plot_class_treatments, path = "~/Documentos/Doctorat/REMEI/results/figures/raw/",
       height = 150,
       width = 180,
       units = 'mm')


###prova sense subset per 1% per plotejant més punts----
reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  #filter(pseudoabundance > 1000) %>%
  ggplot(aes(x=interaction(class_f,phylum_f, sep = ' '), slope_chosen_days, group = class_f))+ #
  geom_point(aes(color = treatment), alpha = 0.8, position = position_jitter())+#
  geom_boxplot(alpha= 0.3)+
  labs(color = "Treatment")+
  #geom_violin(alpha = 0.5, draw_quantiles = c(0.25, 0.75))+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.3,
               colour = "black")+
  labs(x= "Class\n Pylum", y = expression("Growth rate (μ) day"^"-1"))+
  # geom_label(position = position_stack(vjust = 0.99))+
  #geom_text(position = position_dodge(width = 1), aes(x=class, y=-50))+
  #scale_y_continuous(limits = c(0, 10))+
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_treatments_remei)+#palf_large(94)+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 10.4))+
  coord_flip()+
  #coord_flip()+
  #facet_wrap(.~treatment, scales = "fixed")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 4), 
        axis.text.x = element_text(angle = 0, size = 5), 
        axis.title = element_text(size = 7),
        legend.position = "right",
        axis.text.y = element_text(angle = 0, size = 5),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 7),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        legend.margin=margin(0,0,0,0.15))

##plot at order level-----
### el vull al supplementary
###by seasons
plot_order_seasons <- reg_all_slopes_chosen_silva_tax_1perc %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  group_by(order) %>%
  filter(n() > 2) %>%
  #filter(pseudoabundance > 1000) %>%
  ggplot(aes(interaction(order_f, class_f, phylum_f, sep = '\n'), slope_chosen_days))+ #
  geom_point(aes(color = season), alpha = 0.7, position = position_jitter(0.25))+
  geom_boxplot(alpha= 0.3)+
  #geom_violin(alpha = 0.3)+
  labs(x= "Order", y = expression("Growth rate day"^"-1"), color = 'Season')+# (μ) 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 10.4))+
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_seasons_4)+#palf_large(94)
  #facet_wrap(.~class, scales = "free")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(axis.text.x = element_text(angle = 60, size = 6, hjust = 0.9), 
        strip.text.x = element_text(size = 4), 
        axis.title = element_text(size = 7),
        legend.position = "right",
        axis.text.y = element_text(angle = 0, size = 5),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 7),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        legend.margin=margin(0,0,0,0.1))

ggsave('boxplot_order_seasons_1perc.pdf', plot_order_seasons, path = "~/Documentos/Doctorat/REMEI/results/figures/raw/",
       height = 150,
       width = 180,
       units = 'mm')

###by treatments
plot_order_treatments <- reg_all_slopes_chosen_silva_tax_1perc %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  group_by(order) %>%
  filter(n() > 2) %>%
  #filter(pseudoabundance > 1000) %>%
  ggplot(aes(interaction(order_f, class_f,phylum_f, sep = '\n'), slope_chosen_days))+ #
  geom_point(aes(color = treatment), alpha = 0.7, position = position_jitter())+
  geom_boxplot(alpha= 0.3)+
  #geom_violin(alpha = 0.5)+
  labs(x= "Order", y = expression("Growth rate day"^"-1"), color = 'Season')+# (μ) 
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_treatments_remei)+#palf_large(94)
  scale_y_continuous(expand = c(0, 0), limits = c(0, 10.4))+
  #facet_wrap(.~class, scales = "free")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(axis.text.x = element_text(angle = 60, size = 6, hjust = 0.9), 
        strip.text.x = element_text(size = 4), 
        axis.title = element_text(size = 7),
        legend.position = "right",
        axis.text.y = element_text(angle = 0, size = 5),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 7),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        legend.margin=margin(0,0,0,0.1))

ggsave('boxplot_order_treatments_1perc.pdf', plot_order_treatments, path = "~/Documentos/Doctorat/REMEI/results/figures/raw/",
       height = 150,
       width = 180,
       units = 'mm')

##at family level boxplots ----
###by seasons
reg_all_slopes_chosen_silva_tax_1perc %>%
  colnames()
plot_family_seasons <- reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  group_by(family) %>%
  filter(n() > 50) %>%
  #filter(pseudoabundance > 1000) %>%
  ggplot(aes(interaction(family_f, order_f, phylum_f, sep = '\n'), slope_chosen_days))+ #
  geom_point(aes(color = season), alpha = 0.7, position = position_jitter())+
  geom_boxplot(alpha= 0.3)+
  #geom_violin(alpha = 0.5)+
  labs(x= "Order", y = expression("Growth rate (μ) day"^"-1"), color = 'Season')+
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_seasons_4)+#palf_large(94)
  scale_y_continuous(expand = c(0, 0), limits = c(0, 10.4))+
  #facet_wrap(.~class, scales = "free")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(axis.text.x = element_text(angle = 60, size = 6, hjust = 0.9), 
        strip.text.x = element_text(size = 4), 
        axis.title = element_text(size = 7),
        legend.position = "right",
        axis.text.y = element_text(angle = 0, size = 5),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 7),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        legend.margin=margin(0,0,0,0.1))

ggsave('boxplot_family_seasons.pdf', plot_family_seasons, path = "~/Documentos/Doctorat/REMEI/results/figures/raw/",
       height = 150,
       width = 180,
       units = 'mm')

###by treatments
plot_family_treatments <- reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  group_by(family) %>%
  filter(n() > 50) %>%
  #filter(pseudoabundance > 1000) %>%
  ggplot(aes(interaction(family_f, order_f, phylum_f, sep = '\n'), slope_chosen_days))+ #
  geom_point(aes(color = treatment), alpha = 0.7, position = position_jitter())+
  geom_boxplot(alpha= 0.3)+
  #geom_violin(alpha = 0.5)+
  labs(x= "Order", y = expression("Growth rate day"^"-1"), color = 'Treatment')+#(μ) 
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_treatments_remei)+#palf_large(94)
  scale_y_continuous(expand = c(0, 0), limits = c(0, 10.4))+
  #facet_wrap(.~class, scales = "free")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(axis.text.x = element_text(angle = 60, size = 6, hjust = 0.9), 
        strip.text.x = element_text(size = 4), 
        axis.title = element_text(size = 7),
        legend.position = "right",
        axis.text.y = element_text(angle = 0, size = 5),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 7),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        legend.margin=margin(0,0,0,0.1))

ggsave('boxplot_family_treatments.pdf', plot_family_treatments, path = "~/Documentos/Doctorat/REMEI/results/figures/raw/",
       height = 150,
       width = 180,
       units = 'mm')

##exploration a nivell de classes concretes
#"Acidimiccrobiia", "Actinobaccteria", "Alphaproteobacteria", "Bacilli", "Bacteroidia", "Clostridia", "Cyanobacteriia",
#"Gammaproteobacteria",  "Nitrosphaeria", "Phycisphaerae", "Planctomycetes", "Verrucomicrobiae", "Vicinamibacteria"

reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  filter(class == "Actinobacteria") %>%
  ggplot(aes(fct_reorder(order, family), slope_chosen))+ #
  geom_point(aes(color = genus), alpha = 0.7, position = position_jitter())+
  #geom_boxplot(alpha= 0.5)+
  geom_violin(alpha = 0.5)+
  labs(x= "Class", y = expression("Growth rate (μ) day"^"-1"))+
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_phylums)+
  #facet_wrap(.~phylum, scales = "free")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 90), legend.position = "none")

reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  filter(order == "Micrococcales") %>%
  ggplot(aes(fct_reorder(genus, family), slope_chosen))+ #
  geom_point(aes(color = genus), alpha = 0.7, position = position_jitter())+
  geom_boxplot(alpha= 0.5)+
  #geom_violin(alpha = 0.5)+
  labs(x= "Class", y = expression("Growth rate (μ) day"^"-1"))+
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_phylums)+
  #facet_wrap(.~phylum, scales = "free")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 90), legend.position = "none")

reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  filter(class == "Alphaproteobacteria") %>%
  ggplot(aes(fct_infreq(order), slope_chosen_days))+ #
  geom_point(aes(color = season), alpha = 0.9, position = position_jitter())+
  #geom_boxplot(alpha= 0.5)+
  geom_violin(alpha = 0.1)+
  labs(x= "Class", y = expression("Growth rate day"^"-1"))+#(μ) 
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_seasons_4)+
  #facet_wrap(.~season, scales = "fixed")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 12), 
        axis.text.x = element_text(angle = 60, size = 10, hjust = +1), legend.position = "right")

reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  filter(class == "Gammaproteobacteria") %>%
  group_by(order) %>%
  filter(n() > 2) %>%
  ggplot(aes(fct_infreq(order), slope_chosen_days))+ #
  geom_point(aes(color = season), alpha = 0.9, position = position_jitter())+
  #geom_boxplot(alpha= 0.5)+
  geom_violin(alpha = 0.1)+
  labs(x= "Class", y = expression("Growth rate day"^"-1"), color = 'Season')+#(μ) 
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_seasons_4)+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 10.4))+
  #facet_wrap(.~season, scales = "fixed")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 5), 
        axis.text.x = element_text(angle = 60, 
                                   size = 7, hjust = +1), 
        legend.position = "right",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 5),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        legend.margin=margin(0,0,0,0.1))

reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  filter(family == "Hyphomonadaceae") %>%
  ggplot(aes(asv_num, slope_chosen))+ #
  geom_point(aes(color = season), alpha = 0.9, position = position_jitter())+
  geom_boxplot(alpha= 0.5)+
  #geom_violin(alpha = 0.1)+
  labs(x= "ASV", y = expression("Growth rate day"^"-1"))+#(μ)
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_seasons_4)+
  #facet_wrap(.~phylum, scales = "free")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 5), 
        axis.text.x = element_text(angle = 60, 
                                   size = 7, hjust = +1), 
        legend.position = "right",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 5),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        legend.margin=margin(0,0,0,0.1))

reg_all_slopes_chosen_silva_tax_1perc %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  #filter(order == "Rhodobacterales") %>%
  group_by(family) %>%
  filter(n() > 2) %>%
  ggplot(aes(fct_infreq(family), slope_chosen_days))+ #
  geom_point(aes(color = season), alpha = 0.9, position = position_jitter(width=0.2))+
  geom_boxplot(alpha= 0.5)+
  #geom_violin(alpha = 0.1)+
  labs(x= "ASV number", y = expression("Growth rate day"^"-1"))+#(μ) 
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_seasons_4)+
  #facet_wrap(.~season, scales = "free")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 90), legend.position = "bottom")

reg_all_slopes_chosen_silva_tax_1perc %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  filter(order == "SAR11 clade") %>%
  ggplot(aes(asv_num, slope_chosen_days))+ #
  geom_point(aes(color = season), alpha = 0.9, position = position_jitter(width= 0.2))+
  geom_boxplot(alpha= 0.5)+
  #geom_violin(alpha = 0.1, draw_quantiles = c(0.25, 0.5, 0.75))+ #trabar la manera d'afegir la mitja
  labs(x= "Family", y = expression("Growth rate day"^"-1"))+# (μ) 
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_seasons_4)+
  #facet_wrap(.~season, scales = "free")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 90), legend.position = "none")

reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  filter(order == c("SAR11 clade", "Alteromonadales", "Rhodobacterales")) %>%
  ggplot(aes(treatment, slope_chosen_days))+ #
  geom_point(aes(color = season), alpha = 0.9, position = position_jitter())+
  #geom_boxplot(alpha= 0.5)+
  geom_violin(alpha = 0.1, draw_quantiles = c(0.25, 0.75))+ #trabar la manera d'afegir la mitja
  labs(x= "Treatment", y = expression("Growth rate (μ) day"^"-1"))+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.2,
               colour = "black")+
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_seasons)+
  facet_wrap(.~order, scales = "fixed")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 90), legend.position = "bottom")

reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  filter(order == "SAR11 clade") %>%
  ggplot(aes(treatment, slope_chosen_days))+ #
  geom_point(aes(color = season), alpha = 0.9, position = position_jitter())+
  geom_boxplot(alpha= 0.5)+
  #geom_violin(alpha = 0.1, draw_quantiles = c(0.25, 0.75))+ #trabar la manera d'afegir la mitja
  labs(x= "Treatment", y = expression("Growth rate (μ) day"^"-1"))+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.2,
               colour = "black")+
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_seasons_4)+
  facet_wrap(.~family, scales = "fixed")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 90), legend.position = "bottom")

plot_genus_treatments <- reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  #filter(family == "Flavobacteriaceae") %>%
  filter(genus != is.na(genus)) %>%
  group_by(genus) %>%
  filter(n() > 10) %>%
  ggplot(aes(treatment, slope_chosen_days))+ #
  geom_point(aes(color = season), alpha = 0.9, position = position_jitter(0.25), size = 1)+
  geom_boxplot(alpha= 0.5)+
  #geom_violin(alpha = 0.1, draw_quantiles = c(0.25, 0.75))+ #trabar la manera d'afegir la mitja
  labs(x= "Treatment", y = expression("Growth rate day"^"-1"),
       color = 'Season')+# (μ)
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.2,
               colour = "black")+
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_seasons_4)+
  facet_wrap(.~genus_f, scales = "fixed")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 5), 
        axis.text.x = element_text(angle = 90, size = 5),
        legend.position = "bottom",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 5),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        legend.margin=margin(0,0,0,0.1))

ggsave('boxplot_genus_treatments.pdf', plot_genus_treatments, path = "~/Documentos/Doctorat/REMEI/results/figures/raw/",
       height = 200,
       width = 180,
       units = 'mm')

reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  filter(class == "Bacteroidia") %>%
  ggplot(aes(fct_infreq(family), slope_chosen_days))+ #
  geom_point(aes(color = season), alpha = 0.9, position = position_jitter())+
  geom_boxplot(alpha= 0.5)+
  #geom_violin(alpha = 0.1, draw_quantiles = c(0.25, 0.5, 0.75))+ #trabar la manera d'afegir la mitja
  labs(x= "Family", y = expression("Growth rate (μ) day"^"-1"))+
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_seasons_4)+
  #facet_wrap(.~season, scales = "free")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 90), legend.position = "none")

reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  filter(family == "Flavobacteriaceae") %>%
  ggplot(aes(fct_infreq(genus), slope_chosen))+ #
  geom_point(aes(color = season), alpha = 0.9, position = position_jitter())+
  geom_boxplot(alpha= 0.5)+
  #geom_violin(alpha = 0.1, draw_quantiles = c(0.25, 0.5, 0.75))+ #trabar la manera d'afegir la mitja
  labs(x= "Genus", y = expression("Growth rate (μ) day"^"-1"))+
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_seasons_4)+
  #facet_wrap(.~season, scales = "free")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 90), legend.position = "none")

reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  filter(genus == "NS5 marine group") %>%
  ggplot(aes(asv_num, slope_chosen))+ #
  geom_point(aes(color = season), alpha = 0.9, position = position_jitter())+
  geom_boxplot(alpha= 0.5)+
  #geom_violin(alpha = 0.1, draw_quantiles = c(0.25, 0.5, 0.75))+ #trabar la manera d'afegir la mitja
  labs(x= "Season", y = expression("Growth rate (μ) day"^"-1"))+
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_seasons_4)+
  #facet_wrap(.~season, scales = "free")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 90), legend.position = "none")

reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  filter(class == "Gammaproteobacteria") %>%
  ggplot(aes(fct_infreq(family), slope_chosen_days))+ #
  geom_point(aes(color = season), alpha = 0.9, position = position_jitter())+
  geom_boxplot(alpha= 0.5)+
  #geom_violin(alpha = 0.1, draw_quantiles = c(0.25, 0.5, 0.75))+ #trabar la manera d'afegir la mitja
  labs(x= "Family", y = expression("Growth rate (μ) day"^"-1"))+
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_seasons_4)+
  #facet_wrap(.~season, scales = "fixed")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 90), legend.position = "none")

reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  filter(family == "Alteromonadaceae") %>%
  ggplot(aes(asv_num, slope_chosen_days))+ #
  geom_point(aes(color = season), alpha = 0.9, position = position_jitter())+
  geom_boxplot(alpha= 0.1)+
  #geom_violin(alpha = 0.1, draw_quantiles = c(0.25, 0.5, 0.75))+ #trabar la manera d'afegir la mitja
  labs(x= "ASV", y = expression("Growth rate (μ) day"^"-1"))+
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_seasons_4)+
  #facet_wrap(.~season, scales = "free")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 90), legend.position = "none")

reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  filter(family == "Pseudoalteromonadaceae") %>%
  ggplot(aes(asv_num, slope_chosen_days))+ #
  geom_point(aes(color = season), alpha = 0.9, position = position_jitter())+
  geom_boxplot(alpha= 0.2)+
  #geom_violin(alpha = 0.1, draw_quantiles = c(0.25, 0.5, 0.75))+ #trabar la manera d'afegir la mitja
  labs(x= "ASV", y = expression("Growth rate (μ) day"^"-1"))+
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_seasons_4)+
  #facet_wrap(.~season, scales = "free")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 90), legend.position = "bottom")

reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  filter(phylum == "Cyanobacteria") %>%
  ggplot(aes(season, slope_chosen_days))+ #
  geom_point(aes(color = treatment), alpha = 0.9, position = position_jitter())+
  geom_boxplot(alpha= 0.5)+
  #geom_violin(alpha = 0.1, draw_quantiles = c(0.25, 0.5, 0.75))+ #trabar la manera d'afegir la mitja
  labs(x= "Season", y = expression("Growth rate (μ) day"^"-1"))+
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_treatments_remei)+
  #facet_wrap(.~season, scales = "fixed")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 90), legend.position = "none")

reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  #filter(season ==c("Winter", "Spring", "Summer", "Fall")) %>%
  filter(order %in% c("Alteromonadales", "Campylobacter",
                    "Cellvibrionales", "Flavobacteriales", "Flavobact. NS",
                    "Oceanospirillales", "Rhodobacterales", "SAR11 clade", "SAR86",
                    "Vibrionales")) %>%
  ggplot(aes(treatment, slope_chosen_days))+ #
  geom_point(aes(color = season), alpha = 0.9, position = position_jitter())+
  geom_boxplot(alpha= 0.5)+
  #geom_violin(alpha = 0.1, draw_quantiles = c(0.25, 0.75))+ #trabar la manera d'afegir la mitja
  labs(x= "Season", y = expression("Growth rate (μ) day"^"-1"))+
  # stat_summary(fun = "mean",geom = "crossbar", 
  #              width = 0.2, colour = "black")+
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_seasons_4)+
  facet_grid(.~order, scales = "free")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), 
        axis.text.x = element_text(angle = 90), legend.position = "right")

##Effect of light in likely containing proteorhodopsin----
##NOR5, Rhodobacteraeae
reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  filter(genus == "OM60(NOR5) clade", ) %>%
  #filter(family == c("Halieaceae", "Rhodobacteraceae")) %>%
  ggplot(aes(treatment, slope_chosen_days))+ #
  geom_point(aes(color = treatment), alpha = 0.9, position = position_jitter())+
  geom_boxplot(alpha= 0.5)+
  #geom_violin(alpha = 0.1, draw_quantiles = c(0.25, 0.5, 0.75))+ #trabar la manera d'afegir la mitja
  labs(x= "Treatment", y = expression("Growth rate (μ) day"^"-1"))+
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_treatments_remei)+
  facet_wrap(.~season, scales = "fixed")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 90), legend.position = "none")

reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  #filter(genus == "OM60(NOR5) clade", ) %>%
  filter(family == "Rhodobacteraceae") %>%
  ggplot(aes(treatment, slope_chosen_days))+ #
  geom_point(aes(color = treatment), alpha = 0.9, position = position_jitter())+
  geom_boxplot(alpha= 0.5)+
  #geom_violin(alpha = 0.1, draw_quantiles = c(0.25, 0.5, 0.75))+ #trabar la manera d'afegir la mitja
  labs(x= "Treatment", y = expression("Growth rate (μ) day"^"-1"))+
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_treatments_remei)+
  facet_wrap(.~season, scales = "fixed")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 90), legend.position = "none")

reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  #filter(genus == "OM60(NOR5) clade", ) %>%
  filter(order == "Bacteroidales") %>%
  ggplot(aes(treatment, slope_chosen_days))+ #
  geom_point(aes(color = treatment), alpha = 0.9, position = position_jitter())+
  geom_boxplot(alpha= 0.5)+
  #geom_violin(alpha = 0.1, draw_quantiles = c(0.25, 0.5, 0.75))+ #trabar la manera d'afegir la mitja
  labs(x= "Treatment", y = expression("Growth rate (μ) day"^"-1"))+
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_treatments_remei)+
  #facet_wrap(.~family~genus, scales = "fixed")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 90), legend.position = "none")


#busco gèneres amb major taxa de creixement-----
reg_all_slopes_chosen_silva_tax %>%
  #filter(slope_chosen > 0.2) %>%
  filter(pvalue_slope_chosen < 0.05) %$%
  phylum %>%
  unique

# [1] "Pseudoalteromonas" "Polaribacter"      "Dokdonia"          "Microbacterium"    "Glaciecola"       
# [6] "Vibrio"            "Alteromonas"       "TMED14"            "Rubritalea"        NA                 
# [11] "UBA8087"           "AG-422-B15"   

plot.tax.filtered(data = reg_all_t0234_tax,
                  value_slope_pval = 5, 
                  value_intercept_pval = 5, 
                  tax_level = genus, 
                  select_tax_level  = "Sphingomonas")

##per families
reg_all_t0234_tax %>%
  colnames()

reg_all_slopes_chosen_silva_tax %>%
  filter(slope_chosen_days > 0.2) %>%
  filter(slope.pval.x < 0.05 ) %$% #& intercept.pval < 0.05 (NO ES FILTRA NO SERVEIX)
  order %>%
  unique
# 
# # [1] "Alteromonadaceae"   "Flavobacteriaceae"  "Microbacteriaceae"  "Vibrionaceae"       "Schleiferiaceae"   
# # [6] "Akkermansiaceae"    "Rhodobacteraceae"   "SM1A02"             "AG-422-B15"         "Enterobacteriaceae"  
# 
# plot.tax.filtered_colors <- function(data, value_slope_pval, tax_level, select_tax_level){
#   my_var1 <- rlang::enquo(tax_level)
#   data %>%
#     filter(slope.pval < value_slope_pval ) %>% #& intercept.pval < value_intercept_pval
#     filter(!!my_var1 == select_tax_level) %>%
#     ggplot(aes(hours, pseudoabundance, color = tax_ed, shape = season))+
#     geom_point()+
#     labs(x= "Time (h)", y = "log(pseudoabundance)")+
#     geom_smooth(method = "lm")+
#     scale_color_manual(values = palf_large(255))+
#     facet_grid(.~treatment~season, scales = "fixed")+
#     guides(color=guide_legend(ncol = 1))+
#     theme_bw()+
#     theme(strip.text.x = element_text(size = 14))
# }
# 
# plot.tax.filtered.boxplot <- function(data, value_slope_pval){
#   my_var1 <- data %>%
#     filter(slope.pval < value_slope_pval ) %>% #& intercept.pval < value_intercept_pval
#     #filter(!!my_var1 == select_tax_level) %>%
#     ggplot(aes(family, slope, color = family, shape = season))+
#     geom_boxplot()+
#     labs(x= "Time (h)", y = "log(pseudoabundance)")+
#     geom_smooth(method = "lm")+
#     scale_color_manual(values = palf_large(50))+
#     #facet_grid(.~treatment~season, scales = "fixed")+
#     guides(color=guide_legend(ncol = 1))+
#     theme_bw()+
#     theme(strip.text.x = element_text(size = 14))
# }
# 
# reg_all_t0234_tax %>%
#   colnames()
# plot.tax.filtered.boxplot(data = reg_all_t0234_tax,
#                   value_slope_pval = 0.05)
# 
# reg_all_slopes_chosen_silva_tax %>%
#   filter(slope.pval < 0.05 ) %>% #& intercept.pval < value_intercept_pval
#   filter(slope > 0) %>%
#   #filter(pseudoabundance > 500) %>%
#   ggplot(aes(genus, slope, color = order))+
#   geom_boxplot()+
#   geom_point()+
#   labs(x= "Genus", y = "Slope")+
#   #geom_smooth(method = "lm")+
#   scale_color_manual(values = palf_large(153))+
#   facet_wrap(.~class, scales = "free")+
#   #guides(color=guide_legend(ncol = 1))+
#   theme_bw()+
#   theme(strip.text.x = element_text(size = 5), axis.text.x = element_text(angle = 45), legend.position = "none")
# 
# reg_all_slopes_chosen_silva_tax %>%
#   filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
#   filter(slope_chosen > 0) %>%
#   #filter(pseudoabundance > 1000) %>%
#   ggplot(aes(genus, slope, color = family))+
#   geom_boxplot()+
#   geom_point()+
#   labs(x= "Genus", y = "Slope")+
#   #geom_smooth(method = "lm")+
#   scale_color_manual(values = palf_large(140))+
#   #facet_grid(.~treatment~season, scales = "fixed")+
#   #guides(color=guide_legend(ncol = 1))+
#   theme_bw()+
#   theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 90), legend.position = "none")
# 
# reg_all_slopes_chosen_silva_tax %>%
#   filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
#   filter(slope_chosen > 0) %>%
#   #filter(pseudoabundance > 1000) %>%
#   ggplot(aes(genus, slope_chosen, color = order))+
#   geom_boxplot()+
#   geom_point()+
#   labs(x= "Genus", y = "Slope")+
#   #scale_x_discrete()+
#   #geom_smooth(method = "lm")+
#   scale_color_manual(values = palf_large(201))+
#   facet_wrap(.~phylum, scales = "free")+
#   #guides(color=guide_legend(ncol = 1))+
#   theme_bw()+ 
#   theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 90), legend.position = "none")
# 
# reg_all_slopes_chosen_silva_tax %>%
#   colnames()
# 
# reg_all_slopes_chosen_silva_tax %>%
#   filter(pvalue_slope_chosen < 0.05 ) %$% #& intercept.pval < 0.05 (NO ES FILTRA NO SERVEIX)
#   class %>%
#   unique
# # 
# # "Cyanobacteriia"      "Bacteroidia"         "Alphaproteobacteria" "Gammaproteobacteria" "Verrucomicrobiae"   
# # "SAR324"              "Nitrososphaeria"     "Acidimicrobiia"      "Deinococci"          "Poseidoniia"        
# #  "Planctomycetes"      "Bacteriovoracia"     "Actinomycetia"       "UBA2968"             "WGA-4E"             
# #  "UBA1144"             "Clostridia"          "Phycisphaerae"       "UBA8108"             "Dehalococcoidia"    
# #  "Desulfitobacteriia"  "Bacilli"             "Rhodothermia"        "Nitrospiria"         "UBA9160"            
# # "Vicinamibacteria"    "Marinisomatia"       "Blastocatellia"      "Oligoflexia"         "Nitrospinia"        
# #  "Desulfovibrionia"    "Bradimonadia"        "Fusobacteriia"       "Bin61"               "Binatia"            
# # "Vampirovibrionia"    "Lentisphaeria"       "GCA-002687715"       "Zetaproteobacteria"  "Negativicutes"      
# #  "Kiritimatiellae"     "Gemmatimonadetes"    "Desulfobulbia"       "Myxococcia"          "Desulfuromonadia"   
# #  "Anaerolineae"        "UBA1135"             "Coriobacteriia"
# 
# comp_df_tax$treatment <- comp_df_tax$treatment %>% 
#   factor(levels=(c("CL", "CD", "PL", "PD", "DL", "VL")))
# comp_df_tax$season <- comp_df_tax$season %>% 
#   factor(levels=(c("Winter", "Spring", "Summer","Early_fall", "Fall")))
# 
# reg_all_slopes_chosen_silva_tax %>%
#   colnames()
# 
# reg_all_slopes_chosen_silva_tax %>%
#   colnames()
# 
# reg_all_slopes_chosen_silva_tax %>%
#   filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
#   filter(slope_chosen > 0) %>%
#   #subset(class == "Gammaproteobacteria") %>%
#   #filter(pseudoabundance > 1000) %>%
#   ggplot(aes(family, slope_chosen))+ #, color = asv_num
#   geom_boxplot()+
#   geom_point(alpha= 0.6)+
#   labs(x= "Phylogeny", y = "Slope")+
#   scale_x_discrete()+
#   #geom_smooth(method = "lm")+
#   scale_color_manual(values = palf_large(236))+
#   #facet_wrap(.~season, scales = "fixed")+
#   #guides(color=guide_legend(ncol = 1))+
#   theme_bw()+ 
#   theme(strip.text.x = element_text(size = 5), axis.text.x = element_text(angle = 90, size = 5), legend.position = "none")
# 
# SAR11 <- reg_all_t0234_tax %>%
#   filter(slope.pval < 0.05 ) %>% #& intercept.pval < value_intercept_pval
#   filter(slope > 0) %>%
#   subset(order == "SAR11 clade")


##FIGURE 2 (GO TO LINE 4145)-----------------------------
## GROWTH RATES DISTRIBUTION PATTERNS AT DIFFERENT TAXONOMIC LEVELS
##Single ASV growth rates distribution pattern for different phylums
#functions needed
#source("src/growth.rates.distr.tax.ranks.R") ##for all phylums
#source("src/growth.rates.distr.tax.ranks.divided.class.R") ##just for proteobacteria because we separate alpha and gamma
#source("src/growth.rates.distr.tax.ranks.divided.freqpoly.R") #with freqpoly not geom_density
#source("src/growth.rates.distr.tax.ranks.freqpoly.R") 
source('src/growth.rates.distr.tax.ranks.ridges.divided.R')
source('src/growth.rates.distr.tax.ranks.ridges.R')

#packages
library(lubridate)
library(scales)

##proves------
reg_all_slopes_chosen_silva_tax_1perc %>%
  colnames()
reg_all_slopes_chosen_silva_tax_1perc %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  #filter(pseudoabundance > 1000) %>%
  ggplot(aes(x=factor(slope_chosen_days)))+ #
  stat_function(fun = dnorm)+
  #geom_point(aes(color = season), alpha = 0.7, position = position_jitter(0.1))+
  #geom_boxplot(alpha= 0.5)+
  geom_bar(stat = "count")+
  facet_zoom(class == "Actinobaccteria") %>%
  #geom_violin(alpha = 0.0, draw_quantiles = c(0.25, 0.75), )+
  labs(color = "Season")+
  # stat_summary(fun = "mean",
  #              geom = "crossbar", 
  #              width = 0.2,
  #              colour = "black")+
  labs(x= expression("Growth rate (μ) day"^"-1"), y = "GR counts")+
  #scale_x_discrete()+
  #scale_y_continuous(expand = c(0, 0), limits = c(0,10))+ #forçar que comenci al 0 l'eix.
  #geom_smooth(method = "lm")
  scale_color_manual(values = palette_seasons)+
  #facet_wrap(.~phylum, scales = "free")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0), legend.position = "right",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#bi-modal distrubution?
## hi ha algo que no está  bé perquè no pot arribar a 3000 els counts.
reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  mutate(slope_chosen_days_factor = as.factor(slope_chosen_days)) %>%
  #filter(pseudoabundance > 1000) %>%
  group_by(phylum, slope_chosen_days_factor, season, treatment) %>%
  summarise(counts = n()) %>%
  ungroup() %>%
  ggplot(aes(x=as.numeric(slope_chosen_days_factor), y = counts))+ #
  #stat_function(fun = dnorm)+
  #geom_point(aes(color = season), alpha = 0.7, position = position_jitter(0.1))+
  #geom_boxplot(alpha= 0.5)+
  geom_col()+
  #geom_bar()+#stat = "count"
  #facet_zoom(class == "Actinobaccteria") %>%
  #geom_violin(alpha = 0.0, draw_quantiles = c(0.25, 0.75), )+
  labs(color = "Season")+
  # stat_summary(fun = "mean",
  #              geom = "crossbar", 
  #              width = 0.2,
  #              colour = "black")+
  labs(x= expression("Growth rate (μ) day"^"-1"), y = "GR counts")+
  #scale_x_discrete()+
  #scale_y_continuous(expand = c(0, 0), limits = c(0,10))+ #forçar que comenci al 0 l'eix.
  #geom_smooth(method = "lm")
  scale_color_manual(values = palette_seasons)+
  #facet_wrap(.~phylum, scales = "free")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0), legend.position = "right",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

##testing types of graphs-----

reg_all_slopes_chosen_silva_tax %>%
  glimpse()

reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %$%
  slope_chosen_days %>%
  unique()

# test <- reg_all_slopes_chosen_silva_tax %>%
#   filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
#   filter(slope_chosen_days > 0) %>%
#   mutate(slope_chosen_days_n = as.numeric(slope_chosen_days)) 
# test %>%
#   cut_number(test$slope_chosen_days_n, n = NULL, width = 1) %>%
#   mutate(slope_chosen_days_factor = as.factor(slope_chosen_days),
#          interval_for_count = interval(slope_chosen_days, 0.5)) %>%
#   #filter(pseudoabundance > 1000) %>%
#   group_by(phylum_f, slope_chosen_days_factor) %>%
#   summarise(counts = n()) %>%
#   ggplot(aes(x=slope_chosen_days_factor, y = counts))+ #
#   #stat_function(fun = dnorm)+
#   #geom_point(aes(color = season), alpha = 0.7, position = position_jitter(0.1))+
#   #geom_boxplot(alpha= 0.5)+
#   geom_density(adjust = 50)+
#   #geom_col()+
#   #geom_bar(stat = "count")+#
#   #geom_histogram(binwidth = 0.5)+
#   #facet_zoom(class == "Actinobaccteria") %>%
#   #geom_violin(alpha = 0.0, draw_quantiles = c(0.25, 0.75), )+
#   #labs(color = "Season")+
#   # stat_summary(fun = "mean",
#   #              geom = "crossbar", 
#   #              width = 0.2,
#   #              colour = "black")+
#   labs(x= expression("Growth rate (μ) day"^"-1"), y = "GR counts")+
#   #scale_x_discrete()+
#   #scale_y_continuous(expand = c(0, 0), limits = c(0,10))+ #forçar que comenci al 0 l'eix.
#   #geom_smooth(method = "lm")
#   scale_color_manual(values = palette_seasons)+
#   #facet_wrap(.~family_f, scales = "fixed")+
#   #guides(color=guide_legend(ncol = 1))+
#   theme_bw()+ 
#   theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0), legend.position = "right",
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank())

##Density plot----
reg_all_slopes_chosen_silva_tax_1perc %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  #subset(asv_num == 'asv54') %>%
  #mutate(slope_chosen_days_factor = as.factor(slope_chosen_days)) %>%
  #filter(pseudoabundance > 1000) %>%
  #group_by(order, slope_chosen_days_factor) %>%
  #summarise(counts = n()) %>%
  ggplot(aes(x=slope_chosen_days))+ #, colour = phylum
  #geom_density(adjust = 0.25)
  geom_freqpoly(stat = 'count')+
  stat_bin(binwidth = 0.25)+
  #stat_function(fun = dnorm)+
  #geom_freqpoly()
 #geom_histogram()
  #geom_point(aes(color = season), alpha = 0.7, position = position_jitter(0.1))+
  #geom_boxplot(alpha= 0.5)+
  #geom_col()+
  #geom_bar()+#stat = "count"
  #facet_zoom(class == "Actinobaccteria") %>%
  #geom_violin(alpha = 0.0, draw_quantiles = c(0.25, 0.75), )+
  #labs(color = "Season")+
  # stat_summary(fun = "mean",
  #              geom = "crossbar", 
  #              width = 0.2,
  #              colour = "black")+
  labs(x= expression("Growth rate (μ) day"^"-1"), y = "Density")+
  #scale_x_discrete()+
  #scale_y_continuous(expand = c(0, 0), limits = c(0,10))+ #forçar que comenci al 0 l'eix.
  #geom_smooth(method = "lm")
  #scale_color_manual(values = palette_large)+
  #facet_grid(.~treatment, scales = "free")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0), legend.position = "right",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

##different taxonomic ranks for each phylum
reg_all_slopes_chosen_silva_tax_1perc %$%
  phylum %>%
  unique()

reg_all_slopes_chosen_silva_tax_1perc %>%
  subset(phylum == "Proteobacteria") %$%
  order %>%
  unique()

#"Cyanobacteria"  "Proteobacteria" "Bacteroidota" 
reg_all_slopes_chosen_silva_tax_1perc %>%
  colnames()

reg_all_slopes_chosen_silva_tax_1perc %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  subset(phylum == "Proteobacteria") %>%
  #mutate(slope_chosen_days_factor = as.factor(slope_chosen_days)) %>%
  #filter(pseudoabundance > 1000) %>%
  #group_by(order, slope_chosen_days_factor) %>%
  #summarise(counts = n()) %>%
  ggplot(aes(x=slope_chosen_days))+ #, colour = class
  geom_histogram(aes(x=slope_chosen_days, group = phylum, colour = "#B8102A"))+
  # geom_histogram(aes(x=slope_chosen_days, group = order, colour = "#5F9660"))+
  # geom_histogram(aes(x=slope_chosen_days, group = family, colour = "#FFDF26"))+
  # geom_histogram(aes(x=slope_chosen_days, group = genus, colour  = "#6073B6"))+
  geom_line(stat = "density")+ #, group = family, colour = "black"
  labs(x= expression("Growth rate (μ) day"^"-1"), y = "Density")+
  guides(colour = guide_legend(override.aes=list("Class", "Order", "Family", "Species")))+
  #stat_function(fun = dnorm)+
  #geom_point(aes(color = season), alpha = 0.7, position = position_jitter(0.1))+
  #geom_boxplot(alpha= 0.5)+
  #geom_col()+
  #geom_bar()+#stat = "count"
  #facet_zoom(class == "Actinobaccteria") %>%
  #scale_x_discrete()+
  #scale_y_continuous(expand = c(0, 0), limits = c(0,10))+ #forçar que comenci al 0 l'eix.
  #geom_smooth(method = "lm")
  #scale_color_manual(values = palette_large)+
  #facet_grid(.~season, scales = "free")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0), legend.position = "right",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

##Growth rates distribution plots by different taxonomic ranks------
##WHOLE DATASET
##geom_histograms------
##PROTEOBACTERIA
proteobac_class <-  reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  subset(phylum == "Proteobacteria") %>%
  ggplot(aes(x=slope_chosen_days))+ #, colour = class
  geom_histogram(aes(x=slope_chosen_days, group = class), colour = "#5F9660",  fill = "#5F9660", #adjust= 0.5, 
               show.legend = FALSE, alpha = 0.2)+
  labs(x= "", y = "", title ="Proteobacteria")+
  scale_x_continuous(limits = c(0,11))+
  theme_bw()+
  theme(strip.text.x = element_text(size = 0), axis.text.x = element_text(angle = 0, size = 0), 
        axis.ticks = element_blank(), legend.position = "right", plot.title = element_text(size = 10),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin= unit(c(0.2, 1, 0.2, 1), "lines")) #top, right, bottom, left

proteobac_order <- reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  subset(phylum == "Proteobacteria") %>%
  ggplot(aes(x=slope_chosen_days))+ #, colour = class
  geom_histogram(aes(x=slope_chosen_days, group = order), colour = "#B8102A", fill = "#B8102A",
               #adjust= 0.5, 
               show.legend = FALSE, alpha = 0.2)+
  scale_x_continuous(limits = c(0,11))+
  labs(x= "", y = "")+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 0), axis.text.x = element_text(angle = 0, size = 0), 
        axis.ticks = element_blank(), legend.position = "right",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin= unit(c(0.2, 1, 0.2, 1), "lines")) #top, right, bottom, left

proteobac_fam <- reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  subset(phylum == "Proteobacteria") %>%
  ggplot(aes(x=slope_chosen_days))+ #, colour = class
  geom_histogram(aes(x=slope_chosen_days, group = family), colour = "#FFDF26", fill = "#FFDF26", #adjust= 0.5, 
               show.legend = FALSE, alpha = 0.2)+
  scale_x_continuous(limits = c(0,11))+
  labs(x= "", y = "Density")+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 0), axis.text.x = element_text(angle = 0, size = 0), 
        axis.ticks = element_blank(), legend.position = "right",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin= unit(c(0.2, 1, 0.2, 1), "lines")) #top, right, bottom, left

proteobac_genus <-reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  subset(phylum == "Proteobacteria") %>%
  ggplot(aes(x=slope_chosen_days))+ #, colour = class
  geom_histogram(aes(x=slope_chosen_days, group = genus, colour  = "#6073B6", fill = "#6073B6"), #adjust= 0.5, 
               show.legend = FALSE, alpha = 0.2)+
  scale_x_continuous(limits = c(0,11))+
  labs(x= "", y = "")+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 0), axis.text.x = element_text(angle = 0, size = 0), 
        axis.ticks = element_blank(), legend.position = "right",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin= unit(c(0.2, 1, 0.2, 1), "lines")) #top, right, bottom, left

proteobac_sp <-reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  subset(phylum == "Proteobacteria") %>%
  ggplot(aes(x=slope_chosen_days))+ #, colour = class
  geom_histogram(aes(x=slope_chosen_days, group = species), colour  = "#6073B6", fill = "#6073B6", #adjust= 0.5, 
               show.legend = FALSE, alpha = 0.2)+
  labs(x= "", y = "")+
  scale_x_continuous(limits = c(0,11))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 0), axis.text.x = element_text(angle = 0, size = 0), 
        axis.ticks = element_blank(), legend.position = "right",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin= unit(c(0.2, 1, 0.2, 1), "lines")) #top, right, bottom, left

proteobac_asv <-reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  subset(phylum == "Proteobacteria") %>%
  ggplot(aes(x=slope_chosen_days))+ #, colour = class
  geom_histogram(aes(x=slope_chosen_days, group = asv_num), colour  = "#4F507F", 
               fill = "#4F507F", adjust= 0.5, show.legend = FALSE, alpha = 0.2)+
  labs(x= expression("Growth rate (μ) day"^"-1"), y = "")+
  scale_x_continuous(limits = c(0,11))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0), legend.position = "right",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  
        plot.margin= unit(c(0.2, 1, 0.2, 1), "lines"))


# grid.arrange(proteobac_class, proteobac_order, proteobac_fam, proteobac_genus, 
#              proteobac_sp, proteobac_asv, ncol = 1)

##Cyanobacteria
reg_all_slopes_chosen_silva_tax %>%
  str()

cyano_class <-  reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  subset(phylum == "Cyanobacteria") %>%
  ggplot(aes(x=slope_chosen_days))+ #, colour = class
  geom_freqpoly(aes(x=slope_chosen_days, group = class), colour = "#5F9660",  fill = "#5F9660", #adjust= 0.5, 
               show.legend = FALSE, alpha = 0.2)+
  labs(x= "", y = "", title ="Phylum Cyanobacteria")+
 # scale_x_continuous(limits = c(0,11))+
  scale_x_binned()+
  theme_bw()+
  theme(strip.text.x = element_text(size = 0), axis.text.x = element_text(angle = 0, size = 0), 
        axis.ticks = element_blank(), legend.position = "right",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(size = 10),
        plot.margin= unit(c(0.2, 1, 0.2, 1), "lines")) #top, right, bottom, left

cyano_order <- reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  subset(phylum == "Cyanobacteria") %>%
  ggplot(aes(x=slope_chosen_days))+ #, colour = class
  geom_histogram(aes(x=slope_chosen_days, group = order), colour = "#B8102A", fill = "#B8102A",
               #adjust= 0.5, 
               show.legend = FALSE, alpha = 0.2)+
  labs(x= "", y = "")+
  scale_x_continuous(limits = c(0,11))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 0), axis.text.x = element_text(angle = 0, size = 0), 
        axis.ticks = element_blank(), legend.position = "right",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin= unit(c(0.2, 1, 0.2, 1), "lines")) #top, right, bottom, left

cyano_fam <- reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  subset(phylum == "Cyanobacteria") %>%
  ggplot(aes(x=slope_chosen_days))+ #, colour = class
  geom_histogram(aes(x=slope_chosen_days, group = family), colour = "#FFDF26", fill = "#FFDF26", #adjust= 0.5, 
               show.legend = FALSE, alpha = 0.2)+
  labs(x= "", y = "Density")+
  scale_x_continuous(limits = c(0,11))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 0), axis.text.x = element_text(angle = 0, size = 0), 
        axis.ticks = element_blank(), legend.position = "right",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin= unit(c(0.2, 1, 0.2, 1), "lines")) #top, right, bottom, left

cyano_genus <-reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  subset(phylum == "Cyanobacteria") %>%
  ggplot(aes(x=slope_chosen_days))+ #, colour = class
  geom_histogram(aes(x=slope_chosen_days, group = genus, colour  = "#6073B6", fill = "#6073B6"), #adjust= 0.5, 
               show.legend = FALSE, alpha = 0.2)+
  labs(x= "", y = "")+
  scale_x_continuous(limits = c(0,11))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 0), axis.text.x = element_text(angle = 0, size = 0), 
        axis.ticks = element_blank(), legend.position = "right",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin= unit(c(0.2, 1, 0.2, 1), "lines")) #top, right, bottom, left

cyano_sp <-reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  subset(phylum == "Cyanobacteria") %>%
  ggplot(aes(x=slope_chosen_days))+ #, colour = class
  geom_histogram(aes(x=slope_chosen_days, group = species), colour  = "#6073B6", fill = "#6073B6", adjust= 0.5, 
               show.legend = FALSE, alpha = 0.2)+
  labs(x= "", y = "")+
  scale_x_continuous(limits = c(0,11))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 0), axis.text.x = element_text(angle = 0, size = 0), 
        #axis.ticks = element_blank(), legend.position = "right",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin= unit(c(0.2, 1, 0.2, 1), "lines")) #top, right, bottom, left

reg_all_slopes_chosen_silva_tax %>%
  colnames()
#cyano_asv <-
reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  subset(phylum == "Cyanobacteria") %>%
  ungroup() %>%
  ggplot(aes(x=slope_chosen_days))+ #, colour = class
 # geom_histogram(aes(x=slope_chosen_days, group = asv_num), colour  = "#4F507F", fill = "#4F507F", adjust= 0.5, 
 #               show.legend = FALSE, alpha = 0.2)+
  geom_density(adjust = 50, aes(group = asv_num))+
  labs(x= expression("Growth rate (μ) day"^"-1"), y = "")+
  scale_x_continuous(limits = c(0,11))+
  theme_bw()+ 
  theme(legend.position = "bottom", strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  
        plot.margin= unit(c(0.2, 1, 0.2, 1), "lines"))

cyano<- reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  subset(phylum = 'Cyanobacteria') %>%
  group_by(asv_num) %>%
  distinct(asv_num)

test <- reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  subset(phylum == 'Cyanobacteria') %>%
  group_by(asv_num) %>%
  add_count(asv_num, name = 'number of gr')  %>%
  subset(asv_num == 'asv54')

##surt raro el gràfic d'asv_num ##problem based on asv54
reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  subset(phylum == "Cyanobacteria") %>%
  ggplot(aes(x=slope_chosen_days))+ #, colour = class
  # geom_histogram(aes(x=slope_chosen_days, group = asv_num), #adjust= 0.5, 
  #                show.legend = TRUE, alpha = 0.2, binwidth = 0.1)+
  #geom_density(position = 'fill')+
  geom_freqpoly()+
  #labs(x= "", y = "")+
  #scale_x_continuous(limits = c(0,11))+
  facet_wrap(~asv_num, scales = 'free')+
  #scale_y_continuous(limits = c(0,50))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 0), axis.text.x = element_text(angle = 0, size = 0), 
        axis.ticks = element_blank(), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin= unit(c(0.2, 1, 0.2, 1), "lines"))

reg_all_slopes_chosen_silva_tax %>%
  glimpse()

#grid.arrange(cyano_class, cyano_order, cyano_fam, cyano_genus, cyano_sp, cyano_asv, ncol = 1)

##Bacteroidota
bac_class <-  reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  subset(phylum == "Bacteroidota") %>%
  ggplot(aes(x=slope_chosen_days))+ #, colour = class
  geom_histogram(aes(x=slope_chosen_days, group = class), colour = "#5F9660",  fill = "#5F9660", #adjust= 0.5, 
               show.legend = FALSE, alpha = 0.2)+
  labs(x= "", y = "", title = "Phylum Bacteroidota")+
  scale_x_continuous(limits = c(0,11))+
  theme_bw()+
  theme(strip.text.x = element_text(size = 0), axis.text.x = element_text(angle = 0, size = 0), 
        axis.ticks = element_blank(), legend.position = "right", plot.title = element_text(size = 10),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin= unit(c(0.2, 1, 0.2, 1), "lines")) #top, right, bottom, left

bac_order <- reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  subset(phylum == "Bacteroidota") %>%
  ggplot(aes(x=slope_chosen_days))+ #, colour = class
  geom_histogram(aes(x=slope_chosen_days, group = order), colour = "#B8102A", fill = "#B8102A",
               #adjust= 0.5, 
               show.legend = FALSE, alpha = 0.2)+
  labs(x= "", y = "")+
  scale_x_continuous(limits = c(0,11))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 0), axis.text.x = element_text(angle = 0, size = 0), 
        axis.ticks = element_blank(), legend.position = "right",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin= unit(c(0.2, 1, 0.2, 1), "lines")) #top, right, bottom, left

bac_fam <- reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  subset(phylum == "Bacteroidota") %>%
  ggplot(aes(x=slope_chosen_days))+ #, colour = class
  geom_histogram(aes(x=slope_chosen_days, group = family), colour = "#FFDF26", fill = "#FFDF26", #adjust= 0.5, 
               show.legend = FALSE, alpha = 0.2)+
  labs(x= "", y = "Density")+
  scale_x_continuous(limits = c(0,11))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 0), axis.text.x = element_text(angle = 0, size = 0), 
        axis.ticks = element_blank(), legend.position = "right",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin= unit(c(0.2, 1, 0.2, 1), "lines")) #top, right, bottom, left

bac_genus <-reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  subset(phylum == "Bacteroidota") %>%
  ggplot(aes(x=slope_chosen_days))+ #, colour = class
  geom_histogram(aes(x=slope_chosen_days, group = genus, colour  = "#6073B6", fill = "#6073B6"), #adjust= 0.5, 
               show.legend = FALSE, alpha = 0.2)+
  labs(x= "", y = "")+
  scale_x_continuous(limits = c(0,11))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 0), axis.text.x = element_text(angle = 0, size = 0), 
        axis.ticks = element_blank(), legend.position = "right",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin= unit(c(0.2, 1, 0.2, 1), "lines")) #top, right, bottom, left

bac_sp <-reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  subset(phylum == "Bacteroidota") %>%
  ggplot(aes(x=slope_chosen_days))+ #, colour = class
  geom_histogram(aes(x=slope_chosen_days, group = species), colour  = "#6073B6", fill = "#6073B6", #adjust= 0.5, 
               show.legend = FALSE, alpha = 0.2)+
  labs(x= "", y = "")+
  scale_x_continuous(limits = c(0,11))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 0), axis.text.x = element_text(angle = 0, size = 0), 
        axis.ticks = element_blank(), legend.position = "right",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin= unit(c(0.2, 1, 0.2, 1), "lines")) #top, right, bottom, left

bac_asv <-reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  subset(phylum == "Bacteroidota") %>%
  ggplot(aes(x=slope_chosen_days))+ #, colour = class
  geom_histogram(aes(x=slope_chosen_days, group = asv_num), colour  = "#4F507F", fill = "#4F507F", 
               #adjust= 0.5, 
               show.legend = FALSE, alpha = 0.2)+
  scale_x_continuous(limits = c(0,11))+
  labs(x= expression("Growth rate (μ) day"^"-1"), y = "")+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0), legend.position = "right",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  
        plot.margin= unit(c(0.2, 1, 0.2, 1), "lines"))

 grid.arrange(bac_class, bac_order, bac_fam, bac_genus, bac_sp, bac_asv, ncol = 1)

  grid.arrange(
  proteobac_class, cyano_class, bac_class, 
  proteobac_order, cyano_order,  bac_order,
  proteobac_fam, cyano_fam,  bac_fam,
  proteobac_genus, cyano_genus, bac_genus, 
  #proteobac_sp, cyano_sp, bac_sp,
  proteobac_asv, cyano_asv, bac_asv, ncol = 3)

##ALPHAPROTEOBACTERIA
reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  subset(class == "Alphaproteobacteria") %>%
  ggplot(aes(x=slope_chosen_days))+ #, colour = class
  geom_density(aes(x=slope_chosen_days, group = order), colour  = "#4F507F", fill = "#4F507F", 
               #adjust= 0.5,
               show.legend = FALSE, alpha = 0.2)+
  scale_x_continuous(limits = c(0,11))+
  labs(x= expression("Growth rate (μ) day"^"-1"), y = "")+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0), legend.position = "right",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  
        plot.margin= unit(c(0.2, 1, 0.2, 1), "lines"))

alpha_o <- reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  subset(class == "Alphaproteobacteria") %>%
  ggplot(aes(x=slope_chosen_days))+ #, colour = class
  geom_density(aes(x=slope_chosen_days, group = order), colour  = "#4F507F", fill = "#4F507F", 
               #adjust= 0.5,
               show.legend = FALSE, alpha = 0.2)+
  scale_x_continuous(limits = c(0,11))+
  labs(x= expression("Growth rate (μ) day"^"-1"), y = "")+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0), legend.position = "right",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  
        plot.margin= unit(c(0.2, 1, 0.2, 1), "lines"))

alpha_f <- reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  subset(class == "Alphaproteobacteria") %>%
  ggplot(aes(x=slope_chosen_days))+ #, colour = class
  geom_density(aes(x=slope_chosen_days, group = family), colour  = "#4F507F", fill = "#4F507F", 
               #adjust= 0.5,
               show.legend = TRUE, alpha = 0.2)+
  scale_x_continuous(limits = c(0,11))+
  labs(x= expression("Growth rate (μ) day"^"-1"), y = "")+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0), legend.position = "right",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  
        plot.margin= unit(c(0.2, 1, 0.2, 1), "lines"))

alpha_g <- reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  subset(class == "Alphaproteobacteria") %>%
  ggplot(aes(x=slope_chosen_days))+ #, colour = class
  geom_density(aes(x=slope_chosen_days, group = genus), colour  = "#4F507F", fill = "#4F507F", 
               #adjust= 0.5,
               show.legend = FALSE, alpha = 0.2)+
  scale_x_continuous(limits = c(0,11))+
  labs(x= expression("Growth rate (μ) day"^"-1"), y = "")+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0), legend.position = "right",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  
        plot.margin= unit(c(0.2, 1, 0.2, 1), "lines"))

#alpha_asv <- 
  reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  subset(class == "Alphaproteobacteria") %>%
  #group_by(asv_num) %>%
  ggplot(aes(x=slope_chosen_days))+ #, colour = class
  geom_histogram(aes(x=slope_chosen_days, group = asv_num), colour  = "#4F507F", fill = "#4F507F",
               #adjust= 0.5,
               binwidth = 0.1,
               #bins = 90,
               show.legend = FALSE, alpha = 0.2)+
  scale_x_continuous(limits = c(0,11))+
  labs(x= expression("Growth rate (μ) day"^"-1"), y = "")+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0), legend.position = "right",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  
        plot.margin= unit(c(0.2, 1, 0.2, 1), "lines"))

grid.arrange(alpha_o, alpha_f, alpha_g, alpha_asv, ncol = 1)

##GR distribution using functions -------
##Distribution by phylums-----
reg_all_slopes_chosen_silva_tax %$%
  phylum %>%
  unique()

bacteroidota <- growthrates.distribution.tax.rank(data = reg_all_slopes_chosen_silva_tax, 
                                  phylum_to_explore = 'Bacteroidota',
                                  title = 'Phylum Bacteroidota',
                                  bin_width_c = 0.2,
                                  bin_width_o = 0.2,
                                  bin_width_f = 0.2,
                                  bin_width_a = 0.2,
                                  axis_y_title = 'Density') %>%
  as_ggplot()

ciano <- growthrates.distribution.tax.rank(data = reg_all_slopes_chosen_silva_tax, 
                                  phylum_to_explore = 'Cyanobacteria',
                                  title = 'Phylum Cyanobacteria',
                                  bin_width_c = 0.2,
                                  bin_width_o = 0.2,
                                  bin_width_f = 0.2,
                                  bin_width_a = 0.2,
                                  axis_y_title = 'Density') %>%
  as_ggplot()

proteo <- growthrates.distribution.divided.class(data = reg_all_slopes_chosen_silva_tax, 
                                  phylum_to_explore = 'Proteobacteria',
                                  title = 'Phylum Proteobacteria',
                                  c1 = 'Alphaproteobacteria',
                                  c2 = 'Gammaproteobacteria',
                                  bin_width_c = 0.2,
                                  bin_width_o = 0.2,
                                  bin_width_f = 0.2,
                                  bin_width_a = 0.2,
                                  y1 = 0,
                                  y2 = 10,
                                  axis_y_title = 'Density') %>%
  as_ggplot()

acido <- growthrates.distribution.tax.rank(data = reg_all_slopes_chosen_silva_tax, 
                                  phylum_to_explore = 'Acidobacteriota',
                                  title = 'Phylum Acidobacteriota',
                                  bin_width_c = 0.2,
                                  bin_width_o = 0.2,
                                  bin_width_f = 0.2,
                                  bin_width_a = 0.2,
                                  axis_y_title = 'Density') %>%
  as_ggplot()

verruco <- growthrates.distribution.tax.rank(data = reg_all_slopes_chosen_silva_tax, 
                                  phylum_to_explore = 'Verrucomicrobiota',
                                  title = 'Phylum Verrucomicrobiota',
                                  bin_width_c = 0.2,
                                  bin_width_o = 0.2,
                                  bin_width_f = 0.2,
                                  bin_width_a = 0.2,
                                  axis_y_title = 'Density') %>% 
  as_ggplot()
# 
# grid.arrange(bacteroidota, ciano,  proteo, acido, verruco, ncol=5)
#              #arrangeGrob(align = 'v', ncol = 5, rel_heights = c(1/5, 1/5, 2/5, 1/5, 1/5)))
#              #layout_matrix = list(cbind(c(1, 2, 3,3,3, 4, 5))))
#              #layout_matrix = cbind(c(1,1), c(2,2), c(3,3,3,3), c(4,4), c(5,5)))

distribution_gr_plots <- multi_panel_figure(columns = 6, rows = 1, width = 350, height = 150, 
                              column_spacing = 5, unit = 'mm',
                              panel_label_type = 'none')

distribution_gr_plots  %<>%
  fill_panel(bacteroidota, column = 1, row = 1) %<>%
  fill_panel(ciano, column = 2, row = 1) %<>%
  fill_panel(proteo, column = 3:4, row = 1) %<>%
  fill_panel(acido, column = 5, row = 1) %<>%
  fill_panel(verruco, column = 6, row = 1)

ggsave('distribution_gr_plots_geom_density_adjust1.5.pdf', distribution_gr_plots, 
       path = "~/Documentos/Doctorat/REMEI/results/figures/raw/",
       width = 375,
       height = 150,
       units = 'mm')
#> `geom_smooth()` using method = 'loess' and formula 'y ~ x'
#Crec que hauria d'afegir el % que representen de la comunitat i el nº de GR amb les que estem dibuixant.

#el mateix pero amb geom_freqpoly per tenir-ho en tant per 1 i poder comparar millor----
reg_all_slopes_chosen_silva_tax %$%
  phylum %>%
  unique()
library(ggpubr)
reg_all_slopes_chosen_silva_tax %>%
stopifnot(is.data.frame())

reg_all_slopes_chosen_silva_tax %>%
  class()

reg_all_slopes_chosen_silva_tax %>%
  colnames()


library(ggrepel)
reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  #group_by(asv_num) %>%
  ggplot(aes(x=slope_chosen_days,  color = class_f))+ #, colour = class
  # geom_histogram(aes(x=slope_chosen_days, group = asv_num), colour  = "#4F507F", fill = "#4F507F",
  #                #adjust= 0.5,
  #                binwidth = 0.1,
  #                #bins = 90,
  #                show.legend = FALSE, alpha = 0.2)+
  #stat_bin(aes(group = class_f), alpha = 0.2, binwidth = 0.25)+
  #geom_freqpoly(aes(group = class_f), stat = 'bin', binwidth = 0.25)+ 
  geom_density(aes(group = class), adjust = 1.5)+
  #stat_density(adjust = 3)+
  scale_x_continuous(limits = c(0,11))+
  #geom_text(aes(y= 45, label = class_f), check_overlap = TRUE, nudge_x = 1.5, nudge_y = 11)+
  labs(x= expression("Growth rate (μ) day"^"-1"), y = "")+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  
        plot.margin= unit(c(0.2, 1, 0.2, 1), "lines"))

reg_all_slopes_chosen_silva_tax %>%
  colnames()

bacteroidota <- growthrates.distribution.tax.rank.freqpoly(data = reg_all_slopes_chosen_silva_tax, 
                                                  phylum_to_explore = 'Bacteroidota',
                                                  title = 'Bacteroidota',
                                                  bin_width_c = 0.25,
                                                  bin_width_o = 0.25,
                                                  bin_width_f = 0.25,
                                                  bin_width_a = 0.25,
                                                  axis_y_title = 'Frequency',
                                                  y1_ph = 0, y2_ph = 400,
                                                  y1_cl = 0, y2_cl = 400,
                                                  y1_o = 0, y2_o = 400,
                                                  y1_f = 0, y2_f = 400, 
                                                  y1_a = 0, y2_a = 20) %>%
  as_ggplot()

ciano <- growthrates.distribution.tax.rank.freqpoly(data = reg_all_slopes_chosen_silva_tax, 
                                           phylum_to_explore = 'Cyanobacteria',
                                           title = 'Cyanobacteria',
                                           bin_width_c = 0.25,
                                           bin_width_o = 0.25,
                                           bin_width_f = 0.25,
                                           bin_width_a = 0.25,
                                           axis_y_title = 'Frequency',
                                           y1_ph = 0, y2_ph = 30,
                                           y1_cl = 0, y2_cl = 30,
                                           y1_o = 0, y2_o = 30,
                                           y1_f = 0, y2_f = 30, 
                                           y1_a = 0, y2_a = 10) %>%
  as_ggplot()

proteo <- growthrates.distribution.divided.class.freqpoly(data = reg_all_slopes_chosen_silva_tax, 
                                                 phylum_to_explore = 'Proteobacteria',
                                                 title = 'Proteobacteria',
                                                 c1 = 'Alphaproteobacteria',
                                                 c2 = 'Gammaproteobacteria',
                                                 bin_width_c = 0.2,
                                                 bin_width_o = 0.2,
                                                 bin_width_f = 0.2,
                                                 bin_width_a = 0.2,
                                                 y1 = 0,
                                                 y2 = 50,
                                                 y1_a = 0, y2_a = 10,
                                                 axis_y_title = 'Frequency') %>%
  as_ggplot()

acido <- growthrates.distribution.tax.rank.freqpoly(data = reg_all_slopes_chosen_silva_tax, 
                                           phylum_to_explore = 'Acidobacteriota',
                                           title = 'Acidobacteriota',
                                           bin_width_c = 0.2,
                                           bin_width_o = 0.2,
                                           bin_width_f = 0.2,
                                           bin_width_a = 0.2,
                                           axis_y_title = 'Frequency',
                                           y1_ph = 0, y2_ph = 10,
                                           y1_cl = 0, y2_cl = 10,
                                           y1_o = 0, y2_o = 10,
                                           y1_f = 0, y2_f = 10, 
                                           y1_a = 0, y2_a = 5) %>%
  as_ggplot()

verruco <- growthrates.distribution.tax.rank.freqpoly(data = reg_all_slopes_chosen_silva_tax, 
                                             phylum_to_explore = 'Verrucomicrobiota',
                                             title = 'Verrucomicrobiota',
                                             bin_width_c = 0.2,
                                             bin_width_o = 0.2,
                                             bin_width_f = 0.2,
                                             bin_width_a = 0.2,
                                             axis_y_title = 'Frequency',
                                             y1_ph = 0, y2_ph = 15,
                                             y1_cl = 0, y2_cl = 15,
                                             y1_o = 0, y2_o = 15,
                                             y1_f = 0, y2_f = 15, 
                                             y1_a = 0, y2_a = 7.5) %>% 
  as_ggplot()

plancto <- growthrates.distribution.tax.rank.freqpoly(data = reg_all_slopes_chosen_silva_tax, 
                                                      phylum_to_explore = 'Planctomycetota',
                                                      title = 'Planctomycetota',
                                                      bin_width_c = 0.2,
                                                      bin_width_o = 0.2,
                                                      bin_width_f = 0.2,
                                                      bin_width_a = 0.2,
                                                      axis_y_title = 'Frequency', 
                                                      y1_ph = 0, y2_ph = 30,
                                                      y1_cl = 0, y2_cl = 30,
                                                      y1_o = 0, y2_o = 30,
                                                      y1_f = 0, y2_f = 30, 
                                                      y1_a = 0, y2_a = 5) %>% 
  as_ggplot()

crena <- growthrates.distribution.tax.rank.freqpoly(data = reg_all_slopes_chosen_silva_tax, 
                                                      phylum_to_explore = 'Crenarchaeota',
                                                      title = 'Crenarchaeota',
                                                      bin_width_c = 0.2,
                                                      bin_width_o = 0.2,
                                                      bin_width_f = 0.2,
                                                      bin_width_a = 0.2,
                                                      axis_y_title = 'Frequency',
                                                    y1_ph = 0, y2_ph = 15,
                                                    y1_cl = 0, y2_cl = 15,
                                                    y1_o = 0, y2_o = 15,
                                                    y1_f = 0, y2_f = 15, 
                                                    y1_a = 0, y2_a = 3) %>% 
  as_ggplot()

actino <- growthrates.distribution.tax.rank.freqpoly(data = reg_all_slopes_chosen_silva_tax, 
                                                    phylum_to_explore = 'Actinobacteriota',
                                                    title = 'Actinobacteriota',
                                                    bin_width_c = 0.2,
                                                    bin_width_o = 0.2,
                                                    bin_width_f = 0.2,
                                                    bin_width_a = 0.2,
                                                    axis_y_title = 'Frequency',
                                                    y1_ph = 0, y2_ph = 25,
                                                    y1_cl = 0, y2_cl = 25,
                                                    y1_o = 0, y2_o = 25,
                                                    y1_f = 0, y2_f = 25, 
                                                    y1_a = 0, y2_a = 5) %>% 
  as_ggplot()

firmi <- growthrates.distribution.tax.rank.freqpoly(data = reg_all_slopes_chosen_silva_tax, 
                                                     phylum_to_explore = 'Firmicutes',
                                                     title = 'Firmicutes',
                                                     bin_width_c = 0.2,
                                                     bin_width_o = 0.2,
                                                     bin_width_f = 0.2,
                                                     bin_width_a = 0.2,
                                                     axis_y_title = 'Frequency',
                                                    y1_ph = 0, y2_ph = 10,
                                                    y1_cl = 0, y2_cl = 10,
                                                    y1_o = 0, y2_o = 10,
                                                    y1_f = 0, y2_f = 10, 
                                                    y1_a = 0, y2_a = 2.5) %>% 
  as_ggplot()

bdell <- growthrates.distribution.tax.rank.freqpoly(data = reg_all_slopes_chosen_silva_tax, 
                                                    phylum_to_explore = 'Bdellovibrionota',
                                                    title = 'Bdellovibrionota',
                                                    bin_width_c = 0.2,
                                                    bin_width_o = 0.2,
                                                    bin_width_f = 0.2,
                                                    bin_width_a = 0.2,
                                                    axis_y_title = 'Frequency',
                                                    y1_ph = 0, y2_ph = 5,
                                                    y1_cl = 0, y2_cl = 5,
                                                    y1_o = 0, y2_o = 5,
                                                    y1_f = 0, y2_f = 5, 
                                                    y1_a = 0, y2_a = 2.5) %>% 
  as_ggplot()

##for the 5 phylums most relevant in boxplot (with asv at 1% at some point of the experiment)----
distribution_gr_plots_freqplot <- multi_panel_figure(columns = 6, rows = 1, width = 350, height = 150, 
                                            column_spacing = 5, unit = 'mm',
                                            panel_label_type = 'none')

distribution_gr_plots_freqplot  %<>%
  fill_panel(bacteroidota, column = 1, row = 1) %<>%
  fill_panel(ciano, column = 2, row = 1) %<>%
  fill_panel(proteo, column = 3:4, row = 1) %<>%
  fill_panel(acido, column = 5, row = 1) %<>%
  fill_panel(verruco, column = 6, row = 1)

ggsave('distribution_gr_plots_geom_freq0.5.pdf', distribution_gr_plots_freqplot, 
       path = "~/Documentos/Doctorat/REMEI/results/figures/raw/",
       width = 375,
       height = 150,
       units = 'mm')

##all phylums that represent >1% of GR calculated form the whole datset----
distribution_gr_plots_freqplot1perc <- multi_panel_figure(columns = 6, rows = 2, width = 210, height = 297, 
                                                     column_spacing = 1, unit = 'mm',
                                                     panel_label_type = 'none')

distribution_gr_plots_freqplot1perc  %<>%
  fill_panel(proteo, column = 1:2, row = 1) %<>%
  fill_panel(bacteroidota, column = 3, row = 1) %<>%
  fill_panel(plancto, column = 4, row = 1) %<>%
  fill_panel(crena, column = 5, row = 1) %<>%
  fill_panel(actino, column = 6, row = 1) %<>%
  fill_panel(ciano, column = 1, row = 2) %<>%
  fill_panel(verruco, column = 2, row = 2) %<>%
  fill_panel(firmi, column = 3, row = 2) %<>%
  fill_panel(acido, column = 4, row = 2) %<>%
  fill_panel(bdell, column = 5, row = 2) 

ggsave('distribution_gr_plots_geom_freq0.5_1perc.pdf', distribution_gr_plots_freqplot1perc, 
       path = "~/Documentos/Doctorat/REMEI/results/figures/raw/",
       width = 210,
       height = 297,
       units = 'mm')


##GEOM RIDGES FIGURE 2------
##fct_infreq: reorder by number of observations with each level
# reg_all_slopes_chosen_silva_tax %>%
# mutate(rank = factor(rank, levels = c('class', 'order', 'family', 'genus')))

#perquè tots els gràfics estiguin ordenats iguals
reg_all_slopes_chosen_silva_tax <- reg_all_slopes_chosen_silva_tax  %>%
  ungroup() %>%
  mutate(phylum_f = fct_rev(fct_infreq(phylum_f)))

reg_all_slopes_chosen_silva_tax %$% phylum_f %>% levels
palette_phylums_assigned

reg_all_slopes_chosen_silva_tax %>%
  glimpse()

counts <-   reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0 ) %>%
  group_by(phylum_f) %>%
  filter(n() >= 2) %>%
  mutate(counts = n()) %>%
  as_tibble() %>%
  group_by(phylum_f, counts) %>%
  summarize() %>%
  as_tibble()

ridge_ph <- 
  reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0 ) %>%
  group_by(phylum_f) %>%
  filter(n() >= 2) %>%
  mutate(counts = paste('n = ', n())) %>% #
  ggplot(aes(y = fct_rev(fct_infreq(phylum_f)), x = slope_chosen_days, fill = phylum_f,  group = phylum_f, label = counts))+#fct_rev(fct_infreq(phylum_f))
  geom_density_ridges(alpha = 0.8, panel_scaling = TRUE, scale = 1,
                      jittered_points = TRUE,
                      point_shape = 21, point_size = 0.2, point_alpha = 0.0,
                      #position = 
                      quantile_lines = TRUE,
                      quantile_fun = mean)+
  scale_fill_manual(values = palette_phylums_assigned)+
  scale_x_continuous(limits = c(0,11))+
  coord_cartesian(clip = "off") +
  #scale_y_discrete(position = "right") +
  # stat_summary(fun = 'mean',
  #              geom = 'line', color = 'black')+
  #annotate(geom = 'text', label = counts)+
  geom_text(nudge_x = 8.8, nudge_y = 0.3, check_overlap = TRUE, size = 3)+
  #coord_flip()+
  labs(y = '', fill= 'Phylum', x = expression("Growth rate day"^"-1"), title = 'Phylum')+
  theme_bw()+
  theme(legend.position = "none", strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 12),
        plot.margin= unit(c(1, 1, 0.2, 0.2), "lines"), panel.border = element_blank())

# reg_all_slopes_chosen_silva_tax %>%
#   mutate(class_freq = fct_reorder(phylum_f, class_f))

ridge_cl <- reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05,
         slope_chosen_days > 0) %>% #& intercept.pval < value_intercept_pval
  group_by(phylum_f) %>%
  filter(n() >= 2) %>%
  #summarize(n = n())
  # filter(n() > 2) %>% #, .preserve = FALSE
  # ungroup() %>%
  group_by(phylum_f) %>%
  mutate(counts = paste('n = ', n())) %>% #
  ggplot(aes(y = fct_rev(fct_infreq(phylum_f)), x = slope_chosen_days, fill = phylum_f, label = counts))+
  geom_density_ridges(aes(fill = phylum_f, group = class_f), scale = 1, alpha = 0.7,
                      jittered_points = TRUE,
                      point_shape = 21, point_size = 0.2, point_alpha = 0.0,
                      #position = 
                      quantile_lines = TRUE,
                      quantile_fun = mean)+
  scale_fill_manual(values = palette_phylums_assigned)+
  ##facet_grid(vars(class))+
  labs(y = 'Class', fill= 'Phylum', x =expression("Growth rate day"^"-1"), title = 'Class')+
  scale_x_continuous(limits = c(0,11))+
  #facet_grid(rows = vars(class), scales = 'free')+
  geom_text(nudge_x = 7.8, nudge_y = 0.3, check_overlap = TRUE, size = 3)+
  #scale_y_discrete(position = "right") +
  #coord_flip()+
  theme_bw()+
  theme(legend.position = "none", strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 0),
        plot.margin= unit(c(1, 1, 0.2, 0.2), "lines"), axis.title.y = element_blank(), panel.border = element_blank())

ridge_o <- reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  group_by(phylum_f) %>%
  filter(n() >= 2) %>%
  ungroup() %>%
  group_by(phylum_f) %>%
  mutate(counts = paste('n = ', n())) %>% #
  ggplot(aes(y = fct_rev(fct_infreq(phylum_f)), x = slope_chosen_days, fill = phylum_f, label = counts))+#fct_infreq(phylum_f, order_f) per ordenar per freqüència
  geom_density_ridges(aes(fill = phylum_f, group = order_f), scale = 1.0, alpha = 0.6,
                      jittered_points = TRUE,
                      point_shape = 21, point_size = 0.2, point_alpha = 0.0,
                      #position = 
                      quantile_lines = TRUE,
                      quantile_fun = mean)+
  scale_fill_manual(values = palette_phylums_assigned)+
  labs(y = '', fill= 'Phylum', x =expression("Growth rate day"^"-1"), title = 'Order')+
  scale_x_continuous(limits = c(0,11))+
  #facet_grid(rows = vars(class), scales = 'free')+
  #scale_y_discrete(position = "right") +
  geom_text(nudge_x = 7.8, nudge_y = 0.3, check_overlap = TRUE, size = 3)+
  #coord_flip()+
  theme_bw()+
  theme(legend.position = "none", strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 0),
        plot.margin= unit(c(1, 1, 0.2, 0.2), "lines"), panel.border = element_blank())

ridge_f <- reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  group_by(phylum_f) %>%
  filter(n() >= 2) %>%
  # filter(n() > 2) %>%
  # ungroup() %>%
  group_by(phylum_f) %>%
  mutate(counts = paste('n = ', n())) %>% #
  ungroup() %>%
  ggplot(aes(y = fct_rev(fct_infreq(phylum_f)), x = slope_chosen_days, fill = phylum_f, label = counts))+
  geom_density_ridges(aes(fill = phylum_f, group = family_f), scale = 1.0, alpha = 0.6,
                      jittered_points = TRUE,
                      point_shape = 21, point_size = 0.2, point_alpha = 0.0,
                      #position = 
                      quantile_lines = TRUE,
                      quantile_fun = mean)+#fct_rev(fct_infreq(phylum_f))
  scale_fill_manual(values = palette_phylums_assigned)+
  labs(y = '', fill= 'Phylum', x =expression("Growth rate day"^"-1"), title = 'Family')+
  scale_x_continuous(limits = c(0,11))+
  geom_text(nudge_x = 7.8, nudge_y = 0.3, check_overlap = TRUE, size = 3)+
  #facet_grid(rows = vars(class), scales = 'free')+
  theme_bw()+
  theme(legend.position = "none", strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 0),
        plot.margin= unit(c(1, 1, 0.2, 0.2), "lines"), panel.border = element_blank())


##reordeno els phylums perquè al filtrar em canvia l'ordre si ho faig per freqüències
reg_all_slopes_chosen_silva_tax_asv_reorder <- reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  # group_by(asv_f) %>%
  # filter(n() > 2) %>%
  group_by(phylum_f) %>%
  mutate(counts = paste('n = ', n())) %>% #
  ungroup()

reg_all_slopes_chosen_silva_tax_asv_reorder$phylum_f <- reg_all_slopes_chosen_silva_tax_asv_reorder$phylum_f %>%
  factor( levels = c("Proteobacteria",  "Bacteroidota", "Planctomycetota", "Actinobacteriota", 
                    "Crenarchaeota", "Verrucomicrobiota", "Cyanobacteria", "Firmicutes", "Acidobacteriota",
                    "Nitrospinota", "Bdellovibrionota", "Nitrospirota",  "Myxococcota", "Desulfobacterota",
                    "Deinococcota" , "Chloroflexi" , "Fusobacteriota", 
                    "Spirochaetota", "Campilobacterota", "Abditibacteriota", "Latescibacterota",
                    "Methylomirabilota", "Halanaerobiaeota", "Sumerlaeota", "Calditrichota", "Gemmatimonadota"), 
          ordered = TRUE)
# reg_all_slopes_chosen_silva_tax_asv_reorder$phylum_f %>%
#   levels()

ridge_a <- reg_all_slopes_chosen_silva_tax_asv_reorder %>%
  group_by(phylum_f) %>%
  filter(n() >= 2) %>%
  ggplot(aes(y = fct_rev(phylum_f), x = slope_chosen_days, label = counts))+
  geom_density_ridges(aes(fill = phylum_f, group = asv_f), scale = 1.0, alpha = 0.6,
                      jittered_points = TRUE,
                      point_shape = 21, point_size = 0.2, point_alpha = 0.0,
                      #position =
                      quantile_lines = TRUE,
                      quantile_fun = mean)+
  scale_fill_manual(values = palette_phylums_assigned)+
  labs(y = '', fill= 'Phylum', x =expression("Growth rate day"^"-1"), title = 'ASV')+
  scale_x_continuous(limits = c(0,11))+
  geom_text(nudge_x = 7.2, nudge_y = 0.35, check_overlap = TRUE, size = 3)+
  #facet_grid(rows = vars(class), scales = 'free')+
  coord_cartesian(clip = "off") +
  theme_bw()+
  theme(legend.position = "none", strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 0),
        plot.margin= unit(c(1, 1, 0.2, 0.2), "lines"), panel.border = element_blank())


##Percentatge de GR significatives calculades per cada phylum ------ 
#library(forcats)
reg_all_slopes_chosen_silva_tax %$%
  phylum_f %>%
  levels()
perc_gr_phylum <- reg_all_slopes_chosen_silva_tax %>%
  subset(slope_chosen_days != is.na(slope_chosen_days)) %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  group_by(phylum_f) %>%
  filter(n() > 2) %>%
  group_by(phylum_f, slope_chosen_days) %>% 
  summarize(n = n()) %>%
  group_by(phylum_f) %>%
  summarize(n = n()) %>%
  mutate(sum = sum(n),
         perc = n/sum) %>%
  ungroup() %>%
  ggplot(aes(sum, perc, fill = phylum_f))+#group = phylum_f,
  geom_col(aes(fill = fct_reorder(phylum_f, perc, .desc = TRUE)))+ #ct_rev(fct_infreq(phylum_f, ordered = NA)
  #coord_polar(theta='y')+
  #scale_x_continuous(limits = c(14584, 14584))+
  scale_y_continuous(labels = percent_format())+
  scale_fill_manual(values = palette_phylums_assigned)+
  labs(x = 'Phylum', y = 'Percentage of significant\ngrowth rates calculated', fill = 'Phylum')+
  theme_bw()+
  theme(axis.text.x = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), aspect.ratio = 10/2, panel.border = element_blank(),
        axis.ticks = element_blank(), legend.position = 'none', axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14))


# reg_all_slopes_chosen_silva_tax %>%
#   subset(slope_chosen_days != is.na(slope_chosen_days)) %>%
#   group_by(phylum_f, slope_chosen_days) %>% 
#   summarize(n = n()) %>%
#   group_by(phylum_f) %>%
#   summarize(n = n()) %>%
#   mutate(sum = sum(n),
#          perc = n/sum) %>%
#   ggplot(aes(fct_infreq(phylum_f, ordered = NA), perc, color = phylum_f))+ #group = phylum_f,
#   geom_point(aes(color = phylum_f))+
#   scale_y_continuous(labels = percent_format())+
#   scale_color_manual(values = palette_large)+
#   labs(x = 'Phylum', y = 'Percentage of growth\nrates calculated', fill = 'Phylum')+
#   theme_bw()+
#   theme(axis.text.x = element_blank(), panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(), aspect.ratio = 10/3,
#         axis.ticks = element_blank())

##Percentatge de GR significatives representades en cada nivell taxonòmic (+ de 2 grups per poder fer distribució)----
gr_asv_perc <- reg_all_slopes_chosen_silva_tax %>%
  subset(slope_chosen_days != is.na(slope_chosen_days)) %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  group_by(asv_f) %>%
  filter(n() > 2) %>%
  summarize(n = n()) %>%
  mutate(sum = sum(n),
         asv_gr_perc = sum/3749) %>%
  select(asv_gr_perc) %>%
  unique() %>%
  as_tibble()

gr_f_perc <- reg_all_slopes_chosen_silva_tax %>%
  subset(slope_chosen_days != is.na(slope_chosen_days)) %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  group_by(family) %>%
  filter(n() > 2) %>%
  summarize(n = n()) %>%
  mutate(sum = sum(n),
         f_gr_perc = sum/3749) %>%
  select(f_gr_perc) %>%
  unique() %>%
  as_tibble()

gr_o_perc <- reg_all_slopes_chosen_silva_tax %>%
  subset(slope_chosen_days != is.na(slope_chosen_days)) %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  group_by(order) %>%
  filter(n() > 2) %>%
  summarize(n = n()) %>%
  mutate(sum = sum(n),
         order_gr_perc = sum/3749) %>%
  select(order_gr_perc) %>%
  unique() %>%
  as_tibble()

gr_c_perc <- reg_all_slopes_chosen_silva_tax %>%
  subset(slope_chosen_days != is.na(slope_chosen_days)) %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  group_by(class) %>%
  filter(n() > 2) %>%
  summarize(n = n()) %>%
  mutate(sum = sum(n),
         class_gr_perc = sum/3749) %>%
  select(class_gr_perc) %>%
  unique() %>%
  as_tibble()

gr_p_perc <- reg_all_slopes_chosen_silva_tax %>%
  subset(slope_chosen_days != is.na(slope_chosen_days)) %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  group_by(phylum) %>%
  filter(n() > 2) %>%
  summarize(n = n()) %>%
  mutate(sum = sum(n),
         phylum_gr_perc = sum/3749) %>%
  select(phylum_gr_perc) %>%
  unique() %>%
  as_tibble() 


labels_perc_gr_rank <-  as_labeller(c(phylum_gr_perc = 'Phylum',
                                      class_gr_perc = 'Class',
                                      order_gr_perc = 'Order',
                                      f_gr_perc = 'Family',
                                      asv_gr_perc = 'ASV'))

gr_tax_ranks_perc <- cbind(gr_asv_perc, gr_f_perc, gr_o_perc, gr_c_perc, gr_p_perc) %>%
  as_tibble() %>%
  pivot_longer(cols = 1:5) %>%
  ggplot(aes(x = as_factor(name), y = value, group = 1))+
  geom_point()+
  coord_flip()+
  geom_line()+
  labs(y='Percentage', x = 'Growth rates distribution plotted\nat different taxonomic ranks')+
  scale_x_discrete(labels = labels_perc_gr_rank)+
  scale_y_continuous(labels = percent_format())+
  theme_bw()+
  theme(aspect.ratio = 10/3, panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = 'none', axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14))


##Mean GR at different taxonomic levels
labels_mean_rank <- as_labeller(c( mean_phylum = 'Phylum',
                                   mean_class = 'Class',
                                   mean_order = 'Order',
                                   mean_family = 'Family',
                                   mean_asv = 'ASV'))

mean_gr_ranks <- reg_all_slopes_chosen_silva_tax %>%
  subset(slope_chosen_days != is.na(slope_chosen_days)) %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  group_by(phylum_f) %>%
  filter(n() > 2) %>%
  mutate(mean_phylum = mean(slope_chosen_days)) %>%
  ungroup() %>%
  group_by(class) %>%
  mutate(mean_class = mean(slope_chosen_days)) %>%
  ungroup() %>%
  group_by(order) %>%
  mutate(mean_order = mean(slope_chosen_days)) %>%
  ungroup() %>%
  group_by(family) %>%
  mutate(mean_family = mean(slope_chosen_days)) %>%
  ungroup() %>%
  group_by(asv_num) %>%
  mutate(mean_asv = mean(slope_chosen_days)) %>%
  select(starts_with('mean'), 'phylum_f', 'class', 'order', 'family', 'asv_num') %>%
  pivot_longer(cols = starts_with('mean')) %>%
  group_by(phylum_f) %>%
  distinct(value, name, phylum_f) %>%
  ggplot(aes(fct_rev(as_factor(name)), value))+
  geom_point(aes(color = phylum_f), alpha = 0.9, position = position_jitter(0.3))+
  #geom_boxplot(alpha= 0.05)+
  geom_violin(alpha = 0.2, draw_quantiles = c(0.25, 0.75))+
  scale_color_manual(values = palette_phylums_assigned)+
  labs(y = 'Mean growth rate at different\ntaxonomic ranks', x= 'Taxonomic rank', color = 'Phylum')+
  scale_x_discrete(labels = labels_mean_rank)+
  coord_flip()+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = 'none', axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14))


##figure 2 creating pannel design------
distribution_gr_plots_ridge <- multi_panel_figure(columns = 7, rows = 3, width = 410, height = 297, 
                                                          column_spacing = 0.1, unit = 'mm',
                                                          panel_label_type = 'none')

distribution_gr_plots_ridge  %<>%
  fill_panel(perc_gr_phylum, column = 1, row = 1) %<>%
  fill_panel(gr_tax_ranks_perc, column = 1, row = 2) %<>%
  fill_panel(mean_gr_ranks, column = 1, row = 3) %<>%
  fill_panel(ridge_ph, column = 2:3, row = 1:3) %<>%
  fill_panel(ridge_cl, column = 4, row = 1:3) %<>%
  fill_panel(ridge_o, column = 5, row = 1:3) %<>%
  fill_panel(ridge_f, column = 6, row = 1:3) %<>%
  fill_panel(ridge_a, column = 7, row = 1:3)

##CAMBIAR EL NOM ABANS DE GUARDAR!!!!
# ggsave('distribution_gr_plots_ridge_order_freq_perc6.pdf', distribution_gr_plots_ridge, 
#        path = "~/Documentos/Doctorat/REMEI/results/figures/raw/",
#        width = 470,
#        height = 397,
#        units = 'mm')
# palette_phylums

###FIGURE 2 SEPARATED BY TREATMENTS (SUPPLEMENTARY)------
##Geom riges for asv, order and family levels divided by treatments
counts_divided <-   reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0 ) %>%
  group_by(order_f) %>%
  filter(n() >= 10) %>%
  group_by(phylum_f, treatment) %>%
  mutate(counts = n()) %>%
  as_tibble() %>%
  group_by(phylum_f, counts, treatment) %>%
  summarize() %>%
  as_tibble()

write.table(counts_divided, 'results/tables/counts_divided.txt', sep = '\t')

reg_all_slopes_chosen_silva_tax_asv_reorder$phylum_f <- reg_all_slopes_chosen_silva_tax_asv_reorder$phylum_f %>%
  factor( levels = c("Proteobacteria",  "Bacteroidota", "Planctomycetota", "Actinobacteriota", 
                     "Crenarchaeota", "Verrucomicrobiota", "Cyanobacteria", "Firmicutes", "Acidobacteriota",
                     "Nitrospinota", "Bdellovibrionota", "Nitrospirota",  "Myxococcota", "Desulfobacterota",
                     "Deinococcota" , "Chloroflexi" , "Fusobacteriota", 
                     "Spirochaetota", "Campilobacterota", "Abditibacteriota", "Latescibacterota",
                     "Methylomirabilota", "Halanaerobiaeota", "Sumerlaeota", "Calditrichota", "Gemmatimonadota"), 
          ordered = TRUE)
# reg_all_slopes_chosen_silva_tax_asv_reorder$phylum_f %>%
#   levels()

reg_all_slopes_chosen_silva_tax_asv_reorder$treatment <- reg_all_slopes_chosen_silva_tax_asv_reorder$treatment %>%
factor(levels=(c("CL", "CD", "PL", "PD", "DL", "VL")))

ridge_a_divided <- reg_all_slopes_chosen_silva_tax_asv_reorder %>%
  # group_by(asv_num) %>%
  # filter(n() >= 2) %>%
  ggplot(aes(y = fct_rev(phylum_f), x = slope_chosen_days, label = counts))+ #, 
  geom_density_ridges(aes(fill = phylum_f, group = asv_f), scale = 1.0, alpha = 0.6,
                      jittered_points = TRUE,
                      point_shape = 21, point_size = 0.2, point_alpha = 0.0,
                      #position =
                      quantile_lines = TRUE,
                      quantile_fun = mean)+
  scale_fill_manual(values = palette_phylums_assigned)+
  labs(y = '', fill= 'Phylum', x =expression("Growth rate day"^"-1"), title = 'ASV')+
  scale_x_continuous(limits = c(0,11))+
  facet_grid(cols = vars(treatment), scales = 'free')+
  geom_text(nudge_x = 7.2, nudge_y = 0.35, check_overlap = TRUE, size = 3)+
  #annotate('text', label = counts_divided, x=8, y=1)+
  coord_cartesian(clip = "off") +
  theme_bw()+
  theme(legend.position = "none", strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 0),
        plot.margin= unit(c(1, 1, 0.2, 0.2), "lines"), panel.border = element_blank())

ggsave('distribution_gr_plots_ridge_asv_treatment_divided.pdf', ridge_a_divided,
       path = "~/Documentos/Doctorat/REMEI/results/figures/raw/",
       width = 180,
       height = 150,
       units = 'mm')

ridge_family_divided <- reg_all_slopes_chosen_silva_tax_asv_reorder %>%
  group_by(family_f) %>%
  filter(n() >= 10) %>%
  ggplot(aes(y = fct_rev(phylum_f), x = slope_chosen_days, label = counts))+
  geom_density_ridges(aes(fill = phylum_f, group = family_f), scale = 1.0, alpha = 0.6,
                      jittered_points = TRUE,
                      point_shape = 21, point_size = 0.2, point_alpha = 0.0,
                      #position =
                      quantile_lines = TRUE,
                      quantile_fun = mean)+
  scale_fill_manual(values = palette_phylums_assigned)+
  labs(y = '', fill= 'Phylum', x =expression("Growth rate day"^"-1"), title = 'Family')+
  scale_x_continuous(limits = c(0,11))+
  geom_text(nudge_x = 7.2, nudge_y = 0.35, check_overlap = TRUE, size = 3)+
  facet_grid(cols = vars(treatment), scales = 'free')+
  coord_cartesian(clip = "off") +
  theme_bw()+
  theme(legend.position = "none", strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 0),
        plot.margin= unit(c(1, 1, 0.2, 0.2), "lines"), panel.border = element_blank())

ggsave('distribution_gr_plots_ridge_family_treatment_divided.pdf', ridge_family_divided,
       path = "~/Documentos/Doctorat/REMEI/results/figures/raw/",
       width = 180,
       height = 150,
       units = 'mm')

ridge_order_divided <- reg_all_slopes_chosen_silva_tax_asv_reorder %>%
  group_by(order_f) %>%
  filter(n() >= 10) %>%
  ggplot(aes(y = fct_rev(phylum_f), x = slope_chosen_days, label = counts))+
  geom_density_ridges(aes(fill = phylum_f, group = order_f), scale = 1.0, alpha = 0.6,
                      jittered_points = TRUE,
                      point_shape = 21, point_size = 0.2, point_alpha = 0.0,
                      #position =
                      quantile_lines = TRUE,
                      quantile_fun = mean)+
  scale_fill_manual(values = palette_phylums_assigned)+
  labs(y = '', fill= 'Phylum', x =expression("Growth rate day"^"-1"), title = 'Order')+
  scale_x_continuous(limits = c(0,11))+
  geom_text(nudge_x = 7.2, nudge_y = 0.35, check_overlap = TRUE, size = 3)+
  facet_grid(cols = vars(treatment), scales = 'free')+
  coord_cartesian(clip = "off") +
  theme_bw()+
  theme(legend.position = "none", strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 0),
        plot.margin= unit(c(1, 1, 0.2, 0.2), "lines"), panel.border = element_blank())

ggsave('distribution_gr_plots_ridge_order_treatment_divided.pdf', ridge_order_divided,
       path = "~/Documentos/Doctorat/REMEI/results/figures/raw/",
       width = 180,
       height = 150,
       units = 'mm')


##Geom_ridges function
## Proteobacteria    Crenarchaeota     Bacteroidota      Myxococcota       Planctomycetota   Verrucomicrobiota Nitrospirota     
## Actinobacteriota  Acidobacteriota   Nitrospinota      Bdellovibrionota  Deinococcota      Firmicutes        Cyanobacteria   
## Chloroflexi       Desulfobacterota 

##phylums gr distribution separated ordered by frequency
ridge_proteo <- 
  growthrates.distribution.tax.rank.ridges.divided(
  data = reg_all_slopes_chosen_silva_tax_asv_reorder, 
  phylum_to_explore = 'Proteobacteria',
  title = 'Proteobacteria class',
  axis_y_title = '',
  x_c = 10,
  x_o = 10,
  x_f = 10,
  x_a = 8) %>%
  as_ggplot() ##estic afegint fct_rev perquè m'ordeni els gràfics de manera lògica alpha i gamma

ggsave('ridge_proteobacteria_filtered_v4.pdf', ridge_proteo, 
       path = "~/Documentos/Doctorat/REMEI/results/figures/raw/",
       width = 180,
       height = 180,
       units = 'mm')

##proteobacteria specific growth rates------
reg_all_slopes_chosen_silva_tax_specific_gr <- 
  reg_all_slopes_chosen_silva_tax_asv_reorder %>%
  filter(treatment == 'VL')

ridge_proteo <- 
  growthrates.distribution.tax.rank.ridges.divided(
    data = reg_all_slopes_chosen_silva_tax_specific_gr, 
    phylum_to_explore = 'Proteobacteria',
    title = 'Proteobacteria class',
    axis_y_title = '',
    x_c = 2,
    x_o = 4,
    x_f = 4,
    x_a = 4) %>%
  as_ggplot() ##estic afegint fct_rev perquè m'ordeni els gràfics de manera lògica alpha i gamma

ggsave('ridge_proteobacteria_specific_gr_v1.pdf', ridge_proteo, 
       path = "~/Documentos/Doctorat/REMEI/results/figures/raw/",
       width = 180,
       height = 180,
       units = 'mm')

reg_all_slopes_chosen_silva_tax %>%
  colnames()

##canvis en les distribucions entre CL i VL-------
asv_filtered <- reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  # group_by(asv_f) %>%
  # filter(n() > 2) %>%
  group_by(phylum_f) %>%
  mutate(counts = paste('n = ', n())) %>% #
  ungroup()

reg_all_slopes_chosen_silva_tax  %>%
  subset(phylum == 'Proteobacteria') %>%
  #subset(treatment %in% c('CL', 'VL')) %>%
  ggplot(aes(y = fct_rev(order_f), x = slope_chosen_days))+#, label = counts
  geom_density_ridges(aes(fill = treatment, group = order_f), scale = 1.0, alpha = 0.6,
                      jittered_points = TRUE,
                      point_shape = 21, point_size = 0.2, point_alpha = 0.0,
                      #position =
                      quantile_lines = TRUE,
                      quantile_fun = mean)+
  scale_fill_manual(values = palette_treatments_remei)+
  labs(y = '', fill= 'Phylum', x =expression("Growth rate day"^"-1"), title = 'ASV')+
  scale_x_continuous(limits = c(0,11))+
  #geom_text(nudge_x = 7.2, nudge_y = 0.35, check_overlap = TRUE, size = 3)+
  facet_grid(cols = vars(treatment), scales = 'free')+
  coord_cartesian(clip = "off") +
  theme_bw()+
  theme(legend.position = "none", strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 0),
        plot.margin= unit(c(1, 1, 0.2, 0.2), "lines"), panel.border = element_blank())
  
##falten números a asv level
test <- reg_all_slopes_chosen_silva_tax %>%
  subset(phylum == 'Proteobacteria') %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  group_by(asv_num) %>%
  filter(n() > 8) %>%
  mutate(counts = paste('n = ', n())) %>%
  group_by(family) %>%
  filter(n() > 10) %>%
  ungroup() %>%
  ggplot(aes(y = asv_f, x = slope_chosen_days, fill = class_f, label = counts))+
  # geom_histogram(aes(x=slope_chosen_days, group = asv_num), colour  = "#4F507F", fill = "#4F507F", 
  #                #adjust= 0.5, #bins = bin_width_a,
  #                binwidth = bin_width_a,
  #                show.legend = FALSE, alpha = 0.2)+
  geom_density_ridges(aes(fill = class_f), scale = 4, alpha = 0.6,
                      quantile_lines = TRUE,
                      quantile_fun = mean)+
  labs(x= expression("Growth rate day"^"-1"), y = "", title = 'ASV')+
  scale_x_continuous(limits = c(0,11))+
  scale_fill_manual(values = c("#c59a3e","#fcca46", "#806f52"))+
  facet_grid(rows = vars(family_f), scales = 'free', space = 'free', switch = 'y')+
  #scale_y_continuous(limits = c(y1_a,y2_a))+
  #stat_bin(binwidth = bin_width_a)+
  geom_text(nudge_x = 8.8, nudge_y = 0.3,  size = 3)+ #check_overlap = TRUE,
  theme_ridges(center_axis_labels = TRUE)+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #axis.text.y = element_text(size = 0), 
        strip.background = element_blank(), axis.ticks.y = element_blank(),
        strip.text.y = element_text(size = 10),
        plot.margin= unit(c(1.2, 1, 1.2, 1.2), "lines"))
   

reg_all_slopes_chosen_silva_tax %$%
  family_f %>%
  levels()
  
##proves a nivell de classe per endreçar millor els facet grids
 
  reg_all_slopes_chosen_silva_tax %>%
    subset(phylum == 'Proteobacteria') %>%
    filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
    filter(slope_chosen_days > 0) %>%
    group_by(family) %>%
    filter(n() > 10) %>%
    mutate(counts = paste('n = ', n())) %>%
    group_by(order) %>%
    filter(n() > 10) %>%
    ungroup() %>%
    ggplot(aes(y = family, x = slope_chosen_days, fill = phylum_f, label = counts))+ #interaction( family_f, order_f, sep = '\n')
    # geom_histogram(aes(x=slope_chosen_days, group = family), colour = "#FFDF26", fill = "#FFDF26", #adjust= 0.5, 
    #                # bins = bin_width_f,
    #                binwidth = bin_width_f,
    #                show.legend = FALSE, alpha = 0.2)+
    geom_density_ridges(aes(fill = phylum_f), scale = 3, alpha = 0.6,
                        quantile_lines = TRUE,
                        quantile_fun = mean)+
    labs(x= "")+
    #stat_bin(binwidth = bin_width_f)+
    scale_x_continuous(limits = c(0,11))+
    scale_fill_manual(values = palette_phylums)+
    facet_grid(rows = vars(order_f), scales = 'free', space = 'fixed')+
    geom_text(nudge_x = 8.8, nudge_y = 0.3, check_overlap = TRUE, size = 3)+
    # scale_y_continuous(limits = c(y1_f,y2_f))+
    theme_bw()+ 
    theme(strip.text.x = element_text(size = 0), axis.text.x = element_text(angle = 0), 
          axis.ticks = element_blank(), legend.position = "none",axis.title.y = element_text(size = 0),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_blank(),
          strip.text.y = element_blank(),
          plot.margin= unit(c(1, 0.2, 0.2, 0.2), "lines"))
  
##BACTEROIDOTA
ridge_bacterio <- growthrates.distribution.tax.rank.ridges.divided(
  data = reg_all_slopes_chosen_silva_tax, 
  phylum_to_explore = 'Bacteroidota',
  title = 'Bacteroidota \n
  Class',
  axis_y_title = '',
  x_c = 4,
  x_o = 10,
  x_f = 10,
  x_a = 8) %>%
  as_ggplot()

ggsave('ridge_bacteroidota_filtered_v3.pdf', ridge_bacterio, 
       path = "~/Documentos/Doctorat/REMEI/results/figures/raw/",
       width = 180,
       height = 180,
       units = 'mm')


ridge_plancto <- growthrates.distribution.tax.rank.ridges.divided(
  data = reg_all_slopes_chosen_silva_tax, 
  phylum_to_explore = 'Planctomycetota',
  title = 'Planctomycetota',
  axis_y_title = '',
  x_c = 4,
  x_o = 4,
  x_f = 4,
  x_a = 4) %>%
  as_ggplot()

ridge_actino <- growthrates.distribution.tax.rank.ridges.divided(
  data = reg_all_slopes_chosen_silva_tax, 
  phylum_to_explore = 'Actinobacteriota',
  title = 'Actinobacteriota',
  axis_y_title = '',
  x_c = 4,
  x_o = 4,
  x_f = 4,
  x_a = 4) %>%
  as_ggplot()

ridge_crena <- growthrates.distribution.tax.rank.ridges.divided(
  data = reg_all_slopes_chosen_silva_tax, 
  phylum_to_explore = 'Crenarchaeota',
  title = 'Crenarchaeota',
  axis_y_title = '',
  x_c = 4,
  x_o = 4,
  x_f = 4,
  x_a = 4) %>%
  as_ggplot()

ridge_verru <- growthrates.distribution.tax.rank.ridges.divided(
  data = reg_all_slopes_chosen_silva_tax, 
  phylum_to_explore = 'Verrucomicrobiota',
  title = 'Verrucomicrobiota',
  axis_y_title = '',
  x_c = 4,
  x_o = 4,
  x_f = 4,
  x_a = 4) %>%
  as_ggplot()

ridge_firmi <- growthrates.distribution.tax.rank.ridges.divided(
  data = reg_all_slopes_chosen_silva_tax, 
  phylum_to_explore = 'Firmicutes',
  title = 'Firmicutes',
  axis_y_title = '', 
  x_c = 4,
  x_o = 4,
  x_f = 4,
  x_a = 4) %>%
  as_ggplot()


# reg_all_slopes_chosen_silva_tax_filt %$%
#   phylum_f %>%
#   unique()
# 
# distribution_gr_plots_ridge_separated <- multi_panel_figure(columns = 4, rows = 8, width = 410, height = 797, 
#                                                   column_spacing = 0.1, unit = 'mm',
#                                                   row_spacing = 0.1, 
#                                                   panel_label_type = 'none' 
#                                                   )
# 
# distribution_gr_plots_ridge_separated  %<>%
#   #fill_panel(perc_gr_phylum, column = 1, row = 1:2) %<>%
#   fill_panel(ridge_proteo, column = 1:4, row = 1:3)   %<>%
#   fill_panel(ridge_bacterio, column = 1:4, row = 4:5) %<>%
#   fill_panel (ridge_plancto, column = 1:4, row = 6)  %<>% 
#   fill_panel (ridge_actino, column = 1:4, row = 7)  %<>% 
#   fill_panel(ridge_crena, column = 1:4, row = 8)  #%<>%
# 
# 
#   #fill_panel (ridge_firmi, column = 1:4, row = 8)
# 
# ggsave('distribution_gr_plots_ridge_separated_5most_freq.pdf', distribution_gr_plots_ridge_separated, 
#        path = "~/Documentos/Doctorat/REMEI/results/figures/raw/",
#        width = 410,
#        height = 797,
#        units = 'mm')

#plot for phylums with most of the GR calculated Bacteroidota & Proteobacteria--------
# distribution_gr_plots_ridge_top_phylums <- multi_panel_figure(columns = 1, rows = 3, width = 410, height = 597, 
#                                                               column_spacing = 0.1, unit = 'mm',
#                                                               row_spacing = 0.1, 
#                                                               panel_label_type = 'none')
# 
# distribution_gr_plots_ridge_top_phylums  %<>%
#   #fill_panel(perc_gr_phylum, column = 1, row = 1:2) %<>%
#   fill_panel(ridge_proteo, column = 1, row = 1:2)   %<>%
#   fill_panel(ridge_bacterio, column = 1, row = 3)   #%<>%
# 
# 
# #fill_panel (ridge_firmi, column = 1:4, row = 8)
# 
# ggsave('distribution_gr_plots_ridge_top_phylums.pdf', distribution_gr_plots_ridge_top_phylums, 
#        path = "~/Documentos/Doctorat/REMEI/results/figures/raw/",
#        width = 410,
#        height = 597,
#        units = 'mm')

####-----
reg_all_slopes_chosen_silva_tax %>%
  subset(phylum == 'Proteobacteria') %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  group_by(family_f) %>%
  filter(n() > 2) %>%
  ggplot(aes(y = interaction(class_f, family_f), x = slope_chosen_days, fill = phylum_f))+
  geom_density_ridges(aes(fill = phylum_f), scale = 3, alpha = 0.8)+
  scale_fill_manual(values = palette_class_assigned)+
  facet_grid(rows = vars(class), scales = 'free')+
  labs(y = 'Family', fill= 'Phylum', x =expression("Growth rate day"^"-1"))+
  scale_x_continuous(limits = c(0,11))+
  #scale_y_discrete(position = "right") +
  #coord_flip()+
  theme_bw()+
  theme(legend.position = "none", strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  
        plot.margin= unit(c(0.2, 1, 0.2, 1), "lines"))

##All phylums separated without asv level

##ASV level graph
  reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05) %>%
  filter(slope_chosen_days > 0) %>%
  group_by(asv_f) %>%
  filter(n() > 2) %>%
  ggplot(aes(y = asv_f, x = slope_chosen_days, fill = phylum_f))+
  geom_density_ridges(aes(fill = phylum_f), scale = 2, alpha = 0.8)+
  scale_fill_manual(values = palette_phylums)+
  labs(y = 'ASV', fill= 'Phylum', x =expression("Growth rate day"^"-1"))+
  theme_light()+
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 0), 
        plot.margin= unit(c(0.2, 0.2, 0.2, 0.2), "lines"))

##Què passa amb les alphaproteobacteria a asv level?
reg_all_slopes_chosen_silva_tax %>%
  subset(class == 'Alphaproteobacteria',
         slope_chosen_days != is.na(slope_chosen_days))%$%
  family %>%
  unique()

reg_all_slopes_chosen_silva_tax %>%
  subset(class == 'Alphaproteobacteria') %>%
  group_by(order_f) %>%
  ggplot(aes(x=slope_chosen_days, group = order_f, color = order_f))+ #, colour = class
  # geom_histogram(aes(x=slope_chosen_days, group = asv_num), colour  = "#4F507F", fill = "#4F507F", 
  #                #bins = bin_width_a,
  #                binwidth = 0.1,
  #                #adjust= 0.5, 
  #                show.legend = FALSE, alpha = 0.2)+
  #stat_bin(binwidth = bin_width_a)+
  #geom_freqpoly(aes(group = order_f), stat = 'bin', binwidth = 0.25)+
  geom_area(stat = 'bin', binwidth = 0.25, aes(fill = order_f), alpha = 0.6)+
  scale_x_continuous(limits = c(0,11))+
  scale_color_manual(values = palf_large(18))+
  scale_fill_manual(values = palf_large(18))+
  #scale_y_continuous(limits = c(y1,y2))+
  labs(x= expression("Growth rate day"^"-1"), y = "")+
  theme_bw()+ 
  theme(legend.position = "right", strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  
        plot.margin= unit(c(0.2, 1, 0.2, 1), "lines"))

reg_all_slopes_chosen_silva_tax %>%
  subset(class == 'Alphaproteobacteria') %>%
  group_by(family) %>%
  ggplot(aes(x=slope_chosen_days, group = family_f, color = family_f))+ #, colour = class
  # geom_histogram(aes(x=slope_chosen_days, group = asv_num), colour  = "#4F507F", fill = "#4F507F", 
  #                #bins = bin_width_a,
  #                binwidth = 0.1,
  #                #adjust= 0.5, 
  #                show.legend = FALSE, alpha = 0.2)+
  #stat_bin(binwidth = bin_width_a)+
  #geom_freqpoly(aes(group = family_f), stat = 'bin', binwidth = 0.25)+
  geom_area(stat = 'bin', binwidth = 0.25, aes(fill = family_f), alpha = 0.6)+
  scale_x_continuous(limits = c(0,11))+
  scale_color_manual(values = palf_large(48))+
  scale_fill_manual(values = palf_large(48))+
  #scale_y_continuous(limits = c(y1,y2))+
  labs(x= expression("Growth rate day"^"-1"), y = "")+
  theme_bw()+ 
  theme(legend.position = "right", strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  
        plot.margin= unit(c(0.2, 1, 0.2, 1), "lines"))

test <- reg_all_slopes_chosen_silva_tax %>%
  subset(slope_chosen_days != is.na(slope_chosen_days)) %>%
  group_by(phylum) %>%
  add_count(phylum_f, wt = slope_chosen_days)

#prova amb sankey diagram-------
library(viridis)
library(patchwork)
library(hrbrthemes)
library(circlize)
reg_all_slopes_chosen_silva_tax %>%
  subset(slope_chosen_days != is.na(slope_chosen_days)) %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  group_by(phylum_f) %>%
  filter(n() > 2) %>%
  mutate(mean_phylum = mean(slope_chosen_days)) %>%
  ungroup() %>%
  group_by(class) %>%
  mutate(mean_class = mean(slope_chosen_days)) %>%
  ungroup() %>%
  group_by(order) %>%
  mutate(mean_order = mean(slope_chosen_days)) %>%
  ungroup() %>%
  group_by(family) %>%
  mutate(mean_family = mean(slope_chosen_days)) %>%
  ungroup() %>%
  group_by(asv_num) %>%
  mutate(mean_asv = mean(slope_chosen_days)) %>%
  select(starts_with('mean'), 'phylum_f', 'class', 'order', 'family', 'asv_num') %>%
  pivot_longer(cols = starts_with('mean')) %>%
  group_by(phylum_f) %>%
  distinct(value, name, phylum_f) %>%
  ggplot(aes(fct_rev(as_factor(name)), value))+
  geom_point(aes(color = phylum_f), alpha = 0.9, position = position_jitter(0.3))+
  #geom_boxplot(alpha= 0.05)+
  geom_violin(alpha = 0.2, draw_quantiles = c(0.25, 0.75))+
  scale_color_manual(values = palette_phylums)+
  labs(y = 'Mean growth rate at different\ntaxonomic ranks', x= 'Taxonomic rank', color = 'Phylum')+
  scale_x_discrete(labels = labels_mean_rank)+
  coord_flip()+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = 'none', axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14))


#Growing exponentially, MAXIMAL GROWTH RATES CALCULATION------
#T3 i T4 comparison slopes
#Lag at the beginning?
## Example of GR calculation for some ASV (Materials and Methods) -----
#grup asv amb pvalue més significatiu: asv2, asv23, asv3, asv117, asv4 i asv10
pseudoabund_df_wide_fc_silva$season <-  pseudoabund_df_wide_fc_silva$season %>% factor( 
                                               levels = (c('Winter', 'Spring', 'Summer', 'Fall', 'Early_fall')))
pseudoabund_df_wide_fc_silva$treatment <- pseudoabund_df_wide_fc_silva$treatment %>%
  factor(levels=(c("CL", "CD", "PL", "PD", "DL", "VL")))
##4 temps
#asv2 asv2_Bacteroidia_Flavobacteriaceae_Aurantivirga_NA
library(ggpubr)
t4_exemple <- 
  pseudoabund_df_wide_fc_silva %>%
  #select(1:5, 'asv81', 'asv53', 'asv10', 'asv654', 'asv4', 'asv114') %>% # 'asv1',
  #select(1:5, 'asv2', 'asv23', 'asv3', 'asv117', 'asv4', 'asv10') %>%
  select(1:5, 'asv2') %>%
  subset(season != "Early_fall" & 
         treatment == 'DL') %>%
  pivot_longer(cols = starts_with('asv'), names_to = 'asv_num', values_to = 'pseudoabund') %>%
  left_join(reg_all_slopes_chosen_silva_tax, by = 'asv_num') %>%
  ggplot(aes(hours.x, log(pseudoabund), color = season.x, shape = treatment.x))+
  geom_point(shape = 7, size = 3)+
  stat_poly_line(aes(group = season.x), mf.values = TRUE)+
  stat_cor(aes(label = paste(..p.label..)), label.x = 0, p.digits = 2, r.digits = 2)+
  #stat_pvalue_manual(aes(group = season.x), label = 'p.signif')+
  #stat_poly_eq(aes(group = season.x), p.digits = 2)+#label.x = 0
  #stat_regline_equation(aes(group = season.x), label.x = 0.1)+
  scale_y_continuous(labels = scales::scientific, breaks = c(5e+00, 1e+01), limits = c(2e+00, 1.25e+01))+
  #geom_smooth(method = 'lm', group = 'season.x')+ #, linetype = 'longdash', color = '#3d3b42'
  labs(y = 'log(ASV2 pseudoabundance)', x = 'Time (h)')+
  scale_color_manual(values = palette_seasons_4)+
  #facet_grid(asv_num ~., scales = 'free')+
  #facet_grid(vars(season.x), scales = 'fixed', switch = 'y')+
  theme_bw()+
  theme(strip.text = element_text(size = 12), axis.text.x = element_text(angle = 0), legend.position = "none",
        strip.background = element_blank(), strip.placement = 'outside',
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),# aspect.ratio = 1/2,
        plot.margin=unit(c(0.5,0.25,0.5,0.5),"cm")) #, plot.margin=unit(c(0.5,0,0.5,0),"cm") aspect.ratio = 1/2,
  
##3 temps
#library(ggpmisc)
t3_exemple <- 
  pseudoabund_df_wide_fc_silva %>%
  #select(1:5, 'asv81', 'asv53', 'asv10', 'asv654', 'asv4', 'asv114') %>% # 'asv1',
  #select(1:5, 'asv2', 'asv23', 'asv3', 'asv117', 'asv4', 'asv10') %>%
  select(1:5, 'asv2') %>%
  subset(season != "Early_fall" & 
           treatment == 'DL' &
           time != 't4')%>%
  pivot_longer(cols = starts_with('asv'), names_to = 'asv_num', values_to = 'pseudoabund') %>%
  left_join(reg_all_slopes_chosen_silva_tax, by = 'asv_num') %>%
  ggplot(aes(hours.x, log(pseudoabund), color = season.x))+
  geom_point(shape = 7, size = 3)+
  scale_y_continuous(labels = scales::scientific, breaks = c(5e+00, 1e+01), limits = c(2e+00, 1.25e+01))+
  stat_poly_line(aes(group = season.x))+
  #geom_smooth(method = 'lm')+ #, linetype = 'longdash', color = '#3d3b42'
  labs(y = 'log(ASV2 pseudoabundance)', x = 'Time (h)', color = 'Season', shape = 'Treatment')+
  scale_color_manual(values = palette_seasons_4)+
  #facet_grid(asv_num ~., scales = 'free')+
  #facet_grid(vars(treatment.x), scales = 'fixed')+
  scale_x_continuous(limits = c(0,50))+
  stat_cor(aes(label = paste(..p.label..)),label.x = 0, p.digits = 1, r.digits = 2)+
  #stat_poly_eq(aes(group = season.x))+
  scale_shape_manual(name = 'Treatment', labels = 'DL', values = 12)+
  #scale_y_continuous(breaks = c(50000, 10000))+ #breaks = c(1, 250, 500, 750, 1000)
  #geom_cor()+
  #coord_capped_cart(bottom='both')+ #, left='both'
  # stat_fit_glance(method = 'lm',
  #                 method.args = list(formula = formula),
  #                 geom = 'text',
  #                 aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "")))+
  theme_bw()+
  theme(strip.text = element_blank(), axis.text.x = element_text(angle = 0), legend.position = "right",
        axis.title.y = element_blank(), strip.background = element_blank(), axis.text.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        plot.margin=unit(c(0.5,0,0.5,0),"cm")) #aspect.ratio = 1/2, plot.margin=unit(c(0.5,0,0.5,0),"cm")


#no puc posar junts els treatments perquè els GR están calculats per treatments
#grid.arrange(t4_exemple, t3_exemple, ncol = 2 )

example_gr_calculation <- multi_panel_figure(columns = 2, rows = 1, width = 250, height = 150, 
                                            column_spacing = 0.0, unit = 'mm',
                                            panel_label_type = 'none')

example_gr_calculation  %<>%
  fill_panel(t4_exemple, column = 1, row = 1) %<>%
  fill_panel(t3_exemple, column = 2, row = 1)

ggsave('example_gr_calculation2.pdf', example_gr_calculation , 
       path = "~/Documentos/Doctorat/REMEI/results/figures/raw/",
       width = 250,
       height = 150,
       units = 'mm')

  ##provo de fer sobre el mateix gràfic les diferents regressions linials perquè
#quedi millor l'exemple



##Exploració de la phylogenia
reg_all_slopes_chosen_silva_tax %$%
  phylum %>%
  unique()
# 
# [1] "Cyanobacteria"     "Proteobacteria"    "Bacteroidota"      "Crenarchaeota"     "Planctomycetota"  
# [6] "Acidobacteriota"   "Verrucomicrobiota" "Nitrospinota"      "Actinobacteriota"  "Nitrospirota"     
# [11] "Firmicutes"        "Bdellovibrionota"  "Myxococcota"       "Desulfobacterota"  "Deinococcota"     
# [16] "Fusobacteriota"    "Chloroflexi"       "Spirochaetota"     "Campilobacterota"  "Abditibacteriota" 
# [21] "Latescibacterota"  "Methylomirabilota" "Halanaerobiaeota"  "Sumerlaeota"       "Calditrichota"    
# [26] "Gemmatimonadota" 

reg_all_slopes_chosen_silva_tax %>%
  subset(phylum == "Cyanobacteria") %$%
  class %>%
  unique()

reg_all_slopes_chosen_silva_tax %>%
  subset(phylum == "Acidobacteriota") %$%
  class %>%
  unique()

##relació num de treatments vs gr que presentes (classificació com a Bloomer?)

##Question: Do  bloomers respond to general environmental changes or are they specific? Can we classify them?


reg_all_slopes_chosen_silva_tax %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  group_by(phylum_f) %>%
  #summarize(n = n())
  filter(n() > 10) %>%
  ungroup() %>%
  group_by(asv_num, season) %>%
  mutate(n = n()) %>%
  distinct(asv_num, treatment, slope_chosen_days, season, n, phylum_f) %>%
  ggplot(aes(n, slope_chosen_days))+
  geom_point(position = position_jitter(0.35), alpha = 0.8, aes(color = season))+
  geom_violin(aes(group = n, color = season), alpha = 0, draw_quantiles = c(0.25, 0.75))+
  geom_smooth(method = 'loess',  aes(color = season))+#color = 'black',
  scale_color_manual(values = palette_seasons_4)+
  scale_y_continuous(expand = c(0, 0), limits = c(0,10.4))+
  labs(color = 'Season', x = 'Number of times each ASV appeared per treatment', y = expression("Growth rate day"^"-1"))+
  #facet_wrap(vars(treatment))+
  facet_grid(.~season)+#treatment~
  theme_bw()

scale_y_continuous(expand = c(0, 0), limits = c(0,10.4))


####REGRESSIONS FOR COMPLETE DATASET (+IN SITU DATA) (DISCARDED)---------------------------------------
#functions needed
source("src/filter.season.treatment.time.R") ###filtrar les pseudoabund per treatment season i time per calcular regressions
source("src/multiple.linear.regressions.R") ##per fer multiples regressions d'un dataframe
source("src/comparing.reg.R") ##per comparar les taxes de creiement entre 4 i 3 punts (les regressions)

#import data
#asv_tab pseudoabundances (recalculated with pseudoabundances from flow citometry data)
pseudoabund_df_wide_fc_silva <- read.table("data/intermediate_files/regressions_test_datasets/remiau_psueodabund_df_fc_wide_silva.txt", 
                                          header = TRUE, sep = "\t") ##312 mostres (hem filtrat per les que tenen menys de 5000 reads)

#objecte de phyloseq del que traiem la metadata i la tax  
remiau_fc_filt_silva <- readRDS("data/intermediate_files/remiau_phyloseq_silva_fc.rds") 

#substituint silva per GTBD o al revés fas els càlculs d'un o l'altre.
#metadata 
env <- remiau_fc_filt_silva@sam_data
env$treatment <- env$treatment %>% 
  factor(levels=(c("CL", "CD", "PL", "PD", "DL", "VL",'NA')))
env$season <- env$season %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Early_fall", "Fall")))


#asv_tab & metadata
pseudoabund_df_wide_m <- pseudoabund_df_wide_fc_silva %>%
  merge(env, by = c("treatment", "season", "time", "replicate", "hours.x"))

#asv's taxonomy this new column is important for ploting asv with their taxonomic information
rem_fc_filt_silva@tax_table <- rem_fc_filt_silva@tax_table %>% 
  mutate_tax_table(asv_num = str_c( 'asv' , 1:ncol(rem_fc_filt_silva@otu_table)))
tax_table <- rem_fc_filt_silva@tax_table %>%
  mutate_tax_table(tax_ed = (paste(asv_num, class, family, genus, species, sep = "_")))

pseudoabund_df_wide_m %>%
  colnames()

#1. Prepare  datasets for calculating regressions------
###crear una nova variable que combina time and replicate per no tenir problemes amb els temps repetits de les rèpliques.
data <- pseudoabund_df_wide_m %>% #pseudoabundances dataset form fc data in wide format with metadata.
  mutate(t_rep = paste(time, replicate, sep="_")) %>%
  filter(treatment != 'NA') ##eliminem les 4 mostres in situ perquè no serveixen per calcular les regressions.308 mostres

env <- env %>%
  as_tibble() %>%
  mutate(t_rep = paste(time, replicate, sep="_")) %>%
  filter(treatment != 'NA') ##tenim 4 mostres de més perquè no està filtrat per número de reads.


## 2. Filtering for each treatment, season and time
asv_filt_t0234 <- filter.season.treatment.time(data = data,  
                                              treatment_col = treatment, treatment_keep = "VL", 
                                              season_col = season, season_keep = "Winter",
                                              time_col = time, time_keep = c("t0", "t2", "t3", 't4')) %>%
  select(starts_with("asv"), #we clean metadata columns that we don't need now and keep only the column t_rep
         matches(c("t_rep", "hours.x"))) %>%
  mutate(across(everything(), ~replace(., . == 0, "1"))) %>% ##0 values will be substituted by 1 because it is mandatory for transforming pseudoabundances to log
  mutate(across(!t_rep, as.numeric))
 
  # asv_filt_t0234 %>%
  #    glimpse() #las dos ultimas columnas són t_rep i hours.x 
   # asv_filt_t0234 %>%
   #   dim() #5113

## 3. Calculating regressions (CHANGE MATRIX DIM FOR SILVA: 5111) 5113 menys les dues columnes de metadata
regress <- apply(X = asv_filt_t0234[,c(1:5111)], MARGIN = 2, FUN = multiple.linear.regressions, env = asv_filt_t0234$hours.x) ##DIMENSIONS OF THE MATRIX SHOULD BE ADAPTED TO DATA!!!!!
regress <- as.data.frame(t(regress)) %>% # le damos la vuelta y convertimos en data frame
  rownames_to_column(var = "asv_num")

# resultado %T>% #si surten NAs en el resultat és perquè les dimensions a funció d'apPDy están malament.
#   dim() %>%
#   head()

#cambiar nombre  del objeto en función del tratamiento y la estación
#NO ENTENC QUÊ FAIG AQUÍ ara no tinc el dataset asv_m perquè l'he eliminat. Com necessito que sigui per fer els CLots?
##vull un dataset amb tota la informació important el problema és que em quedo sense memòria de l'ordinador haig de simCLificar el procés
asv_m_reg_Winter_VL_t0234 <- filter.season.treatment.time(data = env,  
                                                         treatment_col = treatment, treatment_keep = "VL", 
                                                         season_col = season, season_keep = "Winter",
                                                         time_col = time, time_keep = c("t0", "t2", "t3")) %>%
  select(starts_with("asv"), ".sample", "treatment", "replicate",  "time", "season", "sample_name", "sample_code", "selected_file_name",
         "reads", "light_regime.y", "hours.y", "LNA", "HNA", "All", "BM", "Leu.PB", "SGR", "TD", "t_rep")  %>%
  left_join(asv_filt_t0234, by = "t_rep") %>%
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

##si vull fer totes les regressions amb diferents temps el que hauria de fer és afegir una nova columna que digui 
## de quin temps a quin temps es correspon cada regressió.

##unir tots els datasets
reg_all_t0234_silva <- bind_rows(reg_winter, 
                                 reg_spring,
                                 reg_summer,
                                 reg_early_fall,
                                 reg_fall)
## TINC UN PROBLEMA AMB PL-SPRING PERQUÈ LES REG HAN ESTAT CALCULADES DUES VEGADES I ESTAN DUPLICADES EN EL DATASET
#per això hi ha un pas on les elimino

#regressions 4 punts
write.csv(reg_all_t0234_silva, "data/intermediate_files/remiau_asv_reg_all_t0234_silva.csv")
#regressions 3 punts
write.table(reg_all_t023_silva, "data/intermediate_files/remiau_asv_reg_all_t023_silva.txt", sep = "\t")

##Comparison columna slope 4 temps i 3 temps, decidim quina ens quedem per cada asv ------
##import complete regression datasets 
###GTBD
# reg_all_t0234_gtbd <- read.csv("data/intermediate_files/regressions_test_datasets/asv_reg_all_t0234_gtbd.csv", sep="\t", header = T) %>%
#   as_tibble()
# reg_all_t023_gtbd <- read.csv("data/intermediate_files/regressions_test_datasets/asv_reg_all_t023_gtbd.csv", sep="\t", header = T) %>%
#   as_tibble()
###SILVA
reg_all_t0234_silva <-  read.table("data/intermediate_files/asv_reg_all_t0234_silva.txt", sep="\t", header = T) %>%
  as_tibble()
reg_all_t023_silva <-  read.table("data/intermediate_files/remiau_asv_reg_all_t023_silva.txt", sep="\t", header = T) %>%
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

#df1 ha de ser 3 temps i df2 4 temps pel problema del PL spring (sol·lucionar-ho)
reg_all_slopes_chosen_silva <- comparing.reg(df1 = reg_all_t023_silva, 
                                             df2 = reg_all_t0234_silva,
                                             treatment = treatment,
                                             season = season,
                                             asv_num = asv_num, 
                                             slope = slope, 
                                             slope.pval = slope.pval)


#Afegim taxonomia i fem boxplots generals per tenir un overview de la coherència que hi ha de les ----
##taxes de creiement per grups taxonomics
reg_all_slopes_chosen_silva_tax <- reg_all_slopes_chosen_silva %>%
  left_join(tax_table, by = "asv_num", copy = TRUE) 

write.csv(reg_all_slopes_chosen_silva_tax, "data/intermediate_files/remiau_reg_all_slopes_chosen_silva_tax.csv")

#es podria afegir la rel_abund per filtrar per exemple pels més abundants en general dels experiments.


####GRAPHS FOR COMPLETE DATASET (WITH IN SITU DATA)--------------------------------
## Plotting boxplots for general exploration of the dataset------------
## Import data
### Cambiar GTBD por SILVA en funció del dataset amb el que estigui treballant.
# reg_all_slopes_chosen_gtbd_tax <- read.csv("data/intermediate_files/reg_all_slopes_chosen_gtbd_tax.csv", sep=",") %>%
#   filter(season != "Early_fall")
reg_all_slopes_chosen_silva_tax <- read.csv("data/intermediate_files/remiau_reg_all_slopes_chosen_silva_tax.csv", sep=",") %>%
  filter(season != "Early_fall")

##filter by asv present at least >1% at some treatment or station
#rem_fc_relabun <- read_rds("data/rem_fc_relabun.rds")
remiau_relabun_melt <-  read.table("data/remiau_relabun_melt.txt", sep="\t", header = TRUE) %>%
  filter(season != "Early_fall")

remiau_relabun_melt %$%
  asv_num %>%
  unique() #5111 asv al meu dataset

abundant_asv <- remiau_relabun_melt %>% 
  filter(Abundance > 0.01) %>% #more than 1% of the community at some point
  dplyr::select(asv_num) %>%
  unique() %>%
  as_tibble()

# rem_relabun_melt_1perc <- rem_relabun_melt %>%
#   right_join(abundant_asv, by = "asv_num", copy = TRUE) #per comprovar amb quin % de relabund estic treballant

# reg_all_slopes_chosen_silva_tax_1perc %>%
#   colnames()

reg_all_slopes_chosen_silva_tax_1perc <- reg_all_slopes_chosen_silva_tax %>% 
  right_join(abundant_asv, by = "asv_num", copy = TRUE)

##reorder treatments and seasons for plots
reg_all_slopes_chosen_silva_tax_1perc$treatment <- reg_all_slopes_chosen_silva_tax_1perc$treatment %>% 
  factor(levels=(c("CL", "CD", "PL", "PD", "DL", "VL")))
reg_all_slopes_chosen_silva_tax_1perc$season <- reg_all_slopes_chosen_silva_tax_1perc$season %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))

##transform slopes /hours to /days
reg_all_slopes_chosen_silva_tax_1perc <- reg_all_slopes_chosen_silva_tax_1perc %>%
  mutate(slope_chosen_days = slope_chosen*24)

##Single ASV based growth rates (FIGURE 1)----
reg_all_slopes_chosen_silva_tax %>%
  colnames()

#statistics seasons----
reg_all_slopes_chosen_silva_tax_1perc %>%
  colnames()
anova<-aov(slope_chosen_days ~season, data = reg_all_slopes_chosen_silva_tax_1perc)
summary(anova) #p<0.05 ***
TukeyHSD(anova)#para ver los grupos que son significativamente distintos
aov_residuals <- residuals(object = anova)
# plot(anova, 1)
# plot(anova, 2)
shapiro.test(x = aov_residuals )#From the output, the p-value > 0.05 implying that the distribution of the data are 
#not significantly different from normal distribution. In other words, we can assume the normality.
##p-value < 2.2e-16 NO és normal 
library(FSA)##test the d
dunnTest(slope_chosen_days ~ season, data = reg_all_slopes_chosen_silva_tax_1perc,
         method= 'bonferroni')
results<-dunnTest(slope_chosen_days ~ season, data = reg_all_slopes_chosen_silva_tax_1perc,
                  method='bonferroni')$res
results<-results[results$P.adj<0.05,]
#Differents: fall-spring, fall-summer, fall-winter, spring-winter, summer-winter
X = results$P.adj <= 0.05
names(X) = gsub(" ",  "",  results$Comparison)
multcompLetters(X)

# Fall Spring Summer Winter 
# "a"    "b"    "b"    "c" 
#statistical groups for non-parammetric test
letters_seas <- data.frame(multcompLetters(X)['Letters'])

#plot seasons boxplot----
seas <- reg_all_slopes_chosen_silva_tax_1perc %>%
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
  geom_text(data = letters_seas, aes(y = 9, x = row.names(letters_seas), label = Letters),
            position = position_nudge(x = 0.25), hjust = 0.75, color = "black")+
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

#statistics treatments----
##Treatments
anova<-aov(slope_chosen_days ~ treatment, data = reg_all_slopes_chosen_silva_tax_1perc)
summary(anova) #p<0.05 ***
TukeyHSD(anova)#para ver los grupos que son significativamente distintos
aov_residuals <- residuals(object = anova)
# plot(anova, 1)
# plot(anova, 2)
shapiro.test(x = aov_residuals )#From the output, the p-value > 0.05 implying that the distribution of the data are #not significantly different from normal distribution. In other words, we can assume the normality.
##p-value < 2.2e-16 NO és normal 
#si no es normal:
#hago test no parametrico
kruskal.test(slope_chosen_days ~ treatment, data = reg_all_slopes_chosen_silva_tax_1perc) #
#si sale p<0.05 hago dunn test para ver cuales son significativamente distintos
library(FSA)##test the d
dunnTest(slope_chosen_days ~ treatment, data = reg_all_slopes_chosen_silva_tax_1perc,
         method= 'bonferroni')
results<-dunnTest(slope_chosen_days ~ treatment, data = reg_all_slopes_chosen_silva_tax_1perc,
                  method='bonferroni')$res
#results<-results[results$P.adj<0.05,]
##Diferents: CD-PD, PD-VL
X = results$P.adj <= 0.05

names(X) = gsub(" ",  "",  results$Comparison)
multcompLetters(X)
# CD   CL   DL   PD   PL   VL 
# "a" "ab" "ab"  "b" "ab"  "a"
letters_treat <- data.frame(multcompLetters(X)['Letters'])

treat <- reg_all_slopes_chosen_silva_tax_1perc %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  filter(season != "Early_fall") %>%
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
  geom_text(data = letters_treat, aes(y = 9, x = row.names(letters_treat), label = Letters),
            position = position_nudge(x = 0.25), hjust = 0.78, color = "black")+
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

#library(gridExtra)
gr_asv_dapis <- grid.arrange(treat, gr_dapis_treat, seas, gr_dapis_seas)
# 
ggsave('GR_asv_dapis_stat_remiau.pdf', gr_asv_dapis, path = "~/Documentos/Doctorat/REMEI/results/figures/raw/",
       height = 175,
       width = 250,
       units = 'mm')

##Mean growth rate by season and treatment considering non-growing prokaryotes FIGURE 1. C-----
seas_bulk_com <- reg_all_slopes_chosen_silva_tax_1perc %>%
  filter(season != "Early_fall") %>%
  mutate(slope_chosen_days_all = case_when(is.na(slope_chosen_days) ~ '0',
                                           slope_chosen_days < 0 ~ '0',
                                           pvalue_slope_chosen > 0.05 ~ '0',
                                           TRUE ~ as.character(slope_chosen_days))) %>%
  ggplot(aes(season, as.numeric(slope_chosen_days_all)))+ #
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


treat_bulk_com <- reg_all_slopes_chosen_silva_tax_1perc %>%
  filter(season != "Early_fall") %>%
  mutate(slope_chosen_days_all = case_when(is.na(slope_chosen_days) ~ '0',
                                           slope_chosen_days < 0 ~ '0',
                                           pvalue_slope_chosen > 0.05 ~ '0',
                                           TRUE ~ as.character(slope_chosen_days))) %>%
  ggplot(aes(treatment, as.numeric(slope_chosen_days_all)))+ #
  geom_point(aes(color = season), alpha = 0.7, position = position_jitter(0.25))+
  #geom_boxplot(alpha= 0.5)+
  geom_violin(alpha = 0.2, draw_quantiles = c(0.25, 0.75) )+
  labs(color = "Season")+
  stat_summary(fun = "mean", geom = "crossbar", 
               width = 0.6, colour = "black")+
  # geom_text(data = letters_seas, aes(y = 9, x = row.names(letters_seas), label = Letters),
  #           position = position_nudge(x = 0.2), hjust = 0.7, color = "black")+
  labs(x= "Treatment", y = expression("Growth rate day"^"-1"))+  #(μ) 
  #scale_x_discrete()+
  scale_y_continuous(expand = c(0, 0), limits = c(0,10.4))+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_seasons_4)+
  #facet_wrap(.~phylum, scales = "free")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

##mean DAPIS vs mean ASV (filtered & bulk community) ------
##transform slopes /hours to /days
reg_all_slopes_chosen_silva_tax <- reg_all_slopes_chosen_silva_tax %>%
  mutate(slope_chosen_days = slope_chosen*24)

general_mean_gr_asv <- reg_all_slopes_chosen_silva_tax_1perc %>%
  filter(season != "Early_fall") %>%
  mutate(slope_chosen_days_all = case_when(is.na(slope_chosen_days) ~ '0',
                                           slope_chosen_days < 0 ~ '0',
                                           pvalue_slope_chosen > 0.05 ~ '0',
                                           TRUE ~ as.character(slope_chosen_days))) %>%
  #group_by(treatment, season) %>%
  mutate(mean_asv = mean(as.numeric(slope_chosen_days_all))) %>%
  distinct(mean_asv)

reg_all_slopes_chosen_silva_tax %>%
  colnames()
mean_gr_asv <- reg_all_slopes_chosen_silva_tax %>%
  filter(season != "Early_fall") %>%
  mutate(slope_chosen_days_all = case_when(is.na(slope_chosen_days) ~ '0',
                                           slope_chosen_days < 0 ~ '0',
                                           pvalue_slope_chosen > 0.05 ~ '0',
                                           TRUE ~ as.character(slope_chosen_days))) %>%
  group_by(treatment, season) %>%
  mutate(mean_asv = mean(as.numeric(slope_chosen_days_all))) %>%
  distinct(season, treatment, mean_asv)

##table mean gr general community with 0 (ho faig amb tota la comunitat general o només amb 1perc com en els gràfics)
mean_gr_asv <- reg_all_slopes_chosen_silva_tax_1perc %>%
  filter(season != "Early_fall") %>%
  mutate(slope_chosen_days_all = case_when(is.na(slope_chosen_days) ~ '0',
                                           slope_chosen_days < 0 ~ '0',
                                           pvalue_slope_chosen > 0.05 ~ '0',
                                           TRUE ~ as.character(slope_chosen_days))) %>%
  group_by(season) %>%
  mutate(mean_asv = mean(as.numeric(slope_chosen_days_all)),
         sd_asv = sd(as.numeric(slope_chosen_days_all))) %>%
  distinct(season, mean_asv, sd_asv)


general_mean_gr_dapis <- GR_dapis_OS_filt %>%
  #group_by(treatment, season) %>%
  mutate(mean_dapis = mean(as.numeric(PRO))) %>%
  distinct(mean_dapis)

mean_gr_dapis <- GR_dapis_OS_filt %>%
  group_by(treatment, season) %>%
  mutate(mean_dapis = mean(as.numeric(PRO))) %>%
  distinct(season, treatment, mean_dapis)

mean_gr_asv_dapis <-  mean_gr_asv %>%
  left_join(mean_gr_dapis)

rel_gr_dapis_asv_bulk <- mean_gr_asv_dapis %>%
  as_tibble() %>%
  ggplot(aes(mean_asv, mean_dapis))+
  #geom_smooth(method = 'lm', se = FALSE, color = 'grey')+
  geom_point(aes(shape = treatment, color = season, size = 2))+
  #scale_fill_manual(values = palette_seasons_4)+
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  stat_poly_eq(color = 'black')+
  guides(size = 'none')+
  labs(x = 'Mean growth rates ASV based \n
       (whole community)', y = 'Mean growth rates DAPI based', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (3/4), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

##asv gr mean filtered by negative values, NAs and unsignificant
general_mean_gr_asv_filt <- reg_all_slopes_chosen_silva_tax_1perc %>%
  filter(season != "Early_fall") %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  #group_by(treatment, season) %>%
  mutate(mean_asv = mean(as.numeric(slope_chosen_days))) %>%
  distinct(mean_asv)

mean_gr_asv_filt <- reg_all_slopes_chosen_silva_tax_1perc %>%
  filter(season != "Early_fall") %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  group_by(treatment, season) %>%
  mutate(mean_asv = mean(as.numeric(slope_chosen_days))) %>%
  distinct(season, treatment, mean_asv)

mean_gr_asv_filt_dapis <-  mean_gr_asv_filt %>%
  left_join(mean_gr_dapis)

rel_gr_dapis_asv <- mean_gr_asv_filt_dapis %>%
  as_tibble() %>%
  ggplot(aes(mean_asv, mean_dapis))+
  #geom_smooth(method = 'lm', se = FALSE, color = 'grey')+
  geom_point(aes(shape = treatment, color = season, size = 2))+
  #scale_fill_manual(values = palette_seasons_4)+
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0))+
  scale_x_continuous(limits = c(0, 4), expand = c(0, 0))+
  scale_color_manual(values = palette_seasons_4)+
  stat_poly_line(color = 'black')+
  #stat_poly_eq(color = 'black')+
  stat_cor(aes(label = paste(..p.label..)),label.x = 0.5, p.digits = 1, r.digits = 2)+
  guides(size = 'none')+
  labs(x = 'Mean growth rates ASV based \n
       (growing community)', y = 'Mean growth rates DAPI based', color = 'Season', shape = 'Treatment')+
  theme_bw()+
  theme(aspect.ratio = (4/3), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

##general means graphic
# cbind(general_mean_gr_asv, general_mean_gr_asv_filt, general_mean_gr_dapis) %>%
#   t() %>%
#   as.data.frame() %>%
#   rownames_to_column(var = 'mean_gr') %>%
#   as_tibble() %>%
#   ggplot(aes(mean_gr, V1))+
#   geom_point()+
#   theme_bw()

##FIGURE panel construction----
# gr_asv_community <- grid.arrange(treat_bulk_com, seas_bulk_com)
# gr_asv_dapis_relations <- grid.arrange(rel_gr_dapis_asv, rel_gr_dapis_asv_bulk, nrow = 1)
# 
# gr_asv_dapis_relations <- grid.arrange(treat, seas, gr_dapis_treat, gr_dapis_seas,  treat_bulk_com, seas_bulk_com,
#                                        rel_gr_dapis_asv, rel_gr_dapis_asv_bulk, ncol = 2)

# ggsave('gr_asv_community.pdf', gr_asv_community, path = "~/Documentos/Doctorat/REMEI/results/figures/raw/",
#        height = 200,
#        width = 200,
#        units = 'mm')
# 
# ggsave('gr_asv_dapis_relations2.pdf', gr_asv_dapis_relations, path = "~/Documentos/Doctorat/REMEI/results/figures/raw/",
#        height = 125,
#        width = 200,
#        units = 'mm')

ggsave('gr_asv_dapis_remiau.pdf', rel_gr_dapis_asv, path = "~/Documentos/Doctorat/REMEI/results/figures/raw/",
       height = 125,
       width = 125,
       units = 'mm')

#Box plots----
#Single ASV based mean growth rate plots by different taxonomic ranks (Suppelementary)-------
reg_all_slopes_chosen_silva_tax %>%
  filter(slope.pval.x >0.05 & slope_chosen_days > 0) %>%
  summarize(slope_chosen_days_mean = mean(slope_chosen_days), #0.89
            slope_chosen_days= sd(slope_chosen_days)) #0.78
#mean GR by class
reg_all_slopes_chosen_silva_tax %>%
  filter(slope.pval.x >0.05 & slope_chosen_days > 0) %>% 
  group_by(class) %>%
  summarize(slope_chosen_days_mean = mean(slope_chosen_days), 
            slope_chosen_days= sd(slope_chosen_days))

##Order axis x by higher taxonomy levels (of example: order Classes by Phylum) 1%-----
reg_all_slopes_chosen_silva_tax_1perc %>%
  glimpse()

reg_all_slopes_chosen_silva_tax_1perc <- reg_all_slopes_chosen_silva_tax_1perc %>%
  mutate(domain_f = as_factor(domain),
         phylum_f = as_factor(phylum),
         class_f = as_factor(class),
         order_f = as_factor(order),
         family_f = as_factor(family),
         genus_f = as_factor(genus),
         asv_f = as_factor(asv_num))

reg_all_slopes_chosen_silva_tax_1perc %>%
  colnames()

reg_all_slopes_chosen_silva_tax_1perc$phylum_f <-  factor(reg_all_slopes_chosen_silva_tax_1perc$phylum_f, 
                                                          levels=unique(
                                                            reg_all_slopes_chosen_silva_tax_1perc$phylum_f[
                                                              order(reg_all_slopes_chosen_silva_tax_1perc$domain_f)]), 
                                                          ordered=TRUE)

reg_all_slopes_chosen_silva_tax_1perc$class_f <-  factor(reg_all_slopes_chosen_silva_tax_1perc$class_f, 
                                                         levels=unique(
                                                           reg_all_slopes_chosen_silva_tax_1perc$class_f[
                                                             order(reg_all_slopes_chosen_silva_tax_1perc$domain_f,
                                                                   reg_all_slopes_chosen_silva_tax_1perc$phylum_f)]), 
                                                         ordered=TRUE)
reg_all_slopes_chosen_silva_tax_1perc$order_f <-  factor(reg_all_slopes_chosen_silva_tax_1perc$order_f, 
                                                         levels=unique(
                                                           reg_all_slopes_chosen_silva_tax_1perc$order_f[
                                                             order(reg_all_slopes_chosen_silva_tax_1perc$domain_f,
                                                                   reg_all_slopes_chosen_silva_tax_1perc$phylum_f,
                                                                   reg_all_slopes_chosen_silva_tax_1perc$class_f)]), 
                                                         ordered=TRUE)

reg_all_slopes_chosen_silva_tax_1perc$family_f <-  factor(reg_all_slopes_chosen_silva_tax_1perc$family_f, 
                                                          levels=unique(
                                                            reg_all_slopes_chosen_silva_tax_1perc$family_f[
                                                              order(reg_all_slopes_chosen_silva_tax_1perc$domain_f,
                                                                    reg_all_slopes_chosen_silva_tax_1perc$phylum_f,
                                                                    reg_all_slopes_chosen_silva_tax_1perc$class_f,
                                                                    reg_all_slopes_chosen_silva_tax_1perc$order_f)]), 
                                                          ordered=TRUE)

reg_all_slopes_chosen_silva_tax_1perc$genus_f <-  factor(reg_all_slopes_chosen_silva_tax_1perc$genus_f, 
                                                         levels=unique(
                                                           reg_all_slopes_chosen_silva_tax_1perc$genus_f[
                                                             order(reg_all_slopes_chosen_silva_tax_1perc$domain_f,
                                                                   reg_all_slopes_chosen_silva_tax_1perc$phylum_f,
                                                                   reg_all_slopes_chosen_silva_tax_1perc$class_f,
                                                                   reg_all_slopes_chosen_silva_tax_1perc$order_f,
                                                                   reg_all_slopes_chosen_silva_tax_1perc$family_f)]), 
                                                         ordered=TRUE)
reg_all_slopes_chosen_silva_tax_1perc$asv_f <-  factor(reg_all_slopes_chosen_silva_tax_1perc$asv_f, 
                                                       levels=unique(
                                                         reg_all_slopes_chosen_silva_tax_1perc$asv_f[
                                                           order(reg_all_slopes_chosen_silva_tax_1perc$domain_f,
                                                                 reg_all_slopes_chosen_silva_tax_1perc$phylum_f,
                                                                 reg_all_slopes_chosen_silva_tax_1perc$class_f,
                                                                 reg_all_slopes_chosen_silva_tax_1perc$order_f,
                                                                 reg_all_slopes_chosen_silva_tax_1perc$family_f,
                                                                 reg_all_slopes_chosen_silva_tax_1perc$genus_f)]), 
                                                       ordered=TRUE)

##dataset without filtering by 1%------
reg_all_slopes_chosen_silva_tax <- reg_all_slopes_chosen_silva_tax %>%
  mutate(domain_f = as_factor(domain),
         phylum_f = as_factor(phylum),
         class_f = as_factor(class),
         order_f = as_factor(order),
         family_f = as_factor(family),
         genus_f = as_factor(genus),
         asv_f = as_factor(asv_num))


reg_all_slopes_chosen_silva_tax %>%
  colnames()

reg_all_slopes_chosen_silva_tax$phylum_f <-  factor(reg_all_slopes_chosen_silva_tax$phylum_f, 
                                                    levels=unique(
                                                      reg_all_slopes_chosen_silva_tax$phylum_f[
                                                        order(reg_all_slopes_chosen_silva_tax$domain_f)]), 
                                                    ordered=TRUE)
reg_all_slopes_chosen_silva_tax$phylum_f %>% levels()

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

reg_all_slopes_chosen_silva_tax$phylum_f %>%
  levels()
##exploration at class level---
### el vull al supplementary
#by seasons
reg_all_slopes_chosen_silva_tax_1perc %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  # group_by(phylum_f) %>%
  # filter(n() > 2) %>%
  #filter(pseudoabundance > 1000) %>%
  ggplot(aes(x=interaction(class_f,phylum_f, sep = '\n'), slope_chosen_days, group = class_f))+ #
  geom_point(aes(color = season), alpha = 0.8, position = position_jitter(0.5))+#
  geom_boxplot(alpha= 0.1)+
  labs(color = "Season")+
  #geom_violin(alpha = 0.5, draw_quantiles = c(0.25, 0.75))+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.3,
               colour = "black")+
  labs(x= "Class\n Pylum", y = expression("Growth rate (μ) day"^"-1"))+
  # geom_label(position = position_stack(vjust = 0.99))+
  #geom_text(position = position_dodge(width = 1), aes(x=class, y=-50))+
  #scale_y_continuous(limits = c(0, 10))+
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_seasons_4)+#palf_large(94)+
  #coord_flip()+
  #facet_wrap(.~treatment, scales = "fixed")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 4), 
        axis.text.x = element_text(angle = 0), 
        legend.position = "right",
        axis.text.y = element_text(angle = 0))

#by treatmens
reg_all_slopes_chosen_silva_tax_1perc %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  #filter(pseudoabundance > 1000) %>%
  ggplot(aes(x=interaction(class_f,phylum_f, sep = '\n'), slope_chosen_days, group = class_f))+ #
  geom_point(aes(color = treatment), alpha = 0.8, position = position_jitter())+#
  geom_boxplot(alpha= 0.3)+
  labs(color = "Treatment")+
  #geom_violin(alpha = 0.5, draw_quantiles = c(0.25, 0.75))+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.3,
               colour = "black")+
  labs(x= "Class\n Pylum", y = expression("Growth rate (μ) day"^"-1"))+
  # geom_label(position = position_stack(vjust = 0.99))+
  #geom_text(position = position_dodge(width = 1), aes(x=class, y=-50))+
  #scale_y_continuous(limits = c(0, 10))+
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_treatments_remei)+#palf_large(94)+
  #coord_flip()+
  #facet_wrap(.~treatment, scales = "fixed")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 4), 
        axis.text.x = element_text(angle = 0), 
        legend.position = "right",
        axis.text.y = element_text(angle = 0))


###prova sense subset per 1% per plotejant menys punts----
reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  #filter(pseudoabundance > 1000) %>%
  ggplot(aes(x=interaction(class_f,phylum_f, sep = ' '), slope_chosen_days, group = class_f))+ #
  geom_point(aes(color = treatment), alpha = 0.8, position = position_jitter())+#
  geom_boxplot(alpha= 0.3)+
  labs(color = "Treatment")+
  #geom_violin(alpha = 0.5, draw_quantiles = c(0.25, 0.75))+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.3,
               colour = "black")+
  labs(x= "Class\n Pylum", y = expression("Growth rate (μ) day"^"-1"))+
  # geom_label(position = position_stack(vjust = 0.99))+
  #geom_text(position = position_dodge(width = 1), aes(x=class, y=-50))+
  #scale_y_continuous(limits = c(0, 10))+
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_treatments_remei)+#palf_large(94)+
  coord_flip()+
  #coord_flip()+
  #facet_wrap(.~treatment, scales = "fixed")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 4), 
        axis.text.x = element_text(angle = 0), 
        legend.position = "right",
        axis.text.y = element_text(angle = 0))

##plot at order level-----
### el vull al supplementary
###by seasons
reg_all_slopes_chosen_silva_tax_1perc %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  #filter(pseudoabundance > 1000) %>%
  ggplot(aes(interaction(order_f, class_f, phylum_f, sep = '\n'), slope_chosen_days))+ #
  geom_point(aes(color = season), alpha = 0.7, position = position_jitter())+
  geom_boxplot(alpha= 0.3)+
  #geom_violin(alpha = 0.5)+
  labs(x= "Order", y = expression("Growth rate (μ) day"^"-1"), color = 'Season')+
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_seasons_4)+#palf_large(94)
  #facet_wrap(.~class, scales = "free")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 10.5), 
        axis.text.x = element_text(angle = 60, size = 6, hjust = 0.9), 
        legend.position = "right")
###by treatments
reg_all_slopes_chosen_silva_tax_1perc %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  #filter(pseudoabundance > 1000) %>%
  ggplot(aes(interaction(order_f, class_f,phylum_f, sep = '\n'), slope_chosen_days))+ #
  geom_point(aes(color = treatment), alpha = 0.7, position = position_jitter())+
  geom_boxplot(alpha= 0.3)+
  #geom_violin(alpha = 0.5)+
  labs(x= "Order", y = expression("Growth rate (μ) day"^"-1"), color = 'Season')+
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_treatments_remei)+#palf_large(94)
  #facet_wrap(.~class, scales = "free")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 10.5), 
        axis.text.x = element_text(angle = 60, size = 6, hjust = 0.9), 
        legend.position = "right")

##at family level boxplots ----
###by seasons
reg_all_slopes_chosen_silva_tax_1perc %>%
  colnames()
reg_all_slopes_chosen_silva_tax_1perc %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  #filter(pseudoabundance > 1000) %>%
  ggplot(aes(interaction(family_f, order_f, phylum_f, sep = '\n'), slope_chosen_days))+ #
  geom_point(aes(color = season), alpha = 0.7, position = position_jitter())+
  geom_boxplot(alpha= 0.3)+
  #geom_violin(alpha = 0.5)+
  labs(x= "Order", y = expression("Growth rate (μ) day"^"-1"), color = 'Season')+
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_seasons)+#palf_large(94)
  #facet_wrap(.~class, scales = "free")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 10.5), 
        axis.text.x = element_text(angle = 60, size = 6, hjust = 0.9), 
        legend.position = "right")
###by treatments
reg_all_slopes_chosen_silva_tax_1perc %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  #filter(pseudoabundance > 1000) %>%
  ggplot(aes(interaction(family_f, order_f, phylum_f, sep = '\n'), slope_chosen_days))+ #
  geom_point(aes(color = treatment), alpha = 0.7, position = position_jitter())+
  geom_boxplot(alpha= 0.3)+
  #geom_violin(alpha = 0.5)+
  labs(x= "Order", y = expression("Growth rate (μ) day"^"-1"), color = 'Treatment')+
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_treatments_remei)+#palf_large(94)
  #facet_wrap(.~class, scales = "free")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 10.5), 
        axis.text.x = element_text(angle = 60, size = 6, hjust = 0.9), 
        legend.position = "right")


##exploration a nivell de classes concretes
#"Acidimiccrobiia", "Actinobaccteria", "Alphaproteobacteria", "Bacilli", "Bacteroidia", "Clostridia", "Cyanobacteriia",
#"Gammaproteobacteria",  "Nitrosphaeria", "Phycisphaerae", "Planctomycetes", "Verrucomicrobiae", "Vicinamibacteria"

reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  filter(class == "Actinobacteria") %>%
  ggplot(aes(fct_reorder(order, family), slope_chosen))+ #
  geom_point(aes(color = genus), alpha = 0.7, position = position_jitter())+
  #geom_boxplot(alpha= 0.5)+
  geom_violin(alpha = 0.5)+
  labs(x= "Class", y = expression("Growth rate (μ) day"^"-1"))+
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = new_palette_large)+
  #facet_wrap(.~phylum, scales = "free")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 90), legend.position = "none")

reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  filter(order == "Micrococcales") %>%
  ggplot(aes(fct_reorder(genus, family), slope_chosen))+ #
  geom_point(aes(color = genus), alpha = 0.7, position = position_jitter())+
  geom_boxplot(alpha= 0.5)+
  #geom_violin(alpha = 0.5)+
  labs(x= "Class", y = expression("Growth rate (μ) day"^"-1"))+
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = new_palette_large)+
  #facet_wrap(.~phylum, scales = "free")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 90), legend.position = "none")

reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  filter(class == "Alphaproteobacteria") %>%
  ggplot(aes(order, slope_chosen))+ #
  geom_point(aes(color = season), alpha = 0.9, position = position_jitter())+
  #geom_boxplot(alpha= 0.5)+
  geom_violin(alpha = 0.1)+
  labs(x= "Class", y = expression("Growth rate (μ) day"^"-1"))+
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_seasons_4)+
  #facet_wrap(.~season, scales = "fixed")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 12), axis.text.x = element_text(angle = 60, size = 10, hjust = +1), legend.position = "right")

reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  filter(class == "Gammaproteobacteria") %>%
  ggplot(aes(order, slope_chosen_days))+ #
  geom_point(aes(color = season), alpha = 0.9, position = position_jitter())+
  #geom_boxplot(alpha= 0.5)+
  geom_violin(alpha = 0.1)+
  labs(x= "Class", y = expression("Growth rate (μ) day"^"-1"))+
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_seasons_4)+
  #facet_wrap(.~season, scales = "fixed")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 12), axis.text.x = element_text(angle = 60, size = 10, hjust = +1), legend.position = "right")

reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  filter(family == "Hyphomonadaceae") %>%
  ggplot(aes(asv_num, slope_chosen))+ #
  geom_point(aes(color = season), alpha = 0.9, position = position_jitter())+
  geom_boxplot(alpha= 0.5)+
  #geom_violin(alpha = 0.1)+
  labs(x= "ASV", y = expression("Growth rate (μ) day"^"-1"))+
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_seasons_4)+
  #facet_wrap(.~phylum, scales = "free")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 90), legend.position = "right")

reg_all_slopes_chosen_silva_tax_1perc %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  #filter(order == "Rhodobacterales") %>%
  ggplot(aes(family, slope_chosen_days))+ #
  geom_point(aes(color = season), alpha = 0.9, position = position_jitter(width=0.2))+
  geom_boxplot(alpha= 0.5)+
  #geom_violin(alpha = 0.1)+
  labs(x= "ASV number", y = expression("Growth rate (μ) day"^"-1"))+
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_seasons_4)+
  #facet_wrap(.~season, scales = "free")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 90), legend.position = "bottom")

reg_all_slopes_chosen_silva_tax_1perc %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  filter(order == "SAR11 clade") %>%
  ggplot(aes(asv_num, slope_chosen_days))+ #
  geom_point(aes(color = season), alpha = 0.9, position = position_jitter(width= 0.2))+
  geom_boxplot(alpha= 0.5)+
  #geom_violin(alpha = 0.1, draw_quantiles = c(0.25, 0.5, 0.75))+ #trabar la manera d'afegir la mitja
  labs(x= "Family", y = expression("Growth rate (μ) day"^"-1"))+
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_seasons_4)+
  #facet_wrap(.~season, scales = "free")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 90), legend.position = "none")

reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  filter(order == c("SAR11 clade", "Alteromonadales", "Rhodobacterales")) %>%
  ggplot(aes(treatment, slope_chosen_days))+ #
  geom_point(aes(color = season), alpha = 0.9, position = position_jitter())+
  #geom_boxplot(alpha= 0.5)+
  geom_violin(alpha = 0.1, draw_quantiles = c(0.25, 0.75))+ #trabar la manera d'afegir la mitja
  labs(x= "Treatment", y = expression("Growth rate (μ) day"^"-1"))+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.2,
               colour = "black")+
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_seasons_4)+
  facet_wrap(.~order, scales = "fixed")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 90), legend.position = "bottom")

reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  filter(order == "SAR11 clade") %>%
  ggplot(aes(treatment, slope_chosen_days))+ #
  geom_point(aes(color = season), alpha = 0.9, position = position_jitter())+
  geom_boxplot(alpha= 0.5)+
  #geom_violin(alpha = 0.1, draw_quantiles = c(0.25, 0.75))+ #trabar la manera d'afegir la mitja
  labs(x= "Treatment", y = expression("Growth rate (μ) day"^"-1"))+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.2,
               colour = "black")+
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_seasons_4)+
  facet_wrap(.~family, scales = "fixed")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 90), legend.position = "bottom")

reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  filter(family == "Flavobacteriaceae") %>%
  ggplot(aes(treatment, slope_chosen_days))+ #
  geom_point(aes(color = season), alpha = 0.9, position = position_jitter())+
  geom_boxplot(alpha= 0.5)+
  #geom_violin(alpha = 0.1, draw_quantiles = c(0.25, 0.75))+ #trabar la manera d'afegir la mitja
  labs(x= "Treatment", y = expression("Growth rate (μ) day"^"-1"))+
  stat_summary(fun = "mean",
               geom = "crossbar", 
               width = 0.2,
               colour = "black")+
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_seasons_4)+
  facet_wrap(.~genus, scales = "fixed")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 90), legend.position = "bottom")

reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  filter(class == "Bacteroidia") %>%
  ggplot(aes(family, slope_chosen_days))+ #
  geom_point(aes(color = season), alpha = 0.9, position = position_jitter())+
  geom_boxplot(alpha= 0.5)+
  #geom_violin(alpha = 0.1, draw_quantiles = c(0.25, 0.5, 0.75))+ #trabar la manera d'afegir la mitja
  labs(x= "Family", y = expression("Growth rate (μ) day"^"-1"))+
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_seasons_4)+
  #facet_wrap(.~season, scales = "free")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 90), legend.position = "none")

reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  filter(family == "Flavobacteriaceae") %>%
  ggplot(aes(genus, slope_chosen))+ #
  geom_point(aes(color = season), alpha = 0.9, position = position_jitter())+
  geom_boxplot(alpha= 0.5)+
  #geom_violin(alpha = 0.1, draw_quantiles = c(0.25, 0.5, 0.75))+ #trabar la manera d'afegir la mitja
  labs(x= "Genus", y = expression("Growth rate (μ) day"^"-1"))+
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_seasons_4)+
  #facet_wrap(.~season, scales = "free")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 90), legend.position = "none")

reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  filter(genus == "NS5 marine group") %>%
  ggplot(aes(asv_num, slope_chosen))+ #
  geom_point(aes(color = season), alpha = 0.9, position = position_jitter())+
  geom_boxplot(alpha= 0.5)+
  #geom_violin(alpha = 0.1, draw_quantiles = c(0.25, 0.5, 0.75))+ #trabar la manera d'afegir la mitja
  labs(x= "Season", y = expression("Growth rate (μ) day"^"-1"))+
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_seasons_4)+
  #facet_wrap(.~season, scales = "free")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 90), legend.position = "none")

reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  filter(class == "Gammaproteobacteria") %>%
  ggplot(aes(family, slope_chosen_days))+ #
  geom_point(aes(color = season), alpha = 0.9, position = position_jitter())+
  geom_boxplot(alpha= 0.5)+
  #geom_violin(alpha = 0.1, draw_quantiles = c(0.25, 0.5, 0.75))+ #trabar la manera d'afegir la mitja
  labs(x= "Family", y = expression("Growth rate (μ) day"^"-1"))+
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_seasons_4)+
  #facet_wrap(.~season, scales = "fixed")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 90), legend.position = "none")

reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  filter(family == "Alteromonadaceae") %>%
  ggplot(aes(asv_num, slope_chosen_days))+ #
  geom_point(aes(color = season), alpha = 0.9, position = position_jitter())+
  geom_boxplot(alpha= 0.5)+
  #geom_violin(alpha = 0.1, draw_quantiles = c(0.25, 0.5, 0.75))+ #trabar la manera d'afegir la mitja
  labs(x= "ASV", y = expression("Growth rate (μ) day"^"-1"))+
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_seasons)+
  #facet_wrap(.~season, scales = "free")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 90), legend.position = "none")

reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  filter(family == "Pseudoalteromonadaceae") %>%
  ggplot(aes(asv_num, slope_chosen_days))+ #
  geom_point(aes(color = season), alpha = 0.9, position = position_jitter())+
  geom_boxplot(alpha= 0.5)+
  #geom_violin(alpha = 0.1, draw_quantiles = c(0.25, 0.5, 0.75))+ #trabar la manera d'afegir la mitja
  labs(x= "ASV", y = expression("Growth rate (μ) day"^"-1"))+
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_seasons_4)+
  #facet_wrap(.~season, scales = "free")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 90), legend.position = "bottom")

reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  filter(phylum == "Cyanobacteria") %>%
  ggplot(aes(season, slope_chosen_days))+ #
  geom_point(aes(color = treatment), alpha = 0.9, position = position_jitter())+
  geom_boxplot(alpha= 0.5)+
  #geom_violin(alpha = 0.1, draw_quantiles = c(0.25, 0.5, 0.75))+ #trabar la manera d'afegir la mitja
  labs(x= "Season", y = expression("Growth rate (μ) day"^"-1"))+
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_treatments_remei)+
  #facet_wrap(.~season, scales = "fixed")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 90), legend.position = "none")

reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  #filter(season ==c("Winter", "Spring", "Summer", "Fall")) %>%
  filter(order == c("Alteromonadales", "Campylobacter",
                    "Cellvibrionales", "Flavobacteriales", "Flavobact. NS",
                    "Oceanospirillales", "Rhodobacterales", "SAR11 clade", "SAR86",
                    "Vibrionales")) %>%
  ggplot(aes(treatment, slope_chosen_days))+ #
  geom_point(aes(color = season), alpha = 0.9, position = position_jitter())+
  geom_boxplot(alpha= 0.5)+
  #geom_violin(alpha = 0.1, draw_quantiles = c(0.25, 0.75))+ #trabar la manera d'afegir la mitja
  labs(x= "Season", y = expression("Growth rate (μ) day"^"-1"))+
  # stat_summary(fun = "mean",geom = "crossbar", 
  #              width = 0.2, colour = "black")+
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_seasons_4)+
  facet_grid(.~order, scales = "free")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), 
        axis.text.x = element_text(angle = 90), legend.position = "right")

##Effect of light in likely containing proteorhodopsin
##NOR5, Rhodobacteraeae
reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  filter(genus == "OM60(NOR5) clade", ) %>%
  #filter(family == c("Halieaceae", "Rhodobacteraceae")) %>%
  ggplot(aes(treatment, slope_chosen_days))+ #
  geom_point(aes(color = treatment), alpha = 0.9, position = position_jitter())+
  geom_boxplot(alpha= 0.5)+
  #geom_violin(alpha = 0.1, draw_quantiles = c(0.25, 0.5, 0.75))+ #trabar la manera d'afegir la mitja
  labs(x= "Treatment", y = expression("Growth rate (μ) day"^"-1"))+
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_treatments_remei)+
  facet_wrap(.~season, scales = "fixed")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 90), legend.position = "none")

reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  #filter(genus == "OM60(NOR5) clade", ) %>%
  filter(family == "Rhodobacteraceae") %>%
  ggplot(aes(treatment, slope_chosen_days))+ #
  geom_point(aes(color = treatment), alpha = 0.9, position = position_jitter())+
  geom_boxplot(alpha= 0.5)+
  #geom_violin(alpha = 0.1, draw_quantiles = c(0.25, 0.5, 0.75))+ #trabar la manera d'afegir la mitja
  labs(x= "Treatment", y = expression("Growth rate (μ) day"^"-1"))+
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_treatments_remei)+
  facet_wrap(.~season, scales = "fixed")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 90), legend.position = "none")

reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  #filter(genus == "OM60(NOR5) clade", ) %>%
  filter(order == "Bacteroidales") %>%
  ggplot(aes(treatment, slope_chosen_days))+ #
  geom_point(aes(color = treatment), alpha = 0.9, position = position_jitter())+
  geom_boxplot(alpha= 0.5)+
  #geom_violin(alpha = 0.1, draw_quantiles = c(0.25, 0.5, 0.75))+ #trabar la manera d'afegir la mitja
  labs(x= "Treatment", y = expression("Growth rate (μ) day"^"-1"))+
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palette_treatments_remei)+
  #facet_wrap(.~family~genus, scales = "fixed")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 90), legend.position = "none")


#busco gèneres amb major taxa de creixement
reg_all_slopes_chosen_silva_tax %>%
  #filter(slope_chosen > 0.2) %>%
  filter(pvalue_slope_chosen < 0.05) %$%
  phylum %>%
  unique

# [1] "Pseudoalteromonas" "Polaribacter"      "Dokdonia"          "Microbacterium"    "Glaciecola"       
# [6] "Vibrio"            "Alteromonas"       "TMED14"            "Rubritalea"        NA                 
# [11] "UBA8087"           "AG-422-B15"   

plot.tax.filtered(data = reg_all_t0234_tax,
                  value_slope_pval = 5, 
                  value_intercept_pval = 5, 
                  tax_level = genus, 
                  select_tax_level  = "Sphingomonas")

##per families
reg_all_slopes_chosen_silva_tax %>%
  filter(slope_chosen_days > 0.2) %>%
  filter(pvalue_slope_chosen < 0.05 ) %$% #& intercept.pval < 0.05 (NO ES FILTRA NO SERVEIX)
  order %>%
  unique
# 
# [1] "Alteromonadaceae"   "Flavobacteriaceae"  "Microbacteriaceae"  "Vibrionaceae"       "Schleiferiaceae"   
# [6] "Akkermansiaceae"    "Rhodobacteraceae"   "SM1A02"             "AG-422-B15"         "Enterobacteriaceae"  

plot.tax.filtered_colors <- function(data, value_slope_pval, tax_level, select_tax_level){
  my_var1 <- rlang::enquo(tax_level)
  data %>%
    filter(pvalue_slope_chosen < value_slope_pval ) %>% #& intercept.pval < value_intercept_pval
    filter(!!my_var1 == select_tax_level) %>%
    ggplot(aes(hours, pseudoabundance, color = tax_ed, shape = season))+
    geom_point()+
    labs(x= "Time (h)", y = "log(pseudoabundance)")+
    geom_smooth(method = "lm")+
    scale_color_manual(values = palf_large(255))+
    facet_grid(.~treatment~season, scales = "fixed")+
    guides(color=guide_legend(ncol = 1))+
    theme_bw()+
    theme(strip.text.x = element_text(size = 14))
}

plot.tax.filtered.boxplot <- function(data, value_slope_pval){
  my_var1 <- data %>%
    filter(pvalue_slope_chosen < value_slope_pval ) %>% #& intercept.pval < value_intercept_pval
    #filter(!!my_var1 == select_tax_level) %>%
    ggplot(aes(family, slope, color = family, shape = season))+
    geom_boxplot()+
    labs(x= "Time (h)", y = "log(pseudoabundance)")+
    geom_smooth(method = "lm")+
    scale_color_manual(values = palf_large(50))+
    #facet_grid(.~treatment~season, scales = "fixed")+
    guides(color=guide_legend(ncol = 1))+
    theme_bw()+
    theme(strip.text.x = element_text(size = 14))
}

reg_all_t0234_tax %>%
  colnames()
plot.tax.filtered.boxplot(data = reg_all_t0234_tax,
                          value_slope_pval = 0.05)

reg_all_slopes_chosen_silva_tax %>%
  filter(slope.pval < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope > 0) %>%
  #filter(pseudoabundance > 500) %>%
  ggplot(aes(genus, slope, color = order))+
  geom_boxplot()+
  geom_point()+
  labs(x= "Genus", y = "Slope")+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palf_large(153))+
  facet_wrap(.~class, scales = "free")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+
  theme(strip.text.x = element_text(size = 5), axis.text.x = element_text(angle = 45), legend.position = "none")

reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  #filter(pseudoabundance > 1000) %>%
  ggplot(aes(genus, slope, color = family))+
  geom_boxplot()+
  geom_point()+
  labs(x= "Genus", y = "Slope")+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palf_large(140))+
  #facet_grid(.~treatment~season, scales = "fixed")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 90), legend.position = "none")

reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  #filter(pseudoabundance > 1000) %>%
  ggplot(aes(genus, slope_chosen, color = order))+
  geom_boxplot()+
  geom_point()+
  labs(x= "Genus", y = "Slope")+
  #scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palf_large(201))+
  facet_wrap(.~phylum, scales = "free")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 90), legend.position = "none")

reg_all_slopes_chosen_silva_tax %>%
  colnames()

reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05 ) %$% #& intercept.pval < 0.05 (NO ES FILTRA NO SERVEIX)
  class %>%
  unique
# 
# "Cyanobacteriia"      "Bacteroidia"         "Alphaproteobacteria" "Gammaproteobacteria" "Verrucomicrobiae"   
# "SAR324"              "Nitrososphaeria"     "Acidimicrobiia"      "Deinococci"          "Poseidoniia"        
#  "Planctomycetes"      "Bacteriovoracia"     "Actinomycetia"       "UBA2968"             "WGA-4E"             
#  "UBA1144"             "Clostridia"          "Phycisphaerae"       "UBA8108"             "Dehalococcoidia"    
#  "Desulfitobacteriia"  "Bacilli"             "Rhodothermia"        "Nitrospiria"         "UBA9160"            
# "Vicinamibacteria"    "Marinisomatia"       "Blastocatellia"      "Oligoflexia"         "Nitrospinia"        
#  "Desulfovibrionia"    "Bradimonadia"        "Fusobacteriia"       "Bin61"               "Binatia"            
# "Vampirovibrionia"    "Lentisphaeria"       "GCA-002687715"       "Zetaproteobacteria"  "Negativicutes"      
#  "Kiritimatiellae"     "Gemmatimonadetes"    "Desulfobulbia"       "Myxococcia"          "Desulfuromonadia"   
#  "Anaerolineae"        "UBA1135"             "Coriobacteriia"

comp_df_tax$treatment <- comp_df_tax$treatment %>% 
  factor(levels=(c("CL", "CD", "PL", "PD", "DL", "VL")))
comp_df_tax$season <- comp_df_tax$season %>% 
  factor(levels=(c("Winter", "Spring", "Summer","Early_fall", "Fall")))

reg_all_slopes_chosen_silva_tax %>%
  colnames()

reg_all_slopes_chosen_silva_tax %>%
  colnames()

reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen > 0) %>%
  #subset(class == "Gammaproteobacteria") %>%
  #filter(pseudoabundance > 1000) %>%
  ggplot(aes(family, slope_chosen))+ #, color = asv_num
  geom_boxplot()+
  geom_point(alpha= 0.6)+
  labs(x= "Phylogeny", y = "Slope")+
  scale_x_discrete()+
  #geom_smooth(method = "lm")+
  scale_color_manual(values = palf_large(236))+
  #facet_wrap(.~season, scales = "fixed")+
  #guides(color=guide_legend(ncol = 1))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 5), axis.text.x = element_text(angle = 90, size = 5), legend.position = "none")

SAR11 <- reg_all_t0234_tax %>%
  filter(slope.pval < 0.05 ) %>% #& intercept.pval < value_intercept_pval
  filter(slope > 0) %>%
  subset(order == "SAR11 clade")

##FIGURE 2-----------------------------
## GROWTH RATES DISTRIBUTION PATTERNS AT DIFFERENT TAXONOMIC LEVELS
##Single ASV growth rates distribution pattern for different phylums
#functions needed
source('src/growth.rates.distr.tax.ranks.ridges.divided.R')
source('src/growth.rates.distr.tax.ranks.ridges.R')

##Density plot general
reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  #subset(asv_num == 'asv54') %>%
  #mutate(slope_chosen_days_factor = as.factor(slope_chosen_days)) %>%
  #filter(pseudoabundance > 1000) %>%
  #group_by(order, slope_chosen_days_factor) %>%
  #summarise(counts = n()) %>%
  ggplot(aes(x=slope_chosen_days))+ #, colour = phylum
  #geom_density(adjust = 0.25)
  geom_freqpoly(stat = 'count')+
  stat_bin(binwidth = 0.25)+
  labs(x= expression("Growth rate (μ) day"^"-1"), y = "Density")+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0), legend.position = "right",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

reg_all_slopes_chosen_silva_tax_1perc %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  ggplot(aes(x=slope_chosen_days))+ #, colour = phylum
  #geom_density(adjust = 0.25)
  geom_freqpoly(stat = 'count')+
  stat_bin(binwidth = 0.25)+
  labs(x= expression("Growth rate (μ) day"^"-1"), y = "Density")+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0), legend.position = "right",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

##GEOM RIDGES FIGURE 2----
##fct_infreq: reorder by number of observations with each level
# reg_all_slopes_chosen_silva_tax %>%
#   mutate(rank = factor(rank, levels = c('class', 'order', 'family', 'genus')))

#perquè tots els gràfics estiguin ordenats iguals
reg_all_slopes_chosen_silva_tax <- reg_all_slopes_chosen_silva_tax  %>%
  ungroup() %>%
  mutate(phylum_f = fct_rev(fct_infreq(phylum_f)))

reg_all_slopes_chosen_silva_tax %$% phylum_f %>% levels
reg_all_slopes_chosen_silva_tax %>%
  glimpse()

counts <-   reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0 ) %>%
  group_by(phylum_f) %>%
  filter(n() >= 2) %>%
  mutate(counts = n()) %>%
  as_tibble() %>%
  group_by(phylum_f, counts) %>%
  summarize() %>%
  as_tibble()

reg_all_slopes_chosen_silva_tax %$%
  phylum_f
ridge_ph <- 
  reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0 ) %>%
  group_by(phylum_f) %>%
  filter(n() >= 2) %>%
  mutate(counts = paste('n = ', n())) %>% #
  ggplot(aes(y = fct_rev(fct_infreq(phylum_f)), x = slope_chosen_days, fill = phylum_f,  group = phylum_f, label = counts))+#fct_rev(fct_infreq(phylum_f))
  geom_density_ridges(alpha = 0.8, panel_scaling = TRUE, scale = 1,
                      jittered_points = TRUE,
                      point_shape = 21, point_size = 0.2, point_alpha = 0.0,
                      #position = 
                      quantile_lines = TRUE,
                      quantile_fun = mean)+
  scale_fill_manual(values = palette_phylums_assigned)+
  scale_x_continuous(limits = c(0,11))+
  coord_cartesian(clip = "off") +
  #scale_y_discrete(position = "right") +
  # stat_summary(fun = 'mean',
  #              geom = 'line', color = 'black')+
  #annotate(geom = 'text', label = counts)+
  geom_text(nudge_x = 8.8, nudge_y = 0.3, check_overlap = TRUE, size = 3)+
  #coord_flip()+
  labs(y = '', fill= 'Phylum', x = expression("Growth rate day"^"-1"), title = 'Phylum')+
  theme_bw()+
  theme(legend.position = "none", strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 12),
        plot.margin= unit(c(1, 1, 0.2, 0.2), "lines"), panel.border = element_blank())

# reg_all_slopes_chosen_silva_tax %>%
#   mutate(class_freq = fct_reorder(phylum_f, class_f))

ridge_cl <- reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05,
         slope_chosen_days > 0) %>% #& intercept.pval < value_intercept_pval
  group_by(phylum_f) %>%
  filter(n() >= 2) %>%
  #summarize(n = n())
  # filter(n() > 2) %>% #, .preserve = FALSE
  # ungroup() %>%
  group_by(phylum_f) %>%
  mutate(counts = paste('n = ', n())) %>% #
  ggplot(aes(y = fct_rev(fct_infreq(phylum_f)), x = slope_chosen_days, fill = phylum_f, label = counts))+
  geom_density_ridges(aes(fill = phylum_f, group = class_f), scale = 1, alpha = 0.7,
                      jittered_points = TRUE,
                      point_shape = 21, point_size = 0.2, point_alpha = 0.0,
                      #position = 
                      quantile_lines = TRUE,
                      quantile_fun = mean)+
  scale_fill_manual(values = palette_phylums_assigned)+
  ##facet_grid(vars(class))+
  labs(y = 'Class', fill= 'Phylum', x =expression("Growth rate day"^"-1"), title = 'Class')+
  scale_x_continuous(limits = c(0,11))+
  #facet_grid(rows = vars(class), scales = 'free')+
  geom_text(nudge_x = 7.8, nudge_y = 0.3, check_overlap = TRUE, size = 3)+
  #scale_y_discrete(position = "right") +
  #coord_flip()+
  theme_bw()+
  theme(legend.position = "none", strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 0),
        plot.margin= unit(c(1, 1, 0.2, 0.2), "lines"), axis.title.y = element_blank(), panel.border = element_blank())

ridge_o <- reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  group_by(phylum_f) %>%
  filter(n() >= 2) %>%
  ungroup() %>%
  group_by(phylum_f) %>%
  mutate(counts = paste('n = ', n())) %>% #
  ggplot(aes(y = fct_rev(fct_infreq(phylum_f)), x = slope_chosen_days, fill = phylum_f, label = counts))+#fct_infreq(phylum_f, order_f) per ordenar per freqüència
  geom_density_ridges(aes(fill = phylum_f, group = order_f), scale = 1.0, alpha = 0.6,
                      jittered_points = TRUE,
                      point_shape = 21, point_size = 0.2, point_alpha = 0.0,
                      #position = 
                      quantile_lines = TRUE,
                      quantile_fun = mean)+
  scale_fill_manual(values = palette_phylums_assigned)+
  labs(y = '', fill= 'Phylum', x =expression("Growth rate day"^"-1"), title = 'Order')+
  scale_x_continuous(limits = c(0,11))+
  #facet_grid(rows = vars(class), scales = 'free')+
  #scale_y_discrete(position = "right") +
  geom_text(nudge_x = 7.8, nudge_y = 0.3, check_overlap = TRUE, size = 3)+
  #coord_flip()+
  theme_bw()+
  theme(legend.position = "none", strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 0),
        plot.margin= unit(c(1, 1, 0.2, 0.2), "lines"), panel.border = element_blank())

ridge_f <- reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  group_by(phylum_f) %>%
  filter(n() >= 2) %>%
  # filter(n() > 2) %>%
  # ungroup() %>%
  group_by(phylum_f) %>%
  mutate(counts = paste('n = ', n())) %>% #
  ungroup() %>%
  ggplot(aes(y = fct_rev(fct_infreq(phylum_f)), x = slope_chosen_days, fill = phylum_f, label = counts))+
  geom_density_ridges(aes(fill = phylum_f, group = family_f), scale = 1.0, alpha = 0.6,
                      jittered_points = TRUE,
                      point_shape = 21, point_size = 0.2, point_alpha = 0.0,
                      #position = 
                      quantile_lines = TRUE,
                      quantile_fun = mean)+#fct_rev(fct_infreq(phylum_f))
  scale_fill_manual(values = palette_phylums_assigned)+
  labs(y = '', fill= 'Phylum', x =expression("Growth rate day"^"-1"), title = 'Family')+
  scale_x_continuous(limits = c(0,11))+
  geom_text(nudge_x = 7.8, nudge_y = 0.3, check_overlap = TRUE, size = 3)+
  #facet_grid(rows = vars(class), scales = 'free')+
  theme_bw()+
  theme(legend.position = "none", strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 0),
        plot.margin= unit(c(1, 1, 0.2, 0.2), "lines"), panel.border = element_blank())


##reordeno els phylums perquè al filtrar em canvia l'ordre si ho faig per freqüències
reg_all_slopes_chosen_silva_tax_asv_reorder <- reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  # group_by(asv_f) %>%
  # filter(n() > 2) %>%
  group_by(phylum_f) %>%
  mutate(counts = paste('n = ', n())) %>% #
  ungroup()

reg_all_slopes_chosen_silva_tax_asv_reorder$phylum_f <- reg_all_slopes_chosen_silva_tax_asv_reorder$phylum_f %>%
  factor( levels = c("Proteobacteria",  "Bacteroidota", "Crenarchaeota", "Planctomycetota", "Actinobacteriota", 
                     "Verrucomicrobiota", "Cyanobacteria", "Bdellovibrionota", "Firmicutes", "Nitrospinota", "Nitrospirota", 
                     "Acidobacteriota", "Campilobacterota", "Myxococcota", "Desulfobacterota",
                     "Deinococcota" , "Chloroflexi" , "Fusobacteriota", 
                     "Spirochaetota", "Abditibacteriota", "Latescibacterota",
                     "Methylomirabilota", "Halanaerobiaeota", "Sumerlaeota", "Calditrichota", "Gemmatimonadota"), 
          ordered = TRUE)
# reg_all_slopes_chosen_silva_tax_asv_reorder$phylum_f %>%
#   levels()

ridge_a <- reg_all_slopes_chosen_silva_tax_asv_reorder %>%
  group_by(phylum_f) %>%
  filter(n() >= 2) %>%
  ggplot(aes(y = fct_rev(phylum_f), x = slope_chosen_days, label = counts))+
  geom_density_ridges(aes(fill = phylum_f, group = asv_f), scale = 1.0, alpha = 0.6,
                      jittered_points = TRUE,
                      point_shape = 21, point_size = 0.2, point_alpha = 0.0,
                      #position =
                      quantile_lines = TRUE,
                      quantile_fun = mean)+
  scale_fill_manual(values = palette_phylums_assigned)+
  labs(y = '', fill= 'Phylum', x =expression("Growth rate day"^"-1"), title = 'ASV')+
  scale_x_continuous(limits = c(0,11))+
  geom_text(nudge_x = 7.2, nudge_y = 0.35, check_overlap = TRUE, size = 3)+
  #facet_grid(rows = vars(class), scales = 'free')+
  coord_cartesian(clip = "off") +
  theme_bw()+
  theme(legend.position = "none", strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 0),
        plot.margin= unit(c(1, 1, 0.2, 0.2), "lines"), panel.border = element_blank())


##Percentatge de GR significatives calculades per cada phylum ------ 
#library(scales)
#library(forcats)
reg_all_slopes_chosen_silva_tax %$%
  phylum_f %>%
  levels()
perc_gr_phylum <- reg_all_slopes_chosen_silva_tax %>%
  subset(slope_chosen_days != is.na(slope_chosen_days)) %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  group_by(phylum_f) %>%
  filter(n() > 2) %>%
  group_by(phylum_f, slope_chosen_days) %>% 
  summarize(n = n()) %>%
  group_by(phylum_f) %>%
  summarize(n = n()) %>%
  mutate(sum = sum(n),
         perc = n/sum) %>%
  ungroup() %>%
  ggplot(aes(sum, perc, fill = phylum_f))+#group = phylum_f,
  geom_col(aes(fill = fct_reorder(phylum_f, perc, .desc = TRUE)))+ #ct_rev(fct_infreq(phylum_f, ordered = NA)
  #coord_polar(theta='y')+
  #scale_x_continuous(limits = c(14584, 14584))+
  scale_y_continuous(labels = percent_format())+
  scale_fill_manual(values = palette_phylums_assigned)+
  labs(x = 'Phylum', y = 'Percentage of significant\ngrowth rates calculated', fill = 'Phylum')+
  theme_bw()+
  theme(axis.text.x = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(), aspect.ratio = 10/2, panel.border = element_blank(),
        axis.ticks = element_blank(), legend.position = 'none', axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14))


# reg_all_slopes_chosen_silva_tax %>%
#   subset(slope_chosen_days != is.na(slope_chosen_days)) %>%
#   group_by(phylum_f, slope_chosen_days) %>% 
#   summarize(n = n()) %>%
#   group_by(phylum_f) %>%
#   summarize(n = n()) %>%
#   mutate(sum = sum(n),
#          perc = n/sum) %>%
#   ggplot(aes(fct_infreq(phylum_f, ordered = NA), perc, color = phylum_f))+ #group = phylum_f,
#   geom_point(aes(color = phylum_f))+
#   scale_y_continuous(labels = percent_format())+
#   scale_color_manual(values = new_palette_large)+
#   labs(x = 'Phylum', y = 'Percentage of growth\nrates calculated', fill = 'Phylum')+
#   theme_bw()+
#   theme(axis.text.x = element_blank(), panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(), aspect.ratio = 10/3,
#         axis.ticks = element_blank())

##Percentatge de GR significatives representades en cada nivell taxonòmic (+ de 2 gr per grup per poder fer distribució)----
#calculo el total de gr sense  filtrar per grup més de 2 taxes de creixement per fer el %
reg_all_slopes_chosen_silva_tax %>%
  subset(slope_chosen_days != is.na(slope_chosen_days)) %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  # group_by(asv_f) %>%
  # filter(n() > 2) %>%
  summarize(n = n()) %>%
  summarize(sum = sum(n))

gr_asv_perc <- reg_all_slopes_chosen_silva_tax %>%
  subset(slope_chosen_days != is.na(slope_chosen_days)) %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  group_by(asv_f) %>%
  filter(n() > 2) %>%
  summarize(n = n()) %>%
  mutate(sum = sum(n),
         asv_gr_perc = sum/5543) %>% ###4019 és el número total de gr que he calculat
  select(asv_gr_perc) %>%
  unique() %>%
  as_tibble()

gr_f_perc <- reg_all_slopes_chosen_silva_tax %>%
  subset(slope_chosen_days != is.na(slope_chosen_days)) %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  group_by(family) %>%
  filter(n() > 2) %>%
  summarize(n = n()) %>%
  mutate(sum = sum(n),
         f_gr_perc = sum/5543) %>%
  select(f_gr_perc) %>%
  unique() %>%
  as_tibble()

gr_o_perc <- reg_all_slopes_chosen_silva_tax %>%
  subset(slope_chosen_days != is.na(slope_chosen_days)) %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  group_by(order) %>%
  filter(n() > 2) %>%
  summarize(n = n()) %>%
  mutate(sum = sum(n),
         order_gr_perc = sum/5543) %>%
  select(order_gr_perc) %>%
  unique() %>%
  as_tibble()

gr_c_perc <- reg_all_slopes_chosen_silva_tax %>%
  subset(slope_chosen_days != is.na(slope_chosen_days)) %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  group_by(class) %>%
  filter(n() > 2) %>%
  summarize(n = n()) %>%
  mutate(sum = sum(n),
         class_gr_perc = sum/5543) %>%
  select(class_gr_perc) %>%
  unique() %>%
  as_tibble()

gr_p_perc <- reg_all_slopes_chosen_silva_tax %>%
  subset(slope_chosen_days != is.na(slope_chosen_days)) %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  group_by(phylum) %>%
  filter(n() > 2) %>%
  summarize(n = n()) %>%
  mutate(sum = sum(n),
         phylum_gr_perc = sum/5543) %>%
  select(phylum_gr_perc) %>%
  unique() %>%
  as_tibble() 


labels_perc_gr_rank <-  as_labeller(c(phylum_gr_perc = 'Phylum',
                                      class_gr_perc = 'Class',
                                      order_gr_perc = 'Order',
                                      f_gr_perc = 'Family',
                                      asv_gr_perc = 'ASV'))

gr_tax_ranks_perc <- cbind(gr_asv_perc, gr_f_perc, gr_o_perc, gr_c_perc, gr_p_perc) %>%
  as_tibble() %>%
  pivot_longer(cols = 1:5) %>%
  ggplot(aes(x = as_factor(name), y = value, group = 1))+
  geom_point()+
  coord_flip()+
  geom_line()+
  labs(y='Percentage', x = 'Growth rates distribution plotted\nat different taxonomic ranks')+
  scale_x_discrete(labels = labels_perc_gr_rank)+
  scale_y_continuous(labels = percent_format())+
  theme_bw()+
  theme(aspect.ratio = 10/3, panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = 'none', axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14))


##Mean GR at different taxonomic levels
labels_mean_rank <- as_labeller(c( mean_phylum = 'Phylum',
                                   mean_class = 'Class',
                                   mean_order = 'Order',
                                   mean_family = 'Family',
                                   mean_asv = 'ASV'))

mean_gr_ranks <- reg_all_slopes_chosen_silva_tax %>%
  subset(slope_chosen_days != is.na(slope_chosen_days)) %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  group_by(phylum_f) %>%
  filter(n() > 2) %>%
  mutate(mean_phylum = mean(slope_chosen_days)) %>%
  ungroup() %>%
  group_by(class) %>%
  mutate(mean_class = mean(slope_chosen_days)) %>%
  ungroup() %>%
  group_by(order) %>%
  mutate(mean_order = mean(slope_chosen_days)) %>%
  ungroup() %>%
  group_by(family) %>%
  mutate(mean_family = mean(slope_chosen_days)) %>%
  ungroup() %>%
  group_by(asv_num) %>%
  mutate(mean_asv = mean(slope_chosen_days)) %>%
  select(starts_with('mean'), 'phylum_f', 'class', 'order', 'family', 'asv_num') %>%
  pivot_longer(cols = starts_with('mean')) %>%
  group_by(phylum_f) %>%
  distinct(value, name, phylum_f) %>%
  ggplot(aes(fct_rev(as_factor(name)), value))+
  geom_point(aes(color = phylum_f), alpha = 0.9, position = position_jitter(0.3))+
  #geom_boxplot(alpha= 0.05)+
  geom_violin(alpha = 0.2, draw_quantiles = c(0.25, 0.75))+
  scale_color_manual(values = palette_phylums)+
  labs(y = 'Mean growth rate at different\ntaxonomic ranks', x= 'Taxonomic rank', color = 'Phylum')+
  scale_x_discrete(labels = labels_mean_rank)+
  coord_flip()+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = 'none', axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14))


##figure 2 creating pannel design------
distribution_gr_plots_ridge <- multi_panel_figure(columns = 7, rows = 3, width = 410, height = 297, 
                                                  column_spacing = 0.1, unit = 'mm',
                                                  panel_label_type = 'none')

distribution_gr_plots_ridge  %<>%
  fill_panel(perc_gr_phylum, column = 1, row = 1) %<>%
  fill_panel(gr_tax_ranks_perc, column = 1, row = 2) %<>%
  fill_panel(mean_gr_ranks, column = 1, row = 3) %<>%
  fill_panel(ridge_ph, column = 2:3, row = 1:3) %<>%
  fill_panel(ridge_cl, column = 4, row = 1:3) %<>%
  fill_panel(ridge_o, column = 5, row = 1:3) %<>%
  fill_panel(ridge_f, column = 6, row = 1:3) %<>%
  fill_panel(ridge_a, column = 7, row = 1:3)


##CAMBIAR EL NOM ABANS DE GUARDAR!!!!
ggsave('distribution_gr_plots_ridge_order_freq_perc1_remiau.pdf', distribution_gr_plots_ridge, 
       path = "~/Documentos/Doctorat/REMEI/results/figures/raw/",
       width = 470,
       height = 397,
       units = 'mm')

##Geom_ridges function
## Proteobacteria    Crenarchaeota     Bacteroidota      Myxococcota       Planctomycetota   Verrucomicrobiota Nitrospirota     
## Actinobacteriota  Acidobacteriota   Nitrospinota      Bdellovibrionota  Deinococcota      Firmicutes        Cyanobacteria   
## Chloroflexi       Desulfobacterota 

####Differences in distribution patterns between CL and VL treatments--------------
#perquè tots els gràfics estiguin ordenats iguals
reg_all_slopes_chosen_silva_tax <- reg_all_slopes_chosen_silva_tax  %>%
  ungroup() %>%
  mutate(phylum_f = fct_rev(fct_infreq(phylum_f)))

reg_all_slopes_chosen_silva_tax %$% phylum_f %>% levels
reg_all_slopes_chosen_silva_tax %>%
  glimpse()

counts <-   reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0 ) %>%
  group_by(phylum_f) %>%
  filter(n() >= 2) %>%
  mutate(counts = n()) %>%
  as_tibble() %>%
  group_by(phylum_f, counts) %>%
  summarize() %>%
  as_tibble()

reg_all_slopes_chosen_silva_tax %$%
  phylum_f
ridge_ph <- 
  reg_all_slopes_chosen_silva_tax %>%
  subset(treatment== c('CL', 'VL')) %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0 ) %>%
  group_by(phylum_f) %>%
  filter(n() >= 2) %>%
  mutate(counts = paste('n = ', n())) %>% #
  ggplot(aes(y = fct_rev(fct_infreq(phylum_f)), x = slope_chosen_days, fill = phylum_f,  group = phylum_f, label = counts))+#fct_rev(fct_infreq(phylum_f))
  geom_density_ridges(alpha = 0.8, panel_scaling = TRUE, scale = 1,
                      jittered_points = TRUE,
                      point_shape = 21, point_size = 0.2, point_alpha = 0.0,
                      #position = 
                      quantile_lines = TRUE,
                      quantile_fun = mean)+
  facet_grid(vars(treatment))+
  scale_fill_manual(values = palette_phylums)+
  scale_x_continuous(limits = c(0,11))+
  coord_cartesian(clip = "off") +
  #scale_y_discrete(position = "right") +
  # stat_summary(fun = 'mean',
  #              geom = 'line', color = 'black')+
  #annotate(geom = 'text', label = counts)+
  geom_text(nudge_x = 8.8, nudge_y = 0.3, check_overlap = TRUE, size = 3)+
  #coord_flip()+
  labs(y = '', fill= 'Phylum', x = expression("Growth rate day"^"-1"), title = 'Phylum')+
  theme_bw()+
  theme(legend.position = "none", strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 12),
        plot.margin= unit(c(1, 1, 0.2, 0.2), "lines"), panel.border = element_blank())

# reg_all_slopes_chosen_silva_tax %>%
#   mutate(class_freq = fct_reorder(phylum_f, class_f))

ridge_cl <- reg_all_slopes_chosen_silva_tax %>%
  subset(treatment== c('CL', 'VL')) %>%
  filter(pvalue_slope_chosen < 0.05,
         slope_chosen_days > 0) %>% #& intercept.pval < value_intercept_pval
  group_by(phylum_f) %>%
  filter(n() >= 2) %>%
  #summarize(n = n())
  # filter(n() > 2) %>% #, .preserve = FALSE
  # ungroup() %>%
  group_by(phylum_f) %>%
  mutate(counts = paste('n = ', n())) %>% #
  ggplot(aes(y = fct_rev(fct_infreq(phylum_f)), x = slope_chosen_days, fill = phylum_f, label = counts))+
  geom_density_ridges(aes(fill = phylum_f, group = class_f), scale = 1, alpha = 0.7,
                      jittered_points = TRUE,
                      point_shape = 21, point_size = 0.2, point_alpha = 0.0,
                      #position = 
                      quantile_lines = TRUE,
                      quantile_fun = mean)+
  facet_grid(vars(treatment))+
  scale_fill_manual(values = palette_phylums)+
  ##facet_grid(vars(class))+
  labs(y = 'Class', fill= 'Phylum', x =expression("Growth rate day"^"-1"), title = 'Class')+
  scale_x_continuous(limits = c(0,11))+
  #facet_grid(rows = vars(class), scales = 'free')+
  geom_text(nudge_x = 7.8, nudge_y = 0.3, check_overlap = TRUE, size = 3)+
  #scale_y_discrete(position = "right") +
  #coord_flip()+
  theme_bw()+
  theme(legend.position = "none", strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 0),
        plot.margin= unit(c(1, 1, 0.2, 0.2), "lines"), axis.title.y = element_blank(), panel.border = element_blank())

ridge_o <- reg_all_slopes_chosen_silva_tax %>%
  subset(treatment== c('CL', 'VL')) %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  group_by(phylum_f) %>%
  filter(n() >= 2) %>%
  ungroup() %>%
  group_by(phylum_f) %>%
  mutate(counts = paste('n = ', n())) %>% #
  ggplot(aes(y = fct_rev(fct_infreq(phylum_f)), x = slope_chosen_days, fill = phylum_f, label = counts))+#fct_infreq(phylum_f, order_f) per ordenar per freqüència
  geom_density_ridges(aes(fill = phylum_f, group = order_f), scale = 1.0, alpha = 0.6,
                      jittered_points = TRUE,
                      point_shape = 21, point_size = 0.2, point_alpha = 0.0,
                      #position = 
                      quantile_lines = TRUE,
                      quantile_fun = mean)+
  facet_grid(vars(treatment))+
  scale_fill_manual(values = palette_phylums)+
  labs(y = '', fill= 'Phylum', x =expression("Growth rate day"^"-1"), title = 'Order')+
  scale_x_continuous(limits = c(0,11))+
  #facet_grid(rows = vars(class), scales = 'free')+
  #scale_y_discrete(position = "right") +
  geom_text(nudge_x = 7.8, nudge_y = 0.3, check_overlap = TRUE, size = 3)+
  #coord_flip()+
  theme_bw()+
  theme(legend.position = "none", strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 0),
        plot.margin= unit(c(1, 1, 0.2, 0.2), "lines"), panel.border = element_blank())

ridge_f <- reg_all_slopes_chosen_silva_tax %>%
  subset(treatment== c('CL', 'VL')) %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  group_by(phylum_f) %>%
  filter(n() >= 2) %>%
  # filter(n() > 2) %>%
  # ungroup() %>%
  group_by(phylum_f) %>%
  mutate(counts = paste('n = ', n())) %>% #
  ungroup() %>%
  ggplot(aes(y = fct_rev(fct_infreq(phylum_f)), x = slope_chosen_days, fill = phylum_f, label = counts))+
  geom_density_ridges(aes(fill = phylum_f, group = family_f), scale = 1.0, alpha = 0.6,
                      jittered_points = TRUE,
                      point_shape = 21, point_size = 0.2, point_alpha = 0.0,
                      #position = 
                      quantile_lines = TRUE,
                      quantile_fun = mean)+#fct_rev(fct_infreq(phylum_f))
  facet_grid(vars(treatment))+
  scale_fill_manual(values = palette_phylums)+
  labs(y = '', fill= 'Phylum', x =expression("Growth rate day"^"-1"), title = 'Family')+
  scale_x_continuous(limits = c(0,11))+
  geom_text(nudge_x = 7.8, nudge_y = 0.3, check_overlap = TRUE, size = 3)+
  #facet_grid(rows = vars(class), scales = 'free')+
  theme_bw()+
  theme(legend.position = "none", strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 0),
        plot.margin= unit(c(1, 1, 0.2, 0.2), "lines"), panel.border = element_blank())


##reordeno els phylums perquè al filtrar em canvia l'ordre si ho faig per freqüències
reg_all_slopes_chosen_silva_tax_asv_reorder <- reg_all_slopes_chosen_silva_tax %>%
  subset(treatment== c('CL', 'VL')) %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  # group_by(asv_f) %>%
  # filter(n() > 2) %>%
  group_by(phylum_f) %>%
  mutate(counts = paste('n = ', n())) %>% #
  ungroup()

reg_all_slopes_chosen_silva_tax_asv_reorder$phylum_f <- reg_all_slopes_chosen_silva_tax_asv_reorder$phylum_f %>%
  factor( levels = c("Proteobacteria",  "Bacteroidota", "Crenarchaeota", "Planctomycetota", "Actinobacteriota", 
                     "Verrucomicrobiota", "Cyanobacteria", "Bdellovibrionota", "Firmicutes", "Nitrospinota", "Nitrospirota", 
                     "Acidobacteriota", "Campilobacterota", "Myxococcota", "Desulfobacterota",
                     "Deinococcota" , "Chloroflexi" , "Fusobacteriota", 
                     "Spirochaetota", "Abditibacteriota", "Latescibacterota",
                     "Methylomirabilota", "Halanaerobiaeota", "Sumerlaeota", "Calditrichota", "Gemmatimonadota"), 
          ordered = TRUE)
# reg_all_slopes_chosen_silva_tax_asv_reorder$phylum_f %>%
#   levels()

ridge_a <- reg_all_slopes_chosen_silva_tax_asv_reorder %>%
  subset(treatment== c('CL', 'VL')) %>%
  group_by(phylum_f) %>%
  filter(n() >= 2) %>%
  ggplot(aes(y = fct_rev(phylum_f), x = slope_chosen_days, label = counts))+
  geom_density_ridges(aes(fill = phylum_f, group = asv_f), scale = 1.0, alpha = 0.6,
                      jittered_points = TRUE,
                      point_shape = 21, point_size = 0.2, point_alpha = 0.0,
                      #position =
                      quantile_lines = TRUE,
                      quantile_fun = mean)+
  facet_grid(vars(treatment))+
  scale_fill_manual(values = palette_phylums)+
  labs(y = '', fill= 'Phylum', x =expression("Growth rate day"^"-1"), title = 'ASV')+
  scale_x_continuous(limits = c(0,11))+
  geom_text(nudge_x = 7.2, nudge_y = 0.35, check_overlap = TRUE, size = 3)+
  #facet_grid(rows = vars(class), scales = 'free')+
  coord_cartesian(clip = "off") +
  theme_bw()+
  theme(legend.position = "none", strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 0),
        plot.margin= unit(c(1, 1, 0.2, 0.2), "lines"), panel.border = element_blank())


##figure 2 CL vs VL distributions creating pannel design------
distribution_gr_plots_ridge_CV <- multi_panel_figure(columns = 6, rows = 1, width = 410, height = 297, 
                                                  column_spacing = 0.1, unit = 'mm',
                                                  panel_label_type = 'none')

distribution_gr_plots_ridge_CV  %<>%
  fill_panel(ridge_ph, column = 1:2, row = 1) %<>%
  fill_panel(ridge_cl, column = 3, row = 1) %<>%
  fill_panel(ridge_o, column = 4, row = 1) %<>%
  fill_panel(ridge_f, column = 5, row = 1) %<>%
  fill_panel(ridge_a, column = 6, row = 1)

##CAMBIAR EL NOM ABANS DE GUARDAR!!!!
ggsave('distribution_gr_plots_ridge_order_freq_perc1_remiau_CV.pdf', distribution_gr_plots_ridge_CV, 
       path = "~/Documentos/Doctorat/REMEI/results/figures/raw/",
       width = 470,
       height = 397,
       units = 'mm')



##phylums gr distribution separated ordered by frequency
ridge_proteo <- 
  growthrates.distribution.tax.rank.ridges.divided(
    data = reg_all_slopes_chosen_silva_tax, 
    phylum_to_explore = 'Proteobacteria',
    title = 'Proteobacteria class',
    axis_y_title = '',
    x_c = 4,
    x_o = 10,
    x_f = 10,
    x_a = 8) %>%
  as_ggplot() ##estic afegint fct_rev perquè m'ordeni els gràfics de manera lògica alpha i gamma

reg_all_slopes_chosen_silva_tax %>%
  colnames()

##canvis en les distribucions entre CL i VL-------
asv_filtered <- reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  # group_by(asv_f) %>%
  # filter(n() > 2) %>%
  group_by(phylum_f) %>%
  mutate(counts = paste('n = ', n())) %>% #
  ungroup()

asv_filtered %>%
  subset(phylum == 'Proteobacteria') %>%
  subset(treatment %in% c('CL', 'VL')) %>%
  ggplot(aes(y = fct_rev(order_f), x = slope_chosen_days, label = counts))+
  geom_density_ridges(aes(fill = phylum_f, group = order_f), scale = 1.0, alpha = 0.6,
                      jittered_points = TRUE,
                      point_shape = 21, point_size = 0.2, point_alpha = 0.0,
                      #position =
                      quantile_lines = TRUE,
                      quantile_fun = mean)+
  scale_fill_manual(values = palette_phylums)+
  labs(y = '', fill= 'Phylum', x =expression("Growth rate day"^"-1"), title = 'ASV')+
  scale_x_continuous(limits = c(0,11))+
  geom_text(nudge_x = 7.2, nudge_y = 0.35, check_overlap = TRUE, size = 3)+
  facet_grid(cols = vars(treatment), scales = 'free')+
  coord_cartesian(clip = "off") +
  theme_bw()+
  theme(legend.position = "none", strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 0),
        plot.margin= unit(c(1, 1, 0.2, 0.2), "lines"), panel.border = element_blank())


##falten números a asv level
test <- reg_all_slopes_chosen_silva_tax %>%
  subset(phylum == 'Proteobacteria') %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  group_by(asv_num) %>%
  filter(n() > 8) %>%
  mutate(counts = paste('n = ', n())) %>%
  group_by(family) %>%
  filter(n() > 10) %>%
  ungroup() %>%
  ggplot(aes(y = asv_f, x = slope_chosen_days, fill = class_f, label = counts))+
  # geom_histogram(aes(x=slope_chosen_days, group = asv_num), colour  = "#4F507F", fill = "#4F507F", 
  #                #adjust= 0.5, #bins = bin_width_a,
  #                binwidth = bin_width_a,
  #                show.legend = FALSE, alpha = 0.2)+
  geom_density_ridges(aes(fill = class_f), scale = 4, alpha = 0.6,
                      quantile_lines = TRUE,
                      quantile_fun = mean)+
  labs(x= expression("Growth rate day"^"-1"), y = "", title = 'ASV')+
  scale_x_continuous(limits = c(0,11))+
  scale_fill_manual(values = c("#c59a3e","#fcca46", "#806f52"))+
  facet_grid(rows = vars(family_f), scales = 'free', space = 'free', switch = 'y')+
  #scale_y_continuous(limits = c(y1_a,y2_a))+
  #stat_bin(binwidth = bin_width_a)+
  geom_text(nudge_x = 8.8, nudge_y = 0.3,  size = 3)+ #check_overlap = TRUE,
  theme_ridges(center_axis_labels = TRUE)+ 
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), #axis.text.y = element_text(size = 0), 
        strip.background = element_blank(), axis.ticks.y = element_blank(),
        strip.text.y = element_text(size = 10),
        plot.margin= unit(c(1.2, 1, 1.2, 1.2), "lines"))

ggsave('remiau_ridge_proteobacteria_filtered.pdf', ridge_proteo, 
       path = "~/Documentos/Doctorat/REMEI/results/figures/raw/",
       width = 470,
       height = 397,
       units = 'mm')

reg_all_slopes_chosen_silva_tax %$%
  family_f %>%
  levels()

##proves a nivell de classe per endreçar millor els facet grids

reg_all_slopes_chosen_silva_tax %>%
  subset(phylum == 'Proteobacteria') %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  group_by(family) %>%
  filter(n() > 10) %>%
  mutate(counts = paste('n = ', n())) %>%
  group_by(order) %>%
  filter(n() > 10) %>%
  ungroup() %>%
  ggplot(aes(y = family, x = slope_chosen_days, fill = phylum_f, label = counts))+ #interaction( family_f, order_f, sep = '\n')
  # geom_histogram(aes(x=slope_chosen_days, group = family), colour = "#FFDF26", fill = "#FFDF26", #adjust= 0.5, 
  #                # bins = bin_width_f,
  #                binwidth = bin_width_f,
  #                show.legend = FALSE, alpha = 0.2)+
  geom_density_ridges(aes(fill = phylum_f), scale = 3, alpha = 0.6,
                      quantile_lines = TRUE,
                      quantile_fun = mean)+
  labs(x= "")+
  #stat_bin(binwidth = bin_width_f)+
  scale_x_continuous(limits = c(0,11))+
  scale_fill_manual(values = palette_phylums)+
  facet_grid(rows = vars(order_f), scales = 'free', space = 'fixed')+
  geom_text(nudge_x = 8.8, nudge_y = 0.3, check_overlap = TRUE, size = 3)+
  # scale_y_continuous(limits = c(y1_f,y2_f))+
  theme_bw()+ 
  theme(strip.text.x = element_text(size = 0), axis.text.x = element_text(angle = 0), 
        axis.ticks = element_blank(), legend.position = "none",axis.title.y = element_text(size = 0),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_blank(),
        strip.text.y = element_blank(),
        plot.margin= unit(c(1, 0.2, 0.2, 0.2), "lines"))

##BACTEROIDOTA
ridge_bacterio <- growthrates.distribution.tax.rank.ridges.divided(
  data = reg_all_slopes_chosen_silva_tax, 
  phylum_to_explore = 'Bacteroidota',
  title = 'Bacteroidota \n
  Class',
  axis_y_title = '',
  x_c = 4,
  x_o = 10,
  x_f = 10,
  x_a = 8) %>%
  as_ggplot()

ggsave('remiau_ridge_bacteroidota_filtered_v2.pdf', ridge_bacterio, 
       path = "~/Documentos/Doctorat/REMEI/results/figures/raw/",
       width = 470,
       height = 397,
       units = 'mm')


ridge_plancto <- growthrates.distribution.tax.rank.ridges.divided(
  data = reg_all_slopes_chosen_silva_tax, 
  phylum_to_explore = 'Planctomycetota',
  title = 'Planctomycetota',
  axis_y_title = '') %>%
  as_ggplot()

ridge_actino <- growthrates.distribution.tax.rank.ridges(
  data = reg_all_slopes_chosen_silva_tax, 
  phylum_to_explore = 'Actinobacteriota',
  title = 'Actinobacteriota',
  axis_y_title = '') %>%
  as_ggplot()

ridge_crena <- growthrates.distribution.tax.rank.ridges(
  data = reg_all_slopes_chosen_silva_tax, 
  phylum_to_explore = 'Crenarchaeota',
  title = 'Crenarchaeota',
  axis_y_title = '') %>%
  as_ggplot()

ridge_verru <- growthrates.distribution.tax.rank.ridges(
  data = reg_all_slopes_chosen_silva_tax, 
  phylum_to_explore = 'Verrucomicrobiota',
  title = 'Verrucomicrobiota',
  axis_y_title = '') %>%
  as_ggplot()

ridge_firmi <- growthrates.distribution.tax.rank.ridges(
  data = reg_all_slopes_chosen_silva_tax, 
  phylum_to_explore = 'Firmicutes',
  title = 'Firmicutes',
  axis_y_title = '') %>%
  as_ggplot()


reg_all_slopes_chosen_silva_tax_filt %$%
  phylum_f %>%
  unique()

distribution_gr_plots_ridge_separated <- multi_panel_figure(columns = 4, rows = 8, width = 410, height = 797, 
                                                            column_spacing = 0.1, unit = 'mm',
                                                            row_spacing = 0.1, 
                                                            panel_label_type = 'none' 
)

distribution_gr_plots_ridge_separated  %<>%
  #fill_panel(perc_gr_phylum, column = 1, row = 1:2) %<>%
  fill_panel(ridge_proteo, column = 1:4, row = 1:3)   %<>%
  fill_panel(ridge_bacterio, column = 1:4, row = 4:5) %<>%
  fill_panel (ridge_plancto, column = 1:4, row = 6)  %<>% 
  fill_panel (ridge_actino, column = 1:4, row = 7)  %<>% 
  fill_panel(ridge_crena, column = 1:4, row = 8)  #%<>%


#fill_panel (ridge_firmi, column = 1:4, row = 8)

ggsave('distribution_gr_plots_ridge_separated_5most_freq.pdf', distribution_gr_plots_ridge_separated, 
       path = "~/Documentos/Doctorat/REMEI/results/figures/raw/",
       width = 410,
       height = 797,
       units = 'mm')

#plot for phylums with most of the GR calculated Bacteroidota & Proteobacteria--------
distribution_gr_plots_ridge_top_phylums <- multi_panel_figure(columns = 1, rows = 3, width = 410, height = 597, 
                                                              column_spacing = 0.1, unit = 'mm',
                                                              row_spacing = 0.1, 
                                                              panel_label_type = 'none')

distribution_gr_plots_ridge_top_phylums  %<>%
  #fill_panel(perc_gr_phylum, column = 1, row = 1:2) %<>%
  fill_panel(ridge_proteo, column = 1, row = 1:2)   %<>%
  fill_panel(ridge_bacterio, column = 1, row = 3)   #%<>%


#fill_panel (ridge_firmi, column = 1:4, row = 8)

ggsave('distribution_gr_plots_ridge_top_phylums.pdf', distribution_gr_plots_ridge_top_phylums, 
       path = "~/Documentos/Doctorat/REMEI/results/figures/raw/",
       width = 410,
       height = 597,
       units = 'mm')

####-----
reg_all_slopes_chosen_silva_tax %>%
  subset(phylum == 'Proteobacteria') %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  group_by(family_f) %>%
  filter(n() > 2) %>%
  ggplot(aes(y = family_f, x = slope_chosen_days, fill = phylum_f))+
  geom_density_ridges(aes(fill = phylum_f), scale = 3, alpha = 0.8)+
  scale_fill_manual(values = palette_phylums)+
  facet_grid(rows = vars(class), scales = 'free')+
  labs(y = 'Family', fill= 'Phylum', x =expression("Growth rate day"^"-1"))+
  scale_x_continuous(limits = c(0,11))+
  #scale_y_discrete(position = "right") +
  #coord_flip()+
  theme_bw()+
  theme(legend.position = "none", strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  
        plot.margin= unit(c(0.2, 1, 0.2, 1), "lines"))

##All phylums separated without asv level

##ASV level graph
reg_all_slopes_chosen_silva_tax %>%
  filter(pvalue_slope_chosen < 0.05) %>%
  filter(slope_chosen_days > 0) %>%
  group_by(asv_f) %>%
  filter(n() > 2) %>%
  ggplot(aes(y = asv_f, x = slope_chosen_days, fill = phylum_f))+
  geom_density_ridges(aes(fill = phylum_f), scale = 2, alpha = 0.8)+
  scale_fill_manual(values = palette_phylums)+
  labs(y = 'ASV', fill= 'Phylum', x =expression("Growth rate day"^"-1"))+
  theme_light()+
  theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0), legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 0), 
        plot.margin= unit(c(0.2, 0.2, 0.2, 0.2), "lines"))

##Què passa amb les alphaproteobacteria a asv level?
reg_all_slopes_chosen_silva_tax %>%
  subset(class == 'Alphaproteobacteria',
         slope_chosen_days != is.na(slope_chosen_days))%$%
  family %>%
  unique()

reg_all_slopes_chosen_silva_tax %>%
  subset(class == 'Alphaproteobacteria') %>%
  group_by(order_f) %>%
  ggplot(aes(x=slope_chosen_days, group = order_f, color = order_f))+ #, colour = class
  # geom_histogram(aes(x=slope_chosen_days, group = asv_num), colour  = "#4F507F", fill = "#4F507F", 
  #                #bins = bin_width_a,
  #                binwidth = 0.1,
  #                #adjust= 0.5, 
  #                show.legend = FALSE, alpha = 0.2)+
  #stat_bin(binwidth = bin_width_a)+
  #geom_freqpoly(aes(group = order_f), stat = 'bin', binwidth = 0.25)+
  geom_area(stat = 'bin', binwidth = 0.25, aes(fill = order_f), alpha = 0.6)+
  scale_x_continuous(limits = c(0,11))+
  scale_color_manual(values = palf_large(18))+
  scale_fill_manual(values = palf_large(18))+
  #scale_y_continuous(limits = c(y1,y2))+
  labs(x= expression("Growth rate day"^"-1"), y = "")+
  theme_bw()+ 
  theme(legend.position = "right", strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  
        plot.margin= unit(c(0.2, 1, 0.2, 1), "lines"))

reg_all_slopes_chosen_silva_tax %>%
  subset(class == 'Alphaproteobacteria') %>%
  group_by(family) %>%
  ggplot(aes(x=slope_chosen_days, group = family_f, color = family_f))+ #, colour = class
  # geom_histogram(aes(x=slope_chosen_days, group = asv_num), colour  = "#4F507F", fill = "#4F507F", 
  #                #bins = bin_width_a,
  #                binwidth = 0.1,
  #                #adjust= 0.5, 
  #                show.legend = FALSE, alpha = 0.2)+
  #stat_bin(binwidth = bin_width_a)+
  #geom_freqpoly(aes(group = family_f), stat = 'bin', binwidth = 0.25)+
  geom_area(stat = 'bin', binwidth = 0.25, aes(fill = family_f), alpha = 0.6)+
  scale_x_continuous(limits = c(0,11))+
  scale_color_manual(values = palf_large(48))+
  scale_fill_manual(values = palf_large(48))+
  #scale_y_continuous(limits = c(y1,y2))+
  labs(x= expression("Growth rate day"^"-1"), y = "")+
  theme_bw()+ 
  theme(legend.position = "right", strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  
        plot.margin= unit(c(0.2, 1, 0.2, 1), "lines"))

test <- reg_all_slopes_chosen_silva_tax %>%
  subset(slope_chosen_days != is.na(slope_chosen_days)) %>%
  group_by(phylum) %>%
  add_count(phylum_f, wt = slope_chosen_days)

#prova amb sankey diagram-------
library(viridis)
library(patchwork)
library(hrbrthemes)
library(circlize)
reg_all_slopes_chosen_silva_tax %>%
  subset(slope_chosen_days != is.na(slope_chosen_days)) %>%
  filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
  filter(slope_chosen_days > 0) %>%
  group_by(phylum_f) %>%
  filter(n() > 2) %>%
  mutate(mean_phylum = mean(slope_chosen_days)) %>%
  ungroup() %>%
  group_by(class) %>%
  mutate(mean_class = mean(slope_chosen_days)) %>%
  ungroup() %>%
  group_by(order) %>%
  mutate(mean_order = mean(slope_chosen_days)) %>%
  ungroup() %>%
  group_by(family) %>%
  mutate(mean_family = mean(slope_chosen_days)) %>%
  ungroup() %>%
  group_by(asv_num) %>%
  mutate(mean_asv = mean(slope_chosen_days)) %>%
  select(starts_with('mean'), 'phylum_f', 'class', 'order', 'family', 'asv_num') %>%
  pivot_longer(cols = starts_with('mean')) %>%
  group_by(phylum_f) %>%
  distinct(value, name, phylum_f) %>%
  ggplot(aes(fct_rev(as_factor(name)), value))+
  geom_point(aes(color = phylum_f), alpha = 0.9, position = position_jitter(0.3))+
  #geom_boxplot(alpha= 0.05)+
  geom_violin(alpha = 0.2, draw_quantiles = c(0.25, 0.75))+
  scale_color_manual(values = palette_phylums)+
  labs(y = 'Mean growth rate at different\ntaxonomic ranks', x= 'Taxonomic rank', color = 'Phylum')+
  scale_x_discrete(labels = labels_mean_rank)+
  coord_flip()+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = 'none', axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14))

#Growing exponentially, MAXIMAL GROWTH RATES CALCULATION------
#T3 i T4 comparison slopes
#Lag at the beginning?

## Example of GR calculation for some ASV for the supplementary material-----
reg_all_slopes_chosen_silva_tax %>%
  subset(select)

pseudoabund_df_wide_fc_silva %>%
  colnames()

#grup asv amb pvalue més significatiu: asv2, asv23, asv3, asv117, asv4 i asv10
pseudoabund_df_wide_fc_silva$season <-  pseudoabund_df_wide_fc_silva$season %>% factor( 
  levels = (c('Winter', 'Spring', 'Summer', 'Fall', 'Early_fall')))
pseudoabund_df_wide_fc_silva$treatment <- pseudoabund_df_wide_fc_silva$treatment %>%
  factor(levels=(c("CL", "CD", "PL", "PD", "DL", "VL")))
##4 temps
#asv2 asv2_Bacteroidia_Flavobacteriaceae_Aurantivirga_NA
library(ggpubr)
pseudoabund_df_wide_fc_silva %>%
  colnames()

##trobar un exemple bonic d'asv
reg_all_slopes_chosen_silva_tax_filt2 %>%
  group_by(asv_num) %>%
  summarize(n = n()) %>%
  filter(n >= 15)


reg_all_slopes_chosen_silva_tax_filt2 %>%
  group_by(asv_num, treatment) %>%
  summarize(n = n()) %>%
  filter(n >= 4)

reg_all_slopes_chosen_silva_tax_filt2 <- reg_all_slopes_chosen_silva_tax %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.01) %>%
  group_by(phylum_f) %>%
  #summarize(n = n())
  filter(n() > 2)

reg_all_slopes_chosen_silva_tax_filt %>%
  dim()

t4_exemple <- 
  pseudoabund_df_wide_fc_silva %>%
  #select(1:5, 'asv168', 'asv117',  'asv12', 'asv4') %>% #, 'asv152', 'asv3', 'asv64'
  #select(1:5, 'asv11', 'asv12', 'asv13', 'asv14', 'asv15', 'asv16', 'asv17', 'asv18', 'asv19', 'asv20') %>%
  #select(1:5, 'asv81', 'asv53', 'asv10', 'asv654', 'asv4', 'asv114') %>% # 'asv1',
  #select(1:5, 'asv2', 'asv23', 'asv3', 'asv117', 'asv4', 'asv10') %>%
  select(1:5, 'asv19') %>%
  subset(season != "Early_fall" &
            treatment == 'PL') %>%
  #subset(season != "Early_fall") %>%
  pivot_longer(cols = starts_with('asv'), names_to = 'asv_num', values_to = 'pseudoabund') %>%
  left_join(reg_all_slopes_chosen_silva_tax, by = 'asv_num') %>%
  ggplot(aes(hours.x, log(pseudoabund), color = season.x, shape = treatment.x))+
  geom_point(shape = 7, size = 3)+
  stat_poly_line(aes(group = season.x), mf.values = TRUE)+
  stat_cor(aes(label = paste(..p.label..)), label.x = 0, p.digits = 2, r.digits = 2)+
  #stat_pvalue_manual(aes(group = season.x), label = 'p.signif')+
  #stat_poly_eq(aes(group = season.x), p.digits = 2)+#label.x = 0
  #stat_regline_equation(aes(group = season.x), label.x = 0.1)+
  scale_y_continuous(labels = scales::scientific, breaks = c(7e+00, 1e+01), limits = c(7e+00, 1.1e+01))+
  #geom_smooth(method = 'lm', group = 'season.x')+ #, linetype = 'longdash', color = '#3d3b42'
  labs(y = 'ln(ASV19 pseudoabundance)', x = 'Time (h)')+
  scale_color_manual(values = palette_seasons_4)+
  #facet_grid(asv_num ~., scales = 'free')+
  #facet_grid(vars(season.x), scales = 'fixed', switch = 'y')+
  #facet_wrap(treatment.x~asv_num)+
  theme_bw()+
  theme(strip.text = element_text(size = 12), axis.text.x = element_text(angle = 0), legend.position = "none",
        strip.background = element_blank(), strip.placement = 'outside',
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),# aspect.ratio = 1/2,
        plot.margin=unit(c(0.5,0.25,0.5,0.5),"cm")) #, plot.margin=unit(c(0.5,0,0.5,0),"cm") aspect.ratio = 1/2,

##3 temps
#library(ggpmisc)
t3_exemple <- 
  pseudoabund_df_wide_fc_silva %>%
  #select(1:5, 'asv81', 'asv53', 'asv10', 'asv654', 'asv4', 'asv114') %>% # 'asv1',
  #select(1:5, 'asv2', 'asv23', 'asv3', 'asv117', 'asv4', 'asv10') %>%
  select(1:5, 'asv19') %>%
  subset(season != "Early_fall" & 
           treatment == 'PL' &
           time != 't4')%>%
  pivot_longer(cols = starts_with('asv'), names_to = 'asv_num', values_to = 'pseudoabund') %>%
  left_join(reg_all_slopes_chosen_silva_tax, by = 'asv_num') %>%
  ggplot(aes(hours.x, log(pseudoabund), color = season.x))+
  geom_point(shape = 7, size = 3)+
  scale_y_continuous(labels = scales::scientific, breaks = c(7e+00, 1e+01), limits = c(7e+00, 1.1e+01))+
  stat_poly_line(aes(group = season.x))+
  #geom_smooth(method = 'lm')+ #, linetype = 'longdash', color = '#3d3b42'
  labs(y = 'ln(ASV19 pseudoabundance)', x = 'Time (h)', color = 'Season', shape = 'Treatment')+
  scale_color_manual(values = palette_seasons_4)+
  #facet_grid(asv_num ~., scales = 'free')+
  #facet_grid(vars(treatment.x), scales = 'fixed')+
  scale_x_continuous(limits = c(0,50))+
  stat_cor(aes(label = paste(..p.label..)),label.x = 0, p.digits = 1, r.digits = 2)+
  #stat_poly_eq(aes(group = season.x))+
  scale_shape_manual(name = 'Treatment', labels = 'DL', values = 12)+
  #scale_y_continuous(breaks = c(50000, 10000))+ #breaks = c(1, 250, 500, 750, 1000)
  #geom_cor()+
  #coord_capped_cart(bottom='both')+ #, left='both'
  # stat_fit_glance(method = 'lm',
  #                 method.args = list(formula = formula),
  #                 geom = 'text',
  #                 aes(label = paste("P-value = ", signif(..p.value.., digits = 4), sep = "")))+
  theme_bw()+
  theme(strip.text = element_blank(), axis.text.x = element_text(angle = 0), legend.position = "right",
        axis.title.y = element_blank(), strip.background = element_blank(), axis.text.y = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        plot.margin=unit(c(0.5,0,0.5,0),"cm")) #aspect.ratio = 1/2, plot.margin=unit(c(0.5,0,0.5,0),"cm")


#no puc posar junts els treatments perquè els GR están calculats per treatments
#grid.arrange(t4_exemple, t3_exemple, ncol = 2 )

example_gr_calculation <- multi_panel_figure(columns = 2, rows = 1, width = 250, height = 150, 
                                             column_spacing = 0.0, unit = 'mm',
                                             panel_label_type = 'none')

example_gr_calculation  %<>%
  fill_panel(t4_exemple, column = 1, row = 1) %<>%
  fill_panel(t3_exemple, column = 2, row = 1)

ggsave('remiau_example_gr_calculation3.pdf', example_gr_calculation , 
       path = "~/Documentos/Doctorat/REMEI/results/figures/raw/",
       width = 250,
       height = 150,
       units = 'mm')

##Exploració de la phylogenia
reg_all_slopes_chosen_silva_tax %$%
  phylum %>%
  unique()
# 
# [1] "Cyanobacteria"     "Proteobacteria"    "Bacteroidota"      "Crenarchaeota"     "Planctomycetota"  
# [6] "Acidobacteriota"   "Verrucomicrobiota" "Nitrospinota"      "Actinobacteriota"  "Nitrospirota"     
# [11] "Firmicutes"        "Bdellovibrionota"  "Myxococcota"       "Desulfobacterota"  "Deinococcota"     
# [16] "Fusobacteriota"    "Chloroflexi"       "Spirochaetota"     "Campilobacterota"  "Abditibacteriota" 
# [21] "Latescibacterota"  "Methylomirabilota" "Halanaerobiaeota"  "Sumerlaeota"       "Calditrichota"    
# [26] "Gemmatimonadota" 

reg_all_slopes_chosen_silva_tax %>%
  subset(phylum == "Cyanobacteria") %$%
  class %>%
  unique()

reg_all_slopes_chosen_silva_tax %>%
  subset(phylum == "Acidobacteriota") %$%
  class %>%
  unique()

##relació num de treatments vs gr que presentes (classificació com a Bloomer?)

##Question: Do  bloomers respond to general environmental changes or are they specific? Can we classify them?


reg_all_slopes_chosen_silva_tax %>%
  filter(slope_chosen_days > 0,
         pvalue_slope_chosen < 0.05) %>%
  group_by(phylum_f) %>%
  #summarize(n = n())
  filter(n() > 10) %>%
  ungroup() %>%
  group_by(asv_num, season) %>%
  mutate(n = n()) %>%
  distinct(asv_num, treatment, slope_chosen_days, season, n, phylum_f) %>%
  ggplot(aes(n, slope_chosen_days))+
  geom_point(position = position_jitter(0.35), alpha = 0.8, aes(color = season))+
  geom_violin(aes(group = n, color = season), alpha = 0, draw_quantiles = c(0.25, 0.75))+
  geom_smooth(method = 'loess',  aes(color = season))+#color = 'black',
  scale_color_manual(values = palette_seasons_4)+
  scale_y_continuous(expand = c(0, 0), limits = c(0,10.4))+
  labs(color = 'Season', x = 'Number of times each ASV appeared per treatment', y = expression("Growth rate day"^"-1"))+
  #facet_wrap(vars(treatment))+
  facet_grid(.~season)+#treatment~
  theme_bw()

scale_y_continuous(expand = c(0, 0), limits = c(0,10.4))


