sessionInfo()

library(tidyverse)
library(magrittr)
library(ggplot2)

#functions needed
source("src/filter.season.treatment.time.R") ###filtrar les pseudoabund per treatment season i time per calcular regressions
source("src/multiple.linear.regressions.R") ##per fer multiples regressions d'un dataframe
source("src/comparing.reg.R") ##per comparar les taxes de creiement entre 4 i 3 punts (les regressions)

#import data for recalculating DAPIS GR (l'Olga ja ho tenia fet, ho reaprofitem)------
##flow citometry abundances with t1 to calculate GR
fc <- readxl::read_xlsx('data/envdata/Final_Data_FC_PB_Remei_experiments_ed4.xlsx',
                        col_names = T )

fc %>%
  colnames()

fc$treatment <- fc$treatment %>% 
  factor(levels=(c("CL", "CD", "PL", "PD", "DL", "VL")))
fc$season <- fc$season %>% 
  factor(levels=(c("Winter", "Spring", "Summer","Early_fall", "Fall")))

#check if there is lag time
fc %>%
  subset(time == c("t1", "t2", "t3", "t4")) %>%
  ggplot(aes(hours,log(All), shape = treatment))+
  geom_smooth( colour = "grey", se =FALSE,  span = 0.7)+
  geom_point()+
  #scale_color_manual(values = palette_treatments_remei)+
  #geom_line(color = "grey")+
  facet_grid(.~season)+
  theme_bw()

pseudoabund_df_wide_m %>%
  colnames()

rem_relabun_melt <-  read.table("data/rem_relabun_melt.txt", sep="\t", header = TRUE)

abundant_asv <- rem_relabun_melt %>% 
  filter(Abundance > 0.10) %>% #more than 1% of the community at some point
  select(asv_num) %>%
  unique() %>%
  as_tibble()

# rem_relabun_melt_1perc <- rem_relabun_melt %>%
#   left_join(abundant_asv)
reg_all_slopes_chosen_silva_tax %>%
  colnames()
pseudoabund_df_wide_m %>%
  colnames()
pseudoabund_df_long_m_10perc <- pseudoabund_df_wide_m %>% 
  pivot_longer(cols = starts_with("asv"), names_to = "asv_num") %>%
  right_join(abundant_asv, by = "asv_num", copy = TRUE)

pseudoabund_df_long_m_5perc %>%
  colnames()
pseudoabund_df_long_m_5perc  %>%
  #pivot_longer(cols = starts_with("asv")) %>%
  #subset(time == c("t0", "t2", "t3", "t4")) %>%
  ggplot(aes(hours.x, log(value), color = treatment, group = time))+#
  #geom_point(size = 0.2)+
  geom_boxplot(aes(hours.x, log(value), color = treatment), alpha = 0.7)+
  geom_line()+
  #geom_smooth(aes(hours.x, log(value), color = treatment), se =FALSE)+
  scale_color_manual(values = palette_treatments_remei)+
  #geom_line(color = "grey")+
  facet_grid(.~season~treatment)+
  theme_bw()


###DAPIS GR calculated by Olga Sanchez------
GR_dapis_OS <- read.delim2("data/envdata/GR_DAPIS_OS/GR_REMEI_DAPIS_OS_Ed.csv", sep = ";") %>%
  as_tibble()

GR_dapis_OS_filt <- GR_dapis_OS %>%
  filter(treatment %in% c("CL", "CD", "PL", "PD",  "DL", "VL")) %>%
  mutate(PRO = as.numeric(PRO)) 

  ##reorder treatments and seasons for plots
GR_dapis_OS_filt$treatment <- GR_dapis_OS_filt$treatment %>% 
  factor(levels=(c("CL", "CD", "PL", "PD", "DL", "VL")))
GR_dapis_OS_filt$season <- GR_dapis_OS_filt$season %>% 
  factor(levels=(c("Winter", "Spring", "Summer", "Fall")))

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


##statistics------
###Seasons
library('ggpubr')
library('rstatix')
GR_dapis_OS_filt %>%
  colnames()
anova<-aov(PRO~season, data=GR_dapis_OS_filt)
summary(anova) #p<0.05 ***
tukey <- anova %>%
  TukeyHSD()
# anova <- anova %>%
#   TukeyHSD() %>% #para ver los grupos que son significativamente distintos Summer i Winter no son distintos y Fall i Spring tampoco
#   as.matrix()
# anova[,1] %>%
#   as.data.frame() %>%
#   add_significance('season.p.adj')

aov_residuals <- residuals(object = anova)
plot(anova, 1)
plot(anova, 2)
shapiro.test(x = aov_residuals )#From the output, the p-value > 0.05 implying that the distribution of the data are #not significantly different from normal distribution. In other words, we can assume the normality.
##p-value = 0.218 és normal 
#plot(TukeyHSD(anova, conf.level=.95), las = 2)

#Grups estadístics
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

