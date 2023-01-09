#' Growth rates distribution at different taxonomic ranks
#'
#' Aquesta funció serveix per dibuixar una columna de plots a diferents nivells taxonòmics d'un phylum concret,
#' dibuixem un gràfic d'histogrames amb la distribució dels growth rates. 
#' 
#' @param data tibble with growth rates per day and taxonomy
#' @param phylum_to_explore phylum we want to plot
#' @param bin_width_c with of the bins at class level
#' @param bin_width_o with of the bins at order level
#' @param bin_width_f with of the bins at family level
#' @param bin_width_a with of the bins at asv level
#' @param title Phylum we are plotting
#' @param axis_y_title tittle of y axis
#'
#' @return 5 plots in a column at  different taxonomic levels with the counts of growth rates
#' @export 
#' @import tidyverse
#' @import ggplot2
#'
#' @examples
growthrates.distribution.tax.rank.ridges <- function(data, phylum_to_explore, title,
                                                     axis_y_title){
  data_filt <- data %>%
    subset(phylum == phylum_to_explore) %>%
    filter(pvalue_slope_chosen < 0.05) %>% #& intercept.pval < value_intercept_pval
    filter(slope_chosen_days > 0)
  
  if(!"data.frame" %in% attr(data_filt, "class")) {
    abort("data should be a data frame")
  } 
  
  # plot_phylum <- data_filt %>% 
  #   group_by(phylum_f) %>%
  #   #summarize(n = n())
  #   filter(n() > 2) %>% #, .preserve = FALSE
  #   ggplot(aes(y = phylum_f, x = slope_chosen_days, fill = phylum_f))+
  #   # geom_histogram(aes(x=slope_chosen_days, group = class), colour = "#5F9660",  fill = "#5F9660", #adjust= 0.5,
  #   #                binwidth = bin_width_c,
  #   #                show.legend = FALSE, alpha = 0.2)+# bins = bin_width_c,
  #   geom_density_ridges(aes(fill = phylum_f), scale = 0.5, alpha = 0.6)+
  #   labs(x= "", y = "", title = title)+
  #   #stat_bin(binwidth = bin_width_c)+
  #   scale_x_continuous(limits = c(0,11))+
  #   #scale_y_continuous(limits = c(y1_ph,y2_ph))+
  #   scale_fill_manual(values = palette_phylums)+
  #   theme_bw()+
  #   theme(strip.text.x = element_text(size = 0), axis.text.x = element_text(angle = 0, size = 0), 
  #         axis.ticks = element_blank(), legend.position = "none", plot.title = element_text(size = 10),
  #         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  #         plot.margin= unit(c(0.2, 0.2, 0.2, 0.2), "lines")) #top, right, bottom, left
  
  plot_class <- data_filt %>%
    group_by(class_f) %>%
    filter(n() > 2) %>%
    ggplot(aes(y = class_f, x = slope_chosen_days, fill = phylum_f))+
    # geom_histogram(aes(x=slope_chosen_days, group = class), colour = "#5F9660",  fill = "#5F9660", #adjust= 0.5,
    #                binwidth = bin_width_c,
    #                show.legend = FALSE, alpha = 0.2)+# bins = bin_width_c,
    geom_density_ridges(aes(fill = phylum_f), scale = 1, alpha = 0.6)+
    labs(x= "", y = "")+
    #stat_bin(binwidth = bin_width_c)+
    scale_x_continuous(limits = c(0,11))+
    scale_fill_manual(values = palette_phylums)+
    #scale_y_continuous(limits = c(y1_cl,y2_cl))+
    theme_bw()+
    theme(strip.text.x = element_text(size = 0), axis.text.x = element_text(angle = 0), 
          axis.ticks = element_blank(), legend.position = "none", plot.title = element_text(size = 10),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.margin= unit(c(0.2, 0.2, 0.2, 0.2), "lines")) #top, right, bottom, left
  
  plot_order <- 
    data_filt %>%
    group_by(order) %>%
    filter(n() > 2) %>%
    ggplot(aes(y = interaction(order_f, class_f, sep = '\n'), x = slope_chosen_days, fill = phylum_f))+
    # geom_histogram(aes(x=slope_chosen_days, group = order), colour = "#B8102A", fill = "#B8102A",
    #                #adjust= 0.5, 
    #                binwidth = bin_width_o,
    #                #bins = bin_width_o,
    #                show.legend = FALSE, alpha = 0.2)+
    geom_density_ridges(aes(fill = phylum_f), scale = 2, alpha = 0.6)+
    labs(x= "", y = "")+
    #stat_bin(binwidth = bin_width_o)+
    scale_x_continuous(limits = c(0,11))+
    scale_fill_manual(values = palette_phylums)+
    #scale_y_continuous(limits = c(y1_o,y2_o))+
    theme_bw()+ 
    theme(strip.text.x = element_text(size = 0), axis.text.x = element_text(angle = 0), 
          axis.ticks = element_blank(), legend.position = "none",
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.margin= unit(c(0.2, 0.2, 0.2, 0.2), "lines"))
  
  plot_family <- 
    data_filt %>%
    group_by(family) %>%
    filter(n() > 2) %>%
    ggplot(aes(y = interaction( family_f, order_f, sep = '\n'), x = slope_chosen_days, fill = phylum_f))+
    # geom_histogram(aes(x=slope_chosen_days, group = family), colour = "#FFDF26", fill = "#FFDF26", #adjust= 0.5, 
    #                # bins = bin_width_f,
    #                binwidth = bin_width_f,
    #                show.legend = FALSE, alpha = 0.2)+
    geom_density_ridges(aes(fill = phylum_f), scale = 3, alpha = 0.6)+
    labs(x= "")+
    #stat_bin(binwidth = bin_width_f)+
    scale_x_continuous(limits = c(0,11))+
    scale_fill_manual(values = palette_phylums)+
   # scale_y_continuous(limits = c(y1_f,y2_f))+
    theme_bw()+ 
    theme(strip.text.x = element_text(size = 0), axis.text.x = element_text(angle = 0), 
          axis.ticks = element_blank(), legend.position = "none",axis.title.y = element_text(size = 0),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          plot.margin= unit(c(0.2, 0.2, 0.2, 0.2), "lines"))
  
  plot_asv <- data_filt %>%
    group_by(asv_num) %>%
    filter(n() > 2) %>%
    ggplot(aes(y = asv_f, x = slope_chosen_days, fill = phylum_f))+
    # geom_histogram(aes(x=slope_chosen_days, group = asv_num), colour  = "#4F507F", fill = "#4F507F", 
    #                #adjust= 0.5, #bins = bin_width_a,
    #                binwidth = bin_width_a,
    #                show.legend = FALSE, alpha = 0.2)+
    geom_density_ridges(aes(fill = phylum_f), scale = 4, alpha = 0.6)+
    labs(x= expression("Growth rate day"^"-1"), y = "")+
    scale_x_continuous(limits = c(0,11))+
    scale_fill_manual(values = palette_phylums)+
   # scale_y_continuous(limits = c(y1_a,y2_a))+
    #stat_bin(binwidth = bin_width_a)+
    theme_bw()+ 
    theme(strip.text.x = element_text(size = 6), axis.text.x = element_text(angle = 0), legend.position = "none",
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 0), 
          plot.margin= unit(c(0.2, 0.2, 0.2, 0.2), "lines"))
  
  distrubution_plots_tax_rank <- grid.arrange(
    #plot_phylum,
                                              plot_class, 
                                              plot_order, 
                                              plot_family, 
                                              plot_asv,
                                              nrow = 1)
  print(distrubution_plots_tax_rank)
  return(distrubution_plots_tax_rank)
}
