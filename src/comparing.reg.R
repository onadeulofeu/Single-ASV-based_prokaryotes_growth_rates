#' Comparing regressions
#'
#'This function is used to compare growth rates values between 3 time points or 4 time points
#'and decide to keep the one that is significant and higher
#'
#' @param df1 dataframe with all the regressions 3 points
#' @param df2 dataframe with all the regressions performed with 4 points
#' @param treatment colum from the dataframes varible 1 for deleting repeated slopes values
#' @param season colum from the dataframes varible 2 for deleting repeated slopes values
#' @param asv_num colum from the dataframes varible 3 for deleting repeated slopes values
#' @param slope columnn from the dataframes varible 4 for deleting repeated slopes values
#' @param slope.pval columnn from dataframes varible 5 for deleting repeated slopes values
#'
#' @return a tibble with a new column with the slopes chosen
#' @export
#' @import dplyr
#' @import tidyverse
#'
#' @examples
comparing.reg <- function(df1, df2, treatment, season, asv_num, slope, slope.pval){
  var1 <- rlang::enquo(treatment)
  var2 <- rlang::enquo(season)
  var3 <- rlang::enquo(asv_num)
  var4 <- rlang::enquo(slope)
  var5 <- rlang::enquo(slope.pval)
  df1_filt <- df1 %>% 
    distinct(!!var1, !!var2, !!var3, !!var4, !!var5) #primer eliminem duplicats perquè tinguin les mateixes files els dos datasets
  df2_filt <- df2 %>% 
    distinct(!!var1, !!var2, !!var3, !!var4, !!var5)
  comp_df <- df1_filt %>% 
    left_join(df2_filt, by = c("treatment", "season", "asv_num")) 
  # comp_df <- df1 %>% ##funció per comparar les taxes de creixement (slopes) entre 3 i 4 temps.
  comp_df_ed <- comp_df %>% mutate(slope_chosen = case_when(slope.x > slope.y & slope.pval.x < 0.05 ~ slope.x,
                                                            slope.x < slope.y ~ slope.y),
                                   pvalue_slope_chosen =case_when(slope.x > slope.y & slope.pval.x < 0.05 ~ slope.pval.x,
                                                                  slope.x < slope.y ~ slope.pval.y)) %>%
    as_tibble()
  
  print(comp_df_ed)
  return(comp_df_ed)
}
