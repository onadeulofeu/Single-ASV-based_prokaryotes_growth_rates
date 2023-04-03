#' Calculate exclusive ASVs for a condition 
#'
#' @param data a tibble with a column of ASVs nº of exclusive ASVs and a column of their conditions (season and treatment)
#' @param sample the sample (condition: season and treatment) we want to calculate the number of exclusive ASVs
#'
#' @return a tibble with the number of ASVs exclusive per conditions (treatment and season)
#' @export
#' @import tidyverse
#' @import dplyr
#'
#' @examples
exclusive.asvs <- function(data, sample){
  asvs_present_all_other_samples <- data %>%
    dplyr::filter(seas_treat != sample) %>%
    group_by(asv_num) %>%
    distinct()
  exclusive_asv_sample <- data %>%
    dplyr::filter(seas_treat == sample) %>%
    dplyr::filter(!asv_num %in% asvs_present_all_other_samples$asv_num)
  return(exclusive_asv_sample)
  
}