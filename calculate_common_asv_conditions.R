#' Calculate common ASVs between conditions
#'
#' @param data a long tibble with a column with conditions (season and treatment) and a column with ASV number with all ASVs growing in each condition
#' @param sample1 condition 1 to be compraed with condition two
#' @param sample2 condition 2 to be compared with condition one
#'
#' @return a tibble with the number of ASVs common between the two conditions (sample 1 and sample 2)
#' @export
#' @import dplyr
#' @import tidyverse
#'
#' @examples
common.asvs <- function(data, sample1, sample2){
  
  { 
    data1 <- data %>%
      subset(seas_treat == sample1)
    data2 <- data %>%
      subset(seas_treat == sample2)
  }
  {
    sample12 <- paste0(sample1, '_', sample2) %>%
      as.character() 
    result <-  data1 %>%
      left_join(data2, by = 'asv_num') %>%
      filter(seas_treat.y != is.na(seas_treat.y)) %>%
      group_by(seas_treat.x) %>%
      dplyr::summarise(num_common_asvs = n()) %>%
      dplyr::mutate(seas_treat.x = str_replace(seas_treat.x, !!{sample1}, !!{sample12}))
  }
  return(result)
}