#' Filter season treatment and time in remei tibbles
#'
#' @param data a table with environmental data 
#' @param treatment_col column where treatments are annotated
#' @param treatment_keep treatments we want to keep one or more
#' @param season_col column where seasons are annotated
#' @param season_keep seasons we want to work with one or more
#' @param time_col column where time id are annotated 
#' @param time_keep times we want  to work with one or more
#'
#' @return a dataset filtered
#' @export tidyverse
#' @import dplyr
#'
#' @examples
filter.season.treatment.time <- function(data, treatment_col, treatment_keep, season_col, season_keep, time_col, time_keep){
  my_var1 <- rlang::enquo(treatment_col)
  my_var2 <- rlang::enquo(season_col)
  my_var3 <- rlang::enquo(time_col)
  data_filt <- data %>%
    filter(!!my_var1 %in% treatment_keep & !!my_var2 %in% season_keep & !!my_var3 %in% time_keep)
  return(data_filt)
}