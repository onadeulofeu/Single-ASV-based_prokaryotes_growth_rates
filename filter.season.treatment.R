#' Title
#'
#' @param data 
#' @param treatment_col 
#' @param treatment_keep 
#' @param season_col 
#' @param season_keep 
#'
#' @return
#' @export
#'
#' @examples
filter.season.treatment <- function(data, treatment_col, treatment_keep, season_col, season_keep){
  my_var1 <- rlang::enquo(treatment_col)
  my_var2 <- rlang::enquo(season_col)
  data_filt <- data %>%
    filter(!!my_var1 == treatment_keep & !!my_var2 == season_keep)
  return(data_filt)
}