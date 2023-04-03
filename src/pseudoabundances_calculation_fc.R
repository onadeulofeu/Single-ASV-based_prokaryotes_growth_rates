#' Pseudoabundances calculation
#'
#' @param data a tidy tibble (melted) with relative abundances and DAPI counts per sample as columns
#' @param grp a vector with all the metadata we want to keep for plotting the results and to operate by group
#' @param Abundance relative abundances in 0-1 of each asv
#' @param DAPI counts per each sample
#'
#' @return
#' @export
#' @import tidyverse
#' @import magrittr
#' @import dplyr
#'
#' @examples
pseudoabun_fc <- function(data, grp = NULL, abundance, fc){
  if("package:plyr" %in% search()) detach("package:plyr", unload=TRUE) 
  ##unload library plyr (it doesn't work with this library uploaded) because 
  ##plyr::summari[sz]e on a grouped tbl_df seems to destroy the grouping.
  pseudabund_df <- data %>%
    group_by_at(grp) %>% 
    dplyr::summarize(pseudoabundance = (Abundance*All), .groups = "keep") %>%
    unnest(pseudoabundance) %>%
    ungroup() %>%
    as_tibble()
  print(head(pseudabund_df))
  return(pseudabund_df)
}
