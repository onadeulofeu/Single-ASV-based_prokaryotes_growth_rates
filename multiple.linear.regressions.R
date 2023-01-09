#' Multiple linear regressions
#'
#' @param abund pseudoabundances for each otu in columns
#' @param env.var environmental variable for doing the regressions
#'
#' @return
#' @export table with all the statistics that needs to be transposed afterwords 
#'
#' @examples
multiple.linear.regressions  <- function(abund, env.var){
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
