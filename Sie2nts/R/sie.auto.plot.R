# visulizaed result

source(paste(getwd(), "/R/Sie2nts.Legen.v1.R", sep = ""))
source(paste(getwd(), "/R/Sie2nts.Cheby.v1.R", sep = ""))
source(paste(getwd(), "/R/Sie2nts.Four.v1.R", sep = ""))
source(paste(getwd(), "/R/Sie2nts.Csp.v1.R", sep = ""))
source(paste(getwd(), "/R/Sie2nts.db1-20.v1.R", sep = ""))
source(paste(getwd(), "/R/Sie2nts.plot.v1.R", sep = ""))



#' Plot the Estimate Results by Automatic Fitting
#' @description sie.auto.plot() visualizes the estimation of coefficient, gives the elbow plot and represents the cross validation result.
#' @param res.auto.fit the output from sie.auto.fit() function
#' @param type type indicates which type of basis is used (There are 31 types in this package)
#' @param title give the title for the auto estimate plot
#'
#' @return A list which contains 3 plot related to the estimation of coefficient, Elbow point and cross validation in order
#' @export
#' @examples
#' \donttest{
#' set.seed(137)
#' time.series = c()
#' n = 1024
#' v = 25
#' w = rnorm(n, 0, 1) / v
#' x_ini = runif(1,0,1)

#' for(i in 1:n){
#'   if(i == 1){
#'     time.series[i] = 0.2 + 0.6*cos(2*pi*(i/n))*x_ini  + w[i] #
#'   } else{
#'     time.series[i] = 0.2 + 0.6*cos(2*pi*(i/n))*time.series[i-1] + w[i]
#'   }
#' }
#'
#' res1 = sie.auto.fit(time.series, type = "db10")
#' sie.auto.plot(res1, "db10")
#' res2 = sie.auto.fit(time.series, type = "Legen", method = "LOOCV")
#' sie.auto.plot(res2, "Legen")
#' }




# the plot for auto fit
sie.auto.plot = function(res.auto.fit, type, title = ""){
  aux.res = list()
  aux.res[[1]] = fit.plot.estimate.aux(res.auto.fit$Estimate, res.auto.fit$BC[1], res.auto.fit$BC[2], type, title)
  aux.res[[2]] = fit.plot.elbow(res.auto.fit, res.auto.fit$BC[1], res.auto.fit$BC[2], type, title)
  aux.res[[3]] = fit.plot.cvm(res.auto.fit, type)
  return(aux.res)
}
