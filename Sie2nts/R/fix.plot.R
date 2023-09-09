# visulizaed result


source(paste(getwd(), "/R/Sie2nts.Legen.v1.R", sep = ""))
source(paste(getwd(), "/R/Sie2nts.Cheby.v1.R", sep = ""))
source(paste(getwd(), "/R/Sie2nts.Four.v1.R", sep = ""))
source(paste(getwd(), "/R/Sie2nts.Csp.v1.R", sep = ""))
source(paste(getwd(), "/R/Sie2nts.db1-20.v1.R", sep = ""))
source(paste(getwd(), "/R/Sie2nts.plot.v1.R", sep = ""))

#' Plot Results of Estimating
#' @description fix.plot() visualizes the estimation of coefficient.
#' @param res.fix.fit the output from fix.fit() function
#' @param type type indicates which type of basis is used (There are 31 types in this package)
#' @param title give the title for the fixed estimate plot
#'
#' @return A list which contains 3 plot related to the estimation of coefficient, Elbow point and cross validation in order
#' @export
#' @examples
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
#' res1 = fix.fit(time.series, 5, 1, type = "Legen")
#' fix.plot(res1, "Legen")




# the plot for fix fit
fix.plot = function(res.fix.fit, type, title = ""){
  aux.res = list()

  c = length(res.fix.fit$ols.coef) / length(res.fix.fit$ts.coef)
  b = length(res.fix.fit$ts.coef) - 1
  aux.res[[1]] = fit.plot.estimate(res.fix.fit, c, b, type, title)
  return(aux.res)
}

