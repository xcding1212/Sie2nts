

source(paste(getwd(), "/R/Sie2nts.Legen.v1.R", sep = ""))
source(paste(getwd(), "/R/Sie2nts.Legen.v1.R", sep = ""))
source(paste(getwd(), "/R/Sie2nts.Cheby.v1.R", sep = ""))
source(paste(getwd(), "/R/Sie2nts.Four.v1.R", sep = ""))
source(paste(getwd(), "/R/Sie2nts.Csp.v1.R", sep = ""))
source(paste(getwd(), "/R/Sie2nts.db1-20.v1.R", sep = ""))





#' Predicting Time Series With H Steps
#' @description predict.nts() predicts the time series data basis on the estimation.
#' @param ts The data set which is a time series data typically
#' @param esti.li The output from fix.fit() or sie.auto.fit() function
#' @param h h indicates the number of forecasting points
#'
#' @return A vector which contains h forecasting points
#' @export
#'
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
#' res1.2 = fix.fit(time.series, 5, 1, "Legen")
#' sie.predict(time.series, res1.2, 5)
#' res2.2 = fix.fit(time.series, 2, 1, "db10")
#' sie.predict(time.series, res2.2, 5)
#' }

sie.predict = function(ts, esti.li, h){
  return(predict.legen(ts, esti.li, h))
}
