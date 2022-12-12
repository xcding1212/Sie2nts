# stability test


source(paste(getwd(), "/R/Sie2nts.Legen.v1.R", sep = ""))
source(paste(getwd(), "/R/Sie2nts.Cheby.v1.R", sep = ""))
source(paste(getwd(), "/R/Sie2nts.Four.v1.R", sep = ""))
source(paste(getwd(), "/R/Sie2nts.Csp.v1.R", sep = ""))
source(paste(getwd(), "/R/Sie2nts.db1-20.v1.R", sep = ""))

#' The Test of Stability for AR Approximations With Fixed Parameters
#' @description fix.test() generates a test of Stability for AR Approximations with fixed parameters.
#' @param ts ts is the data set which is a time series data typically
#' @param c c indicates the number of basis used to estimate (For wavelet, the number of basis is 2^c. If
#'     Cspli is chosen, the real number of basis is c-2+or)
#' @param b b is the lag for auto-regressive model
#' @param type type indicates which type of basis is used. There are 31 types in this package
#' @param or indicates the order of spline and only used in Cspli type, default is 4 which indicates cubic spline
#' @param B.s the number of statistics used in multiplier bootstrap, the default value is 1000
#' @param m the number of window size used in multiplier bootstrap, the default value is 0 which uses the
#'     minimum volatility method to determine the number
#'
#' @return p value of the test
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
#' fix.test(time.series, 5, 1, type = "Legen")
#' fix.test(time.series, 3, 1, type = "db10")
#' }

fix.test = function(ts, c, b, type, or=4, B.s = 1000, m = 0){
  wavelet_basis = c("db1", "db2", "db3", "db4", "db5",
                    "db6", "db7", "db8", "db9", "db10",
                    "db11", "db12", "db13", "db14", "db15",
                    "db16", "db17", "db18", "db19", "db20",
                    "cf1", "cf2", "cf3", "cf4", "cf5"
  )

  if(type == "Legen"){
    return(fix.test.legen(ts, c, b, B.s = B.s, m = m))
  } else if (type == "Cheby"){
    return(fix.test.cheby(ts, c, b, B.s = B.s, m = m))
  } else if (type %in% c("tri", "cos", "sin")){
    return(fix.test.four(ts, c, b, ops = type, B.s = B.s, m = m))
  } else if (type == "Cspli"){
    return(fix.test.cspline(ts, c, b, or=or, B.s = B.s, m = m))
  } else if (type %in% wavelet_basis){
    return(fix.test.wavelet(ts, c, b, ops = type, B.s = B.s, m = m))
  } else{
    return(cat("Invalid option!"))
  }
}






