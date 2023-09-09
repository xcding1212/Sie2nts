# Estimate coefficients for time series


source(paste(getwd(), "/R/Sie2nts.Legen.v1.R", sep = ""))
source(paste(getwd(), "/R/Sie2nts.Cheby.v1.R", sep = ""))
source(paste(getwd(), "/R/Sie2nts.Four.v1.R", sep = ""))
source(paste(getwd(), "/R/Sie2nts.Csp.v1.R", sep = ""))
source(paste(getwd(), "/R/Sie2nts.db1-20.v1.R", sep = ""))

#' Estimate the Coefficients of Auto-Regressive (AR) Model by User Specifying
#' @description fix.fit() estimates the coefficients of AR model by sieve methods with user specifying.
#' @param ts ts is the data set which is a time series data typically
#' @param c c indicates the number of basis used to estimate (For wavelet, the real number of basis is 2^c.
#'     For Cubic Spline, the real number of basis is c-2+or)
#' @param b b is the lag for auto-regressive model
#' @param type type indicates which type of basis is used. There are 31 types in this package
#' @param or indicates the order of spline and only used in Cspli type, default is 4 which indicates cubic spline
#' @param m m indicates the number of points of coefficients to estimate
#'
#' @return A list contains 3 objects, the first is a matrix which contains estimates for each basis used
#'     in OLS, the second is a list contains estimates for coefficients in AR model and the last
#'     is a vector contains residuals
#' @export
#'
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
#' res = fix.fit(time.series, c=5, b=1, type = "Legen")
#' cat(res$ols.coef)
#' plot.ts(res$ts.coef[[1]])
#' plot.ts(res$Residuals)



fix.fit = function(ts, c, b, type, or=4, m=500){
  wavelet_basis = c("db1", "db2", "db3", "db4", "db5",
                    "db6", "db7", "db8", "db9", "db10",
                    "db11", "db12", "db13", "db14", "db15",
                    "db16", "db17", "db18", "db19", "db20",
                    "cf1", "cf2", "cf3", "cf4", "cf5"
  )

  if(type == "Legen"){
    return(fix.fit.legen(ts, c, b, m))
  } else if (type == "Cheby"){
    return(fix.fit.cheby(ts, c, b, m))
  } else if (type %in% c("tri", "cos", "sin")){
    return(fix.fit.four(ts, c, b, m, ops = type))
  } else if (type == "Cspli"){
    return(fix.fit.cspline(ts, c, b, or=or, m))
  } else if (type %in% wavelet_basis){
    return(fix.fit.wavelet(ts, c, b, m, ops = type))
  } else{
    return(stop("Invalid option!"))
  }
}


