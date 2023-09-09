# pacf

source(paste(getwd(), "/R/Sie2nts.Legen.v1.R", sep = ""))
source(paste(getwd(), "/R/Sie2nts.Cheby.v1.R", sep = ""))
source(paste(getwd(), "/R/Sie2nts.Four.v1.R", sep = ""))
source(paste(getwd(), "/R/Sie2nts.Csp.v1.R", sep = ""))
source(paste(getwd(), "/R/Sie2nts.db1-20.v1.R", sep = ""))

#' Generate Partial Autocorrelation Function (PACF) by User Specifying
#' @description fix.pacf() generates the PACF with fixed tuning parameters.
#' @param ts ts is the data set which is a time series data typically
#' @param c c indicates the number of basis used to Estimate (For wavelet, the number of basis is 2^c. If
#'     Cspli is chosen, the real number of basis is c-2+or)
#' @param lag lag b is the lag for auto-regressive model
#' @param type type indicates which type of basis is used (There are 31 types in this package)
#' @param or or indicates the order of spline and only used in Cspli type, default is 4 which indicates cubic spline
#' @param m m indicates the number of points of coefficients to Estimate
#'
#' @return A vector which contains the PACF with specific lag
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
#' fix.pacf(time.series, c=5, lag = 1, type = "Legen")



# with what number of lag to generate
fix.pacf = function(ts, c, lag, type, or = 4, m=500){
  wavelet_basis = c("db1", "db2", "db3", "db4", "db5",
                    "db6", "db7", "db8", "db9", "db10",
                    "db11", "db12", "db13", "db14", "db15",
                    "db16", "db17", "db18", "db19", "db20",
                    "cf1", "cf2", "cf3", "cf4", "cf5"
  )

  if(type == "Legen"){
    res = fix.fit.legen(ts, c, lag, m)
    return(res$ts.coef[[lag+1]])
  } else if (type == "Cheby"){
    res = fix.fit.cheby(ts, c, lag, m)
    return(res$ts.coef[[lag+1]])
  } else if (type %in% c("tri", "cos", "sin")){
    res = fix.fit.four(ts, c, lag, m, ops = type)
    return(res$ts.coef[[lag+1]])
  } else if (type == "Cspli"){
    res = fix.fit.cspline(ts, c, lag,or=or, m)
    return(res$ts.coef[[lag+1]])
  } else if (type %in% wavelet_basis){
    res = fix.fit.wavelet(ts, c, lag, m, ops = type)
    return(res$ts.coef[[lag+1]])
  } else{
    return(stop("Invalid option!"))
  }
}

