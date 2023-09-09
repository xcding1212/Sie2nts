# Estimate coefficients of time series by cv

source(paste(getwd(), "/R/Sie2nts.Legen.v1.R", sep = ""))
source(paste(getwd(), "/R/Sie2nts.Cheby.v1.R", sep = ""))
source(paste(getwd(), "/R/Sie2nts.Four.v1.R", sep = ""))
source(paste(getwd(), "/R/Sie2nts.Csp.v1.R", sep = ""))
source(paste(getwd(), "/R/Sie2nts.db1-20.v1.R", sep = ""))


#' Estimate the Coefficients of Auto-Regressive (AR) Model Automatically
#' @description sie.auto.fit() estimates the coefficients of AR model by sieve methods with 2 cross validation methods and elbow method.
#' @param ts ts is the data set which is a time series data typically
#' @param type type indicates which type of basis is used. There are 31 types in this package
#' @param or indicates the order of spline and only used in Cspli type, default is 4 which indicates cubic spline
#' @param method method indicates which method used to choose optimal parameters, 3 methods in this package can be used
#' @param m m indicates the number of points of coefficients to estimate
#' @param threshold threshold determines the bound for Elbow method
#'
#' @return A list contains 4 objects, the first is estimates for coefficients
#'    in AR model, the second is cross validation table, the third is estimates for each basis used
#'    in OLS and the last is optimal parameters
#' @export



sie.auto.fit = function(ts, type, or=4, method = "LOOCV", m = 500, threshold = 0){
  wavelet_basis = c("db1", "db2", "db3", "db4", "db5",
                    "db6", "db7", "db8", "db9", "db10",
                    "db11", "db12", "db13", "db14", "db15",
                    "db16", "db17", "db18", "db19", "db20",
                    "cf1", "cf2", "cf3", "cf4", "cf5"
  )

  if(type == "Legen"){
    return(auto.fit.legen(ts, m = m, method = method, threshold = threshold))
  } else if (type == "Cheby"){
    return(auto.fit.cheby(ts, m = m, method = method, threshold = threshold))
  } else if (type %in% c("tri", "cos", "sin")){
    return(auto.fit.four(ts, m = m, ops = type, method = method, threshold = threshold))
  } else if (type == "Cspli"){
    return(auto.fit.cspline(ts, or=4,  m = m, method = method, threshold = threshold))
  } else if (type %in% wavelet_basis){
    return(auto.fit.wavelet(ts,  m = m, ops = type, method = method, threshold = threshold))
  } else{
    return(stop("Invalid option!"))
  }
}


