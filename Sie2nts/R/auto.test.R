# Stability test by cv


source(paste(getwd(), "/R/Sie2nts.Legen.v1.R", sep = ""))
source(paste(getwd(), "/R/Sie2nts.Cheby.v1.R", sep = ""))
source(paste(getwd(), "/R/Sie2nts.Four.v1.R", sep = ""))
source(paste(getwd(), "/R/Sie2nts.Csp.v1.R", sep = ""))
source(paste(getwd(), "/R/Sie2nts.db1-20.v1.R", sep = ""))
source(paste(getwd(), "/R/sie.auto.fit.R", sep = ""))



#' The Test of Stability for Auto-Regressive (AR) Approximations Automatically
#' @description auto.test() generates a test of Stability for AR Approximations by choosing tuning parameter automatically.
#' @param ts ts is the data set which is a time series data typically
#' @param or indicates the order of spline and only used in Cspli type, default is 4 which indicates cubic spline
#' @param type type indicates which type of basis is used. There are 31 types in this package
#' @param alpha level of the test
#' @param method method indicates which method used to choose optimal parameters, 3 methods in this package can be used
#' @param threshold threshold determines the bound for Elbow method
#' @param B.s the number of statistics used in multiplier bootstrap, the default value is 1000
#'
#' @return p value of the test
#' @export


auto.test = function(ts, or=4, type, alpha = 0.05, method = "LOOCV", threshold = 0, B.s = 1000){
  wavelet_basis = c("db1", "db2", "db3", "db4", "db5",
                    "db6", "db7", "db8", "db9", "db10",
                    "db11", "db12", "db13", "db14", "db15",
                    "db16", "db17", "db18", "db19", "db20",
                    "cf1", "cf2", "cf3", "cf4", "cf5"
  )

  if(type == "Legen"){
    res = sie.auto.fit(ts, type = type, method = method, threshold = threshold)
    cat(paste("The significant level is",alpha, sep = " "))
    return(fix.test.legen(ts, res$BC[1], res$BC[2],  B.s  = B.s, m = 0))
  } else if (type == "Cheby"){
    res = sie.auto.fit(ts, type = type, method = method, threshold = threshold)
    cat(paste("The significant level is",alpha, sep = " "))
    return(fix.test.cheby(ts, res$BC[1], res$BC[2],  B.s  = B.s, m = 0))

  } else if (type %in% c("tri", "cos", "sin")){
    res = sie.auto.fit(ts, type = type, method = method, threshold = threshold)
    cat(paste("The significant level is",alpha, sep = " "))
    return(fix.test.four(ts, res$BC[1], res$BC[2], ops = type,  B.s  = B.s, m = 0))
  } else if (type == "Cspli"){
    res = sie.auto.fit(ts, or =or, type = type, method = method, threshold = threshold)
    cat(paste("The significant level is",alpha, sep = " "))
    return(fix.test.cspline(ts, res$BC[1], res$BC[2], or = or,  B.s  = B.s, m = 0))
  } else if (type %in% wavelet_basis){
    res = sie.auto.fit(ts, type = type, method = method, threshold = threshold)
    cat(paste("The significant level is",alpha, sep = " "))
    return(fix.test.wavelet(ts, res$BC[1], res$BC[2], ops = type,  B.s  = B.s, m = 0))
  } else{
    return(stop("Invalid option!"))
  }
}





