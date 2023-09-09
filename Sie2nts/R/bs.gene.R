# Generating basis

source(paste(getwd(), "/R/Sie2nts.Legen.v1.R", sep = ""))
source(paste(getwd(), "/R/Sie2nts.Cheby.v1.R", sep = ""))
source(paste(getwd(), "/R/Sie2nts.Four.v1.R", sep = ""))
source(paste(getwd(), "/R/Sie2nts.Csp.v1.R", sep = ""))
source(paste(getwd(), "/R/Sie2nts.db1-20.v1.R", sep = ""))

#' Generate Basis
#' @description bs.gene() generates the value of k-th basis function. (The wavelet basis options return the full table)
#' @param type type indicates which type of basis is used. There are 31 types in this package
#' @param k k-th basis function
#' @param point the number of values got from k-th basis function
#' @param ops ops indicates the function uses existing table or theoretical way to generate, the default option is "auto"
#' @param c c only used in Cspli which indicates the total number of knots to generate, the default is 10, c should not be less than k.(for splines, the true
#' number of basis is c-2+or)
#' @param or indicates the order of spline and only used in Cspli type, default is 4 which indicates cubic spline
#' @return A data frame which contains the value of k-th basis function
#' @references \[3] Chen, Xiaohong. “Large Sample Sieve Estimation of Semi-Nonparametric Models.” Handbook of Econometrics, 6(B): 5549–5632,2007.
#' @export
#' 
#' @examples
#' bs.gene("Legen", 2)
#' bs.gene("tri", 2, 300)
 

bs.gene = function(type, k, point = 200, c=10, or = 4, ops = "auto"){
  # library(ggplot2)
  # Bspline and Cspline indicate the k-th basis under the total k+1 basis. Point number is fixed for wavelet.
  wavelet_basis = c("db1", "db2", "db3", "db4", "db5",
                    "db6", "db7", "db8", "db9", "db10",
                    "db11", "db12", "db13", "db14", "db15",
                    "db16", "db17", "db18", "db19", "db20",
                    "cf1", "cf2", "cf3", "cf4", "cf5"
  )

  if(ops == "non-auto"){
    return(cat("general table first"))

  } else{
    if(type == "Legen"){
      return(Legendre_kth_b(k, point))
    } else if (type == "Cheby"){
      return(Chebyshev_kth_b(k, point))
    } else if (type %in% c("tri", "cos", "sin")){
      return(Fourier_kth_b(k, point, ops = type))
    } else if (type == "Cspli"){
      return(Cspline_kth_b(k, point, c, or=or))
    } else if (type %in% wavelet_basis){
      return(wavelet_kth_b(k, ops = type))
    } else{
      return(stop("Invalid option!"))
    }
  }

}

