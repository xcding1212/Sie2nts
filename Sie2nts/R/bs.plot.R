# Plot basis

source(paste(getwd(), "/R/Sie2nts.Legen.v1.R", sep = ""))
source(paste(getwd(), "/R/Sie2nts.Cheby.v1.R", sep = ""))
source(paste(getwd(), "/R/Sie2nts.Four.v1.R", sep = ""))
source(paste(getwd(), "/R/Sie2nts.Csp.v1.R", sep = ""))
source(paste(getwd(), "/R/Sie2nts.db1-20.v1.R", sep = ""))

#' Plots of Basis
#' @description bs.plot() generates the plot of first k basis function.
#' @param type type indicates which type of basis is used (There are 31 types in this package)
#' @param k The k is the number of basis functions represented (If wavelet are chosen, the real number of basis is 2^k. If
#'     Cspli is chosen, the real number of basis is k-2+or)
#' @param title give the title for the basis plot
#' @param or indicates the order of spline and only used in Cspli type, default is 4 which indicates cubic spline
#' @return The plot of 1 to k basis functions
#' @export
#'
#' @examples
#' bs.plot("Legen", 2)
#' bs.plot("tri", 3)


bs.plot = function(type, k, or = 4, title = ""){
  # library(ggplot2)
  # library(splines)
 # for Bspli, the true number of basis function is k+1,  and k+2 for Cspli. For wavelet is 2^k
  wavelet_basis = c("db1", "db2", "db3", "db4", "db5",
                    "db6", "db7", "db8", "db9", "db10",
                    "db11", "db12", "db13", "db14", "db15",
                    "db16", "db17", "db18", "db19", "db20",
                    "cf1", "cf2", "cf3", "cf4", "cf5"
  )

  if(type == "Legen"){
    return(legendre_plot(k, title))
  } else if (type == "Cheby"){
    return(chebyshev_plot(k, title))
  } else if (type %in% c("tri", "cos", "sin")){
    return(fourier_plot(k, ops = type, title))
  } else if (type == "Cspli"){
    return(cspline_plot(k, or = or, title))
  } else if (type %in% wavelet_basis){
    return(dbplot(k, ops = type, title))
  } else{
    return(stop("Invalid option!"))
  }


}
