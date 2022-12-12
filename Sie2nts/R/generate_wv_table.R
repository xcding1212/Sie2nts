

s_phi.0 = function(t){
  return(ifelse(t<1 & t>=0, 1,0))
}


w_find = function(valtable, inter, upper, val){
  if(val < 0 | val >= upper){
    return(0)
  } else{
    return(valtable[which(abs(inter-val) == min(abs(inter-val)))])
  }
}

s_phi = function(coeffi, t){
  aux = sqrt(2)*coeffi[1]*s_phi.0(2*t)
  for(i in 2:length(coeffi)){
    aux = aux + sqrt(2)*coeffi[i]*s_phi.0(2*t - (i-1))
  }
  return(aux)
}

# 1 --- 3
# 2 --- 6
# 3 --- 9
# 4 --- 12
# 5 --- 15

#' Generate wavelet table
#' @description generate_wv_table() generates wavelet table by using Cascade algorithm
#' @param coeffi It contains low-pass filter coefficients to generate wavelet table, the coefficients can get from the website
#' @param len The total number of points
#' @param ite The number of iterations to generate table, more numbers of iterations means more accuracy
#'
#' @return A vector contains the value of scaling wavelet functions
#' @export
#' @references
#' [1] Wasilewski, F. (n.d.). Wavelet browser by pywavelets. Retrieved November 18, 2022, from
#' http://wavelets.pybytes.com/
#'
#' @examples
#' This example generate Daubechies 2 wavelet scaling function
#' coef = c(0.48296291314469025, 0.836516303737469, 0.22414386804185735, -0.12940952255092145)
#' tb = generate_wv_table(coef, 31680, 5)
#' plot(tb)

generate_wv_table = function(coeffi, len, ite){
  t=seq(0, length(coeffi) - 1, length.out = len)
  upper = length(coeffi) - 1
  aux.ite = s_phi(coeffi, t)
  aux.table = c()
  val = 0
  for(i in 1:ite){
    for(ind in 1:len){
      for(i in 0:(length(coeffi)-1)){
        val = val + sqrt(2)*coeffi[i+1]*w_find(aux.ite, t, upper, 2*t[ind] - i)
      }
      aux.table[ind] = val
      val = 0
    }
    aux.ite = aux.table
    aux.table = c()
  }
  return(aux.ite)
}


