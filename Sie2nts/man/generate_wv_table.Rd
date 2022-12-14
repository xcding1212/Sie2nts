% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_wv_table.R
\name{generate_wv_table}
\alias{generate_wv_table}
\title{Generate wavelet table}
\usage{
generate_wv_table(coeffi, len, ite)
}
\arguments{
\item{coeffi}{It contains low-pass filter coefficients to generate wavelet table, the coefficients can get from the website}

\item{len}{The total number of points}

\item{ite}{The number of iterations to generate table, more numbers of iterations means more accuracy}
}
\value{
A vector contains the value of scaling wavelet functions
}
\description{
generate_wv_table() generates wavelet table by using Cascade algorithm
}
\examples{
This example generate Daubechies 2 wavelet scaling function
coef = c(0.48296291314469025, 0.836516303737469, 0.22414386804185735, -0.12940952255092145)
tb = generate_wv_table(coef, 31680, 5)
plot(tb)
}
\references{
[1] Wasilewski, F. (n.d.). Wavelet browser by pywavelets. Retrieved November 18, 2022, from
http://wavelets.pybytes.com/
}
