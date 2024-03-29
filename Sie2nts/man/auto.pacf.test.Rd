% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auto.pacf.test.R
\name{auto.pacf.test}
\alias{auto.pacf.test}
\title{The Test of Lag of Auto-Regressive (AR) Model Automatically}
\usage{
auto.pacf.test(
  ts,
  lag = 3,
  b = 8,
  or = 4,
  type,
  alpha = 0.05,
  method = "LOOCV",
  threshold = 0,
  B.s = 1000
)
}
\arguments{
\item{ts}{ts is the data set which is a time series data typically}

\item{lag}{the lag for auto-regressive model, the default value is 3}

\item{b}{the largest lag for auto-regressive model, the default value is 8, this parameter must be
larger than lag}

\item{or}{or indicates the order of spline and only used in Cspli type, default is 4 which indicates cubic spline}

\item{type}{type indicates which type of basis is used. There are 31 types in this package}

\item{alpha}{level of the test}

\item{method}{method indicates which method used to choose optimal parameters, 3 methods in this package can be used.}

\item{threshold}{threshold determines the bound for Elbow method}

\item{B.s}{the number of statistics used in multiplier bootstrap, the default value is 1000}
}
\value{
p value of the test
}
\description{
auto.pacf.test() generates a test of lag of AR model by choosing tuning parameter automatically.
}
\details{
In the parameter type, this package provides 32 types of basis including options "Legen" for Legendre polynomials, "Cheby" for first kind
Chebyshev polynomials, "tri" for trigonometric polynomials, "cos" for cosine polynomials, "sin" for sine polynomials, "Cspli" for splines means Class of splines
functions, in this option, the first input "c" is knots plus 2 that represent 0 and 1. "or" indicates the order of splines, so the number of basis is number of knots + 2 - 2 plus
the number of order.When functions automatically choose the number of basis for splines, the number is not less than the order of spline. "db1" to "db20" for Daubechies1 wavelet basis to Daubechies20
wavelet basis and "cf1" to "cf5" for Coiflet1 wavelet basis to Coiflet5 wavelet basis. The package provides the wavelet tables are generated by Cascade algorithm using low-pass filter.
If the exact values of wavelet are required, Recursion algorithm should be used.
In the parameter method, it contains 3 options, the default option is "LOOCV", it uses Leave-One-Out Cross-Validation to choose the best tuning parameters.
The second choice is "CV" which uses the Cross-Validation method, it takes 3*log2(n) size as validation set where n is the number of total observations.
The third choice is "Elbow".This method similar as "lOOCV", however, it set the threshold manually. The function will choose the smallest tuning parameters
once the value of LOOCV is less than threshold.
}
\references{
[1] Ding, Xiucai, and Zhou, Zhou. “Estimation and inference for precision matrices of nonstationary time series.” The Annals of Statistics 48(4) (2020): 2455-2477.

[2] Ding, Xiucai, and Zhou, Zhou. “Auto-regressive approximations to non-stationary time series, with inference and applications.” Available online, 2021.
}
