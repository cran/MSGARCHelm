\name{msgarchelm_NG}
\alias{msgarchelm_NG}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Newbold and Granger MS-GARCH-ELM Combination
}
\description{
Computes time series volatility forecast combination of MS-GARCH and ELM model according to the approach by Newbold and Granger (1974) and evaluate its performance.
}
\usage{
msgarchelm_NG(data, stepahead = 10, nlags = 5, modelcomb =
c("sGARCH", "gjrGARCH"), distcomb = c("norm", "std"), freq =
frequency(data), hn = 10, est = c("lm"), rep = 20,
combt = c("mean"))
}

\arguments{
  \item{data}{
Univariate time series data.
}
  \item{stepahead}{
The forecast horizon.
}
  \item{nlags}{
Lags of the data to use as inputs in the Extreme Learning Machine (ELM).
}
  \item{modelcomb}{
Combination of volatility models in two different regimes in the MS-GARCH model. Valid models are "sARCH", "sGARCH", "eGARCH", "gjrGARCH", and "tGARCH". Default: modelcomb=c("sGARCH", "gjrGARCH").
}
  \item{distcomb}{
List with element distribution in the MS-GARCH model. distribution is a character vector (of size 2) of conditional distributions. Valid distributions are "norm", "snorm", "std", "sstd", "ged", and "sged". Default: distcomb=c("norm", "std").
}
  \item{freq}{
Frequency of the time series.
}
  \item{hn}{
Number of hidden nodes in the ELM model.
}
  \item{est}{
Estimation type for output layer weights in the ELM. Can be "lasso" (lasso with CV), "ridge" (ridge regression with CV), "step" (stepwise regression with AIC) or "lm" (linear regression). Default: est=c("lm").
}
  \item{rep}{
Number of networks to train, the result is the ensemble forecast in the ELM.
}
  \item{combt}{
Combination operator for forecasts in the ELM model when rep > 1. Can be "median", "mode" (based on KDE estimation) and "mean". Default: combt=c("mean").
}
}
\details{
It gives the combined volatility forecast of Markov Switching GARCH model and Extreme Learning machine model based on Newbold and Granger (1974) approach. Here MS-GARCH model is restricted to two regime. The Newbold and Granger (1974) approach extracts the combination weights from the estimated mean squared prediction error matrix.
}
\value{
$ fcast_comb:	    Forecasted value of combined model.
$accuracy_combined:      Performance matrices of the combined model.

}
\references{
Hsiao, C., and Wan, S. K. (2014). Is There An Optimal Forecast Combination? Journal of Econometrics, 178(2), 294–309.

Newbold, P., and Granger, C. W. J. (1974). Experience with Forecasting Univariate Time Series and the Combination of Forecasts. Journal of the Royal Statistical Society, Series A, 137(2), 131–165.

}


\examples{

library(MSGARCHelm)
data(ReturnSeries_data)
msgarchelm_NG(ReturnSeries_data)
}
