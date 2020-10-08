
\name{msgarchelm_BG}
\alias{msgarchelm_BG}

\title{
Bates and Granger MS-GARCH-ELM Combination
}
\description{
This function computes combined time series volatility forecast of Markov Switching GARCH (MS-GARCH) and Extreme Learning Machine (ELM) model according to the approach of Bates and Granger (1969).
}
\usage{
msgarchelm_BG(data, stepahead=10, nlags=3, modelcomb=c("sGARCH",
"gjrGARCH"), distcomb=c("norm", "std"), freq = frequency(data),
hn=10, est=c("lm"), rep=20, combt=c("mean"))

}
\arguments{
  \item{data}{
Univariate time series data.
}
  \item{stepahead}{
The forecast horizon
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
Bates and Granger (1969) introduce the idea of combining forecasts. Their approach builds on portfolio diversification theory and uses the diagonal elements of the estimated mean squared prediction error matrix in order to compute combination weights. This function gives the combined volatility forecast of Markov Switching GARCH model and Extreme Learning machine model based on Bates and Granger (1969) approach.
}
\value{
$ fcast_comb:	    Forecasted value of combined model according to Bates and Granger (1969).

$accuracy_combined:      Performance matrices of the combined model.

}
\references{
Bates, J. M., and Granger, C. W. (1969). The Combination of Forecasts. Journal of the Operational Research Society, 20(4), 451–468.

Hsiao, C., and Wan, S. K. (2014). Is There An Optimal Forecast Combination? Journal of Econometrics, 178(2), 294–309.

Timmermann, A. (2006). Forecast Combinations. In: Elliott, G., Granger, C. W. J., and Timmermann, A. (Eds.), Handbook of Economic Forecasting, 1, 135–196.

}

\examples{

library(MSGARCHelm)
data(ReturnSeries_data)
msgarchelm_BG(ReturnSeries_data)
}
