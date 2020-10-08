\name{msgarchelm_OLS}
\alias{msgarchelm_OLS}

\title{
MS-GARCH-ELM combination based on OLS regession
}
\description{
This function computes combined time series volatility forecast of Markov Switching GARCH (MS-GARCH) and Extreme Learning Machine (ELM) model according to weights of ordinary least squares (OLS) regression.
}
\usage{
msgarchelm_OLS(data, stepahead = 10, nlags = 5, modelcomb =
c("sGARCH","gjrGARCH"), distcomb = c("norm", "std"), freq =
frequency(data), hn =10, est = c("lm"), rep = 20, combt =
c("mean"))
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
The OLS combination method (Granger and Ramanathan (1984)) uses ordinary least squares to estimate the weights for the combination. An appealing feature of the method is its bias correction through the intercept – even if one or more of the individual predictors are biased, the resulting combined forecast is unbiased. This function gives the combined volatility forecast of Markov Switching GARCH model and Extreme Learning machine model based on OLS approach.
}
\value{
$fcast_comb:	    Forecasted value of combined model according to ordinary least squares.

$accuracy_combined:      Performance matrices of the combined model.

}
\references{
Granger, C., and Ramanathan, R. (1984). Improved Methods Of Combining Forecasts. Journal of Forecasting, 3(2), 197–204.

Nowotarski, J., Raviv, E., Trueck, S., and Weron, R. (2014). An Empirical Comparison of Alternative Schemes for Combining Electricity Spot Price Forecasts. Energy Economics, 46, 395–412.

Timmermann, A. (2006). Forecast Combinations. In: Elliott, G., Granger, C. W. J., and Timmermann, A. (Eds.), Handbook of Economic Forecasting, 1, 135–196.

}



\examples{

library(MSGARCHelm)
data(ReturnSeries_data)
msgarchelm_OLS(ReturnSeries_data)

}
