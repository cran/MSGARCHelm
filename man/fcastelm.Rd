\name{fcastelm}
\alias{fcastelm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extreme Learning Machine Forecasting
}
\description{
The fcastelm function computes the volatility forecasting performance of Extreme Learning Machine (ELM) model with root mean square error (RMSE), mean absolute error (MAE), MAPE etc.
}
\usage{
fcastelm(data, stepahead=6, nlags=5, freq = frequency(data),
hn=10, est=c("lm"), rep=20, combt=c("mean"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Univariate time series data.
}
  \item{stepahead}{
The forecast horizon.
}
  \item{nlags}{
Lags of the data to use as inputs.
}
  \item{freq}{
Frequency of the time series.
}
  \item{hn}{
Number of hidden nodes.
}
  \item{est}{
Estimation type for output layer weights. Can be "lasso" (lasso with CV), "ridge" (ridge regression with CV), "step" (stepwise regression with AIC) or "lm" (linear regression). Default: est=c("lm").
}
  \item{rep}{
Number of networks to train, the result is the ensemble forecast.
}
  \item{combt}{
Combination operator for forecasts when rep > 1. Can be "median", "mode" (based on KDE estimation) and "mean". Default: combt=c("mean")
}
}
\details{
It helps to find the most appropriate Extreme Learning Machine model for the time series volatility forecasting.
}
\value{
$forecast_elm: 	    Forecasted value of Extreme Learning Machine.

$accuracy_elm:      Performance matrices of ELM model

}
\references{
Engle, R. (1982). Autoregressive conditional heteroscedasticity with estimates of the variance of United Kingdom inflation, Econometrica, 50, 987-1008.

Huang, G.B., Zhu Q.Y., and Siew, C.K. (2006). Extreme learning machine: Theory and applications. Neurocomputing, 70, 489-501.

}


\examples{
library(MSGARCHelm)
data(ReturnSeries_data)
fcastelm(ReturnSeries_data)
}
