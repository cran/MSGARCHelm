##### Function 'fcastelm' ############
requireNamespace("nnfor","MSGARCH", "forecast")
fcastelm <- function(data, stepahead=6, nlags=5, freq = frequency(data), hn=10, est=c("lm"), rep=20, combt=c("mean")){
  fit_arma <- forecast::auto.arima(data)
  resid_arma <- fit_arma$residuals
  data_elm <- (resid_arma)^2
  data_trn_elm <- ts(head(data_elm, round(length(data_elm) - stepahead)))
  data_test <- ts(tail(data_elm, stepahead))
  fit.elm <- nnfor::elm(data_trn_elm, m = freq, hd = hn, type = est, reps = rep, comb = combt,
                        lags = nlags, keep = NULL, difforder = NULL, outplot = c(
                          FALSE), sel.lag = c(FALSE), direct = c(FALSE),
                        allow.det.season = c(FALSE))
  fcast <- predict(fit.elm, h = stepahead)
  fcast_elm <- ts(fcast$mean)
  accuracy_elm <- forecast::accuracy(fcast_elm, data_test)
  All_result<-list(forecast_elm=fcast_elm,
                   accuracy_elm=accuracy_elm)
  return(All_result)
}
