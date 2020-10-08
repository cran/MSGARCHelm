##### Function 'msgarchelm_NG' ###################
requireNamespace("nnfor","MSGARCH", "forecast")
msgarchelm_NG <- function(data, stepahead=10, nlags=5, modelcomb=c("sGARCH", "gjrGARCH"), distcomb=c("norm", "std"), freq = frequency(data), hn=10, est=c("lm"), rep=20, combt=c("mean")){
  fit_arma <- forecast::auto.arima(data)
  resid_arma <- fit_arma$residuals
  data_elm <- (resid_arma)^2
  data_msgarch <- resid_arma
  data_trn_msgarch <- ts(head(data_msgarch, round(length(data_msgarch) - stepahead)))
  data_trn_elm <- ts(head(data_elm, round(length(data_elm) - stepahead)))
  data_test <- ts(tail(data_elm, stepahead))
  fit.elm <- nnfor::elm(data_trn_elm, m = freq, hd = hn, type = est, reps = rep, comb = combt,
                        lags = nlags, keep = NULL, difforder = NULL, outplot = c(
                          FALSE), sel.lag = c(FALSE), direct = c(FALSE),
                        allow.det.season = c(FALSE))
  fcast <- predict(fit.elm, h = stepahead)
  fcast_elm <- ts(fcast$mean)
  msgarch_spec <- MSGARCH::CreateSpec(variance.spec = list(model = modelcomb),
                                      distribution.spec = list(distribution = distcomb),
                                      switch.spec = list(do.mix = FALSE, K = NULL),
                                      constraint.spec = list(fixed = list(), regime.const = NULL),
                                      prior = list(mean = list(), sd = list()))
  fit.ml <- MSGARCH::FitML(spec = msgarch_spec, data = data_trn_msgarch)
  pred <- predict(fit.ml, nahead = stepahead, do.return.draw = TRUE)
  volatality_ms<- pred$vol
  volatality<- ts(volatality_ms)
  volatality_sqr<- (volatality)^2
  fcast_msgarch<- ts(volatality_sqr)
  fcast_matrix <- cbind(fcast_elm, fcast_msgarch)
  error_matrix <- data_test - fcast_matrix
  sample_msqu_pred_error <- (t(error_matrix) %*% error_matrix)/length(data_test)
  e_vec <- as.matrix(rep(1, ncol(fcast_matrix)))
  weights <- as.vector(solve(sample_msqu_pred_error) %*% e_vec)/as.numeric(t(e_vec) %*%
                                                                             solve(sample_msqu_pred_error) %*% e_vec)
  fitted_combined <- as.vector(fcast_matrix %*% weights)
  accuracy_combined <- forecast::accuracy(fitted_combined, data_test)
  All_result<-list(fcast_comb=fitted_combined,
                   accuracy_combined=accuracy_combined)
  return(All_result)
}
