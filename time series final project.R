library(readxl)
library(forecast)
library(tibble)
library(dplyr)
library(ggplot2)
library(tseries)

file_path <- "/Users/msyzdykova/Downloads/master_data.xlsx"

df <- read_excel(file_path)
str(df)

series <- ts(df$MSACSR, start = c(1971, 4),frequency = 12)  # Monthly data assumed

plot(series,
     main = "Monthly Supply of New Houses in the U.S. (Apr 1971-Feb 2025)",
     ylab = "Months' Supply",
     xlab = "Year",
     col = "steelblue",
     lwd = 2)

#----------------------------ACF------------------------------------
acf(ts(series), main = "ACF Plot for Series", lag.max = 40)

#---------------------FIT A LINEAR TREND-----------------------------
t <- 1:length(series)

#MSACSR ~ t
linear_model <- lm(series ~ t)
summary(linear_model)

#fitted values and residuals
fitted_vals <- fitted(linear_model)
residuals_lin <- residuals(linear_model)
adf.test(residuals_lin)

#original series + linear trend
plot(series, 
     main = "Monthly Supply of New Houses with Linear Trend",
     ylab = "Months' Supply",
     xlab = "Year",
     col = "steelblue",
     lwd = 2)
lines(ts(fitted_vals, start = c(1971, 4), frequency = 12), 
      col = "red", 
      lwd = 2)
legend("topright", 
       legend = c("Original Series", "Linear Trend"),
       col = c("steelblue", "red"),
       lwd = 2)

#residuals from linear trend model
plot(residuals_lin, type = "l", col = "black", lwd = 1.5,
     main = "Residuals from Linear Trend Model",
     ylab = "Residuals", xlab = "Time")
abline(h = 0, col = "blue", lty = 2)

#---------------------STATIONARITY-----------------------------
#ADF test on Y
adf.test(series)

#Difference Y to make it stationary
series_diff<-diff(series)
adf.test(series_diff)

#---------------------Modeling CYCLICALITY----------------------
#AIC function
calculate_ar_aic <- function(ts_data, max_p = 12) {
  results <- tibble(p = integer(), AIC = numeric())
  
  for (p in 1:max_p) {
    fit <- Arima(ts_data, order = c(p, 1, 0))  # AR(p)
    results <- add_row(results, p = p, AIC = AIC(fit))
  }
  
  return(results)
}

aic_table <- calculate_ar_aic(series, max_p = 12)
print(aic_table)

#Fit an AR(3) model
ar2_model <- Arima(series, order = c(2, 1, 0))  # AR(2)
summary(ar2_model)

#Fitted values
fitted_ar2 <- fitted(ar2_model)
#fitted_ar5 <- fitted(ar5_model)

plot(series, 
     main = "AR(2) Model: Original vs Fitted", 
     ylab = "MSACSR", 
     xlab = "Time", 
     col = "red", 
     lwd = 1.5)
lines(fitted_ar2, col = "blue", lwd = 2)
legend("topright", legend = c("Original", "Fitted"), 
       col = c("red", "blue"), lwd = 2)

residuals_ar2 <- residuals(ar2_model)
plot(residuals_ar2, type = "l", col = "black", lwd = 1.5,
     main = "Residuals from AR(2) Model",
     ylab = "Residuals", xlab = "Time")
abline(h = 0, col = "red", lty = 2)


#Fit an ARMA model
#no drift
auto_arma_model <- auto.arima(series, d=1,seasonal = FALSE, stepwise = FALSE, approximation = FALSE, ic=c("aic"))
#BEST AIC=1048.28 (but not a big improvement from AR(3))
summary(auto_arma_model)

fitted_arma <- fitted(auto_arma_model)
plot(series, 
     main = "ARMA(2,1,3) Model: Original vs Fitted", 
     ylab = "MSACSR", 
     xlab = "Time", 
     col = "red", 
     lwd = 1.5)
lines(fitted_arma, col = "blue", lwd = 2)
legend("topright", legend = c("Original", "Fitted"), 
       col = c("red", "blue"), lwd = 2)

#time plot of residuals
residuals_arma <- residuals(auto_arma_model)

plot(residuals_arma, 
     main = "Residuals of ARMA(2,1,3) Model", 
     ylab = "Residuals", 
     xlab = "Time")
abline(h = 0, col = "red", lty = 2)


#---------------------INDEPENDENT VARIABLES-----------------------------

#ADF for stationarity
vars <- c("MSACSR","mortgage_rates", "hous_st", "ppi", "hous_pr", "disp_inc")

run_adf <- function(series) {
  result <- adf.test(series)
  return(c("ADF Statistic" = result$statistic, "p-value" = result$p.value))
}

adf_results <- sapply(df[vars], run_adf)

#DIFFERENCE ALL OF THE X's

#mortgage_rates_diff <- diff(df$mortgage_rates)
#unemp_rate_diff     <- diff(df$unemp_rate)
mortgage_rate_diff <- diff(df$mortgage_rates)
hous_st_diff        <- diff(df$hous_st)
ppi_diff            <- diff(df$ppi)
hous_pr_diff        <- diff(df$hous_pr)
disp_inc_diff       <- diff(df$disp_inc)


diff_vars <- list(mortgage_rate_diff,hous_st_diff,
                  ppi_diff, hous_pr_diff, disp_inc_diff 
)

adf_results <- sapply(diff_vars, function(x) {
  test <- adf.test(x)
  c(ADF_Statistic = test$statistic, p_value = test$p.value)
})

#---------------------GRANGER CAUSALITY-----------------------------
library(lmtest)
library(sandwich)
library(car)

lags <- 12 

predictors <- c("mortgage_rate_diff", "hous_st_diff",
                "ppi_diff", "hous_pr_diff", "disp_inc_diff")

wald_results <- list()

for (xvar in predictors) {
  df_temp <- data.frame(
    y = series_diff[(lags + 1):length(series_diff)]
  )

  for (i in 1:lags) {
    df_temp[[paste0("y_lag", i)]] <- series_diff[(lags + 1 - i):(length(series_diff) - i)]
  }
 
  x_data <- get(xvar)
  for (i in 1:lags) {
    df_temp[[paste0("x_lag", i)]] <- x_data[(lags + 1 - i):(length(x_data) - i)]
  }

  model <- lm(y ~ ., data = df_temp)
 
  hyp_string <- paste0("x_lag", 1:lags, " = 0")

  vcov_hc <- vcovHC(model, type = "HC1")
  test_result <- linearHypothesis(model, hyp_string, vcov = vcov_hc)
  
  wald_results[[xvar]] <- test_result
}

wald_results

#Keep only these - mortgage rates, ppi, disp_inc

#---------------------LAG SELECTION FOR ADL----------------------------

mortgage_rates_diff<-ts(mortgage_rate_diff)
ppi_diff<-ts(ppi_diff)
disp_inc_diff<-ts(disp_inc_diff)

Y <- series_diff
X_list <- list(mortgage_rates_diff, ppi_diff, disp_inc_diff)
X_names <- c("mortgage_rates_diff", "ppi_diff", "disp_inc_diff")

max_lag_y <- 12
max_lag_x <- 12

results <- data.frame()

#Grid search
for (p in 1:max_lag_y) {
  for (q1 in 0:max_lag_x) {
    for (q2 in 0:max_lag_x) {
      for (q3 in 0:max_lag_x) {
 
        q_vec <- c(q1, q2, q3)
        max_lag <- max(p, q_vec)
        n <- length(Y)
  
        start_index <- max_lag + 1
        end_index <- n
        
        df <- data.frame(y = Y[start_index:end_index])

        if (p > 0) {
          for (i in 1:p) {
            df[[paste0("y_lag", i)]] <- Y[(start_index - i):(end_index - i)]
          }
        }

        for (j in 1:3) {
          X <- X_list[[j]]
          q <- q_vec[j]
          if (q > 0) {
            for (k in 1:q) {
              df[[paste0(X_names[j], "_lag", k)]] <- X[(start_index - k):(end_index - k)]
            }
          }
        }

        model <- lm(y ~ ., data = df)
        model_aic <- AIC(model)

        results <- rbind(results, data.frame(
          p = p, q1 = q1, q2 = q2, q3 = q3, AIC = model_aic
        ))
      }
    }
  }
}

#best model
best_model <- results[which.min(results$AIC), ]
print(best_model)


#Run the best model
best_p <- best_model$p
best_qs <- c(best_model$q1, best_model$q2, best_model$q3)
max_lag_used <- max(best_p, best_qs)
start_index <- max_lag_used + 1
end_index <- length(Y)

df_best <- data.frame(y = Y[start_index:end_index])

for (i in 1:best_p) {
  df_best[[paste0("y_lag", i)]] <- Y[(start_index - i):(end_index - i)]
}

for (j in 1:3) {
  X <- X_list[[j]]
  q <- best_qs[j]
  
  if (q > 0) {
    for (k in 1:q) {
      df_best[[paste0(X_names[j], "_lag", k)]] <- X[(start_index - k):(end_index - k)]
    }
  }
}

final_model <- lm(y ~ ., data = df_best)

summary(final_model)
AIC(final_model)

res_fin_model <- residuals(final_model)  # actual - predicted

res_ts <- ts(res_fin_model, start = c(1971, 4), frequency = 12)

plot(res_ts, 
     main = "Residuals of Final Model", 
     ylab = "Residuals", 
     xlab = "Time")
abline(h = 0, col = "red", lty = 2)


#---------------------FORECAST USING ADL-----------------------------
# Lags as per your AIC results
p <- 4   # y lags
q1 <- 12  # mortgage_rates lags
q2 <- 8 # ppi lags

X_list <- list(mortgage_rates, ppi)
X_names <- c("mortgage_rates", "ppi")
q_vec <- c(q1, q2)

unrate_f <- rep(NA, 12)
unrate_pi <- matrix(NA, nrow = 12, ncol = 2)

for (h in 1:12) {
  total_lags <- max(p, q1, q2) + h
  
  y_mat <- embed(Y, total_lags)
  y_response <- y_mat[, 1]                      
  y_lags <- y_mat[, (h + 1):(h + p)]           

  X_lags <- list()
  for (j in 1:2) {
    Xj <- X_list[[j]]
    qj <- q_vec[j]
    x_mat <- embed(Xj, total_lags)
    X_lags[[j]] <- x_mat[, (h + 1):(h + qj)]     
  }
  
  df_res <- data.frame(y = y_response, y_lags)
  colnames(df_res)[2:(1 + p)] <- paste0("y_lag", 1:p)
  
  for (j in 1:2) {
    x_df <- as.data.frame(X_lags[[j]])
    colnames(x_df) <- paste0(X_names[j], "_lag", 1:q_vec[j])
    df_res <- cbind(df_res, x_df)
  }

  model_h <- lm(y ~ ., data = df_res)
 
  y_latest <- rev(tail(Y, p))
  x_latest <- unlist(lapply(1:2, function(j) rev(tail(X_list[[j]], q_vec[j]))))
  new_data <- as.data.frame(t(c(y_latest, x_latest)))
  colnames(new_data) <- names(coef(model_h))[-1]  
  
  pred <- predict(model_h, newdata = new_data, interval = "prediction", level = 0.90)
  
  unrate_f[h] <- pred[1]
  unrate_pi[h, ] <- pred[2:3]
}

forecast_result <- data.frame(
  horizon = 1:12,
  forecast = unrate_f,
  lower_90 = unrate_pi[, 1],
  upper_90 = unrate_pi[, 2]
)

print(forecast_result)

#back to lvls
last_value <- tail(series, 1)
level_forecast <- cumsum(c(last_value, forecast_result$forecast))[-1]

level_lower <- cumsum(c(last_value, forecast_result$lower_90))[-1]
level_upper <- cumsum(c(last_value, forecast_result$upper_90))[-1]

forecast_result$level_forecast <- level_forecast
forecast_result$level_lower_90 <- level_lower
forecast_result$level_upper_90 <- level_upper


#plot the forecasts 
freq_Y <- frequency(series)
end_time <- time(series)[length(series)]

time_forecast <- seq(from = end_time + 1/freq_Y, by = 1/freq_Y, length.out = 12)

df_plot <- data.frame(
  time = c(time(series), time_forecast),
  value = c(as.numeric(series), rep(NA, 12)),
  type = "Actual"
)

df_forecast <- data.frame(
  time = time_forecast,
  value = forecast_result$level_forecast,
  lower = forecast_result$level_lower,
  upper = forecast_result$level_upper,
  type = "Forecast"
)

ggplot() +
  geom_line(data = df_plot, aes(x = time, y = value), color = "black", size = 1) +
  geom_line(data = df_forecast, aes(x = time, y = value), color = "blue", size = 1) +
  geom_ribbon(data = df_forecast, aes(x = time, ymin = lower, ymax = upper), alpha = 0.2, fill = "red") +
  labs(title = "Forecast of Y with 90% Prediction Interval",
       x = "Time", y = "Y") +
  theme_minimal()


ggplot() +
  geom_line(data = subset(df_plot, time >= 2020), aes(x = time, y = value), color = "black", size = 1) +
  geom_line(data = subset(df_forecast, time >= 2020), aes(x = time, y = value), color = "blue", size = 1) +
  geom_ribbon(data = subset(df_forecast, time >= 2020), 
              aes(x = time, ymin = lower, ymax = upper), 
              alpha = 0.2, fill = "red") +
  labs(title = "Forecast of monthly supply of new houses in the U.S. (90% prediction interval, from 2020)",
       x = "Time", y = "Monthly supply of new houses") +
  theme_minimal()


#---------------------FORECASTS WITH GARCH-----------------------------
#Shapiro test for the normality of residuals
shapiro.test(as.numeric(res_fin_model))

hist(res_fin_model, breaks = 30, probability=TRUE,main = "Histogram of Residuals from ADL", col = "skyblue", xlab = "Residuals")
curve(dnorm(x, mean = mean(res_fin_model), sd = sd(res_fin_model)),
      col = "red", lwd = 2, add = TRUE)
legend("topright", legend = "Normal Distribution", col = "red", lwd = 2)

#GARCH specification
library(rugarch)

garch_spec <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
  distribution.model = "ged"
)

res_index <- seq(as.Date("1971-04-01"), by = "month", length.out = length(res_fin_model))

res_xts <- xts::xts(res_fin_model, order.by = res_index)

garch_fit <- ugarchfit(spec = garch_spec, data = res_xts)

print(garch_fit)

cond_volatility <- sigma(garch_fit)
std_residuals <- res_xts / cond_volatility

par(mfrow = c(2, 2))
plot(std_residuals, type = "l", main = "Standardized Residuals")
acf(std_residuals, main = "ACF of Standardized Residuals")
acf(std_residuals^2, main = "ACF of Squared Standardized Residuals")
qqnorm(std_residuals); qqline(std_residuals)
par(mfrow = c(1, 1))

weights <- 1/cond_volatility^2

X <- model.matrix(final_model)
y <- final_model$model$y

wls_model <- lm(y ~ X - 1, weights = weights)
summary(wls_model)
AIC(wls_model)

garch_forecast <- ugarchforecast(garch_fit, n.ahead = 12)
sigma_forecast <- sigma(garch_forecast)
print(sigma_forecast)

sigma_forecast <- sigma(garch_forecast)[1:12]  

z_val <- qnorm(0.95)  

point_forecast_diff <- forecast_result$forecast

Y_last <- tail(series, 1)
forecast_levels <- cumsum(c(Y_last, point_forecast_diff))[-1] 

lower_bound <- forecast_levels - z_val * sigma_forecast
upper_bound <- forecast_levels + z_val * sigma_forecast

df_combined <- data.frame(
  Horizon = 1:12,
  Forecast = forecast_levels,
  Lower90 = lower_bound,
  Upper90 = upper_bound
)

print(df_combined)

#PLOT
library(lubridate)

Y_ts <- ts(series, start = c(1971, 4), frequency = 12)  # adjust start date as needed
Y_length <- length(series)
Y_time <- time(Y_ts)

forecast_dates <- seq(as.Date(paste(2025, 03, "01", sep = "-")), by = "month", length.out = 12)

df_plot <- data.frame(
  Date = c(as.Date(as.yearmon(Y_time)), forecast_dates),
  Value = c(as.numeric(series), df_combined$Forecast),
  Type = c(rep("Actual", length(series)), rep("Forecast", 12)),
  Lower = c(rep(NA, length(series)), df_combined$Lower90),
  Upper = c(rep(NA, length(series)), df_combined$Upper90)
)

ggplot(df_plot, aes(x = Date)) +
  geom_line(aes(y = Value, color = Type), size = 1) +
  geom_ribbon(data = subset(df_plot, Type == "Forecast"),
              aes(ymin = Lower, ymax = Upper, x = Date),
              fill = "grey70", alpha = 0.4) +
  labs(
    title = "Forecast of Monthly Supply of New Houses with GARCH Intervals",
    x = "Date", y = "Monthly Supply"
  ) +
  scale_color_manual(values = c("Actual" = "black", "Forecast" = "blue")) +
  theme_minimal() +
  theme(legend.position = "bottom")

#from 2020
df_plot_recent <- subset(df_plot, Date >= as.Date("2020-01-01"))

ggplot(df_plot_recent, aes(x = Date)) +
  geom_line(aes(y = Value, color = Type), size = 1) +
  geom_ribbon(data = subset(df_plot_recent, Type == "Forecast"),
              aes(x = Date, ymin = Lower, ymax = Upper),
              fill = "grey70", alpha = 0.4) +
  labs(
    title = "Forecast of Monthly Supply of New Houses (2020–Present)",
    x = "Date", y = "Monthly Supply"
  ) +
  scale_color_manual(values = c("Actual" = "black", "Forecast" = "blue")) +
  theme_minimal() +
  theme(legend.position = "bottom")

#---------------------OOS PERFORMANCE EVALUATION-----------------------------

Y_last <- tail(series_train, 1)
forecast_levels <- cumsum(c(Y_last, unrate_f))[-1]  # point forecast in levels

file_path <- "/Users/msyzdykova/Downloads/master_data.xlsx"
df <- read_excel(file_path)

series <- ts(df$MSACSR, start = c(1971, 4), frequency = 12)
mortgage_rates <- df$mortgage_rates
ppi <- df$ppi

n_test <- 12
series_train <- head(series, length(series) - n_test)
series_test <- tail(series, n_test)

mortgage_train <- head(mortgage_rates, length(series) - n_test)
ppi_train <- head(ppi, length(series) - n_test)

Y <- diff(series_train)
mortgage_diff <- diff(mortgage_train)
ppi_diff <- diff(ppi_train)

p <- 4   # y lags
q1 <- 12 # mortgage_rates lags
q2 <- 8  # ppi lags

unrate_f <- rep(NA, n_test)
unrate_pi <- matrix(NA, nrow = n_test, ncol = 2)

for (h in 1:n_test) {
  total_lags <- max(p, q1, q2) + h
  
  y_mat <- embed(Y, total_lags)
  y_response <- y_mat[, 1]
  y_lags <- y_mat[, (h + 1):(h + p)]
  
  X_list <- list(mortgage_diff, ppi_diff)
  X_names <- c("mortgage_rates", "ppi")
  q_vec <- c(q1, q2)
  
  X_lags <- list()
  for (j in 1:2) {
    Xj <- X_list[[j]]
    qj <- q_vec[j]
    x_mat <- embed(Xj, total_lags)
    X_lags[[j]] <- x_mat[, (h + 1):(h + qj)]
  }

  df_train <- data.frame(y = y_response, y_lags)
  colnames(df_train)[2:(1 + p)] <- paste0("y_lag", 1:p)
  for (j in 1:2) {
    x_df <- as.data.frame(X_lags[[j]])
    colnames(x_df) <- paste0(X_names[j], "_lag", 1:q_vec[j])
    df_train <- cbind(df_train, x_df)
  }
  
  model_h <- lm(y ~ ., data = df_train)
  
  y_latest <- rev(tail(Y, p))
  x_latest <- unlist(lapply(1:2, function(j) rev(tail(X_list[[j]], q_vec[j]))))
  new_data <- as.data.frame(t(c(y_latest, x_latest)))
  colnames(new_data) <- names(coef(model_h))[-1]
  
  pred <- predict(model_h, newdata = new_data, interval = "prediction", level = 0.90)
  
  unrate_f[h] <- pred[1]
  unrate_pi[h, ] <- pred[2:3]
}


z_val <- qnorm(0.95)
sigma_cum <- sqrt(cumsum(sigma_forecast^2))  # cumulative std devs
lower_bound <- forecast_levels - z_val * sigma_cum
upper_bound <- forecast_levels + z_val * sigma_cum

forecast_errors <- series_test - forecast_levels

rmse <- sqrt(mean(forecast_errors^2))
mae <- mean(abs(forecast_errors))
mape <- mean(abs(forecast_errors / series_test)) * 100
coverage <- mean(series_test >= lower_bound & series_test <= upper_bound) * 100

cat("RMSE:", rmse, "\n")
cat("MAE:", mae, "\n")
cat("MAPE:", mape, "%\n")
cat("Prediction Interval Coverage Rate:", coverage, "%\n")

forecast_dates <- time(tail(series, n_test))
eval_df <- data.frame(
  Date = as.Date(as.yearmon(forecast_dates)),
  Actual = as.numeric(series_test),
  Forecast = as.numeric(forecast_levels),
  Lower = lower_bound,
  Upper = upper_bound
)

ggplot(eval_df, aes(x = Date)) +
  geom_line(aes(y = Actual, color = "Actual"), size = 1) +
  geom_line(aes(y = Forecast, color = "Forecast"), size = 1) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "grey70", alpha = 0.4) +
  labs(
    title = "Out-of-Sample Forecast vs Actual",
    y = "Monthly Supply",
    x = "Date"
  ) +
  scale_color_manual(values = c("Actual" = "black", "Forecast" = "blue")) +
  theme_minimal() +
  theme(legend.position = "bottom")

