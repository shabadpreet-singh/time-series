library(readr)
library(ggplot2)
library(dplyr)
library(forecast)
library(tseries)

# Loading the data and plotting time series
load("C:/Users/DD/Desktop/Time Series Project/Time Series.RData")
plot(dt, ylab = 'Average VWAP')



# Function to check Randomness

turning_point <-function(x)
{
  Q <- 0
  len <- length(x)
  
  for(i in c(2:(len-1)))
  {
    if(x[i] > max(x[i-1], x[i+1])||x[i] < min(x[i-1], x[i+1])){
      Q <- Q + 1
      }
  }
  
  exp_Q <- (len-2) * (2/3)
  var_Q <- (16*len - 29)/90
  
  test_stat <- (Q - exp_Q)/sqrt(var_Q)
  z_alpha <- qnorm(0.025)
  
  if(abs(test_stat) > abs(z_alpha)){
    cat("On the basis of the given data the series is not purely random")
  }
  else{
    cat("On the basis of the given data the series is purely random")
  }
}


# Function to check Trend
rel_ord <- function(data){
  len <- length(data)
  Q = 0
  
  for(i in 1 :(len-1)){
    Q = Q + sum(data[i] > data[(i+1):len])
  }
  
  tau = 1 - 4*Q/(len*(len-1))
  
  var.tau = 2*(2*len+5) / (9*len*(len-1))
  
  z = tau / sqrt(var.tau)
  
  alpha = 0.05
  
  crit.val = qnorm(1-alpha/2)
  
  cat(ifelse(abs(z) > crit.val, 'Trend is present', 'Trend is not present'))
}



##Function to check Seasonality

test_stat <- function(m,c,r){
  x <- 0
  den <- c*r*(r+1)
  num <- c*(r+1)*0.5
  for(i in 1:r){
    x <- x+((m[i]-num) ** 2)
  }
  return((12*x)/den) #den is denominator
}

seasonality_check <- function(y){
  dt_mtrx <- matrix(y,nrow<-12)
  dt_ranked_mtrx <- matrix(0,nrow<-12,ncol<-ncol(dt_mtrx))
  for(i in 1:ncol(dt_mtrx)){
    dt_ranked_mtrx[,i] <- rank(dt_mtrx[,i])
  }
  m_i <- rowSums(dt_ranked_mtrx)
  test_stat <-  test_stat(m_i,ncol(dt_ranked_mtrx),nrow(dt_ranked_mtrx))
  tab_Chi <- qchisq(0.95,nrow(dt_mtrx)-1)
  if(test_stat > tab_Chi){
    cat("We reject the null hypothesis therefore the data shows presence of seansonality")
  }
  else{
    cat("We fail to reject the null hypothesis, therefore may have seasonality")
  }
}


## Test for Randomness
turning_point(dt)
#from above we get to know that our data is not random
#i.e. their exists trend and seasonality


#Trend testing
rel_ord(dt)      #Trend is present

D_dt <- diff(dt) #detrending using differencing
plot(D_dt, ylab = 'detrended data')

rel_ord(D_dt)
#Trend is not present now

#Testing for Seasonality######
seasonality_check(D_dt)
#No seasonality


plot(D_dt, main = "Plot of Detrended and Deseasonalised Data", ylab = "Price") 

#Testing Stationarity#####

adf.test(dt)

adf.test(D_dt)
#Stationary


#Checking White noise by acf and pacf

acf(D_dt) # 1 significant spike suggesting value of q = 1
pacf(D_dt) # 2 significant spike suggesting value of p = 2

#Spikes suggest it is not a white noise process

D_dt <- ts(D_dt[1:239], start = c(2000, 1), frequency = 12)
AIC <- matrix(0, nrow = 5, ncol = 5)
minAIC <- Inf 
best_model <- c(0, 0)

for(i in 0:4)
{
  for(j in 0:4)
  {
    model <- arima(D_dt, c(i, 0, j), method = "ML")
    AIC[i+1, j+1] <- model$aic
    if(minAIC >=  model$aic)
    {
      minAIC = model$aic
      best_model = c(i, j)
    }
  }
}


train_dat <- ts(dt[1:239], start = c(2000, 1), frequency = 12)
test_dat <- ts(dt[240:252], start = c(2019, 12), frequency = 12)
final_model <- arima(train_dat, c(best_model[1], 1, best_model[2]), method = "ML")
min(AIC)

plot(D_dt, type = 'l')


model <- arma(D_dt, c(2, 3))
plot(D_dt, type = "l", ylab = 'det
     rended data')
lines(model$fitted.values, col = 2)
legend ("bottomleft", legend = c('Original', 'Predicted'), fill = c('black', 'red'))

plot(final_model$residuals)

# We see the residuals are scaled, so we scale them down


plot(final_model$residuals, ylab = "residuals")
acf(final_model$residuals)
pacf(final_model$residuals)

#acc to acf and pacf we get to know it is white noise i.e. uncorrelated,
#but in the main plt we can see that the residuals are scaled, therefore before making the qqplot we'll normalize them first

residual <- scale(final_model$residuals, center = T, scale = T)
qqnorm(residual, xlim = c(-2, 3), ylim = c(-2, 3))
abline(0, 1, col = 2)
checkresiduals(final_model) # take histogram from here

#We see that the residuals follow nomal distribution

Gfit <- ts(D_dt[127:139], start = c(2010, 6), frequency = 12)
predicted <- model$fitted.values[127:139]
t <- ts(predicted, start = c(2010, 6), frequency = 12)
kai2 <- data.frame(predicted, Gfit)
View(kai2)

subtract <- numeric(13)
for(i in 1:13)
{
  subtract[i] <- ((kai2[i, 1] - kai2[i, 2])**2)/kai2[i, 1]
}


plot(Gfit)
lines(t)

plot(forecast(final_model, h=15))

lines(test_dat, col = 2)
