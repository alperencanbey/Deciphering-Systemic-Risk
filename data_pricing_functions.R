#This code consists the base functions which are necessary to create the stochastic intensity process,
#simulated CDS prices, and the derivatives of the CDS pricing function wrt the parameters.


#lambda demonstrates the intensity process which follows Cox–Ingersoll–Ross model.
#beta is the speed of the mean reversion, sigma is the volatility, alpha/beta is the mean of the process.
lambda_process <- function(true_parameters,lambda0,days_in_year,year,sample_count){
  
  # --- We initilize values an vectors --- #
  N = days_in_year*year
  T = year
  alpha = true_parameters[1]
  beta = true_parameters[2]
  sigma_sq = true_parameters[3]
  
  lambda_sde = rep(NA,N)
  lambda_sde[1] = lambda0

#the interior function creates the stochastic process through the stochastic differential equation of lambda
  simulation <- function(lambda_sde){
    for (i in 2:N){
      
      delta_b <- rnorm(1,0,sqrt(T/N))
      
      lambda_sde[i] = lambda_sde[i-1] + sqrt(sigma_sq)*sqrt(abs(lambda_sde[i-1]))*sign(lambda_sde[i-1])*(delta_b) + 
        (alpha-beta*lambda_sde[i-1])*(T/N) 
    }
    return(lambda_sde)
  }
  
  
  matrix_sde = matrix(NA, sample_count, N)
  for (i in 1:sample_count) {
    vector_sde <- simulation(lambda_sde)
    matrix_sde[i,1:N] = vector_sde
  }
  
  mean_lambda_sde = colMeans(matrix_sde)
  
  # --- We return our trajectory  --- #
  
  return(mean_lambda_sde)
}


#This is cds pricing formula for the systemic part, derived from Ang&Lonstaff (2013) Appendix
CDS_price <- function(parameters, lambda, days_in_year, maturity){
  #N = T*36
  
  N = days_in_year*maturity
  T = maturity
  alpha = parameters[1]
  beta = parameters[2]
  sigma_sq = parameters[3]
  
  t= seq(0, T-T/N, T/N)
  gamma = 1
  
  x= sqrt(beta^2 + 2*gamma*sigma_sq)
  y= (beta+x)/(beta-x)
  
  A1 = exp(alpha*(beta+x)*t/sigma_sq)*((1-y)/(1-y*exp(x*t)))^(2*alpha/sigma_sq)
  A2 = (beta-x)/sigma_sq + 2*x/(sigma_sq*(1-y*exp(x*t)))
  
  F1 = alpha/x*(exp(x*t)-1)*exp(alpha*(beta+x)*t/sigma_sq)*((1-y)/(1-y*exp(x*t)))^(2*alpha/sigma_sq+1)
  F2 = exp(alpha*(beta+x)*t/sigma_sq+x*t)*((1-y)/(1-y*exp(x*t)))^(2*alpha/sigma_sq+2)
  
  A = A1*exp(A2*lambda)
  F = (F1+F2*lambda)*exp(A2*lambda)
  
  r = seq(0, 0.03*T-0.03*T/N, by=0.03*T/N)
  si = gamma*sum(exp(-r)*F)/sum(exp(-r)*A)
}


#Data creation process depending on parameters and initial lambda.
#First it creates the daily intensity values and use them to create CDS data.
data_creation <- function(true_parameters, lambda0, days_in_year, year, sample_count_BM, random_error){
  

  daily_intensity = lambda_process(true_parameters,lambda0,days_in_year,year,sample_count_BM)
  
  cds_data = matrix(NA, days_in_year*year, 5)
  cds_true_values = matrix(NA, days_in_year*year, 5)
  for (i in 1:length(daily_intensity)) {
    for (T in 1:5) {
      lambda = daily_intensity[i]
      cds_true_values[i,T] = CDS_price(true_parameters,lambda,days_in_year,T)
      cds_data[i,T] = cds_true_values[i,T] + rnorm(1,0,random_error*cds_true_values[i,T])
    }
  }
  
  data <- list("cds_data" = cds_data, "daily_intensity" = daily_intensity, "cds_true_values" = cds_true_values)
  return(data)
}



#This function calculates the difference between the estimated CDS spread and observed (simulated in this case) CDS spread
CDS_price_residual <- function(parameters, lambda, days_in_year, maturity, data_point){
  
  N = days_in_year*maturity
  T = maturity
  alpha = parameters[1]
  beta = parameters[2]
  sigma_sq = parameters[3]
  gamma = 1
  t= seq(0, T-T/N, T/N)
  
  x= sqrt(beta^2 + 2*gamma*sigma_sq)
  y= (beta+x)/(beta-x)
  
  A1 = exp(alpha*(beta+x)*t/sigma_sq)*((1-y)/(1-y*exp(x*t)))^(2*alpha/sigma_sq)
  A2 = (beta-x)/sigma_sq + 2*x/(sigma_sq*(1-y*exp(x*t)))
  
  F1 = alpha/x*(exp(x*t)-1)*exp(alpha*(beta+x)*t/sigma_sq)*((1-y)/(1-y*exp(x*t)))^(2*alpha/sigma_sq+1)
  F2 = exp(alpha*(beta+x)*t/sigma_sq+x*t)*((1-y)/(1-y*exp(x*t)))^(2*alpha/sigma_sq+2)
  
  A = A1*exp(A2*lambda)
  F = (F1+F2*lambda)*exp(A2*lambda)
  
  r = seq(0, 0.03*T-0.03*T/N, by=0.03*T/N)
  si = gamma*sum(exp(-r)*F)/sum(exp(-r)*A)
  
  diff = (si-data_point)
}




#Same as above but squared difference this time.
CDS_price_error <- function(parameters, lambda, days_in_year, maturity, data_point){
  
  N = days_in_year*maturity
  T = maturity
  alpha = parameters[1]
  beta = parameters[2]
  sigma_sq = parameters[3]
  gamma = 1
  t= seq(0, T-T/N, T/N)
  
  x= sqrt(beta^2 + 2*gamma*sigma_sq)
  y= (beta+x)/(beta-x)
  
  A1 = exp(alpha*(beta+x)*t/sigma_sq)*((1-y)/(1-y*exp(x*t)))^(2*alpha/sigma_sq)
  A2 = (beta-x)/sigma_sq + 2*x/(sigma_sq*(1-y*exp(x*t)))
  
  F1 = alpha/x*(exp(x*t)-1)*exp(alpha*(beta+x)*t/sigma_sq)*((1-y)/(1-y*exp(x*t)))^(2*alpha/sigma_sq+1)
  F2 = exp(alpha*(beta+x)*t/sigma_sq+x*t)*((1-y)/(1-y*exp(x*t)))^(2*alpha/sigma_sq+2)
  
  A = A1*exp(A2*lambda)
  F = (F1+F2*lambda)*exp(A2*lambda)
  
  r = seq(0, 0.03*T-0.03*T/N, by=0.03*T/N)
  si = gamma*sum(exp(-r)*F)/sum(exp(-r)*A)
  
  diff = (si-data_point)^2
}



#Below package help us to calculate polynomials to finally estimate the derivatives of the function wrt parameters
library(Deriv)
grad_alpha = Deriv(CDS_price_error, "alpha")
grad_beta = Deriv(CDS_price_error, "beta")
grad_sigma_sq = Deriv(CDS_price_error, "sigma_sq")



grad_alpha_s = Deriv(CDS_price, "alpha")
grad_beta_s = Deriv(CDS_price, "beta")
grad_sigma_sq_s = Deriv(CDS_price, "sigma_sq")

