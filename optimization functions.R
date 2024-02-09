###To recover the parameters, the optimization process has two legs
####combining non-linear least squares for intensity values and gradient descent for parameter optimization.

#This is for the first step. Having the parameters at hand, this function uses NLS to create intensity values.
bootstrap_intensity <- function(parameters, lambda_initial, cds_data){
  
  intensity = rep(NA,nrow(cds_data))
  days_iy = nrow(cds_data)

  for (i in 1:nrow(cds_data)) {
    yy = cds_data[i,]
    m<-nls(yy~ c(CDS_price(parameters, lambda = lambda , days_iy , maturity = 1),
                 CDS_price(parameters, lambda = lambda , days_iy , maturity = 2),
                 CDS_price(parameters, lambda = lambda , days_iy , maturity = 3),
                 CDS_price(parameters, lambda = lambda , days_iy , maturity = 4),
                 CDS_price(parameters, lambda = lambda , days_iy , maturity = 5)), start = list(lambda=lambda_initial)) 
    
    
    intensity[i] = coefficients(m)
  }
  return(intensity)
}


#Gradient Descent algorithm to optimize parameters st minimizing the squared errors
gradient_descent <- function(parameters, intensity, cds_data, step_size){
  
  days_iy = nrow(cds_data)
  
  contribution_alpha = rep(NA,nrow(cds_data)*ncol(cds_data))
  contribution_beta = rep(NA,nrow(cds_data)*ncol(cds_data))
  contribution_sigma_sq = rep(NA,nrow(cds_data)*ncol(cds_data))
  for (i in 1:nrow(cds_data)) {
    for (j in 1:ncol(cds_data)) {
      contribution_alpha[j+ncol(cds_data)*(i-1)] = grad_alpha(parameters, intensity[i], days_iy, j, cds_data[i,j])[[1]]
      contribution_beta[j+ncol(cds_data)*(i-1)] = grad_beta(parameters, intensity[i], days_iy, j, cds_data[i,j])[[1]]
      contribution_sigma_sq[j+ncol(cds_data)*(i-1)] = grad_sigma_sq(parameters, intensity[i], days_iy, j, cds_data[i,j])[[1]]
    }
  }
  contributions = c(
    sum(contribution_alpha),
    sum(contribution_beta),
    sum(contribution_sigma_sq))
  
  parameters[which.max(abs(contributions))] = parameters[which.max(abs(contributions))] - step_size*contributions[which.max(abs(contributions))]
  
  return(parameters)
}



#To calculate the Jacobian Matrix 
jacobian <- function(parameters, intensity, cds_data){
  
  days_iy = nrow(cds_data)
  
  contribution_alpha = matrix(NA, nrow = nrow(cds_data), ncol = ncol(cds_data), byrow = TRUE)
  contribution_beta =  matrix(NA, nrow = nrow(cds_data), ncol = ncol(cds_data), byrow = TRUE)
  contribution_sigma_sq =  matrix(NA, nrow = nrow(cds_data), ncol = ncol(cds_data), byrow = TRUE)
  for (i in 1:nrow(cds_data)) {
    for (j in 1:ncol(cds_data)) {
      contribution_alpha[i,j] = grad_alpha(parameters, intensity[i], days_iy, j, cds_data[i,j])[[1]]
      contribution_beta[i,j] = grad_beta(parameters, intensity[i], days_iy, j, cds_data[i,j])[[1]]
      contribution_sigma_sq[i,j] = grad_sigma_sq(parameters, intensity[i], days_iy, j, cds_data[i,j])[[1]]
    }
  }
  jacobian_matrix = rbind(
    colSums(contribution_alpha),
    colSums(contribution_beta),
    colSums(contribution_sigma_sq))
  
  return(jacobian_matrix)
}


#to calculate the gradient
gradient <- function(parameters, intensity, cds_data){
  
  days_iy = nrow(cds_data)
  
  contribution_alpha = matrix(NA, nrow = nrow(cds_data), ncol = ncol(cds_data), byrow = TRUE)
  contribution_beta =  matrix(NA, nrow = nrow(cds_data), ncol = ncol(cds_data), byrow = TRUE)
  contribution_sigma_sq =  matrix(NA, nrow = nrow(cds_data), ncol = ncol(cds_data), byrow = TRUE)
  for (i in 1:nrow(cds_data)) {
    for (j in 1:ncol(cds_data)) {
      contribution_alpha[i,j] = grad_alpha_s(parameters, intensity[i], days_iy, j)[[1]]
      contribution_beta[i,j] = grad_beta_s(parameters, intensity[i], days_iy, j)[[1]]
      contribution_sigma_sq[i,j] = grad_sigma_sq_s(parameters, intensity[i], days_iy, j)[[1]]
    }
  }
  jacobian_matrix = rbind(
    colSums(contribution_alpha),
    colSums(contribution_beta),
    colSums(contribution_sigma_sq))
  
  return(jacobian_matrix)
}



#This is just to calculate the residuals
error <- function(parameters, intensity, cds_data){
  meansq = rep(NA,nrow(cds_data)*ncol(cds_data))
  days_iy = nrow(cds_data)
  for (i in 1:nrow(cds_data)) {
    for (j in 1:ncol(cds_data)) {
      meansq[j+ncol(cds_data)*(i-1)] = (CDS_price(parameters, lambda = intensity[i], days_iy , maturity = j) - cds_data[i,j])^2
    }
  }
  sum(meansq)/(nrow(cds_data)*ncol(cds_data))
}



#Variance measurement error
var_measurement_error <- function(parameters, intensity, cds_true_values){
  measurement_err = rep(NA,nrow(cds_true_values)*ncol(cds_true_values))
  days_iy = nrow(cds_true_values)
  
  for (i in 1:nrow(cds_true_values)) {
    for (j in 1:ncol(cds_true_values)) {
      measurement_err[j+ncol(cds_true_values)*(i-1)] = cds_true_values[i,j] - CDS_price(parameters, lambda = intensity[i], days_iy , maturity = j)
    }
  }
  var(measurement_err)
}



parameter_errors <- function(final_parameters, true_parameters){
  #avg_deviation = colSums(final_parameters - true_parameters)/100
  #avg_absolute_deviation = colSums(abs(final_parameters - true_parameters))/100
  mean_sq_deviation = colMeans((final_parameters - true_parameters)^2)
  estimated_mean = colMeans(final_parameters)
  variance_errors = var(final_parameters - true_parameters)
  
  output <- list("Estimated Mean" = estimated_mean, "Variance of Errors" = variance_errors, "Mean Squared Error" = mean_sq_deviation)
  return(output) 
}


