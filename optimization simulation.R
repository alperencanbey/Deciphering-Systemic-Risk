#run the functions in these R files before running this optimization simulation:
##data_pricing_functions.R, optimization_functions.R

##days in year * step size = FIXED  = 40 x 0.1

optimization_simulation <- function(true_parameters, initial_parameters, lambda0, days_in_year, year, random_error, step_size, Nsimulation) {
  
  parameters = initial_parameters
  
  simul <- function(parameters) {
    
    #DATA CREATION
    data_created = data_creation(true_parameters = true_parameters, lambda0 = lambda0, days_in_year = days_in_year, year, sample_count_BM = 500, random_error)
    cds_data = data_created[[1]]
    daily_intensity = data_created[[2]]
    cds_true_values = data_created[[3]]
    tryCatch({
      #OPTIMIZATION
      #in case of error termination
      parameters_ex = parameters
      intensity = daily_intensity
      ###################
    error_check = rep(NA,100)
    intensity = bootstrap_intensity(parameters = parameters, lambda_initial = 0.5, cds_data = cds_data)
    error_check[1] =  error(parameters, intensity, cds_data)
 
    parameters = gradient_descent(parameters, intensity, cds_data, step_size)
    error_check[2] = error(parameters, intensity, cds_data)
    intensity = bootstrap_intensity(parameters, mean(intensity), cds_data)
    error_check[3] = error(parameters, intensity, cds_data)
    
    threshold = 0
    i = 2
    #this threshold may vary with random error, sample size and step size 1e-16
    #conditions are to stop if the error improvement is lower than the threshold (error*1e-10)
    while (error_check[i-1]-error_check[i] > threshold & error_check[i] - error_check[i+1] > threshold) {
      
      i = i + 2
      
      parameters_ex = parameters
      parameters = gradient_descent(parameters, intensity, cds_data, step_size)
      error_check[i] = error(parameters, intensity, cds_data)
      
      intensity = bootstrap_intensity(parameters, mean(intensity), cds_data)
      error_check[i+1] =  error(parameters, intensity, cds_data)
   
      threshold = error(parameters, intensity, cds_data) * 1e-10
    }
      
    
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    parameters = parameters_ex
    
#VALUES TO REPORT
mispricing_error = error(parameters, intensity, cds_data)
intensity_error = sum((daily_intensity-intensity)^2)
mispricing_error_2 = error(parameters, intensity, cds_true_values)
variance_measurement_error = var_measurement_error(parameters, intensity, cds_true_values)

output <- list("alpha" = parameters[1], "beta" = parameters[2], "sigma_sq" = parameters[3], "mispricing_error_1" = mispricing_error, "mispricing_error_2" = mispricing_error_2, "variance_measurement_error" = variance_measurement_error, "intensity error" = intensity_error, "iterationN" = i/2)
return(output) 
  }
  
  #CONDUCTING SIMULATIONS
  output_info = data.frame(Doubles=double(),
                           Doubles=double(),
                           Doubles=double(),
                           Doubles=double(),
                           Doubles=double(),
                           Doubles=double(),
                           Doubles=double(),
                           Doubles=double())
  
  colnames(output_info) = c("alpha", "beta", "sigma_sq", "mispricing error_1", "mispricing error_2", "variance_measurement_error", "intensity error", "iterationN")
  
  for (i in 1:Nsimulation) {
    results = simul(parameters)
    output_info[i,] =  data.frame(results)
  }
  
  return(output_info)
  
}


lambda0 = 0.2 
year = 1
days_in_year= 36
random_error = 0.01
step_size = 0.01


#true parameters = c(alpha, beta, sigma_sq)
true_parameters = c(0.4, 0.3, 0.2)

initial_parameters = c(0.4, 0.3, 0.2)

output_info = optimization_simulation(true_parameters, initial_parameters , lambda0 = 0.2, days_in_year = 36, year = 1, random_error = 0.01, step_size = 0.01, Nsimulation=100)

