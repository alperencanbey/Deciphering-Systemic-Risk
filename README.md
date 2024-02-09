# Deciphering-Systemic-Risk

A simulation exercise for an optimization problem to extract parameters from CDS data. Employed a stochastic intensity process (Cox–Ingersoll–Ross model) to model daily default probabilities, integrating this into the CDS pricing formula which can be found in Appendix of Ang and Longstaff 2013. The optimization involves a two-step process combining non-linear least squares for intensity values and gradient descent for parameter optimization, aiming to accurately recover parameters in the stochastic process which are mean, speed of mean reversion, and volatility. 
