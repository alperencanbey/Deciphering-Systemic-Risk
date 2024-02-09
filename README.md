# Deciphering-Systemic-Risk


The aim of the simulation exercises is to measure the performance of the designed optimization algorithm in recovering the parameters that explain the stochastic intensity process and the CDS pricing mechanisms. These exercises are conducted with different sample sizes to demonstrate the asymptotic properties and robustness of the estimator. The parameters are estimated through the minimization process of the objective function which is the sum of squared differences between the estimated and observed CDS prices across the different dates and maturities. 
\begin{equation*}
\sum_{T=1}\sum_{t=1} (\hat{s}_{t,T}-s_{t,T})^2
\end{equation*}
where T is the maturity of the CDS contract and t represents the days where CDS price data observed. 
