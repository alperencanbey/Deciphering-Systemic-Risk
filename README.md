# Deciphering Systemic Risk

## Overview
A simulation exercise for an optimization problem to extract parameters from CDS data. Employed a stochastic intensity process ([Cox–Ingersoll–Ross model](https://en.wikipedia.org/wiki/Cox%E2%80%93Ingersoll%E2%80%93Ross_model)) to model daily default probabilities, integrating this into the CDS pricing formula which can be found in [Appendix of Ang and Longstaff 2013](https://www.sciencedirect.com/science/article/pii/S0304393213000585). The optimization involves a two-step process combining non-linear least squares for intensity values and gradient descent for parameter optimization, aiming to accurately recover parameters in the stochastic process which are mean, speed of mean reversion, and volatility. The study also aims to demonstrate the asymptotic properties and robustness of the estimator across different sample sizes.

## Tools Used
- **R:** Utilized for simulation and parameter estimation.
- **MATLAB:** Utilized for data visualization.

## Methodology
- **Simulation Set-Up:** Parameters are estimated through minimization of the objective function, focusing on systemic risk components in CDS pricing.
- **Data Generation:** Artificial CDS price data is generated to assess the reliability of the model, using predefined parameters and the standard square root model for daily intensity values.
- **Parameter Recovery:** The optimization process involves bootstrapping daily intensity values and employing a gradient descent algorithm to optimize parameters.

## Results and Impact
- Demonstrates the efficacy of the optimization algorithm in parameter recovery.
- Highlights the influence of sample size on the estimator's performance and its asymptotic behavior.

## Challenges and Learnings
- **Challenges:** Overcoming the complexity of modeling and simulating the CDS pricing mechanism.
- **Learnings:** Insights into the robustness of parameter estimation methods and their implications for financial modeling.

## Contact Information
For further discussion or inquiries, please reach out at alperen.canbey@upf.edu.

