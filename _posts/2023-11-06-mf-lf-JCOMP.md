---
layout: post
title:  "Published! Efficient multifidelity likelihood-free Bayesian inference with adaptive computational resource allocation"
categories: publications
---

My collaborative work with [Ruth Baker](https://www.iamruthbaker.com/) and [David Warne](https://twitter.com/davidjwarne)
is now accepted for publication at the [Journal of Computational Physics](https://doi.org/10.1016/j.jcp.2023.112577).

This will have a publication date of 1st January 2024!

The accompanying Julia code is available on [Github](https://github.com/tpprescott/mf-lf).

**Abstract:**
Likelihood-free Bayesian inference algorithms are popular methods for inferring the parameters of complex stochastic models with intractable likelihoods.
These algorithms characteristically rely heavily on repeated model simulations.
However, whenever the computational cost of simulation is even moderately expensive, the significant burden incurred by likelihood-free algorithms leaves them infeasible for many practical applications. 
The multifidelity approach has been introduced in the context of approximate Bayesian computation to reduce the simulation burden of likelihood-free inference without loss of accuracy, by using the information provided by simulating computationally cheap, approximate models in place of the model of interest. 
In this work we demonstrate that multifidelity techniques can be applied in the general likelihood-free Bayesian inference setting. 
Analytical results on the optimal allocation of computational resources to simulations at different levels of fidelity are derived, and subsequently implemented practically. 
We provide an adaptive multifidelity likelihood-free inference algorithm that learns the relationships between models at different fidelities and adapts resource allocation accordingly, and demonstrate that this algorithm produces posterior estimates with near-optimal efficiency.
