---
layout: post
title:  "Efficient multifidelity likelihood-free Bayesian inference with adaptive computational resource allocation"
#date:   2021-06-01 18:08:27 +0100
categories: publications
#talks, research, background, publications
---

The latest collaboration with [Ruth Baker](https://www.iamruthbaker.com/) and [David Warne](https://twitter.com/davidjwarne) is now available to view on [arXiv](https://arxiv.org/abs/2112.11971), with Julia code available on [Github](https://github.com/tpprescott/mf-lf).

**Abstract:**
Likelihood-free Bayesian inference algorithms are popular methods for calibrating the parameters of complex, stochastic models, required when the likelihood of the observed data is intractable.
These algorithms characteristically rely heavily on repeated model simulations.
However, whenever the computational cost of simulation is even moderately expensive, the significant burden incurred by likelihood-free algorithms leaves them unviable in many practical applications.
The multifidelity approach has been introduced (originally in the context of approximate Bayesian computation) to reduce the simulation burden of likelihood-free inference without loss of accuracy, by using the information provided by simulating computationally cheap, approximate models in place of the model of interest.
The first contribution of this work is to demonstrate that multifidelity techniques can be applied in the general likelihood-free Bayesian inference setting.
Analytical results on the optimal allocation of computational resources to simulations at different levels of fidelity are derived, and subsequently implemented practically.
We provide an adaptive multifidelity likelihood-free inference algorithm that learns the relationships between models at different fidelities and adapts resource allocation accordingly, and demonstrate that this algorithm produces posterior estimates with near-optimal efficiency. 
