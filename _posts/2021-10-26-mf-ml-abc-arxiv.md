---
layout: post
title:  "Multifidelity multilevel Monte Carlo to accelerate approximate Bayesian parameter inference for partially observed stochastic processes"
#date:   2021-06-01 18:08:27 +0100
categories: publications
#talks, research, background, publications
---

Over the past few months, [Ruth Baker](https://www.iamruthbaker.com/) and I have been engaged in a productive and interesting collaboration with [David Warne](https://twitter.com/davidjwarne) and [Mat Simpson](http://www.mj-simpson.com) from Queensland University of Technology, including a virtual 'visit' to Oxford from David back in March 2021, supported by the Australian Mathematics Society through a [Liftoff Fellowship](https://austms.org.au/awards-grants/awards/lift-off-fellowships/).
The results of our work are now available to view on [arXiv](https://arxiv.org/abs/2110.14082), with MATLAB code available on [Github](https://github.com/davidwarne/MLMCandMultifidelityForABC).

**Abstract:**
Models of stochastic processes are widely used in almost all fields of science.
Theory validation, parameter estimation, and prediction all require model calibration and statistical inference using data.
However, data are almost always incomplete observations of reality.
This leads to a great challenge for statistical inference because the likelihood function will be intractable for almost all partially observed stochastic processes.
This renders many statistical methods, especially within a Bayesian framework, impossible to implement.
Therefore, computationally expensive likelihood-free approaches are applied that replace likelihood evaluations with realisations of the model and observation process.
For accurate inference, however, likelihood-free techniques may require millions of expensive stochastic simulations.
To address this challenge, we develop a new method based on recent advances in multilevel and multifidelity.
Our approach combines the multilevel Monte Carlo telescoping summation, applied to a sequence of approximate Bayesian posterior targets, with a multifidelity rejection sampler to minimise the number of computationally expensive exact simulations required for accurate inference.
We present the derivation of our new algorithm for likelihood-free Bayesian inference, discuss practical implementation details, and demonstrate substantial performance improvements.
Using examples from systems biology, we demonstrate improvements of more than two orders of magnitude over standard rejection sampling techniques.
Our approach is generally applicable to accelerate other sampling schemes, such as sequential Monte Carlo, to enable feasible Bayesian analysis for realistic practical applications in physics, chemistry, biology, epidemiology, ecology and economics.
