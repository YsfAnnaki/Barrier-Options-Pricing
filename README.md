# Barrier-Options-Pricing
[Youssef Annaki](https://www.linkedin.com/in/youssef-annaki-a91ab5192/)
## Summary

The aim of this work is to build a Monte Carlo pricer (C++) of Up-&-In Calls (UIC) and assess the performance and the accuracy of the pricing for different types of simulations of the diffusion process of the underlying asset (stocks in our case).
For simplicity and error evaluation purposes, we focus on a GBM behavior of the stock price with constant volatility and interest rate paramters (this does not represent a limitation of the pricer as it can easly be generalized to more sophisticated stock price diffusion models: LV, SLV, SR-SLV,...). The focus on the Black-Scholes framework is mainly driven by the fact that we do have analytical formulas for barrier options prices which is very important in order to assess the accuracy of the developed pricing models.


Monte Carlo pricing (C++) of static/moving barrier options using Brownian diffusion Bridges and pre-conditioning/Importance sampling variance reduction techniques with an extension to the multi-step Richardson-Romberg extrapolation.
