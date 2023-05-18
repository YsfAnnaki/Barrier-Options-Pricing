# Barrier-Options-Pricing
[Youssef ANNAKI](https://www.linkedin.com/in/youssef-annaki-a91ab5192/)
## Summary

The aim of this work is to build a Monte Carlo pricer (C++) of Up-&-In Calls (UIC) and assess the performance and the accuracy of the pricing for different types of simulations of the diffusion process of the underlying asset (stocks in our case).
For simplicity and error evaluation purposes, we focus on a GBM behavior of the stock price with constant volatility and interest rate paramters (this does not represent a limitation of the pricer as it can easly be generalized to more sophisticated stock price diffusion models: LV, SLV, SR-SLV,...). The focus on the Black-Scholes framework is mainly driven by the fact that we do have analytical formulas for barrier options prices which is very important in order to assess the accuracy of the developed pricing models.

This work focuses on the following subjects:
- Pricing of an UIC using the discrete Euler Scheme for the stock diffusion process.
- Pricing of an UIC using the continuous Euler Scheme combined with Brownian diffusion bridges that will allow use to introduce a multiplicative correction term in the pricing formula:
  $$ \mathbb{E}((\tilde{S}_T - K)_{+} \mathbbm{1}_{\stackunder{sup}{t \in [0:T]}\tilde{S}_t \leq D}) = \mathbb{E}((\bar{S}_T - K)_{+} \mathbbm{1}_{\stackunder{max}{0 \leq k \leq D} \bar{S}_{t_k} \leq D} \prod_{k = 0}^{N-1}(1 - exp(\frac{-2N}{T} \frac{(\bar{S}_{t_k} - D)(\bar{S}_{t_{k+1}} - D)}{\sigma^2(t_k, \bar{S}_{t_k})})) $$


Monte Carlo pricing (C++) of static/moving barrier options using Brownian diffusion Bridges and pre-conditioning/Importance sampling variance reduction techniques with an extension to the multi-step Richardson-Romberg extrapolation.
