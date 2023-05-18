# Barrier-Options-Pricing
[Youssef ANNAKI](https://www.linkedin.com/in/youssef-annaki-a91ab5192/)
## Summary

The aim of this work is to build a Monte Carlo pricer (C++) of Up-&-In Calls (UIC) and assess the performance and the accuracy of the pricing for different types of simulations of the diffusion process of the underlying asset (stocks in our case).
For simplicity and error evaluation purposes, we focus on a GBM behavior of the stock price with constant volatility and interest rate paramters (this does not represent a limitation of the pricer as it can easly be generalized to more sophisticated stock price diffusion models: LV, SLV, SR-SLV,...). The focus on the Black-Scholes framework is mainly driven by the fact that we do have analytical formulas for barrier options prices which is very important in order to assess the accuracy of the developed pricing models.

This work focuses on the following subjects:
- Pricing an UIC using the price of an UOC and the Put-Call parity applied to barrier options.
- Pricing of an UIC using the discrete Euler Scheme for the stock diffusion process.
- Pricing of an UIC using the continuous Euler Scheme combined with Brownian diffusion bridges that will allow use to introduce a multiplicative correction term in the pricing formula:

```math
\mathbb{E}((\tilde{S}_T - K)_{+} \mathbb{1}_{sup_{t \in [0:T]}(\tilde{S}_t \leq D)}) = \mathbb{E}((\bar{S}_T - K)_{+} \mathbb{1}_{max_{0 \leq k \leq D} (\bar{S}_{t_k} \leq D)} \prod_{k = 0}^{N-1}(1 - exp(\frac{-2N}{T} \frac{(\bar{S}_{t_k} - D)(\bar{S}_{t_{k+1}} - D)}{\sigma^2(t_k, \bar{S}_{t_k})}))
  ```
  - Extension of the previous pricing methodology to moving barrier option ($` barrier_t = \beta * t + \alpha `$).
  - Application of the multi-level technique (Richardson-Romberg extrapolation) to the pricing of an UIC.
  - Comparaison of the pricing performances with the theoretical expected results in various academic papers (see Documents). We mainly focus on the error convergence rate with respect to the Monte Carlo simulation size and the temporal discretization parameter.

## Some Results

### Price Estimators (Static UIC) convergence & Error evolution with respect to the Monte Carlo sample size
![image](https://github.com/YsfAnnaki/Barrier-Options-Pricing/assets/134018406/a50a683c-60e9-40bb-985a-b26a0da677d6)
![image](https://github.com/YsfAnnaki/Barrier-Options-Pricing/assets/134018406/4877c898-2aa1-49a0-8ce1-3b5f8edd5f80)

### Estimators Error as a function of the time discretization parameter
![image](https://github.com/YsfAnnaki/Barrier-Options-Pricing/assets/134018406/0a6600cc-c8e4-4fc9-93e7-5ae6c1b9f7e1)

### Richardson-Romberg extrapolation (3 levels extrapolation)
![image](https://github.com/YsfAnnaki/Barrier-Options-Pricing/assets/134018406/0630ef3d-8bba-4319-b89a-1e983f5e33aa)

##

For any information, feedback or questions, please [Contact me](mailto:annaki.youssef@gmail.com?subject=[GitHub]%20Source%20Han%20Sans)

