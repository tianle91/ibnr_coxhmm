# IBNR with Cox Hidden Markov Model

Model counts of Incurred But Not Reported Claims (IBNR) with a stochastic rate model (Cox Process) whose arrival intensities are conditionally independent given a latent state which follows a Markov Chain (Hidden Markov Model).


# Model

Under certain assumptions on continuous rate intensity parameters, the discretized arrival process follows a Hidden Markov Model with Negative-Binomial marginals under the shape-scale parametrization as follows where `p = 1/(1+s)`.
```
P(X=n|shape=m, scale=s) = binom(n+m-1, m-1) p^m (1-p)^n
```

The latent state transitions are determined by the transition probability matrix `gamma`, so `P(C_{t+1}=j|C_t=i) = gamma[i,j]`. 
Counts of arrivals are conditionally independent given the latent states, so `P(X_t=k|C_t=i) = P(X=k|shape=m_i, scale=s)` which are Negative-Binomial marginals as defined earlier.

## Thinning

If each arrival is observed with probability `q`, the resulting process is also the same Hidden Markov Model with Negative-Binomial marginals, albeit with a scale multiplied by the thinning probability.
Instead of a scale parameter of `s`, as in earlier, we now have a scale parameter of `s*q` and all other parameters as well as the latent state remain unchanged.
For example, if `Y_t` is the thinned version of `X_t` defined earlier with thinning parameters `q_t`, `P(Y_t=k|C_t=i) = P(X=k|shape=m_i, scale=s*q_t)`.

## Reported/Unreported Claims

We only observe reported claims.
If the claims share a common reporting delay distribution (represented by random variable `U`) and occurs at time `t`, the probability of it being reported by `t_val > t` is `P(U < t_val-t)`. 
This is more or less the thinning parameter (i.e. `q_t`) for the Reported Claims process.
The Cox-HMM is fitted on the Reported Claims process in order to estimate the parameters for the True Claims process.
Then samples from the Unreported Claims process are generated to compute an estimate for the number of unreported claims.

## Demo



Check out the demo [here](https://github.com/tianle91/ibnr-coxhmm/blob/master/demo/fitting.ipynb) to fit a synthetically generated Reported Claims process where we can recover the parameters.

# Discussion

1. Stepwise selection of shape parameters is slow
2. Common reporting delay distribution

# Reference

```
Badescu A.L., Lin S., Tang D., A marked Cox model for the number of IBNR claims: Theory, Insurance: Mathematics and Economics, 2016, 69, 29 – 37.
```
