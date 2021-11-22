import nest_asyncio
nest_asyncio.apply()
import stan
import numpy as np
from scipy.optimize import curve_fit
import pandas as pd
import matplotlib.pyplot as plt
data = pd.read_csv("linear_model.txt",delimiter='\t',header=None)
data.columns=['x','y','dy']

def f(x,a,b):
    return a*np.cos(x)+b
popt, pcov = curve_fit(f,data['x'].values,data['y'].values,p0=[2,2],sigma=data['dy'].values)

stan_code= """
data {
  int<lower=0> n;
  vector[n] x;
  vector[n] y;
}
parameters {
  real a;
  real b;
  real<lower=0> sigma;
}
model {
  y ~ normal(b + a * x, sigma);
}
"""
schools_code = """
data {
  int<lower=0> J;         // number of schools
  real y[J];              // estimated treatment effects
  real<lower=0> sigma[J]; // standard error of effect estimates
}
parameters {
  real mu;                // population treatment effect
  real<lower=0> tau;      // standard deviation in treatment effects
  vector[J] eta;          // unscaled deviation from mu by school
}
transformed parameters {
  vector[J] theta = mu + tau * eta;        // school treatment effects
}
model {
  target += normal_lpdf(eta | 0, 1);       // prior log-density
  target += normal_lpdf(y | theta, sigma); // log-likelihood
}
"""
schools_data = {"J": 8,
                "y": [28,  8, -3,  7, -1,  1, 18, 12],
                "sigma": [15, 10, 16, 11,  9, 11, 10, 18]}

posterior = stan.build(schools_code, data=schools_data, random_seed=1)
fit = posterior.sample(num_chains=4, num_samples=1000)

#posterior = stan.build(stan_code, data=stan_data,random_seed=1)

#plt.errorbar(x=data['x'],y=data['y'],yerr=data['dy'],fmt='.')
#yfit=f(data['x'],*popt)
#plt.plot(data['x'],yfit)
#plt.show()