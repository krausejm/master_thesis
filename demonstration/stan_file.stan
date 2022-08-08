data {
	int<lower=0> n; //number of datapoints
	vector[n] x; //predictors
	vector[n] y; //measured quantity
	vector[n] y_Err; //measurement error
}
parameters {
	real a;
	real b;
}
model {
	//likelihood
	y ~ normal(a*x+b, y_err);
	//priors
	a ~ normal(0,1);
	b ~ normal(0,1);
}

