data {
	int<lower=0> n; //number of data points
	vector[n] x; //predictors
	vector[n] y; //measurements with
	vector[n] dy; //corresponding errors
}
parameters {
	real a; //amplitude and offset parameters
	real b;
}

model {
	y ~ normal(b + a * cos(x), dy); //likelihood
	a ~ normal(0,1); // prior
	b ~ normal(0,1); // prior
}
generated quantities {
	//posterior predictive check
	real y_tilde[n] = normal_rng(b + a * cos(x),dy); 
}




