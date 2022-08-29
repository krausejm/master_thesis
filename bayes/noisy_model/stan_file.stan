data {
	int<lower=0> n;
	vector[n] x;
	vector[n] y;
	vector[n] dy;
}
parameters {
	real a;
	real b;
}

model {
	y ~ normal(b + a * cos(x), dy);
	a ~ normal(0,1);
	b ~ normal(0,1);
}
generated quantities {
	real y_tilde[n] = normal_rng(b + a * cos(x),dy);
}




