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
  	a ~ normal(0,100);
  	b ~ normal(0,100);
}
generated quantities {
   array[n] real y_tilde = normal_rng(b + a * cos(x),dy);
}