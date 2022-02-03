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
	y ~ normal(a * cos(pi()/180.*2*(-45-x)), dy);
	a ~ normal(0,1);
	//b ~ normal(0,1);
}
generated quantities {
	real y_tilde[n] = normal_rng(a * cos(pi()/180.*2*(-45-x)),dy);
	vector[n] log_lik;
	vector[n] mu;
	mu=a*cos(pi()/180.*2*(-45-x));
	for (k in 1:n){
    	log_lik[k]=normal_lpdf(y[k]|mu[k], dy[k]);
    }
}
