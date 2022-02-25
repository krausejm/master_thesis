data {
	int<lower=0> n;
	vector[n] x;
	vector[n] y;
	vector[n] dy;
}
parameters {
	real<lower=-1, upper=1> a;
	//real b;
}


model {
	for(k in 1:n){
		//likelihood
		y[k] ~ normal(a * cos(pi()/180.*2*(-45-x[k])), dy[k]);
	}
	//prior
	a ~ normal(0,1) T[-1,1];
	//b ~ normal(0,1);
}
generated quantities {
	real y_tilde[n]; 
	vector[n] log_lik;
	vector[n] mu;
	mu=a*cos(pi()/180.*2*(-45-x));
	for (k in 1:n){
    	log_lik[k]=normal_lpdf(y[k]|mu[k], dy[k]);
		y_tilde[k]=normal_rng(mu[k],dy[k]);
	}
}
