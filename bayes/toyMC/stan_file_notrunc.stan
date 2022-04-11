data {
	int<lower=0> n;
	vector[n] phi;
	vector[n] asym;
	vector[n] dasym;
}
parameters {
	real sigma;
	//real b;
}


model {
	for(k in 1:n){
		//likelihood
		asym[k] ~ normal(sigma * cos(pi()/180.*2*(-45-phi[k])), dasym[k]);
	}
	//prior
	sigma ~ normal(0,1);
	//b ~ normal(0,1);
}
generated quantities {
	real y_tilde[n]; 
	vector[n] log_lik;
	vector[n] mu;
	mu=sigma*cos(pi()/180.*2*(-45-phi));
	for (k in 1:n){
    	log_lik[k]=normal_lpdf(asym[k]|mu[k], dasym[k]);
		y_tilde[k]=normal_rng(mu[k],dasym[k]);
	}
}
