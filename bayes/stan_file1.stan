data {
	int<lower=0> n;
	vector[n] x;
	vector[n] y;
	vector[n] dy;
}
parameters {
	real <lower=0> a;
	real b;
	real<lower=-pi(), upper=pi()> c;
}

model {
	y ~ normal(b + a * cos(x+c), dy);
	a ~ normal(0,100);
	b ~ normal(0,100);
	c ~ normal(0,3.14);
}
generated quantities {
	array[n] real y_tilde = normal_rng(b + a * cos(x+c),dy);
	vector[n] log_lik;
	for (k in 1:n){
    	log_lik[k]=normal_lpdf(y[k]|b + a * cos(x[k]+c), dy[k]);
    }

}
