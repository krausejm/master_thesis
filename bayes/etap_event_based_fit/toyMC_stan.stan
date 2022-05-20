functions{

real mypdf(real phi, real pol, real sigma, vector a, vector b){
    vector[5] m_a;
    vector[5] m_b;
    //the first coefficients are set!
    m_a[1]=0;
    m_b[1]=1;
    //for convenience make vectors of length 5 for easy use in loops
    for(k in 1:4){
        m_a[k+1]=a[k];
        m_b[k+1]=b[k];
    }
    //enumerator and numerator of pdf
    real enmrtr = 1-0.5*m_a[3]*pol*sigma;
    real nmrtr = 1+pol*sigma*cos(2*(-45-phi)*pi()/180.);
    //have to multiply numerator by efficiency function
    real angle=phi*pi()/180.;
    vector[5] sines;
    vector[5] cosines;
    for(k in 0:4){
        sines[k+1]=sin(k*angle);
        cosines[k+1]=cos(k*angle);
    }
    real nnmrtr=nmrtr*(dot_product(sines,m_a)+dot_product(cosines,m_b));
    return nnmrtr/enmrtr;
}
}
data {
    //prmpt peak data
	int<lower=0> N;
	vector[N] phi_prmpt;
	vector[N] pol_prmpt;
    //sideband data
    int<lower=0> M;
    vector[M] phi_side;
    vector[M] pol_side;
    //fraction of bkg events in prmpt peak
    real<lower=0,upper=1> f;
    //fraction of 2pi0 and etap
    real<lower=0,upper=1> f_s;
    real<lower=0,upper=1> f_b;
    //measurement of sigma 2pi0
    real sigma_2pi0_meas;
    real dsigma_2pi0_meas;
}
parameters {
    //to do: dont use vector??
    //18 pars in total

    //sigma and detector coefficients for signal events
	real<lower=-1, upper=1> sigma;
	vector[4] a;

    vector[4] b;

    //sigma and detector coefficients for bkg events
    real<lower=-1, upper=1> sigma_bkg;
	vector[4] a_bkg;
 
    vector[4] b_bkg;

    //true (unknown) value of sigma_2pi0 and prior location
    real sigma_2pi0;
    //real mu_sigma_2pi0;
    //real std_sigma_2pi0;

}
model{
//priors for eff. coefficients, non-informative, broadly around 0
for(k in 1:4){
    a[k]~normal(0,0.1);
    b[k]~normal(0,0.1);
    a_bkg[k]~normal(0,0.1);
    b_bkg[k]~normal(0,0.1);
}
//priors for sigma
sigma ~ normal(0,1) T[-1,1];
sigma_bkg ~ normal(0,1) T[-1,1];
sigma_2pi0 ~uniform(-1,1);
//sigma_2pi0 ~ normal(mu_sigma_2pi0,std_sigma_2pi0);
sigma_2pi0_meas ~ normal(sigma_2pi0,dsigma_2pi0_meas);
//loop over prmpt peak events
for(k in 1:N){
    target+=log(f*mypdf(phi_prmpt[k],pol_prmpt[k],f_s*sigma+f_b*sigma_2pi0_meas,a,b)+(1-f)*mypdf(phi_prmpt[k],pol_prmpt[k],sigma_bkg,a_bkg,b_bkg));
}
//loop over sideband events
for(k in 1:M){
    target+=log(mypdf(phi_side[k],pol_side[k],sigma_bkg,a_bkg,b_bkg));
}


}
