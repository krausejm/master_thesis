functions{

real mylpdf(real phi, real pol, real sigma,vector a, vector b){
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
                        //+m_a[2]*sin(angle)+m_b[2]*cos(angle);
                        //+m_a[3]*sin(2*phi*pi()/180.)+m_b[3]*cos(2*phi*pi()/180.));
                        //+m_a[4]*sin(2*phi*pi()/180.)+m_b[4]*cos(2*phi*pi()/180.));
                        //+m_a[5]*sin(4*phi*pi()/180.)+m_b[5]*cos(4*phi*pi()/180.));
  
    return log(nnmrtr/enmrtr);
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
}
model{

//loop over prmpt peak events
//real mu =my_lpdf(phi_prmpt[1]|pol_prmpt[1],sigma,a,b);
for(k in 1:N){
    target+=(f*mylpdf(phi_prmpt[k],pol_prmpt[k],sigma,a,b)+(1-f)*mylpdf(phi_prmpt[k],pol_prmpt[k],sigma_bkg,a_bkg,b_bkg));
}
//loop over sideband events
for(k in 1:M){
    target+=mylpdf(phi_side[k],pol_side[k],sigma_bkg,a_bkg,b_bkg);
}

//priors for eff. coefficients, non-informative, broadly around 0
for(k in 1:4){
    a[k] ~ normal(0,1);
    b[k] ~ normal(0,1);
    a_bkg[k] ~ normal(0,1);
    b_bkg[k] ~ normal(0,1);
}
sigma ~ normal(0,1) T[-1,1];
sigma_bkg ~ normal(0,1) T[-1,1];
}