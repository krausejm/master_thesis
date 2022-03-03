class pdfObject{
    public:
    Double_t m_p_gamma;
    pdfObject(Double_t p_gamma){
            m_p_gamma=p_gamma;
    }    
    Double_t operator()(Double_t *x, Double_t *p){
        Double_t sigma=p[0]; //beam asymmetry
        //these are the 10 acceptance coefficients
        Double_t a[5];
        Double_t b[5];
        for(int i=0;i<=4;i++){
            a[i]=p[i+1];
            b[i]=p[i+6];
        }
        Double_t nmrtr=1-1./2.*a[2]*m_p_gamma;
        Double_t enmrtr=1+m_p_gamma*sigma*TMath::Cos(2*(-45-x[0])*TMath::Pi()/180.);
        Double_t acc_sum=0;
        for(int i=0;i<=4;i++){
           acc_sum+=a[i]*TMath::Sin(i*x[0]*TMath::Pi()/180.)+b[i]*TMath::Cos(i*x[0]*TMath::Pi()/180.);
        }
        enmrtr*=acc_sum;
        return enmrtr/nmrtr;    
    }

 
};



void toy_MC(){
//std::cout<<"hello world"<<std::endl;
//get pdf for parallel configuration, thus take neg. value
pdfObject n_par(-0.25);
pdfObject n_bot(0.3);
//TF1* test = new TF1("test","(1-0.25*0.5*TMath::Cos(2*(-45-x)*TMath::Pi()/180))*1./10.5*(9.3+0.28*TMath::Cos(x*TMath::Pi()/180.)+0.24*TMath::Sin(3*x*TMath::Pi()/180.))",-180,180);
//test->SetLineColor(kBlue);
//test->SetLineStyle(kDashed);
TF1* par_sig = new TF1("par_sig",n_par,-180,180,11);
//sigma
par_sig->SetParameter(0,0.5);
//a0
par_sig->SetParameter(1,0);
//a1
par_sig->SetParameter(2,0);
//a2
par_sig->SetParameter(3,0);
//a3
par_sig->SetParameter(4,0.24/10.5);
//par_sig->SetParameter(4,0);
//a4
par_sig->SetParameter(5,0);
//b0
par_sig->SetParameter(6,9.3/10.5);
//par_sig->SetParameter(6,1);
//b1
par_sig->SetParameter(7,0.28/10.5);
//par_sig->SetParameter(7,0);
//b2
par_sig->SetParameter(8,0);
//b3
par_sig->SetParameter(9,0);
//b4
par_sig->SetParameter(10,0);

TF1* par_bkg = new TF1("par_bkg",n_par,-180,180,11);
//sigma bkg
par_bkg->SetParameter(0,-0.5);
//a0
par_bkg->SetParameter(1,0);
//a1
par_bkg->SetParameter(2,0);
//a2
par_bkg->SetParameter(3,0);
//a3
par_bkg->SetParameter(4,0.24/10.5);
//a4
par_bkg->SetParameter(5,0);
//b0
par_bkg->SetParameter(6,9.3/10.5);
//b1
par_bkg->SetParameter(7,0.28/10.5);
//b2
par_bkg->SetParameter(8,0);
//b3
par_bkg->SetParameter(9,0);
//b4
par_bkg->SetParameter(10,0);

TF1* bot_sig = new TF1("bot_sig",n_bot,-180,180,11);
//sigma
bot_sig->SetParameter(0,0.5);
//a0
bot_sig->SetParameter(1,0);
//a1
bot_sig->SetParameter(2,0);
//a2
bot_sig->SetParameter(3,0);
//a3
bot_sig->SetParameter(4,0.24/10.5);
//bot_sig->SetParameter(4,0);
//a4
bot_sig->SetParameter(5,0);
//b0
bot_sig->SetParameter(6,9.3/10.5);
//bot_sig->SetParameter(6,1);
//b1
bot_sig->SetParameter(7,0.28/10.5);
//bot_sig->SetParameter(7,0);
//b2
bot_sig->SetParameter(8,0);
//b3
bot_sig->SetParameter(9,0);
//b4
bot_sig->SetParameter(10,0);

TF1* bot_bkg = new TF1("bot_bkg",n_bot,-180,180,11);
//sigma bkg
bot_bkg->SetParameter(0,-0.5);
//a0
bot_bkg->SetParameter(1,0);
//a1
bot_bkg->SetParameter(2,0);
//a2
bot_bkg->SetParameter(3,0);
//a3
bot_bkg->SetParameter(4,0.24/10.5);
//a4
bot_bkg->SetParameter(5,0);
//b0
bot_bkg->SetParameter(6,9.3/10.5);
//b1
bot_bkg->SetParameter(7,0.28/10.5);
//b2
bot_bkg->SetParameter(8,0);
//b3
bot_bkg->SetParameter(9,0);
//b4
bot_bkg->SetParameter(10,0);

TRandom3* random=new TRandom3(0);
//no. of events is poisson distributed
int events_bot=random->Poisson(1000.);
int events_par=random->Poisson(800.);
//prompt peak events
//fraction of signal events is 95 percent (chosen)
int n_sig_bot = std::round(0.95*events_bot);
int n_bkg_bot=events_bot-n_sig_bot;
int n_sig_par = std::round(0.95*events_par);
int n_bkg_par=events_par-n_sig_par;
//sideband events
double weights[7]={15./210.,8./210.,4./210.,10./210.,14./210.,6./210.,11./210.};
int n_sideband_bot[7];
int n_sideband_par[7];
for(int i=0;i<7;i++){
    n_sideband_par[i]=std::round(1./weights[i]*n_bkg_par/7.);
    n_sideband_bot[i]=std::round(1./weights[i]*n_bkg_bot/7.);
}
ofstream myfile;
myfile.open ("test.txt");
myfile << "pol.\tphi\tweight\n";

TH1D* hnbot=new TH1D("hnbot",";phi;n_bot",24,-180,180);
TH1D* hnpar=new TH1D("hnpar",";phi;n_bot",24,-180,180);
TH1D* hn=new TH1D("hn",";phi;n",24,-180,180);
std::cout<<events_bot<<" "<<events_par<<std::endl;

//generate signal events for bot config.
for(int i=0;i<n_sig_bot;i++){
    double phi_0=bot_sig->GetRandom();
    myfile<<n_bot.m_p_gamma<<"\t"<<phi_0<<"\t"<<1<<"\n";
    hnbot->Fill(phi_0);
    hn->Fill(phi_0);
}
//generate bkg events for bot config
for(int i=0;i<n_bkg_bot;i++){
    myfile<<n_bot.m_p_gamma<<"\t"<<bot_bkg->GetRandom()<<"\t"<<1<<"\n";
}
//generate seven times sideband events for bot config
for(int i=0;i<7;i++){
    for(int j=0;j<n_sideband_bot[i];j++){
        myfile<<n_bot.m_p_gamma<<"\t"<<bot_bkg->GetRandom()<<"\t"<<-1*weights[i]<<"\n";    }
}
//generate signal events for par config.
for(int i=0;i<n_sig_par;i++){
    double phi_1=par_sig->GetRandom();
    myfile<<n_par.m_p_gamma<<"\t"<<phi_1<<"\t"<<1<<"\n";
    hnpar->Fill(phi_1);
    hn->Fill(phi_1);
}
//generate bkg events for par config
for(int i=0;i<n_bkg_par;i++){
    myfile<<n_par.m_p_gamma<<"\t"<<par_bkg->GetRandom()<<"\t"<<1<<"\n";
}
//generate seven times sideband events for par config
for(int i=0;i<7;i++){
    for(int j=0;j<n_sideband_par[i];j++){
    myfile<<n_par.m_p_gamma<<"\t"<<par_bkg->GetRandom()<<"\t"<<-1*weights[i]<<"\n";    
    }
}


}