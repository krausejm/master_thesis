class pdfObjectbot{
    public:
    Double_t m_p_gamma;
    pdfObjectbot(Double_t p_gamma){
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
        Double_t nmrtr=1-0.5*a[2]*m_p_gamma*sigma;
        Double_t enmrtr=1+m_p_gamma*sigma*TMath::Cos(2*(-45-x[0])*TMath::Pi()/180.);
        Double_t acc_sum=0;
        for(int i=0;i<=4;i++){
           acc_sum+=a[i]*TMath::Sin(i*x[0]*TMath::Pi()/180.)+b[i]*TMath::Cos(i*x[0]*TMath::Pi()/180.);
        }
        enmrtr*=acc_sum;
        return enmrtr/nmrtr;    
    }

 
};
class pdfObjectpar{
    public:
    Double_t m_p_gamma;
    pdfObjectpar(Double_t p_gamma){
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
        Double_t nmrtr=1+0.5*a[2]*m_p_gamma*sigma;
        Double_t enmrtr=1-m_p_gamma*sigma*TMath::Cos(2*(-45-x[0])*TMath::Pi()/180.);
        Double_t acc_sum=0;
        for(int i=0;i<=4;i++){
           acc_sum+=a[i]*TMath::Sin(i*x[0]*TMath::Pi()/180.)+b[i]*TMath::Cos(i*x[0]*TMath::Pi()/180.);
        }
        enmrtr*=acc_sum;
        return enmrtr;    
    }

 
};



void toy_MC(){
//std::cout<<"hello world"<<std::endl;
//get pdf for parallel configuration, thus take neg. value
pdfObjectpar n_par(0.25);
pdfObjectbot n_bot(0.3);

TF1* par_sig = new TF1("par_sig",n_par,-180,180,11);
TF1* par_bkg = new TF1("par_bkg",n_par,-180,180,11);
TF1* bot_sig = new TF1("bot_sig",n_bot,-180,180,11);
TF1* bot_bkg = new TF1("bot_bkg",n_bot,-180,180,11);

//order of params: sigma,a0-a4,b0-b4
double sig_params[11]={0.5,0,0,0,0.24/10.5,0,9.3/10.5,0.28/10.5,0,0,0};
//double sig_params[11]={0.5,0,0,0,0,0,1,0,0,0,0};
double bkg_params[11]={-0.5,0,0,0,0.24/10.5,0,9.3/10.5,0.28/10.5,0,0,0};

for(int i=0;i<11;i++){
    par_sig->SetParameter(i,sig_params[i]);
    par_bkg->SetParameter(i,bkg_params[i]);
    bot_sig->SetParameter(i,sig_params[i]);
    bot_bkg->SetParameter(i,bkg_params[i]);
}
for(int i=0;i<11;i++){
    std::cout<<par_sig->GetParameter(i)<<std::endl;
}



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
int n_sideband_bot=20*n_bkg_bot;
int n_sideband_par=20*n_bkg_par;
//sideband events
/*double weights[7]={15./210.,8./210.,4./210.,10./210.,14./210.,6./210.,11./210.};
int n_sideband_bot[7];
int n_sideband_par[7];
for(int i=0;i<7;i++){
    n_sideband_par[i]=std::round(1./weights[i]*n_bkg_par/7.);
    n_sideband_bot[i]=std::round(1./weights[i]*n_bkg_bot/7.);
}*/
ofstream myfile;
myfile.open ("./toybins/toybin01.txt");
myfile << "pol.\tphi\tweight\n";

TH1D* hnbot=new TH1D("hnbot",";phi;n_bot",24,-180,180);
TH1D* hnpar=new TH1D("hnpar",";phi;n_bot",24,-180,180);
TH1D* hn=new TH1D("hn",";phi;n",24,-180,180);
TF1* eff = new TF1("f1","[0]*1./10.5*(9.3+0.28*cos(TMath::Pi()/180.*x)+0.24*sin(TMath::Pi()/180.*3*x))",-180,180);
std::cout<<events_bot<<" "<<events_par<<std::endl;
//signal events for bot config
for(int i=0;i<n_sig_bot;i++){
    myfile<<n_bot.m_p_gamma<<"\t"<<bot_sig->GetRandom()<<"\t"<<1<<"\n";
}
//bkg events for bot config
for(int i=0;i<n_bkg_bot;i++){
    myfile<<n_bot.m_p_gamma<<"\t"<<bot_bkg->GetRandom()<<"\t"<<1<<"\n";
}
//sideband events for bot config
for(int i=0;i<n_sideband_bot;i++){
    myfile<<n_bot.m_p_gamma<<"\t"<<bot_bkg->GetRandom()<<"\t"<<-1./20.<<"\n";
}


//signal events for par config
for(int i=0;i<n_sig_par;i++){
    myfile<<n_par.m_p_gamma<<"\t"<<par_sig->GetRandom()<<"\t"<<1<<"\n";
}
//bkg events for bot config
for(int i=0;i<n_bkg_par;i++){
    myfile<<n_par.m_p_gamma<<"\t"<<par_bkg->GetRandom()<<"\t"<<1<<"\n";
}
//sideband events for bot config
for(int i=0;i<n_sideband_par;i++){
    myfile<<n_par.m_p_gamma<<"\t"<<par_bkg->GetRandom()<<"\t"<<-1./20.<<"\n";
}





double scale_bot=1./hnbot->Integral();
double scale_par=1./hnpar->Integral();
//hnbot->Scale(scale_bot);
//hnpar->Scale(scale_par);


//hnbot->Add(hnpar);
//hnbot->Draw("pe");
//hnbot->Fit("f1");



}