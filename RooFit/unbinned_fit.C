//simple 3d fit
double func3(Double_t *x, Double_t *par){
    Double_t r2 = x[0]*x[0] + x[1]*x[1] +x[2]*x[2];
    Double_t arg = (TMath::Sqrt(r2) - par[1])/par[2];
    Double_t val = par[0]/TMath::Power(2*TMath::Pi()*par[2]*par[2],1.5)*TMath::Exp(-0.5*arg*arg);
    return val;
}
//this is the real deal
double pdf(double* x, double* par){
    double phi=x[0];
    double pol=x[1];
    double weight=x[2];


    double sigma = par[0];
    double sigma_bkg = par[9];
    double f = par[18]; //fraction of bkg events in prompt peak
    //eff coefficients
    double a[5];
    double b[5];
    double a_bkg[5];
    double b_bkg[5];
    //first coefficients are set!
    a[0]=0;
    b[0]=1;
    a_bkg[0]=0;
    b_bkg[0]=1;
    //get meaningful var names
    for(int i=1;i<5;i++){
        a[i]=par[i];
        b[i]=par[i+4];
        a_bkg[i]=par[i+9];
        b_bkg[i]=par[i+13];
    }
    //first get efficiency
    double eff=0;
    double eff_bkg=0;
    for(int i=0;i<5;i++){
        eff+=a[i]*sin(i*phi*TMath::DegToRad())+b[i]*cos(i*phi*TMath::DegToRad());
        eff_bkg+=a_bkg[i]*sin(i*phi*TMath::DegToRad())+b_bkg[i]*cos(i*phi*TMath::DegToRad());
    }
    //different scenarios for prompt peak and sideband events
    if(weight<0){//sideband events
        return (1+pol*sigma_bkg*cos(2*(-45-phi)*TMath::DegToRad()))*eff_bkg/(1-0.5*a[2]*pol*sigma_bkg);
    }else if(weight==1){//prompt peak events
        return f*(1+pol*sigma*cos(2*(-45-phi)*TMath::DegToRad()))*eff/(1-0.5*a[2]*pol*sigma)
                +(1-f)*(1+pol*sigma_bkg*cos(2*(-45-phi)*TMath::DegToRad()))*eff_bkg/(1-0.5*a[2]*pol*sigma_bkg);

    }else{
        return 0;
    }

}
double eff(double *x, double* par){
    double scale=par[0];
    double a[5];
    double b[5];
    double phi=x[0];
    
    //first coefficients are set!
    a[0]=0;
    b[0]=1;

    //get meaningful var names
    for(int i=1;i<5;i++){
        a[i]=par[i];
        b[i]=par[i+4];
    }
    //first get efficiency
    double eff=0;
    for(int i=0;i<5;i++){
        eff+=a[i]*sin(i*phi*TMath::DegToRad())+b[i]*cos(i*phi*TMath::DegToRad());
    }
    return scale*eff;

}

void unbinned_fit(){
    TF3* f = new TF3("mypdf",pdf,-180,180,-1,1,-1.2,1.2,19);
    const char * parnames[19]={"sigma","a1","a2","a3","a4","b1","b2","b3","b4","sigma_bkg","a1_bkg","a2_bkg","a3_bkg","a4_bkg","b1_bkg","b2_bkg","b3_bkg","b4_bkg","f"};
    //name parameters to identify later
    for(int i=0;i<19;i++){
        f->SetParName(i,parnames[i]);
    }
    TTree* t = new TTree("t","mytree");
    t->ReadFile("../bayes/etap_event_based_fit/toybins/toybin0000.txt","pol:phi:weight");
    double fr=(t->GetEntries("weight==1")-0.065*t->GetEntries("weight<0"))/t->GetEntries("weight==1");
    std::cout<<fr<<std::endl;
    f->FixParameter(18,fr);
    t->UnbinnedFit("mypdf","phi:pol:weight");
    
    TH1F* hsigma = new TH1F("hsigma",";#Sigma;counts;",100,-1,1);
    TH1F* hsigma_bkg = new TH1F("hsigma_bkg",";#Sigma;counts;",100,-1,1);
    TH1F* hxi = new TH1F("res",";#xi;counts",100,-10,10);
    
    hsigma->Fill(f->GetParameter(0));
    hsigma_bkg->Fill(f->GetParameter(9));
    hxi->Fill((f->GetParameter(0)-0.34)/f->GetParError(0));
   
    TFile* f = new TFile("./toyMC_results.root","RECREATE");
    hsigma->Write();
    hsigma_bkg->Write();
    hxi->Write();
    f->Close(); 
    
    
    /*//check fit of efficiency function

    TH1F* np45 = new TH1F("np45",";#phi / deg;counts",12,-180,180);
    TH1F* nm45 = new TH1F("nm45",";#phi;counts",12,-180,180);
    TH1F* hp45pol = new TH1F("hp45pol",";pol.;counts",100,-1,1);
    TH1F* hm45pol = new TH1F("hm45pol",";pol.;counts",100,-1,1);

    t->Draw("phi>>np45","pol>0&&weight==1","goff");
    t->Draw("phi>>nm45","pol<0&&weight==1","goff");
    //get pol vals
    t->Draw("pol>>hp45pol","pol>0&&weight==1","goff");
    t->Draw("pol>>hm45pol","pol<0&&weight==1","goff");

    double norm_p = np45->Integral();
    double norm_m = nm45->Integral();
    //normalize event yields
    np45->Scale(1./norm_p);
    nm45->Scale(1./norm_m);
    //get mean polarization
    double pol_p45=hp45pol->GetMean();
    double pol_m45=TMath::Abs(hm45pol->GetMean());
    //std::cout<<pol_p45;
    //build weighted sum s.t. only eff function should be visible
    np45->Add(np45,nm45,pol_m45/(pol_p45+pol_m45),pol_p45/(pol_m45+pol_p45));
    np45->Scale(1./np45->GetMaximum());
    //retrieve fit parameters
    TF1* eff_func=new TF1("eff_func",eff,-180,180,9);
    for(int i=1;i<9;i++){
        //write eff coefficients to new func
        eff_func->FixParameter(i,f->GetParameter(i));
    }
    //cosmetics
    eff_func->SetLineColor(kBlue);
    eff_func->SetNpx(1000);
    np45->SetMarkerStyle(kFullSquare);
    np45->SetMarkerColor(kBlack);
    np45->SetLineColor(kBlack);
    np45->GetXaxis()->SetLabelFont(132);
    np45->GetXaxis()->SetTitleFont(132);
    np45->GetYaxis()->SetLabelFont(132);
    np45->GetYaxis()->SetTitleFont(132);
    np45->GetXaxis()->SetTitleSize(.06);
    np45->GetYaxis()->SetTitleSize(.06);
    np45->GetXaxis()->SetLabelSize(.05);
    np45->GetYaxis()->SetLabelSize(.05);

    gStyle->SetOptStat(0);
    np45->Draw("ep");
    gPad->SetBottomMargin(0.2);
    np45->Fit(eff_func);*/
    
    

}
