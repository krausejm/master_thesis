Double_t loggaus(Double_t* x, Double_t* p){
    Double_t amp = p[0];
    Double_t mu = p[1];
    Double_t sigma= p[2];

    Double_t gaus= amp/TMath::Sqrt(2*TMath::Pi()*sigma*sigma)*TMath::Exp(-(x[0]-mu)*(x[0]-mu)/2/sigma/sigma);
    return TMath::Log(gaus);

}
Double_t gaus(Double_t* x, Double_t* p){
    Double_t amp = p[0];
    Double_t mu = p[1];
    Double_t sigma= p[2];

    Double_t gaus= amp/TMath::Sqrt(2*TMath::Pi()*sigma*sigma)*TMath::Exp(-(x[0]-mu)*(x[0]-mu)/2/sigma/sigma);
    return gaus;

}

void toyMC_posteriors(){
    //read sigma posteriors
    TTree* t = new TTree("t","mytree");
    TTree* d = new TTree("d","myothertree");
    t->ReadFile("toy_sigma.csv");
    d->ReadFile("mcse.csv","mcse");
    //set tree branch adresses for posteriors
    const int nbins = 300;
    Float_t currentry[nbins];
    TH1F* posteriors[nbins];
    //get mcse from other tree, mcse is second entry
    Float_t mcse[nbins];
    Float_t dummy;
    d->SetBranchAddress("mcse",&dummy);

    for(int i=0;i<nbins;i++){
        //create corresponding histos
        t->SetBranchAddress(Form("toybin%04d",i),&currentry[i]);
        posteriors[i]=new TH1F(Form("posterior%04d",i),";#Sigma;",100,-1,1);
        //get mcse
        d->GetEntry(i);
        mcse[i]=dummy;
        //std::cout<<mcse[i]<<std::endl;
    }
    //fill histos for each branch
    for(int i=0;i<t->GetEntries();i++){
        t->GetEntry(i);
        for(int j=0;j<nbins;j++){
            posteriors[j]->Fill((currentry[j]));
        }
    }
    //now convert each histo to log scale for better calculations
    for(int i=0;i<nbins;i++){
        for(int j=0;j<posteriors[i]->GetNbinsX();j++){
            //std::cout<<"hussa"<<std::endl;
            int bincontent=posteriors[i]->GetBinContent(j+1);
            if(bincontent!=0) posteriors[i]->SetBinContent(j+1,TMath::Log(bincontent));
        }
    }
    //multiply histos and Divide by prior
    TF1* f = new TF1("f",loggaus,-1,1,3);
    TF1* g = new TF1("f",gaus,-1,1,3);

    double mu = 0;
    double sigma =1;
    double amp = 1.;
    f->SetParameter(0,amp);
    f->SetParameter(1,mu);
    f->SetParameter(2,sigma);
    g->SetParameter(0,amp);
    g->SetParameter(1,mu);
    g->SetParameter(2,sigma);
    TH1F* combinedposterior = new TH1F("combined",";#Sigma;",100,-1,1);
    for(int i=0;i<nbins;i++){
        if(i==0){
            combinedposterior= (TH1F*) posteriors[i]->Clone();
        }else{
            combinedposterior->Add(posteriors[i]);
            //combinedposterior->Add(f,-1);
        }
    }
    combinedposterior->Scale(1./combinedposterior->Integral());
    for(int i=0;i<combinedposterior->GetNbinsX();i++){
        double tmp = TMath::Exp(combinedposterior->GetBinContent(i+1));
        combinedposterior->SetBinContent(i+1,tmp);
    }
    combinedposterior->Draw();
    //posteriors[0]->Draw();
    //std::cout<<g->Integral(-1e3,1e3)<<std::endl;;

}