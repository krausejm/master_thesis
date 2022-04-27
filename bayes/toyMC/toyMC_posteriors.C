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
    t->ReadFile("sigma.csv");
    //d->ReadFile("mcse.csv","mcse");
    //set tree branch adresses for posteriors
    const int nbins = 300;
    Float_t currentry[nbins];
    TH1F* posteriors[nbins];
    //get mcse from other tree, mcse is second entry
    Float_t mcse[nbins];
    Float_t dummy;
    //d->SetBranchAddress("mcse",&dummy);

    for(int i=0;i<nbins;i++){
        //create corresponding histos
        t->SetBranchAddress(Form("toybin%04d",i),&currentry[i]);
        posteriors[i]=new TH1F(Form("posterior%04d",i),";#Sigma;",150,-1,1);
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
            int bincontent=posteriors[i]->GetBinContent(j+1);
            posteriors[i]->SetBinContent(j+1,TMath::Log(bincontent));

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
    TH1F* combinedposterior = new TH1F("combined",";#Sigma;",150,-1,1);
    //multiplication via log scale addition
    for(int i=0;i<nbins;i++){
        if(i==0){
            combinedposterior= (TH1F*) posteriors[i]->Clone();
        }else{
            combinedposterior->Add(posteriors[i]);
            combinedposterior->Add(f,-1);

        }
    }
    TH1F* test = new TH1F("test",";#Sigma;",150,-1,1);
    for(int i=0;i<combinedposterior->GetNbinsX();i++){
        double tmp;
        std::cout<<"lp="<<combinedposterior->GetBinContent(i+1)<<std::endl;
        //convert back to linear scale
        if(combinedposterior->GetBinContent(i+1)>0) tmp = TMath::Exp(combinedposterior->GetBinContent(i+1)-1650);
        else tmp = TMath::Exp(combinedposterior->GetBinContent(i+1));
        std::cout<<"exp(lp)="<<tmp<<std::endl;
        combinedposterior->SetBinContent(i+1,tmp);
        test->SetBinContent(i+1,tmp);

    }
    //divide by prior
    //test->Divide(g,1);
    //fit gaus because
    TF1* fitf = new TF1("fitf","gaus",0.2,0.4);
    fitf->SetNpx(1e5);
    test->Draw("");

    test->Fit(fitf);
    //cosmetics
    gStyle->SetOptStat(0);
    test->GetXaxis()->SetTitleFont(132);
    test->GetXaxis()->SetLabelFont(132);
    test->GetYaxis()->SetTitleFont(132);
    test->GetYaxis()->SetLabelFont(132);
    //test->GetYaxis()->SetLabelSize(0);
    test->GetXaxis()->SetRangeUser(-1,1);
    //test->GetYaxis()->SetNdivisions(1);
    //test->GetYaxis()->SetRangeUser(0,100);
    test->GetYaxis()->SetTitle("#it{p}(#Sigma|D) (arb. units)");
    double fmu, fmu_err, fsigma, fsigma_err;
    fmu=fitf->GetParameter(1);
    fmu_err=fitf->GetParError(1);
    fsigma=fitf->GetParameter(2);
    fsigma_err=fitf->GetParError(2);
    TLatex text;
    text.SetTextAlign(22);
    text.DrawLatex(0.5,95,Form("#font[132]{#color[2]{#mu=%.4f#pm%.4f, #sigma=%.4f#pm %.4f}}",fmu,fmu_err,fsigma,fsigma_err));

}