void toyMC_posteriors(){
    //read sigma posteriors
    TTree* t = new TTree("t","mytree");
    TTree* d = new TTree("d","myothertree");
    t->ReadFile("toy_sigma.csv");
    d->ReadFile("mcse.csv","mcse");
    //set tree branch adresses for posteriors
    const int nbins = 100;
    Float_t currentry[nbins];
    TH1F* posteriors[nbins];
    //get mcse from other tree, mcse is second entry
    Float_t mcse[nbins];
    Float_t dummy;
    d->SetBranchAddress("mcse",&dummy);

    for(int i=0;i<nbins;i++){
        //create corresponding histos
        t->SetBranchAddress(Form("toybin%04d",i),&currentry[i]);
        posteriors[i]=new TH1F(Form("posterior%04d",i),";#Sigma;",50,-1,1);
        //get mcse
        d->GetEntry(i);
        mcse[i]=dummy;
        //std::cout<<mcse[i]<<std::endl;
    }
    //fill histos for each branch
    for(int i=0;i<t->GetEntries();i++){
        t->GetEntry(i);
        for(int j=0;j<nbins;j++){
            posteriors[j]->Fill((0.5-currentry[j]));
        }
    }
    //multiply histos and Divide by prior
    TF1* f = new TF1("f","gaus",-1,1);
    double mu = 0;
    double sigma =1;
    double amp = 1./TMath::Sqrt(2*TMath::Pi()*sigma*sigma);
    f->SetParameter(0,amp);
    f->SetParameter(1,mu);
    f->SetParameter(2,sigma);
    TH1F* combinedposterior = new TH1F("combined",";#Sigma;",50,-1,1);
    for(int i=0;i<nbins;i++){
        if(i==0){
            combinedposterior= (TH1F*) posteriors[i]->Clone();
        }else{
            combinedposterior->Multiply(posteriors[i]);
            combinedposterior->Divide(f,300000*mcse[i]);
        }
    }
    combinedposterior->Draw();
    //posteriors[0]->Draw();

}