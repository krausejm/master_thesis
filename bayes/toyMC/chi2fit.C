void chi2fit(){
    TH1F* sigma = new TH1F("sigma",";#xi;counts",100,-5,5);
    int i=0;
    int nbins=12;
    //now create histos for the two settings
    TH1F* hp45 = new TH1F(Form("hp45_%d",i),";#phi / deg;",nbins,-180,180);
    TH1F* hm45 = new TH1F(Form("hm45_%d",i),";#phi / deg;",nbins,-180,180);
    //copies for error calculation
    TH1F* hp45e = new TH1F(Form("hp45e_%d",i),";#phi / deg;",nbins,-180,180);
    TH1F* hm45e = new TH1F(Form("hm45e_%d",i),";#phi / deg;",nbins,-180,180);
    //enum and nom for asym
    TH1F* nominator = new TH1F(Form("nom_%d",i),";#phi / deg;",nbins,-180,180);
	TH1F* enumerator = new TH1F(Form("enum_%d",i),";#phi / deg;",nbins,-180,180);



    //tree for reading of data
    TTree* t = new TTree(Form("t%d",i),Form("mytree%d",i));

    for(int i=0;i<10000;i++){
        hp45->Reset();
        hm45->Reset();
        hp45e->Reset();
        hm45e->Reset();
        t->Reset();
        t->Refresh();
        t->ReadFile(Form("./toybins/toybin%04d.txt",i),"",'\t');
        float pol, setting, phi;
        t->SetBranchAddress("pol",&pol);
        t->SetBranchAddress("setting",&setting);
        t->SetBranchAddress("phi",&phi);

        for(int j=0;j<t->GetEntries();j++){
            t->GetEntry(j);
            if(setting==+45){
                hp45->Fill(phi);
                hp45e->Fill(phi);
            }else{
                hm45->Fill(phi);
                hm45e->Fill(phi);
            }
        }
        //normalize event yields
		double norm_p_err;
		double norm_m_err;

		double norm_p=hp45->IntegralAndError(0,-1,norm_p_err,"");
		double norm_m=hm45->IntegralAndError(0,-1,norm_m_err,"");

        hp45->Scale(1./norm_p);
        hm45->Scale(1./norm_m);
        //get enumerator and nominator of the asymmetry
		enumerator->Add(hp45,hm45,1,-1);
			
		hp45->Scale(0.25);
		hm45->Scale(0.3);
		nominator->Add(hp45,hm45,1,1);
        enumerator->Divide(nominator);
        //gaussian error propagation
        double final_e[12];
		for(int k=0;k<12;k++){
			double n_bot = hp45e->GetBinContent(k+1);
			double n_par = hm45e->GetBinContent(k+1);
			double n_bot_err=hp45e->GetBinError(k+1);
			double n_par_err=hm45e->GetBinError(k+1);
			double pol_bot = 0.3;
			double pol_par = 0.25;
			double tmp_nom = TMath::Power(pol_par*n_bot*norm_m+pol_bot*n_par*norm_p,4);
			double tmp_enum = TMath::Power(n_bot*n_par*norm_m*(pol_bot+pol_par)*norm_p_err,2)
							+TMath::Power(n_bot*n_par*norm_p*(pol_bot+pol_par)*norm_m_err,2)
							+TMath::Power(n_par*norm_p*norm_m*(pol_bot+pol_par)*n_bot_err,2)
							+TMath::Power(n_bot*norm_p*norm_m*(pol_bot+pol_par)*n_par_err,2);

			final_e[k]=TMath::Sqrt(tmp_enum/tmp_nom);	
			//std::cout<<tmp_nom<<std::endl;
			enumerator->SetBinError(k+1,final_e[k]);
		}
        //enumerator->Draw("ep");
        TF1* f = new TF1(Form("f%d",i),"[0]*cos(2*(-45-x)*TMath::Pi()/180.)",-180,180);
        enumerator->Fit(f,"Q");
        float xi = (f->GetParameter(0)-0.3)/f->GetParError(0);
        sigma->Fill(xi);
        //sigma->Fill(f->GetParameter(0));    

    
    
    }
    sigma->Draw();
    sigma->Fit("gaus");





}