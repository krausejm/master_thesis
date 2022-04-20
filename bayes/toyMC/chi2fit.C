void chi2fit(){
    TH1F* sigma = new TH1F("sigma",";#Sigma;counts",100,-1,1);
    TH1F* sigma_p45 = new TH1F("sigma_p45",";#Sigma;counts",100,-1,1);
    TH1F* sigma_m45 = new TH1F("sigma_m45",";#Sigma;counts",100,-1,1);

    TH1F* res = new TH1F("res",";#xi;counts",100,-5,5);

    int i=0;
    int nbins=24;
    //now create histos for the two settings
    TH1F* hp45 = new TH1F(Form("hp45_%d",i),";#phi / deg;",nbins,-180,180);
    TH1F* hm45 = new TH1F(Form("hm45_%d",i),";#phi / deg;",nbins,-180,180);
    //copies for error calculation
    TH1F* hp45e = new TH1F(Form("hp45e_%d",i),";#phi / deg;",nbins,-180,180);
    TH1F* hm45e = new TH1F(Form("hm45e_%d",i),";#phi / deg;",nbins,-180,180);
    //enum and nom for asym
    TH1F* nominator = new TH1F(Form("nom_%d",i),";#phi / deg;",nbins,-180,180);
	TH1F* enumerator = new TH1F(Form("enum_%d",i),";#phi / deg;",nbins,-180,180);
    //func for fitting event yields
    TF1* fp45 = new TF1("fp45","[0]*(1-0.3*[1]*cos(2*(45-x)*TMath::Pi()/180.))");
    TF1* fm45 = new TF1("fm45","[0]*(1-0.25*[1]*cos(2*(-45-x)*TMath::Pi()/180.))");
    //TF1* fm45 = new TF1("fm45","[0]*(1-0.3*[1]*cos(2*(-45-x)*TMath::Pi()/180.))");


    //tree for reading of data
    TTree* t = new TTree(Form("t%d",i),Form("mytree%d",i));

    for(int i=0;i<10000;i++){
        if(i%100==0) std::cout<<i<<std::endl;
        hp45->Reset();
        hm45->Reset();
        hp45e->Reset();
        hm45e->Reset();
        //fp45->Reset();
        //fm45->Reset();
        t->Reset();
        t->Refresh();
        //t->ReadFile(Form("./toybins/toybin%04d.txt",i),"",'\t');
        //t->ReadFile(Form("./test_toybins/toybin%04d.txt",i),"",'\t');
        t->ReadFile(Form("./pi0_toybins/toybin%04d.txt",i),"",'\t');
        //t->ReadFile(Form("./py_toybins/toybin%04d.txt",i),"pol:setting:phi",',');
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

        hp45e->Fit(fp45,"QN");
        hm45e->Fit(fm45,"QN");
        
        sigma_p45->Fill(fp45->GetParameter(1));
        sigma_m45->Fill(fm45->GetParameter(1));

        //normalize event yields
		double norm_p_err;
		double norm_m_err;

		double norm_p=hp45->IntegralAndError(0,-1,norm_p_err,"");
		double norm_m=hm45->IntegralAndError(0,-1,norm_m_err,"");

        /*norm_p=1;
        norm_m=1;
        norm_p_err=0;
        norm_m_err=0;*/

        hp45->Scale(1./norm_p);
        hm45->Scale(1./norm_m);
        //get enumerator and nominator of the asymmetry
		enumerator->Add(hp45,hm45,1,-1);
			
		hp45->Scale(0.25);
        //hp45->Scale(0.3);
		hm45->Scale(0.3);
		nominator->Add(hp45,hm45,1,1);
        enumerator->Divide(nominator);
        //gaussian error propagation
        double final_e[nbins];
		for(int k=0;k<nbins;k++){
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
        enumerator->Fit(f,"QN");
        float xi = (f->GetParameter(0)-0.3)/f->GetParError(0);
        res->Fill(xi);
        sigma->Fill(f->GetParameter(0));    

    
    
    }
    auto c1= new TCanvas("c1");
    TH1F* hlist[4]={sigma,res,sigma_p45,sigma_m45};
    for(int i=0;i<4;i++){
        hlist[i]->GetXaxis()->SetTitleFont(132);
        hlist[i]->GetXaxis()->SetLabelFont(132);
        hlist[i]->GetXaxis()->SetTitleSize(0.05);
        hlist[i]->GetXaxis()->SetLabelSize(0.04);
        
        hlist[i]->GetYaxis()->SetTitleFont(132);
        hlist[i]->GetYaxis()->SetLabelFont(132);
        hlist[i]->GetYaxis()->SetTitleSize(0.05);
        hlist[i]->GetYaxis()->SetLabelSize(0.04);

        gStyle->SetOptStat(0);
        gStyle->SetOptFit(1);
    }

    c1->Divide(2,2);
    c1->cd(1);
    sigma->Draw();
    sigma->Fit("gaus");
    c1->cd(2);
    res->Draw();
    res->Fit("gaus");
    c1->cd(3);
    sigma_p45->Draw("");
    sigma_p45->Fit("gaus");
    c1->cd(4);
    sigma_m45->Draw("");
    sigma_m45->Fit("gaus");





}