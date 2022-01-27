{
TFile* file = new TFile("eta_histos.root","READ");
TH3F* hp45_0= (TH3F* ) file->Get("yield_p45_0");
TH3F* hp45_1= (TH3F* ) file->Get("yield_p45_1");
TH3F* hp45_2= (TH3F* ) file->Get("yield_p45_2");
TH3F* hp45_3= (TH3F* ) file->Get("yield_p45_3");

TH3F* hm45_0= (TH3F* ) file->Get("yield_m45_0");
TH3F* hm45_1= (TH3F* ) file->Get("yield_m45_1");
TH3F* hm45_2= (TH3F* ) file->Get("yield_m45_2");
TH3F* hm45_3= (TH3F* ) file->Get("yield_m45_3");

TH2F* hp45pol_0= (TH2F* )file->Get("p45_pol_deg_0");
TH2F* hp45pol_1= (TH2F* )file->Get("p45_pol_deg_1");
TH2F* hp45pol_2= (TH2F* )file->Get("p45_pol_deg_2");
TH2F* hp45pol_3= (TH2F* )file->Get("p45_pol_deg_3");

TH2F* hm45pol_0= (TH2F* )file->Get("m45_pol_deg_0");
TH2F* hm45pol_1= (TH2F* )file->Get("m45_pol_deg_1");
TH2F* hm45pol_2= (TH2F* )file->Get("m45_pol_deg_2");
TH2F* hm45pol_3= (TH2F* )file->Get("m45_pol_deg_3");

hp45_0->Add(hp45_1);
hp45_2->Add(hp45_3);

hm45_0->Add(hm45_1);
hm45_2->Add(hm45_3);

hp45pol_0->Add(hp45pol_1);
hp45pol_2->Add(hp45pol_3);

hm45pol_0->Add(hm45pol_1);
hm45pol_2->Add(hm45pol_3);

TH1D* hp45pol_tmp = hp45pol_0->ProjectionY("ppol0_py",9,9);
hp45pol_tmp->Add(hp45pol_2->ProjectionY("ppol2_py",7,7));
TH1D* hm45pol_tmp = hm45pol_0->ProjectionY("mpol0_py",9,9);
hm45pol_tmp->Add(hm45pol_2->ProjectionY("mpol0_py",7,7));

double polp45 = hp45pol_tmp->GetMean();
double polm45 = -1*hm45pol_tmp->GetMean();

auto c = new TCanvas();
c->Divide(4,3);
for (int j=0;j<12;j++){
    TH1D* hp_tmp = hp45_0->ProjectionZ(Form("p45_0_pz%d",j),9,9,j+1,j+1);
    TH1D* hm_tmp = hm45_0->ProjectionZ(Form("m45_0_pz%d",j),9,9,j+1,j+1);
    //add other months
    hp_tmp->Add(hp45_2->ProjectionZ(Form("p45_2_pz%d",j),7,7,j+1,j+1));
    hm_tmp->Add(hm45_2->ProjectionZ(Form("m45_2_pz%d",j),7,7,j+1,j+1));

    hp_tmp->Scale(1./hp_tmp->Integral());
    hm_tmp->Scale(1./hm_tmp->Integral());
    TH1D* enumr = (TH1D*) hp_tmp->Clone();
    enumr->Add(hm_tmp,-1);
    TH1D* nom = (TH1D*) hp_tmp->Clone();
    nom->Add(hp_tmp,hm_tmp,polm45,polp45);
    enumr->Divide(nom);
    c->cd(j+1);
    enumr->Draw();



}



    
}