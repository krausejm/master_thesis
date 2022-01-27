void eta(TTree* t, TH3F* hp45, TH3F* hm45,TH2F* hp45pol, TH2F* hm45pol){
    double phi=0.;
    t->SetBranchAddress("phi",&phi);
    double ebeam=0.;
    t->SetBranchAddress("beam_energy",&ebeam);
    double costheta=0.;
    t->SetBranchAddress("cosinus_theta",&costheta);
    int prompt=0;
    t->SetBranchAddress("prompt",&prompt);
    double setting=0.;
    t->SetBranchAddress("setting",&setting);
    double pol=0.;
    t->SetBranchAddress("linear_polarisation",&pol);
    int whichparticle=0;
    t->SetBranchAddress("which_particle",&whichparticle);
    double w=0.;
    t->SetBranchAddress("weight",&w);
    double edge=0.;
    t->SetBranchAddress("edgecutoff",&edge);
    double peds=0.;
    t->SetBranchAddress("peds",&peds);
    for(int i=0;i<t->GetEntries();i++){
        t->GetEntry(i);
        //std::cout<<phi<<endl;
        if(whichparticle==1){
            if(setting==45){
                hp45->Fill(ebeam,costheta,phi,w);
                hp45pol->Fill(ebeam,pol,w);
            }
            if(setting==-45){
                hm45->Fill(ebeam,costheta,phi,w);
                hm45pol->Fill(ebeam,pol,w);
            }
        }
    }




}
void extract_histos(){
TFile* f= new TFile("farahs_trees.root","READ");
TFile* test = new TFile("./eta_histos.root","RECREATE");
TTree* tj= (TTree *)f->Get("sigma_analysis_july");
TTree* ta= (TTree *)f->Get("sigma_analysis_august");
TTree* ts= (TTree *)f->Get("sigma_analysis_september");
TTree* to= (TTree *)f->Get("sigma_analysis_october");
TTree* t[4]={tj,ta,ts,to};
TH3F* hp45[4];
TH3F* hm45[4];
TH2F* hp45pol[4];
TH2F* hm45pol[4];

for(int i=0;i<4;i++){
    if(i<2){
        hp45[i]=new TH3F(Form("yield_p45_%d",i),";E_{#gamma} / MeV;cos#theta;#phi/deg",9,1130,1670,12,-1,1,12,-180,180);
        hm45[i]=new TH3F(Form("yield_m45_%d",i),";E_{#gamma} / MeV;cos#theta;#phi/deg",9,1130,1670,12,-1,1,12,-180,180);
        hp45pol[i]=new TH2F(Form("p45_pol_deg_%d",i),";E_{#gamma} / MeV; Polarization degree",9,1130,1670,200,0,1);
        hm45pol[i]=new TH2F(Form("m45_pol_deg_%d",i),";E_{#gamma} / MeV; Polarization degree",9,1130,1670,200,-1,0);

    }else{
        hp45[i]=new TH3F(Form("yield_p45_%d",i),";E_{#gamma} / MeV;cos#theta;#phi/deg",9,1250,1790,12,-1,1,12,-180,180);
        hm45[i]=new TH3F(Form("yield_m45_%d",i),";E_{#gamma} / MeV;cos#theta;#phi/deg",9,1250,1790,12,-1,1,12,-180,180);
        hp45pol[i]=new TH2F(Form("p45_pol_deg_%d",i),";E_{#gamma} / MeV; Polarization degree",9,1250,1790,200,0,1);
        hm45pol[i]=new TH2F(Form("m45_pol_deg_%d",i),";E_{#gamma} / MeV; Polarization degree",9,1250,1790,200,-1,0);

    }
}
for(int i=0;i<4;i++){
    eta(t[i],hp45[i],hm45[i],hp45pol[i],hm45pol[i]);
    hp45[i]->Write();
    hm45[i]->Write();
    hp45pol[i]->Write();
    hm45pol[i]->Write();
}

test->Close();
}