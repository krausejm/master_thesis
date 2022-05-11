#include <cstdio>
#include <iostream>




void make_files(){
    for(int i=0;i<11;i++){//energy bins
        for(int j=0;j<12;j++){//costheta bins
            ofstream myfile;
            char buffer[32]; // The filename buffer.
            snprintf(buffer, sizeof(char) * 32, "ebin%02dcostbin%02d.txt", i,j);
            myfile.open(buffer);
            myfile<< "pol\tphi\tweight\n";
            myfile.close();
    }
}

}

void fill_files(TTree* tj,TTree* ta,TTree* ts,TTree* to){
//loop over all trees
TTree* t[4]={tj,ta,ts,to};

for(int i=0;i<4;i++){
    //set branch adresses of tree vars    
    double phi=0.;
    t[i]->SetBranchAddress("phi",&phi);
    double ebeam=0.;
    t[i]->SetBranchAddress("beam_energy",&ebeam);
    double costheta=0.;
    t[i]->SetBranchAddress("cosinus_theta",&costheta);
    int prompt=0;
    t[i]->SetBranchAddress("prompt",&prompt);
    double setting=0.;
    t[i]->SetBranchAddress("setting",&setting);
    double pol=0.;
    t[i]->SetBranchAddress("linear_polarisation",&pol);
    int whichparticle=0;
    t[i]->SetBranchAddress("which_particle",&whichparticle);
    double w=0.;
    t[i]->SetBranchAddress("weight",&w);
    double edge=0.;
    t[i]->SetBranchAddress("edgecutoff",&edge);
    double peds=0.;
    t[i]->SetBranchAddress("peds",&peds);
    for(int k=0;k<11;k++){//energy bins
        if((i==0||i==1)&&k>8){continue;}//respect the coherent edges
        if((i==2||i==3)&&k<2){continue;}
        for(int l=0;l<12;l++){//angle bins
            std::cerr<<"tree: "<<i<<"ebin: "<<k<<" costbin: "<<l<<"\n";
            char buffer[32]; // The filename buffer.
            snprintf(buffer, sizeof(char) * 32, "ebin%02dcostbin%02d.txt", k,l);
            ofstream myfile;
            myfile.open(buffer,std::ios::app);
            
            
            for(int j=0;j<t[i]->GetEntries();j++){//tree entries
                t[i]->GetEntry(j);
                //take only eta events
                if(whichparticle==1){
                    //loop over all possible bins
                    if(1130+k*60<=ebeam&&ebeam<1130+(k+1)*60){
                    if(-1+l*1./6.<=costheta&&costheta<-1+(l+1)*1./6.){
                    myfile<<pol<<"\t"<<phi<<"\t"<<w<<"\n";
                    }  
                    }
                }
            
            }
            myfile.close();
            


        }
    }

}

}


void create_csv(){
TFile* f= new TFile("../realdeal/farahs_trees.root","READ");
TTree* tj= (TTree *)f->Get("sigma_analysis_july");
TTree* ta= (TTree *)f->Get("sigma_analysis_august");
TTree* ts= (TTree *)f->Get("sigma_analysis_september");
TTree* to= (TTree *)f->Get("sigma_analysis_october");
int count=0;
make_files();
fill_files(tj,ta,ts,to);
}

