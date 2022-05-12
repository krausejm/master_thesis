#include <cstdio>
#include <iostream>




void make_files(){
    for(int i=0;i<3;i++){//energy bins
        for(int j=0;j<6;j++){//costheta bins
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
    t[i]->SetBranchAddress("egamma",&ebeam);
    double costheta=0.;
    t[i]->SetBranchAddress("costheta",&costheta);
    double setting=0.;
    t[i]->SetBranchAddress("setting",&setting);
    double pol=0.;
    t[i]->SetBranchAddress("pol",&pol);
    double w=0.;
    t[i]->SetBranchAddress("weight",&w);
    double peds=0.;
    t[i]->SetBranchAddress("PED",&peds);
    for(int k=0;k<3;k++){//energy bins
        if((i==0||i==1)&&k>1){continue;}//respect the coherent edges
        for(int l=0;l<6;l++){//angle bins
            std::cerr<<"tree: "<<i<<"ebin: "<<k<<" costbin: "<<l<<"\n";
            char buffer[32]; // The filename buffer.
            snprintf(buffer, sizeof(char) * 32, "ebin%02dcostbin%02d.txt", k,l);
            ofstream myfile;
            myfile.open(buffer,std::ios::app);
            
            
            for(int j=0;j<t[i]->GetEntries();j++){//tree entries
                t[i]->GetEntry(j);
                //loop over all possible bins
                if(1500+k*100<=ebeam&&ebeam<1500+(k+1)*100){
                if(-1+l*1./3.<=costheta&&costheta<-1+(l+1)*1./3.){
                	myfile<<pol<<"\t"<<phi<<"\t"<<w<<"\n";
                 }  
                 }
             
            
            }
            myfile.close();
            


        }
    }

}

}


void create_csv(){
TFile* f= new TFile("/hiskp3/krause/test/mytrees.root","READ");
TTree* tj= (TTree *)f->Get("mytree_july");
TTree* ta= (TTree *)f->Get("mytree_august");
TTree* ts= (TTree *)f->Get("mytree_september");
TTree* to= (TTree *)f->Get("mytree_october");
int count=0;
make_files();
fill_files(tj,ta,ts,to);
}

