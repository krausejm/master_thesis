using namespace std;
#include <iomanip>      // std::setprecision
#include "TMath.h"
#include "TFile.h"
#include "TF1.h"
#include "TString.h"
#include "TF3.h"
#include "TGraphErrors.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TLine.h"
#include <iostream> 
#include <fstream> 
#include "TPad.h"
#include "TLatex.h"
#include "TTree.h"
#include "TRandom3.h"



double pdf(double* x, double* p){
    return (1-p[0]*p[1]*TMath::Cos(TMath::DegToRad()*2*(p[2]-x[0])));     
}




void make_bin(int k, double p45pol, double m45pol, double sigma, int e1, int e2, int file){
//this is the pdf 
TF1* mypdf=new TF1("mypdf",pdf,-180,180,3);
mypdf->SetNpx(360);
//this is the file data will be written to
ofstream myfile;
if(file==0){myfile.open (Form("./toybins/toybin%04d.txt",k));}
if(file==1){myfile.open (Form("./test_toybins/toybin%04d.txt",k));}
if(file==2){myfile.open (Form("./pi0_toybins/toybin%04d.txt",k));}
if(file==3){myfile.open (Form("./etap_toybins/toybin%04d.txt",k));}

myfile << "pol\tsetting\tphi\n";
TRandom3 r(0);
gRandom->SetSeed(0);
int n1=r.Poisson(e1);
int n2=r.Poisson(e2);
cout<<Form("throwing toy MC experiment no.%d in setting %d",k,file)<<endl;
cout<<n1<<endl;
cout<<n2<<endl;
//generate events for bot config
for(int i=0;i<n1;i++){
    mypdf->FixParameter(0,sigma);//sigma
    mypdf->FixParameter(1,p45pol);//pol
    mypdf->FixParameter(2,+45);//setting
    myfile<<p45pol<<"\t"<<+45<<"\t"<<mypdf->GetRandom()<<"\n";
}
//generate events for parallel config
for(int i=0;i<n2;i++){
    mypdf->FixParameter(0,sigma);//sigma
    mypdf->FixParameter(1,m45pol);//pol
    mypdf->FixParameter(2,-45);//setting
    myfile<<m45pol<<"\t"<<-45<<"\t"<<mypdf->GetRandom()<<"\n";
}
}
void new_toy_MC(){
    for(int k=0;k<10000;k++){
        //eta bins
	make_bin(k,0.3,0.3,0.3,200,200,3);
	//test bins
	//make_bin(k,0.3,0.3,0.3,1000,1000,1);
	//pi0 bins
	//make_bin(k,0.3,0.25,0.3,5000,4000,2);


    }


}
