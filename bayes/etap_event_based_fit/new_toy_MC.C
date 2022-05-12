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
    double eff=1./10.5*(9.3+0.28*cos(x[0]*TMath::DegToRad())+0.24*sin(3*x[0]*TMath::DegToRad()));
    return eff*(1-p[0]*p[1]*TMath::Cos(TMath::DegToRad()*2*(p[2]-x[0])));     
}




void make_bin(int k){
//this is the pdf 
TF1* mypdf=new TF1("mypdf",pdf,-180,180,3);
//this is the file data will be written to
ofstream myfile;
myfile.open (Form("./toybins/toybin%04d.txt",k));
myfile << "pol.\tphi\tweight\n";

TRandom3 r(0);
gRandom->SetSeed(0);
double p45pol=0.300000;
double m45pol=0.300000;
double sigma=0.5;
double sigma_bkg=-0.5;
double weight=0.065;
int n1=r.Poisson(200);//statistics of final state, n1=p45, n2=m45
int n2=r.Poisson(200);
cout<<Form("throwing toy MC experiment no.%d",k)<<endl;
cout<<n1<<endl;
cout<<n2<<endl;
double f = 0.95; // fraction of signal events in prmpt peak

for(int i=0;i<1;i++){// 1 iteration for each time cut
    //signal events for first setting
    for(int j=0;j<int(n1*f);j++){
            mypdf->FixParameter(0,sigma);//sigma
            mypdf->FixParameter(1,p45pol);//pol
            mypdf->FixParameter(2,+45);//setting
            myfile<<p45pol<<"\t"<<mypdf->GetRandom()<<"\t"<<1<<"\n";
    }
    //bkg events for first setting
    for(int j=0;j<int(n1*(1-f));j++){
            mypdf->FixParameter(0,sigma_bkg);//sigma bkg
            mypdf->FixParameter(1,p45pol);//pol
            mypdf->FixParameter(2,+45);//setting
            myfile<<p45pol<<"\t"<<mypdf->GetRandom()<<"\t"<<1<<"\n";
    }
    //sideband events for first setting
    for(int j=0;j<int(n1*(1-f)*1./weight);j++){
            mypdf->FixParameter(0,sigma_bkg);//sigma bkg
            mypdf->FixParameter(1,p45pol);//pol
            mypdf->FixParameter(2,+45);//setting
            myfile<<p45pol<<"\t"<<mypdf->GetRandom()<<"\t"<<-1*weight<<"\n";
    }
    ////////////////////////////////////////////////////////////////////////////
    //signal events for second setting
    for(int j=0;j<int(n2*f);j++){
            mypdf->FixParameter(0,sigma);//sigma
            mypdf->FixParameter(1,m45pol);//pol
            mypdf->FixParameter(2,-45);//setting
            myfile<<-1*m45pol<<"\t"<<mypdf->GetRandom()<<"\t"<<1<<"\n";
    }
    //bkg events for first setting
    for(int j=0;j<int(n2*(1-f));j++){
            mypdf->FixParameter(0,sigma_bkg);//sigma bkg
            mypdf->FixParameter(1,m45pol);//pol
            mypdf->FixParameter(2,-45);//setting
            myfile<<-1*m45pol<<"\t"<<mypdf->GetRandom()<<"\t"<<1<<"\n";
    }
    //sideband events for first setting
    for(int j=0;j<int(n2*(1-f)*1./weight);j++){
            mypdf->FixParameter(0,sigma_bkg);//sigma bkg
            mypdf->FixParameter(1,m45pol);//pol
            mypdf->FixParameter(2,-45);//setting
            myfile<<-1*m45pol<<"\t"<<mypdf->GetRandom()<<"\t"<<-1*weight<<"\n";
    }
}

}
void new_toy_MC(){
    for(int k=0;k<1000;k++){
        make_bin(k);
    }


}