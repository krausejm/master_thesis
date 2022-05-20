//simple 3d fit
double func3(Double_t *x, Double_t *par){
    Double_t r2 = x[0]*x[0] + x[1]*x[1] +x[2]*x[2];
    Double_t arg = (TMath::Sqrt(r2) - par[1])/par[2];
    Double_t val = par[0]*TMath::Exp(-0.5*arg*arg);
    return val;
}



void testing(){
    TF3 *f3 = new TF3("f3",func3,-5,5,-5,5,-5,5,3);
    f3->SetParameters(10,0,1);
    TH3F *h3 = new TH3F("h3","h3",50,-5,5,50,-5,5,50,-5,5);
    Double_t x,y,z;
    for (Int_t i=0;i<10000000;i++) {
        f3->GetRandom3(x,y,z);
        h3->Fill(x,y,z);
    }
    h3->Draw("box");
    h3->Fit(f3);
    //f3->Draw("");
}