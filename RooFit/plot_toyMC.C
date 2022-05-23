
void CanvasPartition(TCanvas *C,const Int_t Nx,const Int_t Ny,
                     Float_t lMargin, Float_t rMargin,
                     Float_t bMargin, Float_t tMargin)
{
   if (!C) return;
 
   // Setup Pad layout:
   Float_t vSpacing = 0.0;
   Float_t vStep  = (1.- bMargin - tMargin - (Ny-1) * vSpacing) / Ny;
 
   Float_t hSpacing = 0.0;
   Float_t hStep  = (1.- lMargin - rMargin - (Nx-1) * hSpacing) / Nx;
 
   Float_t vposd,vposu,vmard,vmaru,vfactor;
   Float_t hposl,hposr,hmarl,hmarr,hfactor;
 
   for (Int_t i=0;i<Nx;i++) {
 
      if (i==0) {
         hposl = 0.0;
         hposr = lMargin + hStep;
         hfactor = hposr-hposl;
         hmarl = lMargin / hfactor;
         hmarr = 0.0;
      } else if (i == Nx-1) {
         hposl = hposr + hSpacing;
         hposr = hposl + hStep + rMargin;
         hfactor = hposr-hposl;
         hmarl = 0.0;
         hmarr = rMargin / (hposr-hposl);
      } else {
         hposl = hposr + hSpacing;
         hposr = hposl + hStep;
         hfactor = hposr-hposl;
         hmarl = 0.0;
         hmarr = 0.0;
      }
 
      for (Int_t j=0;j<Ny;j++) {
 
         if (j==0) {
            vposd = 0.0;
            vposu = bMargin + vStep;
            vfactor = vposu-vposd;
            vmard = bMargin / vfactor;
            vmaru = 0.0;
         } else if (j == Ny-1) {
            vposd = vposu + vSpacing;
            vposu = vposd + vStep + tMargin;
            vfactor = vposu-vposd;
            vmard = 0.0;
            vmaru = tMargin / (vposu-vposd);
         } else {
            vposd = vposu + vSpacing;
            vposu = vposd + vStep;
            vfactor = vposu-vposd;
            vmard = 0.0;
            vmaru = 0.0;
         }
 
         C->cd(0);
 
         char name[16];
         sprintf(name,"pad_%i_%i",i,j);
         TPad *pad = (TPad*) gROOT->FindObject(name);
         if (pad) delete pad;
         pad = new TPad(name,"",hposl,vposd,hposr,vposu);
         pad->SetLeftMargin(hmarl);
         pad->SetRightMargin(hmarr);
         pad->SetBottomMargin(vmard);
         pad->SetTopMargin(vmaru);
 
         pad->SetFrameBorderMode(0);
         pad->SetBorderMode(0);
         pad->SetBorderSize(0);
 
         pad->Draw();
      }
   }


}



void plot_toyMC(){
	TFile* f = new TFile("./toyMC_results_new.root","READ");
	TH1F* hxi=(TH1F*) f->Get("res");
	TH1F* hsigma=(TH1F*) f->Get("hsigma");
	TH1F* hsigma_bkg=(TH1F*) f->Get("hsigma_bkg");

	TH1F* histos[3]={hxi,hsigma,hsigma_bkg};
	//auto res = hsigma->Fit("gaus","S");
	

	gStyle->SetOptStat(0);

   	TCanvas *C = (TCanvas*) gROOT->FindObject("C");
	if (C) delete C;
	C = new TCanvas("C","canvas",1024,640);
	C->SetFillStyle(4000);

	// Number of PADS
	const Int_t Nx = 1;
	const Int_t Ny = 1;

	TPad *pad[Nx][Ny];

	// Margins
	Float_t lMargin = 0.12;
	Float_t rMargin = 0.05;
	Float_t bMargin = 0.15;
	Float_t tMargin = 0.05;

	// Canvas setup
	/*CanvasPartition(C,Nx,Ny,lMargin,rMargin,bMargin,tMargin);
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++){
			C->cd(0);
			char pname[16];
			sprintf(pname,"pad_%i_%i",i,j);
			pad[i][j] = (TPad*) gROOT->FindObject(pname);
			pad[i][j]->Draw();
			pad[i][j]->SetFillStyle(4000);
			pad[i][j]->SetFrameFillStyle(4000);
			pad[i][j]->cd();


         		// Size factors
         		Float_t xFactor = pad[0][0]->GetAbsWNDC()/pad[i][j]->GetAbsWNDC();
         		Float_t yFactor = pad[0][0]->GetAbsHNDC()/pad[i][j]->GetAbsHNDC();

			// y axis range
         		histos[j]->GetYaxis()->SetRangeUser(0.0001,1.2*histos[i]->GetMaximum());
 
         		// Format for y axis
         		histos[j]->GetYaxis()->SetLabelFont(133);
        		histos[j]->GetYaxis()->SetLabelSize(18);
         		histos[j]->GetYaxis()->SetLabelOffset(0.02);
        	 	histos[j]->GetYaxis()->SetTitleFont(133);
         		histos[j]->GetYaxis()->SetTitleSize(20);
         		histos[j]->GetYaxis()->SetTitleOffset(2);
 
         		//histos[i]->GetYaxis()->CenterTitle();
         		histos[j]->GetYaxis()->SetNdivisions(505);
 
         		// TICKS Y Axis
         		histos[j]->GetYaxis()->SetTickLength(xFactor*0.04/yFactor);
 
         		// Format for x axis
         		histos[j]->GetXaxis()->SetLabelFont(133);
         		histos[j]->GetXaxis()->SetLabelSize(18);
         		histos[j]->GetXaxis()->SetLabelOffset(0.02);
         		histos[j]->GetXaxis()->SetTitleFont(133);
         		histos[j]->GetXaxis()->SetTitleSize(20);
         		histos[j]->GetXaxis()->SetTitleOffset(1.5);
        	 	//histos[i]->GetXaxis()->CenterTitle();
         		histos[j]->GetXaxis()->SetNdivisions(505);
			if(i==0)histos[i]->GetXaxis()->SetRangeUser(-5,5);
			//else histos[i]->GetXaxis()->SetRangeUser(-1.2,1.2);
 
         		// TICKS X Axis
         		histos[j]->GetXaxis()->SetTickLength(yFactor*0.06/xFactor);
			histos[j]->Draw();
		
		}
	}*/
	double x,y;
	TLatex t;
	t.SetTextAlign(22);
	const char * names[3] = {"residuals","sigma","sigma_bkg"};
	for(int i=0;i<3;i++){
		auto c = new TCanvas(Form("c%d",i));
		c->cd(0);
		histos[i]->Rebin(2);
		histos[i]->Draw("");
		auto res = histos[i]->Fit("gaus","SQ");
		histos[i]->GetXaxis()->SetLabelFont(133);
		histos[i]->GetXaxis()->SetTitleFont(133);
		histos[i]->GetXaxis()->SetLabelSize(18);
                histos[i]->GetXaxis()->SetTitleSize(20);
		histos[i]->GetYaxis()->SetLabelFont(133);
                histos[i]->GetYaxis()->SetTitleFont(133);
                histos[i]->GetYaxis()->SetLabelSize(18);
                histos[i]->GetYaxis()->SetTitleSize(20);
		histos[i]->GetYaxis()->SetRangeUser(0,1.3*res->Parameter(0));
		if(i==2)histos[i]->GetXaxis()->SetTitle("#Sigma_{bkg}");
		x=0;
		y=1.2*res->Parameter(0);
		//std::cout<<histos[i]->GetMaximum();
		//y=1;
		t.DrawLatex(x,y,Form("#color[2]{#font[132]{#mu=%.4f#pm%.4f, #sigma=%.4f#pm%.4f}}",res->Parameter(1),res->ParError(1),res->Parameter(2),res->ParError(2)));
		c->SaveAs(Form("./plots/%s.pdf",names[i]));







		




	}

}
