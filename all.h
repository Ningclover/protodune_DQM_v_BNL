#include <chrono>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

#include "canvas/Utilities/InputTag.h"
#include "gallery/Event.h"

#include "TFile.h"
#include "TH1F.h"
#include "TInterpreter.h"
#include "TROOT.h"

const int nticks = 5859;
const int nchans = 10240;

        int runNo;
        int evtNo;
        vector<int> evtNo_list;
TString time_string;


TH1I *h_baseline = new TH1I("h_baseline", "Pedestal", nchans, -0.5, nchans - 0.5);
TH2F *h_orig = new TH2F("h_orig", "RawDigits", nchans, -0.5, nchans - 0.5, nticks, 0, nticks);
TH2F *h_orig_sep[4][3];
TH2F *h_orig_ANf = new TH2F("h_orig_ANf", "RawDigits After Nf", nchans, -0.5, nchans - 0.5, nticks, 0, nticks);
TH2F *h_orig_ANf_sep[4][3];
// TH1F *h_rms = new TH1F("h_rms", "raw_rms", nchans, -0.5, nchans - 0.5);
// TH1F *h_rms_ANf = new TH1F("h_rms_ANf", "raw_rms After Nf", nchans, -0.5, nchans - 0.5);
TH1F *h_rms;
TH1F *h_rms_ANf;
TH1F *hfft[4][3];
TH1F *hfft_ANf[4][3];
TH2F *h_cov[4][3];
TH2F *h_cov_ANf[4][3];
// void reset_all_hist(){
//     h_baseline->Reset();

// }
const char* plane_name[3] ={"u","v","w"}; 
int channels[] = {0,800,1600,2560};
int APA_number[] ={1,2,3,4};

double channel_start[4][3]; 
double channel_end[4][3];

TCanvas *c0 = new TCanvas("c0", "c0", 1);

float get_rms(TH1F* h1){

	if(h1->GetEntries()==0) return 0;
	double mean = h1->GetSum()/h1->GetNbinsX();
	double rms = 0;
	int nbin = h1->GetNbinsX();
	for (int j=0;j!=nbin;j++){
		rms += pow(h1->GetBinContent(j+1)-mean,2);
	}
	rms = sqrt(rms/h1->GetNbinsX());
	//   cout<<"rms= "<<rms<<endl;

	TH1F h2("h2","h2",100,mean-10*rms,mean+10*rms);
	for (int j=0;j!=nbin;j++){
		if (fabs(h1->GetBinContent(j+1)-mean)<6*rms)
			h2.Fill(h1->GetBinContent(j+1));
	}
	if(h2.GetEntries()==0) return 0;

	double par[3];
	double xq = 0.5;
	h2.GetQuantiles(1,&par[1],&xq);
	xq = 0.5 + 0.34;
	h2.GetQuantiles(1,&par[0],&xq);
	xq = 0.5 - 0.34;
	h2.GetQuantiles(1,&par[2],&xq);
	mean = par[1];
	rms = std::sqrt((pow(par[0]-par[1],2)+pow(par[2]-par[1],2))/2.);

	// remove signal
	for (int j=0;j!=nbin;j++){
		if(fabs(h1->GetBinContent(j+1) - mean) > 3.5*rms){
			h1->SetBinContent(j+1, mean);
		}
	}
	h2.Delete();
	
	return rms;

}

TH1F *get_rms_hist(TH2F *h2){
	int Nx = h2->GetNbinsX();
	TH1F *h_rms = new TH1F("hrms","RMS",Nx,-0.5,Nx-0.5);
	for(int i = 1; i <= Nx; i++){
		TH1F* h1 = (TH1F*)h2->ProjectionY("proj", i, i);
		double rms = get_rms(h1);
		h_rms->SetBinContent(i,rms);
	}
	return h_rms;
}

void FillRange(vector<int>& v, int start, int end){
	for(int i=start; i<end; i++){
		v.push_back(i);
	}
}



TH1F *CombineFFT(TH2F * hu,int APA = 1, const char* plane="u",TString tag = "ori", int color = 1){

	vector<int> vchns;
	APA = APA-1;
	//FillRange(vchns, 1904, 3072); // collection
	if(std::strcmp(plane, "u") == 0){
		FillRange(vchns, APA*2560, APA*2560+800); // u plane
	}
	else if(std::strcmp(plane, "v") == 0){
		FillRange(vchns, APA*2560+800, APA*2560+1600);
	}
	else if(std::strcmp(plane, "w") == 0){
		FillRange(vchns, APA*2560+1600, APA*2560+2560);
	}

	const int nticks = hu->GetNbinsY();
	const int nchs = hu->GetNbinsX();

	TH1F* hfft= new TH1F(Form("noise_%d%s_",APA+1,plane)+tag,"noise", nticks,0,1./0.512); // MHz, 512ns period
	hfft->SetTitle(Form("FFT -- APA%d %s plane", APA+1,plane));
	hfft->GetXaxis()->SetTitle("Frequency (MHz)");
	int nfft=0; // number of non-empty channels
	// TCanvas *c2 = new TCanvas("c2","",800,600);
	//    c2->SaveAs("/nashome/x/xning/Pictures/test.pdf[");
	for(auto channel: vchns){

		TH1F *h2 = new TH1F(Form("Channel%d", channel),Form("Channel %d", channel),nticks,0,nticks);
		for (int i=0;i!=nticks;i++){
			h2->SetBinContent(i+1,hu->GetBinContent(channel+1,i+1));
		}

		// ExcludeSignal(h2);

		float rms = get_rms(h2);
		TH1F* hmag = (TH1F*)h2->FFT(0, "MAG");
		//cout<<"bin number in hmag : "<<hmag->GetNbinsX()<<endl;
		int dc = hmag->GetBinContent(1);// dc component indicates empty or not
		// if(hmag->GetBinContent(2)<5E3){ // 5E3 avoids the extreme case
		// if(rms<30 and rms>0){
		if(1){
			nfft++;
			for(int i=1; i<nticks; i++){
				double content = hfft->GetBinContent(i+1);
				hfft->SetBinContent(i+1, content + hmag->GetBinContent(i+1));
				//cout<<i+1<<"  "<<content + hmag->GetBinContent(i+1)<<endl;
			}
		}
		//        cout<<"channel = "<<channel<<endl; 
		//    c2->cd();
		//    hmag->Draw();
		//    c2->SaveAs("/nashome/x/xning/Pictures/test.pdf");
		hmag->Delete();
		h2->Delete();
	}
	//c2->SaveAs("/nashome/x/xning/Pictures/test.pdf]");
	//cout << "nfft= " << nfft << endl;

	hfft->Scale(1./nfft);
	hfft->Scale(1400.0/(4096*4)); // convert ADC to mV (TDE)
	hfft->Rebin(27);
	hfft->GetXaxis()->SetRange(0, 4200); // X-axis range limited to 1 MHz
	hfft->GetXaxis()->SetRangeUser(0, 1); // X-axis range limited to 1 MHz
	hfft->SetLineColor(color);
	return hfft;

	}



TH2F* CovMatrix(TH2F *hu, const char* plane, int APA, TString tag = "ori"){   //   how to run ---> root.exe 'CovMatrix(ch_initial, ch_final)'   

	int ch_i, ch_f;
	APA = APA-1;
	if(std::strcmp(plane, "u") == 0){
		ch_i =APA*2560;
		ch_f = APA*2560+800;
	}
	else if(std::strcmp(plane, "v") == 0){
		ch_i =APA*2560+800;
		ch_f = APA*2560+1600;
	}
	else if(std::strcmp(plane, "w") == 0){
		ch_i =APA*2560+1600;
		ch_f = APA*2560+2560;
	}

  gStyle->SetOptStat(0);
  //gStyle->SetPalette(kRainBow);

  int  bins = ch_f - ch_i;

  cout << "total number of channels to analyzer: "<< bins << endl;

  TPrincipal p( bins, "ND");
  Double_t data[bins];

//TH2F* h2d = new TH2F("h2d","ch vs ticks",1000, 0, 1000, 4, 0, 4);



const int nticks = hu->GetNbinsY();
cout<<"n ticks: "<<nticks<<endl;
const int nchs = hu->GetNbinsX();
cout<<"n channels: "<<nchs<<endl;

for (int i=0;i!=nticks;i++){
   for(int channel=ch_i; channel<ch_f; channel++){
      int bin = channel - ch_i;
      //M[i][channel] = hu->GetBinContent(channel+1,i+1);
      data[bin] = hu->GetBinContent(channel+1,i+1);
   }
   p.AddRow(data);
}

const TMatrixD *m=p.GetCovarianceMatrix();
cout<<"n rows = "<<m->GetNrows()<<endl;
TMatrixDSym sym;
sym.Use(m->GetNrows(),m->GetMatrixArray());

cout<<sym(0,0)<<" "<<sym(0,1)<<" "<<sym(1,0)<<" "<<sym(1,1)<<" "<<endl;
//exit(0);
 TH2F* h = new TH2F(Form("CovM_%d%s_",APA+1,plane)+tag,Form("correlation -- APA%d %s plane",APA+1,plane)+tag,bins, ch_i, ch_f, bins, ch_i, ch_f);
 
 for(int i=0; i<bins; i++){

    for(int j=0; j<bins; j++){
        double rho = sym(i,j) / std::sqrt(sym(i,i) * sym(j,j)); // cov(x,y) = rho * sigma_x * sigma_y,
            //cout<<"rho = "<<i<<" "<<j<<" "<<sym(i,j)<<" / "<<std::sqrt(sym(i,i) * sym(j,j))<<" =      "<<rho<<endl;
                if (i != j) {
                  h->SetBinContent(i+1,j+1,rho);
                  h->SetBinContent(j+1,i+1,rho);
        if(rho ==1)  cout<<"rho = "<<i<<" "<<j<<" "<<sym(i,j)<<" / "<<std::sqrt(sym(i,i) * sym(j,j))<<"  = "<<rho<<endl;
                }

    }
  }

  return h;
}

void Init(){
    for(int i=0;i<4;i++)
    	for(int j=0;j<3;j++){
        channel_start[i][j] = (APA_number[i]-1)*2560+channels[j];
        channel_end[i][j] = (APA_number[i]-1)*2560+channels[j+1];
        cout<<"channel range: "<<channel_start[i][j]<<" "<< channel_end[i][j]<<endl;
        cout<<": "<<i<<"  "<<j<<"  "<<APA_number[i]<<" "<<channels[j]<<"   "<<channels[j+1]<<endl;
        h_orig_sep[i][j] = new TH2F(Form("raw_wf_%d_%s",APA_number[i],plane_name[j]),"",channels[j+1]-channels[j],channel_start[i][j]-0.5,channel_end[i][j]-0.5,nticks, 0, nticks);
        h_orig_ANf_sep[i][j] = new TH2F(Form("raw_wf_ANF_%d_%s",APA_number[i],plane_name[j]),"",channels[j+1]-channels[j],channel_start[i][j]-0.5,channel_end[i][j]-0.5,nticks, 0, nticks);
    }
	//exit(0);
}

void Draw_RMS(){
        h_rms = get_rms_hist(h_orig);
        h_rms_ANf = get_rms_hist(h_orig_ANf);
        h_rms->Draw();
        h_rms->GetYaxis()->SetRangeUser(0, 50);
        h_rms->GetXaxis()->SetTitle("channel number");
        h_rms_ANf->SetLineColor(2);
        h_rms_ANf->Draw("same");
        TLegend *leg = new TLegend(0.65, 0.75, 0.87, 0.87);
        leg->SetBorderSize(0);
        leg->AddEntry(h_rms, " Before Noise Filter", "l");
        leg->AddEntry(h_rms_ANf, "After Noise Filter", "l");
        leg->Draw();
    	c0->SaveAs(Form("/nashome/x/xning/Pictures/debug/rms_%d_%d.png", runNo, evtNo));
		h_rms->Delete();
}

void Draw_wf(){
        c0->cd();
        h_orig->Draw("colz");
        h_orig->GetZaxis()->SetRangeUser(-150, 150);
        h_orig->GetXaxis()->SetTitle("channel number");
        h_orig->GetYaxis()->SetTitle("nticks");
        h_orig->SetTitle("Raw waveform at " + time_string);
        c0->SaveAs(Form("/nashome/x/xning/Pictures/debug/raw_%d_%d.png", runNo, evtNo));

        h_orig_ANf->Draw("colz");
        h_orig_ANf->GetZaxis()->SetRangeUser(-150, 150);
        h_orig_ANf->GetXaxis()->SetTitle("channel number");
        h_orig_ANf->GetYaxis()->SetTitle("nticks");
        h_orig_ANf->SetTitle("After Nf at " + time_string);
        c0->SaveAs(Form("/nashome/x/xning/Pictures/debug/raw_ANf_%d_%d.png", runNo, evtNo));
}

void Draw_wf_sep(){
       for(int i=0;i<4;i++)
        for(int j=0;j<2;j++){
            h_orig_sep[i][j]->Draw("colz");
            h_orig_sep[i][j]->GetZaxis()->SetRangeUser(-150, 150);
            h_orig_sep[i][j]->GetXaxis()->SetTitle("channel number");
            h_orig_sep[i][j]->GetYaxis()->SetTitle("nticks");
            h_orig_sep[i][j]->SetTitle(Form("Raw waveform APA%d_%s_%d_%d ",APA_number[i],plane_name[j], runNo, evtNo));
            c0->SaveAs(Form("/nashome/x/xning/Pictures/debug/raw_apa%d_%s_%d_%d.png",APA_number[i],plane_name[j], runNo, evtNo));

            h_orig_ANf_sep[i][j]->Draw("colz");
            h_orig_ANf_sep[i][j]->GetZaxis()->SetRangeUser(-150, 150);
            h_orig_ANf_sep[i][j]->GetXaxis()->SetTitle("channel number");
            h_orig_ANf_sep[i][j]->GetYaxis()->SetTitle("nticks");
            h_orig_ANf_sep[i][j]->SetTitle(Form("Raw waveform ANf APA%d_%s_%d_%d ",APA_number[i],plane_name[j], runNo, evtNo));
            c0->SaveAs(Form("/nashome/x/xning/Pictures/debug/raw_ANf_apa%d_%s_%d_%d.png",APA_number[i],plane_name[j], runNo, evtNo));
        } 
		for(int i=0;i<4;i++)
			for(int j=2;j<3;j++){
			h_orig_sep[i][j]->Draw("colz");
            h_orig_sep[i][j]->GetZaxis()->SetRangeUser(0, 400);
            h_orig_sep[i][j]->GetXaxis()->SetTitle("channel number");
            h_orig_sep[i][j]->GetYaxis()->SetTitle("nticks");
            h_orig_sep[i][j]->SetTitle(Form("Raw waveform APA%d_%s_%d_%d ",APA_number[i],plane_name[j], runNo, evtNo));
            c0->SaveAs(Form("/nashome/x/xning/Pictures/debug/raw_apa%d_%s_%d_%d.png",APA_number[i],plane_name[j], runNo, evtNo));

            h_orig_ANf_sep[i][j]->Draw("colz");
            h_orig_ANf_sep[i][j]->GetZaxis()->SetRangeUser(0, 400);
            h_orig_ANf_sep[i][j]->GetXaxis()->SetTitle("channel number");
            h_orig_ANf_sep[i][j]->GetYaxis()->SetTitle("nticks");
            h_orig_ANf_sep[i][j]->SetTitle(Form("Raw waveform ANf APA%d_%s_%d_%d ",APA_number[i],plane_name[j], runNo, evtNo));
            c0->SaveAs(Form("/nashome/x/xning/Pictures/debug/raw_ANf_apa%d_%s_%d_%d.png",APA_number[i],plane_name[j], runNo, evtNo));
		}
}

void Draw_baseline(){
    h_baseline->Draw();
    c0->SaveAs(Form("/nashome/x/xning/Pictures/debug/baseline_%d_%d.png", runNo, evtNo));
}

void Draw_fft(){
            TLegend *leg = new TLegend(0.65, 0.75, 0.87, 0.87);
            leg->SetBorderSize(0);
    for(int i=0;i<4;i++)
        for(int j=0;j<3;j++){
            hfft[i][j] = CombineFFT(h_orig, APA_number[i], plane_name[j]);
            hfft_ANf[i][j] = CombineFFT(h_orig_ANf, APA_number[i], plane_name[j]);
            hfft[i][j]->Draw("hist");
            hfft_ANf[i][j]->SetLineColor(2);
            hfft_ANf[i][j]->Draw("same hist");
            if(i==0 && j==0){
                leg->AddEntry(hfft[0][0], " Before Noise Filter", "l");
                leg->AddEntry(hfft_ANf[0][0], "After Noise Filter", "l");
            }
            leg->Draw();
            c0->SaveAs(Form("/nashome/x/xning/Pictures/debug/fft%d%s_%d_%d.png",APA_number[i],plane_name[j], runNo, evtNo));
    }
}

void Draw_cov(){
        for(int i=0;i<4;i++)
        for(int j=0;j<3;j++){
            h_cov[i][j] = CovMatrix(h_orig,plane_name[j] , APA_number[i]);
            h_cov[i][j]->GetXaxis()->SetTitle("channel number");
            h_cov[i][j]->GetYaxis()->SetTitle("channel number");
            h_cov[i][j]->GetZaxis()->SetRangeUser(-0.3, 0.5);
            h_cov[i][j]->Draw("colz");
            c0->SaveAs(Form("/nashome/x/xning/Pictures/debug/covMatx_%d%s_%d_%d.png",APA_number[i],plane_name[j], runNo, evtNo));

            h_cov_ANf[i][j] = CovMatrix(h_orig_ANf,plane_name[j] , APA_number[i],"ANf");
            h_cov_ANf[i][j]->GetXaxis()->SetTitle("channel number");
            h_cov_ANf[i][j]->GetYaxis()->SetTitle("channel number");
            h_cov_ANf[i][j]->GetZaxis()->SetRangeUser(-0.3, 0.5);
            h_cov_ANf[i][j]->Draw("colz");
            c0->SaveAs(Form("/nashome/x/xning/Pictures/debug/cov_ANfMatx_%d%s_%d_%d.png",APA_number[i],plane_name[j], runNo, evtNo));

        }
        
}