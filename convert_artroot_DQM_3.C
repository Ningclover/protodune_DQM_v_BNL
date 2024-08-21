#include "all.h"

using namespace art;
using namespace std;
using namespace std::chrono;

void convert_artroot_DQM_3(std::string const &filename = "../ProcessData/np04hd_raw_run027344_0001_dataflow3_datawriter_0_20240621T145342_reco_stage1_reco_stage2_20240731T205640_keepup.root", TString outpath="/exp/dune/data/users/xning/Pictures/DQM/")
{
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kLightTemperature);
    Init();
    // InputTag rawdigits_tag{"tpcrawdecoder:daq"};
    // InputTag AftNoise_tag{"wclsdatahdfilter:raw"};
    InputTag decon_tag{"wclsdatahd:gauss"};
    InputTag gausshit_tag{"gaushit"};


    // Create a vector of length 1, containing the given filename.
    vector<string> filenames(1, filename);
    int n = 0;

    for (gallery::Event ev(filenames); !ev.atEnd(); ev.next())
    {

        auto aux = ev.eventAuxiliary();
        std::cout << "processing event: " << aux.run() << "-" << aux.event() << std::endl;

 
         if (n>15) break;

        // std::cout <<"time  "<<aux.time().value()<<"  "<<aux.time().timeLow()<<"  "<<aux.time().timeHigh()<<std::endl;
        auto t1 = TTimeStamp(aux.time().timeHigh(), aux.time().timeLow());
        time_string = t1.AsString();
        std::cout << t1.AsString() << std::endl;

        runNo = aux.run();
        evtNo = aux.event();
        TString path = outpath+Form("/run_%d/evt_%d/",runNo,evtNo);
        gSystem->mkdir(path.Data(),true);
        outname=path;
        evtNo_list.push_back(evtNo);
        // auto const &rawdigits = *ev.getValidHandle<vector<raw::RawDigit>>(rawdigits_tag);
        // auto const &AftNoise = *ev.getValidHandle<vector<raw::RawDigit>>(AftNoise_tag);
        auto const &deco_wire = *ev.getValidHandle<vector<recob::Wire>>(decon_tag);
        auto const &gaus_hit = *ev.getValidHandle<vector<recob::Hit>>(gausshit_tag);

        if (deco_wire.empty())
        {
            std::cout << "WARNING: no RawDigit found." << std::endl;
            return;
        }

        if (gaus_hit.empty())
        {
            std::cout << "WARNING: no wclsdatahdfilter:raw found." << std::endl;
            return;
        }

        for (auto wire : deco_wire)
        {
            int channel = wire.Channel();
            int nSamples = wire.NSignal();
            std::vector<float> signal = wire.Signal();
            // h_rms->SetBinContent(channel+1, rd.GetSigma());

            for (int j = 0; j < nSamples; j++)
            {
                for(int k=0;k<4;k++)
                    for(int l=0;l<3;l++){
                        if(channel>=channel_start[k][l]&&channel<channel_end[k][l]){
                            h_decon_sep[k][l]->SetBinContent(channel + 1-channel_start[k][l],j + 1, signal[j]);
                            h_decon_sep[k][l]->SetTitle(Form("decon_run_%d_evt_%d_APA%d_%s",runNo,evtNo,APA_number[k],plane_name[l]));
                        }
                }

            }

        } 

           for (auto hit : gaus_hit)
        {
            int channel = hit.Channel();
            float sumADC = hit.SummedADC();
            // float Intgl = hit.Integral();
            h_nHit->Fill(channel);
            h_nHit_ADC->Fill(channel,sumADC);
            // cout<<"channel "<<channel<<" sumADC = "<<sumADC<<" Integral =  "<<Intgl<<endl;
        } 
        Draw_decon_sep();
        Draw_Hit();
        n++;
    }
    // c1->SaveAs("/nashome/x/xning/Pictures/outfile.pdf]");
    //  ofstream fout;
    //  //TString name = Form("/nashome/x/xning/runinfo/runNo%d.dat",runNo);
    //  TString name = Form("/nashome/x/xning/runinfo/new/runNo%d.dat",runNo);
    //  //fout.open(name.Data(),std::ios_base::app);
    //  fout.open(name.Data());
    //  for(auto event : evtNo_list){
    //      fout<<runNo<<"\t"<<event<<endl;
    //  }
    //  fout.close();

}
