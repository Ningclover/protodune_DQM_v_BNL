#include "all.h"

using namespace art;
using namespace std;
using namespace std::chrono;

void convert_artroot_DQM_2(std::string const &filename = "../ProcessData/np04hd_raw_run027980_0663_dataflow3_datawriter_0_20240711T212207_reco_stage1.root")
{
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kLightTemperature);
    Init();
    InputTag rawdigits_tag{"tpcrawdecoder:daq"};
    InputTag AftNoise_tag{"wclsdatahdfilter:raw"};


    // Create a vector of length 1, containing the given filename.
    vector<string> filenames(1, filename);
    int n = 0;

    for (gallery::Event ev(filenames); !ev.atEnd(); ev.next())
    {

        auto aux = ev.eventAuxiliary();
        std::cout << "processing event: " << aux.run() << "-" << aux.event() << std::endl;

         if (n>5) break;

        // std::cout <<"time  "<<aux.time().value()<<"  "<<aux.time().timeLow()<<"  "<<aux.time().timeHigh()<<std::endl;
        auto t1 = TTimeStamp(aux.time().timeHigh(), aux.time().timeLow());
        time_string = t1.AsString();
        std::cout << t1.AsString() << std::endl;

        runNo = aux.run();
        evtNo = aux.event();
        evtNo_list.push_back(evtNo);
        auto const &rawdigits = *ev.getValidHandle<vector<raw::RawDigit>>(rawdigits_tag);
        auto const &AftNoise = *ev.getValidHandle<vector<raw::RawDigit>>(AftNoise_tag);

        if (rawdigits.empty())
        {
            std::cout << "WARNING: no RawDigit found." << std::endl;
            return;
        }

        if (AftNoise.empty())
        {
            std::cout << "WARNING: no wclsdatahdfilter:raw found." << std::endl;
            return;
        }
        double maxchannel = 0;

        for (auto rd : rawdigits)
        {
            if (rd.Channel() == 0)
                std::cout << rd.Samples() << std::endl;
            //       std::cout << "channel: " << rd.Channel()
            //                 << " NADC: " << rd.NADC()
            //                 << " Samples: " << rd.Samples()
            //                 << " Pedestal: " << rd.GetPedestal()
            //                 << " Sigma: " << rd.GetSigma()
            //                 << std::endl;
            if (maxchannel < rd.Channel())
                maxchannel = rd.Channel();

            int channel = rd.Channel();

            h_baseline->SetBinContent(channel + 1, rd.GetPedestal());

            // h_rms->SetBinContent(channel+1, rd.GetSigma());

            int nSamples = rd.Samples();
            for (int j = 0; j < nSamples; j++)
            {
                h_orig->SetBinContent(channel + 1, j + 1, rd.ADC(j) - rd.GetPedestal());
                for(int k=0;k<4;k++)
                    for(int l=0;l<3;l++){
                        if(channel>=channel_start[k][l]&&channel<channel_end[k][l]){
                            h_orig_sep[k][l]->SetBinContent(channel + 1-channel_start[k][l],j + 1, rd.ADC(j) - rd.GetPedestal());
                            h_orig_sep[k][l]->SetTitle(Form("raw_wf_run_%d_evt_%d_APA%d_%s",runNo,evtNo,APA_number[k],plane_name[l]));
                        }
                }


            }

        } // end of rawdigits
        //     cout<<"test content"<<h_orig_sep[1][0]->Integral()<<endl;
           //cout<<"channel_start = "<<channel_start[0][1]<<"  "<<channel_end[0][1]<<" bincontent "<<h_orig_sep[0][1]->GetBinContent(601,5849)<<endl;
        //    cout<<h_orig_sep[0][1]->Integral()<<endl;
        //    h_orig_sep[0][1]->Draw("colz");
        //    h_orig_sep[0][1]->GetZaxis()->SetRangeUser(-150,150);
        //    c0->SaveAs("/nashome/x/xning/Pictures/test.png");
        // exit(0);

        for (auto rd : AftNoise)
        {
            //    std::cout << "channel: " << rd.Channel()
            //              << " NADC: " << rd.NADC()
            //              << " Samples: " << rd.Samples()
            //              << " Pedestal: " << rd.GetPedestal()
            //              << " Sigma: " << rd.GetSigma()
            //              << std::endl;
            if (maxchannel < rd.Channel())
                maxchannel = rd.Channel();

            int channel = rd.Channel();

            int nSamples = rd.Samples();
            for (int j = 0; j < nSamples; j++)
            {
                h_orig_ANf->SetBinContent(channel + 1, j + 1, rd.ADC(j) - rd.GetPedestal());
                    for(int k=0;k<4;k++)
                        for(int l=0;l<3;l++){
                        if(channel>channel_start[k][l]&&channel<channel_end[k][l]){
                            h_orig_ANf_sep[k][l]->SetBinContent(channel + 1-channel_start[k][l],j + 1, rd.ADC(j) - rd.GetPedestal());
                            h_orig_ANf_sep[k][l]->SetTitle(Form("raw_wf_ANf_run_%d_evt_%d_APA%d_%s",runNo,evtNo,APA_number[k],plane_name[l]));
                        }
                    }
            }
        } // end of AftNoise digits
        // cout<<"maxchannel:  "<<maxchannel<<endl;
        Draw_wf();
        Draw_wf_sep();
        Draw_RMS();
        Draw_baseline();
        Draw_fft();
        Draw_cov();
        n++;
    }
    // c1->SaveAs("/nashome/x/xning/Pictures/outfile.pdf]");
    // ofstream fout;
    // TString name = Form("/nashome/x/xning/runinfo/runNo%d.dat",runNo);
    // fout.open(name.Data(),std::ios_base::app);
    // for(auto event : evtNo_list){
    //     fout<<runNo<<"\t"<<event<<endl;
    // }
    // fout.close();

}
