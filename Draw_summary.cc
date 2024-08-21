void Draw_summary(){
    TString path = "/exp/dune/data/users/xning/";

    ifstream fin("/exp/dune/data/users/xning/info/all.dat");
	int run,evt;
    double enc;
    vector<vector<double>> ppp(4,vector<double>(3));
    vector<double> runNo;
    vector<int> evtNo;
    vector<double> rms_mean;
    // vector<vector<vector<double>>> peak;
    vector<double> p1u;
    vector<double>p1v;
    vector<double>p1w;
    vector<double>p2u;
    vector<double>p2v;
    vector<double>p2w;
    vector<double>p3u;
    vector<double>p3v;
    vector<double>p3w;
    vector<double>p4u;
    vector<double>p4v;
    vector<double>p4w;
    int k=0;
    vector<double> count;
    int run_now=27977;
    int evt_count=0;
    while(fin>>run>>evt>>enc>>ppp[0][0]>>ppp[0][1]>>ppp[0][2]>>ppp[1][0]>>ppp[1][1]>>ppp[1][2]>>ppp[2][0]>>ppp[2][1]>>ppp[2][2]>>ppp[3][0]>>ppp[3][1]>>ppp[3][2]){
        evtNo.push_back(evt);
        if(k==0){
        rms_mean.push_back(enc);
        p1u.push_back(ppp[0][0]);
        p1v.push_back(ppp[0][1]);
        p1w.push_back(ppp[0][2]);
        p2u.push_back(ppp[1][0]);
        p2v.push_back(ppp[1][1]);
        p2w.push_back(ppp[1][2]);
        p3u.push_back(ppp[2][0]);
        p3v.push_back(ppp[2][1]);
        p3w.push_back(ppp[2][2]);
        p4u.push_back(ppp[3][0]);
        p4v.push_back(ppp[3][1]);
        p4w.push_back(ppp[3][2]);       
        runNo.push_back(run);
        run_now=run;
        cout<<"init run "<<run_now<<endl;

        }else{
                // cout<<" run "<<run<<" . "<<run_now<<endl;
            if(run==run_now){
                evt_count++;
                rms_mean.back() = rms_mean.back()+enc;
                p1u.back() = p1u.back()+ppp[0][0];
                p1v.back() = p1v.back()+ppp[0][1];
                p1w.back() = p1w.back()+ppp[0][2];
                p2u.back() = p2u.back()+ppp[1][0];
                p2v.back() = p2v.back()+ppp[1][1];
                p2w.back() = p2w.back()+ppp[1][2];
                p3u.back() = p3u.back()+ppp[2][0];
                p3v.back() = p3v.back()+ppp[2][1];
                p3w.back() = p3w.back()+ppp[2][2];
                p4u.back() = p4u.back()+ppp[3][0];
                p4v.back() = p4v.back()+ppp[3][1];
                p4w.back() = p4w.back()+ppp[3][2];
                // cout<<"evtcount = "<<evt_count<<" rms sum = "<<rms_mean.back()<<endl;
            }else{
                cout<<"# evt = "<<evt_count<<endl;
                rms_mean.back() = rms_mean.back()/evt_count;
                p1u.back() = p1u.back()/evt_count;
                p1v.back() = p1v.back()/evt_count;
                p1w.back() = p1w.back()/evt_count;
                p2u.back() = p2u.back()/evt_count;
                p2v.back() = p2v.back()/evt_count;
                p2w.back() = p2w.back()/evt_count;
                p3u.back() = p3u.back()/evt_count;
                p3v.back() = p3v.back()/evt_count;
                p3w.back() = p3w.back()/evt_count;
                p4u.back() = p4u.back()/evt_count;
                p4v.back() = p4v.back()/evt_count;
                p4w.back() = p4w.back()/evt_count;
                cout<<"evtcount = "<<evt_count<<" rms ave = "<<rms_mean.back()<<endl;
                cout<<p1u.back()<<" "
                    <<p1v.back()<<" "
                    <<p1w.back()<<endl;

                evt_count=0;
                run_now=run;
                rms_mean.push_back(enc);
                p1u.push_back(ppp[0][0]);
                p1v.push_back(ppp[0][1]);
                p1w.push_back(ppp[0][2]);
                p2u.push_back(ppp[1][0]);
                p2v.push_back(ppp[1][1]);
                p2w.push_back(ppp[1][2]);
                p3u.push_back(ppp[2][0]);
                p3v.push_back(ppp[2][1]);
                p3w.push_back(ppp[2][2]);
                p4u.push_back(ppp[3][0]);
                p4v.push_back(ppp[3][1]);
                p4w.push_back(ppp[3][2]);       
                runNo.push_back(run);

            }
        }





        count.push_back(k);
        k++;
    }
                    rms_mean.back() = rms_mean.back()/evt_count;
                p1u.back() = p1v.back()/evt_count;
                p1v.back() = p1v.back()/evt_count;
                p1w.back() = p1w.back()/evt_count;
                p2u.back() = p2v.back()/evt_count;
                p2v.back() = p2v.back()/evt_count;
                p2w.back() = p2w.back()/evt_count;
                p3u.back() = p3v.back()/evt_count;
                p3v.back() = p3v.back()/evt_count;
                p3w.back() = p3w.back()/evt_count;
                p4u.back() = p4v.back()/evt_count;
                p4v.back() = p4v.back()/evt_count;
                p4w.back() = p4w.back()/evt_count;

    TGraph *g1 = new TGraph(runNo.size(),&count[0],&rms_mean[0]);
    TGraph *gg1u = new TGraph(runNo.size(),&count[0],&p1u[0]);
    TGraph *gg1v = new TGraph(runNo.size(),&count[0],&p1v[0]);
    TGraph *gg1w = new TGraph(runNo.size(),&count[0],&p1w[0]);
    TGraph *gg2u = new TGraph(runNo.size(),&count[0],&p2u[0]);
    TGraph *gg2v = new TGraph(runNo.size(),&count[0],&p2v[0]);
    TGraph *gg2w = new TGraph(runNo.size(),&count[0],&p2w[0]);
    TGraph *gg3u = new TGraph(runNo.size(),&count[0],&p3u[0]);
    TGraph *gg3v = new TGraph(runNo.size(),&count[0],&p3v[0]);
    TGraph *gg3w = new TGraph(runNo.size(),&count[0],&p3w[0]);
    TGraph *gg4u = new TGraph(runNo.size(),&count[0],&p4u[0]);
    TGraph *gg4v = new TGraph(runNo.size(),&count[0],&p4v[0]);
    TGraph *gg4w = new TGraph(runNo.size(),&count[0],&p4w[0]);
    TCanvas *c1 = new TCanvas("c1","c1",1);

    g1->Draw("AL*");
    g1->GetYaxis()->SetRangeUser(500,600);
    g1->GetYaxis()->SetTitle("Average RMS[ENC]"); 
    c1->SaveAs(path+"/RMS.png");

    gg1u->Draw("AL*");
    gg1u->GetYaxis()->SetRangeUser(2000,3000);
    gg1u->GetYaxis()->SetTitle("peak of noise spectra");
    gg1u->SetLineStyle(1);
    gg1v->Draw("same l *");
    gg1v->SetLineStyle(2);
    gg1w->Draw("same l *");
    gg1w->SetLineStyle(6);

    gg2u->Draw("same l *");
    gg2u->SetLineStyle(1);
    gg2v->Draw("same l *");
    gg2v->SetLineStyle(2);
    gg2w->Draw("same l *");
    gg2w->SetLineStyle(6);
    gg2u->SetLineColor(2);
    gg2v->SetLineColor(2);
    gg2w->SetLineColor(2);

    c1->SaveAs(path+"/peak_1.png");

    gg3u->Draw("AL*");
    gg3u->GetYaxis()->SetRangeUser(2000,3000);
    gg3u->GetYaxis()->SetTitle("peak of noise spectra");
    gg3u->SetLineStyle(1);
    gg3v->Draw("same l *");
    gg3v->SetLineStyle(2);
    gg3w->Draw("same l *");
    gg3w->SetLineStyle(6);

    gg4u->Draw("same l *");
    gg4u->SetLineStyle(1);
    gg4v->Draw("same l *");
    gg4v->SetLineStyle(2);
    gg4w->Draw("same l *");
    gg4w->SetLineStyle(6);
    gg4u->SetLineColor(2);
    gg4v->SetLineColor(2);
    gg4w->SetLineColor(2);

    c1->SaveAs(path+"/peak_2.png");


}