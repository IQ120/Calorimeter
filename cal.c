#include <fstream>
#include <iostream>


#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TString.h"

#include "TStyle.h"




void cal(){
    std::unique_ptr<TFile> myFile (TFile::Open("out_20.root"));
    std::unique_ptr<TFile> myFile2 (TFile::Open("out_40.root"));
    std::unique_ptr<TFile> myFile3 (TFile::Open("out_60.root"));
    std::unique_ptr<TFile> myFile4 (TFile::Open("out_80.root"));
    std::unique_ptr<TFile> myFile5 (TFile::Open("out_100.root"));
    std::unique_ptr<TFile> myFile6 (TFile::Open("calout.root", "RECREATE"));

    TCanvas* canvas = new TCanvas ("canvas", "canvas", 1200, 800);

    myFile->Print();
    TH1D* hist = (TH1D*)myFile->Get("Eabs");
    gStyle->SetOptStat(0);

    hist ->SetLineColor(1);
    hist->SetTitle ("Calorimeter");
    hist->GetXaxis()->SetTitle("bin");
    hist->GetYaxis()->SetTitle("count");
    //hist->GetXaxis()->SetRangeUser(0,100);
    //hist->GetYaxis()->SetRangeUser(0,5000);
    hist->Draw("Hist");

    // fitting gaussian function
    TF1* fit_fce1 = new TF1("fit_fce1", "gaus(0)");
    hist->Fit(fit_fce1);
    fit_fce1->Draw("same");

    //mean and sigma of gauss
    double mean_20GeV= fit_fce1->GetParameter(1);
    double sigma_20Gev = fit_fce1->GetParameter(2);

    //2nd histogram

    TH1D* hist_a = (TH1D*)myFile2->Get("Eabs");

    hist_a ->SetLineColor(2);
    hist_a->Draw("hist same");

    TF1* fit_fce2 = new TF1("fit_fce2", "gaus(0)", 30e3, 40e3);
    fit_fce2->SetParameter(0,750);
    hist_a->Fit(fit_fce2, "R");
    fit_fce2->Draw("same");

    double mean_40GeV= fit_fce2->GetParameter(1);
    double sigma_40Gev = fit_fce2->GetParameter(2);


    //3rd histogram - 60 GeV

    TH1D* hist_b = (TH1D*)myFile3->Get("Eabs");

    hist_b->SetLineColor(3);
    hist_b->Draw("hist same");

    TF1* fit_fce3 = new TF1("fit_fce3", "gaus(0)", 40e3, 55e3);
    hist_b->Fit(fit_fce3, "R");
    fit_fce3->Draw("same");

    double mean_60GeV= fit_fce3->GetParameter(1);
    double sigma_60Gev = fit_fce3->GetParameter(2);

    //4th histogram

    TH1D* hist_c = (TH1D*)myFile4->Get("Eabs");

    hist_c->SetLineColor(4);
    hist_c->Draw("hist same");

    TF1* fit_fce4 = new TF1("fit_fce4", "gaus(0)", 55e3, 75e3);
    hist_c->Fit(fit_fce4, "R");
    fit_fce4->Draw("same");

    double mean_80GeV= fit_fce4->GetParameter(1);
    double sigma_80Gev = fit_fce4->GetParameter(2);

    //5th histogram

    TH1D* hist_d = (TH1D*)myFile5->Get("Eabs");

    hist_d->SetLineColor(5);
    hist_d->Draw("hist same");

    TF1* fit_fce5 = new TF1("fit_fce5", "gaus(0)", 70e3, 90e3);
    hist_d->Fit(fit_fce5, "R");
    fit_fce5->Draw("same");

    double mean_100GeV= fit_fce5->GetParameter(1);
    double sigma_100Gev = fit_fce5->GetParameter(2);

    //legend
    auto legend1 = new TLegend(0.7,0.7,0.9,0.9);
    legend1->AddEntry(hist, "e^- energy 20 GeV", "l");
    legend1->AddEntry(hist_a, "e^- energy 40 GeV", "l");
    legend1->AddEntry(hist_b, "e^- energy 60 GeV", "l");
    legend1->AddEntry(hist_c, "e^- energy 80 GeV", "l");
    legend1->AddEntry(hist_d, "e^- energy 100 GeV", "l");
    legend1->Draw("same");



    //save parameters to file


    TString outfilename = "calout.txt";
    ofstream *infile = 0;
    infile = new ofstream(outfilename.Data());

    (*infile)<<"fitmean_1: " << fit_fce1 -> GetParameter(1) << endl;
    (*infile)<<"fitmeanerror_1: " << fit_fce1-> GetParError(1) << endl;
    (*infile)<<"fitsigma_1: " << fit_fce1 -> GetParameter(2) << endl;
    (*infile)<<"fitsigmaerror_1: " << fit_fce1 -> GetParError(2) << endl;

    (*infile)<<"fitmean_2: " << fit_fce2 -> GetParameter(1) << endl;
    (*infile)<<"fitmeanerror_2: " << fit_fce2-> GetParError(1) << endl;
    (*infile)<<"fitsigma_2: " << fit_fce2 -> GetParameter(2) << endl;
    (*infile)<<"fitsigmaerror_2: " << fit_fce2 -> GetParError(2) << endl;

    (*infile)<<"fitmean_3: " << fit_fce3 -> GetParameter(1) << endl;
    (*infile)<<"fitmeanerror_3: " << fit_fce3-> GetParError(1) << endl;
    (*infile)<<"fitsigma_3: " << fit_fce3 -> GetParameter(2) << endl;
    (*infile)<<"fitsigmaerror_3: " << fit_fce3 -> GetParError(2) << endl;

    (*infile)<<"fitmean_4: " << fit_fce4 -> GetParameter(1) << endl;
    (*infile)<<"fitmeanerror_4: " << fit_fce4-> GetParError(1) << endl;
    (*infile)<<"fitsigma_4: " << fit_fce4 -> GetParameter(2) << endl;
    (*infile)<<"fitsigmaerror_4: " << fit_fce4 -> GetParError(2) << endl;

    (*infile)<<"fitmean_5: " << fit_fce5 -> GetParameter(1) << endl;
    (*infile)<<"fitmeanerror_5: " << fit_fce5-> GetParError(1) << endl;
    (*infile)<<"fitsigma_5: " << fit_fce5 -> GetParameter(2) << endl;
    (*infile)<<"fitsigmaerror_5: " << fit_fce5 -> GetParError(2) << endl;

    if (infile)
        infile->close();

    canvas->SaveAs("histogram.pdf");
    hist -> Write();
    hist_a-> Write();
    hist_b-> Write();
    hist_c-> Write();
    hist_d-> Write();
    canvas->Update();
    canvas->Close();

    double x[100], y[100], ey[100];
    int n = 5;
    for (int i = 1; i < 6; i++){
        x[i-1] = 4*n*i;
    }
    y[0] = fit_fce1 -> GetParameter(1)/1000;
    y[1] = fit_fce2 -> GetParameter(1)/1000;
    y[2] = fit_fce3 -> GetParameter(1)/1000;
    y[3] = fit_fce4 -> GetParameter(1)/1000;
    y[4] = fit_fce5 -> GetParameter(1)/1000;


    ey[0] = fit_fce1 -> GetParError(1)/1000;
    ey[1] = fit_fce2 -> GetParError(1)/1000;
    ey[2] = fit_fce3 -> GetParError(1)/1000;
    ey[3] = fit_fce4 -> GetParError(1)/1000;
    ey[4] = fit_fce5 -> GetParError(1)/1000;

    //defining linear function

    TString name = "linefun";
    TString formula = "x * [0] + [1]";
    TF1* linefun = new TF1(name, formula);
    linefun -> SetParameters(0.9, 0);

    // make a graph with energies (+errorbars) ,fit linear function
    TCanvas* energy_canvas = new TCanvas ("Energy_canvas", "Energy_canvas", 800, 600);
    energy_canvas -> cd(1);


    auto Energy_graph = new TGraphErrors (n, x, y, 0, ey);
    Energy_graph -> SetTitle("Energy absorbed ; e^- energy; absorbed energy");
    Energy_graph -> SetMinimum (0);
    Energy_graph -> SetMaximum (100);
    Energy_graph -> SetMarkerSize(1);
    Energy_graph -> SetMarkerStyle (20);
    Energy_graph -> Draw("AP");
    double xmin = 0, xmax = 100;

    Energy_graph -> Fit(linefun, "", "", xmin, xmax);
    linefun -> Draw("same");

    double yn[100];
    yn[0] = (fit_fce1 -> GetParameter(2) / fit_fce1 -> GetParameter (1));
    yn[1] = (fit_fce2 -> GetParameter(2) / fit_fce2 -> GetParameter (1));
    yn[2] = (fit_fce3 -> GetParameter(2) / fit_fce3 -> GetParameter (1));
    yn[3] = (fit_fce4 -> GetParameter(2) / fit_fce4 -> GetParameter (1));
    yn[4] = (fit_fce5 -> GetParameter(2) / fit_fce5 -> GetParameter (1));

    // normallized graph

    TCanvas* normallized_canvas = new TCanvas ("Energy_canvas", "Energy_canvas", 800, 600);
    normallized_canvas->cd();

    auto Normalized_graph = new TGraph (n, x, yn);
    Normalized_graph -> SetTitle ("Normallized energy graph ; e Energy; sigma/Eabs");

    Normalized_graph -> Draw ("A*");

    TF1* normal_fun = new TF1 ("normal_fun", "[0]*exp(x/[1)");
    Normalized_graph->Fit(normal_fun);
    Normalized_graph -> SetMinimum(0);
    Normalized_graph ->SetMaximum (0.1);

    }
