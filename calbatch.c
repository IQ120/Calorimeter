#include <fstream>
#include <iostream>


#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TString.h"

#include "TStyle.h"




void calbatch(){

    //defining constants, files, functions
    const int nf = 23;
    std::unique_ptr<TFile> myFile (TFile::Open("calout.root", "RECREATE"));
    TString fnames[nf] = {"out_MeV_100.root", "out_MeV_200.root", "out_MeV_300.root", "out_MeV_400.root", "out_MeV_500.root", "out_MeV_600.root", "out_MeV_700.root", "out_MeV_800.root", "out_MeV_900.root", "out_MeV_1000.root", "out_GeV_2.root", "out_GeV_3.root", "out_GeV_4.root", "out_GeV_5.root", "out_GeV_6.root", "out_GeV_7.root", "out_GeV_8.root", "out_GeV_9.root", "out_GeV_10.root", "out_GeV_20.root", "out_GeV_30.root", "out_GeV_40.root", "out_GeV_50.root"};
    TFile* readfiles[nf];
    TH1D* EAbs_hists[nf];
    TH1D* hist_copy[nf];
    TH1D* Egap_hists[nf];
    TF1* fit_fce[nf];
    TString dopt = "";
    double const conversion_constant = 1000;

    TCanvas* canvas = new TCanvas ("histogram", "histogram", 1200, 800);
    canvas -> cd();
    gStyle -> SetOptStat(0);
    double xmin, xmax;

    TString fit_formula = "";

    //for cycle, drawing histogram


/*
    for(int i = 0; i < nf; i++){
        readfiles[i] = new TFile(fnames[i]);
        if (!readfiles[i] || readfiles[i]->IsZombie()){
            std::cerr<<"error opening file"<<endl; //error control
            return;
        }
        EAbs_hists[i] = (TH1D*) readfiles[i] -> Get("Eabs");

        //EAbs_hists[i]->SetLineColor(i+2);
        Egap_hists[i] = (TH1D*) readfiles[i] -> Get("Egap");

        hist_copy[i] = (TH1D*) EAbs_hists[i]->Clone();
        hist_copy[i]->Draw ("hist" + dopt);
        xmin = hist_copy[i]->GetXaxis()->GetXmin();
        xmax = hist_copy[i]->GetXaxis()->GetXmax();
        hist_copy[i]->GetXaxis()->SetLimits(xmin/1000, xmax/1000);
        dopt = "same";
    }
    */
    for(int i = 0; i < nf; i++){
        readfiles[i] = new TFile(fnames[i]);
        if (!readfiles[i] || readfiles[i]->IsZombie()){
            std::cerr<<"error opening file"<<endl; //error control
            return;
        }
        Egap_hists[i] = (TH1D*) readfiles[i] -> Get("Egap");

        hist_copy[i] = (TH1D*) Egap_hists[i]->Clone();
        hist_copy[i]->Draw ("hist" + dopt);
        xmin = hist_copy[i]->GetXaxis()->GetXmin();
        xmax = hist_copy[i]->GetXaxis()->GetXmax();
        hist_copy[i]->GetXaxis()->SetLimits(xmin/1000, xmax/1000);
        dopt = "same";
    }
    //hist_copy[8]->SetLineColor(59);
    hist_copy[0]->SetTitle("Energie deponovana v absorberu");
    hist_copy[0]->GetXaxis()->SetRangeUser (0,4);
    hist_copy[0]->GetXaxis()->SetTitle("E_{dep} (GeV)");
    hist_copy[0]->GetYaxis()->SetTitle("counts");

//legend string, for cycle legend generation
    TString legend[nf] = {"e^{-} energy 10 GeV","e^{-} energy 20 GeV","e^{-} energy 30 GeV","e^{-} energy 40 GeV","e^{-} energy 50 GeV","e^{-} energy 60 GeV","e^{-} energy 70 GeV","e^{-} energy 80 GeV","e^{-} energy 90 GeV","e^{-} energy 100 GeV"};

    auto legend1 = new TLegend(0.7,0.7,0.9,0.9);
    for (int i = 0; i < nf; i++ ){
        legend1->AddEntry(hist_copy[i], legend[i], "l");
    }
    legend1 -> Draw("same");

//saving histogram canvas

    canvas->SaveAs("histogram.pdf");





    //drawing absorbed energy graph + parameter field definitions

    //TCanvas* Energy_canvas = new TCanvas ("Energy_canvas", "Energy_canvas", 1200, 800);
    //Energy_canvas -> cd();
    double Ee[nf], Einf[nf], error_Einf[nf], Sigma[nf], error_Sigma[nf];

    TString line_fun = "[0] * x + [1]";
    TString fun_title = "linear function";

    auto linear_function = new TF1 (fun_title, line_fun, 0, 100);
    linear_function->SetLineWidth(1);

    for (int i = 0; i< 10; i++){
        Ee[i] = (i+1)*0.1;
    }
    for (int i = 0; i< 8; i++){
        Ee[i+10] = (i+2);
    }
    for (int i = 1; i< 6; i++){
        Ee[i+17] = i*10;
    }
    /*for (int i = 0; i < nf; i++){
        Einf[i] = EAbs_hists[i] -> GetMean()/1000;
        error_Einf[i] = EAbs_hists[i] -> GetMeanError()/1000;
        Sigma[i] = EAbs_hists[i] -> GetStdDev()/1000;
        error_Sigma[i] = EAbs_hists[i] -> GetStdDevError()/1000;
    }*/

    /*
    auto Energy_graph = new TGraphErrors(nf, Ee, Einf, 0, error_Einf);
    Energy_graph->SetTitle("Zavislost energie deponovane v absorberu na energii elektronu");
    Energy_graph->GetXaxis()->SetTitle("E e^{-} (GeV)");
    Energy_graph->GetYaxis()->SetTitle("E_{dep} (GeV)");

    Energy_graph -> Draw("A*");
    Energy_graph->Fit(linear_function);



    Energy_canvas -> SaveAs("Energy_graph.pdf");
    */

    //getting gap energy histogram, doing stuff with it
    double gapstddev[nf], gapstddeverror[nf], gapEabs[nf], gapEabs_error[nf], gap_norm_y[nf], gap_norm_y_error[nf];

    for (int i = 0; i<nf; i++){
        gapstddev[i] = Egap_hists[i]->GetStdDev()/1000;
        gapstddeverror[i] = Egap_hists[i]->GetStdDevError()/1000;
        gapEabs[i] = Egap_hists[i] ->GetMean()/1000;
        gapEabs_error[i] = Egap_hists[i]->GetMeanError()/1000;
    }

    //gap stddev graph on electron energy

    TString gapgraph_title = "Gap canvas";
    TCanvas* Gap_canvas = new TCanvas (gapgraph_title, gapgraph_title, 1200, 800);
    Gap_canvas->cd();
    auto gap_graph = new TGraphErrors(nf, Ee, gapstddev, 0, gapstddeverror);
    gap_graph->SetTitle("zavislost smerodatne odchylky energie absorbovane v prostoru na energii elektronu");
    gap_graph->GetXaxis()->SetTitle("E e^{-} (GeV)");
    gap_graph->GetYaxis()->SetTitle("#sigma_{gap} (GeV)");
    gap_graph->GetYaxis()->SetRangeUser(0, 0.1);
    gap_graph->Draw("A*");

    auto gapstddev_function = new TF1 ("gapstddev function", "sqrt(x*[0])", 0, 50);
    gapstddev_function -> SetParameter(0, 0.00015);

    gapstddev_function->SetLineWidth(1);
    gapstddev_function->Draw("same");
    gap_graph->Fit(gapstddev_function);


    auto legend2 = new TLegend(0.1,0.8,0.3,0.9);
    legend2->AddEntry("gap_graph", "data s chybou", "ep");
    legend2->AddEntry("gapstddev_function", "funkce #sqrt{#frac{x}{a}}", "l");
    legend2 -> Draw("same");
    legend2->SetTextSize(0.02);
    Gap_canvas->SaveAs("gap_graph.pdf");

    //Energy dep in scintillator graph

    TCanvas* Gapenergy_canvas = new TCanvas ("gap energy", "gap energy", 1200, 800);
    Gapenergy_canvas->cd();
    auto gapenergy_graph = new TGraphErrors(nf, Ee, gapEabs, 0, gapEabs_error);
    gapenergy_graph->Draw("A*");
    gapenergy_graph->SetTitle("zavislost energie absorbovane scintilatorem na energii elektronu");
    gapenergy_graph->GetXaxis()->SetTitle("E e^{-} (GeV)");
    gapenergy_graph->GetYaxis()->SetTitle("E_{gap} (GeV)");

    auto gapenergy_func = new TF1 ("gapenergy function", "x*[0]+[1]", 0, 50);
    gapenergy_func->SetLineWidth(1);
    gapenergy_graph->Fit(gapenergy_func);

    auto legend3 = new TLegend(0.1,0.8,0.3,0.9);
    legend3->AddEntry("gapenergy_graph", "data s chybou", "ep");
    legend3->AddEntry("gapenergy_func", "funkce ax + b", "l");
    legend3 -> Draw("same");
    legend3->SetTextSize(0.02);
    Gap_canvas->SaveAs("gap_graph.pdf");

    Gapenergy_canvas->SaveAs("gapenergy_graph.pdf");

    //normallized gap stddev graph

    for (int i = 0; i < nf; i++){
        gap_norm_y[i] = gapstddev[i]/gapEabs[i];
        gap_norm_y_error[i] = sqrt((gapstddeverror[i] * gapstddeverror[i])/(gapEabs[i]*gapEabs[i])+(gapEabs_error[i] * gapEabs_error[i])*(gapstddev[i]/(gapEabs[i]*gapEabs[i]) * (gapstddev[i]/(gapEabs[i]*gapEabs[i]))));
    }
    TCanvas* gap_norm_canvas = new TCanvas ("gap normallized_canvas", "gap normallized_canvas", 1200, 800);

    auto gap_norm_graph = new TGraphErrors(nf, Ee, gap_norm_y, 0, gap_norm_y_error);
    gap_norm_graph->SetTitle("Zavislost normovane smerodatne odchylky na energii elektronu");
    gap_norm_graph->GetXaxis()->SetTitle("E e^{-} (GeV)");
    gap_norm_graph->GetYaxis()->SetTitle("#sigma_{gap}/Egap_{dep}");
    gap_norm_graph->GetYaxis()->SetRangeUser(0, 0.8);
    gap_norm_graph -> Draw("A*");




    TString linfrac_fun = "1/sqrt(x*[0])";
    TString linfracfun_title = "linear fraction function";

    auto linfrac_function = new TF1 (linfracfun_title, linfrac_fun, 0, 50);
    linfrac_function -> SetParameter(0, 5);
    linfrac_function->SetLineWidth(1);
    linfrac_function->Draw("same");
    gap_norm_graph->Fit(linfrac_function);
/*
    auto theory_fun2 = new TF1 ("theory_fun2", "1/x", 0, 50);
    theory_fun2->SetLineWidth(1);
    theory_fun2->SetLineColor(1);
    theory_fun2->Draw("same");
*/

    auto legend4 = new TLegend(0.7,0.8,0.9,0.9);
    legend4->AddEntry("gap_norm_graph", "data s chybou", "ep");
    legend4->AddEntry("linfrac_function", "funkce #sqrt{#frac{1}{ax}}", "l");
    legend4 -> Draw("same");
    legend4->SetTextSize(0.02);

    gap_norm_canvas -> SaveAs ("Gap_Normallized_graph.pdf");



    /*
    //sigma graph drawing

    TCanvas* Sigma_canvas = new TCanvas ("sigma_canvas", "sigma_canvas", 1200, 800);
    Sigma_canvas -> cd();

    auto Sigma_graph = new TGraphErrors(nf, Ee, Sigma, 0, error_Sigma);
    Sigma_graph-> SetTitle("zavislost smerodatne odchylky na energii elektronu");
    Sigma_graph->GetXaxis()->SetTitle("E e^{-} (GeV)");
    Sigma_graph->GetYaxis()->SetTitle("#sigma (GeV)");
    Sigma_graph->GetYaxis()->SetRangeUser(0, 3.5);

    Sigma_graph->Fit(linear_function);

    Sigma_graph -> Draw("A*");
    linear_function->Draw("same");
    Sigma_canvas -> SaveAs ("Sigma_graph.pdf");

    //normallized graph drawing

    TCanvas* norm_canvas = new TCanvas ("normallized_canvas", "normallized_canvas", 1200, 800);

    norm_canvas -> cd();

    double norm_y[nf], norm_y_error[nf];

    for (int i = 0; i < nf; i++){
        norm_y[i] = Sigma[i]/Einf[i];
        norm_y_error[i] = sqrt((error_Sigma[i] * error_Sigma[i])/(Einf[i]*Einf[i])+(error_Einf[i] * error_Einf[i])*(Sigma[i]/(Einf[i]*Einf[i]) * (Sigma[i]/(Einf[i]*Einf[i]))));
    }

    auto norm_graph = new TGraphErrors(nf, Ee, norm_y, 0, norm_y_error);
    norm_graph->SetTitle("Depencency of normallized StdDev on electron energy");
    norm_graph->GetXaxis()->SetTitle("E e^{-} (GeV)");
    norm_graph->GetYaxis()->SetTitle("#sigma/E_{dep}");
    norm_graph -> Draw("A*");
    norm_canvas -> SaveAs ("Normallized_graph.pdf");
    */

    }
