#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"

void BetheFitter(Int_t PDGCode = 2212) {

    TFile *f = new TFile("/home/borquez/bachelors/ToyTPC/analysis/run00_ana_v1.root");

    TTree *t = (TTree *)f->Get("Events");

    /* Plot TPC signal */

    t->Draw("dEdx:P_rec>>hist1", Form("TMath::Abs(PDGcode) == %i", PDGCode), "goff");

    TH2F *hist1 = (TH2F *)gDirectory->Get("hist1");

    TCanvas *c1 = new TCanvas("c1", "c1", 800, 800);

    hist1->Draw("COLZ");

    /* Plot dEdx for different ranges of P_rec */

    Int_t nBins = 8;
    Float_t P_rec_min = 0.;
    Float_t P_rec_max = 800.;
    Float_t P_rec_delta = (P_rec_max - P_rec_min) / (Float_t)nBins;

    Float_t P_rec_low, P_rec_high;
    TString current_cut;

    TH1F *hist2;

    std::vector<Float_t> dEdx_mean;
    std::vector<Float_t> dEdx_sigma;

    TCanvas *c2 = new TCanvas("c2", "c2", 1600, 800);
    c2->Divide(4, 2);

    for (Int_t i = 0; i < nBins; i++) {

        c2->cd(i + 1);

        P_rec_low = P_rec_min + i * P_rec_delta;
        P_rec_high = P_rec_min + (i + 1) * P_rec_delta;

        current_cut = Form("TMath::Abs(PDGcode) == %i && P_rec > %f && P_rec < %f", PDGCode, P_rec_low, P_rec_high);

        t->Draw(Form("dEdx>>htemp_%i", i), current_cut, "goff");

        hist2 = (TH1F *)gDirectory->Get(Form("htemp_%i", i));
        hist2->Draw();

        // get mean and sigma of dEdx for each P_rec bin
        dEdx_mean.push_back(hist2->GetMean());
        dEdx_sigma.push_back(hist2->GetRMS());
    }

    /* Plot 1D dEdx distribution vs P_rec */

    TCanvas *c3 = new TCanvas("c3", "c3", 800, 800);

    TH1F *hist3 = new TH1F("hist3", "dEdx vs P_rec", nBins, P_rec_min, P_rec_max);

    for (Int_t i = 0; i < nBins; i++) {
        hist3->SetBinContent(i + 1, dEdx_mean[i]);
        hist3->SetBinError(i + 1, dEdx_sigma[i]);
    }

    hist3->Draw("E");

    /* Fit dEdx vs P_rec with Bethe Bloch formula */

    // TODO
}
