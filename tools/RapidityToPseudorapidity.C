#include <iostream>

#include "TCanvas.h"
#include "TDatabasePDG.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TROOT.h"

/*
 From the particle species PDG code, plot a 2D histogram of Pt vs. Rapidity
 that shows the correct phase space when those two variables give the correct eta range that we are looking for.
 => In our particular case: |eta| < 1.5
 */
void RapidityToPseudorapidity(Int_t part_pdg) {

    /* Transverse Momentum (GeV/c) */

    Float_t min_pt;
    Float_t max_pt;
    Int_t nbins_pt;
    if (part_pdg == 211) {
        min_pt = 0.05;
        max_pt = 0.15;
        nbins_pt = 100;
    }
    if (part_pdg == 321) {
        min_pt = 0.1;
        max_pt = 0.25;
        nbins_pt = 100;
    }
    if (part_pdg == 2212) {
        min_pt = 0.15;
        max_pt = 0.35;
        nbins_pt = 100;
    }
    Float_t step_pt = (max_pt - min_pt) / (Double_t)nbins_pt;

    /* Rapidity */

    Float_t min_y = -1.5;
    Float_t max_y = 1.5;
    Int_t nbins_y = 200;
    Float_t step_y = (max_y - min_y) / (Double_t)nbins_y;

    /* And the rest... */

    TDatabasePDG fPDG;
    Float_t part_mass = fPDG.GetParticle(part_pdg)->Mass();

    Float_t part_phi = 0.;  // independent variable

    Float_t part_px;
    Float_t part_py;
    Float_t part_mt;
    Float_t part_pz;
    Float_t part_p;
    Float_t part_energy;

    TLorentzVector part_tlv;

    TH2I* hist_pt_vs_y = new TH2I("Pt vs Y.", "Pt vs Y.", nbins_y, min_y, max_y, nbins_pt, min_pt, max_pt);

    for (Float_t part_pt = min_pt; part_pt <= max_pt; part_pt += step_pt) {

        std::cout << ">> Pt = " << part_pt << std::endl;

        for (Float_t part_y = min_y; part_y <= max_y; part_y += step_y) {

            std::cout << "   >> Rapidity = " << part_y << std::endl;

            part_px = part_pt * TMath::Cos(part_phi);
            part_py = part_pt * TMath::Sin(part_phi);
            part_mt = TMath::Sqrt(part_mass * part_mass + part_pt * part_pt);
            part_pz = part_mt * TMath::SinH(part_y);
            part_energy = part_mt * TMath::CosH(part_y);
            part_p = TMath::Sqrt(part_px * part_px + part_py * part_py + part_pz * part_pz);

            part_tlv.SetPxPyPzE(part_px, part_py, part_pz, part_energy);

            std::cout << "      >> Eta = " << part_tlv.Eta() << std::endl;

            if (TMath::Abs(part_tlv.Eta()) < 1.5) hist_pt_vs_y->Fill(part_y, part_pt);
        }
    }

    TCanvas* canvas_pt_vs_y = new TCanvas("canvas_pt_vs_y", "canvas_pt_vs_y", 1280, 1280);
    hist_pt_vs_y->Draw("COLZ");
    canvas_pt_vs_y->Draw();
}