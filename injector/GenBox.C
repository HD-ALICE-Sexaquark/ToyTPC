#include <fstream>
#include <iostream>

#include "TDatabasePDG.h"
#include "TF1.h"
#include "TGenPhaseSpace.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TPDGCode.h"
#include "TRandom.h"

/*
 Generate some low-pT charged particles
 */
void GenBox(TString fOutputFilename = "../reconstruction/mc.csv") {

    // parameters
    const Int_t n_times = 50;  // inject `n_times` particles per species
    std::vector<Int_t> pdgVector = {211, 321, 2212};
    std::vector<Float_t> fPtMin = {0.05, 0.1, 0.15};   // (in GeV/c)
    std::vector<Float_t> fPtMax = {0.15, 0.25, 0.35};  // (in GeV/c)
    std::vector<Float_t> fYMin = {-1.25, -0.85, -0.7};
    std::vector<Float_t> fYMax = {1.25, 0.85, 0.7};
    Float_t fPhiMin = 0.;               // (in radians)
    Float_t fPhiMax = 2 * TMath::Pi();  // (in radians)

    // define output file
    std::ofstream fOutputFile;
    fOutputFile.open(fOutputFilename);

    // (debug)
    printf(">> Settings:\n");
    printf("   (for pions)   -> %.2f < Pt < %.2f GeV/c\n", fPtMin[0], fPtMax[0]);
    printf("                 -> %.2f < Rapidity < %.2f\n", fYMin[0], fYMax[0]);
    printf("   (for kaons)   -> %.2f < Pt < %.2f GeV/c\n", fPtMin[1], fPtMax[1]);
    printf("                 -> %.2f < Rapidity < %.2f\n", fYMin[1], fYMax[1]);
    printf("   (for protons) -> %.2f < Pt < %.2f GeV/c\n", fPtMin[2], fPtMax[2]);
    printf("                 -> %.2f < Rapidity < %.2f\n", fYMin[2], fYMax[2]);
    printf("   %.2f < Phi < %.2f rad\n\n", fPhiMin, fPhiMax);

    // init PDG database
    TDatabasePDG pdg;

    // prepare random numbers
    gRandom->SetSeed(0);
    Float_t fRandom[3];

    // declare variables
    Int_t current_particle_pdg;
    Float_t Rapidity, Phi;
    Float_t M, Mt, E;
    Float_t Pt, Px, Py, Pl;
    TLorentzVector Particle;
    TString auxStr;

    for (Int_t sign : {1, -1}) {

        for (Int_t pp = 0; pp < (Int_t)pdgVector.size(); pp++) {

            for (Int_t i = 0; i < n_times; i++) {

                current_particle_pdg = sign * pdgVector[pp];

                gRandom->RndmArray(3, fRandom);

                // set particle properties
                Pt = fPtMin[pp] + fRandom[0] * (fPtMax[pp] - fPtMin[pp]);     // Pt (uniform distribution) (in GeV/c)
                Rapidity = fYMin[pp] + fRandom[1] * (fYMax[pp] - fYMin[pp]);  // rapidity (uniform distribution)
                Phi = fPhiMin + fRandom[2] * (fPhiMax - fPhiMin);             // azimuthal angle (uniform distribution) (in radians)
                M = pdg.GetParticle(current_particle_pdg)->Mass();            // mass (in GeV/c^2)
                Px = Pt * TMath::Cos(Phi);                                    // Px (in GeV/c)
                Py = Pt * TMath::Sin(Phi);                                    // Py (in GeV/c)
                Mt = TMath::Sqrt(M * M + Pt * Pt);  // transverse msass (derived from inv. mass and Pt) (in GeV/c^2)
                Pl = Mt * TMath::SinH(Rapidity);    // longitudinal momentum = Pz (derived from trans. mass and Y) (in GeV/c)
                E = Mt * TMath::CosH(Rapidity);     // energy (derived from trans. mass and Y) (in GeV)

                Particle.SetPxPyPzE(Px, Py, Pl, E);

                // print particle onto CSV file
                auxStr = Form("%i,%.8e,%.8e,%.8e\n", current_particle_pdg, Particle.Px(), Particle.Py(), Particle.Pz());

                // (debug)
                printf("%s", auxStr.Data());

                // (output)
                fOutputFile << auxStr;
            }
        }
    }

    fOutputFile.close();

    // (debug)
    printf("\n>> %s has been generated\n", fOutputFilename.Data());
}
