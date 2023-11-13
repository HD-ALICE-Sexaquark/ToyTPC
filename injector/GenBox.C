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

    TDatabasePDG pdg;

    // define output file
    std::ofstream fOutputFile;
    fOutputFile.open(fOutputFilename);

    // params
    std::vector<Int_t> pdgVector = {211, 321, 2212};
    std::vector<Float_t> fPtMin = {0., 0., 0.};  // (in GeV/c)
    std::vector<Float_t> fPtMax = {0.05, 0.1, 0.15};
    Float_t fThetaMin = 0.439;                    // (in rad) ~ 25.15 rad ~ eta = 1.5
    Float_t fThetaMax = TMath::Pi() - fThetaMin;  // (in rad) ~ 154.84 deg ~ eta = -1.5
    Float_t fPhiMin = 0. * TMath::DegToRad();
    Float_t fPhiMax = 360. * TMath::DegToRad();

    printf(">> Settings:\n");
    printf("   (for pions)   %.2f < Pt < %.2f GeV/c\n", fPtMin[0], fPtMax[0]);
    printf("   (for kaons)   %.2f < Pt < %.2f GeV/c\n", fPtMin[1], fPtMax[1]);
    printf("   (for protons) %.2f < Pt < %.2f GeV/c\n", fPtMin[2], fPtMax[2]);
    printf("   %.2f < Theta < %.2f degrees <=> %.2f < Theta < %.2f rad\n",  //
           fThetaMin * TMath::RadToDeg(), fThetaMax * TMath::RadToDeg(), fThetaMin, fThetaMax);
    printf("   %.2f < Phi < %.2f degrees <=> %.2f < Phi < %.2f rad\n\n",  //
           fPhiMin * TMath::RadToDeg(), fPhiMax * TMath::RadToDeg(), fPhiMin, fPhiMax);

    // prepare random numbers
    gRandom->SetSeed(0);
    Float_t fRandom[3];

    for (Int_t sign : {1, -1}) {

        for (Int_t pp = 0; pp < (Int_t)pdgVector.size(); pp++) {

            Int_t current_particle_pdg = sign * pdgVector[pp];

            gRandom->RndmArray(3, fRandom);

            // set particle properties
            Float_t Pt = fPtMin[pp] + fRandom[0] * (fPtMax[pp] - fPtMin[pp]);  // pt (uniform distribution) (in GeV/c)
            Float_t Theta = fThetaMin + fRandom[1] * (fThetaMax - fThetaMin);  // polar angle (uniform distribution) (in radians)
            Float_t Eta = -TMath::Log(TMath::Tan(Theta / 2));                  // eta
            Float_t Phi = fPhiMin + fRandom[2] * (fPhiMax - fPhiMin);          // azimuthal angle (uniform distribution) (in radians)
            Float_t M = pdg.GetParticle(current_particle_pdg)->Mass();             // mass (in GeV/c^2)

            TLorentzVector Particle;
            Particle.SetPtEtaPhiM(Pt, Eta, Phi, M);

            // print particle onto CSV file
            TString auxStr = Form("%i,%.8e,%.8e,%.8e\n", current_particle_pdg, Particle.Px(), Particle.Py(), Particle.Pz());

            // (debug)
            printf("%s", auxStr.Data());

            // (output)
            fOutputFile << auxStr;
        }
    }

    fOutputFile.close();
    printf("\n>> %s has been generated\n", fOutputFilename.Data());
}
