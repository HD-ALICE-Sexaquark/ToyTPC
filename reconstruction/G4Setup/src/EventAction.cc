//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file G4Setup/src/EventAction.cc
/// \brief Implementation of the B2a::EventAction class

#include "EventAction.hh"
#include "InnerTrackingSystemHit.hh"
#include "TimeProjectionChamberHit.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Trajectory.hh"
#include "G4TrajectoryContainer.hh"
#include "G4ios.hh"

extern std::string output_file;

namespace B2a {

void EventAction::BeginOfEventAction(const G4Event*) {}

void EventAction::EndOfEventAction(const G4Event* event) {
    /*
        G4int fBkgPdgCode = -2112;  // default value: an anti-neutron
        if (bkg_pdg_code) fBkgPdgCode = bkg_pdg_code;

        // get important objects
        auto eventManager = G4EventManager::GetEventManager();

        // get number of stored trajectories
        G4TrajectoryContainer* trajectoryContainer = event->GetTrajectoryContainer();
        G4int n_trajectories = 0;
        if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();

        G4VHitsCollection* hc = event->GetHCofThisEvent()->GetHC(0);
        G4int n_hits = hc->GetSize();

        // periodic printing
        G4int eventID = event->GetEventID();
    */
    /*
    if (eventID < 100 || eventID % 100 == 0) {
        G4cout << ">>> Event: " << eventID << G4endl;
        if (trajectoryContainer) {
            G4cout << "    " << n_trajectories << " trajectories stored" << G4endl;
        }
        G4cout << "    " << n_hits << " hits stored" << G4endl;
    }
    */
    /*
        // declare maps, key = track ID
        std::map<G4int, G4int> NHits;       // number of hits of a track
        std::map<G4int, G4String> Process;  // creation process

        for (size_t i = 0; i < hc->GetSize(); i++) {

            TrackerHit* th = (TrackerHit*)hc->GetHit(i);
            // th->Print();
            Process[th->GetTrackID()] = th->GetProcess();
            NHits[th->GetTrackID()]++;
        }

        // count how many daughters a trajectory had
        std::map<G4int, G4int> NDaughters;
        std::map<G4int, G4int> FirstDaughterID;
        std::map<G4int, G4int> LastDaughterID;
        std::map<G4int, G4int> NNeutralDaughters;
        std::map<G4int, G4int> NChargedDaughters;
        std::map<G4int, std::set<G4int>> DaughtersPDG;

        for (G4int i = 0; i < n_trajectories; i++) {

            G4int trackID = (*trajectoryContainer)[i]->GetTrackID();
            G4int parentID = (*trajectoryContainer)[i]->GetParentID();
            G4int pdg = (*trajectoryContainer)[i]->GetPDGEncoding();
            G4int charge = (G4int)(*trajectoryContainer)[i]->GetCharge();

            NDaughters[parentID]++;

            if (!FirstDaughterID[parentID]) {
                // if not filled, assign current track as the first daughter
                FirstDaughterID[parentID] = trackID;
            } else {
                // when filled,
                LastDaughterID[parentID] = trackID;
                // consequence: when NDaughters == 1, ignore LastDaughterID and get ID with FirstDaughterID
            }

            if (!charge) {
                NNeutralDaughters[parentID]++;
            } else {
                NChargedDaughters[parentID]++;
            }

            DaughtersPDG[parentID].insert(pdg);
        }

        // identify injected particles
        // and check if daughter of signal particles are "nice" reaction channels
        std::map<G4int, G4double> InitialMomentum;
        std::map<G4int, G4int> PDGcode;
        G4int injectedSignalA_ID;
        G4int injectedSignalB_ID;
        G4int injectedBkg_ID;
        */

    /*
for (G4int i = 0; i < n_trajectories; i++) {

    G4int trackID = (*trajectoryContainer)[i]->GetTrackID();
    G4int parentID = (*trajectoryContainer)[i]->GetParentID();

    G4double px = (*trajectoryContainer)[i]->GetInitialMomentum().x();
    G4double py = (*trajectoryContainer)[i]->GetInitialMomentum().y();
    G4double pz = (*trajectoryContainer)[i]->GetInitialMomentum().z();
    G4double x = (*trajectoryContainer)[i]->GetPoint(0)->GetPosition().x();
    G4double y = (*trajectoryContainer)[i]->GetPoint(0)->GetPosition().y();
    G4double z = (*trajectoryContainer)[i]->GetPoint(0)->GetPosition().z();
    G4int pdg = (*trajectoryContainer)[i]->GetPDGEncoding();

    InitialMomentum[trackID] = sqrt(px * px + py * py + pz * pz);
    PDGcode[trackID] = pdg;

    // identify injected particles
    if (!parentID) {
        if (pdg == -3122) injectedSignalA_ID = trackID;
        if (pdg == 310) injectedSignalB_ID = trackID;
        if (pdg == fBkgPdgCode) injectedBkg_ID = trackID;
    }

    G4cout << "    Trajectory = " << trackID << ", Parent = " << parentID << ", PDG = " << pdg  //
           << ", (Px,Py,Pz) = (" << px << ", " << py << ", " << pz                              //
           << ", (x,y,z) = (" << x << ", " << y << ", " << z                                    //
           << "), NHits = " << NHits[trackID] << ", NDaughters = " << NDaughters[trackID] << G4endl;
}
    */

    // trigger condition

    /*
        // (the mentioned probabilities also consider the signal particles decaying into their nice channels)
        G4bool Bkg_V0LikeChannel_AntiNeutron = NDaughters[injectedBkg_ID] > 0 &&                                       //
                                               Process[FirstDaughterID[injectedBkg_ID]] == "anti_neutronInelastic" &&  //
                                               DaughtersPDG[injectedBkg_ID].count(310) &&
                                               DaughtersPDG[injectedBkg_ID].count(-3122);                     // (interaction ~ 2%)
        G4bool Bkg_V0LikeChannel_Neutron = NDaughters[injectedBkg_ID] > 0 &&                                  //
                                           Process[FirstDaughterID[injectedBkg_ID]] == "neutronInelastic" &&  //
                                           DaughtersPDG[injectedBkg_ID].count(310) &&
                                           DaughtersPDG[injectedBkg_ID].count(-3122);  // (interaction ~ 2%)
        G4bool Bkg_V0LikeChannel_AntiProton;                                           // (interaction)
        G4bool Bkg_V0LikeChannel_Proton;                                               // (interaction)
        G4bool Bkg_V0LikeChannel_Gamma;                                                // (conversion, used for testing purposes ~ 10%)
        G4bool Bkg_V0LikeChannel_Pi0;                                                  // (discarded! decay + conversion chance less than
       1%) G4bool Bkg_V0LikeChannel_Eta;                                                  // (discarded! decay + conversion chance less than
       1%) G4bool Bkg_V0LikeChannel_K0Short;                                              // (interaction) G4bool Bkg_V0LikeChannel_K0Long =
       NDaughters[injectedBkg_ID] > 0;              // (interaction ~ 2%) G4bool Bkg_V0LikeChannel_Lambda; // (interaction) G4bool
       Bkg_V0LikeChannel_AntiLambda;                                           // (interaction) G4bool Bkg_V0LikeChannel_EtaPrime; // (decay
       + conversion ~ 3%) G4bool Bkg_V0LikeChannel_Phi;  // (yet to see, for this one, the interaction prob. of the K0Long should be >= 10%)

        // choose which one will be applied based on the input bkg code
        G4bool Bkg_V0LikeChannel = true;  // default, no requirement is needed for the bkg
        if (fBkgPdgCode == -2112) Bkg_V0LikeChannel = Bkg_V0LikeChannel_AntiNeutron;
        if (fBkgPdgCode == 2112) Bkg_V0LikeChannel = Bkg_V0LikeChannel_Neutron;
        if (fBkgPdgCode == -2212) Bkg_V0LikeChannel = Bkg_V0LikeChannel_AntiProton;
        if (fBkgPdgCode == 2212) Bkg_V0LikeChannel = Bkg_V0LikeChannel_Proton;
        if (fBkgPdgCode == 22) Bkg_V0LikeChannel = Bkg_V0LikeChannel_Gamma;
        if (fBkgPdgCode == 111) Bkg_V0LikeChannel = Bkg_V0LikeChannel_Pi0;
        if (fBkgPdgCode == 221) Bkg_V0LikeChannel = Bkg_V0LikeChannel_Eta;
        if (fBkgPdgCode == 310) Bkg_V0LikeChannel = Bkg_V0LikeChannel_K0Short;
        if (fBkgPdgCode == 130) Bkg_V0LikeChannel = Bkg_V0LikeChannel_K0Long;
        if (fBkgPdgCode == 3122) Bkg_V0LikeChannel = Bkg_V0LikeChannel_Lambda;
        if (fBkgPdgCode == -3122) Bkg_V0LikeChannel = Bkg_V0LikeChannel_AntiLambda;
        if (fBkgPdgCode == 331) Bkg_V0LikeChannel = Bkg_V0LikeChannel_EtaPrime;
        if (fBkgPdgCode == 333) Bkg_V0LikeChannel = Bkg_V0LikeChannel_Phi;

        G4bool OnlyBkg = signal_file == "0";
        G4bool SignalA_NiceChannel = true;
        G4bool SignalB_NiceChannel = true;
        if (!OnlyBkg) {
            SignalA_NiceChannel = NDaughters[injectedSignalA_ID] == 2 &&  //
                                  DaughtersPDG[injectedSignalA_ID].count(-2212) && DaughtersPDG[injectedSignalA_ID].count(211);
            SignalB_NiceChannel = NDaughters[injectedSignalB_ID] == 2 &&  //
                                  DaughtersPDG[injectedSignalB_ID].count(211) && DaughtersPDG[injectedSignalB_ID].count(-211);
        }
        if (SignalA_NiceChannel && SignalB_NiceChannel && Bkg_V0LikeChannel) {
            StoreEvent(event);
            eventManager->KeepTheCurrentEvent();
        }
        */
}

void EventAction::StoreEvent(const G4Event* event) {
    //
    //
    //

    G4String fOutputFilename = "../output_e" + std::to_string(event->GetEventID()) + ".csv";  // default test value
    if (output_file != "") fOutputFilename = output_file;

    std::ofstream fOutputFile;
    fOutputFile.open(fOutputFilename);

    // get trajectories
    G4TrajectoryContainer* trajectoryContainer = event->GetTrajectoryContainer();
    G4int n_trajectories = 0;
    if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();

    // map trajectories info to trackID
    std::map<G4int, G4int> PdgCode;
    std::map<G4int, G4ThreeVector> InitialMomentum;
    std::map<G4int, G4ThreeVector> InitialPosition;
    std::map<G4int, G4int> MotherID;
    std::map<G4int, G4bool> IsPrimary;
    std::map<G4int, G4bool> IsSignal;

    for (G4int i = 0; i < n_trajectories; i++) {

        G4int trackID = (*trajectoryContainer)[i]->GetTrackID();
        PdgCode[trackID] = (*trajectoryContainer)[i]->GetPDGEncoding();
        InitialPosition[trackID] = (*trajectoryContainer)[i]->GetPoint(0)->GetPosition();
        InitialMomentum[trackID] = (*trajectoryContainer)[i]->GetInitialMomentum();
        MotherID[trackID] = (*trajectoryContainer)[i]->GetParentID();
        IsPrimary[trackID] = MotherID[trackID] == 0;
        G4bool is_secondary = InitialPosition[trackID].x() != 0. &&  //
                              InitialPosition[trackID].y() != 0. &&  //
                              InitialPosition[trackID].z() != 0.;
        G4bool is_signal = IsPrimary[trackID] &&                                      //
                           (PdgCode[trackID] == 310 || PdgCode[trackID] == -3122) &&  //
                           is_secondary;
        G4bool is_mother_secondary = InitialPosition[MotherID[trackID]].x() != 0. &&  //
                                     InitialPosition[MotherID[trackID]].y() != 0. &&  //
                                     InitialPosition[MotherID[trackID]].z() != 0.;
        G4bool is_mother_signal = IsPrimary[MotherID[trackID]] &&                                                //
                                  (PdgCode[MotherID[trackID]] == 310 || PdgCode[MotherID[trackID]] == -3122) &&  //
                                  is_mother_secondary;
        IsSignal[trackID] = is_signal || is_mother_signal;
    }

    // declare columns
    G4int csv_eventID;
    G4int csv_trackID;
    G4int csv_chamberNb;
    G4double csv_PDGcode;
    G4double csv_x;
    G4double csv_y;
    G4double csv_z;
    G4double csv_px;
    G4double csv_py;
    G4double csv_pz;
    G4double csv_x_ini;
    G4double csv_y_ini;
    G4double csv_z_ini;
    G4double csv_px_ini;
    G4double csv_py_ini;
    G4double csv_pz_ini;
    G4double csv_Edep;
    G4String csv_process;
    G4bool csv_issignal;

    G4int csv_motherID;
    G4int csv_mother_PDGcode;
    G4bool csv_mother_issignal;
    G4double csv_mother_x;
    G4double csv_mother_y;
    G4double csv_mother_z;
    G4double csv_mother_px;
    G4double csv_mother_py;
    G4double csv_mother_pz;

    // loop over hits
    /*
        G4VHitsCollection* hc = event->GetHCofThisEvent()->GetHC(0);
        G4int n_hits = hc->GetSize();

        for (G4int i = 0; i < n_hits; i++) {

            TrackerHit* th = (TrackerHit*)hc->GetHit(i);

            csv_eventID = event->GetEventID();
            csv_trackID = th->GetTrackID();
            csv_chamberNb = th->GetChamberNb();

            csv_PDGcode = PdgCode[csv_trackID];
            csv_x = th->GetPosition().x();
            csv_y = th->GetPosition().y();
            csv_z = th->GetPosition().z();
            csv_px = th->GetMomentum().x();
            csv_py = th->GetMomentum().y();
            csv_pz = th->GetMomentum().z();
            csv_x_ini = InitialPosition[csv_trackID].x();
            csv_y_ini = InitialPosition[csv_trackID].y();
            csv_z_ini = InitialPosition[csv_trackID].z();
            csv_px_ini = InitialMomentum[csv_trackID].x();
            csv_py_ini = InitialMomentum[csv_trackID].y();
            csv_pz_ini = InitialMomentum[csv_trackID].z();
            csv_Edep = th->GetEdep();
            csv_process = th->GetProcess();
            csv_issignal = IsSignal[csv_trackID];

            csv_motherID = MotherID[csv_trackID];
            csv_mother_PDGcode = PdgCode[csv_motherID];
            csv_mother_issignal = IsSignal[csv_motherID];
            csv_mother_x = InitialPosition[csv_motherID].x();
            csv_mother_y = InitialPosition[csv_motherID].y();
            csv_mother_z = InitialPosition[csv_motherID].z();
            csv_mother_px = InitialMomentum[csv_motherID].x();
            csv_mother_py = InitialMomentum[csv_motherID].y();
            csv_mother_pz = InitialMomentum[csv_motherID].z();

            // (output)
            fOutputFile << csv_eventID << "," << csv_trackID << "," << csv_chamberNb << ","                               //
                        << (G4long)csv_PDGcode << "," << csv_x / cm << "," << csv_y / cm << "," << csv_z / cm << ","      //
                        << csv_px / GeV << "," << csv_py / GeV << "," << csv_pz / GeV << ","                              //
                        << csv_x_ini / cm << "," << csv_y_ini / cm << "," << csv_z_ini / cm << ","                        //
                        << csv_px_ini / GeV << "," << csv_py_ini / GeV << "," << csv_pz_ini / GeV << ","                  //
                        << csv_Edep / GeV << "," << csv_process << "," << (G4int)csv_issignal << ","                      //
                        << csv_motherID << "," << (G4long)csv_mother_PDGcode << "," << (G4int)csv_mother_issignal << ","  //
                        << csv_mother_x / cm << "," << csv_mother_y / cm << "," << csv_mother_z / cm << ","               //
                        << csv_mother_px / GeV << "," << csv_mother_py / GeV << "," << csv_mother_pz / GeV << G4endl;
        }
         */
    fOutputFile.close();
}

}  // namespace B2a
