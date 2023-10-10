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
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Trajectory.hh"
#include "G4TrajectoryContainer.hh"
#include "G4ios.hh"

extern std::string output_file;

namespace {

/*
 Utility function which finds a hit collection with the given Id
 and print warnings if not found
 [copied from examples/B5/src/EventAction.cc]
*/
G4VHitsCollection* GetHC(const G4Event* event, G4int collId) {

    G4HCofThisEvent* hce = event->GetHCofThisEvent();
    if (!hce) {
        G4ExceptionDescription msg;
        msg << "No hits collection of this event found." << G4endl;
        G4Exception("EventAction::EndOfEventAction()", "Code001", JustWarning, msg);
        return nullptr;
    }

    G4VHitsCollection* hc = hce->GetHC(collId);
    if (!hc) {
        G4ExceptionDescription msg;
        msg << "Hits collection " << collId << " of this event not found." << G4endl;
        G4Exception("EventAction::EndOfEventAction()", "Code001", JustWarning, msg);
    }
    return hc;
}

}  // namespace

namespace B2a {

/*
 Find hit collections IDs by names (just once)
 and save them in the data members of this class
*/
void EventAction::BeginOfEventAction(const G4Event*) {

    auto sdManager = G4SDManager::GetSDMpointer();

    // hit collections names
    G4String itsHCname = "ITS/ITS_hitsCollection";
    G4String tpcHCname = "TPC/TPC_hitsCollection";

    // hit collections IDs
    itsHC_id = sdManager->GetCollectionID(itsHCname);
    tpcHC_id = sdManager->GetCollectionID(tpcHCname);

    G4cout << "itsHC_id = " << itsHC_id << G4endl;
    G4cout << "tpcHC_id = " << tpcHC_id << G4endl;
}

void EventAction::EndOfEventAction(const G4Event* event) {

    // get important objects
    auto eventManager = G4EventManager::GetEventManager();

    G4int eventID = event->GetEventID();
    G4cout << "> Event " << eventID << G4endl;

    /* Trajectories */

    G4TrajectoryContainer* trajectoryContainer = event->GetTrajectoryContainer();
    G4int n_trajectories = 0;
    if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
    G4cout << "  >> N Trajectories: " << n_trajectories << G4endl;

    // debug
    for (G4int i = 0; i < n_trajectories; i++) {
        G4cout << "     >> Trajectory " << i                               //
               << ", Track " << (*trajectoryContainer)[i]->GetTrackID()    //
               << ", Parent " << (*trajectoryContainer)[i]->GetParentID()  //
               << ", PDG " << (*trajectoryContainer)[i]->GetPDGEncoding()  //
               << ", Charge " << (G4int)(*trajectoryContainer)[i]->GetCharge() << G4endl;
    }

    /* ITS Hits */

    G4VHitsCollection* itsHC = GetHC(event, itsHC_id);
    G4int n_its_hits = 0;
    if (itsHC) n_its_hits = itsHC->GetSize();
    G4cout << "  >> N ITS Hits: " << n_its_hits << G4endl;

    // declare maps, key = track ID
    // std::map<G4int, G4int> NHits;       // number of hits of a track
    // std::map<G4int, G4String> Process;  // creation process

    InnerTrackingSystemHit* its_h = nullptr;
    for (G4int i = 0; i < n_its_hits; i++) {
        its_h = (InnerTrackingSystemHit*)itsHC->GetHit(i);
        G4cout << "     >> Hit " << i                //
               << ", Track " << its_h->GetTrackID()  //
               << ", Layer " << its_h->GetChamberNb() << G4endl;
    }

    /* TPC Hits */

    G4VHitsCollection* tpcHC = GetHC(event, tpcHC_id);
    G4int n_tpc_hits = 0;
    if (tpcHC) n_tpc_hits = tpcHC->GetSize();
    G4cout << "  >> N TPC Hits: " << n_tpc_hits << G4endl;

    TimeProjectionChamberHit* tpc_h = nullptr;
    for (G4int i = 0; i < n_tpc_hits; i++) {
        tpc_h = (TimeProjectionChamberHit*)tpcHC->GetHit(i);
        G4cout << "     >> Hit " << i                //
               << ", Layer " << tpc_h->GetLayerID()  //
               << ", Time " << tpc_h->GetTime() << G4endl;
    }

    // there's no trigger condition
    // StoreEvent(event);
    eventManager->KeepTheCurrentEvent();
}

/*
 Store event info into a .CSV file
 */
void EventAction::StoreEvent(const G4Event* event) {

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

    G4VHitsCollection* hc = event->GetHCofThisEvent()->GetHC(0);
    G4int n_hits = hc->GetSize();

    for (G4int i = 0; i < n_hits; i++) {

        InnerTrackingSystemHit* th = (InnerTrackingSystemHit*)hc->GetHit(i);

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

    fOutputFile.close();
}

}  // namespace B2a
