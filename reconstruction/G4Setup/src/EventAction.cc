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

extern std::string traj_file;
extern std::string its_file;
extern std::string tpc_file;

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

    /* Debug --Trajectories */

    G4TrajectoryContainer* trajectoryContainer = event->GetTrajectoryContainer();
    G4int n_trajectories = 0;
    if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
    G4cout << "  >> N Trajectories: " << n_trajectories << G4endl;
    /*
        for (G4int i = 0; i < n_trajectories; i++) {
            G4cout << "     >> Trajectory " << i                               //
                   << ", Track " << (*trajectoryContainer)[i]->GetTrackID()    //
                   << ", Parent " << (*trajectoryContainer)[i]->GetParentID()  //
                   << ", PDG " << (*trajectoryContainer)[i]->GetPDGEncoding()  //
                   << ", Charge " << (G4int)(*trajectoryContainer)[i]->GetCharge() << G4endl;
        }
    */

    /* Debug -- ITS Hits */

    G4VHitsCollection* itsHC = GetHC(event, itsHC_id);
    G4int n_its_hits = 0;
    if (itsHC) n_its_hits = itsHC->GetSize();
    G4cout << "  >> N ITS Hits: " << n_its_hits << G4endl;
    /*
        InnerTrackingSystemHit* its_hit = nullptr;
        for (G4int i = 0; i < n_its_hits; i++) {
            its_hit = (InnerTrackingSystemHit*)itsHC->GetHit(i);
            G4cout << "     >> Hit " << i                  //
                   << ", Track " << its_hit->GetTrackID()  //
                   << ", Layer " << its_hit->GetLayerNb() << G4endl;
        }
    */

    /* Debug -- TPC Hits */

    G4VHitsCollection* tpcHC = GetHC(event, tpcHC_id);
    G4int n_tpc_hits = 0;
    if (tpcHC) n_tpc_hits = tpcHC->GetSize();
    G4cout << "  >> N TPC Hits: " << n_tpc_hits << G4endl;
    /*
        TimeProjectionChamberHit* tpc_hit = nullptr;
        for (G4int i = 0; i < n_tpc_hits; i++) {
            tpc_hit = (TimeProjectionChamberHit*)tpcHC->GetHit(i);
            G4cout << "     >> Hit " << i                  //
                   << ", Track " << tpc_hit->GetTrackID()  //
                   << ", Time " << tpc_hit->GetTime() << G4endl;
        }
    */

    // there's no trigger condition
    StoreEvent(event);
    eventManager->KeepTheCurrentEvent();
}

/*
 Store event info into .CSV files
 */
void EventAction::StoreEvent(const G4Event* event) {

    G4int eventID = event->GetEventID();

    /* Trajectories */

    G4String fTrajectoriesFilename = "../event" + std::to_string(eventID) + "_traj.csv";  // default test value
    if (traj_file != "") fTrajectoriesFilename = traj_file;

    std::ofstream fTrajectoriesFile;
    fTrajectoriesFile.open(fTrajectoriesFilename);

    G4TrajectoryContainer* trajectoryContainer = event->GetTrajectoryContainer();
    G4int n_trajectories = 0;
    if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();

    G4int traj_trackID;
    G4long traj_PDGcode;
    G4ThreeVector traj_position;
    G4double traj_x_ini, traj_y_ini, traj_z_ini;
    G4ThreeVector traj_momentum;
    G4double traj_px_ini, traj_py_ini, traj_pz_ini;
    G4int traj_parentID;
    G4int traj_charge;

    for (G4int i = 0; i < n_trajectories; i++) {

        traj_trackID = (*trajectoryContainer)[i]->GetTrackID();
        traj_PDGcode = (G4long)(*trajectoryContainer)[i]->GetPDGEncoding();
        traj_position = (*trajectoryContainer)[i]->GetPoint(0)->GetPosition();
        traj_x_ini = traj_position.x();
        traj_y_ini = traj_position.y();
        traj_z_ini = traj_position.z();
        traj_momentum = (*trajectoryContainer)[i]->GetInitialMomentum();
        traj_px_ini = traj_momentum.x();
        traj_py_ini = traj_momentum.y();
        traj_pz_ini = traj_momentum.z();
        traj_parentID = (*trajectoryContainer)[i]->GetParentID();
        traj_charge = (G4int)(*trajectoryContainer)[i]->GetCharge();

        fTrajectoriesFile << eventID << "," << traj_trackID << "," << traj_PDGcode << ","                      //
                          << traj_x_ini / cm << "," << traj_y_ini / cm << "," << traj_z_ini / cm << ","        //
                          << traj_px_ini / MeV << "," << traj_py_ini / MeV << "," << traj_pz_ini / MeV << ","  //
                          << traj_parentID << "," << traj_charge << G4endl;
    }

    fTrajectoriesFile.close();

    /* ITS Hits */

    G4String fITSFilename = "../event" + std::to_string(eventID) + "_its.csv";  // default test value
    if (its_file != "") fITSFilename = its_file;

    std::ofstream fITSFile;
    fITSFile.open(fITSFilename);

    G4VHitsCollection* itsHC = GetHC(event, itsHC_id);
    G4int n_its_hits = 0;
    if (itsHC) n_its_hits = itsHC->GetSize();

    InnerTrackingSystemHit* its_hit = nullptr;
    G4int its_trackID;
    G4int its_layer_n;
    G4ThreeVector its_position;
    G4double its_x, its_y, its_z;
    G4double its_dedx;
    G4String its_process;

    for (G4int i = 0; i < n_its_hits; i++) {

        its_hit = (InnerTrackingSystemHit*)itsHC->GetHit(i);

        its_trackID = its_hit->GetTrackID();
        its_layer_n = its_hit->GetLayerNb();
        its_position = its_hit->GetPosition();
        its_x = its_position.x();
        its_y = its_position.y();
        its_z = its_position.z();
        its_dedx = its_hit->GetEdep();
        its_process = its_hit->GetProcess();

        fITSFile << eventID << "," << its_trackID << "," << its_layer_n << ","   //
                 << its_x / cm << "," << its_y / cm << "," << its_z / cm << ","  //
                 << its_dedx / MeV << "," << its_process << G4endl;
    }

    fITSFile.close();

    /* TPC Hits */

    G4String fTPCFilename = "../event" + std::to_string(eventID) + "_tpc.csv";  // default test value
    if (tpc_file != "") fTPCFilename = tpc_file;

    std::ofstream fTPCFile;
    fTPCFile.open(fTPCFilename);

    G4VHitsCollection* tpcHC = GetHC(event, tpcHC_id);
    G4int n_tpc_hits = 0;
    if (tpcHC) n_tpc_hits = tpcHC->GetSize();

    TimeProjectionChamberHit* tpc_hit = nullptr;

    G4int tpc_trackID;
    G4ThreeVector tpc_position;
    G4double tpc_x, tpc_y, tpc_z;
    G4ThreeVector tpc_momentum;
    G4double tpc_px, tpc_py, tpc_pz;
    G4double tpc_time;
    G4double tpc_dedx;
    G4String tpc_process;

    for (G4int i = 0; i < n_tpc_hits; i++) {

        tpc_hit = (TimeProjectionChamberHit*)tpcHC->GetHit(i);

        tpc_trackID = tpc_hit->GetTrackID();
        tpc_position = tpc_hit->GetWorldPos();
        tpc_x = tpc_position.x();
        tpc_y = tpc_position.y();
        tpc_z = tpc_position.z();
        tpc_momentum = tpc_hit->GetMomentum();
        tpc_px = tpc_momentum.x();
        tpc_py = tpc_momentum.y();
        tpc_pz = tpc_momentum.z();
        tpc_time = tpc_hit->GetTime();
        tpc_dedx = tpc_hit->GetEdep();
        tpc_process = tpc_hit->GetProcess();

        fTPCFile << eventID << "," << tpc_trackID << ","                               //
                 << tpc_x / cm << "," << tpc_y / cm << "," << tpc_z / cm << ","        //
                 << tpc_px / MeV << "," << tpc_py / MeV << "," << tpc_pz / MeV << ","  //
                 << tpc_time / nanosecond << "," << tpc_dedx / MeV << "," << tpc_process << G4endl;
    }

    fTPCFile.close();
}

}  // namespace B2a
