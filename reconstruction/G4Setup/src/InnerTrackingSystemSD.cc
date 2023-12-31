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
/// \file G4Setup/src/InnerTrackingSystemSD.cc
/// \brief Implementation of the B2a::InnerTrackingSystemSD class

#include "InnerTrackingSystemSD.hh"

#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4VProcess.hh"
#include "G4ios.hh"

namespace B2a {

InnerTrackingSystemSD::InnerTrackingSystemSD(G4String name)  //
    : G4VSensitiveDetector(name) {
    collectionName.insert("ITS_hitsCollection");
}

InnerTrackingSystemSD::~InnerTrackingSystemSD() {}

void InnerTrackingSystemSD::Initialize(G4HCofThisEvent* hce) {
    // Create hits collection

    fHitsCollection = new InnerTrackingSystemHitsCollection(SensitiveDetectorName, collectionName[0]);

    // Add this collection in hce

    G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
    hce->AddHitsCollection(hcID, fHitsCollection);
}

G4bool InnerTrackingSystemSD::ProcessHits(G4Step* aStep, G4TouchableHistory*) {
    // energy deposit

    G4double edep = aStep->GetTotalEnergyDeposit();
    if (edep == 0.) return false;

    InnerTrackingSystemHit* newHit = new InnerTrackingSystemHit();

    newHit->SetTrackID(aStep->GetTrack()->GetTrackID());
    newHit->SetLayerNb(aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber());
    newHit->SetEdep(edep);
    newHit->SetMomentum(aStep->GetPostStepPoint()->GetMomentum());
    newHit->SetPosition(aStep->GetPostStepPoint()->GetPosition());

    const G4VProcess* process = aStep->GetTrack()->GetCreatorProcess();
    if (process) {
        newHit->SetProcess(process->GetProcessName());
    } else {
        newHit->SetProcess("Unknown");
    }

    fHitsCollection->insert(newHit);

    return true;
}

void InnerTrackingSystemSD::EndOfEvent(G4HCofThisEvent*) {
    if (verboseLevel > 1) {
        G4int nofHits = fHitsCollection->entries();
        // G4cout << G4endl << "-------->Hits Collection: in this event they are " << nofHits << " hits in the tracker chambers: " <<
        // G4endl; for (G4int i = 0; i < nofHits; i++) (*fHitsCollection)[i]->Print();
    }
}

}  // namespace B2a
