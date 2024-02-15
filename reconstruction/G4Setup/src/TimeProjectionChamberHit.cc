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
/// \file G4Setup/src/TimeProjectionChamberHit.cc
/// \brief Implementation of the B2a::TimeProjectionChamberHit class
//         [copied from example B5]

#include "TimeProjectionChamberHit.hh"

#include "G4AttDef.hh"
#include "G4AttDefStore.hh"
#include "G4AttValue.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4SystemOfUnits.hh"
#include "G4UIcommand.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"
#include "G4ios.hh"

namespace B2a {

G4ThreadLocal G4Allocator<TimeProjectionChamberHit>* TimeProjectionChamberHitAllocator;

TimeProjectionChamberHit::TimeProjectionChamberHit() {}

TimeProjectionChamberHit::TimeProjectionChamberHit(G4int layerID) : fLayerID(layerID) {}

TimeProjectionChamberHit::~TimeProjectionChamberHit() {}

// TimeProjectionChamberHit::TimeProjectionChamberHit(const TimeProjectionChamberHit &right)
// : G4VHit(),
//   fLayerID(right.fLayerID),
//   fTime(right.fTime),
//   fLocalPos(right.fLocalPos),
//   fWorldPos(right.fWorldPos)
// {}

// const TimeProjectionChamberHit& TimeProjectionChamberHit::operator=(const TimeProjectionChamberHit &right)
// {
//   fLayerID = right.fLayerID;
//   fTime = right.fTime;
//   fLocalPos = right.fLocalPos;
//   fWorldPos = right.fWorldPos;
//   return *this;
// }

G4bool TimeProjectionChamberHit::operator==(const TimeProjectionChamberHit& /*right*/) const { return false; }

void TimeProjectionChamberHit::Draw() {
    auto visManager = G4VVisManager::GetConcreteInstance();
    if (!visManager) return;

    G4Circle circle(fWorldPos);
    circle.SetScreenSize(2);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(1., 0.5, 0.);  // orange
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    visManager->Draw(circle);
}

const std::map<G4String, G4AttDef>* TimeProjectionChamberHit::GetAttDefs() const {
    G4bool isNew;
    auto store = G4AttDefStore::GetInstance("TimeProjectionChamberHit", isNew);

    if (isNew) {
        (*store)["HitType"] = G4AttDef("HitType", "Hit Type", "Physics", "", "G4String");
        (*store)["ID"] = G4AttDef("ID", "ID", "Physics", "", "G4int");
        (*store)["Time"] = G4AttDef("Time", "Time", "Physics", "G4BestUnit", "G4double");
        (*store)["Pos"] = G4AttDef("Pos", "Position", "Physics", "G4BestUnit", "G4ThreeVector");
    }

    return store;
}

std::vector<G4AttValue>* TimeProjectionChamberHit::CreateAttValues() const {
    auto values = new std::vector<G4AttValue>;

    values->push_back(G4AttValue("HitType", "TimeProjectionChamberHit", ""));
    values->push_back(G4AttValue("ID", G4UIcommand::ConvertToString(fLayerID), ""));
    values->push_back(G4AttValue("Time", G4BestUnit(fTime, "Time"), ""));
    values->push_back(G4AttValue("Pos", G4BestUnit(fWorldPos, "Length"), ""));

    return values;
}

void TimeProjectionChamberHit::Print() {
    G4cout << "  Layer[" << fLayerID << "] : time " << fTime / ns << " (nsec) --- local (x,y) " << fLocalPos.x() << ", " << fLocalPos.y()
           << G4endl;
}

}  // namespace B2a
