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
/// \file G4Setup/include/TimeProjectionChamberHit.hh
/// \brief Definition of the B2a::TimeProjectionChamberHit class
//         [copied from example B5]

#ifndef B2aTimeProjectionChamberHit_hh
#define B2aTimeProjectionChamberHit_hh 1

#include "G4Allocator.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4VHit.hh"

class G4AttDef;
class G4AttValue;

namespace B2a {

/*
 TPC hit
 It records:
 - the layer ID
 - the particle time
 - the particle local and global positions
*/
class TimeProjectionChamberHit : public G4VHit {
   public:
    TimeProjectionChamberHit();
    TimeProjectionChamberHit(G4int layerID);
    TimeProjectionChamberHit(const TimeProjectionChamberHit& right) = default;
    ~TimeProjectionChamberHit() override;

    TimeProjectionChamberHit& operator=(const TimeProjectionChamberHit& right) = default;
    G4bool operator==(const TimeProjectionChamberHit& right) const;

    inline void* operator new(size_t);
    inline void operator delete(void* aHit);

    void Draw() override;
    const std::map<G4String, G4AttDef>* GetAttDefs() const override;
    std::vector<G4AttValue>* CreateAttValues() const override;
    void Print() override;

    void SetTrackID(G4int track) { fTrackID = track; };
    G4int GetTrackID() const { return fTrackID; };

    void SetMomentum(G4ThreeVector pxpypz) { fMomentum = pxpypz; };
    G4ThreeVector GetMomentum() const { return fMomentum; };

    void SetEdep(G4double energy) { fEdep = energy; };
    G4double GetEdep() const { return fEdep; };

    void SetProcess(G4String process) { fProcess = process; };
    G4String GetProcess() const { return fProcess; };

    void SetLayerID(G4int z) { fLayerID = z; }
    G4int GetLayerID() const { return fLayerID; }

    void SetTime(G4double t) { fTime = t; }
    G4double GetTime() const { return fTime; }

    void SetLocalPos(G4ThreeVector xyz) { fLocalPos = xyz; }
    G4ThreeVector GetLocalPos() const { return fLocalPos; }

    void SetWorldPos(G4ThreeVector xyz) { fWorldPos = xyz; }
    G4ThreeVector GetWorldPos() const { return fWorldPos; }

   private:
    G4int fTrackID;
    G4ThreeVector fMomentum;
    G4double fEdep;
    G4String fProcess;
    G4int fLayerID = -1;
    G4double fTime = 0.;
    G4ThreeVector fLocalPos;
    G4ThreeVector fWorldPos;
};

using TimeProjectionChamberHitsCollection = G4THitsCollection<TimeProjectionChamberHit>;

extern G4ThreadLocal G4Allocator<TimeProjectionChamberHit>* TimeProjectionChamberHitAllocator;

inline void* TimeProjectionChamberHit::operator new(size_t) {
    if (!TimeProjectionChamberHitAllocator) {
        TimeProjectionChamberHitAllocator = new G4Allocator<TimeProjectionChamberHit>;
    }
    return (void*)TimeProjectionChamberHitAllocator->MallocSingle();
}

inline void TimeProjectionChamberHit::operator delete(void* aHit) {
    TimeProjectionChamberHitAllocator->FreeSingle((TimeProjectionChamberHit*)aHit);
}

}  // namespace B2a

#endif
