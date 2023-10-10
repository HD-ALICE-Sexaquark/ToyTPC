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
/// \file G4Setup/include/InnerTrackingSystemSD.hh
/// \brief Definition of the B2a::InnerTrackingSystemSD class

#ifndef B2aInnerTrackingSystemSD_hh
#define B2aInnerTrackingSystemSD_hh 1

#include "G4VSensitiveDetector.hh"

#include "InnerTrackingSystemHit.hh"

#include <vector>

class G4Step;
class G4HCofThisEvent;

namespace B2a {

/*
 Tracker sensitive detector class
 The hits are accounted in hits in ProcessHits() function which is called
 by Geant4 kernel at each step. A hit is created with each step with non zero
 energy deposit.
*/
class InnerTrackingSystemSD : public G4VSensitiveDetector {
   public:
    InnerTrackingSystemSD(G4String name);
    virtual ~InnerTrackingSystemSD();

    // methods from base class
    virtual void Initialize(G4HCofThisEvent* hitCollection);
    virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);
    virtual void EndOfEvent(G4HCofThisEvent* hitCollection);

   private:
    InnerTrackingSystemHitsCollection* fHitsCollection;
};

}  // namespace B2a

#endif
