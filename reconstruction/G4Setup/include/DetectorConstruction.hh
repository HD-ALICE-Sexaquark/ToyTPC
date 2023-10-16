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
/// \file G4Setup/include/DetectorConstruction.hh
/// \brief Definition of the B2a::DetectorConstruction class

#ifndef B2aDetectorConstruction_hh
#define B2aDetectorConstruction_hh 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "tls.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;
class G4UserLimits;
class G4GlobalMagFieldMessenger;

namespace B2a {

class DetectorMessenger;

/// Detector construction class to define materials, geometry
/// and global uniform magnetic field.

class DetectorConstruction : public G4VUserDetectorConstruction {
   public:
    DetectorConstruction();
    ~DetectorConstruction() override;

   public:
    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;

    // Set methods
    void SetMaxStep(G4double);
    void SetCheckOverlaps(G4bool);

   private:
    void DefineMaterials();
    G4VPhysicalVolume* DefineVolumes();

    static G4ThreadLocal G4GlobalMagFieldMessenger* fMagFieldMessenger;

    G4UserLimits* fStepLimit_ITS = nullptr;  // pointer to user step limits
    G4UserLimits* fStepLimit_TPC = nullptr;  // pointer to user step limits

    DetectorMessenger* fMessenger = nullptr;  // messenger

    G4bool fCheckOverlaps = true;  // option to activate checking of volumes overlaps
};

}  // namespace B2a

#endif
