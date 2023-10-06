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
/// \file G4Setup/src/DetectorConstruction.cc
/// \brief Implementation of the B2a::DetectorConstruction class

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "InnerTrackingSystemSD.hh"
#include "TimeProjectionChamberSD.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"

#include "G4AutoDelete.hh"
#include "G4Box.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Polycone.hh"
#include "G4Tubs.hh"

#include "G4GeometryManager.hh"
#include "G4GeometryTolerance.hh"

#include "G4UserLimits.hh"

#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include "G4SystemOfUnits.hh"

namespace B2a {

G4ThreadLocal G4GlobalMagFieldMessenger* DetectorConstruction::fMagFieldMessenger = nullptr;

DetectorConstruction::DetectorConstruction() { fMessenger = new DetectorMessenger(this); }

DetectorConstruction::~DetectorConstruction() {
    delete fStepLimit;
    delete fMessenger;
}

G4VPhysicalVolume* DetectorConstruction::Construct() {
    DefineMaterials();
    return DefineVolumes();
}

void DetectorConstruction::DefineMaterials() {
    // Material definition

    G4NistManager* nistManager = G4NistManager::Instance();

    // Air defined using NIST Manager
    nistManager->FindOrBuildMaterial("G4_AIR");
    nistManager->FindOrBuildMaterial("G4_Ne");
    nistManager->FindOrBuildMaterial("G4_CARBON_DIOXIDE");
    nistManager->FindOrBuildMaterial("G4_N");
    nistManager->FindOrBuildMaterial("G4_Si");

    /*
    // (optional) increase Silicon density
    G4Material* RegularSi = nistManager->FindOrBuildMaterial("G4_Si");
    fChamberMaterial = new G4Material("DenseSi",                         //
                                      RegularSi->GetElement(0)->GetZ(),  //
                                      RegularSi->GetElement(0)->GetA(),  //
                                      10. * RegularSi->GetDensity());
    */

    // Print materials
    G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

G4VPhysicalVolume* DetectorConstruction::DefineVolumes() {

    /* Define materials */
    // (pending: could relocate)

    G4Material* Air = G4Material::GetMaterial("G4_AIR");
    G4Material* Ne = G4Material::GetMaterial("G4_Ne");
    G4Material* CO2 = G4Material::GetMaterial("G4_CARBON_DIOXIDE");
    G4Material* N2 = G4Material::GetMaterial("G4_N");

    G4double Ne_comp = 90. / 105.;
    G4double CO2_comp = 10. / 105.;
    G4double N2_comp = 5. / 105.;
    G4double DriftGas_density = Ne_comp * Ne->GetDensity() + CO2_comp * CO2->GetDensity() + N2_comp * N2->GetDensity();

    G4Material* DriftGas = new G4Material("Ne+CO2+N2", DriftGas_density, 3);
    DriftGas->AddMaterial(Ne, Ne_comp * perCent);
    DriftGas->AddMaterial(CO2, CO2_comp * perCent);
    DriftGas->AddMaterial(N2, N2_comp * perCent);

    G4Material* ITSMaterial = G4Material::GetMaterial("G4_Si");

    /*** World ***/

    G4double worldLengthX = 600. * cm;
    G4double worldLengthY = 600. * cm;
    G4double worldLengthZ = 600. * cm;

    G4GeometryManager::GetInstance()->SetWorldMaximumExtent(worldLengthZ);

    G4cout << "Computed tolerance = " << G4GeometryTolerance::GetInstance()->GetSurfaceTolerance() / mm << " mm" << G4endl;

    G4Box* World_S = new G4Box("World_S", 0.5 * worldLengthX, 0.5 * worldLengthY, 0.5 * worldLengthZ);
    G4LogicalVolume* World_LV = new G4LogicalVolume(World_S, Air, "World_LV");

    G4VPhysicalVolume* World_PV = new G4PVPlacement(nullptr,          // no rotation
                                                    G4ThreeVector(),  // must be at (0,0,0)
                                                    World_LV,         // its logical volume
                                                    "World_PV",       // its name
                                                    nullptr,          // its mother volume
                                                    false,            // no boolean operations
                                                    0,                // copy number
                                                    fCheckOverlaps);  // checking overlaps

    /*** ALICE ITS ***/

    // define an chamber of air inside the TPC, where the ITS will be located
    /*
    G4double chamberLength = 500. * cm;
    G4double chamberRadius = 78.8 * cm;

    G4Tubs* Chamber_S = new G4Tubs("Chamber", 0., chamberRadius, 0.5 * chamberLength, 0. * deg, 360. * deg);
    G4LogicalVolume* Chamber_LV = new G4LogicalVolume(Chamber_S, Air, "Chamber", nullptr, nullptr, nullptr);

    new G4PVPlacement(nullptr,                    // no rotation
                      G4ThreeVector(0., 0., 0.),  // (x,y,z) w.r.t. mother's volume center
                      Chamber_LV,                 // its logical volume
                      "Chamber",                  // its name
                      World_LV,                   // its mother volume
                      false,                      // no boolean operations
                      0,                          // copy number
                      fCheckOverlaps);            // checking overlaps
    */

    // [taken from Table 1.1 of the TDR of the ITS Upgrade]
    const G4int NLayersITS = 7;
    G4double layerLength[NLayersITS] = {271. * mm, 271. * mm, 271. * mm, 843. * mm, 843. * mm, 1475. * mm, 1475. * mm};
    G4double layerInnerRadius[NLayersITS] = {22.4 * mm, 30.1 * mm, 37.8 * mm, 194.4 * mm, 243.9 * mm, 342.3 * mm};
    G4double layerOuterRadius[NLayersITS] = {26.7 * mm, 34.6 * mm, 42.1 * mm, 197.7 * mm, 247. * mm, 345.4 * mm, 394.9 * mm};
    G4double layerMidRadius[NLayersITS];
    for (G4int layerNo = 0; layerNo < NLayersITS; layerNo++)
        layerMidRadius[layerNo] = (layerInnerRadius[layerNo] + layerOuterRadius[layerNo]) / 2.;
    // calculated from the material budget of 0.3% per layer
    // -- radiation length (Si) = 21.82 g cm^-2
    // -- density (Si) = 2.3 g cm^-3
    // => RL / D = 9.49 cm ~ 10 cm
    // => (0.3 / 100) x (RL / D ) = 0.03 cm = 0.3 mm
    G4double layerThickness = 0.3 * mm;

    G4Tubs* Layer_S[NLayersITS];
    G4LogicalVolume* Layer_LV[NLayersITS];

    for (G4int layerNo = 0; layerNo < NLayersITS; layerNo++) {

        Layer_S[layerNo] = new G4Tubs("Layer_S",                                       // name
                                      layerMidRadius[layerNo] - 0.5 * layerThickness,  // inner radius
                                      layerMidRadius[layerNo] + 0.5 * layerThickness,  // outer radius
                                      0.5 * layerLength[layerNo],                      // length in z
                                      0. * deg,                                        // minimum phi
                                      360. * deg);                                     // maximum phi
        Layer_LV[layerNo] = new G4LogicalVolume(Layer_S[layerNo], ITSMaterial, "Layer_LV", nullptr, nullptr, nullptr);

        new G4PVPlacement(nullptr,                                    // no rotation
                          G4ThreeVector(0., 0., 0.),                  // (x,y,z) w.r.t. mother's volume center
                          Layer_LV[layerNo],                          // its logical volume
                          "Layer" + std::to_string(layerNo) + "_PV",  // its name
                          World_LV,                                   // its mother volume
                          false,                                      // no boolean operations
                          layerNo,                                    // copy number
                          fCheckOverlaps);                            // checking overlaps
    }

    /*** ALICE TPC ***/

    // [taken from Table 1 of the TPC paper]
    G4double tpcLength = 500. * cm;
    G4double tpcInnerRadius = 78.8 * cm;
    G4double tpcOuterRadius = 258. * cm;

    G4Tubs* TPC_S = new G4Tubs("TPC_S", tpcInnerRadius, tpcOuterRadius, 0.5 * tpcLength, 0. * deg, 360. * deg);
    G4LogicalVolume* TPC_LV = new G4LogicalVolume(TPC_S, DriftGas, "TPC_LV", nullptr, nullptr, nullptr);

    new G4PVPlacement(nullptr,                    // no rotation
                      G4ThreeVector(0., 0., 0.),  // (x,y,z) w.r.t. mother's volume center
                      TPC_LV,                     // its logical volume
                      "TPC_PV",                   // its name
                      World_LV,                   // its mother volume
                      false,                      // no boolean operations
                      0,                          // copy number
                      fCheckOverlaps);            // checking overlaps

    /*** Visualization Attributes ***/

    G4VisAttributes* boxVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));    // white
    G4VisAttributes* layerVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));  // yellow
    G4VisAttributes* tpcVisAtt = new G4VisAttributes(G4Colour(1.0, 0.647, 0.0));  // orange

    World_LV->SetVisAttributes(boxVisAtt);
    for (G4int layerNo = 0; layerNo < NLayersITS; layerNo++) Layer_LV[layerNo]->SetVisAttributes(layerVisAtt);
    TPC_LV->SetVisAttributes(tpcVisAtt);

    // Below is an example of how to set tracking constraints in a given logical volume
    // Sets a max step length in the tracker region, with G4StepLimiter

    G4double maxStep = 0.1;  // borquez: should it be the layer thickness?
    fStepLimit = new G4UserLimits(maxStep);
    TPC_LV->SetUserLimits(fStepLimit);

    /// Set additional contraints on the track, with G4UserSpecialCuts
    ///
    /// G4double maxLength = 2*trackerLength, maxTime = 0.1*ns, minEkin = 10*MeV;
    /// trackerLV->SetUserLimits(new G4UserLimits(maxStep,
    ///                                           maxLength,
    ///                                           maxTime,
    ///                                           minEkin));

    // always return the physical world
    return World_PV;
}

/*
 Construct sensitive detectors and magnetic field
*/
void DetectorConstruction::ConstructSDandField() {

    /** ITS **/

    InnerTrackingSystemSD* ITS_SD = new InnerTrackingSystemSD("/ITS_SD", "InnerTrackingSystemHitsCollection");
    G4SDManager::GetSDMpointer()->AddNewDetector(ITS_SD);
    SetSensitiveDetector("Layer_LV", ITS_SD, true);  // apply to all LVs named "Layer_LV"

    /** TPC **/

    TimeProjectionChamberSD* TPC_SD = new TimeProjectionChamberSD("/TPC_SD");
    G4SDManager::GetSDMpointer()->AddNewDetector(TPC_SD);
    SetSensitiveDetector("TPC_LV", TPC_SD);

    /** Magnetic Field **/

    // uniform magnetic field is then created automatically if the field value is not zero
    G4ThreeVector fieldValue = G4ThreeVector(0., 0., 0.);  // set to zero
    fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
    fMagFieldMessenger->SetVerboseLevel(1);

    G4AutoDelete::Register(fMagFieldMessenger);
}

void DetectorConstruction::SetMaxStep(G4double maxStep) {
    if ((fStepLimit) && (maxStep > 0.)) fStepLimit->SetMaxAllowedStep(maxStep);
}

void DetectorConstruction::SetCheckOverlaps(G4bool checkOverlaps) { fCheckOverlaps = checkOverlaps; }

}  // namespace B2a
