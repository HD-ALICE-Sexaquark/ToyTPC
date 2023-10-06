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
/// \file G4Setup/src/DetectorMessenger.cc
/// \brief Implementation of the B2a::DetectorMessenger class

#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"

#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"

extern std::string input_file;
extern std::string output_file;

namespace B2a {

DetectorMessenger::DetectorMessenger(DetectorConstruction* det) : fDetectorConstruction(det) {

    fALICE = new G4UIdirectory("/ALICE/");
    fALICE->SetGuidance("File management options");

    fInputFileCmd = new G4UIcmdWithAString("/ALICE/input_file", this);
    fInputFileCmd->SetGuidance("Select input file");
    fInputFileCmd->SetParameterName("filename", false);
    fInputFileCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fOutputFileCmd = new G4UIcmdWithAString("/ALICE/output_file", this);
    fOutputFileCmd->SetGuidance("Select output file");
    fOutputFileCmd->SetParameterName("filename", false);
    fOutputFileCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

DetectorMessenger::~DetectorMessenger() {
    delete fInputFileCmd;
    delete fOutputFileCmd;
    delete fALICE;
}

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue) {

    if (command == fInputFileCmd) {
        input_file = newValue;
    }
    if (command == fOutputFileCmd) {
        output_file = newValue;
    }
}

}  // namespace B2a
