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

#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIdirectory.hh"

extern std::string input_file;
extern std::string traj_file;
extern std::string its_file;
extern std::string tpc_file;
extern G4int pdg_single_part;
extern G4float py_single_part;

namespace B2a {

DetectorMessenger::DetectorMessenger(DetectorConstruction* det) : fDetectorConstruction(det) {

    fALICE = new G4UIdirectory("/ALICE/");
    fALICE->SetGuidance("File management options");

    fInputFileCmd = new G4UIcmdWithAString("/ALICE/input_file", this);
    fInputFileCmd->SetGuidance("Select input file");
    fInputFileCmd->SetParameterName("filename", false);
    fInputFileCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fTrajFileCmd = new G4UIcmdWithAString("/ALICE/traj_file", this);
    fTrajFileCmd->SetGuidance("Select output filename for trajectories info");
    fTrajFileCmd->SetParameterName("filename", false);
    fTrajFileCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fITSFileCmd = new G4UIcmdWithAString("/ALICE/its_file", this);
    fITSFileCmd->SetGuidance("Select output filename for ITS info");
    fITSFileCmd->SetParameterName("filename", false);
    fITSFileCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fTPCFileCmd = new G4UIcmdWithAString("/ALICE/tpc_file", this);
    fTPCFileCmd->SetGuidance("Select output filename for TPC info");
    fTPCFileCmd->SetParameterName("filename", false);
    fTPCFileCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fMagFieldCmd = new G4UIcmdWithADoubleAndUnit("/ALICE/mag_field", this);
    fMagFieldCmd->SetGuidance("Change the z-component of the magnetic field");
    fMagFieldCmd->SetParameterName("mag_field", false);
    fMagFieldCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
    fMagFieldCmd->SetUnitCandidates("tesla");

    fSinglePartPDGCmd = new G4UIcmdWithAnInteger("/ALICE/pdg_single_part", this);
    fSinglePartPDGCmd->SetGuidance("Choose the PDG code of the single particle to inject");
    fSinglePartPDGCmd->SetParameterName("pdg_code", false);
    fSinglePartPDGCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    fSinglePartPyCmd = new G4UIcmdWithADoubleAndUnit("/ALICE/py_single_part", this);
    fSinglePartPyCmd->SetGuidance("Choose the momentum's y-component of the single particle to inject");
    fSinglePartPyCmd->SetParameterName("py", false);
    fSinglePartPyCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
    fSinglePartPyCmd->SetUnitCandidates("MeV GeV");
}

DetectorMessenger::~DetectorMessenger() {
    delete fInputFileCmd;
    delete fTrajFileCmd;
    delete fITSFileCmd;
    delete fTPCFileCmd;
    delete fMagFieldCmd;
    delete fSinglePartPDGCmd;
    delete fSinglePartPyCmd;
    delete fALICE;
}

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue) {

    if (command == fInputFileCmd) {
        input_file = newValue;
    }
    if (command == fTrajFileCmd) {
        traj_file = newValue;
    }
    if (command == fITSFileCmd) {
        its_file = newValue;
    }
    if (command == fTPCFileCmd) {
        tpc_file = newValue;
    }
    if (command == fMagFieldCmd) {
        fDetectorConstruction->SetMagneticField(fMagFieldCmd->GetNewDoubleValue(newValue));
    }
    if (command == fSinglePartPDGCmd) {
        pdg_single_part = fSinglePartPDGCmd->GetNewIntValue(newValue);
    }
    if (command == fSinglePartPyCmd) {
        py_single_part = fSinglePartPyCmd->GetNewDoubleValue(newValue);
    }
}

}  // namespace B2a
