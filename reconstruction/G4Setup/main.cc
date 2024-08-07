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
//
/// \file main.cc
/// \brief Main program of the B2a example

#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"

#include "FTFP_BERT_HP.hh"

#include "G4Run.hh"
#include "G4RunManagerFactory.hh"
#include "G4StateManager.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4SteppingVerbose.hh"
#include "G4SystemOfUnits.hh"
#include "G4UImanager.hh"

#include "Randomize.hh"

#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"

// global variables
// (bad practice, but who cares at this point)
std::string input_file = "";
std::string traj_file = "";
std::string its_file = "";
std::string tpc_file = "";
G4int pdg_single_part = 0;
G4float py_single_part = 0.;
bool cli_trigger_condition = false;

int main(int argc, char** argv) {

    /* Process input command-line arguments */

    G4UIExecutive* ui = nullptr;

    G4String input_filename;
    G4String output_filename_traj;
    G4String output_filename_ITS;
    G4String output_filename_TPC;
    G4float magfield_z;

    if (argc == 1) {
        ui = new G4UIExecutive(argc, argv);
    } else if (argc == 7) {
        input_filename = argv[1];
        output_filename_traj = argv[2];
        output_filename_ITS = argv[3];
        output_filename_TPC = argv[4];
        cli_trigger_condition = std::stoi(argv[5]);
        magfield_z = std::stof(argv[6]) * tesla;
    } else {
        // clang-format off
        G4cerr << "main.cc :: ERROR: incorrect number of arguments." << G4endl;
        G4cerr << "main.cc ::        -> for command-line mode, you need exactly four arguments:" << G4endl;
        G4cerr << "main.cc ::           ./main <input_filename> <output_filename_traj> <output_filename_ITS> <output_filename_TPC> <trigger_condition 1/0> <magfield_z in tesla>" << G4endl;
        G4cerr << "main.cc ::        -> for graphic-interactive mode, you need no arguments:" << G4endl;
        G4cerr << "main.cc ::           ./main" << G4endl;
        return 1;
        // clang-format on
    }

    // (debug)
    G4cout << "main.cc :: initiating..." << G4endl;
    G4cout << "main.cc :: >> input_filename        = " << input_filename << G4endl;
    G4cout << "main.cc :: >> output_filename_traj  = " << output_filename_traj << G4endl;
    G4cout << "main.cc :: >> output_filename_ITS   = " << output_filename_ITS << G4endl;
    G4cout << "main.cc :: >> output_filename_TPC   = " << output_filename_TPC << G4endl;
    G4cout << "main.cc :: >> cli_trigger_condition = " << cli_trigger_condition << G4endl;
    G4cout << "main.cc :: >> magfield_z            = " << magfield_z / tesla << " tesla" << G4endl;

    // optional: choose a different random engine...
    // G4Random::setTheEngine(new CLHEP::MTwistEngine);

    // use G4SteppingVerboseWithUnits
    G4int precision = 4;
    G4SteppingVerbose::UseBestUnit(precision);

    // Construct the default run manager
    auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);

    // Set mandatory initialization classes
    runManager->SetUserInitialization(new B2a::DetectorConstruction());

    // Select a Physics List and augment it
    G4VModularPhysicsList* physicsList = new FTFP_BERT_HP;
    physicsList->RegisterPhysics(new G4StepLimiterPhysics());

    // Set user action classes
    runManager->SetUserInitialization(physicsList);
    runManager->SetUserInitialization(new B2a::ActionInitialization());

    // Initialize visualization
    G4VisManager* visManager = new G4VisExecutive;
    // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
    // G4VisManager* visManager = new G4VisExecutive("Quiet");
    visManager->Initialize();

    // Get the pointer to the User Interface manager
    G4UImanager* UImanager = G4UImanager::GetUIpointer();

    // Useful objects to check for Trigger Condition
    const G4Run* run = nullptr;
    const std::vector<const G4Event*>* events = nullptr;
    G4int nKeptEvents;
    G4String nProcesses = getenv("GEANT4_NUM_THREADS") ? getenv("GEANT4_NUM_THREADS") : "1";  // left as a string on purpose
    G4int nAttempts = int(1000. / std::stof(nProcesses));                                     // hardcoded value, so job doesn't run forever

    // Process macro or start UI session
    if (!ui) {
        // batch mode
        UImanager->ApplyCommand("/run/numberOfThreads " + nProcesses);
        UImanager->ApplyCommand("/run/initialize");  // G4RunManager::Initialize().
        UImanager->ApplyCommand("/tracking/verbose 0");
        UImanager->ApplyCommand("/tracking/storeTrajectory 2");  // IMPORTANT!!
        UImanager->ApplyCommand("/ALICE/input_file " + input_filename);
        UImanager->ApplyCommand("/ALICE/traj_file " + output_filename_traj);
        UImanager->ApplyCommand("/ALICE/its_file " + output_filename_ITS);
        UImanager->ApplyCommand("/ALICE/tpc_file " + output_filename_TPC);
        UImanager->ApplyCommand("/globalField/setValue 0 0 " + std::to_string(magfield_z / tesla) + " tesla");
        G4int nAttempt_i = 0;
        do {
            UImanager->ApplyCommand("/run/beamOn " + nProcesses);
            run = runManager ? runManager->GetCurrentRun() : nullptr;
            events = run ? run->GetEventVector() : nullptr;
            nKeptEvents = events ? (G4int)events->size() : 0;
            nAttempt_i++;
        } while (nAttempt_i < nAttempts && nKeptEvents < 1 && cli_trigger_condition);
    } else {
        // graphical mode
        UImanager->ApplyCommand("/control/execute init_vis.mac");
        if (ui->IsGUI()) {
            UImanager->ApplyCommand("/control/execute gui.mac");
        }
        ui->SessionStart();
        delete ui;
    }

    // Job termination
    // Free the store: user actions, physics_list and detector_description are
    // owned and deleted by the run manager, so they should not be deleted
    // in the main() program !
    delete visManager;
    delete runManager;
}
