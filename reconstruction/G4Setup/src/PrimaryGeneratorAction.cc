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
/// \file G4Setup/src/PrimaryGeneratorAction.cc
/// \brief Implementation of the B2a::PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"

#include "G4Box.hh"
#include "G4Event.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

extern std::string input_file;
extern G4int pdg_single_part;
extern G4float py_single_part;

namespace B2a {

PrimaryGeneratorAction::PrimaryGeneratorAction() {
    G4int nofParticles = 1;
    fParticlesGun = new G4ParticleGun(nofParticles);
}

PrimaryGeneratorAction::~PrimaryGeneratorAction() { delete fParticlesGun; }

void PrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent) {

    if (pdg_single_part != 0 && py_single_part != 0.) {

        /* Let Geant4 create the particles */

        fParticlesGun = new G4ParticleGun(1);
        fParticlesGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle(pdg_single_part));
        fParticlesGun->SetParticleMomentum(G4ThreeVector(0., py_single_part, 0.));
        fParticlesGun->SetParticlePosition(G4ThreeVector(0., 0., 0.));
        fParticlesGun->GeneratePrimaryVertex(anEvent);

    } else {

        /* Create particles from CSV file */

        std::vector<int> mcStatus, mcPdgCode, mcMotherPdgCode;
        std::vector<double> mcPx, mcPy, mcPz;

        G4String input_filename = "../mc.csv";  // default test value
        if (input_file != "") input_filename = input_file;

        std::ifstream mcFile(input_filename);
        if (!mcFile.is_open()) {
            G4cerr << "ERROR: input file " << input_filename << " not found." << G4endl;
            return;
        }

        std::string line;

        while (std::getline(mcFile, line)) {

            // protection
            if (line == "") continue;

            std::istringstream iss(line);
            std::string token;

            // Read each column separated by commas
            std::getline(iss, token, ',');

            // (debug)
            // std::cout << "TOKEN: " << token << std::endl;

            mcStatus.push_back(std::stoi(token));

            std::getline(iss, token, ',');
            mcPdgCode.push_back(std::stoi(token));

            std::getline(iss, token, ',');
            mcPx.push_back(std::stod(token));

            std::getline(iss, token, ',');
            mcPy.push_back(std::stod(token));

            std::getline(iss, token, ',');
            mcPz.push_back(std::stod(token));
        }

        mcFile.close();

        /* Gun */

        fParticlesGun = new G4ParticleGun(1);

        for (int i = 0; i < (int)mcStatus.size(); i++) {

            fParticlesGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle(mcPdgCode[i]));
            fParticlesGun->SetParticleMomentum(G4ThreeVector(mcPx[i] * GeV, mcPy[i] * GeV, mcPz[i] * GeV));
            fParticlesGun->SetParticlePosition(G4ThreeVector(0., 0., 0.));

            fParticlesGun->GeneratePrimaryVertex(anEvent);
        }
    }
}
}  // namespace B2a
