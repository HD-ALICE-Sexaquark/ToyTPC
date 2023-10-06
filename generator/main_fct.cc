// main01.cc is a part of the PYTHIA event generator.
// Copyright (C) 2023 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

#include <fstream>
#include <iostream>
#include <sstream>  // stringstream, useful to convert char to str
#include <string>

#include "Pythia8/HeavyIons.h"
#include "Pythia8/Pythia.h"

using namespace Pythia8;

int main(int argc, char* argv[]) {
    // usage: after compilation, ./main_hi main_hi.cfg

    // declare generator
    Pythia pythia;

    // read config file
    pythia.readFile(argv[1]);
    int nEvent = pythia.mode("Main:numberOfEvents");

    // initialize
    pythia.init();

    // prepare output file
    std::ofstream output_file;

    // string and auxiliary types to help with formatting
    std::string output_filename;
    std::stringstream aux_ss;
    char buffer_a[48];
    char buffer_b[96];

    // aux. variables to hold mother's pid
    int mother1;
    int mother2;
    int id_mother;

    // counter of valid particles
    int n_particles;

    // event loop
    for (int iEvent = 0; iEvent < nEvent; iEvent++) {

        // generate event, skip if error
        if (!pythia.next()) {
            continue;
        }

        // prepare output filename
        if (nEvent == 1) {
            // test value
            output_filename = "../reco/bkg.csv";
        } else {
            aux_ss.clear();
            sprintf(buffer_a, "event%03d_bkg.csv", iEvent);  // format string, insert int as a 3-digit number
            aux_ss << buffer_a;
            aux_ss >> output_filename;
        }
        output_file.open(output_filename);

        // reset counter
        n_particles = 0;

        // (debug)
        std::cout << "main_fct :: >> Event " << iEvent << std::endl;

        // particle loop -- print particle info
        for (int i = 0; i < pythia.event.size(); i++) {

            // (cut) on status
            if (!pythia.event[i].isFinal()) {
                continue;
            }

            // (cut) on theta
            if (pythia.event[i].theta() < 3.04 || pythia.event[i].theta() > 3.13) {
                continue;
            }

            // get mother indices
            mother1 = pythia.event[i].mother1();
            mother2 = pythia.event[i].mother2();
            // assign mother's pid
            // verify the "normal" mother case, where it is meaningful to speak of one single mother to several products
            // otherwise, set PID of mother to 0
            id_mother = (mother1 > 0 && mother2 == 0) ? pythia.event[mother1].id() : 0;

            sprintf(buffer_b, "%i,%i,%i,%.8e,%.8e,%.8e\n",                      //
                    pythia.event[i].status(), pythia.event[i].id(), id_mother,  //
                    pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz());

            // (debug)
            std::cout << "main_fct :: " << buffer_b;

            // (output)
            output_file << buffer_b;

            // update counter
            n_particles++;
        }

        output_file.close();

        if (n_particles) {
            // (debug)
            std::cout << "main_fct :: >> " << n_particles << " particles have been generated and stored in " << output_filename
                      << std::endl;
        } else {
            // (debug)
            std::cout << "main_fct :: >> " << n_particles << " particles have been generated. Running it again..." << std::endl;
            iEvent--;
        }

        // (debug) empty line
        std::cout << std::endl;

    }  // end event loop

    // print statistics
    pythia.stat();

    return 0;
}
