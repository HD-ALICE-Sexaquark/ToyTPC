// main01.cc is a part of the PYTHIA event generator.
// Copyright (C) 2023 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

/*
 USAGE: after compilation, ./main_GenCollision config_pp.cmnd <n_events> <output_dir>
*/

#include <fstream>
#include <iostream>
#include <sstream>  // stringstream, useful to convert char to str
#include <string>

#include "Pythia8/HeavyIons.h"
#include "Pythia8/Pythia.h"

#define DEBUG 1

using namespace Pythia8;

void print_help() {
    std::cout << "GenCollision :: USAGE: ./main_GenCollision [COMMANDS]" << std::endl;
    std::cout << "GenCollision ::        where [COMMANDS] could be:" << std::endl;
    std::cout << "GenCollision ::        --help                              : print this message, and then exit" << std::endl;
    std::cout << "GenCollision ::        --n <n>                             : number of events" << std::endl;
    std::cout << "GenCollision ::        --config <config_file>  (mandatory) : input configuration file" << std::endl;
    std::cout << "GenCollision ::        --output <output_dir>               : output directory" << std::endl;
}

int main(int argc, char* argv[]) {

    /* Parse command-line options */

    std::string config_file;
    std::string output_dir;
    int nEvent = 0;
    for (int i = 0; i < argc; i++) {
        if ((std::string)argv[i] == "--help") {
            print_help();
            return 0;
        } else if ((std::string)argv[i] == "--n") {
            nEvent = atoi(argv[i + 1]);
        } else if ((std::string)argv[i] == "--config") {
            config_file = argv[i + 1];
        } else if ((std::string)argv[i] == "--output") {
            output_dir = argv[i + 1];
        }
    }
    if (output_dir == "") output_dir = ".";
    if (config_file == "") {
        std::cerr << "GenCollision :: ERROR: an input configuration file must be provided." << std::endl;
        print_help();
        return 1;
    }

    std::cout << "GenCollision :: Input Options" << std::endl;
    std::cout << "GenCollision :: =============" << std::endl;
    std::cout << "GenCollision :: >> n_events    = " << nEvent << std::endl;
    std::cout << "GenCollision :: >> config_file = " << config_file << std::endl;
    std::cout << "GenCollision :: >> output_dir  = " << output_dir << std::endl;

    // declare generator
    Pythia pythia;

    // read config file
    pythia.readFile(config_file);
    if (!nEvent) nEvent = pythia.mode("Main:numberOfEvents");

    // initialize
    pythia.init();

    // prepare output file
    std::ofstream output_file;

    // string and auxiliary types to help with formatting
    std::string output_filename;
    std::stringstream aux_ss;
    char buffer_a[48];
    char buffer_b[96];

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
            output_filename = "mc.csv";
        } else {
            aux_ss.clear();
            sprintf(buffer_a, "%s/event%03d_mc.csv", output_dir.c_str(), iEvent);  // format string, insert int as a 3-digit number
            aux_ss << buffer_a;
            aux_ss >> output_filename;
        }
        output_file.open(output_filename);

        // reset counter
        n_particles = 0;

#if DEBUG
        std::cout << "GenCollision :: Event " << iEvent << std::endl;
#endif

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

            sprintf(buffer_b, "%i,%.8e,%.8e,%.8e\n",                      //
                    pythia.event[i].id(), pythia.event[i].px(), pythia.event[i].py(), pythia.event[i].pz());

#if DEBUG
            std::cout << "GenCollision :: " << buffer_b;
#endif

            // (output)
            output_file << buffer_b;

            // update counter
            n_particles++;
        }

        output_file.close();

#if DEBUG
        if (n_particles) {
            std::cout << "GenCollision :: " << n_particles << " particles have been generated." << std::endl;
            std::cout << "GenCollision :: File " << output_filename << " has been created." << std::endl;
        } else {
            std::cout << "GenCollision :: No particles have been generated. Running it again..." << std::endl;
        }
        std::cout << std::endl;
#endif

        if (!n_particles) iEvent--;
    }  // end event loop

    // print statistics
    pythia.stat();

    return 0;
}
