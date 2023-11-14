#include <fstream>
#include <iostream>
#include <map>
#include <sstream>

#include "TF1.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TString.h"
#include "TTree.h"
#include "TVector3.h"

/*
 Container of ITS hit information.
 */
struct ITSHit_tt {
    Int_t layerNb;
    Float_t x;
    Float_t y;
    Float_t z;
    Float_t edep;
};

/*
 Container of TPC hit information.
 */
struct TPCHit_tt {
    Float_t x;
    Float_t y;
    Float_t z;
    Float_t px;
    Float_t py;
    Float_t pz;
    Float_t time;
    Float_t edep;
};

/*
 Container for particle information.
 */
struct Particle_tt {
    Int_t trackID;
    Int_t PDGcode;
    Float_t x_ini;
    Float_t y_ini;
    Float_t z_ini;
    Float_t px_ini;
    Float_t py_ini;
    Float_t pz_ini;
    Int_t parentID;
    Int_t n_daughters;
    Bool_t is_primary;
    Int_t charge;
    std::vector<ITSHit_tt> its_hits;
    std::vector<TPCHit_tt> tpc_hits;
};

/*
 Container for event information.
 */
struct Event_tt {
    Int_t eventID;
    std::vector<Particle_tt> particles;
};

void ParseCSVFiles(Int_t run_n = 0) {

    // set input/output filenames
    TString input_traj_file = Form("run%03d_traj.csv", run_n);
    TString input_its_file = Form("run%03d_its.csv", run_n);
    TString input_tpc_file = Form("run%03d_tpc.csv", run_n);
    TString output_filename = Form("run%03d_ana.root", run_n);

    /*** Part 1: Read and Store Input ***/

    // init particles -- main container
    std::vector<Event_tt> Events;

    // vector and map to link between eventID, trackID and main vector indices
    std::vector<std::map<Int_t, Int_t>> map_event_index;  // map[eventID][trackID] = index
    std::map<Int_t, Int_t> map_track_index;               // key = trackID, value = index

    /** Part 1a: Trajectories **/

    // load file into memory
    std::fstream traj_file;
    traj_file.open(input_traj_file);
    if (!traj_file) {
        std::cerr << "ParseCSVFiles.C :: ERROR :: File " << input_traj_file << " not found." << std::endl;
        return;
    }

    /* (Auxiliary variables) */

    Particle_tt aux_particle;
    Event_tt aux_event;
    Int_t current_eventID;
    Int_t prev_eventID = -1;

    std::vector<std::map<Int_t, Int_t>> map_event_ndaughters;  // map[eventID][trackID] = n_daughters
    std::map<Int_t, Int_t> n_daughters;                        // key = trackID, value = n_daughters

    // read line by line ~ loop over particles
    TObjArray *token = nullptr;
    std::string line;
    while (std::getline(traj_file, line)) {

        // (protection)
        if (line == "") {
            std::cerr << "ParseCSVFiles.C :: WARNING :: Line empty, skipping..." << std::endl;
            continue;
        }

        // convert the line into a TString, then split
        token = ((TString)line).Tokenize(",");

        /*  [0] eventID     */ current_eventID = ((TString)(token->At(0)->GetName())).Atoi();

        if (current_eventID != prev_eventID) {
            // if the event value changed, we're in a different event now
            if (prev_eventID != -1) {
                // then, store past event info into main vector, clean particles vector
                aux_event.eventID = prev_eventID;
                Events.push_back(aux_event);
                aux_event.particles.clear();
                // assign map of tracks->indices to corresponding event, clean
                map_event_index.push_back(map_track_index);
                map_track_index.clear();
                // assign map of tracks->n_daughters to corresponding event, clean
                map_event_ndaughters.push_back(n_daughters);
                n_daughters.clear();
            }
            // now, we're inside another event
            prev_eventID = current_eventID;
        }

        /*  [1] trackID     */ aux_particle.trackID = ((TString)(token->At(1)->GetName())).Atoi();
        /*  [2] PDGcode     */ aux_particle.PDGcode = ((TString)(token->At(2)->GetName())).Atoi();
        /*  [3] x_ini       */ aux_particle.x_ini = ((TString)(token->At(3)->GetName())).Atof();
        /*  [4] y_ini       */ aux_particle.y_ini = ((TString)(token->At(4)->GetName())).Atof();
        /*  [5] z_ini       */ aux_particle.z_ini = ((TString)(token->At(5)->GetName())).Atof();
        /*  [6] px_ini      */ aux_particle.px_ini = ((TString)(token->At(6)->GetName())).Atof();
        /*  [7] py_ini      */ aux_particle.py_ini = ((TString)(token->At(7)->GetName())).Atof();
        /*  [8] pz_ini      */ aux_particle.pz_ini = ((TString)(token->At(8)->GetName())).Atof();
        /*  [9] parentID    */ aux_particle.parentID = ((TString)(token->At(9)->GetName())).Atoi();
        /*  [-] is_primary  */ aux_particle.is_primary = aux_particle.parentID == 0;
        /* [10] charge      */ aux_particle.charge = ((TString)(token->At(10)->GetName())).Atoi();

        aux_event.particles.push_back(aux_particle);

        // update maps
        map_track_index[aux_particle.trackID] = (Int_t)aux_event.particles.size() - 1;
        n_daughters[aux_particle.parentID]++;
    }  // end of reading file

    traj_file.close();

    /* Posteriori variable: N Daughters */

    for (Event_tt &evt : Events) {
        for (Particle_tt &part : evt.particles) {
            part.n_daughters = map_event_ndaughters[evt.eventID][part.trackID];
        }
    }

    /** Part 1b: ITS Hits **/

    /* (Auxiliary variables) */

    ITSHit_tt aux_its_hit;
    Int_t aux_track_id;
    Int_t aux_index;

    // load file into memory
    std::fstream its_file;
    its_file.open(input_its_file);
    if (!its_file) {
        std::cerr << "ParseCSVFiles.C :: ERROR :: File " << input_its_file << " not found." << std::endl;
        return;
    }

    // read line by line ~ loop over hits
    line.clear();
    while (std::getline(its_file, line)) {

        // (protection)
        if (line == "") {
            std::cerr << "ParseCSVFiles.C :: WARNING :: Line empty, skipping..." << std::endl;
            continue;
        }

        // convert the line into a TString, then split
        token = ((TString)line).Tokenize(",");

        /* [0] eventID  */ current_eventID = ((TString)(token->At(0)->GetName())).Atoi();
        /* [1] trackID  */ aux_track_id = ((TString)(token->At(1)->GetName())).Atoi();
        /* [2] layerNb  */ aux_its_hit.layerNb = ((TString)(token->At(2)->GetName())).Atoi();
        /* [3] x        */ aux_its_hit.x = ((TString)(token->At(3)->GetName())).Atof();
        /* [4] y        */ aux_its_hit.y = ((TString)(token->At(4)->GetName())).Atof();
        /* [5] z        */ aux_its_hit.z = ((TString)(token->At(5)->GetName())).Atof();
        /* [6] edep     */ aux_its_hit.edep = ((TString)(token->At(6)->GetName())).Atof();
        /* [7] process  */  // not used

        // get index from the track ID
        aux_index = map_event_index[current_eventID][aux_track_id];

        // store current hit info on corresponding index
        Events.at(current_eventID).particles.at(aux_index).its_hits.push_back(aux_its_hit);
    }  // end of reading file

    its_file.close();

    /** Part 1c: TPC Hits **/

    /* (Auxiliary variables) */

    TPCHit_tt aux_tpc_hit;

    // load file into memory
    std::fstream tpc_file;
    tpc_file.open(input_tpc_file);
    if (!tpc_file) {
        std::cerr << "ParseCSVFiles.C :: ERROR :: File " << input_tpc_file << " not found." << std::endl;
        return;
    }

    // read line by line ~ loop over hits
    line.clear();
    while (std::getline(tpc_file, line)) {

        // (protection)
        if (line == "") {
            std::cerr << "ParseCSVFiles.C :: WARNING :: Line empty, skipping..." << std::endl;
            continue;
        }

        // convert the line into a TString, then split
        token = ((TString)line).Tokenize(",");

        /*  [0] eventID  */ current_eventID = ((TString)(token->At(0)->GetName())).Atoi();
        /*  [1] trackID  */ aux_track_id = ((TString)(token->At(1)->GetName())).Atoi();
        /*  [2] x        */ aux_tpc_hit.x = ((TString)(token->At(2)->GetName())).Atof();
        /*  [3] y        */ aux_tpc_hit.y = ((TString)(token->At(3)->GetName())).Atof();
        /*  [4] z        */ aux_tpc_hit.z = ((TString)(token->At(4)->GetName())).Atof();
        /*  [5] px       */ aux_tpc_hit.px = ((TString)(token->At(5)->GetName())).Atof();
        /*  [6] py       */ aux_tpc_hit.py = ((TString)(token->At(6)->GetName())).Atof();
        /*  [7] pz       */ aux_tpc_hit.pz = ((TString)(token->At(7)->GetName())).Atof();
        /*  [8] time     */ aux_tpc_hit.time = ((TString)(token->At(8)->GetName())).Atof();
        /*  [9] edep     */ aux_tpc_hit.edep = ((TString)(token->At(9)->GetName())).Atof();
        /* [10] process  */  // not used

        // get index from the track ID
        aux_index = map_event_index[current_eventID][aux_track_id];

        // store current hit info on corresponding index
        Events.at(current_eventID).particles.at(aux_index).tpc_hits.push_back(aux_tpc_hit);
    }  // end of reading file

    tpc_file.close();

    /*** Debug: Check content of Events vector ***/

    for (Event_tt &evt : Events) {

        std::cout << "ParseCSVFiles.C :: Printing Event #" << evt.eventID << std::endl;
        std::cout << "ParseCSVFiles.C :: (n particles = " << (Int_t)evt.particles.size() << ")" << std::endl;
        std::cout << "ParseCSVFiles.C :: " << std::endl;

        for (Particle_tt &part : evt.particles) {

            std::cout << "ParseCSVFiles.C :: part " << part.trackID << std::endl;
            std::cout << "ParseCSVFiles.C :: >> pdg,is_primary " << part.PDGcode << ", " << part.is_primary << std::endl;
            std::cout << "ParseCSVFiles.C :: >> parent,n_daughters " << part.parentID << ", " << part.n_daughters << std::endl;
            std::cout << "ParseCSVFiles.C :: >> x,y,z " << part.x_ini << ", " << part.y_ini << ", " << part.z_ini << ", " << std::endl;
            std::cout << "ParseCSVFiles.C :: >> px,py,pz " << part.px_ini << ", " << part.py_ini << ", " << part.pz_ini << ", "
                      << std::endl;
            std::cout << "ParseCSVFiles.C :: >> n_its_hits " << (Int_t)part.its_hits.size() << std::endl;
            std::cout << "ParseCSVFiles.C :: >> n_tpc_hits " << (Int_t)part.tpc_hits.size() << std::endl;

            /*
            for (ITSHit_tt &its_hit : part.its_hits) {
                std::cout << "ParseCSVFiles.C ::    >> layer " << its_hit.layerNb << std::endl;
                std::cout << "ParseCSVFiles.C ::    >> edep " << its_hit.edep << std::endl;
                std::cout << "ParseCSVFiles.C ::    >> x,y,z " << its_hit.x << ", " << its_hit.y << ", " << its_hit.z << ", " << std::endl;
            }
            for (TPCHit_tt &tpc_hit : part.tpc_hits) {
                std::cout << "ParseCSVFiles.C ::    >> edep " << tpc_hit.edep << std::endl;
                std::cout << "ParseCSVFiles.C ::    >> x,y,z " << tpc_hit.x << ", " << tpc_hit.y << ", " << tpc_hit.z << ", " << std::endl;
                std::cout << "ParseCSVFiles.C ::    >> px,py,pz " << tpc_hit.x << ", " << tpc_hit.y << ", " << tpc_hit.z << ", "
                          << std::endl;
            }
            */
        }
        std::cout << "ParseCSVFiles.C :: " << std::endl;
    }

    /*** Part 2: Trees ***/

}  // end of macro
