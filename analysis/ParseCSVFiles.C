#include <fstream>
#include <iostream>
#include <map>
#include <sstream>

#include "TDatabasePDG.h"
#include "TF1.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TRandom3.h"
#include "TString.h"
#include "TTree.h"
#include "TVector3.h"

/*
 Container of ITS hit information.
 */
struct ITSHit_t {
    Int_t layerNb;
    Float_t x;
    Float_t y;
    Float_t z;
    Float_t dedx;
};

/*
 Container of TPC hit information.
 */
struct TPCHit_t {
    Float_t x;
    Float_t y;
    Float_t z;
    Float_t px;
    Float_t py;
    Float_t pz;
    Float_t time;
    Float_t dedx;
};

/*
 Container for particle information.
 */
struct Particle_t {
    Int_t eventID;
    Int_t trackID;
    Int_t PDGcode;
    Float_t x_ini;
    Float_t y_ini;
    Float_t z_ini;
    Float_t px_ini;
    Float_t py_ini;
    Float_t pz_ini;
    Int_t parentID;
    Bool_t is_primary;
    Int_t charge;
    TString process;
    std::vector<ITSHit_t> its_hits;
    std::vector<TPCHit_t> tpc_hits;
};

void ParseCSVFiles(Int_t event_n = 0) {

    // set input/output filenames
    TString input_traj_file = Form("event%i_traj.csv", event_n);  // %03i
    TString input_its_file = Form("event%i_its.csv", event_n);    // %03i
    TString input_tpc_file = Form("event%i_tpc.csv", event_n);    // %03i
    TString output_filename = Form("event%03i_ana.root", event_n);

    /*** Part 1: Read and Store Input ***/

    // init particles -- main container
    std::vector<Particle_t> Particles;

    // map to link between track ID and main vector indices
    std::map<Int_t, Int_t> map_track_index;  // map[trackID] = index

    // init PDG database object
    TDatabasePDG pdg;
    const Float_t kMassPion = pdg.GetParticle(211)->Mass();
    const Float_t kMassKaon = pdg.GetParticle(321)->Mass();
    const Float_t kMassProton = pdg.GetParticle(2212)->Mass();

    /** Part 1a: Trajectories **/

    // load file into memory
    std::fstream traj_file;
    traj_file.open(input_traj_file);
    if (!traj_file) {
        std::cerr << "ParseCSVFiles.C :: ERROR :: File " << input_traj_file << " not found." << std::endl;
        return;
    }

    // auxiliary object
    Particle_t aux_particle;

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

        /*  [0] eventID     */ aux_particle.eventID = ((TString)(token->At(0)->GetName())).Atoi();
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
        /* [10] charge      */ aux_particle.charge = ((TString)(token->At(0)->GetName())).Atoi();

        // update map
        map_track_index[aux_particle.trackID] = (Int_t)Particles.size();

        // store current particle info into main vector
        Particles.push_back(aux_particle);
    }  // end of reading file

    traj_file.close();

    /* (Auxiliary variables) */

    Int_t aux_track_id;
    Int_t aux_index;
    TString aux_process;

    /** Part 1b: ITS Hits **/

    // load file into memory
    std::fstream its_file;
    its_file.open(input_its_file);
    if (!its_file) {
        std::cerr << "ParseCSVFiles.C :: ERROR :: File " << input_its_file << " not found." << std::endl;
        return;
    }

    // auxiliary object
    ITSHit_t aux_its_hit;

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

        /* [0] eventID  */  // not used
        /* [1] trackID  */ aux_track_id = ((TString)(token->At(1)->GetName())).Atoi();
        /* [2] layerNb  */ aux_its_hit.layerNb = ((TString)(token->At(2)->GetName())).Atoi();
        /* [3] x        */ aux_its_hit.x = ((TString)(token->At(3)->GetName())).Atof();
        /* [4] y        */ aux_its_hit.y = ((TString)(token->At(4)->GetName())).Atof();
        /* [5] z        */ aux_its_hit.z = ((TString)(token->At(5)->GetName())).Atof();
        /* [6] Edep     */ aux_its_hit.dedx = ((TString)(token->At(6)->GetName())).Atof();
        /* [7] process  */ aux_process = (TString)(token->At(7)->GetName());

        // get index from the track ID
        aux_index = map_track_index[aux_track_id];

        // if not filled, fill particle's process
        if (!Particles.at(aux_index).process) {
            Particles.at(aux_index).process = aux_process;
        }

        // store current hit info on corresponding index
        Particles.at(aux_index).its_hits.push_back(aux_its_hit);
    }  // end of reading file

    its_file.close();

    /** Part 1c: TPC Hits **/

    // load file into memory
    std::fstream tpc_file;
    tpc_file.open(input_tpc_file);
    if (!tpc_file) {
        std::cerr << "ParseCSVFiles.C :: ERROR :: File " << input_tpc_file << " not found." << std::endl;
        return;
    }

    // auxiliary object
    TPCHit_t aux_tpc_hit;

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

        /*  [0] eventID  */  // not used
        /*  [1] trackID  */ aux_track_id = ((TString)(token->At(1)->GetName())).Atoi();
        /*  [2] x        */ aux_tpc_hit.x = ((TString)(token->At(2)->GetName())).Atof();
        /*  [3] y        */ aux_tpc_hit.y = ((TString)(token->At(3)->GetName())).Atof();
        /*  [4] z        */ aux_tpc_hit.z = ((TString)(token->At(4)->GetName())).Atof();
        /*  [5] px       */ aux_tpc_hit.px = ((TString)(token->At(5)->GetName())).Atof();
        /*  [6] py       */ aux_tpc_hit.py = ((TString)(token->At(6)->GetName())).Atof();
        /*  [7] pz       */ aux_tpc_hit.pz = ((TString)(token->At(7)->GetName())).Atof();
        /*  [8] time     */ aux_tpc_hit.time = ((TString)(token->At(8)->GetName())).Atof();
        /*  [9] dedx     */ aux_tpc_hit.dedx = ((TString)(token->At(9)->GetName())).Atof();
        /* [10] process  */ aux_process = (TString)(token->At(10)->GetName());

        // get index from the track ID
        aux_index = map_track_index[aux_track_id];

        // if not filled, fill particle's process
        if (!Particles.at(aux_index).process) {
            Particles.at(aux_index).process = aux_process;
        }

        // store current hit info on corresponding index
        Particles.at(aux_index).tpc_hits.push_back(aux_tpc_hit);
    }  // end of reading file

    tpc_file.close();

    /*** Debug: Check content of Particles vector ***/

    std::cout << "ParseCSVFiles.C :: >> N_Particles = " << (Int_t)Particles.size() << std::endl;
    std::cout << "ParseCSVFiles.C :: " << std::endl;

    for (Particle_t part : Particles) {
        std::cout << "ParseCSVFiles.C :: part " << part.trackID << std::endl;
        std::cout << "ParseCSVFiles.C :: >> pdg,is_primary,process " << part.PDGcode << ", " << part.is_primary << ", " << part.process
                  << std::endl;
        std::cout << "ParseCSVFiles.C :: >> x,y,z " << part.x_ini << ", " << part.y_ini << ", " << part.z_ini << ", " << std::endl;
        std::cout << "ParseCSVFiles.C :: >> px,py,pz " << part.px_ini << ", " << part.py_ini << ", " << part.pz_ini << ", " << std::endl;
        std::cout << "ParseCSVFiles.C :: >> n_its_hits " << (Int_t)part.its_hits.size() << std::endl;
        for (ITSHit_t &its_hit : part.its_hits) {
            std::cout << "ParseCSVFiles.C ::    >> layer " << its_hit.layerNb << std::endl;
            std::cout << "ParseCSVFiles.C ::    >> edep " << its_hit.dedx << std::endl;
            std::cout << "ParseCSVFiles.C ::    >> x,y,z " << its_hit.x << ", " << its_hit.y << ", " << its_hit.z << ", " << std::endl;
        }
        std::cout << "ParseCSVFiles.C :: >> n_tpc_hits " << (Int_t)part.tpc_hits.size() << std::endl;
        for (TPCHit_t &tpc_hit : part.tpc_hits) {
            std::cout << "ParseCSVFiles.C ::    >> edep " << tpc_hit.dedx << std::endl;
            std::cout << "ParseCSVFiles.C ::    >> x,y,z " << tpc_hit.x << ", " << tpc_hit.y << ", " << tpc_hit.z << ", " << std::endl;
            std::cout << "ParseCSVFiles.C ::    >> px,py,pz " << tpc_hit.x << ", " << tpc_hit.y << ", " << tpc_hit.z << ", " << std::endl;
        }
    }

    /*** Part 2: Histograms! ***/

}  // end of macro
