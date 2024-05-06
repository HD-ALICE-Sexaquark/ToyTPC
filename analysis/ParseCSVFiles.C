#include <fstream>
#include <iostream>
#include <map>
#include <sstream>

#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TString.h"
#include "TTree.h"
#include "TVector3.h"

#include "EnergyLoss.h"
#include "HelixFit.h"
#include "LeastSquaresCircleFit.h"

/*** STRUCTURES ***/

/*
 Container of ITS hit information.
*/
struct ITSHit_tt {
    Int_t layerNb;
    Float_t x;
    Float_t y;
    Float_t z;
    Float_t px;
    Float_t py;
    Float_t pz;
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
    Float_t edep;
};

/*
 Container for particle information.
*/
struct Particle_tt {
    // true info
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
    // ITS
    std::vector<ITSHit_tt> its_hits;
    Float_t its_pt_rec;
    Float_t its_pz_rec;
    Float_t its_p_rec;
    Float_t its_charge_rec;
    Float_t its_chi2_circle;
    Float_t its_chi2_helix;
    // -- first hit info, for debugging
    Float_t its_fh_px;
    Float_t its_fh_py;
    Float_t its_fh_pt;
    Float_t its_fh_pz;
    Float_t its_fh_p;
    // -- last hit info
    Float_t its_lh_px;
    Float_t its_lh_py;
    Float_t its_lh_pt;
    Float_t its_lh_pz;
    Float_t its_lh_p;
    // TPC
    std::vector<TPCHit_tt> tpc_hits;
    Float_t tpc_pt_rec;
    Float_t tpc_pz_rec;
    Float_t tpc_p_rec;
    Float_t tpc_charge_rec;
    Float_t tpc_chi2_circle;
    Float_t tpc_chi2_helix;
    Float_t tpc_dx;  // average distance between two hits in the TPC
    Float_t tpc_dE_dx;
    // -- first hit info
    Float_t tpc_fh_px;
    Float_t tpc_fh_py;
    Float_t tpc_fh_pt;
    Float_t tpc_fh_pz;
    Float_t tpc_fh_p;
};

/*
 Container for event information.
*/
struct Event_tt {
    Int_t eventID;
    std::vector<Particle_tt> particles;
};

/*** MAIN ***/

void ParseCSVFiles(Int_t job_n = 0, Int_t run_n = 0) {

    // set input/output filenames -- hardcoded (not anymore :))

    TString simdir = "/misc/alidata130/alice_u/kleine/software/ToyTPC/output/";

    TString input_traj_file = simdir + Form("%02d/", job_n) + Form("run%02d/", run_n) + Form("run%02d_traj.csv", run_n);
    TString input_its_file = simdir + Form("%02d/", job_n) + Form("run%02d/", run_n) + Form("run%02d_its.csv", run_n);
    TString input_tpc_file = simdir + Form("%02d/", job_n) + Form("run%02d/", run_n) + Form("run%02d_tpc.csv", run_n);

    TString output_filename = Form("run%02d_ana.root", run_n);

    // temporary test snippet
    /*
    void ParseCSVFiles(Int_t run_n = 0) {

        TString input_traj_file = Form("run%02d_traj_itest.csv", run_n);
        TString input_its_file = Form("run%02d_its_itest.csv", run_n);
        TString input_tpc_file = Form("run%02d_tpc_itest.csv", run_n);
        TString output_filename = Form("run%02d_ana_itest.root", run_n);
    */

    // conversion factors
    const Double_t MeVToGeV = 1E-3;
    const Double_t GeVToMeV = 1E3;

    const Double_t cmTom = 1E-2;
    const Double_t mTocm = 1E2;

    /************************************/
    /*** Part 1: Read and Store Input ***/

    // init events -- main container
    std::vector<Event_tt> Events;

    // vector and map to link between eventID, trackID and main vector indices
    std::vector<std::map<Int_t, Int_t>> part_index;       // [eventID][trackID]
    std::vector<std::map<Int_t, Int_t>> part_ndaughters;  // [eventID][trackID]

    /** Part 1a: Trajectories **/

    // load file into memory
    std::fstream traj_file;
    traj_file.open(input_traj_file);
    if (!traj_file) {
        std::cerr << "ParseCSVFiles.C :: ERROR :: File " << input_traj_file << " not found." << std::endl;
        return;
    }

    Particle_tt aux_particle;
    Int_t current_eventID;
    Int_t prev_eventID = -1;

    // read line by line ~ loop over particles
    TObjArray *token = nullptr;
    std::string line;
    while (std::getline(traj_file, line)) {

        // (protection)
        if (line == "") {
            std::cerr << "ParseCSVFiles.C :: WARNING :: Skipping empty line..." << std::endl;
            continue;
        }

        // convert the line into a TString, then split
        token = ((TString)line).Tokenize(",");

        current_eventID = ((TString)(token->At(0)->GetName())).Atoi();  // [0] eventID

        if (current_eventID != prev_eventID) {
            // resize vectors, because eventID starts at 0
            Events.resize(current_eventID + 1);
            part_index.resize(current_eventID + 1);
            part_ndaughters.resize(current_eventID + 1);
            Events.at(current_eventID).eventID = current_eventID;
            prev_eventID = current_eventID;  // we're in a different event now
        }

        aux_particle.trackID = ((TString)(token->At(1)->GetName())).Atoi();            //     [1] trackID
        aux_particle.PDGcode = ((TString)(token->At(2)->GetName())).Atoi();            //     [2] PDGcode
        aux_particle.x_ini = ((TString)(token->At(3)->GetName())).Atof();              //     [3] x_ini
        aux_particle.y_ini = ((TString)(token->At(4)->GetName())).Atof();              //     [4] y_ini
        aux_particle.z_ini = ((TString)(token->At(5)->GetName())).Atof();              //     [5] z_ini
        aux_particle.px_ini = ((TString)(token->At(6)->GetName())).Atof() * MeVToGeV;  //     [6] px_ini
        aux_particle.py_ini = ((TString)(token->At(7)->GetName())).Atof() * MeVToGeV;  //     [7] py_ini
        aux_particle.pz_ini = ((TString)(token->At(8)->GetName())).Atof() * MeVToGeV;  //     [8] pz_ini
        aux_particle.parentID = ((TString)(token->At(9)->GetName())).Atoi();           //     [9] parentID
        aux_particle.is_primary = aux_particle.parentID == 0;                          //     [-] is_primary
        aux_particle.charge = ((TString)(token->At(10)->GetName())).Atoi();            //    [10] charge

        Events.at(current_eventID).particles.push_back(aux_particle);

        // update maps
        part_index[current_eventID][aux_particle.trackID] = (Int_t)Events.at(current_eventID).particles.size() - 1;
        part_ndaughters[current_eventID][aux_particle.parentID]++;
    }  // end of reading file

    traj_file.close();

    /* Posteriori variable: N Daughters */
    /* (now that all particles are loaded) */

    for (Event_tt &evt : Events) {
        for (Particle_tt &part : evt.particles) {
            part.n_daughters = part_ndaughters[evt.eventID][part.trackID];
        }
    }

    /** Part 1b: ITS Hits **/

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
            std::cerr << "ParseCSVFiles.C :: WARNING :: Skipping empty line..." << std::endl;
            continue;
        }

        // convert the line into a TString, then split
        token = ((TString)line).Tokenize(",");

        current_eventID = ((TString)(token->At(0)->GetName())).Atoi();            // [0] eventID
        aux_track_id = ((TString)(token->At(1)->GetName())).Atoi();               // [1] trackID
        aux_its_hit.layerNb = ((TString)(token->At(2)->GetName())).Atoi();        // [2] layerNb
        aux_its_hit.x = ((TString)(token->At(3)->GetName())).Atof();              // [3] x
        aux_its_hit.y = ((TString)(token->At(4)->GetName())).Atof();              // [4] y
        aux_its_hit.z = ((TString)(token->At(5)->GetName())).Atof();              // [5] z
        aux_its_hit.px = ((TString)(token->At(6)->GetName())).Atof() * MeVToGeV;  // [6] px
        aux_its_hit.py = ((TString)(token->At(7)->GetName())).Atof() * MeVToGeV;  // [7] py
        aux_its_hit.pz = ((TString)(token->At(8)->GetName())).Atof() * MeVToGeV;  // [8] pz

        // get index from the track ID
        aux_index = part_index[current_eventID][aux_track_id];

        // store current hit info on corresponding index
        Events.at(current_eventID).particles.at(aux_index).its_hits.push_back(aux_its_hit);
    }  // end of reading file

    its_file.close();

    /** Part 1c: TPC Hits **/

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
            std::cerr << "ParseCSVFiles.C :: WARNING :: Skipping empty line..." << std::endl;
            continue;
        }

        // convert the line into a TString, then split
        token = ((TString)line).Tokenize(",");

        current_eventID = ((TString)(token->At(0)->GetName())).Atoi();            //   [0] eventID
        aux_track_id = ((TString)(token->At(1)->GetName())).Atoi();               //   [1] trackID
        aux_tpc_hit.x = ((TString)(token->At(2)->GetName())).Atof();              //   [2] x
        aux_tpc_hit.y = ((TString)(token->At(3)->GetName())).Atof();              //   [3] y
        aux_tpc_hit.z = ((TString)(token->At(4)->GetName())).Atof();              //   [4] z
        aux_tpc_hit.px = ((TString)(token->At(5)->GetName())).Atof() * MeVToGeV;  //   [5] px
        aux_tpc_hit.py = ((TString)(token->At(6)->GetName())).Atof() * MeVToGeV;  //   [6] py
        aux_tpc_hit.pz = ((TString)(token->At(7)->GetName())).Atof() * MeVToGeV;  //   [7] pz
        aux_tpc_hit.edep = ((TString)(token->At(9)->GetName())).Atof();           //   [9] edep

        // get index from the track ID
        aux_index = part_index[current_eventID][aux_track_id];

        // store current hit info on corresponding index
        Events.at(current_eventID).particles.at(aux_index).tpc_hits.push_back(aux_tpc_hit);
    }  // end of reading file

    tpc_file.close();

    /** Tracking **/
    /* (now that all hits are loaded) */

    Double_t p_ini;

    Double_t x_c, y_c, radius, chi2_circle;
    Bool_t circlefit_state;

    Double_t angle, charge, chi2_helix;
    Int_t direction;
    Bool_t helixfit_state;

    Int_t n_its_hits, n_tpc_hits;

    for (Event_tt &evt : Events) {
        for (Particle_tt &part : evt.particles) {

            /* Protection */

            p_ini = TMath::Sqrt(part.px_ini * part.px_ini + part.py_ini * part.py_ini + part.pz_ini * part.pz_ini);
            if (p_ini < 0.01) continue;

            /** ITS Tracking **/

            n_its_hits = (Int_t)part.its_hits.size();

            if (n_its_hits > 3) {

                Double_t its_x[n_its_hits];
                Double_t its_y[n_its_hits];
                Double_t its_z[n_its_hits];

                for (Int_t hit_idx = 0; hit_idx < n_its_hits; hit_idx++) {
                    its_x[hit_idx] = part.its_hits[hit_idx].x;
                    its_y[hit_idx] = part.its_hits[hit_idx].y;
                    its_z[hit_idx] = part.its_hits[hit_idx].z;
                }

                /* Fit ITS hits to a circle */

                circlefit_state = LeastSquaresCircleFit(n_its_hits, its_x, its_y, x_c, y_c, radius, chi2_circle);

                /* Fit ITS hits to helix */

                helixfit_state = HelixFit(n_its_hits, its_x, its_y, its_z, x_c, y_c, radius, angle, charge, direction, chi2_helix);

                /* Store results */

                part.its_fh_px = part.its_hits[0].px;
                part.its_fh_py = part.its_hits[0].py;
                part.its_fh_pt = TMath::Sqrt(part.its_fh_px * part.its_fh_px + part.its_fh_py * part.its_fh_py);
                part.its_fh_pz = part.its_hits[0].pz;
                part.its_fh_p = TMath::Sqrt(part.its_fh_pt * part.its_fh_pt + part.its_fh_pz * part.its_fh_pz);

                part.its_lh_px = part.its_hits[n_its_hits - 1].px;
                part.its_lh_py = part.its_hits[n_its_hits - 1].py;
                part.its_lh_pt = TMath::Sqrt(part.its_lh_px * part.its_lh_px + part.its_lh_py * part.its_lh_py);
                part.its_lh_pz = part.its_hits[n_its_hits - 1].pz;
                part.its_lh_p = TMath::Sqrt(part.its_lh_pt * part.its_lh_pt + part.its_lh_pz * part.its_lh_pz);

                part.its_pt_rec = 0.3 * 0.2 * radius * cmTom;
                part.its_pz_rec = direction * part.its_pt_rec * TMath::Abs(TMath::Tan(angle)) * GeVToMeV;
                part.its_p_rec = TMath::Sqrt(part.its_pt_rec * part.its_pt_rec + part.its_pz_rec * part.its_pz_rec);
                part.its_charge_rec = charge;
                part.its_chi2_circle = chi2_circle;
                part.its_chi2_helix = chi2_helix;
            }

            /** TPC Tracking **/

            n_tpc_hits = (Int_t)part.tpc_hits.size();

            if (n_tpc_hits >= 100) {

                Double_t tpc_x[n_tpc_hits];
                Double_t tpc_y[n_tpc_hits];
                Double_t tpc_z[n_tpc_hits];
                Double_t tpc_edep[n_tpc_hits];

                for (Int_t hit_idx = 0; hit_idx < n_tpc_hits; hit_idx++) {
                    tpc_x[hit_idx] = part.tpc_hits[hit_idx].x;
                    tpc_y[hit_idx] = part.tpc_hits[hit_idx].y;
                    tpc_z[hit_idx] = part.tpc_hits[hit_idx].z;
                    tpc_edep[hit_idx] = part.tpc_hits[hit_idx].edep;
                }

                /* Fit TPC hits to a circle */

                circlefit_state = LeastSquaresCircleFit(n_tpc_hits, tpc_x, tpc_y, x_c, y_c, radius, chi2_circle);

                /* Fit TPC hits to helix */

                helixfit_state = HelixFit(n_tpc_hits, tpc_x, tpc_y, tpc_z, x_c, y_c, radius, angle, charge, direction, chi2_helix);

                /* Compute the differential energy loss of the particle */

                Double_t dx, dE_dx;
                EnergyLoss(n_tpc_hits, tpc_edep, tpc_x, tpc_y, tpc_z, dE_dx, dx);

                /* Store results */

                part.tpc_fh_px = part.tpc_hits[0].px;
                part.tpc_fh_py = part.tpc_hits[0].py;
                part.tpc_fh_pt = TMath::Sqrt(part.tpc_fh_px * part.tpc_fh_px + part.tpc_fh_py * part.tpc_fh_py);
                part.tpc_fh_pz = part.tpc_hits[0].pz;
                part.tpc_fh_p = TMath::Sqrt(part.tpc_fh_pt * part.tpc_fh_pt + part.tpc_fh_pz * part.tpc_fh_pz);

                part.tpc_pt_rec = 0.3 * 0.2 * radius * cmTom;
                part.tpc_pz_rec = direction * part.tpc_pt_rec * TMath::Abs(TMath::Tan(angle)) * GeVToMeV;
                part.tpc_p_rec = TMath::Sqrt(part.tpc_pt_rec * part.tpc_pt_rec + part.tpc_pz_rec * part.tpc_pz_rec);
                part.tpc_charge_rec = charge;
                part.tpc_chi2_circle = chi2_circle;
                part.tpc_chi2_helix = chi2_helix;

                part.tpc_dx = dx;
                part.tpc_dE_dx = dE_dx;
            }
        }
    }

    /*********************/
    /*** Part 2: Trees ***/

    // prepare output
    TFile *output_file = new TFile(output_filename, "RECREATE");
    TTree *output_tree = new TTree("Events", "Events");

    // set aux objects
    Int_t aux_eventID;
    std::vector<Int_t> aux_trackID;
    std::vector<Int_t> aux_PDGcode;
    std::vector<Float_t> aux_Px_ini;
    std::vector<Float_t> aux_Py_ini;
    std::vector<Float_t> aux_Pz_ini;
    std::vector<Float_t> aux_Pt_ini;
    std::vector<Float_t> aux_P_ini;
    std::vector<Int_t> aux_parentID;
    std::vector<Int_t> aux_nDaughters;
    std::vector<Bool_t> aux_isPrimary;
    std::vector<Float_t> aux_charge;
    std::vector<Int_t> aux_nITShits;
    std::vector<Float_t> aux_ITS_Pt_rec;
    std::vector<Float_t> aux_ITS_Pz_rec;
    std::vector<Float_t> aux_ITS_P_rec;
    std::vector<Float_t> aux_ITS_charge_rec;
    std::vector<Float_t> aux_ITS_chi2_circle;
    std::vector<Float_t> aux_ITS_chi2_helix;
    std::vector<Float_t> aux_ITS_FH_Px;
    std::vector<Float_t> aux_ITS_FH_Py;
    std::vector<Float_t> aux_ITS_FH_Pt;
    std::vector<Float_t> aux_ITS_FH_Pz;
    std::vector<Float_t> aux_ITS_FH_P;
    std::vector<Float_t> aux_ITS_LH_Px;
    std::vector<Float_t> aux_ITS_LH_Py;
    std::vector<Float_t> aux_ITS_LH_Pt;
    std::vector<Float_t> aux_ITS_LH_Pz;
    std::vector<Float_t> aux_ITS_LH_P;
    std::vector<Int_t> aux_nTPChits;
    std::vector<Float_t> aux_TPC_Pt_rec;
    std::vector<Float_t> aux_TPC_Pz_rec;
    std::vector<Float_t> aux_TPC_P_rec;
    std::vector<Float_t> aux_TPC_charge_rec;
    std::vector<Float_t> aux_TPC_chi2_circle;
    std::vector<Float_t> aux_TPC_chi2_helix;
    std::vector<Float_t> aux_TPC_dx;
    std::vector<Float_t> aux_TPC_dEdx;
    std::vector<Float_t> aux_TPC_FH_Px;
    std::vector<Float_t> aux_TPC_FH_Py;
    std::vector<Float_t> aux_TPC_FH_Pt;
    std::vector<Float_t> aux_TPC_FH_Pz;
    std::vector<Float_t> aux_TPC_FH_P;

    // set tree branches
    output_tree->Branch("eventID", &aux_eventID);
    output_tree->Branch("trackID", &aux_trackID);
    output_tree->Branch("PDGcode", &aux_PDGcode);
    output_tree->Branch("Px_ini", &aux_Px_ini);
    output_tree->Branch("Py_ini", &aux_Py_ini);
    output_tree->Branch("Pz_ini", &aux_Pz_ini);
    output_tree->Branch("Pt_ini", &aux_Pt_ini);
    output_tree->Branch("P_ini", &aux_P_ini);
    output_tree->Branch("parentID", &aux_parentID);
    output_tree->Branch("nDaughters", &aux_nDaughters);
    output_tree->Branch("isPrimary", &aux_isPrimary);
    output_tree->Branch("charge", &aux_charge);
    output_tree->Branch("nITShits", &aux_nITShits);
    output_tree->Branch("ITS_Pt_rec", &aux_ITS_Pt_rec);
    output_tree->Branch("ITS_Pz_rec", &aux_ITS_Pz_rec);
    output_tree->Branch("ITS_P_rec", &aux_ITS_P_rec);
    output_tree->Branch("ITS_charge_rec", &aux_ITS_charge_rec);
    output_tree->Branch("ITS_chi2_circle", &aux_ITS_chi2_circle);
    output_tree->Branch("ITS_chi2_helix", &aux_ITS_chi2_helix);
    output_tree->Branch("ITS_FH_Px", &aux_ITS_FH_Px);
    output_tree->Branch("ITS_FH_Py", &aux_ITS_FH_Py);
    output_tree->Branch("ITS_FH_Pt", &aux_ITS_FH_Pt);
    output_tree->Branch("ITS_FH_Pz", &aux_ITS_FH_Pz);
    output_tree->Branch("ITS_FH_P", &aux_ITS_FH_P);
    output_tree->Branch("ITS_LH_Px", &aux_ITS_LH_Px);
    output_tree->Branch("ITS_LH_Py", &aux_ITS_LH_Py);
    output_tree->Branch("ITS_LH_Pt", &aux_ITS_LH_Pt);
    output_tree->Branch("ITS_LH_Pz", &aux_ITS_LH_Pz);
    output_tree->Branch("ITS_LH_P", &aux_ITS_LH_P);
    output_tree->Branch("nTPChits", &aux_nTPChits);
    output_tree->Branch("TPC_Pt_rec", &aux_TPC_Pt_rec);
    output_tree->Branch("TPC_Pz_rec", &aux_TPC_Pz_rec);
    output_tree->Branch("TPC_P_rec", &aux_TPC_P_rec);
    output_tree->Branch("TPC_charge_rec", &aux_TPC_charge_rec);
    output_tree->Branch("TPC_chi2_circle", &aux_TPC_chi2_circle);
    output_tree->Branch("TPC_chi2_helix", &aux_TPC_chi2_helix);
    output_tree->Branch("TPC_dx", &aux_TPC_dx);
    output_tree->Branch("TPC_dEdx", &aux_TPC_dEdx);
    output_tree->Branch("TPC_FH_Px", &aux_TPC_FH_Px);
    output_tree->Branch("TPC_FH_Py", &aux_TPC_FH_Py);
    output_tree->Branch("TPC_FH_Pt", &aux_TPC_FH_Pt);
    output_tree->Branch("TPC_FH_Pz", &aux_TPC_FH_Pz);
    output_tree->Branch("TPC_FH_P", &aux_TPC_FH_P);

    for (Event_tt &evt : Events) {

        aux_eventID = evt.eventID;

        for (Particle_tt &part : evt.particles) {

            /* Cuts */

            // if ((Int_t)part.tpc_hits.size() < 100) continue;

            p_ini = TMath::Sqrt(part.px_ini * part.px_ini + part.py_ini * part.py_ini + part.pz_ini * part.pz_ini);
            if (p_ini < 0.01) continue;

            // if (part.tpc_chi2_helix > 100) continue;

            // if (part.tpc_chi2_circle > 5) continue;

            aux_trackID.push_back(part.trackID);
            aux_PDGcode.push_back(part.PDGcode);
            aux_Px_ini.push_back(part.px_ini);
            aux_Py_ini.push_back(part.py_ini);
            aux_Pz_ini.push_back(part.pz_ini);
            aux_Pt_ini.push_back(TMath::Sqrt(part.px_ini * part.px_ini + part.py_ini * part.py_ini));
            aux_P_ini.push_back(p_ini);
            aux_parentID.push_back(part.parentID);
            aux_nDaughters.push_back(part.n_daughters);
            aux_isPrimary.push_back(part.is_primary);
            aux_charge.push_back(part.charge);
            aux_nITShits.push_back((Int_t)part.its_hits.size());
            aux_ITS_Pt_rec.push_back(part.its_pt_rec);
            aux_ITS_Pz_rec.push_back(part.its_pz_rec);
            aux_ITS_P_rec.push_back(part.its_p_rec);
            aux_ITS_charge_rec.push_back(part.its_charge_rec);
            aux_ITS_chi2_circle.push_back(part.its_chi2_circle);
            aux_ITS_chi2_helix.push_back(part.its_chi2_helix);
            aux_ITS_FH_Px.push_back(part.its_fh_px);
            aux_ITS_FH_Py.push_back(part.its_fh_py);
            aux_ITS_FH_Pt.push_back(part.its_fh_pt);
            aux_ITS_FH_Pz.push_back(part.its_fh_pz);
            aux_ITS_FH_P.push_back(part.its_fh_p);
            aux_ITS_LH_Px.push_back(part.its_lh_px);
            aux_ITS_LH_Py.push_back(part.its_lh_py);
            aux_ITS_LH_Pt.push_back(part.its_lh_pt);
            aux_ITS_LH_Pz.push_back(part.its_lh_pz);
            aux_ITS_LH_P.push_back(part.its_lh_p);
            aux_nTPChits.push_back((Int_t)part.tpc_hits.size());
            aux_TPC_Pt_rec.push_back(part.tpc_pt_rec);
            aux_TPC_Pz_rec.push_back(part.tpc_pz_rec);
            aux_TPC_P_rec.push_back(part.tpc_p_rec);
            aux_TPC_charge_rec.push_back(part.tpc_charge_rec);
            aux_TPC_chi2_circle.push_back(part.tpc_chi2_circle);
            aux_TPC_chi2_helix.push_back(part.tpc_chi2_helix);
            aux_TPC_dx.push_back(part.tpc_dx);
            aux_TPC_dEdx.push_back(part.tpc_dE_dx);
            aux_TPC_FH_Px.push_back(part.tpc_fh_px);
            aux_TPC_FH_Py.push_back(part.tpc_fh_py);
            aux_TPC_FH_Pt.push_back(part.tpc_fh_pt);
            aux_TPC_FH_Pz.push_back(part.tpc_fh_pz);
            aux_TPC_FH_P.push_back(part.tpc_fh_p);
        }

        // at the end of the event
        if ((Int_t)aux_trackID.size()) output_tree->Fill();

        // clear vectors
        aux_trackID.clear();
        aux_PDGcode.clear();
        aux_Px_ini.clear();
        aux_Py_ini.clear();
        aux_Pz_ini.clear();
        aux_Pt_ini.clear();
        aux_P_ini.clear();
        aux_parentID.clear();
        aux_nDaughters.clear();
        aux_isPrimary.clear();
        aux_charge.clear();
        aux_nITShits.clear();
        aux_ITS_Pt_rec.clear();
        aux_ITS_Pz_rec.clear();
        aux_ITS_P_rec.clear();
        aux_ITS_charge_rec.clear();
        aux_ITS_chi2_circle.clear();
        aux_ITS_chi2_helix.clear();
        aux_ITS_FH_Px.clear();
        aux_ITS_FH_Py.clear();
        aux_ITS_FH_Pt.clear();
        aux_ITS_FH_Pz.clear();
        aux_ITS_FH_P.clear();
        aux_ITS_LH_Px.clear();
        aux_ITS_LH_Py.clear();
        aux_ITS_LH_Pt.clear();
        aux_ITS_LH_Pz.clear();
        aux_ITS_LH_P.clear();
        aux_nTPChits.clear();
        aux_TPC_Pt_rec.clear();
        aux_TPC_Pz_rec.clear();
        aux_TPC_P_rec.clear();
        aux_TPC_charge_rec.clear();
        aux_TPC_chi2_circle.clear();
        aux_TPC_chi2_helix.clear();
        aux_TPC_dx.clear();
        aux_TPC_dEdx.clear();
        aux_TPC_FH_Px.clear();
        aux_TPC_FH_Py.clear();
        aux_TPC_FH_Pt.clear();
        aux_TPC_FH_Pz.clear();
        aux_TPC_FH_P.clear();
    }

    output_tree->Write();

    output_file->Save();
    output_file->Close();
}  // end of macro
