#include <fstream>
#include <iostream>
#include <map>
#include <sstream>

#include "Fit/Fitter.h"
#include "Math/Functor.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TRandom3.h"
#include "TString.h"
#include "TTree.h"
#include "TVector3.h"

/*** STRUCTURES ***/

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

/*** FUNCTIONS ***/

/*
 Smear px, py or pz
 The parameters for sigma_p were obtained from the ALICE performance paper
 - after extracting the points from Fig. 23
 - using Eq. 14 to get sigma_pt/pt
 - then made a linear fit on those values
 (reference: Int. J. Mod. Phys. A 2014.29)
 NOTE: it says resolution of pt, but I'm generalizing it for px,py,pz...
 - Units: in GeV/c
 */
Double_t SmearMomentum(Double_t p, TRandom3 *rnd) {

    std::cerr << "ParseCVSFiles.C :: SmearMomentum :: >> p = " << p << std::endl;

    Double_t sigma_p = 0.0053 * p * p + 0.0032 * TMath::Abs(p);
    std::cerr << "ParseCVSFiles.C :: SmearMomentum :: >> sigma_p = " << sigma_p << std::endl;

    // protection
    if (sigma_p < 1E-4) {
        return p;
    }

    TF1 function_p("function_p", "gaus", 0., 10.);
    function_p.SetParameter(0, 1.);
    function_p.SetParameter(1, p);
    function_p.SetParameter(2, sigma_p);

    Double_t smeared_p = function_p.GetRandom(rnd);
    std::cerr << "ParseCVSFiles.C :: SmearMomentum :: >> smeared_p = " << smeared_p << std::endl;

    return smeared_p;
}

/*
 Try to fit a circle from the TPC hits...
 NOTE: currently not working...
 */
Bool_t FitCircle(const Int_t nnn, Double_t x[], Double_t y[], Double_t &radius) {

    TGraph graph(nnn, x, y);

    auto chi2Function = [&](const Double_t *par) {
        // minimisation function computing the sum of squares of residuals looping at the graph points
        Int_t np = graph.GetN();
        Double_t f = 0;
        Double_t *x = graph.GetX();
        Double_t *y = graph.GetY();
        for (Int_t i = 0; i < np; i++) {
            Double_t u = x[i] - par[0];
            Double_t v = y[i] - par[1];
            Double_t dr = par[2] - TMath::Sqrt(u * u + v * v);
            f += dr * dr;
        }
        return f;
    };

    // wrap chi2 function in a function object for the fit
    // 3 is the number of fit parameters (size of array par)
    ROOT::Math::Functor fcn(chi2Function, 3);
    ROOT::Fit::Fitter fitter;
    Double_t pStart[3] = {0, 0, 1};
    fitter.SetFCN(fcn, pStart);
    fitter.Config().ParSettings(0).SetName("x0");
    fitter.Config().ParSettings(1).SetName("y0");
    fitter.Config().ParSettings(2).SetName("R");

    // do the fit
    Bool_t ok = fitter.FitFCN();
    if (!ok) {
        Error("line3Dfit", "Line3D Fit failed");
    }
    const ROOT::Fit::FitResult &result = fitter.Result();
    // result.Print(std::cout);

    radius = result.Parameter(2);
    return ok;
}

/*** MAIN ***/

void ParseCSVFiles(Int_t run_n = 0) {

    // set input/output filenames -- hardcoded
    TString input_traj_file = Form("run%03d_traj.csv", run_n);
    TString input_its_file = Form("run%03d_its.csv", run_n);
    TString input_tpc_file = Form("run%03d_tpc.csv", run_n);
    TString output_filename = Form("run%03d_ana.root", run_n);

    // prepare TRandom object, necessary later for smearing
    TRandom3 *rnd = new TRandom3(0.);

    // conversion factors
    const Double_t MeVToGeV = 1E-3;
    const Double_t GeVToMeV = 1E3;

    /*** Part 1: Read and Store Input ***/

    // init particles -- main container
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

    /* (Auxiliary variables) */

    Particle_tt aux_particle;
    Int_t current_eventID;
    Int_t prev_eventID = -1;

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
            // resize vectors, because eventID starts at 0
            Events.resize(current_eventID + 1);
            part_index.resize(current_eventID + 1);
            part_ndaughters.resize(current_eventID + 1);
            Events.at(current_eventID).eventID = current_eventID;
            prev_eventID = current_eventID;  // we're in a different event now
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

        Events.at(current_eventID).particles.push_back(aux_particle);

        // update maps
        part_index[current_eventID][aux_particle.trackID] = (Int_t)Events.at(current_eventID).particles.size() - 1;
        part_ndaughters[current_eventID][aux_particle.parentID]++;
    }  // end of reading file

    traj_file.close();

    /* Posteriori variable: N Daughters */

    for (Event_tt &evt : Events) {
        for (Particle_tt &part : evt.particles) {
            part.n_daughters = part_ndaughters[evt.eventID][part.trackID];
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
        aux_index = part_index[current_eventID][aux_track_id];

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
        aux_index = part_index[current_eventID][aux_track_id];

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

            if ((Int_t)part.tpc_hits.size() >= 100) {
                std::cout << "ParseCSVFiles.C :: event,part " << evt.eventID << ", " << part.trackID << std::endl;
                std::cout << "ParseCSVFiles.C :: >> pdg,is_primary " << part.PDGcode << ", " << part.is_primary << std::endl;
                std::cout << "ParseCSVFiles.C :: >> parent,n_daughters " << part.parentID << ", " << part.n_daughters << std::endl;
                std::cout << "ParseCSVFiles.C :: >> x,y,z " << part.x_ini << ", " << part.y_ini << ", " << part.z_ini << ", " << std::endl;
                std::cout << "ParseCSVFiles.C :: >> px,py,pz " << part.px_ini << ", " << part.py_ini << ", " << part.pz_ini << ", "
                          << std::endl;
                std::cout << "ParseCSVFiles.C :: >> n_its_hits " << (Int_t)part.its_hits.size() << std::endl;
                std::cout << "ParseCSVFiles.C :: >> n_tpc_hits " << (Int_t)part.tpc_hits.size() << std::endl;

                Double_t tpc_x[(Int_t)part.tpc_hits.size()];
                Double_t tpc_y[(Int_t)part.tpc_hits.size()];
                for (Int_t hit_idx = 0; hit_idx < (Int_t)part.tpc_hits.size(); hit_idx++) {
                    tpc_x[hit_idx] = part.tpc_hits[hit_idx].x;
                    tpc_y[hit_idx] = part.tpc_hits[hit_idx].y;
                }
                Double_t radius;
                FitCircle((Int_t)part.tpc_hits.size(), tpc_x, tpc_y, radius);
                std::cout << "ParseCSVFiles.C :: >> radius " << radius << std::endl;
                std::cout << "ParseCSVFiles.C :: >> pt_ini " << TMath::Sqrt(TMath::Power(part.px_ini, 2) + TMath::Power(part.py_ini, 2))
                          << std::endl;
                std::cout << "ParseCSVFiles.C :: >> pt (from radius) " << 0.3 * 0.2 * 0.1 * radius << std::endl;
            }

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

    // prepare output
    TFile *output_file = new TFile(output_filename, "RECREATE");
    TTree *output_tree = new TTree("Events", "Events");

    // set aux objects
    Int_t aux_eventID;
    std::vector<Int_t> aux_trackID;
    std::vector<Int_t> aux_PDGcode;
    std::vector<Bool_t> aux_isPrimary;
    std::vector<Int_t> aux_parentID;
    std::vector<Int_t> aux_nDaughters;
    std::vector<Int_t> aux_nITShits;
    std::vector<Int_t> aux_nTPChits;
    std::vector<Float_t> aux_Px_ini;
    std::vector<Float_t> aux_Py_ini;
    std::vector<Float_t> aux_Pz_ini;
    std::vector<Float_t> aux_Pt_ini;
    std::vector<Float_t> aux_Pt_sme;
    std::vector<Float_t> aux_Pt_rec;
    std::vector<Float_t> aux_dEdx;

    // set tree branches
    output_tree->Branch("eventID", &aux_eventID);
    output_tree->Branch("trackID", &aux_trackID);
    output_tree->Branch("PDGcode", &aux_PDGcode);
    output_tree->Branch("isPrimary", &aux_isPrimary);
    output_tree->Branch("parentID", &aux_parentID);
    output_tree->Branch("nDaughters", &aux_nDaughters);
    output_tree->Branch("nITShits", &aux_nITShits);
    output_tree->Branch("nTPChits", &aux_nTPChits);
    output_tree->Branch("Px_ini", &aux_Px_ini);
    output_tree->Branch("Py_ini", &aux_Py_ini);
    output_tree->Branch("Pz_ini", &aux_Pz_ini);
    output_tree->Branch("Pt_ini", &aux_Pt_ini);
    output_tree->Branch("Pt_sme", &aux_Pt_sme);
    output_tree->Branch("Pt_rec", &aux_Pt_rec);
    output_tree->Branch("dEdx", &aux_dEdx);

    for (Event_tt &evt : Events) {

        aux_eventID = evt.eventID;

        for (Particle_tt &part : evt.particles) {

            // (cut)
            if ((Int_t)part.tpc_hits.size() < 100) {
                continue;
            }

            aux_trackID.push_back(part.trackID);
            aux_PDGcode.push_back(part.PDGcode);
            aux_isPrimary.push_back(part.is_primary);
            aux_parentID.push_back(part.parentID);
            aux_nDaughters.push_back(part.n_daughters);
            aux_nITShits.push_back((Int_t)part.its_hits.size());
            aux_nTPChits.push_back((Int_t)part.tpc_hits.size());
            aux_Px_ini.push_back(part.px_ini);
            aux_Py_ini.push_back(part.py_ini);
            aux_Pz_ini.push_back(part.pz_ini);
            aux_Pt_ini.push_back(TMath::Sqrt(part.px_ini * part.px_ini + part.py_ini * part.py_ini));
            aux_Pt_sme.push_back(SmearMomentum(aux_Pt_ini.back() * MeVToGeV, rnd) * GeVToMeV);
            aux_Pt_rec.push_back(0.);  // (pending)
            aux_dEdx.push_back(0.);    // (pending)
        }

        // at the end of the event
        if ((Int_t)aux_trackID.size()) output_tree->Fill();

        // clear vectors
        aux_trackID.clear();
        aux_PDGcode.clear();
        aux_isPrimary.clear();
        aux_parentID.clear();
        aux_nDaughters.clear();
        aux_nITShits.clear();
        aux_nTPChits.clear();
        aux_Px_ini.clear();
        aux_Py_ini.clear();
        aux_Pz_ini.clear();
        aux_Pt_ini.clear();
        aux_Pt_sme.clear();
        aux_Pt_rec.clear();
        aux_dEdx.clear();
    }

    output_tree->Write();

    output_file->Save();
    output_file->Close();
}  // end of macro
