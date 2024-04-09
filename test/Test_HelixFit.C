#include <fstream>
#include <iostream>
#include <map>
#include <sstream>

#include "TArc.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TGraph2D.h"
#include "TH2.h"
#include "TObjArray.h"
#include "TString.h"

#include "Fit/Fitter.h"
#include "Math/Functor.h"

#include "../analysis/LeastSquaresCircleFit.h"
// #include "../analysis/jonas_04.04/Helixfit_root.h"

struct MCParticle_tt {
    Double_t pid;
    Double_t px, py, pz;
};

struct Particle_tt {
    std::vector<Double_t> x_TPC_hits;
    std::vector<Double_t> y_TPC_hits;
    std::vector<Double_t> z_TPC_hits;
};

// Fit z -> (rcos(w*z+phi),rsin(w*z+phi), z) with slope r*w and angle arctan(r*w)
void HelixFit(Int_t N, Double_t *x, Double_t *y, Double_t *z,  // hits
              Double_t x_c, Double_t y_c, Double_t r,          // initial guess
              Double_t &w, Double_t &phi                       // output
) {

    auto chi2Function = [&](const Double_t *par) {
        // minimisation function computing the sum of squares of residuals looping at the graph points
        Double_t f = 0;
        for (Int_t i = 0; i < N; i++) {
            Double_t u = x[i] - x_c;
            Double_t v = y[i] - y_c;
            Double_t dx = u - r * TMath::Cos(par[0] * z[i] + par[1]);
            Double_t dy = v - r * TMath::Sin(par[0] * z[i] + par[1]);
            f += dx * dx + dy * dy;
        }
        return f;
    };

    // wrap chi2 function in a function object for the fit
    // 2 is the number of fit parameters (size of array par)
    ROOT::Math::Functor fcn(chi2Function, 2);
    ROOT::Fit::Fitter fitter;
    Double_t pStart[2] = {0, 0};
    fitter.SetFCN(fcn, pStart);
    fitter.Config().ParSettings(0).SetName("Frequency");
    fitter.Config().ParSettings(1).SetName("phase");

    // do the fit
    Bool_t ok = fitter.FitFCN();
    if (!ok) {
        Error("line3Dfit", "Line3D Fit failed");
    }
    const ROOT::Fit::FitResult &result = fitter.Result();

    w = result.Parameter(0);
    phi = result.Parameter(1);
}

void HelixFit_Wiki(Int_t N, Double_t *x, Double_t *y, Double_t *z,             // hits
                   Double_t x_c, Double_t y_c, Double_t r,                     // initial guess
                   Double_t &omega, Double_t &phi, Double_t &b, Double_t &z_0  // output
) {

    auto chi2Function = [&](const Double_t *par) {
        // minimisation function computing the sum of squares of residuals looping at the graph points
        Double_t f = 0;
        for (Int_t i = 0; i < N; i++) {
            Double_t t = (Double_t)i;
            Double_t dx = x[i] - x_c - r * TMath::Cos(par[0] * t + par[1]);
            Double_t dy = y[i] - y_c - r * TMath::Sin(par[0] * t + par[1]);
            Double_t dz = z[i] - par[2] * t - par[3];
            f += dx * dx + dy * dy + dz * dz;
        }
        return f;
    };

    // wrap chi2 function in a function object for the fit
    // 2 is the number of fit parameters (size of array par)
    ROOT::Math::Functor fcn(chi2Function, 4);
    ROOT::Fit::Fitter fitter;
    Double_t pStart[4] = {0, 0, 0, 0};
    fitter.SetFCN(fcn, pStart);
    fitter.Config().ParSettings(0).SetName("Frequency");
    fitter.Config().ParSettings(1).SetName("phase");

    // do the fit
    Bool_t ok = fitter.FitFCN();
    if (!ok) {
        Error("line3Dfit", "Line3D Fit failed");
    }
    const ROOT::Fit::FitResult &result = fitter.Result();

    omega = result.Parameter(0);
    phi = result.Parameter(1);
    b = result.Parameter(2);
    z_0 = result.Parameter(3);
}

void Test_HelixFit() {

    std::map<Int_t, Particle_tt> particles;
    std::map<Int_t, MCParticle_tt> true_particles;
    Int_t aux_track_id;

    /* Read file containing True info */

    TString input_traj_file = "run00_traj.csv";
    std::fstream traj_file;
    traj_file.open(input_traj_file);
    if (!traj_file) {
        std::cerr << "ParseCSVFiles.C :: ERROR :: File " << input_traj_file << " not found." << std::endl;
        return;
    }

    // read line by line ~ loop over hits
    TObjArray *token = nullptr;
    std::string line;
    line.clear();
    while (std::getline(traj_file, line)) {

        // (protection)
        if (line == "") {
            std::cerr << "ParseCSVFiles.C :: WARNING :: Line empty, skipping..." << std::endl;
            continue;
        }

        // convert the line into a TString, then split
        token = ((TString)line).Tokenize(",");

        aux_track_id = ((TString)(token->At(1)->GetName())).Atoi();                      //   [1] trackID
        true_particles[aux_track_id].pid = ((TString)(token->At(2)->GetName())).Atoi();  //   [2] pid
        true_particles[aux_track_id].px = ((TString)(token->At(6)->GetName())).Atof();   //   [6] px
        true_particles[aux_track_id].py = ((TString)(token->At(7)->GetName())).Atof();   //   [7] py
        true_particles[aux_track_id].pz = ((TString)(token->At(8)->GetName())).Atof();   //   [8] pz
    }                                                                                    // end of reading file

    traj_file.close();

    /* Read file containing TPC hits */

    TString input_tpc_file = "run00_tpc_cropped.csv";
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

        aux_track_id = ((TString)(token->At(1)->GetName())).Atoi();                                 //   [1] trackID
        particles[aux_track_id].x_TPC_hits.push_back(((TString)(token->At(2)->GetName())).Atof());  //   [2] x
        particles[aux_track_id].y_TPC_hits.push_back(((TString)(token->At(3)->GetName())).Atof());  //   [3] y
        particles[aux_track_id].z_TPC_hits.push_back(((TString)(token->At(4)->GetName())).Atof());  //   [4] z
    }                                                                                               // end of reading file

    tpc_file.close();

    /* choose particle */

    Int_t particle_to_test = 12;
    Int_t N = particles[particle_to_test].x_TPC_hits.size();

    std::cout << "Chosen particle: " << particle_to_test << std::endl;
    std::cout << "# True Info:" << std::endl;
    std::cout << "  -> pid = " << true_particles[particle_to_test].pid << std::endl;
    std::cout << "  -> px  = " << true_particles[particle_to_test].px * 1E-3 << std::endl;
    std::cout << "  -> py  = " << true_particles[particle_to_test].py * 1E-3 << std::endl;
    std::cout << "  -> pz  = " << true_particles[particle_to_test].pz * 1E-3 << std::endl;
    std::cout << "  => pt  = "
              << TMath::Sqrt(true_particles[particle_to_test].px * true_particles[particle_to_test].px +
                             true_particles[particle_to_test].py * true_particles[particle_to_test].py) *
                     1E-3
              << std::endl;
    std::cout << "Number of hits: " << N << std::endl;

    /* fit circle */

    Double_t x_c, y_c, r;
    LeastSquaresCircleFit(N, &particles[particle_to_test].x_TPC_hits[0], &particles[particle_to_test].y_TPC_hits[0], x_c, y_c, r);

    std::cout << "# From Circle Fit:" << std::endl;
    std::cout << "  -> x_c = " << x_c << std::endl;
    std::cout << "  -> y_c = " << y_c << std::endl;
    std::cout << "  -> r   = " << r << std::endl;
    std::cout << "  => pt  = " << 0.3 * 0.2 * r * 1E-2 << std::endl;

    /* retrieve circle from fit values */

    Double_t aux_y;
    std::vector<Double_t> x_fitted;
    std::vector<Double_t> y_fitted;
    std::vector<Double_t> z_fitted;

    for (Int_t i = 0; i < N; i++) {

        aux_y = y_c +
                TMath::Sqrt(r * r - (particles[particle_to_test].x_TPC_hits[i] - x_c) * (particles[particle_to_test].x_TPC_hits[i] - x_c));

        // std::cout << particles[particle_to_test].x_TPC_hits[i] << " " << aux_y << " " << 0 << std::endl;

        x_fitted.push_back(particles[particle_to_test].x_TPC_hits[i]);
        y_fitted.push_back(aux_y);
        z_fitted.push_back(1);
    }

    /* fit helix */

    /* Double_t w, phi;
    HelixFit(N, &particles[particle_to_test].x_TPC_hits[0], &particles[particle_to_test].y_TPC_hits[0],
             &particles[particle_to_test].z_TPC_hits[0], x_c, y_c, r, w, phi); */
    Double_t omega, phi, b, z_0;
    HelixFit_Wiki(N, &particles[particle_to_test].x_TPC_hits[0], &particles[particle_to_test].y_TPC_hits[0],
                  &particles[particle_to_test].z_TPC_hits[0], x_c, y_c, r, omega, phi, b, z_0);

    /* retrieve helix from fit values */

    Double_t aux_x;
    Double_t aux_z;
    std::vector<Double_t> x_h_fitted;
    std::vector<Double_t> y_h_fitted;
    std::vector<Double_t> z_h_fitted;

    for (Int_t i = 0; i < N; i++) {

        /*         aux_x = x_c + r * TMath::Cos(omega * particles[particle_to_test].z_TPC_hits[i] + phi);
                aux_y = y_c + r * TMath::Sin(omega * particles[particle_to_test].z_TPC_hits[i] + phi);
                aux_z = particles[particle_to_test].z_TPC_hits[i];
         */

        aux_x = x_c + r * TMath::Cos(omega * i + phi);
        aux_y = y_c + r * TMath::Sin(omega * i + phi);
        aux_z = b * i + z_0;

        x_h_fitted.push_back(aux_x);
        y_h_fitted.push_back(aux_y);
        z_h_fitted.push_back(aux_z + 1);
    }

    std::cout << "# From Helix Fit:" << std::endl;
    std::cout << "  -> omega = " << omega << std::endl;
    std::cout << "  -> phi   = " << phi << std::endl;
    std::cout << "  -> b     = " << b << std::endl;
    std::cout << "  -> z_0   = " << z_0 << std::endl;
    std::cout << "  => a / b = " << r / b << std::endl;
    std::cout << "  => pz    = " << TMath::ATan(omega * r / b) << std::endl;

    /* plot */

    TGraph2D *graph = new TGraph2D(N, &particles[particle_to_test].x_TPC_hits[0], &particles[particle_to_test].y_TPC_hits[0],
                                   &particles[particle_to_test].z_TPC_hits[0]);

    graph->SetMarkerStyle(7);
    graph->SetMarkerSize(2);
    graph->SetMarkerColor(kBlue);

    graph->SetMinimum(-200);
    graph->SetMaximum(0);

    TGraph2D *fitted_circle = new TGraph2D(N, &x_fitted[0], &y_fitted[0], &z_fitted[0]);

    fitted_circle->SetMarkerStyle(7);
    fitted_circle->SetMarkerSize(2);
    fitted_circle->SetMarkerColor(kRed);

    TGraph2D *fitted_helix = new TGraph2D(N, &x_h_fitted[0], &y_h_fitted[0], &z_h_fitted[0]);

    fitted_helix->SetMarkerStyle(7);
    fitted_helix->SetMarkerSize(2);
    fitted_helix->SetMarkerColor(kMagenta);

    /* canvas */

    TCanvas *c1 = new TCanvas("c1", "c1", 800, 800);

    graph->Draw("P");

    graph->GetXaxis()->SetLimits(-200, 0);
    graph->GetYaxis()->SetLimits(-260, -60);

    fitted_circle->Draw("P SAME");
    fitted_helix->Draw("P SAME");
}
