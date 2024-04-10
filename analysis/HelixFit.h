#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>

#include "Fit/Fitter.h"
#include "Math/Functor.h"
#include "TMath.h"

Double_t costFunc(Int_t N, Double_t *x, Double_t *y, Double_t *z, Double_t x_c, Double_t y_c, Double_t r, Double_t omega, Double_t phi) {
    // minimisation function computing the sum of squares of residuals looping at the graph points
    Double_t f = 0;
    for (Int_t i = 0; i < N; i++) {
        Double_t u = x[i] - x_c;
        Double_t v = y[i] - y_c;
        Double_t dx = u - r * TMath::Cos(omega * z[i] + phi);
        Double_t dy = v - r * TMath::Sin(omega * z[i] + phi);
        f += dx * dx + dy * dy;
    }
    return f;
}

// Fit z -> (rcos(w*z+phi),rsin(w*z+phi), z) with slope r*w and angle arctan(r*w)
Bool_t HelixFit(Int_t N, Double_t *x, Double_t *y, Double_t *z, Double_t x_c, Double_t y_c, Double_t r, Double_t &angle, Double_t &charge,
                Double_t &chi2, Int_t &direction) {

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

    Double_t omega = result.Parameter(0);
    Double_t phi = result.Parameter(1);

    angle = TMath::ATan((1 / (r * omega)));
    direction = (Int_t)TMath::Sign(1, z[1] - z[0]);
    charge = -1 * TMath::Sign(1.0, omega) * direction;

    chi2 = costFunc(N, x, y, z, x_c, y_c, r, omega, phi);

    return ok;
}