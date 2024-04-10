#include <iostream>

#include "Fit/Fitter.h"
#include "Math/Functor.h"
#include "TMath.h"

/*
  Minimisation function computing the sum of squares of residuals looping at the graph points.
*/
Double_t HelixMinFcn(Int_t N, Double_t *x, Double_t *y, Double_t *z, Double_t x_c, Double_t y_c, Double_t r, Double_t omega, Double_t phi) {
    Double_t dx, dy;
    Double_t res = 0.;
    for (Int_t i = 0; i < N; i++) {
        dx = x[i] - x_c - r * TMath::Cos(omega * z[i] + phi);
        dy = y[i] - y_c - r * TMath::Sin(omega * z[i] + phi);
        res += dx * dx + dy * dy;
    }
    return res;
}

/*
  From a collection of `N` points, and using `x_c`, `y_c` and `r` from a previous circle fit, fit a helix with the
  following parametrization:
  `(x, y, z) = (r*cos(omega*z + phi) + x_c, r*sin(omega*z + phi) + y_c, z)`
  - Input: `x`, `y`, `z` arrays of length `N`, `x_c`, `y_c`, `r` initial guesses
  - Output: `angle`, `charge`, `chi2`, `direction`
  - Return: `true` if successful, `false` if fit failed
  Reference: https://www.geometrictools.com/Documentation/HelixFitting.pdf
*/
Bool_t HelixFit(Int_t N, Double_t *x, Double_t *y, Double_t *z, Double_t x_c, Double_t y_c, Double_t r, Double_t &angle, Double_t &charge,
                Int_t &direction, Double_t &chi2) {

    auto chi2Function = [&](const Double_t *par) { return HelixMinFcn(N, x, y, z, x_c, y_c, r, par[0], par[1]); };

    // wrap chi2 function in a function object for the fit
    ROOT::Math::Functor fcn(chi2Function, 2);
    ROOT::Fit::Fitter fitter;
    Double_t pStart[2] = {0, 0};
    fitter.SetFCN(fcn, pStart);
    fitter.Config().ParSettings(0).SetName("Frequency");
    fitter.Config().ParSettings(1).SetName("Phase");

    // do the fit
    Bool_t fit_status = fitter.FitFCN();
    if (!fit_status) {
        std::cerr << "Fit failed!" << std::endl;
        return fit_status;
    }
    const ROOT::Fit::FitResult &result = fitter.Result();

    Double_t omega = result.Parameter(0);
    Double_t phi = result.Parameter(1);

    angle = TMath::ATan((1 / (r * omega)));
    direction = (Int_t)TMath::Sign(1, z[1] - z[0]);
    charge = -1 * TMath::Sign(1.0, omega) * direction;

    chi2 = result.MinFcnValue();

    return fit_status;
}
