#include "TMath.h"
#include "TVector3.h"

#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"

/*
 Given the parametric form of a helix, find the distance between it and a point.
 - Variable: `t`, the main parameter of the helix.
 - Parameters: `point`, the point to which we want to calculate the distance.
 - Parameters: `linePoint`, `lineDir`, the position and direction of the line.
 - Return: the square of the distance between the point and the line.
*/
Double_t SquaredDistancePointToHelix(const Double_t* t, Double_t point[], Double_t x_c, Double_t y_c, Double_t r, Double_t omega,
                                     Double_t phi) {
    return TMath::Power(point[0] - r * TMath::Cos(omega * t[0] + phi) - x_c, 2) +  //
           TMath::Power(point[1] - r * TMath::Sin(omega * t[0] + phi) - y_c, 2) +  //
           TMath::Power(point[2] - omega * t[0], 2);
}

/*
 Calculate the distance of closest approach between a point and a helix.
 - Input: `point`, the point to which we want to calculate the distance
 - Output: `PCA`, the point of closest approach
 - Return: the distance of closest approach
 */
Double_t CalculatePointToHelixDCA(Double_t point[], Double_t helix_params[], TVector3& PCA) {

    Double_t ref[3] = {point[0], point[1], point[2]};
    Double_t x_c = helix_params[0];
    Double_t y_c = helix_params[1];
    Double_t r = helix_params[2];
    Double_t omega = helix_params[3];
    Double_t phi = helix_params[4];

    auto func = [&ref, &x_c, &y_c, &r, &omega, &phi](const Double_t* t) {
        return SquaredDistancePointToHelix(t, ref, x_c, y_c, r, omega, phi);
    };
    ROOT::Math::Functor f(func, 1);

    ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit", "Migrad");
    min->SetFunction(f);
    min->SetVariable(0, "t", -999., 0.01);
    min->Minimize();

    const Double_t* xs = min->X();

    PCA.SetXYZ(r * TMath::Cos(omega * xs[0] + phi) + x_c, r * TMath::Sin(omega * xs[0] + phi) + y_c, xs[0]);

    return TMath::Sqrt(SquaredDistancePointToHelix(xs, ref, x_c, y_c, r, omega, phi));
}
