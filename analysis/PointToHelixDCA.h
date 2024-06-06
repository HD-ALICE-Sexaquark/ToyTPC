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
Double_t SquaredDistancePointToHelix(const Double_t* t, Double_t point[], Double_t helix[]) {
    return TMath::Power(point[0] - helix[2] * TMath::Cos(helix[3] * t[0] + helix[4]) - helix[0], 2) +  //
           TMath::Power(point[1] - helix[2] * TMath::Sin(helix[3] * t[0] + helix[4]) - helix[1], 2) +  //
           TMath::Power(point[2] - t[0], 2);
}

/*
 Calculate the distance of closest approach between a point and a helix.
 - Input: `point`, the point to which we want to calculate the distance
 - Output: `PCA`, the point of closest approach
 - Return: the distance of closest approach
*/
Double_t CalculatePointToHelixDCA(Double_t point[], Double_t helix_params[], Double_t& sol, TVector3& PCA) {

    Double_t ref[3] = {point[0], point[1], point[2]};
    Double_t helix[5] = {helix_params[0], helix_params[1], helix_params[2], helix_params[3], helix_params[4]};

    auto func = [&ref, &helix](const Double_t* t) { return SquaredDistancePointToHelix(t, ref, helix); };
    ROOT::Math::Functor f(func, 1);

    ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit", "Migrad");
    min->SetFunction(f);
    min->SetVariable(0, "t", 0., 0.01);
    min->Minimize();

    sol = min->X()[0];
    PCA.SetXYZ(helix[2] * TMath::Cos(helix[3] * sol + helix[4]) + helix[0], helix[2] * TMath::Sin(helix[3] * sol + helix[4]) + helix[1],
               sol);

    return TMath::Sqrt(SquaredDistancePointToHelix(min->X(), ref, helix));
}
