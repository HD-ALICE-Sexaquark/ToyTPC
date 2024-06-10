#include "TMath.h"
#include "TVector3.h"

#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"

/**
 * @brief Given the parametric form of a helix, evaluate it at a given parameter `t`,
 * and return the square of its distance w.r.t. an external point.
 *
 * @param t Evaluation parameter.
 * @param point External point coordinates. Format: `{x, y, z}`.
 * @param helix Helix parameters. Format: `{x_c, y_c, R, omega, phi}`.
 *
 * @returns The square of the distance between the point and the helix.
 */
Double_t SquaredDistancePointToHelix(const Double_t* t, Double_t point[], Double_t helix[]) {
    return TMath::Power(point[0] - helix[2] * TMath::Cos(helix[3] * t[0] + helix[4]) - helix[0], 2) +
           TMath::Power(point[1] - helix[2] * TMath::Sin(helix[3] * t[0] + helix[4]) - helix[1], 2) +  //
           TMath::Power(point[2] - t[0], 2);
}

/**
 * @brief Calculate the distance of closest approach between a point and a helix.
 *
 * @param point External point coordinates. Format: `{x, y, z}`.
 * @param helix_params Helix parameters. Format: `{x_c, y_c, R, omega, phi}`.
 * @param sol (set by reference) Solution of the minimization.
 * @param PCA (set by reference) Point of Closest Approach.
 *
 * @returns Distance of Closest Approach (DCA) between a point and a helix.
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

/**
 * Overload of `CalculatePointToHelixDCA()`, without the need of returning the solution and PCA.
 */
Double_t CalculatePointToHelixDCA(Double_t point[], Double_t helix_params[]) {
    Double_t sol;
    TVector3 PCA;
    return CalculatePointToHelixDCA(point, helix_params, sol, PCA);
}
