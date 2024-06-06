#include "TMath.h"
#include "TVector3.h"

#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"

Double_t SquaredDistancePointToLine(const Double_t* t, Double_t point[], Double_t linePoint[], Double_t lineDir[]) {
    return TMath::Power(point[0] - linePoint[0] - t[0] * lineDir[0], 2) +  //
           TMath::Power(point[1] - linePoint[1] - t[0] * lineDir[1], 2) +  //
           TMath::Power(point[2] - linePoint[2] - t[0] * lineDir[2], 2);
}

Double_t CalculatePointToLineDCA(Double_t point[], Double_t line_pos[], Double_t line_dir[], Double_t& sol, TVector3& PCA) {

    Double_t ref[3] = {point[0], point[1], point[2]};
    Double_t pos[3] = {line_pos[0], line_pos[1], line_pos[2]};
    Double_t dir[3] = {line_dir[0], line_dir[1], line_dir[2]};

    // lambda function
    auto func = [&ref, &pos, &dir](const Double_t* t) { return SquaredDistancePointToLine(t, ref, pos, dir); };
    ROOT::Math::Functor f(func, 1);

    // initialize minimizer
    ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit", "Migrad");
    min->SetFunction(f);
    min->SetVariable(0, "t", 0., 0.01);
    min->Minimize();

    const Double_t* xs = min->X();
    PCA.SetXYZ(pos[0] + xs[0] * dir[0], pos[1] + xs[0] * dir[1], pos[2] + xs[0] * dir[2]);

    return TMath::Sqrt(SquaredDistancePointToLine(xs, ref, pos, dir));
}

void SimpleMinimizer() {

    Double_t pv[3] = {0., 0., 0.};
    Double_t line_pos[3] = {1., 1., 0.};
    Double_t line_dir[3] = {1., 0., 0.};

    Double_t sol;
    TVector3 PCA;

    Double_t dca = CalculatePointToLineDCA(pv, line_pos, line_dir, sol, PCA);

    std::cout << "DCA: " << dca << std::endl;
    std::cout << "PCA: (" << PCA.X() << ", " << PCA.Y() << ", " << PCA.Z() << ")" << std::endl;
    std::cout << "Solution: " << sol << std::endl;
}
