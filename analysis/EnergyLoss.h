#include <iostream>

#include "TMath.h"

/*
  Determine the differential energy loss of a particle with a certain momentum.
*/
Int_t EnergyLoss(Int_t N, Double_t *depositedEnergy, Double_t *x, Double_t *y, Double_t *z, Double_t &dE_dx, Double_t &dxaverage) {

    Double_t dx[N - 1];
    Double_t diffEnergyLoss[N - 1];
    for (Int_t i = 1; i < N; i++) {
        // calculate the distance between two hits to get dx, but start at 1 because we need an initial distance
        // dE is given by the deposited energy in MeV
        dx[i - 1] = (x[i] - x[i - 1]) * (x[i] - x[i - 1]) + (y[i] - y[i - 1]) * (y[i] - y[i - 1]) + (z[i] - z[i - 1]) * (z[i] - z[i - 1]);
        dx[i - 1] = TMath::Sqrt(dx[i - 1]);
        // calculate all differential energy losses
        diffEnergyLoss[i - 1] = depositedEnergy[i] / dx[i - 1];
    }

    // sort diffEnergyLoss
    Int_t ind[N - 1];
    TMath::Sort(N - 1, diffEnergyLoss, ind);

    // Cut out 50 highest and 50 lowest energy losses in sum
    // now take the average of the remaining diffEnergyLoss array
    Double_t sum = 0;
    Double_t sum_dx = 0;
    for (Int_t i = 50; i < N - 51; i++) {
        sum += diffEnergyLoss[ind[i]];
        sum_dx += dx[ind[i]];
    }

    // take average to obtain [dE/dx] in MeV/cm
    Double_t average = sum / (N - 101);
    Double_t average_dx = sum_dx / (N - 101);
    dE_dx = average;
    dxaverage = average_dx * 0.01;

    return 0;
}
