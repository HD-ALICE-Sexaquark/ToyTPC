#include <iostream>

#include "TMath.h"

Double_t CHI_squared(Int_t N, Double_t *x, Double_t *y, Double_t x_c, Double_t y_c, Double_t r) {

    // get the smallest distance of all points from the circle

    Double_t dist[N];

    Double_t res = 0;  // residual

    // compute distance from loop center and compare to radius of the loop and residual
    for (Int_t i = 0; i < N; i++) {
        dist[i] = r - TMath::Sqrt((x[i] - x_c) * (x[i] - x_c) + (y[i] - y_c) * (y[i] - y_c));

        res += dist[i] * dist[i];
    }

    return res;
}

/*
 A "new" method of making a least-squares circle fit to a set of points.
 - Input: `x`, `y` arrays of length `N`
 - Output: `x_c`, `y_c`, `r`
 - Return: `0` if successful, `1` if determinant is too small
 * Reference: https://dtcenter.org/community-code/model-evaluation-tools-met/documentation
*/
Int_t LeastSquaresCircleFit(Int_t N, Double_t *x, Double_t *y, Double_t &x_c, Double_t &y_c, Double_t &r, Double_t &chi2) {

    /* Get average */

    Double_t avg_x = 0.;
    Double_t avg_y = 0.;

    for (Int_t i = 0; i < N; i++) {
        avg_x += x[i];
        avg_y += y[i];
    }

    avg_x /= N;
    avg_y /= N;

    /* Define: u = x - avg_x, v = y - avg_y */

    Double_t u[N], v[N];

    for (Int_t i = 0; i < N; i++) {
        u[i] = x[i] - avg_x;
        v[i] = y[i] - avg_y;
    }

    /* Define sums */

    Double_t Suu = 0.;
    Double_t Suv = 0.;
    Double_t Svv = 0.;
    Double_t Suuu = 0.;
    Double_t Suvv = 0.;
    Double_t Svvv = 0.;
    Double_t Svuu = 0.;

    for (Int_t i = 0; i < N; i++) {
        Suu += u[i] * u[i];
        Suv += u[i] * v[i];
        Svv += v[i] * v[i];
        Suuu += u[i] * u[i] * u[i];
        Suvv += u[i] * v[i] * v[i];
        Svvv += v[i] * v[i] * v[i];
        Svuu += v[i] * u[i] * u[i];
    }

    /* Calculate determinant */

    Double_t det = Suu * Svv - Suv * Suv;
    if (det < 1E-6) {
        std::cerr << "The determinant is too small!" << std::endl;
        return 1;
    }

    /* Calculate the solution */

    Double_t u_c = (0.5 * (Suuu + Suvv) * Svv - 0.5 * (Svvv + Svuu) * Suv) / det;
    Double_t v_c = (0.5 * (Svvv + Svuu) * Suu - 0.5 * (Suuu + Suvv) * Suv) / det;

    Double_t alpha = u_c * u_c + v_c * v_c + (Suu + Svv) / N;

    /* Define circle properties */

    x_c = u_c + avg_x;
    y_c = v_c + avg_y;
    r = TMath::Sqrt(alpha);
    chi2 = CHI_squared(N, x, y, x_c, y_c, r);

    return 0;
}
