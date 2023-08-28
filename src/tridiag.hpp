#ifndef TRIDIAG_HPP
#define TRIDIAG_HPP

#include <vector>

void tridiag(double a, double b, double c, std::vector<double> &r, std::vector<double> &u)
{
    int j;
    double bet;

    int n = r.size();
    std::vector<double> gam(n);

    u[0] = r[0] / (bet = b);

    for (j = 1; j < n; j++)
    {
        gam[j] = c / bet;
        bet = b - a * gam[j];
        u[j] = (r[j] - a * u[j - 1]) / bet;
    }
    for (j = (n - 2); j >= 0; j--)
        u[j] -= gam[j + 1] * u[j + 1];
}

#endif