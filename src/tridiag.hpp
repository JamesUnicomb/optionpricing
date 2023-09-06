#ifndef TRIDIAG_HPP
#define TRIDIAG_HPP

#define MAX_MESH_SIZE 65535

#include <cmath>
#include <vector>
#include <cassert>
#include <iostream>

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

class TriDiagSolver
{
private:
    int n, nn, lnn;
    double *x, *y, *aDiag, *bDiag, *cDiag, *aDiagCpy, *bDiagCpy, *cDiagCpy;
    double a, b, c;

public:
    TriDiagSolver(int n, double a, double b, double c)
        : n(n), a(a), b(b), c(c)
    {
        nn = 0;
        while (nn < n && nn < MAX_MESH_SIZE)
        {
            nn = 2 * (nn + 1) - 1;
        }

        if (nn < n)
        {
            throw std::invalid_argument("mesh size too large");
        }

        lnn = log2(nn + 1) - 1;

        x = nn > 0 ? new double[nn] : nullptr;
        y = nn > 0 ? new double[nn] : nullptr;

        aDiag = nn > 0 ? new double[nn] : nullptr;
        bDiag = nn > 0 ? new double[nn] : nullptr;
        cDiag = nn > 0 ? new double[nn] : nullptr;

        aDiagCpy = nn > 0 ? new double[nn] : nullptr;
        bDiagCpy = nn > 0 ? new double[nn] : nullptr;
        cDiagCpy = nn > 0 ? new double[nn] : nullptr;

        int i;
        for (i = 0; i < n; i++)
        {
            y[i] = 0.0;
            bDiag[i] = bDiagCpy[i] = b;
            aDiag[i] = aDiagCpy[i] = a;
            cDiag[i] = cDiagCpy[i] = c;
        }

        for (i = n; i < nn; i++)
        {
            y[i] = 1.0;
            bDiag[i] = bDiagCpy[i] = 1.0;
            aDiag[i] = aDiagCpy[i] = 0.0;
            cDiag[i] = cDiagCpy[i] = 0.0;
        }

        // fix diagonal endpoints
        aDiag[0] = aDiagCpy[0] = 0.0;
        cDiag[n - 1] = cDiagCpy[n - 1] = 0.0;
    }

    void solve_parallel()
    {
        int i, j, index1, index2, offset;
        double alpha, gamma;

#pragma omp parallel
        for (i = 0; i < nn; i++)
        {
            x[i] = 0.0;
            bDiag[i] = bDiagCpy[i];
            aDiag[i] = aDiagCpy[i];
            cDiag[i] = cDiagCpy[i];
        }

        for (i = 0; i < lnn; i++)
        {
            int step = pow(2, i + 1);
#pragma omp parallel shared(aDiag, cDiag, bDiag, y, nn) private(j, index1, index2, alpha, gamma)
            {
#pragma omp for
                for (j = pow(2, i + 1) - 1; j < nn; j = j + step)
                {
                    index1 = j - pow(2, i);
                    index2 = j + pow(2, i);

                    alpha = aDiag[j] / bDiag[index1];
                    gamma = cDiag[j] / bDiag[index2];

                    aDiag[j] = -aDiag[index1] * (alpha);
                    bDiag[j] = bDiag[j] - cDiag[index1] * alpha - aDiag[index2] * gamma;
                    cDiag[j] = -cDiag[index2] * (gamma);

                    y[j] = y[j] - y[index1] * alpha - y[index2] * gamma;
                }
            }
        }

        int index = (nn - 1) / 2;
        x[index] = y[index] / bDiag[index];

        for (i = log2(nn + 1) - 2; i >= 0; i--)
        {
            int step = pow(2, i + 1);
#pragma omp parallel shared(x, y, aDiag, cDiag, bDiag, nn) private(j, index1, index2, alpha, gamma)

            {
#pragma omp for
                for (j = pow(2, i + 1) - 1; j < nn; j = j + step)
                {
                    offset = pow(2, i);
                    index1 = j - offset;
                    index2 = j + offset;

                    if (index1 - offset < 0)
                    {

                        x[index1] = (y[index1] - cDiag[index1] * x[index1 + offset]) / bDiag[index1];
                    }
                    else
                    {
                        x[index1] = (y[index1] - aDiag[index1] * x[index1 - offset] - cDiag[index1] * x[index1 + offset]) / bDiag[index1];
                    }

                    if (index2 + offset >= nn)
                    {
                        x[index2] = (y[index2] - aDiag[index2] * x[index2 - offset]) / bDiag[index2];
                    }
                    else
                    {
                        x[index2] = (y[index2] - aDiag[index2] * x[index2 - offset] - cDiag[index2] * x[index2 + offset]) / bDiag[index2];
                    }
                }
            }
        }
    }

    void solve_serial()
    {
        int j;
        double bet;

        double gam[n];

        x[0] = y[0] / (bet = bDiag[0]);

        for (j = 1; j < n; j++)
        {
            gam[j] = cDiag[j] / bet;
            bet = bDiag[j] - aDiag[j] * gam[j];
            x[j] = (y[j] - aDiag[j] * x[j - 1]) / bet;
        }
        for (j = (n - 2); j >= 0; j--)
            x[j] -= gam[j + 1] * x[j + 1];
    }

    void sety(std::vector<double> ynew)
    {
        for (int i = 0; i < n; i++)
            y[i] = ynew[i];
        for (int i = n; i < nn; i++)
            y[i] = 1.0;
    }

    std::vector<double> getx()
    {
        std::vector<double> xout;
        for (int i = 0; i < n; i++)
            xout.push_back(x[i]);

        return xout;
    }
    ~TriDiagSolver()
    {
        if (x)
            delete[] x;
        if (y)
            delete[] y;
        if (aDiag)
            delete[] aDiag;
        if (bDiag)
            delete[] bDiag;
        if (cDiag)
            delete[] cDiag;
        if (aDiagCpy)
            delete[] aDiagCpy;
        if (bDiagCpy)
            delete[] bDiagCpy;
        if (cDiagCpy)
            delete[] cDiagCpy;
    }
};

#endif