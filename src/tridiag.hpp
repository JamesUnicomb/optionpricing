#ifndef TRIDIAG_HPP
#define TRIDIAG_HPP

#include <cmath>
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

class TriDiagSolver
{
private:
    int n;
    double *x, *d, *y;
    double a, b, c;

public:
    TriDiagSolver(int n, double a, double b, double c)
        : n(n), a(a), b(b), c(c),
          x(n > 0 ? new double[n] : nullptr),
          y(n > 0 ? new double[n] : nullptr),
          d(n > 0 ? new double[n] : nullptr)
    {
    }
    void solve()
    {
        double alpha, beta, gamma;
        double an, bn, cn;
        double r = 1.0 / (b - a * c);

        d0 = y[0] / b;
        for (int i = 1; i < n; i++)
        {
            d[i] =
        }

        alpha = a / b;
        beta = c / b;
        gamma = c * beta;

        // an = -a * alpha;
        // bn = b - c * gamma;
        // cn = -c * gamma;

        int id = (n - 1) / 2;
        x[id] = y[id] / bn;

        for (int i = log2(n + 1) - 2; i >= 0; i--)
        {
#pragma omp parallel
            {
                int step = pow(2, i + 1);
                int offset, id1, id2;

                offset = pow(2, i);
#pragma omp for
                for (int j = pow(2, i + 1) - 1; j < n; j = j + step)
                {
                    id1 = j - offset;
                    id2 = j + offset;

                    if (id1 - offset < 0)
                    {
                        x[id1] = (y[id1] - cn * x[id1 + offset]) / bn;
                    }
                    else
                    {
                        x[id1] = (y[id1] - an * x[id1 - offset] - cn * x[id1 + offset]) / bn;
                    }

                    if (id2 + offset >= n)
                    {
                        x[id2] = (y[id2] - an * x[id2 - offset]) / bn;
                    }
                    else
                    {
                        x[id2] = (y[id2] - an * x[id2 - offset] - cn * x[id2 + offset]) / bn;
                    }
                }
            }
        }
    }

    void sety(std::vector<double> ynew)
    {
        for (int i = 0; i < n; i++)
            y[i] = ynew[i];
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
    }
};

#endif