#ifndef BINOMIAL_HPP
#define BINOMIAL_HPP

#include <omp.h>

template <class T>
class BinomialMesh
{
private:
    int n;
    T **v;

public:
    BinomialMesh(int n) : n(n), v(n > 0 ? new T *[n + 1] : nullptr)
    {
        int i, nel = n * (n + 1) / 2;

        if (v)
            v[0] = nel > 0 ? new T[nel] : nullptr;

        for (i = 1; i <= n; i++)
            v[i] = v[i - 1] + i;
    }
    BinomialMesh(int n, const T &a) : n(n), v(n > 0 ? new T *[n + 1] : nullptr)
    {
        int i, j, nel = n * (n + 1) / 2;
        if (v)
            v[0] = nel > 0 ? new T[nel] : nullptr;

        for (i = 1; i < n; i++)
            v[i] = v[i - 1] + i;

        for (i = 0; i < n; i++)
            for (j = 0; j <= i; j++)
                v[i][j] = a;
    }

    double calculate_serial()
    {
        for (int i = n - 1; i >= 0; i--)
        {
            for (int j = 0; j < i; j++)
            {
                v[i - 1][j] = 0.5 * (v[i][j] + v[i][j + 1]);
            }
        }
        return v[0][0];
    }

    double calculate_parallel1()
    {
        for (int i = n - 1; i >= 0; i--)
        {
#pragma omp parallel for
            for (int j = 0; j < i; j++)
            {
                v[i - 1][j] = 0.5 * (v[i][j] + v[i][j + 1]);
            }
        }
        return v[0][0];
    }

    double calculate_parallel2()
    {
        for (int i = n - 1; i >= 0; i--)
        {
#pragma omp parallel for
            for (int j = 0; j < i; j = j + 2)
            {
                v[i - 1][j] = 0.5 * (v[i][j] + v[i][j + 1]);
            }
#pragma omp parallel for
            for (int j = 1; j < i; j = j + 2)
            {
                v[i - 1][j] = 0.5 * (v[i][j] + v[i][j + 1]);
            }
        }
        return v[0][0];
    }

    void set_initial_condition()
    {
        for (int i = 0; i <= n; i++)
        {
            v[n - 1][i] = (i < n / 2) ? 1.0 : 0.0;
        }
    }

    void print()
    {
        for (int i = 0; i < n; i++)
        {
            printf("%d ", i);
            for (int j = 0; j <= i; j++)
            {
                printf("%f ", v[i][j]);
            }
            printf("\n");
        }
    }
};

class BinomialPricing
{
private:
    const int n;

public:
    BinomialPricing(int n) : n(n) {}
    int sum_parallel()
    {
        int sum = 0;
#pragma omp parallel for reduction(+ : sum)
        for (int i = 0; i <= n; ++i)
        {
            sum += i;
        }
        return sum;
    }

    int sum_serial()
    {
        int sum = 0;
        for (int i = 0; i <= n; ++i)
        {
            sum += i;
        }
        return sum;
    }
};

#endif