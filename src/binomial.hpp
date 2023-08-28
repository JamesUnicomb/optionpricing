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

    T get(int i, int j)
    {
        return v[i][j];
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