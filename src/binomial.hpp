#ifndef BINOMIAL_HPP
#define BINOMIAL_HPP

#include <omp.h>

template <class T>
inline T SQR(const T a) { return a * a; }

template <class T>
class BinomialMesh
{
private:
    const int n;
    const int nel;
    T **v;

public:
    BinomialMesh(int n);
    BinomialMesh(int n, const T &a);
    T get(int i, int j);
    double calculate_serial();
    double calculate_parallel();
    void set_initial_condition();
    void print();
    ~BinomialMesh();
};

template <class T>
BinomialMesh<T>::BinomialMesh(int n) : n(n), nel(n * (n + 1) / 2), v(n > 0 ? new T *[n] : nullptr)
{
    int i;

    if (v)
        v[0] = nel > 0 ? new T[nel] : nullptr;

    for (i = 1; i < n; i++)
        v[i] = v[i - 1] + i;
}

template <class T>
BinomialMesh<T>::BinomialMesh(int n, const T &a) : n(n), nel(n * (n + 1) / 2), v(n > 0 ? new T *[n] : nullptr)
{
    int i, j;
    if (v)
        v[0] = nel > 0 ? new T[nel] : nullptr;

    for (i = 1; i < n; i++)
        v[i] = v[i - 1] + i;

    for (i = 0; i < n; i++)
        for (j = 0; j <= i; j++)
            v[i][j] = a;
}

template <class T>
T BinomialMesh<T>::get(int i, int j)
{
    return v[i][j];
}

template <class T>
double BinomialMesh<T>::calculate_serial()
{
    for (int i = n - 1; i > 0; i--)
    {
        for (int j = i - 1; j >= 0; j--)
        {
            v[i - 1][j] = 0.5 * (v[i][j] + v[i][j + 1]);
        }
    }
    return v[0][0];
}

template <class T>
double BinomialMesh<T>::calculate_parallel()
{
    int num = 2;
    int level = n - 1;
    int blocksize = level / num + 1;

    while (blocksize > 1)
    {
        for (int idx = 0; idx < level / blocksize + 1; idx++)
        {
            int jstart, jstop, jdiff;

            jstart = idx * blocksize;
            jstop = std::min((idx + 1) * blocksize, level);
            jdiff = jstop - jstart;

            int ii, jj;
            int step = 1;
            for (int i = 0; i < blocksize; i++)
            {
                ii = level - i;
                for (int j = 0; j < blocksize - i - 1; j++)
                {
                    jj = jstart + j;
                    printf("loop1: v[%d][%d] = 0.5 * (v[%d][%d] + v[%d][%d]) = 0.5 * (%f + %f)\n", ii - 1, jj, ii, jj, ii, jj + 1, v[ii][j], v[ii][jj + 1]);
                    v[ii - 1][jj] = 0.5 * (v[ii][jj] + v[ii][jj + 1]);
                }
                step++;
            }

            // #pragma omp barrier

            // printf("level = %d\n", level);

            // // step = blocksize - 1;
            // for (int i = level - 1; i > level - blocksize; i--)
            // {
            //     printf("i = %d, jstart, jstop = %d, %d\n", i, jstart, jstop);
            //     for (int j = jstop - 1; j > jstart; j--)
            //     {
            //         printf("loop2: v[%d][%d] = 0.5 * (v[%d][%d] + v[%d][%d]) = 0.5 * (%f + %f)\n", i - 1, j, i, j, i, j + 1, v[i][j], v[i][j + 1]);
            //         v[i - 1][j] = 0.5 * (v[i][j] + v[i][j + 1]);
            //     }
            //     step--;
            // }
        }

        level -= blocksize - 1;
        blocksize = level / num + 1;
    }

    return v[0][0];
}

template <class T>
void BinomialMesh<T>::set_initial_condition()
{
    for (int j = 0; j < n; j++)
    {
        v[n - 1][j] = (j < n / 2) ? 1.0 : 0.0;
    }
}

template <class T>
void BinomialMesh<T>::print()
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

template <class T>
BinomialMesh<T>::~BinomialMesh()
{
    if (v != nullptr)
    {
        delete[] (v[0]);
        delete[] (v);
    }
}

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