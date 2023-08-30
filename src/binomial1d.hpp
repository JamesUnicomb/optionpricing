#ifndef BINOMIAL1D_HPP
#define BINOMIAL1D_HPP

#include <omp.h>

template <class T>
class BinomialMesh1D
{
private:
    const int n;
    T *v;

public:
    BinomialMesh1D(int n);
    BinomialMesh1D(int n, const T &a);
    double calculate_serial();
    double calculate_parallel();
    void set_initial_condition();
    ~BinomialMesh1D();
};

template <class T>
BinomialMesh1D<T>::BinomialMesh1D(int n) : n(n), v(n > 0 ? new T[n] : nullptr) {}

template <class T>
BinomialMesh1D<T>::BinomialMesh1D(int n, const T &a) : n(n), v(n > 0 ? new T[n] : nullptr)
{
    for (int i = 0; i < n; i++)
        v[i] = a;
}

template <class T>
double BinomialMesh1D<T>::calculate_serial()
{
    for (int i = 0; i < n - 1; i++)
    {
        for (int j = 0; j < n - i - 1; j++)
        {
            printf("i = %02d, j = %02d : v[%02d] = 0.5 * (v[%02d] + v[%02d]) = 0.5 * (%f + %f) = %f\n", i, j, j, j, j + 1, v[j], v[j + 1], 0.5 * (v[j] + v[j + 1]));
            v[j] = 0.5 * (v[j] + v[j + 1]);
        }
    }
    return v[0];
}

template <class T>
double BinomialMesh1D<T>::calculate_parallel()
{
    // #pragma omp parallel
    //     {
    //         int level = 0;
    //         int num = omp_get_num_threads();
    //         int idx = omp_get_thread_num();

    //         int blocksize = n / num + 1;

    //         while (level < n)
    //         {
    //             int istart, jstart;

    //             jstart = idx * blocksize;
    //             istart = level;

    //             int ii, jj;

    //             for (int i = 0; i < blocksize - 1; i++)
    //             {
    //                 ii = istart - i;
    //                 for (int j = 0; j < blocksize - i; j++)
    //                 {
    //                     jj = jstart + j;
    //                     if (jj < ii)
    //                     {
    //                         printf("i = %02d, j = %02d : v[%02d] = 0.5 * (v[%02d] + v[%02d]) = 0.5 * (%f + %f) = %f\n", ii, jj, jj, jj, jj + 1, v[jj], v[jj + 1], 0.5 * (v[jj] + v[jj + 1]));
    //                         v[jj] = 0.5 * (v[jj] + v[jj + 1]);
    //                     }
    //                 }
    //             }

    // #pragma omp barrier

    //             for (int i = 0; i < blocksize - 1; i++)
    //             {
    //                 ii = istart - i;
    //                 for (int j = 0; j < blocksize - i - 1; j++)
    //                 {
    //                     jj = jstart + j;
    //                     if (jj < ii)
    //                     {
    //                         printf("i = %02d, j = %02d : v[%02d] = 0.5 * (v[%02d] + v[%02d]) = 0.5 * (%f + %f) = %f\n", ii, jj, jj, jj, jj + 1, v[jj], v[jj + 1], 0.5 * (v[jj] + v[jj + 1]));
    //                         v[jj] = 0.5 * (v[jj] + v[jj + 1]);
    //                     }
    //                 }
    //             }

    // #pragma omp barrier

    //             level += blocksize - 1;
    //             blocksize = (n - level) / num + 1;
    //         }

    //         // if (idx == 0)
    //         // {
    //         //     for (int i = 0; i < level; i++)
    //         //     {
    //         //         for (int j = 0; j < n - i - 1; j++)
    //         //         {
    //         //             printf("i = %02d, j = %02d : v[%02d] = 0.5 * (v[%02d] + v[%02d]) = 0.5 * (%f + %f) = %f\n", i, j, j, j, j + 1, v[j], v[j + 1], 0.5 * (v[j] + v[j + 1]));
    //         //             v[j] = 0.5 * (v[j] + v[j + 1]);
    //         //         }
    //         //     }
    //         // }
    //     }

    return v[0];
}

template <class T>
void BinomialMesh1D<T>::set_initial_condition()
{
    for (int j = 0; j < n; j++)
    {
        v[j] = (j < n / 2) ? 1.0 : 0.0;
    }
}

template <class T>
BinomialMesh1D<T>::~BinomialMesh1D()
{
    if (v != nullptr)
    {
        delete[] (v);
    }
}

#endif