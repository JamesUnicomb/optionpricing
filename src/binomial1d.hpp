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
            printf("v[%02d] = 0.5 * (v[%02d] + v[%02d]) = 0.5 * (%f + %f) = %f\n", j, j, j + 1, v[j], v[j + 1], 0.5 * (v[j] + v[j + 1]));
            v[j] = 0.5 * (v[j] + v[j + 1]);
        }
    }
    return v[0];
}

template <class T>
double BinomialMesh1D<T>::calculate_parallel()
{
#pragma omp parallel num_threads(4)
    {
        int level, blocksize;
        int num = omp_get_num_threads();
        int idx = omp_get_thread_num();

        for (level = n - 1, blocksize = level / num + 1; level > 0; level -= blocksize, blocksize = level / num + 1)
        {
            int i, j, k, h;
            for (i = 0; i < blocksize; i++)
            {
                for (j = 0; j < blocksize - i; j++)
                {
                    h = level - i;
                    k = j + idx * blocksize;

                    if (k < h && k >= 0)
                    {
                        printf("v[%02d] = 0.5 * (v[%02d] + v[%02d]) = 0.5 * (%f + %f) = %f\n", k, k, k + 1, v[k], v[k + 1], 0.5 * (v[k] + v[k + 1]));
                        // printf("set 1 level = %02d, h = %02d, k = %02d, idx = %02d, i = %02d, j = %02d, idx = %02d, blocksize = %02d\n", level, h, k, idx, i, j, idx, blocksize);
                        v[k] = 0.5 * (v[k] + v[k + 1]);
                    }
                }
            }

#pragma omp barrier
            for (i = 1; i < blocksize; i++)
            {
                for (j = 0; j < i; j++)
                {
                    h = level - i;
                    k = j + idx * blocksize + (blocksize - i);

                    if (k < h && k >= 0)
                    {
                        printf("v[%02d] = 0.5 * (v[%02d] + v[%02d]) = 0.5 * (%f + %f) = %f\n", k, k, k + 1, v[k], v[k + 1], 0.5 * (v[k] + v[k + 1]));
                        // printf("set 2 level = %02d, h = %02d, k = %02d, idx = %02d, i = %02d, j = %02d, idx = %02d, blocksize = %02d\n", level, h, k, idx, i, j, idx, blocksize);
                        v[k] = 0.5 * (v[k] + v[k + 1]);
                    }
                }
            }
#pragma omp barrier
        }
    }

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