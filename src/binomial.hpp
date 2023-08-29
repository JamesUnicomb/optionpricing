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
            // printf("v[%02d][%02d] = 0.5 * (v[%02d][%02d] + v[%02d][%02d]) = 0.5 * (%f + %f) = %f\n", i - 1, j, i, j, i, j + 1, v[i][j], v[i][j + 1], 0.5 * (v[i][j] + v[i][j + 1]));
            v[i - 1][j] = 0.5 * (v[i][j] + v[i][j + 1]);
        }
    }
    return v[0][0];
}

template <class T>
double BinomialMesh<T>::calculate_parallel()
{
#pragma omp parallel
    {
        int level = n - 1;
        int num = omp_get_num_threads();
        int idx = omp_get_thread_num();

        int blocksize = level / num + 1;

        while (blocksize > 1)
        {
            int istart, jstart;

            jstart = idx * blocksize;
            istart = level;

            int ii, jj;

            for (int i = 0; i < blocksize - 1; i++)
            {
                ii = istart - i;
                for (int j = 0; j < blocksize - i; j++)
                {
                    jj = jstart + j;
                    if (jj < ii)
                    {
                        // printf("v[%02d][%02d] = 0.5 * (v[%02d][%02d] + v[%02d][%02d]) = 0.5 * (%f + %f) = %f\n", ii - 1, jj, ii, jj, ii, jj + 1, v[ii][jj], v[ii][jj + 1], 0.5 * (v[ii][jj] + v[ii][jj + 1]));
                        v[ii - 1][jj] = 0.5 * (v[ii][jj] + v[ii][jj + 1]);
                    }
                }
            }

#pragma omp barrier

            for (int i = 0; i < blocksize - 1; i++)
            {
                ii = istart - i;
                for (int j = blocksize - 1; j > blocksize - i - 1; j--)
                {
                    jj = jstart + j;
                    if (jj < ii)
                    {
                        // printf("v[%02d][%02d] = 0.5 * (v[%02d][%02d] + v[%02d][%02d]) = 0.5 * (%f + %f) = %f\n", ii - 1, jj, ii, jj, ii, jj + 1, v[ii][jj], v[ii][jj + 1], 0.5 * (v[ii][jj] + v[ii][jj + 1]));
                        v[ii - 1][jj] = 0.5 * (v[ii][jj] + v[ii][jj + 1]);
                    }
                }
            }

#pragma omp barrier

            level -= blocksize - 1;
            blocksize = level / num + 1;
        }

        if (idx == 0)
        {
            for (int i = level; i > 0; i--)
            {
                for (int j = i - 1; j >= 0; j--)
                {
                    // printf("v[%02d][%02d] = 0.5 * (v[%02d][%02d] + v[%02d][%02d]) = 0.5 * (%f + %f) = %f\n", i - 1, j, i, j, i, j + 1, v[i][j], v[i][j + 1], 0.5 * (v[i][j] + v[i][j + 1]));
                    v[i - 1][j] = 0.5 * (v[i][j] + v[i][j + 1]);
                }
            }
        }
    }

    v[0][0] = 0.5 * (v[1][0] + v[1][1]);

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