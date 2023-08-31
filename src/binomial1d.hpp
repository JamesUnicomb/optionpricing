#ifndef BINOMIAL1D_HPP
#define BINOMIAL1D_HPP

#include <omp.h>

class BinomialMesh1D
{
private:
    const int n;
    double *v;

public:
    BinomialMesh1D(int n);
    BinomialMesh1D(int n, const double &a);
    double getValue(
        double &p, double &pn);
    double getValue(
        double &p, double &pn, int current_step, int num_up_movements);
    double getExerciseValue(
        int current_step, int num_up_movements);
    void parallelStencilTriangle(
        double *p, double *past_edge_points, int jstart, int m, int triangle_size, int level);
    void parallelStencilRhombus(
        double *p, double *past_edge_points, double *curr_edge_points, int jstart, int m, int width, int level);
    double calculate_serial();
    double calculate_parallel();
    void set_initial_condition();
    ~BinomialMesh1D();
};

BinomialMesh1D::BinomialMesh1D(int n) : n(n), v(n > 0 ? new double[n] : nullptr) {}

BinomialMesh1D::BinomialMesh1D(int n, const double &a) : n(n), v(n > 0 ? new double[n] : nullptr)
{
    for (int i = 0; i < n; i++)
        v[i] = a;
}

double BinomialMesh1D::getValue(
    double &p, double &pn)
{
    return 0.5 * (p + pn);
}

inline double BinomialMesh1D::getValue(
    double &p, double &pn, int current_step, int num_up_movements)
{
    return 0.5 * (p + pn);
}

inline double BinomialMesh1D::getExerciseValue(
    int current_step, int num_up_movements)
{
    return 0.0;
}

void BinomialMesh1D::parallelStencilTriangle(
    double *p, double *past_edge_points, int jstart, int m, int triangle_size, int level)
{
    int i, j;
    double value;

    for (i = 0; i < triangle_size - 1; i++)
    {
        for (j = 0; j < triangle_size - i - 1; j++)
        {
<<<<<<< HEAD
            value = getValue(p[jstart + j], p[jstart + j + 1]);
            p[jstart + j] = value;

            if (jstart != 0 && j == 0)
            {
                past_edge_points[jstart - m + i + 1] = value;
            }
=======
            printf("v[%02d] = 0.5 * (v[%02d] + v[%02d]) = 0.5 * (%f + %f) = %f\n", j, j, j + 1, v[j], v[j + 1], 0.5 * (v[j] + v[j + 1]));
            v[j] = 0.5 * (v[j] + v[j + 1]);
>>>>>>> dfbfbd09be07cd77661fe3b16671fae3ae68f0e3
        }
    }
}

void BinomialMesh1D::parallelStencilRhombus(
    double *p, double *past_edge_points, double *curr_edge_points, int jstart, int m, int width, int level)
{
<<<<<<< HEAD
    int i, j, row_length;
    double value;

    for (i = 0; i < m; i++)
    {
        row_length = std::min(i + 1, width);
        for (j = 0; j < row_length; j++)
        {
            if (j == i)
            {
                value = getValue(
                    p[jstart + (m - i - 1) + j],
                    past_edge_points[jstart + i]);
                p[jstart + (m - i - 1) + j] = value;
            }
            else
=======
#pragma omp parallel num_threads(4)
    {
        int level, blocksize;
        int num = omp_get_num_threads();
        int idx = omp_get_thread_num();

        for (level = n - 1, blocksize = level / num + 1; level > 0; level -= blocksize, blocksize = level / num + 1)
        {
            int i, j, k, h;
            for (i = 0; i < blocksize; i++)
>>>>>>> dfbfbd09be07cd77661fe3b16671fae3ae68f0e3
            {
                value = getValue(
                    p[jstart + (m - i - 1) + j],
                    p[jstart + (m - i - 1) + j + 1]);
                p[jstart + (m - i - 1) + j] = value;

                if (jstart != 0 && i == m - 1 && j == 0)
                {
<<<<<<< HEAD
                    curr_edge_points[jstart = m] = value;
                }
            }
        }
    }
    for (i = 0; i < width - 1; i++)
    {
        row_length = std::min(m - i - 1, width - i - 1);
        for (j = 0; j < row_length; j++)
        {
            value = getValue(
                p[jstart + j],
                p[jstart + j + 1]);
            p[jstart + j] = value;
            if (jstart != 0 && j == 0)
            {
                curr_edge_points[jstart - m + i + 1] = value;
            }
        }
    }
}

double BinomialMesh1D::calculate_serial()
{
    for (int i = 0; i < n - 1; i++)
    {
        for (int j = 0; j < n - i - 1; j++)
        {
            v[j] = getValue(v[j], v[j + 1]);
=======
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
>>>>>>> dfbfbd09be07cd77661fe3b16671fae3ae68f0e3
        }
    }
    return v[0];
}

double BinomialMesh1D::calculate_parallel()
{
    int numblocks = 12;
    int blocksize = n / numblocks;
    int edgeblocksize = n % blocksize;

<<<<<<< HEAD
    double *past_edge_points = new double[numblocks * blocksize + edgeblocksize];

    for (int i = 0; i < numblocks; i++)
    {
        parallelStencilTriangle(
            v,
            past_edge_points,
            i * blocksize,
            blocksize,
            blocksize,
            n);
    }
    if (edgeblocksize > 0)
    {
        parallelStencilTriangle(
            v,
            past_edge_points,
            numblocks * blocksize,
            blocksize,
            edgeblocksize,
            n);
    }
=======
>>>>>>> dfbfbd09be07cd77661fe3b16671fae3ae68f0e3
    return v[0];
}

void BinomialMesh1D::set_initial_condition()
{
    for (int j = 0; j < n; j++)
    {
        v[j] = (j < n / 2) ? 1.0 : 0.0;
    }
}

BinomialMesh1D::~BinomialMesh1D()
{
    if (v != nullptr)
    {
        delete[] (v);
    }
}

#endif