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
            value = getValue(p[jstart + j], p[jstart + j + 1]);
            p[jstart + j] = value;

            if (jstart != 0 && j == 0)
            {
                past_edge_points[jstart - m + i + 1] = value;
            }
        }
    }
}

void BinomialMesh1D::parallelStencilRhombus(
    double *p, double *past_edge_points, double *curr_edge_points, int jstart, int m, int width, int level)
{
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
            {
                value = getValue(
                    p[jstart + (m - i - 1) + j],
                    p[jstart + (m - i - 1) + j + 1]);
                p[jstart + (m - i - 1) + j] = value;

                if (jstart != 0 && i == m - 1 && j == 0)
                {
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
        }
    }
    return v[0];
}

double BinomialMesh1D::calculate_parallel()
{
    int numblocks = 12;
    int blocksize = n / numblocks;
    int edgeblocksize = n % blocksize;

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