#include <iostream>
#include <iomanip>
#include <math.h>
#include <omp.h>
#include <time.h>

using namespace std;

int n, sizeN, nn;
int maxMatsize = 25;

int main()
{

    int count = 0;
    bool flag = false;
    cout << "Enter the nn of the system : ";
    cin >> n;

    double *p = new double[maxMatsize];

    for (int k = 0; k < maxMatsize; k++)
    {

        p[k] = pow(2, k) - 1;

        if (p[k] == n)
        {
            nn = p[k];
            flag = true;

            cout << "sizeN = " << sizeN;
            cout << endl;
        }
        else
        {
            if (p[k] >= n)
            {
                count = count + 1;
            }
        }
    }

    sizeN = p[maxMatsize - count];

    cout << nn << endl;
    cout << sizeN << endl;
    cout << n << endl;
    for (int i = 0; i < maxMatsize; i++)
    {
        cout << p[i] << " ";
    }
    cout << endl;

    if (flag == false)
    {
        nn = sizeN;
    }
    cout << nn << "  " << n;
    cout << endl;

    clock_t start = clock();

    int i, j;
    int index1, index2, offset;
    double alpha, gamma;

    double *x = new double[nn];
    for (i = 0; i < nn; i++)
    {
        x[i] = 0.0;
    }

    double *y = new double[nn];
    double *bDiag = new double[nn];
    double *aDiag = new double[nn];
    double *cDiag = new double[nn];

    if (flag == false)
    {
        // #pragma omp parallel for
        for (i = 0; i < n; i++)
        {

            y[i] = 1.0;
            bDiag[i] = 2.0;
            aDiag[i] = -1.0;
            cDiag[i] = -1.0;
        }

        // #pragma omp parallel for
        for (i = n; i < nn; i++)
        {

            y[i] = 1.0;
            bDiag[i] = 1.0;
            aDiag[i] = 0.0;
            cDiag[i] = 0.0;
        }
    }

    else
    {

        for (i = 0; i < nn; i++)
        {

            y[i] = 1.0;
            bDiag[i] = 2.0;
            aDiag[i] = -1.0;
            cDiag[i] = -1.0;
        }
    }

    aDiag[0] = 0.0;
    cDiag[n - 1] = 0.0;

    int lnn = log2(nn + 1) - 1;

    /// Cyclic Reduction Step 1
    for (i = 0; i < lnn; i++)
    {
        int step = pow(2, i + 1);
#pragma omp parallel shared(aDiag, cDiag, bDiag, y, nn) private(j, index1, index2, alpha, gamma)
        {
#pragma omp for
            for (j = pow(2, i + 1) - 1; j < nn; j = j + step)
            {

                // offset = pow(2,i);
                index1 = j - pow(2, i);
                index2 = j + pow(2, i);

                alpha = aDiag[j] / bDiag[index1];
                gamma = cDiag[j] / bDiag[index2];

                // #pragma omp atomic capture
                aDiag[j] = -aDiag[index1] * (alpha);
                bDiag[j] = bDiag[j] - cDiag[index1] * alpha - aDiag[index2] * gamma;
                cDiag[j] = -cDiag[index2] * (gamma);
                y[j] = y[j] - y[index1] * alpha - y[index2] * gamma;
            }
        }
    }

    for (i = 0; i < nn; i++)
    {
        std::cout << i << "," << aDiag[i] << "," << bDiag[i] << "," << cDiag[i] << "," << x[i] << "," << y[i] << std::endl;
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

                // printf("Executed by %d \n", omp_get_thread_num());
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

    for (i = 0; i < n; i++)
    {
        cout << i << " = " << x[i] << endl;
    }

    // Stop measuring time and calculate the elapsed time
    clock_t end = clock();
    double elapsed = double(end - start) / CLOCKS_PER_SEC;

    //    printf("Time measured: %.3f seconds.\n", elapsed);
    cout << "The Time measured is : " << elapsed << endl;
    return 0; // return 0 to the OS.
}