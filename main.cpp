#include <iostream>
#include <cmath>

#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>

#include "Consts.hpp"

// Многочлены Чебышева первого рода (рекуррентное соотношение)
double T(int n, double x) {
    if (n == 0)
        return 1.0;
    else if (n == 1)
        return x;
    else
        return 2.0 * x * T(n - 1, x) - T(n - 2, x);
}

int main()
{
    double A[ADAAi::N+1][ADAAi::N+1] = { 0.0 }; // Коэффициенты системы
    double B[ADAAi::N+1] = { 0.0 }; // Вектор свободных членов системы
    for (int k = 0; k < ADAAi::N; k++) {
        A[k][k] = -1.0;
        for (int n = k+1; n <= ADAAi::N; n++) {
            if ((n % 2) == 0) {
               if ((k % 2) == 0)
                  A[k][n] = 0.0;
               else if (k == 1)
                  A[1][n] = 2*n;
               else
                  A[k][n] = 2*n;
            } else
               if ((k % 2) != 0)
                  A[k][n] = 0.0;
               else if (k == 0)
                  A[0][n] = n;
               else
                  A[k][n] = 2*n;
        }
    }
    A[ADAAi::N][0] = 1.0;
    for (int n = 2; n <= ADAAi::N; n+=2)
        A[ADAAi::N][n] = -A[ADAAi::N][n-2];
    B[ADAAi::N] = 1.0;

    gsl_matrix_view m = gsl_matrix_view_array ((double *)A, ADAAi::N+1, ADAAi::N+1);
    gsl_vector_view b = gsl_vector_view_array (B, ADAAi::N+1);
    gsl_vector * x = gsl_vector_alloc (ADAAi::N+1);
    int s;
    gsl_permutation * p = gsl_permutation_alloc (ADAAi::N+1);
    gsl_linalg_LU_decomp (&m.matrix, p, &s);
    gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);

    printf("c = ");
    gsl_vector_fprintf (stdout, x, "%g");
    gsl_permutation_free (p);

    const int NP = 200;
    double maxERR = 0.0;
    for (int i = 0; i < NP; i++) {
        double xx = -1.0 + 2.0*i/(NP-1);
        double e = 0.0;
        for (int i = 0; i <= ADAAi::N; i++)
            e += x->data[i]*T(i, xx);
        double err = std::abs(1.0 - e/std::exp(xx));
        if (err > maxERR)
           maxERR = err;
    }

    std::cout << "Max Relative ERR = " << maxERR << std::endl;

    gsl_vector_free (x);
}

