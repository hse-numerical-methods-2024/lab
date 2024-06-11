#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

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
    double a_an[ADAAi::N + 1]; // Аналитическое решение
    for (int k = 0; k <= ADAAi::N; k++)
        a_an[k] = 2.0 * (std::exp(M_PI) * k * std::sin(M_PI * k) + std::exp(M_PI) * std::cos(M_PI * k) - 1.0) / (M_PI * (k * k + 1));
    std::cout << "Analytically: ";
    for (int k = 0; k <= ADAAi::N; k++)
        std::cout << a_an[k] << " ";
    std::cout << std::endl;
    // Решение через квадратуру Гаусса-ЧЕбышева
    double a_gc[ADAAi::N + 1] = { 0.0 };
    double xi[ADAAi::N + 2];
    for (int i = 1; i <= ADAAi::N + 1; i++)
        xi[i] = std::cos(M_PI*(2*i-1)/(2*(ADAAi::N+1)));
    for (int k = 0; k <= ADAAi::N; k++) {
        for (int i = 1; i <= ADAAi::N + 1; i++)
            a_gc[k] += std::exp(std::acos(xi[i])) * T(k, xi[i]);
        a_gc[k] *= 2.0 / (ADAAi::N + 1);
    }
    std::cout << "Gauss-Chebyshev: ";
    for (int k = 0; k <= ADAAi::N; k++)
        std::cout << a_gc[k] << " ";
    std::cout << std::endl;

    double x[(ADAAi::N+1)*2] = { 0.0 };
    gsl_complex_packed_array data = x;
    gsl_fft_complex_wavetable * wavetable = gsl_fft_complex_wavetable_alloc (ADAAi::N+1);
    gsl_fft_complex_workspace * workspace = gsl_fft_complex_workspace_alloc (ADAAi::N+1);
    data[0] = 0.5*a_gc[0];
    for (int k = 1; k <= ADAAi::N; k++) {
        data[k*2] = a_gc[k];
    }
    gsl_fft_complex_forward(data, 1, ADAAi::N+1, wavetable, workspace);
    std::cout << "Max relative EPS<approx with koeffs calculated by Gauss=Chebyshev>:" << std::endl;
    double maxEPS = 0.0;
    for (int k = 0; k <= ADAAi::N/2; k++) {
        double xxi = M_PI*k/(ADAAi::N+1);
        double EPS = std::abs(1.0 -  data[k*2]/std::exp(2.0*xxi));
        if (EPS > maxEPS)
           maxEPS = EPS;
    }
    std::cout << 100*maxEPS << "%" << std::endl;

    gsl_fft_complex_wavetable_free(wavetable);
    gsl_fft_complex_workspace_free(workspace);
    return 0;
}

