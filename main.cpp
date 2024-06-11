#define _USE_MATH_DEFINES

#include <type_traits>
#include <cmath>
#include <limits>

#include <iomanip>
#include <iostream>

#include "Consts.hpp"
#include "Exp.hpp"

template<typename F>
constexpr void experiment(F & maxEPS)
{
    constexpr int NP = 1000;
    for (int i = 0; i < NP; i++) {
        F X = -50.0 + 100.0*i/(NP-1);
        F EPS = std::abs((1.0 - std::abs(ADAAi::exp(X)/std::exp(X))) / ADAAi::Eps<F>);
        if (EPS > maxEPS)
           maxEPS = EPS;
    }
}

int main() {
    double EPS = 0.0;
    experiment(EPS);
    long double EPSL = 0.0L;
    experiment(EPSL);
    float EPSF = 0.0f;
    experiment(EPSF);
    std::cout << std::setprecision(22) << "float<EPS=" << ADAAi::Eps<float> << ">: " << EPSF << std::endl;
    std::cout << std::setprecision(22) << "double<EPS=" << ADAAi::Eps<double> << ">: " << EPS << std::endl;
    std::cout << std::setprecision(22) << "long double<EPS=" << ADAAi::Eps<long double> << ">: " << EPSL << std::endl;
}

