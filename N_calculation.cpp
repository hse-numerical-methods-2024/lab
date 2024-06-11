#include <cmath>
#include <cassert>

constexpr double constexpr_exp(double x) {
    double sum = 1.0;
    double term = 1.0;
    for (int i = 1; i < 100; ++i) {
        term *= x / i;
        sum += term;
        if (term < std::numeric_limits<double>::epsilon()) {
            break;
        }
    }
    return sum;
}

constexpr double constexpr_pow(double base, int exp) {
    double result = 1.0;
    for (int i = 0; i < exp; ++i) {
        result *= base;
    }
    return result;
}

namespace ADAII {
    template<typename F>
    constexpr inline int MKExpTaylorOrder() {
        const F SQRT2 = M_SQRT2;
        F factorial_for_g = 1;
        for (int g = 1; g < 1000; g++) {
            factorial_for_g *= g;
            if (constexpr_exp(SQRT2)* constexpr_pow(1,g) / factorial_for_g  < std::numeric_limits<F>::epsilon())
                return g - 1;
        }
        assert(false);
    }

}