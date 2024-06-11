#pragma once

#define _USE_MATH_DEFINES

#include <type_traits>
#include <cmath>
#include <limits>

#include "Consts.hpp"

namespace ADAAi {
    template<typename F>
    constexpr F exp(F x) {
        static_assert(std::is_floating_point_v<F>);
        F y00;
        F y01 = std::modf(x/Ln2<F>, &y00);
        if (y00 < std::numeric_limits<int>::min())
            return 0.0;
        else if (y00 > std::numeric_limits<int>::max())
            return std::numeric_limits<F>::infinity();
        int n = int(y00);
        F y1 = y01;
        if (y01 <= (F) -0.5L) {
            y1 = (F)1.0L + y01;
            n--;
        }
        else if (y01 >= (F)0.5L) {
            y1 = y01 - (F)1.0L;
            n++;
        }

        F x1 = y1 * Ln2<F>;
        F f1 = (F)0.0L;
        F delta = 10.0L * Eps<F>;
        F C = x1 > 0.0L ? Sqrt2<F> : (F)1.0L;

        F addition_p = (F)1.0L;
        F addition_q = (F)1.0L;
        int k = 1;
     
        while (C*std::abs(addition_p/addition_q) >= delta) {
            f1 += addition_p/addition_q;
            addition_p *= x1;
            addition_q *= k;
            k++;
        }

        return std::ldexp(f1, n);
    }
}