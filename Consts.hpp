#pragma once

#define _USE_MATH_DEFINES

#include <type_traits>
#include <cmath>
#include <limits>

namespace ADAAi {
    template<typename F>
    constexpr inline F Ln2;

    template<>
    constexpr inline float Ln2<float> = 1.0f/M_LOG2Ef;

    template<>
    constexpr inline double Ln2<double> = 1.0/M_LOG2E;

    template<>
    constexpr inline long double Ln2<long double> = 1.0L/M_LOG2El;

    template<typename F>
    constexpr inline F Sqrt2;

    template<>
    constexpr inline float Sqrt2<float> = M_SQRT2f;

    template<>
    constexpr inline double Sqrt2<double> = M_SQRT2;

    template<>
    constexpr inline long double Sqrt2<long double> = M_SQRT2l;

    template<typename F>
    constexpr inline F Eps = std::numeric_limits<F>::epsilon();
}