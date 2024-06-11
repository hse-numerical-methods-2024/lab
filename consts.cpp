#include <cmath>
#include <cfloat>
#include <climits>
#include <type_traits>

namespace ADAAI{
    template <typename F>
    constexpr inline F LOG2E;
    //specialization
    template <>
    constexpr inline float LOG2E <float> = (float)(M_LOG2E);
    template <>
    constexpr inline double LOG2E <double> = (double)(M_LOG2E);
    template <>
    constexpr inline long double LOG2E <long double> = (long double)(M_LOG2E);

    template <typename F>
    constexpr inline F SQRT2;
    //specialization
    template <>
    constexpr inline float SQRT2 <float> = (float)(M_SQRT2);
    template <>
    constexpr inline double SQRT2 <double> = (double)(M_SQRT2);
    template <>
    constexpr inline long double SQRT2 <long double> = (long double)(M_SQRT2);

    //TODO
}