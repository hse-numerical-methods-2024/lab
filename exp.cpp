#include <cmath>
#include <climits>
#include <vector>

#include "N_calculation.cpp"

#include "pade_calculation.cpp"

namespace ADAII {
    enum class MethodE: int {
        Taylor,
        Pade
    };
    template<MethodE M = MethodE::Pade, typename F>

    [[maybe_unused]] constexpr F Exp(F x) {
        const F LOG2E = M_LOG2E;
        const F LN2 = M_LN2;
        const F SQRT2 = M_SQRT2;

        if constexpr (M == MethodE::Taylor) {
        if (x - 1 < INT_MIN)
            return 0.0;
        if (x + 1 > INT_MAX)
            return INFINITY;

        F y = x * LOG2E;
        F y0;
        F y1 = std::modf(y, &y0); //y = y0 + y1, y0 is int

        //making borders
        if (y1 > 0.5) {
            y0++;
            y1--;
        } else if (y1 < -0.5) {
            y1++;
            y0--;
        }

            F x1;
            F accuracy;
            F Tailor = 0.0;
            y1 /= LOG2E;
            int n = 0;
            int res = y0;
            long long int factorial = 1;

            //making accuracy for Tailor series
            if (std::abs(y1) > std::abs(y1 - (LN2 / 2))) {
                accuracy = SQRT2;
                x1 = LN2 / 2;
            } else if (std::abs(y1) > std::abs(y1 + (LN2 / 2))) {
                accuracy = SQRT2 / 2;
                x1 = -LN2 / 2;
            } else {
                accuracy = 1;
                x1 = 0;
            }

            constexpr int N = MKExpTaylorOrder<F>();
            //Tailor series
            for (int i = 0; i < N; i++) {
                Tailor += accuracy * pow((y1 - x1), n) / factorial;
                n++;
                factorial *= n;
            }

            return std::ldexp(Tailor, res);
        }

        else
            if constexpr (M == MethodE::Pade) {

                std::vector<F> quotient;
                std::vector<F> remainder;
                if (std::is_same<F, float>::value) {
                    std::vector<F> P = {1.0/120.0, 1.0/10.0, 1.0/2.0, 1.0};
                    std::vector<F> Q = {-1.0/120.0, 1.0/10.0, -1.0/2.0, 1.0};
                    F divide_answer = dividePolynomials<F>(P, Q, x, quotient, remainder);
                    return divide_answer;
                }
                else
                if (std::is_same<F, double>::value) {
                    std::vector<F> P = {1.0/332640.0, 1.0/9240.0, 1.0/528.0, 2.0/99.0, 3.0/22.0, 6.0/11.0, 1.0};
                    std::vector<F> Q = {-1.0/55440.0, 1.0/1584.0, -1.0/99.0, 1.0/11.0, -5.0/11.0, 1.0};
                    F divide_answer = dividePolynomials<F>(P, Q, x, quotient, remainder);
                    return divide_answer;
                }
                else
                if (std::is_same<F, long double>::value) {
                    std::vector<F> P = {1.0/8648640.0, 7.0/1235520.0, 7.0/51480.0, 7.0/3432.0, 35.0/1716.0, 7.0/52.0, 7.0/13.0, 1.0};
                    std::vector<F> Q = {1.0/1235520.0, -1.0/25740.0, 1.0/1144.0, -5.0/429.0, 5.0/52.0, -6.0/13.0, 1.0};
                    F divide_answer = dividePolynomials<F>(P, Q, x, quotient, remainder);
                    return divide_answer;
                }
            }
    }
}