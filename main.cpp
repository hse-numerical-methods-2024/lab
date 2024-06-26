#include <cmath>
//#include "AAD.cpp"

const double h = 0.0001;
const double n = 2.0;

enum class DiffMethod : int {
    Stencil3,
    Stencil3Extra,
    Stencil5,
    Stencil5Extra,
    FwdAAD
};

enum class WhichD : int {
    X,
    Y,
    XX,
    XY,
    YY
};

double F(double x, double y) {
//    some func
//    return func;
}

template<WhichD W, DiffMethod M, typename Callable>
double Differentiator(Callable const &F, double x, double y, double h1) {
    if constexpr (M == DiffMethod::Stencil3) {
        if (W == WhichD::X) {
            if (x >= 1.0 && x <= -1.0) h1 = h1 * fabs(x);
            double func_dif = (F(x + h1, y) - F(x - h1, y)) / 2.0 / h1;
            return func_dif;
        } else if (W == WhichD::Y) {
            if (y >= 1.0 && y <= -1.0) h1 = h1 * fabs(y);
            double func_dif = (F(x, y + h1) - F(x, y - h1)) / 2.0 / h1;
            return func_dif;
        } else if (W == WhichD::XX) {
            if (x >= 1.0 && x <= -1.0) h1 = h1 * fabs(x);
            double func_dif = (F(x + h1, y) - 2.0 * F(x, y) + F(x - h1, y)) / h1 / h1;
            return func_dif;
        } else if (W == WhichD::YY) {
            if (y >= 1.0 && y <= -1.0) h1 = h1 * fabs(y);
            double func_dif = (F(x, y + h1) - 2.0 * F(x, y) + F(x, y - h1)) / h1 / h1;
            return func_dif;
        } else if (W == WhichD::XY) {
            double h2 = h1;
            if (y >= 1.0 && y <= -1.0) h1 = h1 * fabs(y);
            if (x >= 1.0 && x <= -1.0) h2 = h1 * fabs(x);
            double func_dif =
                    (F(x + h2, y + h1) - F(x - h2, y + h1) - F(x + h2, y - h1) + F(x - h2, y - h1)) / 4.0 / h1 / h2;
            return func_dif;
        }
    } else if constexpr (M == DiffMethod::Stencil3Extra) {
        double d3_1 = Differentiator<W, DiffMethod::Stencil3>(F, x, y, h);
        double d3_2 = Differentiator<W, DiffMethod::Stencil3>(F, x, y, h / n);
        double func_dif = (n * n * d3_2 - d3_1) / (n * n - 1.0);
        return func_dif;
    } else if constexpr (M == DiffMethod::Stencil5) {
        if (W == WhichD::X) {
            if (x >= 1.0 && x <= -1.0) h1 = h1 * fabs(x);
            double func_dif = F(x - 2.0 * h1, y) / 12.0 - 2.0 / 3.0 * F(x - h1, y) + 2.0 / 3.0 * F(x + h1, y) -
                              F(x + 2.0 * h1, y) / 12.0;
            return func_dif;
        } else if (W == WhichD::Y) {
            if (y >= 1 && y <= -1) h1 = h1 * fabs(y);
            double func_dif = F(x, y - 2.0 * h1) / 12.0 - 2.0 / 3.0 * F(x, y - h1) + 2.0 / 3.0 * F(x, y + h1) -
                              F(x, y + 2.0 * h1) / 12.0;
            return func_dif;
        } else if (W == WhichD::XX) {
            if (x >= 1 && x <= -1) h1 = h1 * fabs(x);
            double func_dif = -F(x - 2.0 * h1, y) / 12.0 + 4.0 / 3.0 * F(x - h1, y) - 5.0 / 2.0 * F(x, y) +
                              4.0 / 3.0 * F(x + h1, y) - F(x + 2.0 * h1, y) / 12.0;
            return func_dif;
        } else if (W == WhichD::YY) {
            if (y >= 1 && y <= -1) h1 = h1 * fabs(y);
            double func_dif = -F(x, y - 2.0 * h1) / 12.0 + 4.0 / 3.0 * F(x, y - h1) - 5.0 / 2.0 * F(x, y) +
                              4.0 / 3.0 * F(x, y + h1) - F(x, y + 2.0 * h1) / 12.0;
            return func_dif;
        } else if (W == WhichD::XY) {
            double h2 = h1;
            if (y >= 1 && y <= -1) h1 = h1 * fabs(y);
            if (x >= 1 && x <= -1) h2 = h1 * fabs(x);
            double func_dif =
                    (F(x - 2.0 * h2, y - 2.0 * h1) / 12.0 - 2.0 / 3.0 * F(x - h2, y - 2.0 * h1) +
                     2.0 / 3.0 * F(x + h2, y - 2.0 * h1) - F(x + 2.0 * h2, y - 2.0 * h1) / 12.0) / 12.0 -
                    2.0 / 3.0 * (F(x - 2.0 * h2, y - h1) / 12.0 - 2.0 / 3.0 * F(x - h2, y - h1) + 2.0 / 3.0 *
                                                                                                  F(x + h2, y - h1) -
                                 F(x + 2.0 * h2, y - h1) / 12.0) +
                    2.0 / 3.0 * (F(x - 2.0 * h2, y + h1) / 12.0 - 2.0 / 3.0 * F(x - h2, y + h1) + 2.0 / 3.0 *
                                                                                                  F(x + h2, y + h1) -
                                 F(x + 2.0 * h2, y + h1) / 12.0) -
                    (F(x - 2.0 * h2, y + 2.0 * h1) / 12.0 - 2.0 / 3.0 * F(x - h2, y + 2.0 * h1) +
                     2.0 / 3.0 * F(x + h2, y + 2.0 * h1) - F(x + 2.0 * h2, y + 2.0 * h1) / 12.0) / 12.0;
            return func_dif;
        }
    } else if constexpr (M == DiffMethod::Stencil5Extra) {
        double d5_1 = Differentiator<W, DiffMethod::Stencil5>(F, x, y, h);
        double d5_2 = Differentiator<W, DiffMethod::Stencil5>(F, x, y, h / n);
        double func_dif = (n * n * n * n * d5_2 - d5_1) / (n * n * n * n - 1); //not sure about this
        return func_dif;
    } else if constexpr (M == DiffMethod::FwdAAD) {
        //func for AAD method
    }
    }
