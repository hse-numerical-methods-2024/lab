#include <iostream>
#include <cstdlib>
#include <ctime>
#include <limits>
#include <valarray>

#include "exp.cpp"

using namespace ADAII;
template <typename T>
T getRandom(T min, T max) {
    T randomFraction = static_cast<T>(std::rand()) / RAND_MAX;
    return min + randomFraction * (max - min);
}

template <typename T>
void checkError(int c) {
    T tolerance_factor = static_cast<T>(600);
    T err;
    for (int i = 1; i <= 100; i++) {
        std::cout << "Test " << i << "/100 ";
        T x = getRandom(static_cast<T>(-0.35), static_cast<T>(0.35));
        T res = Exp<MethodE::Taylor>(x);
        if (c == 1) {
            T res = Exp<MethodE::Pade>(x);
        }
        T stdRes = std::exp(x);
        if (x <= static_cast<T>(0)) {
            err = std::abs(res - stdRes);
        } else {
            err = std::abs(res / stdRes - static_cast<T>(1.0));
        }
        T eps = tolerance_factor * std::numeric_limits<T>::epsilon();
        if (err <= eps) {
            std::cout << "success" << std::endl;
            std::cout << "x = " << x << " std::exp = " << stdRes << " Exp<T> = " << res << " error = " << err << std::endl;
        } else {
            std::cout << "failed" << std::endl;
            std::cout << "x = " << x << " std::exp = " << stdRes << " Exp<T> = " << res << " error = " << err << std::endl;
            return;
        }
        std::cout << std::endl;
    }
}

int main() {
    std::cout<< "float"<<std::endl;
    checkError<float>(0);
    std::cout<< "double"<<std::endl;
    checkError<double>(0);
    std::cout<< "long double"<<std::endl;
    checkError<long double>(0);

    std::cout<< "Pade\nfloat"<<std::endl;
    checkError<float>(1);
    std::cout<< "double"<<std::endl;
    checkError<double>(1);
    std::cout<< "long double"<<std::endl;
    checkError<long double>(1);
    return 0;
}