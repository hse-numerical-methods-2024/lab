#include <iostream>
#include "calculate_param_func.cpp"
#include "C_d_calc.cpp"

#define pi M_PI

float const D = 0.216; // calibre of a Paris gun
float const Mass = 106; // weight


float Q_calc(float h, float v_in_2, float d) {
    float S = pi * d * d / 4;
    float M = sqrt(v_in_2) / sqrt(pressure(h) / density(h));
    float Q = C_d(M) * S * v_in_2 / 2 * density(h);
    return Q;
}


//void test(std::vector<float> vec){
//    std::vector<float> cmp(4); //
//    cmp[0] = 940.665;
//    cmp[1] = -0.3;
//    cmp[2] = 699.335;
//    cmp[3] = -10.03;
//    for (int i = 0; i < 4; i++)
//    std::cout << abs(cmp[i] - vec[i]) << "\n";
//}

int main() {
    std::vector<float> u(4);
    std::vector<float> u_der(4);

//    std::cin >> u[0] >> u[1] >> u[2] >> u[3];
//    u[0] = 1;
//    u[1] = 940.665;
//    u[2] = 1.428;
//    u[3] = 699.335;

    float v_in_2 = u[1] * u[1] + u[3] * u[3];
    float Q = Q_calc(u[2], v_in_2, D);

    u_der[0] = u[1];
    u_der[2] = u[3];
    u_der[1] = -Q * u[1] / sqrt(v_in_2) / Mass;
    u_der[3] = -Q * u[3] / sqrt(v_in_2) / Mass - G;
//    test(u_der);

//    std::cout << u_der[0] << std::endl << u_der[1] << std::endl << u_der[2] << std::endl << u_der[3] << std::endl;
    return 0;
}
