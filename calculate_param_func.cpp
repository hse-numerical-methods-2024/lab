#include <assert.h>
#include <cmath>

#include "Consts.cpp"

int h_bouder(float h) {
    if (h >= h_0 && h < h_1)
        return 0;
    if (h >= h_1 && h < h_2)
        return 1;
    if (h >= h_2 && h < h_3)
        return 2;
    if (h >= h_3 && h < h_4)
        return 3;
    assert(false);
}

float tempreture(float h) {
    float T = T_vec[h_bouder(h)] - r_vec[h_bouder(h)] * (h - h_vec[h_bouder(h)]);
    return T;
}

float pressure(float h) {
    int layer = h_bouder(h);
    float P = 0;
    if (r_vec[layer] == 0) {
        P = p_vec[layer] * exp(-G * (h - h_vec[layer]) / (R * T_vec[layer]));
    } else {
        P = p_vec[layer] * exp(-G / (R * r_vec[layer])) * (1 - r_vec[layer]*(h-h_vec[layer])/T_vec[layer]);
    }
    return P;
}

float density(float h) {
    float ro = pressure(h) / (R * tempreture(h));
    return ro;
}