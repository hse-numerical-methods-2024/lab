#include "C_d_consts.cpp"

float C_d(float M) {
    int i = 1;
    while (M > C_d_value[i].first) {
        i++;
        if (i == 16) {
            assert(false);
        }
    }
    float C_d = C_d_value[i - 1].second + (C_d_value[i].second - C_d_value[i - 1].second) *
                                          (M - C_d_value[i - 1].first) / (C_d_value[i].first - C_d_value[i - 1].first);
    return C_d;
}