#include <vector>

float const G = 9.80655;
float const R = 287.0528;

float const h_0 = 0;
float const h_1 = 11000;
float const h_2 = 20000;
float const h_3 = 32000;
float const h_4 = 47000;
std::vector<float> const h_vec={h_0, h_1, h_2, h_3, h_4};

float const r_0 = 0.0065;
float const r_1 = 0;
float const r_2 = -0.001;
float const r_3 = -0.0028;
std::vector<float> const r_vec={r_0, r_1, r_2, r_3};

float const T_0 = 288.15;
float const T_1 =216.8;
float const T_2 = 216.7;
float const T_3 = 228.5;
std::vector<float> const T_vec={T_0, T_1, T_2, T_3};

float const p_0 = 101325;
float const p_1 = 22700;
float const p_2 = 5529;
float const p_3 = 889;
std::vector<float> p_vec={p_0, p_1, p_2, p_3};
