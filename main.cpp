#include <iostream>
#include <fstream>
#include <cmath>

const int k = 10; // Порядок метода Everhart
const int N = 20; // Число итераций на одном шаге метода
const int p = 3; // Размерность пространства

const double mu = 398600.4; // km3/s2
const double J2 = 1.0827E-3;
const double Re = 6378.137; // km
const double a = 7500.0; // km
const double v0 = std::sqrt(mu / a); // km/s

const double tau = 1.0; // Шаг по времени

const double T = 80000; // Конечное время

typedef void (*F)(double t, double y[], double ys[], double f[]);

// Рекурсивная функция вычисления коэффициентов B
void Calculate_B(int q, double B[p], double F[k+1][p], double hh, int s, int m, int v, double r) {
    if (s < v) {
        Calculate_B(q, B, F, hh, s + 1, m, v, r);
        Calculate_B(q, B, F, hh, s + 1, m - 1, v, -r * s * hh);
    } else
        B[m] += F[v][q] * r;
}

// Шаг метода Everhart
void EverhartStep(F f, double h, double t0, double y0[p], double ys0[p], double & t1, double y1[p], double ys1[p]) {
    double f0[p];
    f(t0, y0, ys0, f0);

    double hh = h / k;
    double Fi[k + 1][k + 1][p] = {0.0};
    double yi[k + 1][p] = {0.0};
    double ysi[k + 1][p] = {0.0};
    for (int i = 0; i < p; i++) {
        Fi[0][0][i] = f0[i];
        yi[0][i] = y0[i];
        ysi[0][i] = ys0[i];
    }
    for (int q = 0; q < p; q++)
        for (int i = 1; i <= k; i++) {
            yi[i][q] = yi[0][q] + ysi[0][q] * i * hh + Fi[0][0][q] * (i * hh * i * hh) / 2.0;
            ysi[i][q] = ysi[0][q] + Fi[0][0][q] * (i * hh);
        }
    for (int n = 0; n < N; n++) { // Итерации уточнения
        for (int i = 1; i <= k; i++) {
            f(t0 + i * hh, yi[i], ysi[i], Fi[i][0]);
        }
        for (int i = 1; i <= k; i++) {
            for (int m = 0; m <= k - i; m++)
                for (int q = 0; q < p; q++)
                    Fi[m][i][q] = (Fi[m + 1][i - 1][q] - Fi[m][i - 1][q]) / (i * hh);
        }
        double B[k + 1][p][k + 1] = {0.0};
        for (int m = 0; m <= k; m++) {
            for (int q = 0; q < p; q++) {
                B[m][q][0] = Fi[0][0][q];
                for (int i = 1; i <= k; i++)
                    Calculate_B(q, B[m][q], Fi[0], hh, 0, i, i, 1.0);
            }
        }
        for (int i = 1; i <= k; i++) {
            for (int q = 0; q < p; q++) {
                ysi[i][q] = ysi[0][q];
                yi[i][q] = yi[0][q] + ysi[0][q] * (i * hh);
                for (int j = 0; j <= k; j++) {
                    ysi[i][q] += B[i][q][j] * std::pow(i * hh, j + 1) / (j + 1);
                    yi[i][q] += B[i][q][j] * std::pow(i * hh, j + 2) / ((j + 1) * (j + 2));
                }
            }
        }
    }
    t1 = t0 + h;
    for (int i = 0; i < p; i++) {
        y1[i] = yi[k][i];
        ys1[i] = ysi[k][i];
    }
}

void f(double t, double yy[p], double ys[p], double ff[p]) { // Правая часть системы = grad(U)
    double x = yy[0];
    double y = yy[1];
    double z = yy[2];

    double r2 = x * x + y * y + z * z;
    double r = std::sqrt(r2);
    double r5 = std::pow(r, 5.0);
    double r7 = r5 * r2;

    ff[0] =
        mu * x * (2 * (4 * r2 + J2 * Re * Re) / r5 - 5 * (2 * r2 * r2 + J2 * Re * Re * (x * x + y * y - 2 * z * z)) / r7) / 2;
    ff[1] =
        mu * y * (2 * (4 * r2 + J2 * Re * Re) / r5 - 5 * (2 * r2 * r2 + J2 * Re * Re * (x * x + y * y - 2 * z * z)) / r7) / 2;
    ff[2] =
        mu * z * (4 * (2 * r2 - J2 * Re * Re) / r5 - 5 * (2 * r2 * r2 + J2 * Re * Re * (x * x + y * y - 2 * z * z)) / r7) / 2;
}

int main()
{
    double t1 = 0.0;
    double r[p] = { 0.0, 0.0, a };
    double rs[p] = { v0 };
    double r1[p] = { 0.0 };
    double rs1[p] = { 0.0 };

    std::ofstream out("sattelite.dat"); // Данные сохраняются в файл sattelite.dat

    while (t1 < T) {
        EverhartStep(f, tau, t1, r, rs, t1, r1, rs1);
        for (int i = 0; i < p; i++) {
            r[i] = r1[i];
            rs[i] = rs1[i];
        }
        if (t1 - 1000 * (int)(t1 / 1000) < 1.0) {
            std::cout << "Time: " << t1 << std::endl;
        }
        if (t1 - 50 * (int)(t1 / 50) < 1.0) {
            if (out)
                out << t1 << " " << r[0] << " " << r[1] << " " << r[2] << std::endl;
        }
    }

    out.close();

    return 0;
}
