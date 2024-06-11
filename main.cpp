#include <iostream>
#include <cmath>

const int n = 100;

const double Smax = 100;
const double K = 80;
const double rTau = 1.0;
const double sigmaTau = 1.5;
const double sigmaTau2 = sigmaTau * sigmaTau;

const double h = Smax / n;

const double T = 10.0; // Максимальное время
const double tau = 0.0000001; // Шаг интегрирования по времени
// Граничное условие сверху
double phiU(double t, double S) {
	return Smax - K;
}
// Граничное условие снизу
double phiL(double t, double S) {
	return 0.0;
}
// Правые части системы
void rhs(double t, double C[n + 1], double F[n + 1]) {
	F[1] = rTau * h * (-25*C[1]+48*C[2]-36*C[3]+16*C[4]-3*C[5]) / (12 * h) +
		0.5*sigmaTau2*h*h*(35*C[1]-104*C[2]+114*C[3]-56*C[4]+11*C[5]) / (12 * h * h) - rTau*C[1];
	F[n-1] = rTau * (n-1) * h * (3 * C[n-5] - 16 * C[n-4] + 36 * C[n-3] - 48 * C[n-2] + 25 * C[n-1]) / (12 * h) +
		0.5 * sigmaTau2 * (n-1) * (n-1) * h * h * (11 * C[n-5] - 56 * C[n-4] + 114 * C[n-3] - 104 * C[n-2] + 35 * C[n-1]) / (12 * h * h) - rTau * C[n-1];
	for (int i = 2; i < n-1; i++)
		F[i] = rTau * i * h * (-C[i + 2] + 8 * C[i + 1] - 8 * C[i - 1] + C[i - 2]) / (12 * h) +
			0.5 * sigmaTau2 * i * i * h * h * (-C[i + 2] + 16 * C[i + 1] - 30 * C[i] + 16 * C[i - 1] - C[i - 2]) / (12 * h * h) - rTau * C[i];
}

int main()
{
	double C[n + 1] = { 0.0 };
	double C1[n + 1];

	double a[n + 1];
	double b[n + 1];
	double c[n + 1];
	double d[n + 1];

	double t = 0.0;
	C[0] = phiL(t, 0.0);
	C[n] = phiU(t, Smax);
	for (int i = 1; i < n; i++)
		C[i] = std::max(0.0, i*h - K);
	while (t < T) {
		// Метод Рунге-Кутта-4
		rhs(t, C, a);
		for (int i = 1; i < n; i++) {
			a[i] *= tau;
			C1[i] = C[i] + 0.5 * a[i];
		}
		C1[0] = phiL(t, 0.0);
		C1[n] = phiU(t, Smax);
		rhs(t + 0.5 * tau, C1, b);
		for (int i = 1; i < n; i++) {
			b[i] *= tau;
			C1[i] = C[i] + 0.5 * b[i];
		}
		C1[0] = phiL(t, 0.0);
		C1[n] = phiU(t, Smax);
		rhs(t + 0.5 * tau, C1, c);
		for (int i = 1; i < n; i++) {
			c[i] *= tau;
			C1[i] = C[i] + c[i];
		}
		C1[0] = phiL(t, 0.0);
		C1[n] = phiU(t, Smax);
		rhs(t + tau, C1, d);
		for (int i = 1; i < n; i++) {
			d[i] *= tau;
		}
		for (int i = 1; i < n; i++) {
			C[i] += (a[i] + 2*b[i] + 2*c[i] + d[i]) / 6.0;
		}
		C[0] = phiL(t, 0.0);
		C[n] = phiU(t, Smax);

		if (t - (int)t < tau) {
			std::cout << "t = " << t << " : C = [ ";
			for (int i = 0; i <= n; i++)
				std::cout << C[i] << " ";
			std::cout << std::endl << std::endl;
		}

		t += tau;
	}
}
