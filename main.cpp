#include <iostream>
#include <cmath>

const int n = 100;

const double Smax = 100;
const double K = 80;
const double rTau = 1.0;
const double sigmaTau = 1.5;
const double sigmaTau2 = sigmaTau * sigmaTau;

const double h = Smax / n;

const double T = 20.0; // Максимальное время
const double tau = 0.001; // Шаг интегрирования по времени
// Граничное условие сверху
double phiU(double t, double S) {
	return Smax - K;
}
// Граничное условие снизу
double phiL(double t, double S) {
	return 0.0;
}

int main()
{
	double C[n + 1] = { 0.0 };
	double A[n + 1];
	double B[n + 1];
	double D[n + 1];
	double F[n + 1];
	double BS[n + 1];
	double FS[n + 1];

	double t = 0.0;
	C[0] = phiL(t, 0.0);
	C[n] = phiU(t, Smax);
	for (int i = 1; i < n; i++)
		C[i] = std::max(0.0, i*h - K);
	while (t < T) {
		// Неявный метод Эйлера
		for (int i = 0; i <= n; i++) {
			A[i] = 0.5*sigmaTau2*i*i - 0.5*rTau*i;
			B[i] = -rTau - sigmaTau2*i*i - 1/tau;
			D[i] = 0.5*sigmaTau2*i*i + 0.5*rTau*i;
			F[i] = -C[i]/tau;
		}

		BS[0] = 1.0;
		FS[0] = phiL(t, 0.0);
		for (int i = 1; i < n; i++) {
			BS[i] = B[i] - D[i-1]*A[i]/BS[i-1];
			FS[i] = F[i] - FS[i-1]*A[i]/BS[i-1];
		}

		C[n] = phiU(t, Smax);

		for (int i = n-1; i >= 0; i--) {
			C[i] = (FS[i] - D[i]*C[i+1])/BS[i];
		}

		if (t - (int)t < tau) {
			std::cout << "t = " << t << " : C = [ ";
			for (int i = 0; i <= n; i++)
				std::cout << C[i] << " ";
			std::cout << std::endl << std::endl;
		}

		t += tau;
	}
}
