/**
 * Runge-Kutta 4th order
 * Equation -- Y + sqrt(pow(X,2) + pow(X,2)) - XY' = 0
 * Cauchy -- Y(Xo) = -0.5
 * Limits -- [0.0 ; 1.0]
 * Solution -- Y(X) = (pow(x,2) - 1) / 2
 */

/**
 * Possible testing equation system:
 * y1 = sinx;
 * y2 = cosx;
 * ----------
 * y1' = cosx = y2
 * y2' = -sinx = -y1;
 */

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <iomanip>

// Prototypes
double inline solution(int neq, double x) { return (pow(x, 2) - 1) / 2; }
double inline derivative(int neq, double x, double *y, double yp) { return (y[0] + sqrt(pow(x, 2) + pow(y[0], 2))) / x; }

void Output(double **yApproximate, double **yAnalytical, double *xk, int neq, int nOutIter);

// Method Function
double RungeKutta(int neq, double Xk, double Xk1, double h, double *y, double yp)
{
	// coefficients for Runge Kutta 4th order
	double *k1 = new double[neq];
	double *k2 = new double[neq];
	double *k3 = new double[neq];
	double *k0 = new double[neq];

	double x;

	double finalValue;

	int numOfIterations = (int)((Xk1 - Xk) / h); // Xk1 - Xk >> h --> nt is quite big
	if (numOfIterations == 0) numOfIterations = 1; // check to avoid division by zero
	h = (Xk1 - Xk) / numOfIterations; // make h more precise

	for (int it = 1; it <= numOfIterations; it++)
	{
		double a = Xk + (it - 1) * h; // left side (reserved, unlike x)

		for (int i = 0; i < neq; i++)
		{
			k0[i] = y[i]; // make copies of an array of initial values
		}

		yp = derivative(neq, a, y, yp); // calculate derivative Yp = F(Xk, Yk)

		for (int i = 0; i < neq; i++)
		{
			k1[i] = h * yp;
			y[i] = k0[i] + k1[i] / 2.0; // set next Yi for next Ki
		}

		x = a + h / 2.0; // set next Xi for next Ki
		yp = derivative(neq, x, y, yp);

		for (int i = 0; i < neq; i++)
		{
			k2[i] = h * yp;
			y[i] = k0[i] + k2[i] / 2.0;
		}

		x = a + h / 2.0;
		yp = derivative(neq, x, y, yp);

		for (int i = 0; i < neq; i++)
		{
			k3[i] = h * yp;
			y[i] = k0[i] + k3[i];
		}

		x = a + h;
		yp = derivative(neq, x, y, yp);

		for (int i = 0; i < neq; i++)
		{
			y[i] = k0[i] + 1.0 / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + yp * h); // k4 = yp * h
		}
	}
	return y[0];
}

int main(void)
{
	// limits of an interval
	double a = pow(10, -10);
	double b = 1.0;

	constexpr int neq = 1; // number of equations

	// integrational steps for different accuracy of calculations
	constexpr double h1 = 0.2;			// h1 = 1/5
	constexpr double h2 = h1 / 5.0;		// h2 = 1/25
	constexpr double h3 = h1 / 25.0;	// h3 = 1/125
	constexpr double h = h3;

	// output preparation
	double outputStep = 0.1; // bigger step for output purpose
	// correction (not a must-have)
	int nOutIter = (int)((b - a) / outputStep);  // number of output iterations
	outputStep = (b - a) / nOutIter;

	// arrays to store OUTPUT values
	// TODO remake the whole idea to use these arrays w/o using array y (and maybe yp)
	double **yApproximate = new double*[neq];
	double **yAnalytical = new double*[neq];
	for (int i = 0; i < neq; ++i)
	{
		yApproximate[i] = new double[nOutIter + 1];
		yAnalytical[i] = new double[nOutIter + 1];
	}
	// array to store arguments
	double *xk = new double[32];

	// array to store CALCULATED values
	double *y = new double[neq];

	// Cauchy condition
	y[0] = -0.5; // TODO get rid of when remaking the whole idea
	yApproximate[0][0] = -0.5;
	// function derivative
	double yp = 0;

	for (int it = 0; it < nOutIter; it++)
	{
		// [Xki; Xki+1] -- intervals we iterate over [from a to b]
		double Xk = a + outputStep * (it);
		double Xk1 = a + outputStep * (it + 1);
		xk[it] = Xk;

		for (int i = 0; i < neq; i++)
		{
			yApproximate[i][it + 1] = RungeKutta(neq, Xk, Xk1, h, y, yp);
			y[i] = yApproximate[0][it + 1];

			yAnalytical[i][it] = solution(neq, Xk);
		}

	}

	Output(yApproximate, yAnalytical, xk, neq, nOutIter);
	
	std::cin.get();
	for (int i = 0; i < neq; i++)
	{
		delete[] yApproximate[i];
		delete[] yAnalytical[i];
	}
	delete[] yApproximate;
	delete[] yAnalytical;
	delete[] xk;
	delete[] y;

	return 0;
}

void Output(double **yApproximate, double **yAnalytical, double *xk, int neq, int nOutIter)
{
	std::cout << std::setprecision(12);
	// header
	std::cout << " Xk              ";
	for (int i = 0; i < neq; i++)
	{
		std::cout << " | Y(Xk)           "
			<< " | Y(X)            " << " | E(x)            "
			<< " | 100*E(x)/Yk" << std::left << std::endl;
	}
	for (int i = 0; i < 19 + 76 * neq; i++) std::cout << "-"; // header line
	std::cout << std::endl;

	// body
	for (int j = 0; j < nOutIter; j++)
	{
		std::cout << " " << std::setw(16) << xk[j];
		for (int i = 0; i < neq; i++)
		{
			std::cout << " | " << std::setw(16) << yApproximate[i][j]
				<< " | " << std::setw(16) << yAnalytical[i][j] << " | " << std::setw(16) << yAnalytical[i][j] - yApproximate[i][j]
				<< " | " << std::setw(16) << (yApproximate[i][j] - yAnalytical[i][j]) * 100.0 / yAnalytical[i][j] << std::fixed << std::endl;
		}
	}
}
