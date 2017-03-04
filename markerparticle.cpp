//Marker particle method simulation to model one-dimensional interface 
//Initial interface is defined by ro(s) = (1 + 0.5cos(4s))(cos(s), sin(s)) for s in [0, 2*pi]
//Silmulation with 128 marker particles and time steps of 10^-3, 5*10^-4 and 2.5*10^-4
//No reparametrization included

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <cassert>
#include <vector>
#include <fstream>

using namespace std;


int main(int argc, char* argv[])
{
	//Defining number of marker particles
	int n = 128;

	//defining variables for mean curvature and force function
	double kbar, k_sum = 0.0;
	std::vector<double> k(n, 0.0);
	std::vector<double> F(n, 0.0);

	//Defining the time steps
	double dt = 2.5*pow(10, -4);

	//Defining the step size to accomodate 128 points on the initial interface
	double ds = (2 * M_PI) / n;

	//Defining the points on the interface with step size ds
	std::vector<double> s(n, 0.0);

	for (int i = 1; i < s.size(); i++)
	{
		s[i] = s[i - 1] + ds;
	}

	std::vector<double> a(s.size(), 0.0);
	for (int i = 0; i < n; i++)
	{
		a[i] = (1 + 0.5*cos(4 * s[i]));
	}

	std::vector<double> X(s.size(), 0.0);
	std::vector<double> Y(s.size(), 0.0);
	std::vector<double> xt(n, 0.0);
	std::vector<double> yt(n, 0.0);

	for (int i = 0; i < s.size(); i++)
	{
		X[i] = a[i] * cos(s[i]);
		Y[i] = a[i] * sin(s[i]);
	}

	//Defining the first and second difference parameters dX, dY, d2X, d2Y and initializing them to zero
	std::vector<double> dX(s.size(), 0.0);
	std::vector<double> d2X(s.size(), 0.0);
	std::vector<double> dY(s.size(), 0.0);
	std::vector<double> d2Y(s.size(), 0.0);

	//Setting the time intervals and obtaining different time dependin on time steps.
	//The simulation is to be run to t = 0.5,
	double tmax = 0.5;
	double tmin = 0.0;
	double tn = (tmax - tmin) / dt;
	std::vector<double> t(tn, 0.0);

	//Marching in time
	for (int i = 0; i < tn; i++)
	{
		//Implementing Central Difference scheme to calculate the mean curvature
		for (int i = 1; i < n - 1; i++)
		{
			dX[i] = (X[i + 1] - X[i - 1]) / (2 * ds);
			dY[i] = (Y[i + 1] - Y[i - 1]) / (2 * ds);
			d2X[i] = (X[i + 1] - 2 * X[i] + X[i - 1]) / (ds*ds);
			d2Y[i] = (Y[i + 1] - 2 * Y[i] + Y[i - 1]) / (ds*ds);
		}

		//Calculating 1st and last point considering the periodic boundary condition
		dX[0] = (X[1] - X[n - 1]) / (2 * ds);
		dX[n - 1] = (X[0] - X[n - 2]) / (2 * ds);
		dY[0] = (Y[1] - Y[n - 1]) / (2 * ds);
		dY[n - 1] = (Y[0] - Y[n - 2]) / (2 * ds);

		d2X[0] = (X[1] - 2 * X[0] + X[n - 1]) / (ds*ds);
		d2X[n - 1] = (X[0] - 2 * X[n - 1] + X[n - 2]) / (ds*ds);
		d2Y[0] = (Y[1] - 2 * Y[0] + Y[n - 1]) / (ds*ds);
		d2Y[n - 1] = (Y[0] - 2 * Y[n - 1] + Y[n - 2]) / (ds*ds);

		//Calculating the mean curvature
		for (int i = 0; i < n; i++)
		{
			k[i] = (dX[i] * d2Y[i] - dY[i] * d2X[i]) / (pow((pow(dX[i], 2) + pow(dY[i], 2)), 1.5));
		}

		for (int i = 0; i < k.size(); i++)
		{
			k_sum += k[i];
		}

		kbar = k_sum / k.size();

		//Calculating the value of F
		for (int i = 0; i < n; i++)
		{
			F[i] = kbar - k[i];
		}

		//Calculating the velocity functional in x and y directions
		for (int i = 0; i < n; i++)
		{
			xt[i] = X[i] + dt*((F[i] * dY[i]) / (pow((pow(dX[i], 2) + pow(dY[i], 2)), 0.5)));
			yt[i] = Y[i] + dt*((F[i] * dX[i]) / (pow((pow(dX[i], 2) + pow(dY[i], 2)), 0.5)));
		}

		//Updating the values of X and Y
		for (int i = 0; i < n; i++)
		{
			X[i] = xt[i];
			Y[i] = yt[i];
		}
		//X[n] = X[0];
		//Y[n] = Y[0];
	}	

	return 0;
}