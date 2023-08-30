#include "solver.h"
#include <math.h>
using namespace std;
void displayProgressBar(int progress, int total, int barWidth = 60);

/*update data*/
void Euler::refresh()
{
	for (int i = 0; i <= size - 1; i++)
	{
		data[i].rou = U[i][rou_r];
		data[i].u = U[i][rou_u] / U[i][rou_r];
		data[i].E = U[i][rou_e] / U[i][rou_r];
		data[i].p = (data[i].E * data[i].rou - 0.5 * data[i].rou * pow(data[i].u, 2)) * (gamma - 1);
		data[i].T = data[i].p / data[i].rou * gamma * pow(data[i].u / data[i].c, 2);
		data[i].c = pow(gamma * data[i].p / data[i].rou, 0.5);
	}
}
/*calculate new U*/
void Euler::calculation(vector *R)
{
	for (int i = 1; i <= size - 1; i++)
		U[i] = U[i] + dt * R[i];
	refresh();
}
void Euler::calculation(vector *u0, vector *R, double ration)
{
	for (int i = 1; i <= size - 1; i++)
		U[i] = (1 - ration) * u0[i] + ration * U[i] + (ration * dt) * R[i];
	refresh();
}
/*calculate F using U*/
vector Euler::U2F(vector U1)
{
	vector F;
	double density, velocity, energy, pressure;
	density = U1[rou_r];
	velocity = U1[rou_u] / density;
	energy = U1[rou_e] / density;
	pressure = (energy * density - 0.5 * density * pow(velocity, 2)) * (gamma - 1);
	F[0] = density * velocity;
	F[1] = density * pow(velocity, 2) + pressure;
	F[2] = (density * energy + pressure) * velocity;
	return F;
}
/*rou average process*/
inline void rouAverage(vector L, vector R, double &uu, double &hh, double &cc)
{
	double roul = L[rou_r], rour = R[rou_r];
	double ul = L[rou_u] / roul, ur = R[rou_u] / rour;
	double el = L[rou_e] / roul, er = R[rou_e] / rour;
	double hl, hr, pl, pr;
	double gamma = 1.4;
	pl = (el * roul - 0.5 * roul * pow(ul, 2)) * (gamma - 1);
	pr = (er * rour - 0.5 * rour * pow(ur, 2)) * (gamma - 1);
	hl = el + pl / roul;
	hr = er + pr / rour;
	uu = (pow(roul, 0.5) * ul + pow(rour, 0.5) * ur) / (pow(roul, 0.5) + pow(rour, 0.5));
	hh = (pow(roul, 0.5) * hl + pow(rour, 0.5) * hr) / (pow(roul, 0.5) + pow(rour, 0.5));
	cc = (gamma - 1) * (hh - 0.5 * pow(uu, 2));
	cc = pow(cc, 0.5);
}
/*A*/
inline matrix A_cal(double uu, double cc)
{
	matrix tmp(3, 3);
	for (int i = 0; i <= 2; i++)
		for (int j = 0; j <= 2; j++)
			tmp.m[i][j] = 0.0;
	tmp.m[0][0] = abs(uu - cc);
	tmp.m[1][1] = abs(uu);
	tmp.m[2][2] = abs(uu + cc);
	return tmp;
}
/*R*/
inline matrix R_cal(double uu, double cc, double H)
{
	matrix tmp(3, 3);
	for (int i = 0; i <= 2; i++)
		for (int j = 0; j <= 2; j++)
			tmp.m[i][j] = 0.0;
	tmp.m[0][0] = tmp.m[0][1] = tmp.m[0][2] = 1;
	tmp.m[1][0] = uu - cc;
	tmp.m[1][1] = uu;
	tmp.m[1][2] = uu + cc;
	tmp.m[2][0] = H - uu * cc;
	tmp.m[2][1] = pow(uu, 2) / 2;
	tmp.m[2][2] = H + uu * cc;
	return tmp;
}
/*L*/
inline matrix L_cal(double uu, double cc)
{
	matrix tmp(3, 3);
	double gamma = 1.4;
	for (int i = 0; i <= 2; i++)
		for (int j = 0; j <= 2; j++)
			tmp.m[i][j] = 0.0;
	tmp.m[0][0] = 0.5 * (0.5 * (gamma - 1) * pow(uu, 2) / pow(cc, 2) + uu / cc);
	tmp.m[0][1] = -0.5 * ((gamma - 1) * uu / pow(cc, 2) + 1 / cc);
	tmp.m[0][2] = 0.5 * (gamma - 1) / pow(cc, 2);
	tmp.m[1][0] = 1 - 0.5 * (gamma - 1) * pow(uu, 2) / pow(cc, 2);
	tmp.m[1][1] = (gamma - 1) * uu / pow(cc, 2);
	tmp.m[1][2] = -(gamma - 1) / pow(cc, 2);
	tmp.m[2][0] = 0.5 * (0.5 * (gamma - 1) * pow(uu, 2) / pow(cc, 2) - uu / cc);
	tmp.m[2][1] = -0.5 * ((gamma - 1) * uu / pow(cc, 2) - 1 / cc);
	tmp.m[2][2] = 0.5 * (gamma - 1) / pow(cc, 2);
	return tmp;
}
//----------------------------------------------------------------------------------------------------------//
Euler::Euler(int size, double t)
{
	this->size = size;
	this->t = t;
	data = new status[size];
	U = new vector[size];
	dx = 1.0 / size;
	dt = 0.0001;
	for (int i = 0; i <= size - 1; i++)
	{
		data[i].x = (i + 1) * dx;
		if (data[i].x < 0.5)
		{
			data[i].u = 0;
			data[i].rou = 1;
			data[i].p = 1;
		}
		else
		{
			data[i].u = 0;
			data[i].rou = 0.125;
			data[i].p = 0.1;
		}
		data[i].E = data[i].p / (gamma - 1) / data[i].rou + 0.5 * pow(data[i].u, 2);
		data[i].T = data[i].p / data[i].rou * gamma * pow(data[i].u / data[i].c, 2);
		data[i].c = pow(gamma * data[i].p / data[i].rou, 0.5);
		U[i][rou_r] = data[i].rou;
		U[i][rou_u] = data[i].rou * data[i].u;
		U[i][rou_e] = data[i].rou * data[i].E;
	}
}
/*rou formate*/
void Euler::Rou(vector *&R)
{
	vector Ul, Ur;
	vector *F_2;
	F_2 = new vector[size];
	double uu, hh, cc; // aveage vaule;
	matrix mA, mR, mL;
	for (int i = 0; i <= size - 1; i++)
	{
		if (i < size - 1)
			Ur = U[i + 1];
		else
			Ur = U[i];
		Ul = U[i];
		rouAverage(Ul, Ur, uu, hh, cc);
		mA = A_cal(uu, cc);
		mR = R_cal(uu, cc, hh);
		mL = L_cal(uu, cc);
		F_2[i] = 0.5 * (U2F(Ul) + U2F(Ur));
		F_2[i] = F_2[i] - 0.5 * (mR * mA * mL) * (Ur - Ul);
	}

	for (int i = 0; i <= size - 1; i++)
	{
		if (i > 0)
			R[i] = (-1.0) * (F_2[i] + (-1.0) * F_2[i - 1]) / dx;
		else
			R[i] = F_2[i] - F_2[i];
	}
	delete[] F_2;
}
/*the process of runge-kutta*/
void Euler::Runge_Kutta()
{
	vector *U0, *R;
	U0 = new vector[size];
	R = new vector[size];
	for (int i = 0; i <= size - 1; i++)
		U0[i] = U[i];
	Rou(R);
	calculation(R);
	Rou(R);
	calculation(U0, R, 1.0 / 4);
	Rou(R);
	calculation(U0, R, 2.0 / 3);
	delete[] U0, R;
}
/*solve the equation*/
void Euler::solver()
{
	double tt = 0;
	int count = 0;
	while (tt < t)
	{
		count++;
		displayProgressBar(count, t / dt);
		Runge_Kutta();
		tt = tt + dt;
	}
	cout << endl;
}
/*delete the memory*/
Euler::~Euler()
{
	if (data != NULL)
		delete[] data;
	if (U != NULL)
		delete[] U;
}
/*output result*/
void Euler::output(string name) const
{
	name.append(".dat");
	ofstream out(name);
	for (int i = 0; i <= size - 1; i++)
	{
		// out << data[i].x << '\t' << data[i].p << endl;
		out << data[i].x << '\t' << data[i].u << '\t' << data[i].p
			<< '\t' << data[i].rou << '\t' << data[i].T << endl;
	}
	out.close();
}
