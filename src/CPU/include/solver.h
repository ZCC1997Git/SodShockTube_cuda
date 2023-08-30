#ifndef solver_h
#define solver_h
#include"prepare.h"
#include<iostream>
#include<fstream>
using namespace std;
/*define some variables which respond to the physical status */
struct status
{
	double x;                 /*position*/
	double rou;               /*density;*/
	double u;                 /*velocity;*/
	double p;                 /*pressure;*/
	double E;                 /*total energy;*/
	double T;                 /*tmpture*/
	double c;                 /*the speed of sound*/
};
/*help to specific the physical meaning of U and F for different index*/
typedef enum
{
	rou_r,                     /*density*/
	rou_u,                     /*density x velocity*/
	rou_e,                     /*density x energy*/
} U;
//------------------------------------------------------//
class Euler
{
private:
	const double R0=287.0;
	const double gamma=1.4;
	int size;
	double dt,dx,t;
	status* data;
	vector* U;
	/*update data*/
	void refresh();
	/*calculate new U*/
	void calculation(vector* R);
	void calculation(vector* u0,vector* R,double ration);
	/*calculate F using U*/
	vector U2F(vector U1);
public:
	Euler(int size,double t);
	~Euler();
	void Rou(vector*& R);
	void Runge_Kutta();
	void solver();
	void output(string name)const;
};
//------------------------------------------------------//
#endif 
