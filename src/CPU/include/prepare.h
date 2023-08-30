#ifndef prepare_h
#define prepare_h
#include<cstddef>
#include<math.h>
#include<iostream>
#include<fstream>
using namespace std;
/*define what is a one-dimension matrix*/
class matrix;
class vector
{
	friend matrix;
	friend vector operator+(const vector& a,const vector& b);
	friend vector operator-(const vector& a,const vector& b);
	friend vector operator*(const double ration, const vector& a);
	friend vector operator/(const vector& a,double ration);
	friend vector operator*(const matrix& a,const vector & b);
	private:
		int size;
		double* v;
	public:
		/* construction function*/
		vector(int size=3);
		/*copy & construction function*/
		vector(const vector& tmp);
		/*deconstruction function*/
		~vector();
		/*find value according to index*/
		double& operator[](int index);
		/*define v=v*/
		vector &operator=(const vector& a);
		/*get size*/
		int Size() const {return size;}
};
//---------------------------------------------------------------------------------------------------------//
class matrix
{
	friend matrix operator*(double ration,const matrix& a);
	friend matrix operator*(const matrix& a,const matrix& b);
	friend vector operator*(const matrix& a,const vector & b);
	private:
		int size1,size2;
	public:
		vector * m;
		/*construction function*/
		matrix(int size1=3,int size2=3);
		/*copy & constuction function*/
		matrix(const matrix& tmp);
		/*deconstruction function*/
		~matrix();
		/*return colum*/
		vector& operator[](int index);
		/*M=M*/
		matrix& operator=(const matrix& a);
		/*get size1*/
		int Size1() const {return size1;};
		/*get size2*/
		int Size2() const {return size2;};
};
#endif
