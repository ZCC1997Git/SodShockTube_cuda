#include"prepare.h"
/*******************************************************/
/*   define the friend function of class vector        */
/*******************************************************/
/*define v+v*/
vector operator+(const vector& a,const vector& b)
{
    if(a.Size()!=b.Size()) cerr<<"the Size() of vectors is not same"<<endl;
		vector tmp(a.Size());
		for(int i=0;i<=a.Size()-1;i++)
			tmp.v[i]=a.v[i]+b.v[i];
		return tmp;
}
/*define v-v*/
vector operator-(const vector& a,const vector& b)
{
    if(a.Size()!=b.Size()) cerr<<"the Size() of vectors is not same"<<endl;
    vector tmp(a.Size());
    for(int i=0;i<=a.Size()-1;i++)
        tmp.v[i]=a.v[i]-b.v[i];
    return tmp;
}
/*define a*v*/
vector operator*(const double ration, const vector& a)
{
    vector tmp(a.Size());
    for(int i=0;i<=a.Size()-1;i++)
        tmp.v[i]=ration*a.v[i];	
    return tmp;
}
/*define V/a*/
vector operator/(const vector& a,double ration)
{
    vector tmp(a.Size());
    for(int i=0;i<=a.Size()-1;i++)
        tmp.v[i]=a.v[i]/ration;
    return tmp;
}
/*******************************************************/
/*         define the function of class vector         */
/*******************************************************/
/*construction function*/
vector::vector(int size)
{
    this->size=size;
    v=new double[size];
}
/*copy & construction function*/
vector::vector(const vector& tmp)
{
    size=tmp.Size();
    v=new double[size];
    for(int i=0;i<=size-1;i++)
        v[i]=tmp.v[i];
}
/*deconstruction function*/
vector::~vector()
{
    if(v!=NULL) delete[] v;
}
/*find value according to index*/
double& vector::operator[](int index) 
{
    return v[index];
}
/*define v=v*/
vector& vector::operator=(const vector& a)
{
    if(this==&a) return *this;
    delete[] v;
    size=a.Size();
    v=new double[a.Size()];
    for(int i=0;i<=a.Size()-1;i++)
        v[i]=a.v[i];
    return *this;
}
/*******************************************************/
/*   define the friend function of class matrix        */
/*******************************************************/
/*define a*M */
matrix operator*(double ration,const matrix& a)
{
    matrix tmp(a.Size1(),a.Size2());
    for(int i=0;i<=a.Size1()-1;i++)
        tmp.m[i]=ration*a.m[i];
    return tmp;
} 
/*define M*M*/
matrix operator*(const matrix& a,const matrix& b)
{
    if(a.Size2()!=b.Size1()) cerr<<"the size of matrixes is not same"<<endl;
    matrix tmp(a.Size1(),b.Size2());
    for(int i=0;i<=a.Size1()-1;i++)
    {
        for(int j=0;j<=b.Size2()-1;j++)
        {
            tmp.m[i][j]=0;
            for(int k=0;k<=a.Size2()-1;k++)
                tmp.m[i][j]+=a.m[i][k]*b.m[k][j];
        }
    }
    return tmp;
}
/*define M*v*/
vector operator*(const matrix& a,const vector & b)
{
    if(a.Size2()!=b.size) cerr<<"the size of matrixes is not same"<<endl;
    vector tmp(a.Size1());
    for(int i=0;i<=a.Size1()-1;i++)
        tmp.v[i]=0;
    for(int i=0;i<=a.Size1()-1;i++)
    {
        for(int k=0;k<=a.Size2()-1;k++)
            tmp[i]+=a.m[i][k]*b.v[k];
    }
    return tmp;
}
/*******************************************************/
/*         define the function of class Matrix        */
/*******************************************************/
/*construction function*/
matrix::matrix(int size1,int size2)
{
    this->size1=size1;
    this->size2=size2;
    m=new vector[size1];
    for(int i=0;i<=size1-1;i++)
    {
        m[i].size=size2;
        delete[] m[i].v;
        m[i].v=new double[size2];
    }
}
/*copy & constuction function*/
matrix::matrix(const matrix& tmp)
{
    size1=tmp.Size1();
    size2=tmp.Size2();
    m=new vector[size1];
    for(int i=0;i<=size1-1;i++)
    {
        m[i].size=size2;
        delete[] m[i].v;
        m[i].v=new double[size2];
    }
    for(int i=0;i<=size1-1;i++)
        for(int j=0;j<=size2-1;j++)
            m[i][j]=tmp.m[i][j];
}
/*deconstruction function*/
matrix::~matrix()
{
    if(m!=NULL) delete[] m;
}
/*return colum*/
vector& matrix::operator[](int index)
{
    return m[index];
}
/*M=M*/
matrix& matrix::operator=(const matrix& a)
{
    if(this==&a) return *this;
    size1=a.Size1();
    size2=a.Size2();
    delete[] m;
    m=new vector[a.Size1()];
    for(int i=0;i<=size1-1;i++)
        m[i]=a.m[i];
    return *this;
}