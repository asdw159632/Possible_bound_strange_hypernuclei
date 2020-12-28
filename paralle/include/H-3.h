#include "mass.h"

//H_3, 1 for p, 2 for n, 3 for n
const double m[3]={mp,mn,mn};
const double M=(m[0]+m[1]+m[2]);
const double Q=1/sqrt(M/(m[0]*m[1]*m[2]));
const int two_T=1;
const int two_J=1;
const int two_s1=1;
const int two_s2=1;
const int two_s3=1;
const int two_t1=1;
const int two_t2=1;
const int two_t3=1;
const int tjkmin=0;
const int tjkmax=1;
const int Ptotal=1;
char nuclear[]="H_3";
#define H_3

double r12(double r, double alpha)
{
	double res=sqrt(Q*(m[0]+m[1])/(m[0]*m[1]))*r*cos(alpha);
	return res;
}

double r13(double r, double alpha)
{
	double res=sqrt(Q*(m[0]+m[2])/(m[0]*m[2]))*r*cos(alpha);
	return res;
}

double r23(double r, double alpha)
{
	double res=sqrt(Q*(m[1]+m[2])/(m[1]*m[2]))*r*cos(alpha);
	return res;
}

