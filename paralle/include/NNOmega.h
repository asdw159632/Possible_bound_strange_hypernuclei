#include "mass.h"

//NNOmega, 1 for Omega, 2 for N, 3 for N
#define NNOmega

#define pnOmega
//#define nnOmega
//#define ppOmega

#ifdef pnOmega
double m[3]={mOmega,mp,mn};
char nuclear[]="pnOmega";
const int two_T=0;
const int two_J=5;
const int tjkmin=0;
const int tjkmax=0;
#endif

#ifdef nnOmega
double m[3]={mOmega,mn,mn};
char nuclear[]="nnOmega";
const int two_T=2;
const int two_J=3;
const int tjkmin=1;
const int tjkmax=1;
#endif

#ifdef ppOmega
double m[3]={mOmega,mp,mp};
char nuclear[]="ppOmega";
const int two_T=2;
const int two_J=3;
const int tjkmin=1;
const int tjkmax=1;
#endif

double M=(m[0]+m[1]+m[2]);
double Q=1/sqrt(M/(m[0]*m[1]*m[2]));
const int two_s1=3;
const int two_s2=1;
const int two_s3=1;
const int two_t1=0;
const int two_t2=1;
const int two_t3=1;
const int Ptotal=1;

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
