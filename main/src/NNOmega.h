#include <mass.h>

//NNOmega, 1 for Omega, 2 for N, 3 for N

#define pnOmega
//#define nnOmega
//#define ppOmega

#ifdef pnOmega
double m[3]={mOmega,mp,mn};
char nuclear[]="pnOmega";
#endif

#ifdef nnOmega
double m[3]={mOmega,mn,mn};
char nuclear[]="nnOmega";
#endif

#ifdef ppOmega
double m[3]={mOmega,mp,mp};
char nuclear[]="ppOmega";
#endif

double M=(m[1]+m[2]+m[3]);
double Q=1/sqrt(M/(m[1]*m[2]*m[3]));

double r12(double r, double alpha)
{
	double res=sqrt(Q*(m[1]+m[2])/(m[1]*m[2]))*r*cos(alpha);
	return res;
}

double r13(double r, double alpha)
{
	double res=sqrt(Q*(m[1]+m[3])/(m[1]*m[3]))*r*cos(alpha);
	return res;
}

double r23(double r, double alpha)
{
	double res=sqrt(Q*(m[2]+m[3])/(m[2]*m[3]))*r*cos(alpha);
	return res;
}
