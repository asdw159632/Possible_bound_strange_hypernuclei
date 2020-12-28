#include <mass.h>

//He_3, 1 for p, 2 for p, 3 for n
double m[3]={mp,mp,mn};
double M=(m[1]+m[2]+m[3]);
double Q=1/sqrt(M/(m[1]*m[2]*m[3]));
char nuclear[]="He_3";

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
