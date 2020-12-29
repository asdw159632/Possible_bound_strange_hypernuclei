#include "mass.h"

//He_3, 1 for n, 2 for n, 3 for p
double m[3]={mn,mp,mp};
double M=(m[0]+m[1]+m[2]);
double Q=1/sqrt(M/(m[0]*m[1]*m[2]));
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

char nuclear[]="He_3";

#define ABB

//neutron
double Tkin_p1 = 0.1116;//GeV
double rho_0_p1 = 0.98;
const int num_p1 = 35;//_arry=40;//int((gRandom->Rndm()-0.5)*10.)+38;//

//proton
double Tkin_p2 = 0.1116;//GeV
double rho_0_p2 = 0.98;
const int num_p2 = 35;//_arry=40;//int((gRandom->Rndm()-0.5)*10.)+35;//38;//

/*proton2
double Tkin_p3 = 0.1116;//GeV
double rho_0_p3 = 0.98;
const int num_p3 = 35;//_arry=40;//int((gRandom->Rndm()-0.5)*10.)+35;//38;*/

char ti_p1_pT_Dst[]="neutron_pT_Dst";
char ti_p2_pT_Dst[]="proton_pT_Dst";
//char ti_p3_pT_Dst[]="proton2_pT_Dst";

char ti_func_pT_p1[]="func_pT_neutron";
char ti_func_pT_p2[]="func_pT_proton";
//char ti_func_pT_p3[]="func_pT_proton2";

char ti_func_coordinates_p1[]="func_coordinates_neutron";
char ti_func_coordinates_p2[]="func_coordinates_proton";
//char ti_func_coordinates_p3[]="func_coordinates_proton2";

Double_t GA = 1./3.;//GA=g0_A+2*g0_12->D*g0_D3
//g0=(SA*2+1)/((s1*2+1)*...*(sn*2+1))/NI;
//NI counts the iospin state
//A. Polleri, R. Mattiello, I.N. Mishustin, J.P. Bondorf, Nucl. Phys. A 661 (1999) 452.

#define He_3

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
