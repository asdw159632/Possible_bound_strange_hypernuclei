#include "mass.h"
#include "blw_para_particle.h"

//NNOmega, 1 for Omega, 2 for N, 3 for N
#define NNOmega

#define pnOmega
//#define nnOmega
//#define ppOmega

//Omega
double Tkin_p1 = Tkin_Omega;//GeV
double rho_0_p1 = rho_0_Omega;
const int num_p1 = num_Omega;//need to divide 100*1.8

char ti_p1_pT_Dst[]="Omega_pT_Dst";
char ti_func_pT_p1[]="func_pT_Omega";
char ti_func_coordinates_p1[]="func_coordinates_Omega";

#ifdef pnOmega

#define coulomb12 1

double m[3]={mOmega,mp,mn};
char nuclear[]="pnOmega";
const int two_T=0;
const int two_J=5;
const int tjkmin=0;
const int tjkmax=0;

//proton
double Tkin_p2 = Tkin_nucleon;//GeV
double rho_0_p2 = rho_0_nucleon;
const int num_p2 = num_nucleon;//_arry=40;//int((gRandom->Rndm()-0.5)*10.)+35;//38;//

//neutron
double Tkin_p3 = Tkin_nucleon;//GeV
double rho_0_p3 = rho_0_nucleon;
const int num_p3 = num_nucleon;//_arry=40;//int((gRandom->Rndm()-0.5)*10.)+35;//38;//

char ti_p2_pT_Dst[]="proton_pT_Dst";
char ti_p3_pT_Dst[]="neutron_pT_Dst";

char ti_func_pT_p2[]="func_pT_proton";
char ti_func_pT_p3[]="func_pT_neutron";

char ti_func_coordinates_p2[]="func_coordinates_proton";
char ti_func_coordinates_p3[]="func_coordinates_neutron";

Double_t GA = 3./8.;
//GA=(SA*2+1)/((s1*2+1)*...*(sn*2+1))/NI;
//NI counts the iospin state; take it as 1
//A. Polleri, R. Mattiello, I.N. Mishustin, J.P. Bondorf, Nucl. Phys. A 661 (1999) 452.

#endif

#ifdef nnOmega
#define ABB
double m[3]={mOmega,mn,mn};
char nuclear[]="nnOmega";
const int two_T=2;
const int two_J=3;
const int tjkmin=1;
const int tjkmax=1;

//neutron
double Tkin_p2 = Tkin_nucleon;//GeV
double rho_0_p2 = rho_0_nucleon;
const int num_p2 = num_nucleon;//_arry=40;//int((gRandom->Rndm()-0.5)*10.)+35;//38;//

char ti_p2_pT_Dst[]="neutron_pT_Dst";

char ti_func_pT_p2[]="func_pT_neutron";

char ti_func_coordinates_p2[]="func_coordinates_neutron";

Double_t GA = 1./4.;
//GA=(SA*2+1)/((s1*2+1)*...*(sn*2+1))/NI;
//NI counts the iospin state; take it as 1
//A. Polleri, R. Mattiello, I.N. Mishustin, J.P. Bondorf, Nucl. Phys. A 661 (1999) 452.

#endif

#ifdef ppOmega
#define ABB

#define coulomb12 1
#define coulomb23 -1
#define coulomb13 1

double m[3]={mOmega,mp,mp};
char nuclear[]="ppOmega";
const int two_T=2;
const int two_J=3;
const int tjkmin=1;
const int tjkmax=1;

//proton
double Tkin_p2 = Tkin_nucleon;//GeV
double rho_0_p2 = rho_0_nucleon;
const int num_p2 = num_nucleon;//_arry=40;//int((gRandom->Rndm()-0.5)*10.)+35;//38;//

char ti_p2_pT_Dst[]="proton_pT_Dst";

char ti_func_pT_p2[]="func_pT_proton";

char ti_func_coordinates_p2[]="func_coordinates_proton";

Double_t GA = 1./4.;
//GA=(SA*2+1)/((s1*2+1)*...*(sn*2+1))/NI;
//NI counts the iospin state; take it as 1
//A. Polleri, R. Mattiello, I.N. Mishustin, J.P. Bondorf, Nucl. Phys. A 661 (1999) 452.

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
