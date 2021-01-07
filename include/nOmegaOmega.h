#include "mass.h"
#include "blw_para_particle.h"

//NOmegaOmega, 1 for N, 2 for Omega, 3 for Omega
#define two_NOmega_interaction

#define ABB
double m[3]={mn,mOmega,mOmega};
double M=(m[0]+m[1]+m[2]);
double Q=1/sqrt(M/(m[0]*m[1]*m[2]));
const int two_s1=1;
const int two_s2=3;
const int two_s3=3;
const int two_t1=1;
const int two_t2=0;
const int two_t3=0;
const int Ptotal=1;
const int two_T=1;
const int two_J=1;
const int sjkmin=0;
const int sjkmax=0;
const int tjkmin=0;
const int tjkmax=0;

char nuclear[]="nOmegaOmega";

//neutron
double Tkin_p1 = Tkin_nucleon;//GeV
double rho_0_p1 = rho_0_nucleon;
const int num_p1 = num_nucleon;//_arry=40;//int((gRandom->Rndm()-0.5)*10.)+35;//38;//

//Omega
double Tkin_p2 = Tkin_Omega;//GeV
double rho_0_p2 = rho_0_Omega;
const int num_p2 = num_Omega;//need to divide 100*1.8

char ti_p1_pT_Dst[]="neutron_pT_Dst";
char ti_p2_pT_Dst[]="Omega_pT_Dst";

char ti_func_pT_p1[]="func_pT_neutron";
char ti_func_pT_p2[]="func_pT_Omega";

char ti_func_coordinates_p1[]="func_coordinates_neutron";
char ti_func_coordinates_p2[]="func_coordinates_Omega";

Double_t GA = 1./16.;
//GA=(SA*2+1)/((s1*2+1)*...*(sn*2+1))/NI;
//NI counts the iospin state; take it as 1
//A. Polleri, R. Mattiello, I.N. Mishustin, J.P. Bondorf, Nucl. Phys. A 661 (1999) 452.
