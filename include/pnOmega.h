#include "mass.h"
#include "blw_para_particle.h"

//NNOmega, 1 for Omega, 2 for N, 3 for N
#define two_NOmega_interaction

#define coulomb12 1

double m[3]={mOmega,mp,mn};
double M=(m[0]+m[1]+m[2]);
double Q=1/sqrt(M/(m[0]*m[1]*m[2]));
const int two_s1=3;
const int two_s2=1;
const int two_s3=1;
const int two_t1=0;
const int two_t2=1;
const int two_t3=1;
const int Ptotal=1;
const int two_T=0;
const int two_J=5;
const int sjkmin=0;
const int sjkmax=1;
const int tjkmin=0;
const int tjkmax=0;

char nuclear[]="pnOmega";

//Omega
double Tkin_p1 = Tkin_Omega;//GeV
double rho_0_p1 = rho_0_Omega;
const int num_p1 = num_Omega;//need to divide 100*1.8

//proton
double Tkin_p2 = Tkin_nucleon;//GeV
double rho_0_p2 = rho_0_nucleon;
const int num_p2 = num_nucleon;//_arry=40;//int((gRandom->Rndm()-0.5)*10.)+35;//38;//

//neutron
double Tkin_p3 = Tkin_nucleon;//GeV
double rho_0_p3 = rho_0_nucleon;
const int num_p3 = num_nucleon;//_arry=40;//int((gRandom->Rndm()-0.5)*10.)+35;//38;//

char ti_p1_pT_Dst[]="Omega_pT_Dst";
char ti_p2_pT_Dst[]="proton_pT_Dst";
char ti_p3_pT_Dst[]="neutron_pT_Dst";

char ti_func_pT_p1[]="func_pT_Omega";
char ti_func_pT_p2[]="func_pT_proton";
char ti_func_pT_p3[]="func_pT_neutron";

char ti_func_coordinates_p1[]="func_coordinates_Omega";
char ti_func_coordinates_p2[]="func_coordinates_proton";
char ti_func_coordinates_p3[]="func_coordinates_neutron";

Double_t GA = 3./8.;
//GA=(SA*2+1)/((s1*2+1)*...*(sn*2+1))/NI;
//NI counts the iospin state; take it as 1
//A. Polleri, R. Mattiello, I.N. Mishustin, J.P. Bondorf, Nucl. Phys. A 661 (1999) 452.
