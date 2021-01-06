#include "mass.h"
#include "blw_para_particle.h"

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
const int Ptotal=1;
const int sjkmin=0;
const int sjkmax=1;
const int tjkmin=0;
const int tjkmax=1;

char nuclear[]="H_3";

#define ABB

//proton
double Tkin_p1 = Tkin_nucleon;//GeV
double rho_0_p1 = rho_0_nucleon;
const int num_p1 = num_nucleon;//_arry=40;//int((gRandom->Rndm()-0.5)*10.)+35;//38;//

//neutron
double Tkin_p2 = Tkin_nucleon;//GeV
double rho_0_p2 = rho_0_nucleon;
const int num_p2 = num_nucleon;//_arry=40;//int((gRandom->Rndm()-0.5)*10.)+35;//38;//

char ti_p1_pT_Dst[]="proton_pT_Dst";
char ti_p2_pT_Dst[]="neutron_pT_Dst";

char ti_func_pT_p1[]="func_pT_proton";
char ti_func_pT_p2[]="func_pT_neutron";

char ti_func_coordinates_p1[]="func_coordinates_proton";
char ti_func_coordinates_p2[]="func_coordinates_neutron";

Double_t GA = 1./4.;
//g0=(SA*2+1)/((s1*2+1)*...*(sn*2+1))/NI;
//NI counts the iospin state; take it as 1 in this case;
//A. Polleri, R. Mattiello, I.N. Mishustin, J.P. Bondorf, Nucl. Phys. A 661 (1999) 452.
