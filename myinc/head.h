#include <fstream>
#include <iostream>
#include <ostream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <ctime>
#include <cstring>
#include <string.h>
#include <math.h>
#include <TMath.h>
#include <TMatrixDfwd.h>
#include <TMatrixDSymfwd.h>
#include <TMatrixDBase.h>
#include <TMatrixDSymEigen.h>
#include <TNtuple.h>
#include <TH2.h>
#include <THn.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TFrame.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TRandom3.h>
#include <TBenchmark.h>
#include <TInterpreter.h>
#include <stdio.h>
#include <stdlib.h>
#include <TComplex.h>
#include <TROOT.h>
#include <TFile.h>
#include <TH3.h>
#include <TF3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMinuit.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include "TVector3.h"
#include "TVector2.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TSystem.h"
#include "Math/SpecFuncMathMore.h"


double Pi=TMath::Pi();
double fmmev=0.00506;

double wigner_3j(int two_j1, int two_j2, int two_j3, int two_m1, int two_m2, int two_m3);
	//Calculates the Wigner 3j coupling coefficients:
	// (j1, j2, j3
	//  m1, m2, m3)
	//where j1,m1,...etc are integers or half integers. The function takes as input arguments only integers which corresponds to half integer units, e.g two_j1 = 2 * j1
	
double wigner_6j(int two_j1, int two_j2, int two_j3, int two_j4, int two_j5, int two_j6);
	//Calculates the Wigner 6j coupling coefficients: 
	// (j1, j2, j3
	//  j4, j5, j6)
	
double wigner_9j(int two_j1, int two_j2, int two_j3, int two_j4, int two_j5, int two_j6,  int two_j7, int two_j8, int two_j9);
	//Calculates the Wigner 6j coupling coefficients: 
	// (j1, j2, j3
	//  j4, j5, j6
	//  j7, j8, j9)
	
double CGcoeff(int two_j1, int two_j2, int two_j3, int two_m1, int two_m2, int two_m3);
	//Calculate the ClebschGordan coefficient (j1,m1;j2,m2|j1,j2;j3,m3)

double Wcoeff(int two_j1, int two_j2, int two_j3, int two_j23, int two_j12, int two_J);
	//Calculate the Racah W coefficient (j1,j2j3(j23);J|j1j2(j12),j3;J)

double assoc_laguerre(unsigned n, double m, double x);

/////////////////////////////////////////////


double wigner_3j(int two_j1, int two_j2, int two_j3, int two_m1, int two_m2, int two_m3)
{
	return ROOT::Math::wigner_3j(two_j1, two_j2, two_j3, two_m1, two_m2, two_m3);
}

double wigner_6j(int two_j1, int two_j2, int two_j3, int two_j4, int two_j5, int two_j6)
{
	return ROOT::Math::wigner_6j(two_j1, two_j2, two_j3, two_j4, two_j5, two_j6);
}

double wigner_9j(int two_j1, int two_j2, int two_j3, int two_j4, int two_j5, int two_j6,  int two_j7, int two_j8, int two_j9)
{
	return ROOT::Math::wigner_9j(two_j1, two_j2, two_j3, two_j4, two_j5, two_j6, two_j7, two_j8, two_j9);
}

double CGcoeff(int two_j1, int two_j2, int two_j3, int two_m1, int two_m2, int two_m3)
{
  double W3j=wigner_3j(two_j1, two_j2, two_j3, two_m1, two_m2, -two_m3);
	int negpow=(two_m3+two_j1-two_j2)/2;
	int neg=pow(-1,negpow);
	return neg*sqrt(two_j3+1)*W3j;
}

double Wcoeff(int two_j1, int two_j2, int two_j12, int two_j3, int two_J, int two_j23)
{
	int negpow=(two_j1+two_j2+two_j3+two_J)/2;
	int neg=pow(-1,negpow);
	return neg*sqrt((two_j12+1)*(two_j23+1))*wigner_6j(two_j1,two_j2,two_j12,two_j3,two_J,two_j23);
}

double assoc_laguerre(unsigned n, double m, double x)
{
	return ROOT::Math::assoc_laguerre(n,m,x);
}
