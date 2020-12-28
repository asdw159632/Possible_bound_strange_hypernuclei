//Initialization of potential
//The NN potential parameter 
double C1_10 = -513.968*fmmev;//C1_10, the 10 represent the isospin-spin,the C1_01, VNN_10, VNN_01 is the same
double C1_01 = -626.885*fmmev;
double mu1 = 1.55;
double C2 = 1438.72*fmmev;
double mu2 = 3.11;

//The NOmega potential parameter
double mpi = 146;//use the effective mass as HAL QCD result
double b[] = {-313.0*fmmev, 81.7, -252*fmmev, 0.85};//choose the parameters of t/a=12

//The Coulomb potential perameter
double alpha = 1/137.036;

double VNN_10(double r)//The potential between NN, with (I,S)=(1,0), like nn pp
{
	double result_VNN_10 = C1_10 * exp(-mu1 * r)/r + C2 * exp(-mu2 * r)/r;
	return result_VNN_10;
}

double VNN_01(double r)//The potential between NN, with (I,S)=(0,1), like pn
{
	double result_VNN_01 = C1_01 * exp(-mu1 * r)/r + C2 * exp(-mu2 * r)/r;
	return result_VNN_01;
}

double VNOmega(double r)//The potential between NOmega
{
	double result_VNOme = b[0] * exp(-b[1] * r*r) + b[2] * (1 - exp(-b[3] * r*r))*TMath::Power(exp(-mpi * r)/r,2);
	return result_VNOme;
}

double Vc(double r)
{
	double result_Vc=-alpha/r;
	return result_Vc;
}

