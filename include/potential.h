//Initialization of potential
//The NN potential parameter 
double C1_10 = -513.968*fmmev;//C1_10, the 10 represent the isospin-spin,the C1_01, VNN_10, VNN_01 is the same
double C1_01 = -626.885*fmmev;
double mu1 = 1.55;
double C2 = 1438.72*fmmev;
double mu2 = 3.11;

//The NOmega potential parameter
double mpi = 146*fmmev;//use the effective mass as HAL QCD result
double b[] = {-313.0*fmmev, 81.7, -252*fmmev, 0.85};//choose the parameters of t/a=12

//Th3 OmegaOmega potential parameter
double cj[] = {914, 305, -112};
double dj[] = {0.143*fmmev, 0.305*fmmev, 0.949*fmmev};

//The Coulomb potential perameter
double alpha = 1/137.036;

double VNN_10(double r)//The potential between NN, with (I,S)=(1,0), like nn pp
{
	return C1_10 * exp(-mu1 * r)/r + C2 * exp(-mu2 * r)/r;
}

double VNN_01(double r)//The potential between NN, with (I,S)=(0,1), like pn
{
	return C1_01 * exp(-mu1 * r)/r + C2 * exp(-mu2 * r)/r;
}

double VNOmega(double r)//The potential between NOmega
{
	return b[0] * exp(-b[1] * r*r) + b[2] * (1 - exp(-b[3] * r*r))*TMath::Power(exp(-mpi * r)/r,2);
}

double VOmegaOmega(double r)//The potential between OmegaOmega
{
	double result_VOmeOme=0;
	for(int j=0;j<3;j++)result_VOmeOme+=cj[j]*exp(-r*r/dj[j]/dj[j]);
	return result_VOmeOme;
}

double Vc(double r)
{
	return -alpha/r;
}

