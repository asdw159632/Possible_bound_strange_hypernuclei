#include <head.h>

#define fmmev 0.00506
//Initialization of mass
double mp=938.27*fmmev;
double mn=939.57*fmmev;
double mOmega=1672.45*fmmev;

//Initialization of mass done

/////////////////////////

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

//The Coulomb potential perameter
double alpha = 1/137.036;

double c=30;
int l=1;
double Pi=TMath::Pi();
double Q=sqrt((mp+mn+mn)/(mp*mn*mn));
char nuclear[]="H_3";

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

double V(double r, double alpha)//V(r)=Sum_{i<j} {V_ij(TMath::Sqrt(Q(mi+mj)/(mi*mj))*r*cos(alpha))} do not integrated by alpha, the hypersphercial angle.
{
	double result_V;
	/*TMath::Sqrt(Q(mi+mj)/(mi*mj))*r*cos(alpha) is actually the distance between i and j in different particle index set,
	  r is indepent of particle index set, but alpha do, so we need to do an approximate to integrate alpha.
	  The integrate of alpha is also independ of alpha, so in the different particle index set the alpha can be integrated in one integrate.*/
	double r12=TMath::Sqrt(Q*(mp+mn)/(mp*mn))*r*cos(alpha);
	double r13=TMath::Sqrt(Q*(mp+mn)/(mp*mn))*r*cos(alpha);
	double r23=TMath::Sqrt(Q*(mn+mn)/(mn*mn))*r*cos(alpha);
	result_V=VNN_01(r12)+VNN_10(r23)+VNN_01(r13);//alpha integrate toghter with r later
	return result_V;
}

double orthbasis(double r, int n)//orbital angular momentum number l and variation parameter c is defined as global.
{
	if(n<l+1)return 0;
	double result_basis=TMath::Sqrt(TMath::Power(2*c/n,3)*TMath::Factorial(-1 - l + n)/(2*n*TMath::Factorial(l + n)))*exp(-c*r/n)*TMath::Power(2*c*r/n, l)
	      *assoc_laguerre(-1-l+n, 1+2*l, 2*c*r/n)*r;
		/*The bases are given by hydrogen, it can be write as |nl> in the annotation.
	      orthbasis[r, n, l] = <r|nl> *r is the reduced radial wave function. 
	      The unit of reduced radial wave function is always [r]^(-1/2) in n-dim, that's why choose it as base.
	      Attention : the bases with different l are not orthogonal.
          And the c is the variational parameter.*/
	return result_basis;
}

double fun_orthbasis(double *x, double *par)//orbital angular momentum number l and variation parameter c is defined as global.
{
	double r=x[0];
	int n=(int)par[0];
	double result_basis=orthbasis(r,n);
	return result_basis;
}

double integral_H(double *x, double *Par)
{
	double r=x[0];//hyper spherical radius
	double alpha=x[1];//hyper spherical angle
	int n=Par[0];//n is matrix row index of H_nm
	int m=Par[1];//m is matrix column index of H_nm
	TF1 *basis_m=new TF1("basis_m",fun_orthbasis,0,50,1);
	basis_m->SetParameter(0,m);
	double d2basis_m=basis_m->Derivative2(r);
	double res_integrand_H=orthbasis(r,n)*(-d2basis_m/(2*Q)+(V(r,alpha)+15/(4*r*r*2*Q))
		                       *orthbasis(r,m))*cos(alpha)*cos(alpha)*sin(alpha)*sin(alpha)*16/Pi;
	  /*15/(4*x^2*2Q) come from the derivative of u, the reduce radial wave function. 
	    orthbasis=r*R_nl as the reduce radial wave function, R_nl is 3-dim radial wave function of hydrogen.
	    In this case u=r^(5/2)S_nl, S_nl is 6-dim radial wave function.
	    No matter in any dimensions the unit of u is always [r]^-1/2, so using orthbasis as the base of u in this case.
	    In 6-dim, d^2 u/dr^2=d^2 S/dr^2+5/r*dS/dr+15/(4*r^2)*S and the kinetic part is -(d^2 S/dr^2+5/r*dS/dr)/2m=-d^2 u/dr^2*1/2m+15/(4*r^2*2m)*S
	    The potential part add a extra potential from derivative.
	    And cos(alpha)*cos(alpha)*sin(alpha)*sin(alpha) is the volume element of alpha integral.
	    16/Pi is the normalization factor of alpha integral.*/
	delete basis_m;
	return res_integrand_H;
}

int test()
{
  TF2 *test=new TF2("test",integral_H,0,50,0,TMath::Pi()/2,2);
  double n=1;
  double m=1;
  l=0;
  test->SetParameters(n,m);
  test->Draw();
  delete test;
}
