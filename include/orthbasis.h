
double orthbasis_radial(double r, int n);
double fun_orthbasis_radial(double *x, int *par);
double orthbasis_alpha(int q, int lx, int ly, double alpha);
double fun_orthbasis(double *x, double *par);





///////////////////////////////////////////////////////////////////////

double orthbasis_radial(double r, int n)//number l and variation parameter c is defined as global.
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

double orthbasis_alpha(int q, int lx, int ly, double alpha)
{
	//Kievsky, A., Viviani, M., & Rosati, S. (1993). Nuclear Physics A, 551(2) 241;
	//Kovalchuk, V. I. (2014). International Journal of Modern Physics E, 23(11), 1450069.
	//Eigen function of hyper angle momentum with wigen value K=2*q+lx+ly;
	double Nk=sqrt(2*TMath::Factorial(q)*(2*q+lx+ly+2)*TMath::Factorial(q+lx+ly+1)/(tgamma(q+lx+1.5)*tgamma(q+ly+1.5)));
	return Nk*pow(cos(alpha),lx)*pow(sin(alpha),ly)*JacobiP(q,ly+0.5,lx+0.5,cos(2*alpha));
}

double fun_orthbasis_radial(double *x, double *par)
{
	double r=x[0];
	int n=par[0];
	return orthbasis_radial(r,n);
}

double fun_orthbasis(double *x, double *par)
{
	double r=x[0];
	double alpha=x[1];
	int n=par[0];
	int q=par[1];
	int lx=par[2];
	int ly=par[3];
	return orthbasis_radial(r,n)*orthbasis_alpha(q,lx,ly,alpha);
}

