#include <head.h>

using namespace std;

#include <jacobi_polynomial.hpp>

///////////////////////////////////////////////////////////////
//
//HyperSpherical method solving three-body nuclear wave function
//
///////////////////////////////////////////////////////////////
//Written by Zhang Liang, 2020-12
//The unit is all in fm
///////////////////////////////////////////////////////////////

///////////////////////////////////////
//Initialization
///////////////////////////////////////

//#include "../../include/mass.h"
#include "../../include/potential.h"

//choice a system, choice H_3, He_3 or NNOmega
#include "../../include/H-3.h"
//#include "../../include/He-3.h"
//#include "../../include/NNOmega.h"//choose different NN-pair in NNOmega.h

double V23(int sjk, double r, double alpha);
double V13(int sjk, double r, double alpha);
double V12(int sjk, double r, double alpha);
//Initialization of potential done

/////////////////////////

//Initialization of test wave function

#define nmax 5//nmax represent the number of used orthogonal bases.
double Y=1/TMath::Sqrt(1/4 * Pi);/*Spherical harmonic function: sph_legendre(l, m, theta)*exp(i*m*phi). 
    						For S-wave, l = m = 0, Y = 1/TMath::Sqrt(4 Pi).*/
double c;//variation parameter in orthbasis
int l=1;//number in orthbasis_radial. The bases with different l are not orthogonal, so set it as global.
#include "../../include/readRR.h"

#include "../../include/orthbasis.h"

/*double orthbasis_radial(double r, int n);
double fun_orthbasis_radial(double *x, int *par);
double orthbasis_alpha(int q, int lx, int ly, double alpha);
double fun_orthbasis(double *x, double *par);*/

//Initialization of test wave function done

/////////////////////////

//Initialization of variation
double c0=2000;//c0 is the begin variational parameter of every variational iteration loop
double c00=TMath::Power(10,(int) (floor(TMath::Log10(c0))-1));//A parameter which determine the precision/step-size of every variational iteration loop
#define Ncmax 20/*The span of one variational iteration loop is given by Ncmax*c00. 
			  It means, in this submit, c varies from c0-c00*(Ncmax/2) to c0+c00*(Ncmax/2)*/
void Hmin_reset();

//Initialization of variation done

///////////////////////////////////////
//Initialization done
///////////////////////////////////////

//////////////////////////////////////////////////
//////////////////////////////////////////////////

double integrand_H(double *x, double *Par);

int main(int argc, char *argv[])
{
	//ofstream save("/home/zhangliang/Possible_bound_strange_hypernuclei/code/bin/save.txt");
	//if(!save)return 0;

	cout<<"///////////////////////////////////////////////////////////////"<<endl;
	cout<<"//"<<endl;
	cout<<"//Hartree method solving three-body nuclear wave function"<<endl;
	cout<<"//"<<endl;
	cout<<"///////////////////////////////////////////////////////////////"<<endl;
	cout<<endl;
	cout<<"///////////////////////////////////////"<<endl;
	cout<<"//Variation"<<endl;
	cout<<"///////////////////////////////////////"<<endl;
	cout<<endl;

	if(readRR(0))
	{
		cerr<<"ERROR: READ RR_MATRIX FILED"<<endl;
		return 0;
	}

	int qmax=10;
	int lmax=4;

	cout<<"Nc  q   lx  ly  L   sjk two_Sa tjk"<<endl;

	int anglemomentlist[100][7];
	int aNc=0;
	for(int q=0;q<=qmax;q++)
	{
		for(int lx=0;lx<=lmax;lx++)
		{
			for(int ly=0;ly<=lmax;ly++)
			{
				if(pow(-1,lx+ly)!=Ptotal)continue;
				for(int sjk=0;sjk<=1;sjk++)
				{
					for(int tjk=tjkmin;tjk<=tjkmax;tjk++)
					{
						if(pow(-1,sjk+tjk+lx)!=-1)continue;
						for(int two_Sa=abs(sjk*2-two_s1);two_Sa<=sjk*2+two_s1;two_Sa+=2)
						{
							for(int L=abs(lx-ly);L<=lx+ly;L++)
							{
								if(abs(L*2-two_Sa)>two_J || 2*L+two_Sa<two_J || L>0)continue;
								cout<<aNc+1<<"   "<<q<<"   "<<lx<<"   "<<ly<<"   "<<L<<"   "<<sjk<<"   "<<two_Sa<<"   "<<tjk<<endl;
								anglemomentlist[aNc][0]=q;
								anglemomentlist[aNc][1]=lx;
								anglemomentlist[aNc][2]=ly;
								anglemomentlist[aNc][3]=L;
								anglemomentlist[aNc][4]=sjk;
								anglemomentlist[aNc][5]=two_Sa;
								anglemomentlist[aNc][6]=tjk;

								aNc++;
							}
						}
					}
				}
			}
		}
	}
	const int Hdim=aNc*nmax;
	cout<<"Hdim: "<<Hdim<<endl;

	//Variation to solve Schrodinger equation and get the parameter c and solution A

	int H_row=0;//Point to the row of H
	int H_col=0;//Point to the column of H
	double Hmin;
	double Hmin_Vec[Hdim];
  memset(Hmin_Vec,0,sizeof(Hmin_Vec));
	
	int nstart=2;
	int nstop=nstart+nmax-1;
	
	int Nc=atoi(argv[1]);
	TMatrixDSym H(Hdim);
	//H is the matrix form of Hamiltonian
	c = c0 + c00*(Nc - Ncmax/2 + 1);//Giving c value in this precision

	if(c <= 0)
	{
		cerr<<"In Nc="<<Nc<<", c <= 0, c should be postive!"<<endl;
		return 1;
	}
	else if(c>10000)
	{
		cerr<<"In Nc="<<Nc<<", c > 10000, out of range"<<endl;
		return 1;
	}
	c*=fmmev;
	int process=0;
	for(int base_n=nstart; base_n<=nstop; base_n++)
	{
		if(process<(base_n-nstart+1)*10/nmax)
		{
			process=(base_n-nstart+1)*10/nmax;
			cout<<"process: "<<process*10<<"%"<<endl;
		}
		for(int iaNc1=0;iaNc1<aNc;iaNc1++)
		{
			int q1=anglemomentlist[iaNc1][0];
			int lx1=anglemomentlist[iaNc1][1];
			int ly1=anglemomentlist[iaNc1][2];
			int L1=anglemomentlist[iaNc1][3];
			int sjk1=anglemomentlist[iaNc1][4];
			int two_Sa1=anglemomentlist[iaNc1][5];
			int tjk1=anglemomentlist[iaNc1][6];

			for(int base_m=nstart; base_m<=nstop; base_m++)
			{
				for(int iaNc2=0;iaNc2<aNc;iaNc2++)
				{
					int q2=anglemomentlist[iaNc2][0];
					int lx2=anglemomentlist[iaNc2][1];
					int ly2=anglemomentlist[iaNc2][2];
					int L2=anglemomentlist[iaNc2][3];
					int sjk2=anglemomentlist[iaNc2][4];
					int two_Sa2=anglemomentlist[iaNc2][5];
					int tjk2=anglemomentlist[iaNc2][6];

					if(H_row<H_col)continue;
					//For the H matrix is a symmetric matrix, so the matrix elements just need to be calculated until diagonal element

					//c=5.06;


					TF2 *fun_integrand_H=new TF2("fun_integrand_H",integrand_H,0,50,0,Pi/2,16);
					double par[]={base_n,q1,lx1,ly1,L1,sjk1,two_Sa1,tjk1,base_m,q2,lx2,ly2,L2,sjk2,two_Sa2,tjk2};

					//double xtest[]={2,0.67};
					//cout<<"c="<<c<<", integrand_H({2,0.67},par)="<<integrand_H(xtest,par)<<endl;
					//return 0;

					fun_integrand_H->SetParameters(par);
					H(H_row,H_col)=fun_integrand_H->TF2::Integral(0.,50,0.,Pi*0.5,1.e-300);//calculate the matrix elements of H

					if(H_row!=H_col) H(H_col,H_row) = H(H_row,H_col);
					delete fun_integrand_H;

					//cout<<q1<<" "<<lx1<<" "<<ly1<<" "<<q2<<" "<<lx2<<" "<<ly2<<" "<<endl;
					//cout<<"H("<<H_row<<","<<H_col<<"): "<<H(H_row,H_col)/fmmev<<endl;
					//return 0;
					H_col++;
				}//iaNc2 loop;
			}//base_m loop, index of orthbasis
			H_row++;
			H_col=0;//reset H_col.
		}//iaNc1 loop;
	}//base_n loop, index of orthbasis

	//H.Print();
	//the H matrix is complete for this c
	H_row=0;//reset H_row.

	TMatrixDSymEigen HEigen(H);//TMatrixDSym need to transform into TMatrixDSymEigen to calculate the Eigen values and vectors
	TVectorD HVa=HEigen.GetEigenValues();//HVa is the vector which consists of eigen values
	TMatrixD HVe=HEigen.GetEigenVectors();
	//HVe is the matrix which consists of eigen vectors, and HVe(*,i)(the ith column of HVe) corresponds to HVa(i)
	Hmin=HVa(Hdim-1);
	for(int i=0;i<Hdim;i++)Hmin_Vec[i]=TMatrixDColumn(HVe,Hdim-1)[i];

	cout<<endl;
	cout<<"///////////////////////////////////////"<<endl;
	cout<<"//Variation done"<<endl;
	cout<<"///////////////////////////////////////"<<endl;
	cout<<endl;
	cout<<"///////////////////////////////////////"<<endl;
	cout<<"//Result"<<endl;
	cout<<"///////////////////////////////////////"<<endl;
	cout<<endl;
	cout<<nuclear<<endl;

	cout<<"Variation parameter c="<<c/fmmev<<" MeV, base state energy: "<<Hmin/fmmev<<" MeV"<<endl;
	cout<<setiosflags(ios::fixed)<<setprecision(5)<<"Eigen vector: "<<endl;
	double normlize=0;
	for(int i=0;i<Hdim;i++)normlize+=Hmin_Vec[i]*Hmin_Vec[i];
	for(int i=0;i<Hdim;i++){Hmin_Vec[i]/=normlize;cout<<Hmin_Vec[i]<<" "<<endl;}
	cout<<endl;
	cout<<"///////////////////////////////////////////////////////////////"<<endl;
	
	char txtsave[100];
	sprintf(txtsave,"%s/c_%.4f_energy_%.6f.txt",nuclear,c/fmmev,Hmin/fmmev);
	ofstream f;
	f.open(txtsave);
	for(int i=0;i<Hdim;i++){f<<Hmin_Vec[i]<<endl;}
	f.close();


	//HEigen.Delete();
	//HVa.Print();

	return 0;
}





//////////////////////////////////////////////////////////////

double V23(int sjk, double r, double alpha)
{
	double V=0;
	if(sjk==1)V=VNN_01(r23(r,alpha));
	if(sjk==0)V=VNN_10(r23(r,alpha));
	if(sjk==2)V=VNOmega(r23(r,alpha));
	return V;
}

double V12(int sjk, double r, double alpha)
{
	double V=0;
	if(sjk==1)V=VNN_01(r12(r,alpha));
	if(sjk==0)V=VNN_10(r12(r,alpha));
	if(sjk==2)V=VNOmega(r12(r,alpha));
	return V;
}

double V13(int sjk, double r, double alpha)
{
	double V=0;
	if(sjk==1)V=VNN_01(r13(r,alpha));
	if(sjk==0)V=VNN_10(r13(r,alpha));
	if(sjk==2)V=VNOmega(r13(r,alpha));
	return V;
}

double Vc12(double r, double alpha)
{
	return Vc(r12(r,alpha));
}

double Vc13(double r, double alpha)
{
	return Vc(r13(r,alpha));
}

double Vc23(double r, double alpha)
{
	return Vc(r23(r,alpha));
}

/*double orthbasis_radial(double r, int n)//number l and variation parameter c is defined as global.
	{
	if(n<l+1)return 0;
	double result_basis=TMath::Sqrt(TMath::Power(2*c/n,3)*TMath::Factorial(-1 - l + n)/(2*n*TMath::Factorial(l + n)))*exp(-c*r/n)*TMath::Power(2*c*r/n, l)
 *assoc_laguerre(-1-l+n, 1+2*l, 2*c*r/n)*r;
//The bases are given by hydrogen, it can be write as |nl> in the annotation.
// orthbasis[r, n, l] = <r|nl> *r is the reduced radial wave function. 
// The unit of reduced radial wave function is always [r]^(-1/2) in n-dim, that's why choose it as base.
//  Attention : the bases with different l are not orthogonal.
//  And the c is the variational parameter.
return result_basis;
}*/

/*double orthbasis_alpha(int q, int lx, int ly, double alpha)
	{
//Kievsky, A., Viviani, M., & Rosati, S. (1993). Nuclear Physics A, 551(2) 241;
//Kovalchuk, V. I. (2014). International Journal of Modern Physics E, 23(11), 1450069.
//Eigen function of hyper angle momentum with wigen value K=2*q+lx+ly;
double Nk=sqrt(2*TMath::Factorial(q)*(2*q+lx+ly+2)*TMath::Factorial(q+lx+ly+1)/(tgamma(q+lx+1.5)*tgamma(q+ly+1.5)));
return Nk*pow(cos(alpha),lx)*pow(sin(alpha),ly)*JacobiP(q,ly+0.5,lx+0.5,cos(2*alpha));
}*/

/*double fun_orthbasis_radial(double *x, double *par)
	{
	double r=x[0];
	int n=par[0];
	return orthbasis_radial(r,n);
	}*/

/*double fun_orthbasis(double *x, double *par)
	{
	double r=x[0];
	double alpha=x[1];
	int n=par[0];
	int q=par[1];
	int lx=par[2];
	int ly=par[3];
	return orthbasis_radial(r,n)*orthbasis_alpha(q,lx,ly,alpha);
}*/

double integrand_H(double *x, double *Par)
{
	double r=x[0];//hyper spherical radius
	double alpha=x[1];//hyper spherical angle

	int base_n=Par[0];//n is index of R_n (radial basis)
	int q1=Par[1];//q is index of alpha basis
	int lx1=Par[2];//lx is the orbital angle momentum of Jacobian coordinate axix x;
	int ly1=Par[3];//ly is the orbital angle momentum of Jacobian coordinate axix y;
	int L1=Par[4];//Total orbital momentum;
	int sjk1=Par[5];//sjk is the spin of pair jk;
	int two_Sa1=Par[6];//Sa is the total spin of three particle; two_Sa1=2*Sa1;
	int tjk1=Par[7];//tjk is the iospin of pair jk;
	int base_m=Par[8];
	int q2=Par[9];
	int lx2=Par[10];
	int ly2=Par[11];
	int L2=Par[12];
	int sjk2=Par[13];
	int two_Sa2=Par[14];
	int tjk2=Par[15];
	//the matrix H is expand as
	//  <n q1 lx1 ly1 L1 sjk1 Sa1 J tjk1 T|H|m q2 lx2 ly2 L2 sjk2 Saj J tjk2 T> 
	//i and j is the coordinate index.
	
	if(two_Sa1!=two_Sa2 || sjk1!=sjk2 || tjk1!=tjk2 || L1!=L2)return 0;//Hamilation is diag of spin and iospin;

	int K1=q1*2+lx1+ly1;//hyper angle momentum of coordinate i;
	int K2=q2*2+lx2+ly2;//hyper angle momentum of coordinate j;

  double elem_volume=pow(cos(alpha)*sin(alpha),2);

	double base1=0;
  if(lx1==lx2 && ly1==ly2 && q1==q2)
	{
		TF1 *basis_radial_m=new TF1("basis_m",fun_orthbasis_radial,0,50,1);
		basis_radial_m->SetParameter(0,base_m);
		double d2basis_radial_m=basis_radial_m->Derivative2(r);
		base1=orthbasis_radial(r,base_n)*(-d2basis_radial_m/(2.*Q)+(K2+1.5)*(K2+2.5)/(r*r*2.*Q)*orthbasis_radial(r,base_m))*16./Pi;
		delete basis_radial_m;
	}

  double base2=0;
	if(lx1==lx2 && ly1==ly2)
	{
		base2=V23(sjk1,r,alpha)*orthbasis_radial(r,base_n)*orthbasis_radial(r,base_m)*orthbasis_alpha(q1,lx1,ly1,alpha)*orthbasis_alpha(q2,lx2,ly2,alpha);
#ifdef coulomb23
		base2+=coulomb23*Vc23(r,alpha)*orthbasis_radial(r,base_n)*orthbasis_radial(r,base_m)*orthbasis_alpha(q1,lx1,ly1,alpha)*orthbasis_alpha(q2,lx2,ly2,alpha);
#endif
	}

	double base3=0;
	for(int lxs=0;lxs<=min(K1,K2);lxs++)
	{
		for(int lys=0;lys<=min(K1,K2);lys++)
		{
			if(lxs+lys<L1 || abs(lxs-lys)>L1 || lxs+lys>K1 || lxs+lys>K2 || (K1-lxs-lys)%2!=0 || (K2-lxs-lys)%2!=0)continue;
			int qs1=(K1-lxs-lys)/2;
			int qs2=(K2-lxs-lys)/2;

#ifdef coulomb12
			base3+=RRM(1,2,q1,lx1,ly1,qs1,lxs)*RRM(2,1,qs2,lxs,lys,q2,lx2)*Vc12(r,alpha)*coulomb12/2.*
				       orthbasis_radial(r,base_n)*orthbasis_radial(r,base_m)*
								 orthbasis_alpha(qs1,lxs,lys,alpha)*orthbasis(qs2,lxs,lys,alpha);
#endif
#ifdef coulomb23
			base3+=RRM(1,3,q1,lx1,ly1,qs1,lxs)*RRM(3,1,qs2,lxs,lys,q2,lx2)*Vc13(r,alpha)*coulomb23/2.*
				       orthbasis_radial(r,base_n)*orthbasis_radial(r,base_m)*
								 orthbasis_alpha(qs1,lxs,lys,alpha)*orthbasis(qs2,lxs,lys,alpha);
#endif

#ifdef NNOmega
			base3+=RRM(1,2,q1,lx1,ly1,qs1,lxs)*RRM(2,1,qs2,lxs,lys,q2,lx2)*V12(2,r,alpha)*
			//base3+=RRM(1,2,q1,lx1,ly1,qs1,lxs,lys)*RRM(1,2,q2,lx2,ly2,qs2,lxs,lys)*V12(2,r,alpha)*
				       orthbasis_radial(r,base_n)*orthbasis_radial(r,base_m)*
								 orthbasis_alpha(qs1,lxs,lys,alpha)*orthbasis(qs2,lxs,lys,alpha);
#else
			for(int two_sjks=abs(two_s1-two_s3);two_sjks<=(two_s1+two_s3);two_sjks+=2)
			{
				if(two_sjks+two_s2<two_Sa1 || abs(two_sjks-two_s2)>two_Sa1)continue;
				base3+=RRM(1,2,q1,lx1,ly1,qs1,lxs)*RRM(2,1,qs2,lxs,lys,q2,lx2)*V12(two_sjks/2,r,alpha)*
				//base3+=RRM(1,2,q1,lx1,ly1,qs1,lxs)*RRM(1,2,q2,lx2,ly2,qs2,lxs)*V12(two_sjks/2,r,alpha)*
								Wcoeff(two_s1,two_s3,two_sjks,two_s2,two_Sa1,2*sjk1)*Wcoeff(two_s2,two_s3,2*sjk2,two_s1,two_Sa2,two_sjks)*
				          orthbasis_radial(r,base_n)*orthbasis_radial(r,base_m)*
								    orthbasis_alpha(qs1,lxs,lys,alpha)*orthbasis_alpha(qs2,lxs,lys,alpha);
			}
#endif
		}
	}

	/*double base4=0;
	for(int lxs=0;lxs<=min(K1,K2);lxs++)
	{
		for(int lys=0;lys<=min(K1,K2);lys++)
		{
			if(lxs+lys<L1 || abs(lxs-lys)>L1 || lxs+lys>K1 || lxs+lys>K2 || (K1-lxs-lys)%2!=0 || (K2-lxs-lys)%2!=0)continue;
			int qs1=(K1-lxs-lys)/2;
			int qs2=(K2-lxs-lys)/2;
#ifdef NNOmega
			if(sjk1!=sjk2)continue;
			//base4+=RRM(1,3,q1,lx1,ly1,qs1,lxs)*RRM(3,1,qs2,lxs,lys,q2,lx2)*V13(2,r,alpha)*
			base4+=RRM(1,3,q1,lx1,ly1,qs1,lxs)*RRM(1,3,q2,lx2,ly2,qs2,lxs)*V13(2,r,alpha)*
				       orthbasis_radial(r,base_n)*orthbasis_radial(r,base_m)*
								 orthbasis_alpha(qs1,lxs,lys,alpha)*orthbasis_alpha(qs2,lxs,lys,alpha);
#else
			for(int two_sjks=abs(two_s1-two_s2);two_sjks<=(two_s1+two_s2);two_sjks+=2)
			{
				if(two_sjks+two_s3<two_Sa1 || abs(two_sjks-two_s3)>two_Sa1)continue;
				//base4+=RRM(1,3,q1,lx1,ly1,qs1,lxs)*RRM(3,1,qs2,lxs,lys,q2,lx2)*V13(two_sjks/2,r,alpha)*
				base4+=RRM(1,3,q1,lx1,ly1,qs1,lxs)*RRM(1,3,q2,lx2,ly2,qs2,lxs)*V13(two_sjks/2,r,alpha)*
								Wcoeff(two_s1,two_s2,two_sjks,two_s3,two_Sa1,2*sjk2)*Wcoeff(two_s3,two_s2,2*sjk1,two_s1,two_Sa1,two_sjks)*
				          orthbasis_radial(r,base_n)*orthbasis_radial(r,base_m)*
								    orthbasis_alpha(qs1,lxs,lys,alpha)*orthbasis_alpha(qs2,lxs,lys,alpha);
			}
#endif
		}
	}*/
	return (2*base3+base1+base2)*elem_volume;
}

