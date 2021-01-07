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

//#define E200GeV
#define E2_76TeV

//choice a system, choice H_3, He_3 or NNOmega
//#include "../../include/H-3.h"
//#include "../../include/He-3.h"
#include "../../include/pnOmega.h"
//#include "../../include/nnOmega.h"
//#include "../../include/ppOmega.h"
//#include "../../include/nOmegaOmega.h"
//#include "../../include/pOmegaOmega.h"

double r23(double r, double alpha);
double r13(double r, double alpha);
double r12(double r, double alpha);
double V23(int sjk, double r, double alpha);
double V13(int sjk, double r, double alpha);
double V12(int sjk, double r, double alpha);
double Vc23(double r, double alpha);
double Vc13(double r, double alpha);
double Vc12(double r, double alpha);
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
			  It means, in this loop, c varies from c0-c00/10^loop*(Ncmax/2) to c0+c00/10^loop*(Ncmax/2)*/
int eprecision=0;
double precision=pow(10,-eprecision);// determine the final precision of c.
int Nprecisionloop=TMath::Log10(c00/precision);/*The Nprecisionloop determine loop times to reach the precision defined before*/
void Hmin_reset();

//Initialization of variation done

///////////////////////////////////////
//Initialization done
///////////////////////////////////////

//////////////////////////////////////////////////
//////////////////////////////////////////////////

double integrand_H(double *x, double *Par);

int main()
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

	int qmax=8;
	int lmax=4;

	cout<<"Nc  q   lx  ly  L   sjk Sa  tjk"<<endl;

	int anglemomentlist[100][7];
	int aNc=0;
	for(int q=0;q<=qmax;q++)
	{
		for(int lx=0;lx<=lmax;lx++)
		{
			for(int ly=0;ly<=lmax;ly++)
			{
				if(pow(-1,lx+ly)!=Ptotal)continue;
				for(int sjk=sjkmin;sjk<=sjkmax;sjk++)
				{
					for(int tjk=tjkmin;tjk<=tjkmax;tjk++)
					{
						if(pow(-1,sjk+tjk+lx)!=-1)continue;
						for(int two_Sa=abs(sjk*2-two_s1);two_Sa<=sjk*2+two_s1;two_Sa+=2)
						{
							for(int L=abs(lx-ly);L<=lx+ly;L++)
							{
								if(abs(L*2-two_Sa)>two_J || 2*L+two_Sa<two_J || L>0)continue;
								cout<<aNc+1<<"   "<<q<<"   "<<lx<<"   "<<ly<<"   "<<L<<"   "<<sjk<<"   "<<two_Sa<<"/2"<<"   "<<tjk<<endl;
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
	double Hmin[2][Ncmax];
	double Hmin_Vec[Ncmax][Hdim];
  memset(Hmin,0,sizeof(Hmin));
  memset(Hmin_Vec,0,sizeof(Hmin_Vec));
  /*An array record the value of miniest eigenvalue in this precisionloop and its eigenvector and c. 
    The order of data is: c, eigenvalue, Nc(position of the elements). And Hmin_Vec store the eigenvectors*/
	
	int nstart=2;
	int nstop=nstart+nmax-1;
	
	for(int precisionloop=0; precisionloop<=Nprecisionloop; precisionloop++)//The precision, in one loop, is c00/10^precisionloop
	{
		for(int Nc=0; Nc<Ncmax; Nc++)//variation of c in this precisionloop
		{
			TMatrixDSym H(Hdim);
				//H is the matrix form of Hamiltonian
			c = c0 + c00/TMath::Power(10,precisionloop)*(Nc - Ncmax/2 + 1);//Giving c value in this precision
			if(c >= 0 && c<=10000)
			{
				Hmin[0][Nc]=c;//Record the c value in this precisionloop
				c*=fmmev;
				for(int base_n=nstart; base_n<=nstop; base_n++)
				{
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
				Hmin[1][Nc]=HVa(Hdim-1);
				for(int i=0;i<Hdim;i++){
					Hmin_Vec[Nc][i]=TMatrixDColumn(HVe,Hdim-1)[i];
				}

				cout<<"Variation parameter c="<<c/fmmev<<" MeV, base state energy: "<<Hmin[1][Nc]/fmmev<<" MeV"<<endl;

				//HEigen.Delete();
				HVa.Clear();
				HVe.Clear();
				//HVa.Print();
				//return 0;

			}//This condition is c must lager than 0
			H.Clear();
		}//Nc loop. This loop is same as variation of c, but find optimal solution by scan in every precision.
		double min_Hmin=99999;
		int min_Nc=-1;
		for(int i=0;i<Ncmax;i++)
		{
			if(Hmin[1][i]<min_Hmin)
			{
				min_Hmin=Hmin[1][i];
				min_Nc=i;
			}
		}
		c0=Hmin[0][min_Nc];
		double process=100*(precisionloop+1)/(1+Nprecisionloop);
		cout<<setiosflags(ios::fixed)<<setprecision(2)<<"Variation process: "<<process<<"%"<<endl;;
	}//Precisionloop is used to give c a precision, the precision of c is c00/(10^Nprecisionloop) in every loop.

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

	double min_Hmin=99999;
	int min_Nc=-1;
	for(int i=0;i<Ncmax;i++)
	{
		if(Hmin[1][i]<min_Hmin)
		{
			min_Hmin=Hmin[1][i];
			min_Nc=i;
		}
	}
	c=Hmin[0][min_Nc];
	cout<<setiosflags(ios::fixed)<<setprecision(eprecision)<<"Variation parameter c = "<<c<<" fm^-1"<<endl;
	cout<<setiosflags(ios::fixed)<<setprecision(5)<<"Eigen vector: ";
	double normlize=0;
	for(int i=0;i<Hdim;i++)normlize+=Hmin_Vec[min_Nc][i]*Hmin_Vec[min_Nc][i];
	for(int i=0;i<Hdim;i++){Hmin_Vec[min_Nc][i]/=normlize;cout<<Hmin_Vec[min_Nc][i]<<" ";}
	cout<<endl;
	cout<<"Binding energy: "<<min_Hmin<<" MeV"<<endl;
	/*TF1 *probability=new TF1("probability",prob,0,50,nmax);
	probability->SetParameters(Hmin_Vec[min_Nc]);
	double max_x=probability->GetMaximumX(0,10);
	delete probability;
	cout<<"Most probable radius: "<<max_x<<" fm"<<endl;*/
	cout<<endl;
	cout<<"///////////////////////////////////////////////////////////////"<<endl;

	//save.close();
	return 0;
}





//////////////////////////////////////////////////////////////

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
								 orthbasis_alpha(qs1,lxs,lys,alpha)*orthbasis_alpha(qs2,lxs,lys,alpha);
#endif
#ifdef coulomb13
			base3+=RRM(1,3,q1,lx1,ly1,qs1,lxs)*RRM(3,1,qs2,lxs,lys,q2,lx2)*Vc13(r,alpha)*coulomb13/2.*
				       orthbasis_radial(r,base_n)*orthbasis_radial(r,base_m)*
								 orthbasis_alpha(qs1,lxs,lys,alpha)*orthbasis_alpha(qs2,lxs,lys,alpha);
#endif

#ifdef two_NOmega_interaction
			base3+=RRM(1,2,q1,lx1,ly1,qs1,lxs)*RRM(2,1,qs2,lxs,lys,q2,lx2)*V12(2,r,alpha)*
				       orthbasis_radial(r,base_n)*orthbasis_radial(r,base_m)*
								 orthbasis_alpha(qs1,lxs,lys,alpha)*orthbasis_alpha(qs2,lxs,lys,alpha);
#else
			for(int two_sjks=abs(two_s1-two_s3);two_sjks<=(two_s1+two_s3);two_sjks+=2)
			{
				if(two_sjks+two_s2<two_Sa1 || abs(two_sjks-two_s2)>two_Sa1)continue;
				base3+=RRM(1,2,q1,lx1,ly1,qs1,lxs)*RRM(2,1,qs2,lxs,lys,q2,lx2)*V12(two_sjks/2,r,alpha)*
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
#ifdef two_NOmega_interaction
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

