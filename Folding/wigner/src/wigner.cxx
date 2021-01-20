#include <head.h>

using namespace std;

//infomation
double c;
int l;
int nmax;
int nstart;
int nstop;
double par[100]={0};

//plot infomation
#define scale
#ifdef scale
double aBorhr=1;//fm
int Nr=50;
double rmin=-8;
double rmax=2;
double dr=(rmax-rmin)/Nr;

double ku=1/aBorhr/fmmev;//MeV
int Np=50;
double pmin=-8;
double pmax=1;
double dp=(pmax-pmin)/Np;
#else
int Nr=50;
double rmin=0;//fm
double rmax=10;//fm
double dr=(rmax-rmin)/Nr;//fm

int Np=50;
double pmin=0;//MeV
double pmax=1000;//MeV
double dp=(pmax-pmin)/Np;//MeV
#endif
int Ntheta=24;
double thetamin=0;
double thetamax=Pi;
double dtheta=(thetamax-thetamin)/Ntheta;

//functions
#include <jacobi_polynomial.hpp>
#include "../../include/orthbasis.h"
double wignerdst_Integrand(double *x, double *par);

int main(int argc, char *argv[])
{
	char *file=argv[1];
	char filepath[500];
	sprintf(filepath,"./stateinfo/%s.txt",file);
	ifstream read;
	read.open(filepath);
	if(!(read.is_open()))
	{
		cerr<<"open "<<filepath<<" fail."<<endl;
		return 1;
	}
	read>>c;
	read>>l;
	read>>nmax;
	read>>nstart;
	nstop=nstart+nmax-1;
	for(int i=nstart;i<=nstop;i++)read>>par[i];

	/*double paratest[]={2,2};
		double xtest[]={2,0.2,0.8};
		cout<<wignerdst_Integrand(xtest,paratest)<<endl;
		return 0;*/

	TF2 *wigdstInte=new TF2("wignerdst_Integrand",wignerdst_Integrand,0,50,0,Pi,3);

	TH2D wigdst("wigdst","wigner density plot",Nr,rmin,rmax,Np,pmin,pmax);

	double process=0;

	//double parr[]={2.,2.*fmmev};
	//double x[]={2.,3.,1.};
	//cout<<wignerdst_Integrand(x,parr)<<endl;
	//return 0;

	int thetal=atoi(argv[2]);
	double theta=thetal*dtheta+thetamin;

	for(int i=1;i<=Nr;i++)
	{
#ifdef scale
		double ri=exp(rmin+i*dr)*aBorhr;
#else
		double ri=rmin+i*dr;
#endif

		if(floor(i*10/Nr)>process)
		{
			process=floor(i*10/Nr);
			cout<<"precess: "<<process*10<<"%"<<endl;
		}
		//double resMax=0;

		for(int j=1;j<=Nr;j++)
		{
			//clock_t start,end;
			//start=clock();
#ifdef scale
			double pj=exp(pmin+dp*j)*ku;
#else
			double pj=pmin+dp*j;
#endif
			pj*=fmmev;

			wigdstInte->SetParameters(ri,pj,theta);
			double res=wigdstInte->Integral(0.,50.,0.,Pi,1.e-20);

			wigdst.SetBinContent(i,j,l,res);
		}
	}

	char savepath[100];
#ifdef scale
	sprintf(savepath,"./wigdst/%s/%s_thetal_%d-scaled.root",file,file,thetal);
#else
	sprintf(savepath,"./wigdst/%s/%s_thetal_%d-scaled.root",file,file,thetal);
#endif
	TFile save(savepath,"recreate");
	if(save.IsZombie())
	{
		cerr<<"file "<<savepath<<" create error"<<endl;
		return 1;
	}
	wigdst.Write();
	save.Close();
	return 0;
}

double wignerdst_Integrand(double *x,double *par)//no angle relate
{
	double y=x[0];
	double theta_y=x[1];//the angle between y and r

	double r=par[0];
	double p=par[1];
	double theta=par[2];//the angle between r and p

	double rp=sqrt(r*r+y*y/4.-r*y*cos(theta_y));//|vec_r+vec_y/2|
	double rn=sqrt(r*r+y*y/4.+r*y*cos(theta_y));//|vec_r-vec_y/2|

	double res=0;
	for(int n1=nstart;n1<=nstop;n1++)
	{
		for(int n2=nstart;n2<nstop;n2++)
		{
			double res1=par[n1]*par[n2]*orthbasis_radial(rp,n1)/pow(rp,1./2.)*orthbasis_radial(rn,n2)/pow(rn,1./2.)*1./(4.*Pi);
			res+=res1;
		}
	}

	double res2=TMath::BesselJ0(p*y*sin(theta)*sin(theta_y))*2*Pi;
	//equal to Integreat[exp(-I*p*sin(theta)*y*sin(theta_y)*cos(phi_y)),{phi_y,0,2*Pi}]
	res*=res2;

	double res3=cos(p*cos(theta)*y*cos(theta_y));//real part of exp(-I*p*cos(theta)*y*cos(theta_y))
	res*=res3;

	res*=pow(y,2)*sin(theta_y)*2*Pi;//volume element of hyper-spherical coordinate (y,theta_y,phi_y);

	return res;
}
