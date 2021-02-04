#include <head.h>

using namespace std;

//infomation storage
TH1D *stateinfo;
TH1D *function_parameter;
TH2I *anglemomentlist;

//infomation
int qmax;
int lmax;
int Lmax;
int angnum;
int nmax;
int nstart;
double c;

int readinfo(TFile *read,char *file);

int Nselect;
double state_select[100][2];
double norm=0;
int select();

//plot infomation
#define scale
#ifdef scale
double aBorhr=1;//fm
const int Nr=50;
double rmin=-8;//fm
double rmax=2;//fm
double dr=(rmax-rmin)/Nr;

double ku=1/aBorhr/fmmev;//MeV
const int Np=50;
double pmin=-8;//MeV
double pmax=1;//MeV
double dp=(pmax-pmin)/Np;

#else
const int Nr=50;
double rmin=0;//fm
double rmax=10;//fm
double dr=(rmax-rmin)/Nr;//fm

const int Np=50;
double pmin=0;//MeV
double pmax=1000;//MeV
double dp=(pmax-pmin)/Np;//MeV
#endif
const int Ntheta=24;
double thetamin=0;
double thetamax=Pi;
double dtheta=(thetamax-thetamin)/Ntheta;

//functions
int l=1;
#include <jacobi_polynomial.hpp>
#include "../../include/orthbasis.h"
double wignerdst_Integrand(double *x, double *par);

int main(int argc, char *argv[])
{
	char *file=argv[1];
	char filepath[520];
	sprintf(filepath,"./stateinfo/%s.root",file);
	TFile *read=new TFile(filepath,"read");
	if(readinfo(read,filepath))return 1;

	Nselect=select();
	if(Nselect==-1)
	{
		cout<<"angle relate"<<endl;
		return 0;
	}
	if(Nselect==-2)
	{
		cout<<"select is not enough(norm<0.94)"<<endl;
		return 0;
	}
	
	/*double paratest[]={2,2};
	double xtest[]={2,0.2,0.8};
	cout<<wignerdst_Integrand(xtest,paratest)<<endl;
	return 0;*/

	TF2 *wigdstInte=new TF2("wignerdst_Integrand",wignerdst_Integrand,0,50,0,Pi,3);

	int tl=atof(argv[2])+1;
	double thetal=thetamin+dtheta*(tl-0.5);
	
	char plotname[500];
	sprintf(plotname,"wigdst_thetal_%d",tl);
	TH2D wigdst(plotname,"wigner density plot",Nr,rmin,rmax,Np,pmin,pmax);

	double process=0;

	//double parr[]={2.,2.*fmmev};
	//double x[]={2.,3.,1.};
	//cout<<wignerdst_Integrand(x,parr)<<endl;
	//return 0;

	for(int i=1;i<=Nr;i++)
	{
		if(floor(i*10/Nr)>process)
		{
			process=floor(i*10/Nr);
			cout<<"precess: "<<process*10<<"%"<<endl;
		}

		double ri=wigdst.GetXaxis()->GetBinCenter(i);
#ifdef scale
		ri=exp(ri)*aBorhr;
#endif

		//double resMax=0;

		for(int j=1;j<=Nr;j++)
		{
			double pj=wigdst.GetYaxis()->GetBinCenter(j);
#ifdef scale
			pj=exp(pj)*ku;
#endif
			pj*=fmmev;

			/*double paratest[]={ 0.301194,0.261846,thetal};
				double xtest[]={2,0.2};
				cout<<wignerdst_Integrand(xtest,paratest)<<endl;
				return 0;*/

			wigdstInte->SetParameters(ri,pj,thetal);
			double res=wigdstInte->Integral(0.,50.,0.,Pi,1.e-20);

			//if(i>30&&i<40&&j>30 &&j<40 && res<0.0001)cout<<i<<" "<<j<<" "<<ri<<" "<<pj<<res<<endl;

			/*double absres=abs(res);
				if(absres<resMax/1000)
				{
				break;
				}
				else if(absres>resMax)
				{
				resMax=absres;
				}*/

			wigdst.SetBinContent(i,j,res);
		}
	}

	char savepath[500];
#ifdef scale
	sprintf(savepath,"./wigdst/%s/%s_thetal_%d-scaled.root",file,file,tl);
#else
	sprintf(savepath,"./wigdst/%s/%s_thetal_%d.root",file,file,tl);
#endif
	TFile save(savepath,"recreate");
	if(save.IsZombie())
	{
		cerr<<"file "<<file<<" already exist"<<endl;
		return 1;
	}
	wigdst.Write();
	save.Close();
	delete read;
	delete wigdstInte;
	return 0;
}

int readinfo(TFile *read,char *file)
{
	if(read->IsZombie())
	{
		cerr<<"file "<<file<<" open error"<<endl;
		return 1;
	}
	stateinfo=(TH1D*)read->Get("stateinfo");
	function_parameter=(TH1D*)read->Get("function_parameter");
	anglemomentlist=(TH2I*)read->Get("anglemomentlist");

	qmax=(int)stateinfo->GetBinContent(1);
	lmax=(int)stateinfo->GetBinContent(2);
	Lmax=(int)stateinfo->GetBinContent(3);
	angnum=(int)stateinfo->GetBinContent(4);
	nmax=(int)stateinfo->GetBinContent(5);
	nstart=(int)stateinfo->GetBinContent(6);
	c=stateinfo->GetBinContent(7);
	c*=fmmev;

	return 0;
}

int select()
{
	int i=0;
	for(int iNc1=1;iNc1<=angnum;iNc1++)
	{
		int q1=anglemomentlist->GetBinContent(iNc1,1);
		int lx1=anglemomentlist->GetBinContent(iNc1,2);
		int ly1=anglemomentlist->GetBinContent(iNc1,3);
		int L1=anglemomentlist->GetBinContent(iNc1,4);
		int sjk1=anglemomentlist->GetBinContent(iNc1,5);
		int two_Sa1=anglemomentlist->GetBinContent(iNc1,6);
		int tjk1=anglemomentlist->GetBinContent(iNc1,7);

		for(int n1=nstart;n1<=nstart+nmax;n1++)
		{
			double par1=function_parameter->GetBinContent(iNc1+(n1-nstart)*angnum);
			if(abs(par1)<=0.11)continue;
			//cout<<i<<" "<<n1<<" "<<q1<<" "<<lx1<<" "<<ly1<<endl;
			if(q1!=0||lx1!=0||ly1!=0)return -1;
			state_select[i][0]=n1;
			state_select[i][1]=par1;
			//cout<<par1<<endl;
			i++;
			norm+=par1*par1;
		}
	}
	//cout<<norm<<endl;
	if(norm<0.94)return -2;
	return i;
}

double wignerdst_Integrand(double *x,double *par)//no angle relate
{
	double y=x[0];
	double theta_y1=x[1];//the angle between y and r

	double r=par[0];
	double p=par[1];
	double theta=par[2];//the angle between r and p

	double rp=sqrt(r*r+y*y/4.-r*y*cos(theta_y1));//|vec_r+vec_y/2|
	double rn=sqrt(r*r+y*y/4.+r*y*cos(theta_y1));//|vec_r-vec_y/2|
	
	double res=0;
	for(int i1=0;i1<Nselect;i1++)
	{
		int n1=(int)state_select[i1][0];
		double par1=state_select[i1][1];

		for(int i2=0;i2<Nselect;i2++)
		{
			int n2=(int)state_select[i2][0];
			double par2=state_select[i2][1];

			double res1=par1*par2*orthbasis_radial(rp,n1)/pow(rp,5./2.)*orthbasis_radial(rn,n2)/pow(rn,5./2.)*16./Pi*1./(4.*Pi)*1./(4.*Pi);
			//u(rp)*u(rp)*Y_00(theta_rp_J1,phi_rp_J1)*Y_00(theta_rp_J2,phi_rp_J2)*P_0,0+1/2,0+1/2(alpha_rp)*conjugate_Y_00(theta_rn_J1,phi_rn_J1)*conjugate_Y_00(theta_rn_J2,phi_rn_J2)*P_0,0+1/2,0+1/2(alpha_rn)
			//Ji represent the jacobi-transformed coordinate

			//cout<<"res1: "<<res1<<endl;
			res+=res1;

		}
	}
	//psi(vec_r+vec_y/2)*conjugate_psi(vec_r-vec_y/2)
	//cout<<norm<<endl;
	res/=norm;

	//exp(-I*vec_p*vec_y)=exp(-I*p*sin(theta)*y*sin(theta_y1)*cos(theta_y2)-I*p*cos(theta)*y*cos(theta_y1))
	double res2=4*(sin(p*sin(theta)*y*sin(theta_y1))-p*sin(theta)*y*sin(theta_y1)*cos(p*sin(theta)*y*sin(theta_y1)))/pow(p*sin(theta)*y*sin(theta_y1),3);
	//equal to Integrate[exp(-I*p*sin(theta)*y*sin(theta_y1)*cos(theta_y2))*pow(sin(theta_y2),3),{theta_y2,0,Pi}]
	//theta_y2 is the angle between y and rp-plane
	//cout<<"res2: "<<res2<<endl;
	res*=res2;
	
	double res3=cos(p*cos(theta)*y*cos(theta_y1));//real part of exp(-I*p*cos(theta)*y*cos(theta_y1))
	//cout<<"res3: "<<res3<<endl;
	res*=res3;

	res*=pow(y,5)*pow(sin(theta_y1),4)*2*Pi*Pi;//volume element of hyper-spherical coordinate (y,theta_y1,theta_y2,theta_y3,theta_y4,phi_y);

	res/=pow(2.*Pi,6);//fourie transform normalization coefficient

	//res*=Pi*Pi*Pi*pow(r,5);//volume element of vec_r;
	//res*=8./3.*Pi*Pi*pow(p,5)*pow(sin(theta),4);//volume element of vec_p;

	return res;
}

/*double wignerdst_Integrand(double *x,double *par)
{
	double y=x[0];
	double theta_y1=x[1];//the angle between y and r

	double r=par[0];
	double p=par[1];
	double theta=x[2];//the angle between r and p

	double rp=sqrt(r*r+y*y/4.-r*y*cos(theta_y1));//|vec_r+vec_y/2|
	double rn=sqrt(r*r+y*y/4.+r*y*cos(theta_y1));//|vec_r-vec_y/2|
	
	double norm=0;
	double res=0;
	for(int iNc1=1;iNc1<=angnum;iNc1++)
	{
		int q1=anglemomentlist->GetBinContent(iNc1,1);
		int lx1=anglemomentlist->GetBinContent(iNc1,2);
		int ly1=anglemomentlist->GetBinContent(iNc1,3);
		int L1=anglemomentlist->GetBinContent(iNc1,4);
		int sjk1=anglemomentlist->GetBinContent(iNc1,5);
		int two_Sa1=anglemomentlist->GetBinContent(iNc1,6);
		int tjk1=anglemomentlist->GetBinContent(iNc1,7);

		for(int iNc2=1;iNc2<=angnum;iNc2++)
		{
			int q2=anglemomentlist->GetBinContent(iNc2,1);
			int lx2=anglemomentlist->GetBinContent(iNc2,2);
			int ly2=anglemomentlist->GetBinContent(iNc2,3);
			int L2=anglemomentlist->GetBinContent(iNc2,4);
			int sjk2=anglemomentlist->GetBinContent(iNc2,5);
			int two_Sa2=anglemomentlist->GetBinContent(iNc2,6);
			int tjk2=anglemomentlist->GetBinContent(iNc2,7);

			for(int i1=0;n1<=Nselect;i1++)
			{
				int n1=(int)state_select[i1][0];
				double par1=state_select[i1][1];

				for(int n2=nstart;n2<=nstart+nmax;n2++)
				{
					double par2=function_parameter->GetBinContent(iNc2+(n2-nstart)*angnum);
					if(sjk1!=sjk2 || two_Sa1!=two_Sa2 || tjk1!=tjk2 || abs(par1)<=0.1 || abs(par2)<=0.1)continue;
					double res1=par1*par2*orthbasis_radial(rp,n1)/pow(rp,5/2)*orthbasis_radial(rn,n2)/pow(rn,5/2)*16/Pi*1/(4*Pi)*1/(4*Pi);
					//u(rp)*u(rp)*Y_00(theta_rp_J1,phi_rp_J1)*Y_00(theta_rp_J2,phi_rp_J2)*P_0,0+1/2,0+1/2(alpha_rp)*conjugate_Y_00(theta_rn_J1,phi_rn_J1)*conjugate_Y_00(theta_rn_J2,phi_rn_J2)*P_0,0+1/2,0+1/2(alpha_rn)
					//Ji represent the jacobi-transformed coordinate

					//cout<<"res1: "<<res1<<endl;
					res+=res1;

					if(n1==n2 && iNc1==iNc2)norm+=par1*par1;
				}
			}
		}
	}
	//psi(vec_r+vec_y/2)*conjugate_psi(vec_r-vec_y/2)
	//cout<<norm<<endl;
	res/=norm;

	//exp(-I*vec_p*vec_y)=exp(-I*p*sin(theta)*y*sin(theta_y1)*cos(theta_y2)-I*p*cos(theta)*y*cos(theta_y1))
	double res2=4*(sin(p*sin(theta)*y*sin(theta_y1))-p*sin(theta)*y*sin(theta_y1)*cos(p*sin(theta)*y*cos(theta_y1)))/pow(p*sin(theta)*y*sin(theta_y1),3);
	//equal to Integreat[exp(-I*p*sin(theta)*y*sin(theta_y1)*cos(theta_y2))*pow(sin(theta_y2),3),{theta_y2,0,Pi}]
	//theta_y2 is the angle between y and rp-plane
	//cout<<"res2: "<<res2<<endl;
	res*=res2;
	
	double res3=cos(p*cos(theta)*y*sin(theta_y1));//real part of exp(-I*p*cos(theta)*y*cos(theta_y1))
	//cout<<"res3: "<<res3<<endl;
	res*=res3;

	res*=pow(y,5)*pow(sin(theta_y1),4)*2*Pi*Pi;//volume element of hyper-spherical coordinate (y,theta_y1,theta_y2,theta_y3,theta_y4,phi_y);

	res/=pow(2*Pi,6);//fourie transform normalization coefficient

	res*=Pi*Pi*Pi*pow(r,5);//volume element of vec_r;
	res*=8/3*Pi*Pi*pow(p,5)*pow(sin(theta),4);//volume element of vec_p;

	return res;
}*/
