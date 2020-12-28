#include "head.h"
//#include "../../include/H-3.h"
//#include "../../include/He-3.h"
#include "../../include/NNOmega.h"

using namespace std;

double phiki[3][3];
int assignmentphiki();

double Cabc(int a, int b, int c);
double fabc(int a, int b, int c);
double RR(int i, int k, int qi, int lxi, int lyi, int qk, int lxk, int lyk, int L);

int main(int argc, char *argv[])
{
  assignmentphiki();

	int L=atoi(argv[1]);
	const int qmax=12;
	const int lmax=6;

	char pathroot[100];
	sprintf(pathroot,"../RR_Matrix/RR@{qmax=%d,lmax=%d,L=%d,%s}.root",qmax,lmax,L,nuclear);
	TFile RRsaveroot(pathroot,"create");
	if(RRsaveroot.IsZombie())
	{
		cout<<pathroot<<" has been created"<<endl;
		return 0;
	}

  const int dim=qmax+lmax+1;
	//double RRM[3][3][dim+1][dim*2+1][dim*2+1][dim+1][dim*2+1];

	const int bins[]={3,3,dim,dim,dim,dim,dim};
	const double xmin[]={1,1,0,0,0,0,0};
	const double xmax[]={3,3,dim,dim,dim,dim,dim};
	THnD *RRMn=new THnD("RRcoeff","RRMatrix",7,bins,xmin,xmax);

	for(int i=1;i<=3;i++)
	{
		for(int j=2;j<=3;j++)
		{
			for(int qi=0;qi<=qmax;qi++)
			{
				for(int lxi=0;lxi<=lmax;lxi++)
				{
					for(int lyi=abs(lxi-L);lyi<=min(lxi+L,lmax);lyi++)
					{
						for(int qj=0;qj<=(qi*2+lxi+lyi)/2;qj++)
						{
							for(int lxj=0;lxj<=qi*2+lxi+lyi-qj*2;lxj++)
							{
								int lyj=qi*2+lxi+lyi-qj*2-lxj;
								if(lyj<0||lyj+lxj<L||abs(lyj-lxj)>L)continue;

								int pi=lyi+dim*(lxi+dim*(i));
								int pj=lyj+dim*(lxj+dim*(j));
								if(pi>pj)continue;
								double RR_point=RR(i,j,qi,lxi,lyi,qj,lxj,lyj,L);
								int idx[]={i,j,qi,lxi,lyi,qj,lxj};
								RRMn->SetBinContent(idx,RR_point);
								int jdx[]={j,i,qj,lxj,lyj,qi,lxi};
								RRMn->SetBinContent(jdx,RR_point);

								//RRM[i][j][qi][lxi][lyi][qj][lxj]=RR_point;
								//RRM[j][i][qj][lxj][lyj][qi][lyi]=RR_point;
							}
						}
					}
				}
			}
		}
	}

	RRMn->Write();
	RRsaveroot.Close();
	delete RRMn;

	/*char path[100];
	sprintf(path,"../RR_Matrix/RR@{qmax=%d,lmax=%d,L=%d,%s}.txt",qmax,lmax,L,nuclear);
	ofstream RRsave(path);

	for(int i=0;i<3;i++)
	{
		for(int j=0;j<3;j++)
		{
			for(int qi=0;qi<=qmax+lmax;qi++)
			{
				for(int lxi=0;lxi<=lmax+lmax+qmax*2;lxi++)
				{
					for(int lyi=0;lyi<=lmax*2+qmax*2;lyi++)
					{
						for(int qj=0;qj<=(qi*2+lxi+lyi)/2;qj++)
						{
							for(int lxj=0;lxj<=qi*2+lxi+lyi-qj*2;lxj++)
							{
								RRsave<<RRM[i][j][qi][lxi][lyi][qj][lxj]<<endl;
							}
						}
					}
				}
			}
		}
	}
	RRsave.close();*/

  return 1;
}

double Cabc(int a, int b, int c)
{
	double den=TMath::Gamma(2*a+b+c+2);
	double num=TMath::Gamma(a+b+1.5)*TMath::Gamma(a+c+1.5)*TMath::Gamma(a+1)*TMath::Gamma(a+b+c+2);
	return den/num;
}

double fabc(int a,int b,int c)
{
	return sqrt((2*a+1)*(2*b+1))*CGcoeff(2*a,2*b,2*c,0,0,0);
}

double RR(int i, int k, int qi, int lxi, int lyi, int qk, int lxk, int lyk, int L)
{
	if(i==k)
	{
		if(qi==qk && lxi==lxk &&lyi==lyk)
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}
	i-=1;
	k-=1;
	if(qi*2+lxi+lyi!=2*qk+lxk+lyk)return 0;
	double sumRe=0;
	double sumIm=0;
	for(int lambda1=0;lambda1<=qi*2+lxi+lyi;lambda1++)
	{
		for(int lambda2=0;lambda2<=qi*2+lxi+lyi-lambda1;lambda2++)
		{
			for(int lambda3=0;lambda3<=qi*2+lxi+lyi-lambda1-lambda2;lambda3++)
			{
				for(int lambda4=0;lambda4<=qi*2+lxi+lyi-lambda1-lambda2-lambda3;lambda4++)
				{
					double suminside=0;
					for(int mu=0;mu<=(qi*2+lxi+lyi-lambda1-lambda2-lambda3-lambda4)/2;mu++)
					{
						int nv=(qi*2+lxi+lyi-lambda1-lambda2-lambda3-lambda4)/2-mu;
						suminside+=pow(-1,mu)*Cabc(mu,lambda3,lambda4)*Cabc(nv,lambda1,lambda2)
						             *pow(cos(phiki[i][k]),2*nv+lambda1+lambda2)*pow(-sin(phiki[i][k]),2*mu+lambda3+lambda4);
					}
					double Im=lambda3+lambda4+lyi-lyk;
					if(Im/2-(int)(Im/2)==0)
					{
						sumRe+=pow(-1,lambda1+lambda2+(int)Im/2)*fabc(lambda1,lambda3,lxk)*fabc(lambda4,lambda2,lyk)
					        *fabc(lambda1,lambda4,lxi)*fabc(lambda3,lambda2,lyi)*wigner_9j(2*lambda1,2*lambda3,2*lxk,2*lambda4,2*lambda2,2*lyk,2*lxi,2*lyi,2*L)*suminside;
					}
					else
					{
						sumIm+=pow(-1,lambda1+lambda2+(Im-1)/2)*fabc(lambda1,lambda3,lxk)*fabc(lambda4,lambda2,lyk)
					        *fabc(lambda1,lambda4,lxi)*fabc(lambda3,lambda2,lyi)*wigner_9j(2*lambda1,2*lambda3,2*lxk,2*lambda4,2*lambda2,2*lyk,2*lxi,2*lyi,2*L)*suminside;
					}
				}
			}
		}
	}
	sumRe*=Pi/4/sqrt(Cabc(qi,lxi,lyi)*Cabc(qk,lxk,lyk));
	sumIm*=Pi/4/sqrt(Cabc(qi,lxi,lyi)*Cabc(qk,lxk,lyk));
	if(sumIm>1.e-6)
	{
		cout<<"Im="<<sumIm<<", "<<(sumIm>1.e-6)<<endl;
	}
	return sumRe;
}

int assignmentphiki()
{
	phiki[0][0]=0;
	phiki[1][1]=0;
	phiki[2][2]=0;

	phiki[0][1]=atan(-sqrt(M*m[2]/(m[1]*m[0])));
	phiki[0][2]=atan(sqrt(M*m[1]/(m[2]*m[0])));
	phiki[1][0]=atan(sqrt(M*m[2]/(m[0]*m[1])));
	phiki[1][2]=atan(-sqrt(M*m[0]/(m[1]*m[2])));
	phiki[2][0]=atan(-sqrt(M*m[1]/(m[0]*m[2])));
	phiki[2][1]=atan(sqrt(M*m[0]/(m[1]*m[2])));

	return 0;
}
