#include <head.h>

using namespace std;

#include "../../include/H-3.h"
//#include "../../include/He-3.h"
//#include "../../include/NNOmega.h"

TH3D *rho_density;

#ifdef E200GeV
double R0 = 12;//fm
double tau0 = 9;//fm/c
double delta_tau = 3.5;//fm/c  */
#elif defined E2_76TeV
double R0 = 19.7;//fm
double tau0 = 15.5;//fm/c
double delta_tau = 1;//fm/c  */
#endif

double hbar = 197.3269788;// MeV fm/c
double aBorhr=1;//fm
double ku=1/aBorhr/fmmev;//MeV


double nDst(TLorentzVector vmP, double vTkin, double vrho_0);
Double_t dN_momentum_BLW(double *x, double *par);
Double_t dN_coordinates(double *x, double *par);
Double_t dN_tau(double *x, double *par);

void coal_p1p2p3(TLorentzVector *vmP_p1, TLorentzVector *vmR_p1, TLorentzVector *vmP_p2, TLorentzVector *vmR_p2, TLorentzVector *vmP_p3, TLorentzVector *vmR_p3, int vnum_p1, int vnum_p2, int vnum_p3);
double rho_wigner(TLorentzVector *vcompP, TLorentzVector *vcompR);

TH1D *cluster_pT_Dst = new TH1D("cluster_pT_Dst","cluster_pT_Dst",100,0,5);

Long_t nseed_cmd();
void executeCMD(const char *cmd, char *result);

int main(int argc, char **argv)
{
	cout<<nuclear<<endl;

  Long_t nseedc = nseed_cmd();
  cout<<"global nseedc="<<nseedc<<endl;
  //gRandom->SetSeed(time(0));
  gRandom->SetSeed(nseedc);

  TList *myList = new TList();
  myList->SetName("myList");

  char *workDir = argv[1];
  char *jobID = argv[2];

  char inRhoDir[512];
	char *target=argv[3];
  sprintf(inRhoDir,"%s/wigdst/%s.root",workDir,target);
  TFile fin(inRhoDir);
  rho_density = (TH3D *) fin.Get("wigdst");

  //init par

  //aBorhr = sqrt(hbar*hbar/(Q/fmmev*V0Scale));// fm for length unit
  //ku = hbar/aBorhr;// MeV/c  for momentum unit
  //init par end

  TH1D *p1_pT_Dst = new TH1D(ti_p1_pT_Dst,ti_p1_pT_Dst,50,0,5);
  myList->Add(p1_pT_Dst);

  TH1D *p2_pT_Dst = new TH1D(ti_p2_pT_Dst,ti_p2_pT_Dst,50,0,5);
  myList->Add(p2_pT_Dst);

#ifndef ABB
  TH1D *p3_pT_Dst = new TH1D(ti_p3_pT_Dst,ti_p3_pT_Dst,50,0,5);
  myList->Add(p3_pT_Dst);
#endif

  myList->Add(cluster_pT_Dst);
  
  double dpT = p1_pT_Dst->GetBinWidth(1);

  const int numEvt=100000;

  TLorentzVector mP_p1[num_p1];
  TLorentzVector mR_p1[num_p1];
  TLorentzVector mP_p2[num_p2];
  TLorentzVector mR_p2[num_p2];
#ifndef ABB
  TLorentzVector mP_p3[num_p3];
  TLorentzVector mR_p3[num_p3];
#endif

  double pT;
  double phi_p;
  double eta_p;

  double tau;
  double r_rho;
  double eta_s;
  double phi_s;

  TF1 *func_pT_p1 = new TF1(ti_func_pT_p1,dN_momentum_BLW,0.,5.,3);
  func_pT_p1->SetParameters(m[0]/fmmev/1000/*GeV*/, Tkin_p1,rho_0_p1);

  TF2 *func_coordinates_p1 = new TF2(ti_func_coordinates_p1,dN_coordinates,0,R0, -0.5,0.5, 7);//

  cout<<"func_p1 created!!!"<<endl;

  TF1 *func_pT_p2 = new TF1(ti_func_pT_p2,dN_momentum_BLW,0.,5,3);
  func_pT_p2->SetParameters(m[1]/fmmev/1000/*GeV*/, Tkin_p2,rho_0_p2);

  TF2 *func_coordinates_p2 = new TF2(ti_func_coordinates_p2,dN_coordinates,0,R0, -0.5,0.5, 7);//

  cout<<"func_p2 created!!!"<<endl;

#ifndef ABB
  TF1 *func_pT_p3 = new TF1(ti_func_pT_p3,dN_momentum_BLW,0.,5,3);
  func_pT_p3->SetParameters(m[2]/fmmev/1000/*GeV*/, Tkin_p3,rho_0_p3);

  TF2 *func_coordinates_p3 = new TF2(ti_func_coordinates_p3,dN_coordinates,0,R0, -0.5,0.5, 7);//

  cout<<"func_p3 created!!!"<<endl;
#endif

  TF1 *func_tau = new TF1("func_tau",dN_tau,tau0-0.5*delta_tau, tau0+0.5*delta_tau, 2);
  func_tau->SetParameters(tau0, delta_tau);

  TH1D *x_dst = new TH1D("x_dst","x_dst",30,-5,5);
  TH1D *y_dst = new TH1D("y_dst","y_dst",30,-5,5);
  TH1D *z_dst = new TH1D("z_dst","z_dst",30,-5,5);

  TH1D *r_dst = new TH1D("r_dst","r_dst",40,0,20);
  TH1D *rho_dst = new TH1D("rho_dst","r_dst",40,0,20);

  myList->Add(x_dst);
  myList->Add(y_dst);
  myList->Add(z_dst);
  myList->Add(r_dst);
  myList->Add(rho_dst);

  TH1F *h1_event = new TH1F("h1_event","h1_event",1,0,1);
  myList->Add(h1_event);

  double rMag, mT;

  for(int ievt=0;ievt<numEvt; ievt++)
  {
    h1_event->Fill(0.);

    //p1
    for(int ipart=0; ipart<num_p1; ipart++)
		{
			pT = func_pT_p1->GetRandom();

			p1_pT_Dst->Fill(pT, 1./(2. * Pi * pT * dpT));

			phi_p = gRandom->Rndm() * 2. * Pi;
			eta_p = gRandom->Rndm() - 0.5;

			phi_s = gRandom->Rndm() * 2. * Pi;
			func_coordinates_p1->SetParameters(m[0]/fmmev/1000/*GeV*/, Tkin_p1, rho_0_p1,  pT, eta_p, phi_p, phi_s);

			func_coordinates_p1->GetRandom2(r_rho, eta_s);
			tau = func_tau->GetRandom();

			x_dst->Fill(r_rho*cos(phi_s));
			y_dst->Fill(r_rho*sin(phi_s));
			z_dst->Fill(tau*sinh(eta_s));

			rMag = sqrt(r_rho*r_rho + pow(tau*sinh(eta_s),2));
			if(rMag>1e-7) r_dst->Fill(rMag, 1/rMag/rMag);

			rho_dst->Fill(r_rho, 1./r_rho);

			mT = sqrt(pT*pT+m[0]/fmmev/1000/*GeV*/*m[0]/fmmev/1000/*GeV*/);
			mP_p1[ipart].SetPxPyPzE(pT*cos(phi_p), pT*sin(phi_p), mT*sinh(eta_p), mT*cosh(eta_p));
			mR_p1[ipart].SetXYZT(r_rho*cos(phi_s), r_rho*sin(phi_s), tau*sinh(eta_s), tau*cosh(eta_s));
		}

		//p2
		for(int ipart=0; ipart<num_p2; ipart++)
		{
			pT = func_pT_p2->GetRandom();

			p2_pT_Dst->Fill(pT, 1./(2. * Pi * pT * dpT));

			phi_p = gRandom->Rndm() * 2. * Pi;
			eta_p = gRandom->Rndm() - 0.5;

			phi_s = gRandom->Rndm() * 2. * Pi;
			func_coordinates_p2->SetParameters(m[1]/fmmev/1000/*GeV*/, Tkin_p2, rho_0_p2,  pT, eta_p, phi_p, phi_s);

			func_coordinates_p2->GetRandom2(r_rho, eta_s);
			tau = func_tau->GetRandom();

			mT = sqrt(pT*pT+m[1]/fmmev/1000/*GeV*/*m[1]/fmmev/1000/*GeV*/);
			mP_p2[ipart].SetPxPyPzE(pT*cos(phi_p), pT*sin(phi_p), mT*sinh(eta_p), mT*cosh(eta_p));
			mR_p2[ipart].SetXYZT(r_rho*cos(phi_s), r_rho*sin(phi_s), tau*sinh(eta_s), tau*cosh(eta_s));
		}

#ifndef ABB
		//p3
		for(int ipart=0; ipart<num_p3; ipart++)
		{
			pT = func_pT_p3->GetRandom();

			p3_pT_Dst->Fill(pT, 1./(2. * Pi * pT * dpT));

			phi_p = gRandom->Rndm() * 2. * Pi;
			eta_p = gRandom->Rndm() - 0.5;

			phi_s = gRandom->Rndm() * 2. * Pi;
			func_coordinates_p3->SetParameters(m[2]/fmmev/1000/*GeV*/, Tkin_p3, rho_0_p3,  pT, eta_p, phi_p, phi_s);

			func_coordinates_p3->GetRandom2(r_rho, eta_s);
			tau = func_tau->GetRandom();

			mT = sqrt(pT*pT+m[2]/fmmev/1000/*GeV*/*m[2]/fmmev/1000/*GeV*/);
			mP_p3[ipart].SetPxPyPzE(pT*cos(phi_p), pT*sin(phi_p), mT*sinh(eta_p), mT*cosh(eta_p));
			mR_p3[ipart].SetXYZT(r_rho*cos(phi_s), r_rho*sin(phi_s), tau*sinh(eta_s), tau*cosh(eta_s));
		}
#endif

#ifdef ABB
		coal_p1p2p3(mP_p1, mR_p1, mP_p2, mR_p2, mP_p2, mR_p2, num_p1, num_p2, num_p2);
#else
		coal_p1p2p3(mP_p1, mR_p1, mP_p2, mR_p2, mP_p3, mR_p3, num_p1, num_p2, num_p3);
#endif

		cout<<ievt+1<<" events finished"<<endl;
	}

	char outName[512];
	sprintf(outName,"%s/results/result-%s.root",workDir,jobID);

	TFile *fout = new TFile(outName,"RECREATE");

	myList->Write();
	fout->Write();

	delete fout;
}

Double_t dN_momentum_BLW(double *x, double *par)
{
	double vpT = x[0];

	double mass = par[0];
	double Tkin = par[1];
	double rho_0 = par[2];

	//double mass_p = 1.67245;

	TLorentzVector mP;

	//double pT, phi, eta;
	double phi, eta;

	int Nphi_p=24;//24
	double dphi = 2.*Pi/Nphi_p;

	//set eta=0;
	eta=0;

	double sumN=0;
	for(int iphi = 0; iphi<Nphi_p; iphi++)
	{
		phi = dphi * iphi;

		mP.SetPtEtaPhiM(vpT, eta, phi, mass);
		sumN += nDst(mP, Tkin, rho_0) * dphi;
	}

	return sumN * vpT;
}

double nDst(TLorentzVector vmP, double vTkin, double vrho_0)
{
	double Tkin = vTkin;
	double rho_0 = vrho_0;

	double mT = vmP.Mt();
	double pT = vmP.Pt();
	double phi_p = vmP.Phi();
	double eta_p = vmP.Rapidity();

	//int it=0;//tau
	int ies=0;//eta_s
	int ir=0;//r
	int iphis=0;//phi_s

	//int Ntau = 10;//50
	//double dtau = 6.*delta_tau / Ntau;//from -5sigma to 5sigma, here delta_tau is sigma of proper time distribution

	int Neta_s = 20;//50
	double deta_s = 10.0 / Neta_s;//from -5 to 5

	int Nr = 50;//100
	double dr = R0 / Nr;

	int Nphi_s=24;//24
	double dphi_s = 2.*Pi / Nphi_s;

	double eta_s, r, phi_s;

	double rho;
	double pmu_umu;
	//  double J_tau;

	double sum = 0;
	//for(it=0; it<Ntau; it++)
	//{
	//tau = tau0 - 3.*delta_tau + dtau * it;

	for(ies=0; ies<Neta_s; ies++)
	{
		eta_s = -5. + deta_s * ies;

		for(ir=0; ir<Nr; ir++)
		{
			r = 0. + dr *ir;

			for(iphis=0; iphis<Nphi_s; iphis++)
			{
				phi_s = 0. + dphi_s * iphis;

				rho = rho_0 * r / R0;
				pmu_umu = mT * cosh(rho) *cosh(eta_s - eta_p) - pT * sinh(rho) * cos(phi_p - phi_s);
				//	    J_tau = 1. / (delta_tau * sqrt(2. * Pi)) * exp( - pow( (tau - tau0) / delta_tau, 2) / 2.);
				//sum += mT * cosh(eta_s - eta_p) / (exp(pmu_umu/Tkin) + 1.) * J_tau * tau * dtau * deta_s * r * dr * dphi_s;
				sum += mT * cosh(eta_s - eta_p) / (exp(pmu_umu/Tkin) + 1.) * deta_s * r * dr * dphi_s;
			}
		}
	}
	//}

	return sum;
}

Double_t dN_coordinates(double *x, double *par)
{
	double r = x[0];
	double eta_s = x[1];

	double mass  = par[0];
	double Tkin  = par[1];
	double rho_0 = par[2];
	double pT    = par[3];
	double eta_p = par[4];
	double phi_p = par[5];
	double phi_s = par[6];

	double mT = sqrt(pT*pT + mass * mass);

	double rho = rho_0 * r / R0;
	double pmu_umu = mT * cosh(rho) *cosh(eta_s - eta_p) - pT * sinh(rho) * cos(phi_p - phi_s);

	double theVal = mT * cosh(eta_s - eta_p) / (exp(pmu_umu/Tkin) + 1.) * r;

	return theVal;
}

Double_t dN_tau(double *x, double *par)
{
	double tau = x[0];

	double vtau0 = par[0];
	double vdelta_tau = par[1];

	double Pi = TMath::Pi();

	double J_tau = 1. / (vdelta_tau * sqrt(2. * Pi)) * exp( - pow( (tau - vtau0) / vdelta_tau, 2) / 2.);
	double theVal = J_tau * tau;

	return theVal;
}


void coal_p1p2p3(TLorentzVector *vmP_p1, TLorentzVector *vmR_p1, TLorentzVector *vmP_p2, TLorentzVector *vmR_p2, TLorentzVector *vmP_p3, TLorentzVector *vmR_p3, int vnum_p1, int vnum_p2, int vnum_p3)
{
  double rho_W=0;

  TLorentzVector mP[3];
  TLorentzVector mR[3];
  TLorentzVector mP_p1p2p3;
  
  double dpT = cluster_pT_Dst->GetXaxis()->GetBinWidth(1);
  for(int i=0; i<vnum_p1; i++)
  {
    mP[0] = vmP_p1[i];
    mR[0] = vmR_p1[i];
    for(int j=0; j<vnum_p2; j++)
		{
			mP[1] = vmP_p2[j];
			mR[1] = vmR_p2[j];
			for(int l=j+1; l<vnum_p3; l++)
			{
				mP[2] = vmP_p3[l];
				mR[2] = vmR_p3[l];

				mP_p1p2p3.SetVectM((mP[0] + mP[1] + mP[2]).Vect(), M/fmmev/1000/*GeV*/);
				//if(mP_p1p2p3.Pt()<1.e-7) continue;
				if(fabs(mP_p1p2p3.Rapidity())>0.5) continue;

				rho_W = rho_wigner(mP, mR);
				//if(rho_W<1e-20) continue;
				//if(rho_W>1e-5)cout<<"Fill: "<<mP_p1p2p3.Pt()<<" "<<1./(2.*Pi*mP_p1p2p3.Pt()*dpT) * rho_W * GA<<endl;

				cluster_pT_Dst->Fill(mP_p1p2p3.Pt(), 1./(2.*Pi*mP_p1p2p3.Pt()*dpT) * rho_W * GA);
			}
		}
	}
}

double rho_wigner(TLorentzVector *vcompP, TLorentzVector *vcompR)
{
  double ETotalLab; //total energy of cluster
  TVector3 Pcmtolab; //3d-momenta of cluster
  TVector3 betalabtocm; //velocity of lab relative to center of cluster
  //TLorentzRotation LRltoc;

  TLorentzVector pcpt1Rcm; //coordinates of particle in center of mass frame
  TLorentzVector pcpt2Rcm;
  TLorentzVector pcpt3Rcm;

  double pcpt1Timecm;//time of particle in center of mass frame
  double pcpt2Timecm;
  double pcpt3Timecm;

  TVector3 pcpt1R3cm;//3d-coordinates of particle in center of mass frame
  TVector3 pcpt2R3cm;
  TVector3 pcpt3R3cm;

  TLorentzVector pcpt1Pcm;//momenta of particle in center of mass frame
  TLorentzVector pcpt2Pcm;
  TLorentzVector pcpt3Pcm;

  TVector3 pcpt1P3cm;//3d-momenta of particle in center of mass frame
  TVector3 pcpt2P3cm;
  TVector3 pcpt3P3cm;

  TVector3 CMR;//3d-coordinates of center of mass frame
  //TVector3 rho;//3d-relative-coordinates
  TVector3 CMK;//3d-momenta of center of mass frame
  //TVector3 krho;//3d-relative-momenta

  double rhoWH=0.; //possibility

  TLorentzVector current_mP; //momenta of cluster

  ETotalLab = vcompP[0].Energy() + vcompP[1].Energy() + vcompP[2].Energy();
  //      Pcmtolab = pcpt1[i1].particleP.Vect() + pcpt2[i2].particleP.Vect();

  //here skip some cluster

  current_mP = vcompP[0]+vcompP[1]+vcompP[2];
  Pcmtolab = current_mP.Vect();

  betalabtocm = -1./ETotalLab*Pcmtolab;
  TLorentzRotation LRltoc(betalabtocm);

  pcpt1Rcm = LRltoc*vcompR[0];//boost R(4d) from lb to cm frame
  pcpt2Rcm = LRltoc*vcompR[1];
  pcpt3Rcm = LRltoc*vcompR[2];

  pcpt1Timecm = pcpt1Rcm.T();
  pcpt2Timecm = pcpt2Rcm.T();
  pcpt3Timecm = pcpt3Rcm.T();

  //if(fabs(pcpt1Timecm-pcpt2Timecm)>1.0) return rhoWH;

  pcpt1R3cm = pcpt1Rcm.Vect();//extract R(3d) from R(4d)
  pcpt2R3cm = pcpt2Rcm.Vect();
  pcpt3R3cm = pcpt3Rcm.Vect();

  pcpt1Pcm = LRltoc*vcompP[0];//boost P(4d) from lb to cm frame
  pcpt2Pcm = LRltoc*vcompP[1];
  pcpt3Pcm = LRltoc*vcompP[2];

  pcpt1P3cm = pcpt1Pcm.Vect();//extract P(3d) from P(4d)
  pcpt2P3cm = pcpt2Pcm.Vect();
  pcpt3P3cm = pcpt3Pcm.Vect();


  double m1 = vcompP[0].M();
  double m2 = vcompP[1].M();
  double m3 = vcompP[2].M();
	double M_in = m1+m2+m3;
	double Q_in = sqrt(m1*m2*m3/M_in);
	//J is the Jacobian; r = ( r1, r2, r3 )^T; p = ( p1, p2, p3)^T
	//
	//M=m1+m2+m3;
	//Q=sqrt(m1 m2 m3/M); 
	//both of M and Q is defined in *nuclear*.h
	//
	//J = {          m1/M                         m2/M                          m3/M            }
	//    {           0                   -sqrt(m2 m3/((m2+m3)Q))      sqrt(m2 m3/((m2+m3)Q))   }
	//    { -sqrt(m1 (m2+m3)/(M Q))      sqrt(m2^2 m1/(M(m2+m3)Q))    sqrt(m3^2 m1/(M(m2+m3)Q)) }
	//
	//(J^-1)^T = {          1                         1                         1            }
	//           {          0              -sqrt(m3 Q/(m2 (m2+m3)))   sqrt(m2 Q/(m3(m2+m3))) }
	//					 { -sqrt((m2+m3)Q/(M m1))    sqrt(m1 Q/(M(m2+m3)))    sqrt(m1 Q/(M(m2+m3)))  }
	//
	//R=J.r; K=(J^-1)^T.p
  CMR = 1./M_in * (m1*pcpt1R3cm + m2*pcpt2R3cm + m3*pcpt3R3cm);
  CMK = pcpt1P3cm + pcpt2P3cm +pcpt3P3cm;

  //convert to unit used in phi-wave function and rho-densoty

	TVector3 rho1=-sqrt( m2*m3/((m2+m3)*Q_in) )*(pcpt2R3cm - pcpt3R3cm);
	TVector3 rho2=-sqrt( (m2+m3)*m1/(M_in*Q_in) )*pcpt1R3cm + sqrt( m1/(M_in*Q_in*(m2+m3)) )*(m2*pcpt2R3cm + m3*pcpt3R3cm);
	double rho_mag=sqrt(rho1.Mag2()+rho2.Mag2());//sqrt(rho1^2+rho2^2); fm
  //cout<<"rho1: "<<rho1.Mag()<<" rho2:"<<rho2.Mag()<<endl;

//  cout<<"rho_mag="<<rho_mag<<endl;

	TVector3 krho1=-sqrt( Q_in/(m2+m3) )*( sqrt(m3/m2)*pcpt2P3cm - sqrt(m2/m3)*pcpt3P3cm);
	TVector3 krho2=-sqrt( (m2+m3)*Q_in/(M_in*m1) )*pcpt1P3cm + sqrt( m1*Q_in/(M_in*(m2+m3)) )*(pcpt2P3cm + pcpt3P3cm);
  double krho_mag = sqrt(krho1.Mag2()+krho2.Mag2());//sqrt(krho1^2+krho2^2); GeV
  //cout<<"krho1: "<<krho1.Mag()<<" krho2:"<<krho2.Mag()<<endl;

	double theta=acos(( rho1.Dot(krho1)+rho2.Dot(krho2) )/(rho_mag*krho_mag));
	//cout<<"theta:"<<theta<<" k:"<<krho_mag<<" rho:"<<rho_mag<<" rhoWH:";

//scale
	rho_mag/=aBorhr;
	rho_mag=log(rho_mag);

	krho_mag=krho_mag*1000/ku;//GeV->MeV, then scale it
	krho_mag=log(krho_mag);

  int rBin = rho_density->GetXaxis()->FindBin(rho_mag);
  int kBin = rho_density->GetYaxis()->FindBin(krho_mag);
	int thetaBin = rho_density->GetZaxis()->FindBin(theta);

  rhoWH = rho_density->GetBinContent(rBin, kBin, thetaBin); 
	if(rhoWH>0)cout<<"("<<rBin<<" "<<kBin<<" "<<thetaBin<<"): "<<rhoWH<<endl;

  if(std::isnan(rhoWH)||std::isinf(rhoWH)) return 0;

  return rhoWH;
}

Long_t nseed_cmd()
{
  //cout<<time(0)+seed()<<endl;
//  char thecmd[64]={"date +%S"};//MacOs
  char thecmd[64]={"date +%N"};//Linux
  char thetns[64]={""};
  executeCMD(thecmd,thetns);

  //    cout<<"thetns is "<<thetns<<endl;

  Long_t thens = Long_t(atof(thetns));//system("date +%N");
  Long_t ts = time(0);

  //    cout<<"thens="<<thens<<endl;
  //    cout<<"the ts="<<ts<<endl;
  //    cout<<"nseed="<<ts+thens<<endl;

  return ts+thens;
  //    sprintf(thetns,"%s","");
}
void executeCMD(const char *cmd, char *result)
{
  char buf_ps[1024];
  char ps[1024]={0};
  FILE *ptr;
  strcpy(ps, cmd);
  if((ptr=popen(ps, "r"))!=NULL)
  {
    while(fgets(buf_ps, 1024, ptr)!=NULL)
    {
      strcat(result, buf_ps);
      if(strlen(result)>1024)
        break;
    }
    pclose(ptr);
    ptr = NULL;
  }
  else
  {
    printf("popen %s error\n", ps);
  }
}
