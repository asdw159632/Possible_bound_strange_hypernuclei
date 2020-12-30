void normlize_test()
{
	double fmmev=5.06/1000;
	double aBohr=1.;//fm
	double ku=1/aBohr/fmmev;//MeV
	double Pi=TMath::Pi();

	cout<<"test which file: ";
	char file[100];
	cin>>file;
	TFile read(file,"read");
	TH3D *wigdst=(TH3D*)read.Get("wigdst");
	
	double norm=0;
	int Nrmin=wigdst->GetXaxis()->GetFirst();
	int Nrmax=wigdst->GetXaxis()->GetLast();
	double dr=wigdst->GetXaxis()->GetBinWidth(1);
	int Npmin=wigdst->GetYaxis()->GetFirst();
	int Npmax=wigdst->GetYaxis()->GetLast();
	double dp=wigdst->GetYaxis()->GetBinWidth(1);
	int Nthetamin=wigdst->GetZaxis()->GetFirst();
	int Nthetamax=wigdst->GetZaxis()->GetLast();
	double dtheta=wigdst->GetZaxis()->GetBinWidth(1);
	for(int nr=Nrmin;nr<=Nrmax;nr++)
	{
		for(int np=Npmin;np<=Npmax;np++)
		{
			for(int ntheta=Nthetamin;ntheta<=Nthetamax;ntheta++)
			{
				double wigner=wigdst->GetBinContent(nr,np,ntheta);
				double r=wigdst->GetXaxis()->GetBinUpEdge(nr);
				r=exp(r)*aBohr;
				double p=wigdst->GetYaxis()->GetBinUpEdge(np);
				p=exp(p)*ku*fmmev;
				double theta=wigdst->GetZaxis()->GetBinUpEdge(ntheta);
				norm+=wigner*pow(Pi,5)*pow(r,6)*8./3.*pow(p,6)*pow(sin(theta),4)*dr*dp*dtheta;
			}
		}
	}
	cout<<norm<<endl;
}
