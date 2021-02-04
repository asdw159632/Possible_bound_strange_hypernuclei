void Dnormlize_test()
{
	double hbar = 197.3269788;// MeV fm/c
	double aBorhr;
	double ku;
	aBorhr = 1;// fm for length unit
	ku = hbar/aBorhr;// MeV/c  for momentum unit
	double Pi=TMath::Pi();

	cout<<"test which file: ";
	char file[500];
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
	int Ntmin=wigdst->GetZaxis()->GetFirst();
	int Ntmax=wigdst->GetZaxis()->GetLast();
	double dt=wigdst->GetZaxis()->GetBinWidth(1);
	for(int nt=Ntmin;nt<=Ntmax;nt++)
	{
		double th=wigdst->GetZaxis()->GetBinCenter(nt);
		double one=0;
		for(int nr=Nrmin;nr<=Nrmax;nr++)
		{
			double r=wigdst->GetXaxis()->GetBinCenter(nr);
			if(r>2)continue;
			r=exp(r)*aBorhr;
			for(int np=Npmin;np<=Npmax;np++)
			{
				double p=wigdst->GetYaxis()->GetBinCenter(np);
				if(p>1)continue;
				p=exp(p)*ku;
				p/=hbar;
				double wigner=wigdst->GetBinContent(nr,np,nt);
				one+=wigner*8./3.*pow(Pi,5)*pow(r,6)*pow(p,6)*dr*dp/*/(pow(2*Pi,3))*/;
				norm+=wigner*8./3.*pow(Pi,5)*pow(r,6)*pow(p,6)*dr*dp*pow(sin(th),4)*dt/*/(pow(2*Pi,3))*/;
			}
		}
		//cout<<one<<endl;
	}
	cout<<norm<<endl;
}
