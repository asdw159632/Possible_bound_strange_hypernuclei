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
	int Npmin=wigdst->GetYaxis()->GetFirst();
	int Npmax=wigdst->GetYaxis()->GetLast();
	int Ntmin=wigdst->GetZaxis()->GetFirst();
	int Ntmax=wigdst->GetZaxis()->GetLast();
	for(int nt=Ntmin;nt<=Ntmax;nt++)
	{
		double dt=wigdst->GetZaxis()->GetBinWidth(nt);
		double th=wigdst->GetZaxis()->GetBinCenter(nt);
		th-=dt;
		for(int nr=Nrmin;nr<=Nrmax;nr++)
		{
			double dr=wigdst->GetXaxis()->GetBinWidth(nr);
			double r=wigdst->GetXaxis()->GetBinCenter(nr);
			r=exp(r)*aBorhr;
			for(int np=Npmin;np<=Npmax;np++)
			{
				double dp=wigdst->GetYaxis()->GetBinWidth(np);
				double p=wigdst->GetYaxis()->GetBinCenter(np);
				p=exp(p)*ku;
				p/=hbar;
				double wigner=wigdst->GetBinContent(nr,np,nt);
				norm+=wigner*4*Pi*2*Pi*pow(r,3)*pow(p,3)*dr*dp*sin(th)*dt/(pow(2*Pi,3));
			}
		}
	}
	cout<<norm<<endl;
}
