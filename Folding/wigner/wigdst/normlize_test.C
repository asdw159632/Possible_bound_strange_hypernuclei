void normlize_test()
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
	TH2D *wigdst=(TH2D*)read.Get("wigdst");
	
	double norm=0;
	int Nrmin=wigdst->GetXaxis()->GetFirst();
	int Nrmax=wigdst->GetXaxis()->GetLast();
	double dr=wigdst->GetXaxis()->GetBinWidth(1);
	int Npmin=wigdst->GetYaxis()->GetFirst();
	int Npmax=wigdst->GetYaxis()->GetLast();
	double dp=wigdst->GetYaxis()->GetBinWidth(1);
	for(int nr=Nrmin;nr<=Nrmax;nr++)
	{
		double r=wigdst->GetXaxis()->GetBinCenter(nr);
		r=exp(r)*aBorhr;
		for(int np=Npmin;np<=Npmax;np++)
		{
			double p=wigdst->GetYaxis()->GetBinCenter(np);
			p=exp(p)*ku;
			p/=hbar;

			double wigner=wigdst->GetBinContent(nr,np);
			norm+=wigner*8.*pow(Pi,2)*pow(r,3)*pow(p,3)*dr*dp/(pow(2*Pi,3));
		}
	}
	cout<<norm<<endl;
}
