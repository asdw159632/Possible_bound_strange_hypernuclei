void transform()
{
	char path[100];
	char num[10];
	double logsrmax=2;
	double logsrmin=-8;
	int Nr=50;
	double logspmax=1;
	double logspmin=-8;
	int Np=50;

	TH2D wigdst_test("wigtest","wigtest",Nr,logsrmin,logsrmax,Np,logspmin,logspmax);
	double dp=(logspmax-logspmin)/Np;
	for(int i=1; i<=Nr; i++)
	{
		double ri=logsrmin+(logsrmax-logsrmin)/Nr*i;
		sprintf(num,"%.8f",ri);
		if(ri>=1.)
		{
			sprintf(path,"%.4s.txt",num);
		}else if(ri<1. && ri>0)
		{
			sprintf(path,"%.5s.txt",num);
		}else if(ri==0.)
		{
			sprintf(path,"%.1s.txt",num);
		}else if(ri<0. && ri>-1.)
		{
			sprintf(path,"%.6s.txt",num);
		}else if(ri<=-1.)
		{
			sprintf(path,"%.5s.txt",num);
		}

		ifstream f;
		f.open(path);
		if(!f)
		{
			cout<<path<<" open error"<<endl;
			return;
		}
		double rin;
		double pin;
		double wdstin;
		int count=0;
		while(f>>rin)
		{
			count++;
			f>>pin>>wdstin;
			int pbin=floor((pin+0.001-logspmin)/dp);
			//if(i==1)cout<<count<<" "<<pbin<<" "<<wigdst_test.GetYaxis()->FindBin(pin-0.001)<<endl;
			wigdst_test.SetBinContent(i,pbin,wdstin);
			if(wdstin<0.3&&i<40&&pbin<40)cout<<i<<" "<<pbin<<endl;
		}
		f.close();
	}
	TFile *save = new TFile("pnOmega_qmax3_lmax2_Lmax0_angnum8_nmax5_nstart2_2600-test-scaled.root","recreate");
	wigdst_test.Write();
	save->Close();
	delete save;
}
