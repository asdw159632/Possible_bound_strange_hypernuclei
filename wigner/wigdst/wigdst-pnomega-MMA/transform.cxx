void transform()
{
	char path[100];
	char num[10];
	double rmax=10;
	double rmin=0;
	int Nr=50;
	double pmax=1000;
	double pmin=0;
	int Np=50;

	TH2D wigdst_test("wigtest","wigtest",Nr,rmin,rmax,Np,pmin,pmax);
	double dp=(pmax-pmin)/Np;
	for(int i=1; i<=Nr; i++)
	{
		double ri=rmin+(rmax-rmin)/Nr*i;
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
			int pbin=floor((pin+0.001-pmin)/dp);
			wigdst_test.SetBinContent(i,pbin,wdstin);
		}
		f.close();
	}
	TFile *save = new TFile("pnOmega_qmax3_lmax2_Lmax0_angnum8_nmax5_nstart2_2600-test-scaled.root","recreate");
	wigdst_test.Write();
	save->Close();
	delete save;
}
