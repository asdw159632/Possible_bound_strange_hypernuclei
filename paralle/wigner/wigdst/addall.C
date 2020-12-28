void addall()
{
	double fmmev=5.06/1000;
	double Pi=TMath::Pi();
	//plot infomation
	int scale=1;
	if(scale=1)
	{
		double aBorhr=1;
		const int Nr=50;
		double rmin=-8;
		double rmax=2;
		double dr=(rmax-rmin)/Nr;

		double ku=1/aBorhr/fmmev;
		const int Np=50;
		double pmin=-8;
		double pmax=1;
		double dp=(pmax-pmin)/Np;
	}
	else
	{
		const int Nr=50;
		double rmin=0;//fm
		double rmax=10;//fm
		double dr=(rmax-rmin)/Nr;//fm

		const int Np=50;
		double pmin=0;//MeV
		double pmax=1000;//MeV
		double dp=(pmax-pmin)/Np;//MeV
	}

	const int Ntheta=24;
	double thetamin=0;
	double thetamax=Pi;
	double dtheta=(thetamax-thetamin)/Ntheta;

	char savepath[100];
	char workdir[100];

	cout<<"WorkDir: "<<endl;
	cin>>workdir;
	sprintf(savepath,"%s/%s.root",workdir,workdir);
	TH3D *wigdst=new TH3D("wigdst","wigdst",Nr,rmin,rmax,Np,pmin,pmax,Ntheta,thetamin,thetamax);
	int process=0;
	for(int l=1;l<=Ntheta;l++)
	{
		if(process<(l*10/Ntheta))
		{
			process=l*10/Ntheta;
			cout<<"process: "<<process*10<<"%"<<endl;
		}
		char readpath[100];
		sprintf(readpath,"%s/%s_thetal_%d.root",workdir,workdir,l);
		char wigname[100];
		sprintf(wigname,"wigdst_thetal_%d",l);
		TFile read(readpath,"read");
		TH2D *readth=(TH2D*)read.Get(wigname);
		for(int i=1;i<=Nr;i++)
		{
			for(int j=1;j<Np;j++)
			{
				double dst=readth->GetBinContent(i,j);
				wigdst->SetBinContent(i,j,l,dst);
			}
		}
		read.Close();
	}

	TFile create(savepath,"recreate");
	wigdst->Write();
	cout<<"Program complete"<<endl;
}
