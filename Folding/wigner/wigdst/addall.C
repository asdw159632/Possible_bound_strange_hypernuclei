void addall()
{
	double fmmev=5.06/1000;
	double Pi=TMath::Pi();
	//plot infomation
	int scale=1;

	double aBorhr;
	const int Nr=50;
	double rmin;
	double rmax;
	double dr;

	double ku;
	const int Np=50;
	double pmin;
	double pmax;
	double dp;

	if(scale==1)
	{
		aBorhr=1;
		rmin=-8;
		rmax=2;
		dr=(rmax-rmin)/Nr;

		ku=1/aBorhr/fmmev;
		pmin=-8;
		pmax=1;
		dp=(pmax-pmin)/Np;
	}
	else
	{
		rmin=0;//fm
		rmax=10;//fm
		dr=(rmax-rmin)/Nr;//fm

		pmin=0;//MeV
		pmax=1000;//MeV
		dp=(pmax-pmin)/Np;//MeV
	}

	const int Ntheta=24;
	double thetamin=0;
	double thetamax=Pi;
	double dtheta=(thetamax-thetamin)/Ntheta;

	char savepath[500];
	char workdir[500];

	cout<<"WorkDir: "<<endl;
	cin>>workdir;

	if(scale==1)
	{
		sprintf(savepath,"%s-scaled.root",workdir);
	}else{
		sprintf(savepath,"%s.root",workdir);
	}
	TH3D *wigdst=new TH3D("wigdst","wigdst",Nr,rmin,rmax,Np,pmin,pmax,Ntheta,thetamin,thetamax);
	
	int process=0;
	for(int l=1;l<=Ntheta;l++)
	{
		if(process<(l*10/Ntheta))
		{
			process=l*10/Ntheta;
			cout<<"process: "<<process*10<<"%"<<endl;
		}
		char readpath[520];
		if(scale==1)
		{
			sprintf(readpath,"%s/%s_thetal_%d-scaled.root",workdir,workdir,l);
		}else{
			sprintf(readpath,"%s/%s_thetal_%d.root",workdir,workdir,l);
		}
		char wigname[100];
		sprintf(wigname,"wigdst_thetal_%d",l);
		TFile read(readpath,"read");
		TH2D *readth=(TH2D*)read.Get(wigname);
		for(int i=1;i<=Nr;i++)
		{
			for(int j=1;j<=Np;j++)
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
