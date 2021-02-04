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
		ku=1/aBorhr/fmmev;
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

	char readpath0[520];
	if(scale==1)
	{
		sprintf(readpath0,"%s/%s_thetal_0-scaled.root",workdir,workdir);
	}else{
		sprintf(readpath0,"%s/%s_thetal_0.root",workdir,workdir);
	}
	TFile read0(readpath0,"read");
	TH2D *readth0=(TH2D*)read0.Get("wigdst");
	rmin=readth0->GetXaxis()->GetBinLowEdge(readth0->GetXaxis()->GetFirst());
	rmax=readth0->GetXaxis()->GetBinUpEdge(readth0->GetXaxis()->GetLast());
	dr=readth0->GetXaxis()->GetBinWidth(1);

	pmin=readth0->GetYaxis()->GetBinLowEdge(readth0->GetYaxis()->GetFirst());
	pmax=readth0->GetYaxis()->GetBinUpEdge(readth0->GetYaxis()->GetLast());
	dp=readth0->GetYaxis()->GetBinWidth(1);
	read0.Close();

	TH3D *wigdst=new TH3D("wigdst","wigdst",Nr,rmin,rmax,Np,pmin,pmax,Ntheta,thetamin,thetamax);
	
	int process=0;
	for(int l=0;l<Ntheta;l++)
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
		TFile read(readpath,"read");
		TH2D *readth=(TH2D*)read.Get("wigdst");
		for(int i=1;i<=Nr;i++)
		{
			for(int j=1;j<=Np;j++)
			{
				double dst=readth->GetBinContent(i,j);
				wigdst->SetBinContent(i,j,l+1,dst);
			}
		}
		read.Close();
	}

	TFile create(savepath,"recreate");
	wigdst->Write();
	cout<<"Program complete"<<endl;
}
