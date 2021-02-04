void normtest()
{
	TFile f("nOmega.root","read");
	TH3D *rho=(TH3D*)f.Get("wigdst");

	int imax=rho->GetNbinsX();
	int jmax=rho->GetNbinsY();
	int lmax=rho->GetNbinsZ();
	
	double dpi=2*3.1415926535897932384626;

	double integ=0;
	for(int i=1;i<=imax;i++)
	{
		double r=rho->GetXaxis()->GetBinCenter(i);
		double dr=rho->GetXaxis()->GetBinWidth(i);
		for(int j=1;j<=jmax;j++)
		{
			double k=rho->GetYaxis()->GetBinCenter(j);
			double dk=rho->GetYaxis()->GetBinWidth(j);
			for(int l=1;l<=lmax;l++)
			{
				double th=rho->GetZaxis()->GetBinCenter(l);
				double dt=rho->GetZaxis()->GetBinWidth(l);

				double rho_rk=rho->GetBinContent(i,j,l);
				integ+=rho_rk*2*pow(dpi,2)*r*r*k*k*dr*dk*sin(th)*dt/pow(dpi,3)*pow(5.06/1000,3);
			}
		}
	}
	cout<<integ<<endl;
}
