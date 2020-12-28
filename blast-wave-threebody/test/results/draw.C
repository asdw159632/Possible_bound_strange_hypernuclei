void draw()
{

	TFile res("result-all.root","read");
	TH1D *proton=(TH1D*)res->Get("proton_pT_Dst");
	TH1D *neutron=(TH1D*)res->Get("neutron_pT_Dst");
	TH1D *Omega=(TH1D*)res->Get("Omega_pT_Dst");
	TH1D *cluster=(TH1D*)res->Get("cluster_pT_Dst");

	TH2D *blank=new TH2D("blank","blank",50,0,5,100,1.e-3,1.e8);

	proton->SetMarkerStyle(kOpenCircle);
	proton->SetMarkerColor(9);
	neutron->SetMarkerStyle(kOpenSquare);
	neutron->SetMarkerColor(8);
	Omega->SetMarkerStyle(kOpenTriangleUp);
	Omega->SetMarkerColor(4);
	cluster->SetMarkerStyle(kOpenTriangleDown);
	cluster->SetMarkerColor(2);

	TCanvas *c1 = new TCanvas("c1","c1",10,10,800,600);
	gPad->SetLogy();

	blank->GetYaxis()->SetNdivisions(9,5,9);
	//proton->GetYaxis()->SetNdivisions(10,5,10);

	blank->Draw();
	proton->Draw("SAME E0 X0");
	neutron->Draw("SAME E0 X0");
	Omega->Draw("SAME E0 X0");
	cluster->Draw("SAME E0 X0");
	TFile *save = new TFile("picture.root","recreate");
	c1->Write();
	save->Close();
}
