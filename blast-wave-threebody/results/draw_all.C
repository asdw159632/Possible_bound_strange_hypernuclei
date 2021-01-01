void draw_all()
{
	//adjustable para.
	int rebin=5;

	double Pi=TMath::Pi();

	double exppT_p[100];
	double expproduct_p[100];
	double exppT_Omega[100];
	double expproduct_Omega[100];
	double exppT_t[100];
	double expproduct_t[100];
	double pT, product;
	int num_t=0;
	int num_p=0;
	int num_Omega=0;
	ifstream readdata;

	readdata.open("ALICE-PbPb-p.txt");
	//readdata.open("RHIC-AUAU-p.txt");
	while(readdata>>pT>>product)
	{
		exppT_p[num_p]=pT;
		expproduct_p[num_p]=product;
		num_p++;
	}
	readdata.close();
	TGraph *expdata_p=new TGraph(num_p,exppT_p,expproduct_p);

	expdata_p->SetMarkerStyle(kFullCircle);
	expdata_p->SetMarkerColor(1);
	expdata_p->SetMarkerSize(1);

	readdata.open("ALICE-PbPb-Omega.txt");
	while(readdata>>pT>>product)
	{
		exppT_Omega[num_Omega]=pT;
		expproduct_Omega[num_Omega]=product;
		num_Omega++;
	}
	readdata.close();
	TGraph *expdata_Omega=new TGraph(num_Omega,exppT_Omega,expproduct_Omega);

	expdata_Omega->SetMarkerStyle(kFullTriangleDown);
	expdata_Omega->SetMarkerColor(46);
	expdata_Omega->SetMarkerSize(1);

	readdata.open("ALICE-PbPb-He3.txt");
	//readdata.open("RHIC-AUAU-triton.txt");
	while(readdata>>pT>>product)
	{
		exppT_t[num_t]=pT;
		expproduct_t[num_t]=product;
		num_t++;
	}
	readdata.close();
	TGraph *expdata_t=new TGraph(num_t,exppT_t,expproduct_t);

	expdata_t->SetMarkerStyle(kFullSquare);
	expdata_t->SetMarkerColor(9);
	expdata_t->SetMarkerSize(1);

	TFile readroot1("PbPb2_76-pnOmega/result-all.root","read");
	TH1D *model_p=(TH1D*)readroot1.Get("proton_pT_Dst");
	TH1D *model_Omega=(TH1D*)readroot1.Get("Omega_pT_Dst");
	TH1D *model_pnOmega=(TH1D*)readroot1.Get("cluster_pT_Dst");
	TH1F *hevt1=(TH1F*)readroot1.Get("h1_event");

	TFile readroot2("PbPb2_76-H3/result-all.root","read");
	TH1D *model_t=(TH1D*)readroot2.Get("cluster_pT_Dst");
	TH1F *hevt2=(TH1F*)readroot2.Get("h1_event");

	model_p->Rebin(rebin);
	model_t->Rebin(rebin);
	model_Omega->Rebin(rebin);
	model_pnOmega->Rebin(rebin);

	int Nevent1=hevt1->GetBinContent(1);
	int Nevent2=hevt2->GetBinContent(1);
	TF1 *divevt=new TF1("divevt","1",0,5);

	model_p->Divide(divevt,Nevent1*rebin);
	model_Omega->Divide(divevt,Nevent1*100*1.9785995*rebin);//Omega at 2.76TeV
	model_t->Divide(divevt,Nevent2/pow(2*Pi,6)*rebin);//triton
	//model_pnOmega->Divide(divevt,Nevent1/pow(2*Pi,6)*100*1.8);//NNOmega at 200GeV
	model_pnOmega->Divide(divevt,Nevent1/pow(2*Pi,6)*100*1.9785995*rebin);//NNOmega at 2.76TeV

	model_p->SetMarkerStyle(kOpenCircle);
	model_p->SetMarkerColor(1);
	model_p->SetMarkerSize(1);

	model_t->SetMarkerStyle(kOpenSquare);
	model_t->SetMarkerColor(9);
	model_t->SetMarkerSize(1);

	model_Omega->SetMarkerStyle(kOpenTriangleDown);
	model_Omega->SetMarkerColor(46);
	model_Omega->SetMarkerSize(1);

	model_pnOmega->SetMarkerStyle(kOpenTriangleUp);
	model_pnOmega->SetMarkerColor(30);
	model_pnOmega->SetMarkerSize(1);

	TH2D *blank=new TH2D("backgrand","title",10,0.,5.,10,1.e-10,1.e2);
	
	TCanvas c1("c1","c1",20,10,800,600);
	gPad->SetLogy();

	blank->Draw();

	model_t->Draw("same e0 x0");
	model_p->Draw("same e0 x0");
	model_Omega->Draw("same e0 x0");
	model_pnOmega->Draw("same e0 x0");
	expdata_p->Draw("same P");
	expdata_Omega->Draw("same P");
	expdata_t->Draw("same P");

	TLegend *legend=new TLegend(0.7,0.85,0.6,0.9);
	legend->SetHeader("title");
	legend->AddEntry(model_p,"blwc proton");
	legend->AddEntry(expdata_p,"exp proton");
	legend->AddEntry(model_Omega,"blwc Omega");
	legend->AddEntry(expdata_Omega,"exp Omega");
	legend->AddEntry(model_t,"blwc triton");
	legend->AddEntry(expdata_t,"exp He3");
	legend->AddEntry(model_pnOmega,"blwc pnOmega");
	legend->Draw();

	TFile save("picture.root","recreate");
	c1.Write();
	save.Close();
}

