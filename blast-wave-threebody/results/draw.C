void draw()
{
	double Pi=TMath::Pi();

	double exppT_p[100];
	double expproduct_p[100];
	double exppT_t[100];
	double expproduct_t[100];
	double pT, product;
	int num_t=0;
	int num_p=0;
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
/*
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
	expdata_t->SetMarkerSize(1);*/

	TFile readroot("result-all.root","read");
	TH1D *model_p=(TH1D*)readroot.Get("proton_pT_Dst");
	TH1D *model_t=(TH1D*)readroot.Get("cluster_pT_Dst");
	TH1F *hevt=(TH1F*)readroot.Get("h1_event");

	int Nevent=hevt->GetBinContent(1);
	TF1 *divevt=new TF1("divevt","1",0,5);

	model_p->Divide(divevt,Nevent);
	//model_t->Divide(divevt,Nevent/pow(2*Pi,6));//triton
	//model_t->Divide(divevt,Nevent/pow(2*Pi,6)*100*1.8);//NNOmega at 200GeV
	model_t->Divide(divevt,Nevent/pow(2*Pi,6)*100*1.9785995);//NNOmega at 2.76TeV

	model_p->SetMarkerStyle(kOpenCircle);
	model_p->SetMarkerColor(1);
	model_p->SetMarkerSize(1);

	model_t->SetMarkerStyle(kOpenSquare);
	model_t->SetMarkerColor(9);
	model_t->SetMarkerSize(1);

	TH2D *blank=new TH2D("backgrand","title",10,0.,5.,10,1.e-10,1.e2);
	
	TCanvas c1("c1","c1",20,10,800,600);
	gPad->SetLogy();

	blank->Draw();

	model_t->Draw("same e0 x0");
	model_p->Draw("same e0 x0");
	expdata_p->Draw("same P");
	//expdata_t->Draw("same P");

	TLegend *legend=new TLegend(0.7,0.85,0.6,0.9);
	legend->SetHeader("title");
	legend->AddEntry(model_p,"blwc proton");
	legend->AddEntry(expdata_p,"exp proton");
	//legend->AddEntry(model_t,"blwc triton");
	legend->AddEntry(model_t,"blwc pnOmega");
	//legend->AddEntry(expdata_t,"exp He3");
	legend->Draw();

	TFile save("picture.root","recreate");
	c1.Write();
	save.Close();
}

