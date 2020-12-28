{
  gROOT->Reset();
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetFillColor(0);
  gStyle->SetTitleColor(1);
  gStyle->SetStatColor(0);

  gStyle->SetTextFont(132);
  gStyle->SetLabelFont(132,"X");
  gStyle->SetTitleFont(132,"X");
  gStyle->SetLabelFont(132,"Y");
  gStyle->SetTitleFont(132,"Y");
  gStyle->SetStatFont(132);

  //ALICE_Data/PRC_88_044910_2013_Fig4C_C0005.txt
  Double_t p_pTC0005[42];
  Double_t p_pTC0005Err[42]={0};
  Double_t p_dstC0005[42];
  Double_t p_dstC0005Err[42];

  char skip[512];
  int ipoint=0;
  Double_t xx, yy;
  Double_t Errs[2];
  ifstream pinC0005("../ampt-spectra/plot/ALICE_Data/PRC_88_044910_2013_Fig4C_C0005.txt");
  for(int i=0; i<10; i++) pinC0005.getline(skip,512);
  cout<<skip<<endl;
  while(pinC0005.good())
  {
    pinC0005>>xx>>skip>>skip>>yy>>Errs[0]>>skip>>Errs[1]>>skip>>skip>>skip>>skip>>skip>>skip;
    if(!pinC0005.good()) continue;

    p_pTC0005[ipoint] = xx;
    p_dstC0005[ipoint] = yy;
    p_dstC0005Err[ipoint] = sqrt(Errs[0]*Errs[0] + Errs[1]*Errs[1]);

    ipoint++;
  }
  pinC0005.close();

  TGraphErrors *gre_p_pTC0005 = new TGraphErrors(ipoint, p_pTC0005, p_dstC0005, p_pTC0005Err, p_dstC0005Err);
  gre_p_pTC0005->SetMarkerStyle(20);

  ipoint=0;
  Double_t p_pTC0510[42];
  Double_t p_pTC0510Err[42]={0};
  Double_t p_dstC0510[42];
  Double_t p_dstC0510Err[42];

  ifstream pinC0510("../ampt-spectra/plot/ALICE_Data/PRC_88_044910_2013_Fig4C_C0510.txt");
  for(int i=0; i<10; i++) pinC0510.getline(skip,512);
  cout<<skip<<endl;
  while(pinC0510.good())
  {
    pinC0510>>xx>>skip>>skip>>yy>>Errs[0]>>skip>>Errs[1]>>skip>>skip>>skip>>skip>>skip>>skip;
    if(!pinC0510.good()) continue;

    p_pTC0510[ipoint] = xx;
    p_dstC0510[ipoint] = yy;
    p_dstC0510Err[ipoint] = sqrt(Errs[0]*Errs[0] + Errs[1]*Errs[1]);

    ipoint++;
  }
  pinC0510.close();

  TGraphErrors *gre_p_pTC0510 = new TGraphErrors(ipoint, p_pTC0510, p_dstC0510, p_pTC0510Err, p_dstC0510Err);
  gre_p_pTC0510->SetMarkerStyle(21);

  Double_t p_pTC0010[42];
  Double_t p_pTC0010Err[42]={0};
  Double_t p_dstC0010[42];
  Double_t p_dstC0010Err[42];

  for(int i=0; i<ipoint; i++) 
  {
    p_pTC0010[i] = p_pTC0510[i];
    p_dstC0010[i] = p_dstC0005[i] + p_dstC0510[i];
    p_dstC0010Err[i] = sqrt(p_dstC0005Err[i]*p_dstC0005Err[i] + p_dstC0510Err[i]*p_dstC0510Err[i]);
  }
  TGraphErrors *gre_p_pTC0010 = new TGraphErrors(ipoint, p_pTC0010, p_dstC0010, p_pTC0010Err, p_dstC0010Err);
  gre_p_pTC0010->SetMarkerStyle(22);

  //ALICE_Data/PLB728_216_2014_Figure3_BD_C0010.txt
    ipoint=0;
  Double_t Omega_pTC0010[42];
  Double_t Omega_pTC0010Err[42]={0};
  Double_t Omega_dstC0010[42];
  Double_t Omega_dstC0010Err[42];

  ifstream OmegainC0010("../ampt-spectra/plot/ALICE_Data/PLB728_216_2014_Figure3_BD_C0010.txt");
  for(int i=0; i<10; i++) OmegainC0010.getline(skip,512);
  cout<<skip<<endl;
  while(OmegainC0010.good())
  {
    OmegainC0010>>xx>>skip>>skip>>yy>>Errs[0]>>skip>>Errs[1]>>skip>>skip>>skip>>skip>>skip>>skip;
    if(!OmegainC0010.good()) continue;

    Omega_pTC0010[ipoint] = xx;
    Omega_dstC0010[ipoint] = yy/(2. * TMath::Pi() * xx);
    Omega_dstC0010Err[ipoint] = sqrt(Errs[0]*Errs[0] + Errs[1]*Errs[1])/(2. * TMath::Pi() * xx);

    ipoint++;
  }
  OmegainC0010.close();

  TGraphErrors *gre_Omega_pTC0010 = new TGraphErrors(ipoint, Omega_pTC0010, Omega_dstC0010, Omega_pTC0010Err, Omega_dstC0010Err);
  gre_Omega_pTC0010->SetMarkerStyle(23);


  //ALICE_Data/PRC_93_024917_2016_Fig9.csv
  ipoint=0;
  Double_t deuteron_pTC0010[42];
  Double_t deuteron_pTC0010Err[42]={0};
  Double_t deuteron_dstC0010[42];
  Double_t deuteron_dstC0010Err[42];

  ifstream deuteroninC0010("../ampt-spectra/plot/ALICE_Data/PRC_93_024917_2016_Fig9.csv");
  for(int i=0; i<14; i++) deuteroninC0010.getline(skip,512);
  cout<<skip<<endl;
  while(deuteroninC0010.good() && ipoint<21)
  {
    deuteroninC0010>>xx>>skip>>skip>>yy>>Errs[0]>>skip>>Errs[1]>>skip;
    if(!deuteroninC0010.good()) continue;

    deuteron_pTC0010[ipoint] = xx;
    deuteron_dstC0010[ipoint] = yy/(1e3);
    deuteron_dstC0010Err[ipoint] = sqrt(Errs[0]*Errs[0] + Errs[1]*Errs[1])/(1e3);

    ipoint++;
  }
  deuteroninC0010.close();

  TGraphErrors *gre_deuteron_pTC0010 = new TGraphErrors(ipoint, deuteron_pTC0010, deuteron_dstC0010, deuteron_pTC0010Err, deuteron_dstC0010Err);
  gre_deuteron_pTC0010->SetMarkerStyle(24);


  //blast-wave
  TFile finO("test.root");
  TF1 *constant = new TF1("constant","1.",0,100);

  /*
  TH1D *h1_Omega = (TH1D *) finO.Get("Omega_pT_Dst");

  numEvt = 10000;
  h1_Omega->Divide(constant, numEvt);

  h1_Omega->SetLineColor(kRed);
*/

  TH1D *h1_proton = (TH1D *) finO.Get("proton_pT_Dst");

  numEvt = 10000.;
  h1_proton->Divide(constant, numEvt);

  h1_proton->SetLineColor(kBlue);

  TH1D *h1_deuteron = (TH1D *) finO.Get("cluster_pT_Dst");

  numEvt = 10000;
  h1_deuteron->Divide(constant, numEvt);

  h1_deuteron->SetMarkerStyle(20);
  h1_deuteron->SetLineColor(kBlue);
  h1_deuteron->SetMarkerColor(kBlue);

  TH2F *back;
  //back = new TH2F("back","back",1000,14,1100,10,-0.06,0.29);
  back = new TH2F("back","back",10,0,5,100,1e-6,50);
  back->GetXaxis()->SetNdivisions(505);
  back->GetYaxis()->SetNdivisions(505);
  back->GetXaxis()->CenterTitle(true);
  back->GetYaxis()->CenterTitle(true);
  back->GetXaxis()->SetTitleOffset(1.2);
  back->GetYaxis()->SetTitleOffset(1.2);
  back->GetXaxis()->SetLabelOffset(0.015);
  back->GetYaxis()->SetLabelOffset(0.015);
  back->GetYaxis()->SetTitle("dN^{2}/(2#pip_{T}dp_{T}dy");
  back->GetXaxis()->SetTitle("p_{T} GeV/c");
  back->GetXaxis()->SetTitleFont(132);
  back->GetYaxis()->SetTitleFont(132);
  back->GetXaxis()->SetLabelSize(0.07);
  back->GetYaxis()->SetLabelSize(0.07);
  back->GetXaxis()->SetTitleSize(0.07);
  back->GetYaxis()->SetTitleSize(0.07);
  //back->GetXaxis()->SetTickLength(0.0);
  //back->GetYaxis()->SetTickLength(0.0);

  TLegend *leg = new TLegend(0.2, 0.65, 0.59, 0.95,NULL,"brNDC");
  leg->SetName("leg");
  leg->SetBorderSize(0);
  leg->SetTextSize(0.07);
  leg->SetTextFont(132);
  leg->SetLineColor(0);
  leg->SetLineStyle(0);
  leg->SetLineWidth(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetHeader("6.73 TeV, ^{12}C+^{12}C, b=0fm");
  //leg->AddEntry(gre_C12_v2_pT_2PC,"Woods-Saxon","p");
  //leg->AddEntry(gre_C12Triangle_v2_pT_2PC,"Triangle","p");
  //leg->AddEntry(gre_C12Chain_v2_pT_2PC,"Chain","p");


  TCanvas *c1 = new TCanvas("c1", "c1",50,0,700,500);
  //c1->SetLogx(1);
  //c1->SetLogy(1);
  //  c1->SetGrid(1,1);
  //  c1->SetFillStyle(4000);
  //  c1->SetFrameFillStyle(4000);
  c1->SetBorderMode(0);
  c1->SetFillColor(10);
  c1->SetFrameFillColor(0);
  c1->SetFrameBorderMode(0);
  c1->SetLeftMargin(0.15);
  c1->SetTopMargin(0.05);
  c1->SetRightMargin(0.02);
  c1->SetBottomMargin(0.19);
  //c1->Divide(2,2,0.01,0.01);
  //c1->Divide(2,1,0.0,0.0);


  //Double_t padX[5] = {0.0,0.3,0.6,0.9,1.0};
  Double_t padX[2] = {0.0,1.};
  Double_t padY[2] = {0.0,1.0};

  TPad *pad[1][1];
  char padName[32];
  int xN=0;
  int yN=0;
  for(int i=0; i<1; i++)
  {
    for(int j=0; j<1; j++)
    {
	sprintf(padName,"pad%d%d",i,j);

	pad[i][j] = new TPad(padName,padName,padX[i],padY[j],padX[i+1],padY[j+1]);

	pad[i][j]->Draw();
	pad[i][j]->cd();
	pad[i][j]->Range(0,0,1,1);
	pad[i][j]->SetTickx(1);
	pad[i][j]->SetTicky(1);
	pad[i][j]->SetLogy(1);

	pad[i][j]->SetBottomMargin(0.18);
	pad[i][j]->SetLeftMargin(0.18);
	pad[i][j]->SetTopMargin(0.01);
	pad[i][j]->SetRightMargin(0.01);
//	if(j==0) pad[i][j]->SetTopMargin(0.);
//	if(j==1) pad[i][j]->SetBottomMargin(0.);

	back->Draw("c");

	gre_p_pTC0005->Draw("ep");
	gre_p_pTC0510->Draw("ep");
	gre_p_pTC0010->Draw("ep");

	gre_Omega_pTC0010->Draw("ep");

	//gre_p_blw->Draw("Line");

	//gre_Omega_blw->Draw("Line");

	gre_deuteron_pTC0010->Draw("ep");

//h1_Omega->Draw("same");
h1_proton->Draw("same");
h1_deuteron->Draw("same");

	pad[i][j]->Modified();
	c1->cd();
    }
  }

  //char saveName[128];
  //sprintf(saveName,"./figs/fig-C12-vn-pT.eps");
  //c1->SaveAs(saveName);
}
