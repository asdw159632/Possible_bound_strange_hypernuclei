#include <fstream>
#include <istream>

using namespace std;

void transform()
// a program transform function_parameter.txt into function_parameter.root
{
	char filename[100];
	cout<<"Input *nuclear_qmax_lmax_Lmax_angnum_nmax_c*: ";
	cin>>filename;

	char partxtname[100];
	sprintf(partxtname,"./funpar-MMA/%s.txt",filename);
	ifstream fpar;
	fpar.open(partxtname);
	if(!fpar)
	{
		cerr<<"file "<<partxtname<<" not exist"<<endl;
		return ;
	}

	char listtxtname[100];
	sprintf(listtxtname,"./anglemomentlist-MMA/%s.txt",filename);
	ifstream flist;
	flist.open(listtxtname);
	if(!flist)
	{
		cerr<<"file "<<listtxtname<<" not exist"<<endl;
		return ;
	}

	char rootname[100];
	sprintf(rootname,"%s.root",filename);
	TFile save(rootname,"recreate");
	if(save.IsZombie())
	{
		cerr<<"file exist"<<endl;
		return;
	}

	//input state infomation
	double qmax;
	double lmax;
	double Lmax;
	double angnum;
	double nmax;
	double nstart;
	double c;
	cout<<"state information: "<<endl;
	cout<<"qmax=";
	cin>>qmax;
	cout<<"lmax=";
	cin>>lmax;
	cout<<"Lmax=";
	cin>>Lmax;
	cout<<"anglestate number=";
	cin>>angnum;
	cout<<"nmax=";
	cin>>nmax;
	cout<<"nstart=";
	cin>>nstart;
	cout<<"c=";
	cin>>c;
	cout<<endl;

	TH1D *stateinfo=new TH1D("stateinfo","state infomation",7,0,7);//record state infomation
	stateinfo->SetXTitle("qmax; lmax; Lmax; anglestate number; nmax; nstart; c");
	stateinfo->SetBinContent(1,qmax);//qmax;
	stateinfo->SetBinContent(2,lmax);//lmax;
	stateinfo->SetBinContent(3,Lmax);//Lmax;
	stateinfo->SetBinContent(4,angnum);//anglestate number;
	stateinfo->SetBinContent(5,nmax);//nmax;
	stateinfo->SetBinContent(6,nstart);//nstart;
	stateinfo->SetBinContent(7,c);//c;
	stateinfo->Write();
	delete stateinfo;

	//read function parameter
	double par;
	TH1D *funpar=new TH1D("function_parameter","function_parameter",angnum*nmax,0,1);
	int i=1;
	while(fpar>>par)
	{
		funpar->SetBinContent(i,par);
		i++;
	}
	funpar->Write();
	delete funpar;

	fpar.close();

	//read anglemomentlist
	int Nc;
	int q;
	int lx;
	int ly;
	int L;
	int sjk;
	int two_Sa;
	int tjk;

	char listhead[100];
	flist.getline(listhead,30);
	cout<<listhead<<endl;

	TH2I *anglemoment=new TH2I("anglemomentlist","anglemomentlist",(int)angnum,0,1,7,0,7);//x-axis for different state(Nc); y-axis for quantum number;
	anglemoment->SetXTitle("Nc");
	anglemoment->SetYTitle("quantum number");
	while(flist>>Nc>>q>>lx>>ly>>L>>sjk>>two_Sa>>tjk)
	{
		cout<<Nc<<" "<<q<<" "<<lx<<" "<<ly<<" "<<L<<" "<<sjk<<" "<<two_Sa<<" "<<tjk<<endl;;
		anglemoment->SetBinContent(Nc,1,q);
		anglemoment->SetBinContent(Nc,2,lx);
		anglemoment->SetBinContent(Nc,3,ly);
		anglemoment->SetBinContent(Nc,4,L);
		anglemoment->SetBinContent(Nc,5,sjk);
		anglemoment->SetBinContent(Nc,6,two_Sa);
		anglemoment->SetBinContent(Nc,7,tjk);
	}
	anglemoment->Write();
	delete anglemoment;

	flist.close();

	save.Close();
}
