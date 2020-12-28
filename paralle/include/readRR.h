THnD *RRMTHn;

const int readRR(int L)
{
	char path[100];
	sprintf(path,"/home/zhangliang/Possible_bound_strange_hypernuclei/RRCoeff/RR_Matrix/RR@{qmax=12,lmax=6,L=%d,%s}.root",L,nuclear);
	TFile f(path,"read");
	if(f.IsZombie())
	{
		cout<<"Zombie"<<endl;
		return 1;
	}
	RRMTHn=(THnD*)f.Get("RRcoeff");
	f.Close();
	return 0;
}

double RRM(int i, int j, int qi, int lxi, int lyi, int qj, int lxj)
{
	int idx[]={i,j,qi,lxi,lyi,qj,lxj};
	return RRMTHn->GetBinContent(idx);
}
