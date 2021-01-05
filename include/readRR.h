THnD *RRMTHn;

const int readRR(int L)
{
	char path[100];
	sprintf(path,"~/Possible_bound_strange_hypernuclei/RRCoeff/RR_Matrix/RR_qmax12_lmax6_L%d_%s.root",L,nuclear);
	TFile f(path,"read");
	if(f.IsZombie())
	{
		cout<<path<<"Zombie"<<endl;
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
