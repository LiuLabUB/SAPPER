#include "analysis.hpp"

void CalModel_Homo(const vector<double> &top1qualprob,const vector<double> &inputtop1qualprob,const vector<double> &top2qualprob,const vector<double> &inputtop2qualprob,
		double &lnL,double &BIC)
{
	int i;
	lnL=0;

	for(i=0;i<top1qualprob.size();i++) lnL+=log(1-top1qualprob[i]);
	for(i=0;i<inputtop1qualprob.size();i++) lnL+=log(1-inputtop1qualprob[i]);

	for(i=0;i<top2qualprob.size();i++) lnL+=log(top2qualprob[i]);
	for(i=0;i<inputtop2qualprob.size();i++) lnL+=log(inputtop2qualprob[i]);

	BIC=-2*lnL;
}

void CalModel_Heter_noAS(const vector<double> &top1qualprob,const vector<double> &inputtop1qualprob,const vector<double> &top2qualprob,const vector<double> &inputtop2qualprob,
		double &lnL,double &BIC,int &kc,int &ki)
{
        //int i;

	//for kc
	int tnc=top1qualprob.size()+top2qualprob.size();
	double da;

	if(tnc==0) {da=0;kc=-1;}
	else da=GreedyMaxFunctionNoAS(top1qualprob.size(),top2qualprob.size(),tnc,top1qualprob,top2qualprob,kc);

	//for ki
	int tni=inputtop1qualprob.size()+inputtop2qualprob.size();
	double db;

	if(tni==0) {db=0;ki=-1;}
	else db=GreedyMaxFunctionNoAS(inputtop1qualprob.size(),inputtop2qualprob.size(),tni,inputtop1qualprob,inputtop2qualprob,ki);

	if(tnc==0) BIC=-2*da-2*db+log(double(tni));
	else if(tni==0) BIC=-2*da-2*db+log(double(tnc));
	else BIC=-2*da-2*db+log(double(tnc))+log(double(tni));

	lnL=da+db;
}

void CalModel_Heter_AS(const vector<double> &top1qualprob,const vector<double> &inputtop1qualprob,const vector<double> &top2qualprob,const vector<double> &inputtop2qualprob,
		double &lnL,double &BIC,int &kc,int &ki,double &AS_alleleratio)
{
        //int i;

	//for kc
	int tnc=top1qualprob.size()+top2qualprob.size();
	double da;

	if(tnc==0) {da=0;kc=-1;}
	else da=GreedyMaxFunctionAS(top1qualprob.size(),top2qualprob.size(),tnc,top1qualprob,top2qualprob,AS_alleleratio,kc);//AS_alleleratio=(double)kc/tnc

	//for ki
	int tni=inputtop1qualprob.size()+inputtop2qualprob.size();
	double db;

	if(tni==0) {db=0;ki=-1;}
	else db=GreedyMaxFunctionNoAS(inputtop1qualprob.size(),inputtop2qualprob.size(),tni,inputtop1qualprob,inputtop2qualprob,ki);

	if(tnc==0) BIC=-2*da-2*db+log(double(tni));
	else if(tni==0) BIC=-2*da-2*db+2*log(double(tnc));
	else BIC=-2*da-2*db+2*log(double(tnc))+log(double(tni));

	lnL=da+db;
}

void CalModel_Homo(const vector<double> &top1qualprob,const vector<double> &top2qualprob,
		double &lnL,double &BIC)
{
	int i;
	lnL=0;

	for(i=0;i<top1qualprob.size();i++) lnL+=log(1-top1qualprob[i]);
	for(i=0;i<top2qualprob.size();i++) lnL+=log(top2qualprob[i]);

	BIC=-2*lnL;
}

void CalModel_Heter_noAS(const vector<double> &top1qualprob,const vector<double> &top2qualprob,
		double &lnL,double &BIC,int &kc)
{
        //int i;

	//for kc
	int tnc=top1qualprob.size()+top2qualprob.size();
	double da;

	if(tnc==0) {da=0;kc=-1;}
	else da=GreedyMaxFunctionNoAS(top1qualprob.size(),top2qualprob.size(),tnc,top1qualprob,top2qualprob,kc);

	if(tnc==0) BIC=-2*da;
	else BIC=-2*da+log(double(tnc));

	lnL=da;
}

void CalModel_Heter_AS(const vector<double> &top1qualprob,const vector<double> &top2qualprob,
		double &lnL,double &BIC,int &kc,double &AS_alleleratio)
{
        //int i;

	//for kc
	int tnc=top1qualprob.size()+top2qualprob.size();
	double da;

	if(tnc==0) {da=0;kc=-1;}
	else da=GreedyMaxFunctionAS(top1qualprob.size(),top2qualprob.size(),tnc,top1qualprob,top2qualprob,AS_alleleratio,kc);//AS_alleleratio=(double)kc/tnc

	if(tnc==0) BIC=-2*da;
	else BIC=-2*da+2*log(double(tnc));

	lnL=da;
}

double GreedyMaxFunctionAS(const int m,const int n,const int tn,const vector<double> &me,const vector<double> &ne,double &ar,int &k)
{
	double dnew,dold,rold,rnew;
	int kold,knew;
	bool btemp=false;
	int k0;

	if(tn==1)
	{
		double dl=CallnFunction(m,n,tn,me,ne,0,0);
		double dr=CallnFunction(m,n,tn,me,ne,1,1);

		if(dl>dr) {k=0;return dl;}
		else {k=1;return dr;}
	}
	else if(m==0) k0=m+1;
	else if(m==tn) k0=m-1;
	else k0=m;

	double d0=CallnFunction(m,n,tn,me,ne,(double)k0/tn,k0);
	double d1l=CallnFunction(m,n,tn,me,ne,(double)(k0-1)/tn,k0-1);
	double d1r=CallnFunction(m,n,tn,me,ne,(double)(k0+1)/tn,k0+1);

	if(d0>d1l-1e-8 && d0>d1r-1e-8) {k=k0;ar=(double)k0/tn;return d0;}
	else if(d1l>d0)
	{
		dold=d1l;
		kold=k0-1;
		rold=(double)(k0-1)/tn;
		while(kold>=1)//when kold=1 still run, than knew=0 is the final run
		{
			knew=kold-1;

			rnew=(double)(knew)/tn;

			dnew=CallnFunction(m,n,tn,me,ne,rnew,knew);

			if(dnew-1e-8<dold) {btemp=true;break;}
			kold=knew;
			dold=dnew;
			rold=rnew;
		};

		if(btemp)//maximum L value is in [1,m-1];
		{
			k=kold;ar=rold;return dold;
		}
		else//L(k=0) is the max for [0,m-1]
		{
			k=kold;ar=rold;return dold;
		}
	}
	else if(d1r>d0)
	{
		dold=d1r;
		kold=k0+1;
		rold=(double)(k0+1)/tn;
		while(kold<=tn-1)//when kold=tn-1 still run, than knew=tn is the final run
		{
			knew=kold+1;

			rnew=(double)(knew)/tn;

			dnew=CallnFunction(m,n,tn,me,ne,rnew,knew);

			if(dnew-1e-8<dold) {btemp=true;break;}
			kold=knew;
			dold=dnew;
			rold=rnew;
		};

		if(btemp)//maximum L value is in [m+1,tn-1]
		{
			k=kold;ar=rold;return dold;
		}
		else//L(k=tn) is the max for [m+1,tn]
		{
			k=kold;ar=rold;return dold;
		}
	}
	else {cout<<"error in GreedyMaxFunctionAS"<<endl;exit(0);}
}

double GreedyMaxFunctionNoAS(const int m,const int n,const int tn,const vector<double> &me,const vector<double> &ne,int &k)
{
	double dnew,dold;
	int kold,knew;
	bool btemp=false;
	int k0;

	double bg_r=0.5;

	if(tn==1)
	{
		double dl=CallnFunction(m,n,tn,me,ne,bg_r,0);
		double dr=CallnFunction(m,n,tn,me,ne,bg_r,1);

		if(dl>dr) {k=0;return dl;}
		else {k=1;return dr;}
	}
	else if(m==0) k0=m+1;
	else if(m==tn) k0=m-1;
	else k0=m;

	double d0=CallnFunction(m,n,tn,me,ne,bg_r,k0);
	double d1l=CallnFunction(m,n,tn,me,ne,bg_r,k0-1);
	double d1r=CallnFunction(m,n,tn,me,ne,bg_r,k0+1);

	if(d0>d1l-1e-8 && d0>d1r-1e-8) {k=k0;return d0;}
	else if(d1l>d0)
	{
		dold=d1l;
		kold=k0-1;
		while(kold>=1)//when kold=1 still run, than knew=0 is the final run
		{
			knew=kold-1;
			dnew=CallnFunction(m,n,tn,me,ne,bg_r,knew);

			if(dnew-1e-8<dold) {btemp=true;break;}
			kold=knew;
			dold=dnew;
		};

		if(btemp)//maximum L value is in [1,m-1];
		{
			k=kold;return dold;
		}
		else//L(k=0) is the max for [0,m-1]
		{
			k=kold;return dold;
		}
	}
	else if(d1r>d0)
	{
		dold=d1r;
		kold=k0+1;
		while(kold<=tn-1)//when kold=tn-1 still run, than knew=tn is the final run
		{
			knew=kold+1;
			dnew=CallnFunction(m,n,tn,me,ne,bg_r,knew);

			if(dnew-1e-8<dold) {btemp=true;break;}
			kold=knew;
			dold=dnew;
		};

		if(btemp)//maximum L value is in [m+1,tn-1]
		{
			k=kold;return dold;
		}
		else//L(k=tn) is the max for [m+1,tn]
		{
			k=kold;return dold;
		}
	}
	else {cout<<"error in GreedyMaxFunctionNoAS"<<endl;exit(0);}
}

double CallnFunction(const int m,const int n,const int tn,const vector<double> &me,const vector<double> &ne,const double r,const int &k)
{
	int i;
	double da=0;

	if(r<1e-6 || r>1-1e-6) da+=0;
	else da+=k*log(r)+(tn-k)*log(1-r);

	if(k==0 || k==tn) da+=0;
	else if(k<=tn/2)
	{
		for(i=0;i<k;i++) da+=log(double(tn-i)/(k-i));
	}
	else
	{
		for(i=0;i<tn-k;i++) da+=log(double(tn-i)/(tn-k-i));
	}

	for(i=0;i<m;i++)
	{
		da+=log((1-me[i])*((double)k/tn)+me[i]*(1-(double)k/tn));
	}
	for(i=0;i<n;i++)
	{
		da+=log((1-ne[i])*(1-(double)k/tn)+ne[i]*((double)k/tn));
	}

	return da;
}

int CalGQscore(const double lnL1,const double lnL2,const double lnL3)
{
	long double L1=powl(2.7182818,lnL1-lnL1);
	long double L2=powl(2.7182818,lnL2-lnL1);
	long double L3=powl(2.7182818,lnL3-lnL1);

	if(L1>1) L1=1;
	if(L2>1) L2=1;
	if(L3>1) L3=1;
	if(L1<1e-110) L1=1e-110;
	if(L2<1e-110) L2=1e-110;
	if(L3<1e-110) L3=1e-110;

	long double sum=L1+L2+L3;
	long double tmp=(L2+L3)/sum;
	int GQ_score;
	if(tmp>1e-110) GQ_score=(int)(-4.34294*logl(tmp));
	else GQ_score=255;

	return GQ_score;
}

int CalGQ_heterASsig_score(const double lnL1,const double lnL2)
{
	long double L1=powl(2.7182818,lnL1-lnL1);
	long double L2=powl(2.7182818,lnL2-lnL1);

	if(L1>1) L1=1;
	if(L2>1) L2=1;
	if(L1<1e-110) L1=1e-110;
	if(L2<1e-110) L2=1e-110;

	long double sum=L1+L2;
	long double tmp=L2/sum;
	int ASsig_score;
	if(tmp>1e-110) ASsig_score=(int)(-4.34294*logl(tmp));
	else ASsig_score=255;

	return ASsig_score;
}


void CalTop2NT(const vector<vector<double> > &Qual_set,const vector<vector<double> > &InputQual_set,int &top1ntindex,int &top2ntindex,double &prob_top1top2)
{
	int i;
	int itotal=0;
	vector<int> vitemp(5);
	for(i=0;i<5;i++) {vitemp[i]=Qual_set[i].size()+InputQual_set[i].size();itotal+=vitemp[i];}

	top1ntindex=(int)(max_element(vitemp.begin(),vitemp.end())-vitemp.begin());

	vitemp[top1ntindex]=-1;

	top2ntindex=(int)(max_element(vitemp.begin(),vitemp.end())-vitemp.begin());

	prob_top1top2=double(Qual_set[top1ntindex].size()+InputQual_set[top1ntindex].size()+Qual_set[top2ntindex].size()+InputQual_set[top2ntindex].size())/itotal;
}

void CalTop2NT(const vector<vector<double> > &Qual_set,int &top1ntindex,int &top2ntindex,double &prob_top1top2)
{
	int i;
	int itotal=0;
	vector<int> vitemp(5);
	for(i=0;i<5;i++) {vitemp[i]=Qual_set[i].size();itotal+=vitemp[i];}

	top1ntindex=(int)(max_element(vitemp.begin(),vitemp.end())-vitemp.begin());

	vitemp[top1ntindex]=-1;

	top2ntindex=(int)(max_element(vitemp.begin(),vitemp.end())-vitemp.begin());

	prob_top1top2=double(Qual_set[top1ntindex].size()+Qual_set[top2ntindex].size())/itotal;
}


void CalLikelihoodGTprob(PosReadsInfor &myReadInfor)
{
	double prob_top1top2;

	CalTop2NT(myReadInfor.Qual_set,myReadInfor.InputQual_set,myReadInfor.top1ntindex,myReadInfor.top2ntindex,prob_top1top2);

	if(myReadInfor.Qual_set[myReadInfor.top1ntindex].size()+myReadInfor.InputQual_set[myReadInfor.top1ntindex].size()==0) {myReadInfor.filterout=true;return;}
	if(prob_top1top2<0.8 || NTindex2char(myReadInfor.top1ntindex)=='N') {myReadInfor.filterout=true;return;}
	if(find(myReadInfor.fermiNTs.begin(),myReadInfor.fermiNTs.end(),NTindex2char(myReadInfor.top1ntindex))==myReadInfor.fermiNTs.end() && find(myReadInfor.fermiNTs.begin(),myReadInfor.fermiNTs.end(),NTindex2char(myReadInfor.top2ntindex))==myReadInfor.fermiNTs.end()) {myReadInfor.hasfermiinfor=false;return;}

	myReadInfor.hasfermiinfor=true;

	vector<double> top1qualprob;
	vector<double> top2qualprob;
	vector<double> inputtop1qualprob;
	vector<double> inputtop2qualprob;

	if(find(myReadInfor.fermiNTs.begin(),myReadInfor.fermiNTs.end(),NTindex2char(myReadInfor.top1ntindex))!=myReadInfor.fermiNTs.end() && find(myReadInfor.fermiNTs.begin(),myReadInfor.fermiNTs.end(),NTindex2char(myReadInfor.top2ntindex))!=myReadInfor.fermiNTs.end())
	{
		top1qualprob=myReadInfor.Qual_set[myReadInfor.top1ntindex];
		top2qualprob=myReadInfor.Qual_set[myReadInfor.top2ntindex];
		inputtop1qualprob=myReadInfor.InputQual_set[myReadInfor.top1ntindex];
		inputtop2qualprob=myReadInfor.InputQual_set[myReadInfor.top2ntindex];
	}
	else if(find(myReadInfor.fermiNTs.begin(),myReadInfor.fermiNTs.end(),NTindex2char(myReadInfor.top1ntindex))!=myReadInfor.fermiNTs.end())
	{
		top1qualprob=myReadInfor.Qual_set[myReadInfor.top1ntindex];
		inputtop1qualprob=myReadInfor.InputQual_set[myReadInfor.top1ntindex];
	}
	else if(find(myReadInfor.fermiNTs.begin(),myReadInfor.fermiNTs.end(),NTindex2char(myReadInfor.top2ntindex))!=myReadInfor.fermiNTs.end())
	{
		if(myReadInfor.Qual_set[myReadInfor.top2ntindex].size()+myReadInfor.InputQual_set[myReadInfor.top2ntindex].size()==0) {myReadInfor.filterout=true;return;}

		int itemp=myReadInfor.top1ntindex;
		myReadInfor.top1ntindex=myReadInfor.top2ntindex;
		myReadInfor.top2ntindex=itemp;
		top1qualprob=myReadInfor.Qual_set[myReadInfor.top2ntindex];
		inputtop1qualprob=myReadInfor.InputQual_set[myReadInfor.top2ntindex];
	}

	CalModel_Homo(top1qualprob,inputtop1qualprob,top2qualprob,inputtop2qualprob,
			myReadInfor.lnL_homo_majar,myReadInfor.BIC_homo_majar);

	CalModel_Homo(top2qualprob,inputtop2qualprob,top1qualprob,inputtop1qualprob,
			myReadInfor.lnL_homo_minor,myReadInfor.BIC_homo_minor);

	CalModel_Heter_noAS(top1qualprob,inputtop1qualprob,top2qualprob,inputtop2qualprob,
			myReadInfor.lnL_heter_noAS,myReadInfor.BIC_heter_noAS,myReadInfor.heter_noAS_kc,myReadInfor.heter_noAS_ki);

	CalModel_Heter_AS(top1qualprob,inputtop1qualprob,top2qualprob,inputtop2qualprob,
			myReadInfor.lnL_heter_AS,myReadInfor.BIC_heter_AS,myReadInfor.heter_AS_kc,myReadInfor.heter_AS_ki,myReadInfor.heter_AS_alleleratio);

	//
	myReadInfor.GQ_homo_majar=0;
	myReadInfor.GQ_heter_noAS=0;
	myReadInfor.GQ_heter_AS=0;
	myReadInfor.GQ_heter_ASsig=0;



	if(myReadInfor.ref!=NTindex2char(myReadInfor.top1ntindex) && myReadInfor.BIC_homo_majar<myReadInfor.BIC_homo_minor && myReadInfor.BIC_homo_majar<myReadInfor.BIC_heter_noAS && myReadInfor.BIC_homo_majar<myReadInfor.BIC_heter_AS)
	{
		myReadInfor.type="homo";
		myReadInfor.GQ_homo_majar=CalGQscore(myReadInfor.lnL_homo_majar,myReadInfor.lnL_homo_minor,myReadInfor.lnL_heter_noAS);

	}
	else if(myReadInfor.BIC_heter_noAS<myReadInfor.BIC_homo_majar && myReadInfor.BIC_heter_noAS<myReadInfor.BIC_homo_minor && myReadInfor.BIC_heter_noAS<myReadInfor.BIC_heter_AS+1e-8)
	{
		myReadInfor.type="heter_noAS";
		myReadInfor.GQ_heter_noAS=CalGQscore(myReadInfor.lnL_heter_noAS,myReadInfor.lnL_homo_majar,myReadInfor.lnL_homo_minor);
	}
	else if(myReadInfor.BIC_heter_AS<myReadInfor.BIC_homo_majar && myReadInfor.BIC_heter_AS<myReadInfor.BIC_homo_minor && myReadInfor.BIC_heter_AS<myReadInfor.BIC_heter_noAS)
	{
		myReadInfor.type="heter_AS";
		myReadInfor.GQ_heter_AS=CalGQscore(myReadInfor.lnL_heter_AS,myReadInfor.lnL_homo_majar,myReadInfor.lnL_homo_minor);
		myReadInfor.GQ_heter_ASsig=CalGQ_heterASsig_score(myReadInfor.lnL_heter_AS,myReadInfor.lnL_heter_noAS);
	}
	else if(myReadInfor.ref==NTindex2char(myReadInfor.top1ntindex) && myReadInfor.BIC_homo_majar<myReadInfor.BIC_homo_minor && myReadInfor.BIC_homo_majar<myReadInfor.BIC_heter_noAS && myReadInfor.BIC_homo_majar<myReadInfor.BIC_heter_AS)
	{
		myReadInfor.type="homo_ref";

	}
	else myReadInfor.type="unsure";
}

void CalSNVAS(map<int,PosReadsInfor> &pos2Readsinfo)
{
	//calculate likelihood function and probability of different genotype
	map<int,PosReadsInfor>::iterator miPtmp;
	for(miPtmp=pos2Readsinfo.begin();miPtmp!=pos2Readsinfo.end();++miPtmp)//for each position
	{
		CalLikelihoodGTprob(miPtmp->second);
	}
}
