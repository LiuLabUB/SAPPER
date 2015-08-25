#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <map>
#include <sstream>
#include <algorithm>

using namespace std;

struct SNVASinfor
{
	string chr;
	int pos;
	char ref;
	int ChIPseq_depth,input_depth;
	string top2NT;
	vector<int> top2No;
	string predicttype;
	int GQ;//genotype quality for homozygous or heterozygous
	int GQ_ASsig;//genotype quality for heterozygous with AS compared to no AS, only for predicttype="heter_AS"
	double heter_AS_alleleratio;
	int top1plus,top2plus,top1minus,top2minus;
};

void ExtractInfor(const string INFO,string &predicttype,string &top2NT,vector<int> &top2No,int &GQ,int &GQ_ASsig,int &top1plus,int &top2plus,int &top1minus,int &top2minus,double &heter_AS_alleleratio,int &ChIPseq_depth,int &input_depth)
{
	int i,start=0;

	vector<string> vstemp;
	for(i=0;i<INFO.size();i++)
	{
		if(INFO[i]==';')
		{
			vstemp.push_back(INFO.substr(start,i-start));
			start=i+1;
		}
	}
	vstemp.push_back(INFO.substr(start,INFO.size()-start));

	//for(i=0;i<vstemp.size();i++) cout<<vstemp[i]<<endl;

	//
	if(vstemp[0].substr(0,12)!="MinBIC_model") {cout<<"wrong INFO begin: "<<vstemp[0]<<endl;exit(0);}
	predicttype=vstemp[0].substr(13);

	int GQ_homo=0,GQ_heter_noAS=0,GQ_heter_AS=0;
	for(i=0;i<vstemp.size();i++)
	{
		if(vstemp[i].substr(0,3)=="DP_")
		{
			if(vstemp[i].substr(0,8)=="DP_ChIP=") ChIPseq_depth=atoi(vstemp[i].substr(8).c_str());
			else if(vstemp[i].substr(0,9)=="DP_input=") input_depth=atoi(vstemp[i].substr(9).c_str());
		}
		else if(vstemp[i].substr(0,5)=="top1=") {top2NT.push_back(vstemp[i][vstemp[i].size()-1]);top2No.push_back(atoi(vstemp[i].substr(5,vstemp[i].size()-6).c_str()));}
		else if(vstemp[i].substr(0,5)=="top2=") {top2NT.push_back(vstemp[i][vstemp[i].size()-1]);top2No.push_back(atoi(vstemp[i].substr(5,vstemp[i].size()-6).c_str()));}
		else if(vstemp[i].substr(0,4)=="GQ_h")
		{
			if(vstemp[i].substr(0,7)=="GQ_homo") GQ_homo=atof(vstemp[i].substr(8).c_str());
			else if(vstemp[i].size()>=11 && vstemp[i].substr(0,11)=="GQ_heter_no") GQ_heter_noAS=atoi(vstemp[i].substr(14).c_str());
			else if(vstemp[i].size()>=12 && vstemp[i].substr(0,12)=="GQ_heter_AS=") GQ_heter_AS=atoi(vstemp[i].substr(12).c_str());
			else if(vstemp[i].size()>=14 && vstemp[i].substr(0,14)=="GQ_heter_ASsig") GQ_ASsig=atoi(vstemp[i].substr(15).c_str());
		}
		else if(vstemp[i].substr(0,9)=="top1plus=") top1plus=atoi(vstemp[i].substr(9).c_str());
		else if(vstemp[i].substr(0,9)=="top2plus=") top2plus=atoi(vstemp[i].substr(9).c_str());
		else if(vstemp[i].substr(0,10)=="top1minus=") top1minus=atoi(vstemp[i].substr(10).c_str());
		else if(vstemp[i].substr(0,10)=="top2minus=") top2minus=atoi(vstemp[i].substr(10).c_str());
		else if(vstemp[i].substr(0,22)=="Allele_ratio_heter_AS=") heter_AS_alleleratio=atof(vstemp[i].substr(22).c_str());
	}



	if(predicttype=="homo") GQ=GQ_homo;
	else if(predicttype=="heter_noAS") GQ=GQ_heter_noAS;
	else if(predicttype=="heter_AS") GQ=GQ_heter_AS;
	else {cout<<"wrong predicted type: "<<predicttype<<endl;exit(0);}
}

int CalPeakRegionAndLength(const string PeakRegionFile,map<string,vector<int> > &chr2start,map<string,vector<int> > &chr2end)
{
	string sbuf;
	int itemp=0;

	ifstream is0(PeakRegionFile.c_str());
	if(!is0) {cout<<"error "<<PeakRegionFile<<endl;exit(0);}

	do
	{
		sbuf.clear();
		is0>>sbuf;
		if(sbuf.empty()) break;

		int ia,ib;
		is0>>ia>>ib;
		if(chr2start.find(sbuf)==chr2start.end())
		{
			vector<int> via,vib;
			via.push_back(ia+1);
			vib.push_back(ib);
			chr2start[sbuf]=via;
			chr2end[sbuf]=vib;
		}
		else
		{
			chr2start[sbuf].push_back(ia+1);
			chr2end[sbuf].push_back(ib);
		}

		itemp+=ib-ia;

		getline(is0,sbuf);

	}while(!is0.eof());
	is0.close();

	return itemp;
}


string CalTsOrTv_heter(const char ref,const char top1NT,const char top2NT,const string position)
{
	char alter;

	if(ref!='A' && ref!='T' && ref!='C' && ref!='G') {cout<<"wrong ref "<<position<<endl;exit(0);}

	if(top1NT==ref)
	{
		if(top2NT!='N') alter=top2NT;
		else {cout<<"top1NT is N "<<position<<endl;exit(0);}
	}
	else if(top2NT==ref)
	{
		if(top1NT!='N') alter=top1NT;
		else {cout<<"top2NT is N "<<position<<endl;exit(0);}
	}
	else return "skip";

	if(ref=='A' && alter=='G') return "ts";
	if(ref=='G' && alter=='A') return "ts";
	if(ref=='C' && alter=='T') return "ts";
	if(ref=='T' && alter=='C') return "ts";

	return "tv";
}

string CalTsOrTv_homo(const char ref,const char top1NT,const string position)
{
	char alter;

	if(ref!='A' && ref!='T' && ref!='C' && ref!='G') {cout<<"wrong ref "<<position<<endl;exit(0);}

	if(top1NT==ref) {cout<<"top1NT is ref "<<position<<endl;exit(0);}
	if(top1NT=='N') {cout<<"top1NT is N "<<position<<endl;exit(0);}

	alter=top1NT;

	if(ref=='A' && alter=='G') return "ts";
	if(ref=='G' && alter=='A') return "ts";
	if(ref=='C' && alter=='T') return "ts";
	if(ref=='T' && alter=='C') return "ts";

	return "tv";
}

bool WhetherIn(const int pos,const string chr,map<string,vector<int> > &chr2start,map<string,vector<int> > &chr2end)
{
	bool btemp=false;

	int i;

	if(chr2start.find(chr)==chr2start.end()) {cout<<"error!chr "<<chr<<endl;exit(0);}

	vector<int> via(chr2start[chr]);
	vector<int> vib(chr2end[chr]);

	for(i=0;i<via.size();i++)
	{
		if(pos>=via[i] && pos<=vib[i]) {btemp=true;break;}
	}

	return btemp;
}

void ReadSNVAS_result(const string filename,const int depthcutoff,map<string,vector<int> > &chr2start,map<string,vector<int> > &chr2end,vector<string> &pred01set,vector<string> &pred11set,vector<double> &score01set,vector<double> &score11set,vector<string> &tsortv01set,vector<string> &tsortv11set)
{
	string sbuf,INFO;
	string position;

	ifstream is(filename.c_str());
	if(!is) {cout<<"can not open "<<filename<<endl;exit(0);}

	do
	{
		sbuf.clear();
		is>>sbuf;
		if(sbuf.empty()) break;

		if(sbuf[0]=='#') {getline(is,sbuf);continue;}

		SNVASinfor snvastemp;

		snvastemp.chr=sbuf;
		is>>snvastemp.pos>>sbuf>>snvastemp.ref>>sbuf>>sbuf>>sbuf>>INFO;getline(is,sbuf);

		if(!WhetherIn(snvastemp.pos,snvastemp.chr,chr2start,chr2end)) continue;

		stringstream ss;
		ss<<snvastemp.chr<<":"<<snvastemp.pos;
		position=ss.str();

		ExtractInfor(INFO,snvastemp.predicttype,snvastemp.top2NT,snvastemp.top2No,snvastemp.GQ,snvastemp.GQ_ASsig,snvastemp.top1plus,snvastemp.top2plus,snvastemp.top1minus,snvastemp.top2minus,snvastemp.heter_AS_alleleratio,snvastemp.ChIPseq_depth,snvastemp.input_depth);

		if(snvastemp.ChIPseq_depth+snvastemp.input_depth<depthcutoff) continue;

		if(snvastemp.predicttype=="homo")
		{
			pred11set.push_back(position);
			score11set.push_back(snvastemp.GQ);
			tsortv11set.push_back(CalTsOrTv_homo(snvastemp.ref,snvastemp.top2NT[0],position));
		}
		else if(snvastemp.predicttype=="heter_noAS" || snvastemp.predicttype=="heter_AS")
		{
			pred01set.push_back(position);
			score01set.push_back(snvastemp.GQ);
			tsortv01set.push_back(CalTsOrTv_heter(snvastemp.ref,snvastemp.top2NT[0],snvastemp.top2NT[1],position));
		}

	}while(!is.eof());
	is.close();
}

void Output(char *outputfile,const vector<double> &score_set,const vector<string> &tsortv_set,const int peak_length)
{
	ofstream os(outputfile,ofstream::app);
	if(!os) {cout<<"can not open "<<outputfile<<endl;exit(0);}

	vector<double> vdtemp(score_set);

	sort(vdtemp.begin(),vdtemp.end());
	vdtemp.erase( unique( vdtemp.begin(), vdtemp.end() ), vdtemp.end() );

	int i,j;
	for(i=0;i<=vdtemp.size();i++)
	{
		if(i!=vdtemp.size())
		{
			double cutoff=vdtemp[i];
			int ts=0,tv=0,skip=0,all=0;

			for(j=0;j<score_set.size();j++)
			{
				if(score_set[j]>cutoff-1e-8)
				{
					if(tsortv_set[j]=="ts") ts++;
					else if(tsortv_set[j]=="tv") tv++;
					else if(tsortv_set[j]=="skip") skip++;
					else {cout<<"error!! tsortv "<<tsortv_set[j]<<endl;exit(0);}

					all++;
				}
			}

			os<<cutoff<<"\t"<<(double)all*1000.0/peak_length<<"\t"<<(double)(ts+0.001)/(tv+0.001)<<endl;
		}
	}
}

int main(int argc,char *argv[])
{
	if ( argc < 6 )
	{
		cout<<"Program: SNVAS_cutoff_statistics (Statistics of genotype quality cutoff from predicted SNVs by SNVAS)\n";
		cout<<"Version: 0.1\n";
		cout<<"Contacts: Liqing Tian <liqingti@buffalo.edu> & Tao Liu <tliu4@buffalo.edu>\n";
		cout<<"Usage: SNVAS_cutoff_statistics <snvas.vcf> <depth_cutoff> <peaks.bed> <hetero_cutoff_statistics.txt> <homo_cutoff_statistics.txt>\n\n";
		cout<<"Options: <snvas.vcf>                            The raw output vcf file of SNV calling from SNVAS\n";
		cout<<"         <depth_cutoff>                         Only show SNVs with read depth >=depth_cutoff (recommend:20)\n";
		cout<<"         <peaks.bed>                            The BED file for peak regions, which is used in SNV calling by SNVAS\n";
		cout<<"         <hetero_cutoff_statistics.txt>         The output cutoff statistics file for predicted heterozygous SNVs\n"
		    <<"                                                the 1st column: genotype quality cutoff\n"
		    <<"                                                the 2nd column: density of predicted heterozygous SNVs per kbp\n"
		    <<"                                                the 3rd column: ts/tv ratio of predicted heterozygous SNVs\n\n";
		cout<<"         <homo_cutoff_statistics.txt>           The output cutoff statistics file for predicted homozygous SNVs\n"
		    <<"                                                the 1st column: genotype quality cutoff\n"
		    <<"                                                the 2nd column: density of predicted homozygous SNVs per kbp\n"
		    <<"                                                the 3rd column: ts/tv ratio of predicted homozygous SNVs\n\n";
		return(1);
	}

	//
	const int depthcutoff=atoi(argv[2]);
	if(depthcutoff<0) {cout<<"wrong argv[2], which is an integer >=0>"<<endl;return(1);}

	//read peak region bed file
	map<string,vector<int> > chr2start;
	map<string,vector<int> > chr2end;

	const string PeakRegionFile(argv[3]);
	int peak_length=CalPeakRegionAndLength(PeakRegionFile,chr2start,chr2end);

	//read my predicted vcf
	vector<string> pred01set,pred11set;//01 means heter, include 0/1 1/2 types
	vector<double> score01set,score11set;
	vector<string> tsortv01set,tsortv11set;//ts,tv,skip

	const string vcffile(argv[1]);
	ReadSNVAS_result(vcffile,depthcutoff,chr2start,chr2end,pred01set,pred11set,score01set,score11set,tsortv01set,tsortv11set);

	//
	ofstream os_heter_statistics(argv[4]);
	if(!os_heter_statistics) {cout<<"can not open "<<argv[4]<<endl;exit(0);}
	os_heter_statistics<<"GQCutoff\tSNVsPerKb\tTsTv"<<endl;

	Output(argv[4],score01set,tsortv01set,peak_length);

	//
	ofstream os_homo_statistics(argv[5]);
	if(!os_homo_statistics) {cout<<"can not open "<<argv[5]<<endl;exit(0);}
	os_homo_statistics<<"GQCutoff\tSNVsPerKb\tTsTv"<<endl;

	Output(argv[5],score11set,tsortv11set,peak_length);
}

