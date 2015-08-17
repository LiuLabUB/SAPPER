#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>

using namespace std;

struct SNVASinfor
{
	string chr;
	int pos;
	char ref;
	string top2NT;
	vector<int> top2No;
	string predicttype;
	int GQ;//genotype quality for homozygous or heterozygous
	int GQ_ASsig;//genotype quality for heterozygous with AS compared to no AS, only for predicttype="heter_AS"
	double heter_AS_alleleratio;
	int top1plus,top2plus,top1minus,top2minus;
};

void ExtractInfor(const string INFO,string &predicttype,string &top2NT,vector<int> &top2No,int &GQ,int &GQ_ASsig,int &top1plus,int &top2plus,int &top1minus,int &top2minus,double &heter_AS_alleleratio)
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
		if(vstemp[i].substr(0,5)=="top1:") {top2NT.push_back(vstemp[i][vstemp[i].size()-1]);top2No.push_back(atoi(vstemp[i].substr(5,vstemp[i].size()-6).c_str()));}
		else if(vstemp[i].substr(0,5)=="top2:") {top2NT.push_back(vstemp[i][vstemp[i].size()-1]);top2No.push_back(atoi(vstemp[i].substr(5,vstemp[i].size()-6).c_str()));}
		else if(vstemp[i].substr(0,4)=="GQ_h")
		{
			if(vstemp[i].substr(0,7)=="GQ_homo") GQ_homo=atof(vstemp[i].substr(8).c_str());
			else if(vstemp[i].size()>=11 && vstemp[i].substr(0,11)=="GQ_heter_no") GQ_heter_noAS=atoi(vstemp[i].substr(14).c_str());
			else if(vstemp[i].size()>=12 && vstemp[i].substr(0,12)=="GQ_heter_AS:") GQ_heter_AS=atoi(vstemp[i].substr(12).c_str());
			else if(vstemp[i].size()>=14 && vstemp[i].substr(0,14)=="GQ_heter_ASsig") GQ_ASsig=atoi(vstemp[i].substr(15).c_str());
		}
		else if(vstemp[i].substr(0,9)=="top1plus:") top1plus=atoi(vstemp[i].substr(9).c_str());
		else if(vstemp[i].substr(0,9)=="top2plus:") top2plus=atoi(vstemp[i].substr(9).c_str());
		else if(vstemp[i].substr(0,10)=="top1minus:") top1minus=atoi(vstemp[i].substr(10).c_str());
		else if(vstemp[i].substr(0,10)=="top2minus:") top2minus=atoi(vstemp[i].substr(10).c_str());
		else if(vstemp[i].substr(0,22)=="Allele_ratio_heter_AS:") heter_AS_alleleratio=atof(vstemp[i].substr(22).c_str());
	}

	if(predicttype=="homo") GQ=GQ_homo;
	else if(predicttype=="heter_noAS") GQ=GQ_heter_noAS;
	else if(predicttype=="heter_AS") GQ=GQ_heter_AS;
	else {cout<<"wrong predicted type: "<<predicttype<<endl;exit(0);}
}

int main(int argc,char *argv[])
{
	if ( argc < 5 )
	{
		cout<<"Program: SNVAS_filter (filter vcf file of SNV calling from SNVAS)\n";
		cout<<"Version: 0.1\n";
		cout<<"Contacts: Liqing Tian <liqingti@buffalo.edu> & Tao Liu <tliu4@buffalo.edu>\n";
		cout<<"Usage: SNVAS_filter <snv.vcf> <homo/hete/heterAS/heterNonAS> <cutoff> <snv_afterfilter.vcf>\n\n";
		cout<<"Options: <snv.vcf>                            raw output vcf file of SNV calling from SNVAS\n";
		cout<<"         <homo/hete/heter_AS/heter_noAS>      homo: get homozygous SNV after filtering\n"
			<<"                                              hete: get heterozygous SNV after filtering\n"
			<<"                                              heter_AS: get heterozygous SNV with allele-specific binding after filtering\n"
			<<"                                              heter_noAS: get heterozygous SNV with non allele-specific binding after filtering\n";
		cout<<"         <cutoff>                             genotype quality cutoff for SNV filtering, \n"
			<<"                                              which should be integer more than 0\n";
		cout<<"         <snv_afterfilter.vcf>                output vcf file after filtering\n\n";
		return(1);
	}

	//
	string sbuf,INFO,s1,s2,s3,s4,s5,s6,s7,s8,s9;

	const string type(argv[2]);
	if(type!="homo" && type!="hete" && type!="heter_AS" && type!="heter_noAS")
	{
		cout<<"wrong argv[2], which should be <homo/hete/heter_AS/heter_noAS>"<<endl;return(1);
	}

	const int GQcutoff=atoi(argv[3]);

	ifstream is(argv[1]);
	if(!is) {cout<<"can not open "<<argv[1]<<endl;exit(0);}

	ofstream os(argv[4]);
	if(!os) {cout<<"can not open "<<argv[4]<<endl;exit(0);}

	do
	{
		sbuf.clear();
		is>>sbuf;
		if(sbuf.empty()) break;

		if(sbuf[0]=='#') {os<<sbuf;getline(is,sbuf);os<<sbuf<<endl;continue;}

		s1=sbuf;
		is>>s2>>s3>>s4>>s5>>s6>>s7>>INFO>>s8>>s9;

		SNVASinfor snvastemp;
		ExtractInfor(INFO,snvastemp.predicttype,snvastemp.top2NT,snvastemp.top2No,snvastemp.GQ,snvastemp.GQ_ASsig,snvastemp.top1plus,snvastemp.top2plus,snvastemp.top1minus,snvastemp.top2minus,snvastemp.heter_AS_alleleratio);

		if(type=="homo")
		{
			if(snvastemp.predicttype=="homo")
			{
				if(snvastemp.GQ>=GQcutoff) os<<s1<<"\t"<<s2<<"\t"<<s3<<"\t"<<s4<<"\t"<<s5<<"\t"<<s6<<"\t"<<s7<<"\t"<<INFO<<"\t"<<s8<<"\t"<<s9<<endl;
			}
		}
		else if(type=="hete")
		{
			if(snvastemp.predicttype=="heter_noAS" || snvastemp.predicttype=="heter_AS")
			{
				if(snvastemp.GQ>=GQcutoff) os<<s1<<"\t"<<s2<<"\t"<<s3<<"\t"<<s4<<"\t"<<s5<<"\t"<<s6<<"\t"<<s7<<"\t"<<INFO<<"\t"<<s8<<"\t"<<s9<<endl;
			}
		}
		else if(type=="heter_AS")
		{
			if(snvastemp.predicttype=="heter_AS")
			{
				if(snvastemp.GQ>=GQcutoff) os<<s1<<"\t"<<s2<<"\t"<<s3<<"\t"<<s4<<"\t"<<s5<<"\t"<<s6<<"\t"<<s7<<"\t"<<INFO<<"\t"<<s8<<"\t"<<s9<<endl;
			}
		}
		else if(type=="heter_noAS")
		{
			if(snvastemp.predicttype=="heter_noAS")
			{
				if(snvastemp.GQ>=GQcutoff) os<<s1<<"\t"<<s2<<"\t"<<s3<<"\t"<<s4<<"\t"<<s5<<"\t"<<s6<<"\t"<<s7<<"\t"<<INFO<<"\t"<<s8<<"\t"<<s9<<endl;
			}
		}

	}while(!is.eof());
	is.close();
}

