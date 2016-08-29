#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <cstdlib>
#include "SAPPER_common.hpp"

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

	//
	if(vstemp[0].substr(0,12)!="MinBIC_model") {cout<<"wrong INFO begin: "<<vstemp[0]<<endl;exit(0);}
	predicttype=vstemp[0].substr(13);

	int GQ_homo=0,GQ_heter_noAS=0,GQ_heter_AS=0;
	for(i=0;i<vstemp.size();i++)
	{
		if(vstemp[i].substr(0,4)=="raw_")
                {
                        if(vstemp[i].substr(0,15)=="raw_depth_ChIP=") ChIPseq_depth=atoi(vstemp[i].substr(15).c_str());
                        else if(vstemp[i].substr(0,16)=="raw_depth_input=") input_depth=atoi(vstemp[i].substr(16).c_str());
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

