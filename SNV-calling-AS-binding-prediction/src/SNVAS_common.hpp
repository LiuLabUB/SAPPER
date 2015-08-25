#ifndef SNVAS_COMMON_HPP
#define SNVAS_COMMON_HPP

#include <string>
#include <vector>

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

void ExtractInfor(const string INFO,string &predicttype,string &top2NT,vector<int> &top2No,int &GQ,int &GQ_ASsig,int &top1plus,int &top2plus,int &top1minus,int &top2minus,double &heter_AS_alleleratio,int &ChIPseq_depth,int &input_depth);

#endif
