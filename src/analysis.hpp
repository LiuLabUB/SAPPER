/* The MIT License
   Copyright (c) 2016 Liqing Tian <liqingti@buffalo.edu>
   
   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:
   
   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.
   
   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

#ifndef ANALYSIS_HPP
#define ANALYSIS_HPP

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <sstream>
#include <algorithm>
#include <cstdlib>
#include <cmath>

using namespace std;

struct PosReadsInfor
{
	char ref;
	bool filterout;//if true, do not output

	vector<vector<double> > Qual_set;//A C G T N, only stat for the reads consistent with any contig
	vector<vector<double> > InputQual_set;//A C G T N, only stat for the reads consistent with any contig
	vector<int> ChIP_rawNo;//A C G T N
	vector<int> Input_rawNo;//A C G T N
	//int raw_read_depth;//ChIP_rawNo[0].size()+ChIP_rawNo[1].size()+...+Input_rawNo[0].size()+Input_rawNo[1].size()+...

	//
	int top1ntindex,top2ntindex;
	double lnL_homo_majar,lnL_heter_AS,lnL_heter_noAS,lnL_homo_minor;
	double BIC_homo_majar,BIC_heter_AS,BIC_heter_noAS,BIC_homo_minor;
	int heter_noAS_kc,heter_noAS_ki;
	int heter_AS_kc,heter_AS_ki;
	double heter_AS_alleleratio;

	//
	int GQ_homo_majar,GQ_heter_noAS,GQ_heter_AS;//phred scale of prob by standard formular
	int GQ_heter_ASsig;//phred scale of prob, to measure the difference between AS and noAS

	//
	string type;

	//
	bool hasfermiinfor;//if no fermi bam overlap in the position, false; if fermi bam in the position GT: N, false; if anyone of top2NT is not in fermi GT NTs, false;
	vector<char> fermiNTs;
};


inline char NTindex2char(const int ntindex)
{
	if(ntindex==0) return 'A';
	else if(ntindex==1) return 'C';
	else if(ntindex==2) return 'G';
	else if(ntindex==3) return 'T';
	else if(ntindex==4) return 'N';
	else {cout<<"wrong nucleotide index: "<<ntindex<<endl;exit(0);}
}


void CalSNVAS(map<int,PosReadsInfor> &pos2Readsinfo,const double top2nt_minpercent);

//Calculate GQscore
void CalLikelihoodGTprob(PosReadsInfor &myReadInfor,const double top2nt_minpercent);

int CalGQscore(const double lnL1,const double lnL2,const double lnL3);
int CalGQ_heterASsig_score(const double lnL1,const double lnL2);


//has input
void CalModel_Homo(const vector<double> &top1qualprob,const vector<double> &inputtop1qualprob,const vector<double> &top2qualprob,const vector<double> &inputtop2qualprob,double &lnL,double &BIC);
void CalModel_Heter_noAS(const vector<double> &top1qualprob,const vector<double> &inputtop1qualprob,const vector<double> &top2qualprob,const vector<double> &inputtop2qualprob,double &lnL,double &BIC,int &kc,int &ki);
void CalModel_Heter_AS(const vector<double> &top1qualprob,const vector<double> &inputtop1qualprob,const vector<double> &top2qualprob,const vector<double> &inputtop2qualprob,double &lnL,double &BIC,int &kc,int &ki,double &AS_alleleratio);

void CalTop2NT(const vector<vector<double> > &Qual_set,const vector<vector<double> > &InputQual_set,int &top1ntindex,int &top2ntindex,double &prob_top1top2);

//no input
void CalModel_Homo(const vector<double> &top1qualprob,const vector<double> &top2qualprob,double &lnL,double &BIC);
void CalModel_Heter_noAS(const vector<double> &top1qualprob,const vector<double> &top2qualprob,double &lnL,double &BIC,int &kc);
void CalModel_Heter_AS(const vector<double> &top1qualprob,const vector<double> &top2qualprob,double &lnL,double &BIC,int &kc,double &AS_alleleratio);

void CalTop2NT(const vector<vector<double> > &Qual_set,int &top1ntindex,int &top2ntindex,double &prob_top1top2);

//core algorithm
double GreedyMaxFunctionAS(const int m,const int n,const int tn,const vector<double> &me,const vector<double> &ne,double &ar,int &k);
double GreedyMaxFunctionNoAS(const int m,const int n,const int tn,const vector<double> &me,const vector<double> &ne,int &k);
double CallnFunction(const int m,const int n,const int tn,const vector<double> &me,const vector<double> &ne,const double r,const int &k);


#endif
