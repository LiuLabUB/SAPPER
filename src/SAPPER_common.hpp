/* The MIT License
   Copyright (c) 2016 Liqing Tian <tianlq@gmail.com>
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

#ifndef SAPPER_COMMON_HPP
#define SAPPER_COMMON_HPP

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
