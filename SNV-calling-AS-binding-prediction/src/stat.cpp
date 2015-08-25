#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <map>
#include <unistd.h>
#include "SNVAS_common.hpp"

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

int Write_stat(string outputfile,const vector<double> &score_set,const vector<string> &tsortv_set,const int peak_length)
{

  ofstream os(outputfile);
  if(!os) {cerr<<"Can not open "<<outputfile<<endl;return 1;}
  
  os<<"GQCutoff\tSNVsPerKb\tTsTv"<<endl;

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
		  else {cerr<<"error!! tsortv "<<tsortv_set[j]<<endl;exit(0);}

		  all++;
		}
	    }

	  os<<cutoff<<"\t"<<(double)all*1000.0/peak_length<<"\t"<<(double)(ts+0.001)/(tv+0.001)<<endl;
	}
    }
  os.close();
  return 0;
}

int main_stat(int argc,char *argv[])
{
  int c;
  int depthcutoff = 20;
  string vcffile, heterofile, homofile, PeakRegionFile;

  while ((c = getopt(argc, argv, "i:b:e:o:d:")) >= 0) {
    switch (c) {
    case 'i': vcffile.assign(optarg); break;
    case 'd': depthcutoff = atoi(optarg); break;
    case 'b': PeakRegionFile.assign(optarg); break;
    case 'e': heterofile.assign(optarg); break;
    case 'o': homofile.assign(optarg); break;
    }
  }

  if ( vcffile=="" or PeakRegionFile=="" or heterofile=="" or homofile=="" ) {
    cerr<<"Program: SNVAS stat -- Statistics of genotype quality cutoff from predicted SNVs by SNVAS\n";
    cerr<<"Version: 0.1\n";
    cerr<<"Contacts: Liqing Tian <liqingti@buffalo.edu> & Tao Liu <tliu4@buffalo.edu>\n";
    cerr<<"Usage: SNVAS stat <-i SNVAS.vcf> <-b peaks.bed> <-e hetero_stat.txt> <-o homo_stat.txt >[-d depthCutoff]\n\n";
    cerr<<"Required arguments:\n"
	<<"    <-i SNVAS.vcf>          The raw output VCF file of SNV calling from SNVAS\n"
	<<"    <-b peaks.bed>          The BED/narrowPeak/BroadPeak file for peak regions, which is used in SNV calling by SNVAS\n"
	<<"    <-e hetero_stat.txt>    The output cutoff statistics file for predicted heterozygous SNVs\n"
	<<"                              the 1st column: genotype quality cutoff\n"
	<<"                              the 2nd column: density of predicted heterozygous SNVs per kbp\n"
	<<"                              the 3rd column: ts/tv ratio of predicted heterozygous SNVs\n"
	<<"    <-o homo_stat.txt>      The output cutoff statistics file for predicted homozygous SNVs\n"
	<<"                              the 1st column: genotype quality cutoff\n"
	<<"                              the 2nd column: density of predicted heterozygous SNVs per kbp\n"
	<<"                              the 3rd column: ts/tv ratio of predicted heterozygous SNVs\n\n"
	<<"Options:\n"
	<<"    [-d depthCutoff]        Keep the SNVs with read depth >= depthCutoff (Default:20). Must be a positive integer\n";
    return 1;
  }

  if(depthcutoff<0) {cerr<<"Wrong depthCutoff, which must be an integer >=0>"<<endl;return 1;}

  //read peak region bed file
  map<string,vector<int> > chr2start;
  map<string,vector<int> > chr2end;

  int peak_length=CalPeakRegionAndLength(PeakRegionFile,chr2start,chr2end);

  //read my predicted vcf
  vector<string> pred01set,pred11set;//01 means heter, include 0/1 1/2 types
  vector<double> score01set,score11set;
  vector<string> tsortv01set,tsortv11set;//ts,tv,skip

  ReadSNVAS_result(vcffile,depthcutoff,chr2start,chr2end,pred01set,pred11set,score01set,score11set,tsortv01set,tsortv11set);

  //
  Write_stat(heterofile,score01set,tsortv01set,peak_length);

  //
  Write_stat(homofile,score11set,tsortv11set,peak_length);

  return 0;
}

