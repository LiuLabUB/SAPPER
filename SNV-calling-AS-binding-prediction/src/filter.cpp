#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <unistd.h>
#include "SNVAS_common.hpp"

using namespace std;

int main_filter(int argc,char *argv[])
{
  int c;
  int depthcutoff = 20;
  int GQcutoff = 50;
  string type, ifile, ofile;

  while ((c = getopt(argc, argv, "i:d:q:t:o:")) >= 0) {
    switch (c) {
    case 'i': ifile.assign(optarg); break;
    case 'd': depthcutoff = atoi(optarg); break;
    case 'q': GQcutoff = atoi(optarg); break;
    case 't': type.assign(optarg); break;
    case 'o': ofile.assign(optarg); break;
    }
  }

  if ( ifile=="" or ofile=="" or type=="" ) {
    cerr<<"Program: SNVAS filter -- apply cutoff to filter VCF file from SNVAS call\n";
    cerr<<"Version: 0.1\n";
    cerr<<"Contacts: Liqing Tian <liqingti@buffalo.edu> & Tao Liu <tliu4@buffalo.edu>\n";
    cerr<<"Usage: SNVAS filter <-i SNVAS.vcf> <-o output.vcf> <-t genotype> [-d depthCutoff] [-q GQCutoff]\n\n";
    cerr<<"Required arguments:\n"
	<<"    <-i SNVAS.vcf>          The raw output VCF file of SNV calling from SNVAS\n"
	<<"    <-o output.vcf>         The VCF file containing SNVs after filtering\n"
	<<"    <-t homo|hetero|hetero_AS|hetero_nonAS>   Genotypes chosen from:\n"
	<<"                                    homo: homozygous SNVs\n"
	<<"                                    hetero: heterozygous SNVs\n"
	<<"                                    hetero_AS: heterozygous SNVs with allele-specific binding\n"
	<<"                                    hetero_nonAS: heterozygous SNV with non allele-specific binding\n"
	<<"Options:\n"
	<<"    [-d depthCutoff]        Keep the SNVs with read depth >= depthCutoff (Default:20). Must be a positive integer\n" 
	<<"    [-q GQCutoff]           Genotype quality cutoff. (Recommend: 50 for heterozygous SNVs and 10 for\n"
	<<"                            homozygous SNVs. Default: 50) Must be a positive integer.\n";
    return 1;
  }
  string sbuf,INFO,s1,s2,s3,s4,s5,s6,s7,s8,s9;

  if(type!="homo" && type!="hetero" && type!="hetero_AS" && type!="hetero_nonAS") {
    cerr<<"Wrong genotype! It should be either homo|hetero|hetero_AS|hetero_nonAS."<<endl;return(1);
  }

  if(depthcutoff<0) {
    cerr<<"Wrong depthCutoff, which must be an integer >=0>"<<endl;return(1);
  }
  if(GQcutoff<0){
    cerr<<"Wrong GQCutoff, which must be an integer >=0>"<<endl;return(1);
  }
  
  ifstream is(ifile.c_str());
  if(!is) {cerr<<"Can not open "<<ifile<<endl;return 1;}

  ofstream os(ofile.c_str());
  if(!os) {cerr<<"Can not open "<<ofile<<endl;return 1;}

  do
    {
      sbuf.clear();
      is>>sbuf;
      if(sbuf.empty()) break;
      
      if(sbuf[0]=='#') {os<<sbuf;getline(is,sbuf);os<<sbuf<<endl;continue;}
      
      s1=sbuf;
      is>>s2>>s3>>s4>>s5>>s6>>s7>>INFO>>s8>>s9;
      
      SNVASinfor snvastemp;
      ExtractInfor(INFO,snvastemp.predicttype,snvastemp.top2NT,snvastemp.top2No,snvastemp.GQ,snvastemp.GQ_ASsig,snvastemp.top1plus,snvastemp.top2plus,snvastemp.top1minus,snvastemp.top2minus,snvastemp.heter_AS_alleleratio,snvastemp.ChIPseq_depth,snvastemp.input_depth);
      
      //cout<<s1<<":"<<s2<<"\t"<<snvastemp.raw_ChIPseq_depth<<" "<<snvastemp.raw_input_depth<<" "<<depthcutoff<<endl;
      if(snvastemp.ChIPseq_depth+snvastemp.input_depth<depthcutoff) continue;
      
      if(type=="homo") {
	if(snvastemp.predicttype=="homo") {
	  if(snvastemp.GQ>=GQcutoff) os<<s1<<"\t"<<s2<<"\t"<<s3<<"\t"<<s4<<"\t"<<s5<<"\t"<<s6<<"\t"<<s7<<"\t"<<INFO<<"\t"<<s8<<"\t"<<s9<<endl;
	}
      }
      else if(type=="hetero") {
	if(snvastemp.predicttype=="heter_noAS" || snvastemp.predicttype=="heter_AS") {
	  if(snvastemp.GQ>=GQcutoff) os<<s1<<"\t"<<s2<<"\t"<<s3<<"\t"<<s4<<"\t"<<s5<<"\t"<<s6<<"\t"<<s7<<"\t"<<INFO<<"\t"<<s8<<"\t"<<s9<<endl;
	}
      }
      else if(type=="hetero_AS") {
	if(snvastemp.predicttype=="heter_AS") {
	  if(snvastemp.GQ>=GQcutoff) os<<s1<<"\t"<<s2<<"\t"<<s3<<"\t"<<s4<<"\t"<<s5<<"\t"<<s6<<"\t"<<s7<<"\t"<<INFO<<"\t"<<s8<<"\t"<<s9<<endl;
	}
      }
      else if(type=="hetero_nonAS") {
	if(snvastemp.predicttype=="heter_noAS") {
	  if(snvastemp.GQ>=GQcutoff) os<<s1<<"\t"<<s2<<"\t"<<s3<<"\t"<<s4<<"\t"<<s5<<"\t"<<s6<<"\t"<<s7<<"\t"<<INFO<<"\t"<<s8<<"\t"<<s9<<endl;
	}
      }
      
    } while(!is.eof());
  is.close();
  return 0;
}

