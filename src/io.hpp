#ifndef IO_HPP
#define IO_HPP

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <algorithm>
#include <sstream>
#include <cmath>
#include <cstring>
#include <zlib.h>
#include <stdint.h>
#include <sys/stat.h>
#include <cstdlib>
#include <ctime>
#include "swalign.hpp"
#include "analysis.hpp"

using namespace std;

/****
sam format flag
1	read paired
2	read mapped in proper pair
4	read unmapped
8	mate unmapped
16	read reverse strand
32	mate reverse strand
64	first in pair
128	second in pair
256	not primary alignment
512	read fails platform/vendor quality checks
1024	read is PCR or optical duplicate
2048	supplementary alignment
Usaually, 99/147 [3 + (-1) | 3 - (-2)]  83/163 [3 - (-1) | 3 + (-2)]
***/

struct BedRegion
{
	string chr;
	int start;//1-base
	int end;
};

struct BamInfor
{
	string readname;
	string chr;
	unsigned long start;//1-base
	unsigned long end;
	unsigned short firstsegment;//0: SE, 1: PE /1, 2 PE /2
	bool reversestrand;//true: -, false: +
	unsigned int mapq;
	string cigar_old;
	string cigar;//revised e.g. 20M81S --> 20M
	string seq;
	string bq;
	unsigned int nm;
	string md;
	string refseq;
};

int GetReadLengthFromBamFile(const string Bamfile);

void ReadPeakBedFile(const string Peakbedfile,const string Bamfile,vector<BedRegion> &peakbedregion_set);
void ReadInputBamfile(const vector<BedRegion> &peakbedregion_set,const string InputBamfile,vector<vector<BamInfor> > &AllPeakInputBamInfor);
void ReadBamfile(const double top2nt_minpercent,const double Fermi_overlap_minpercent,const int ReadLength,const vector<BedRegion> &peakbedregion_set,const string Bamfile,const vector<vector<BamInfor> > &AllPeakInputBamInfor,const string fermi_location,const string tmpfilefolder,const string OutputVcffile);
void GetReverseComplementary(const string &input,string &output);
void AssembleAndSNVAS(const double top2nt_minpercent,const double Fermi_overlap_minpercent,const int ReadLength,const string fermi_location,const string tmpfilefolder,const string OutputVcffile,const int PeakIndex,const string regionchr,const vector<BamInfor> &PeakBamInfor,const vector<BamInfor> &PeakInputBamInfor);
void MyFree(vector<BamInfor> &PeakBamInfor);

void FillChIPraw(const int snvpos,const char readnt,map<int,PosReadsInfor> &pos2Readsinfo);
void FillControlraw(const int snvpos,const char readnt,map<int,PosReadsInfor> &pos2Readsinfo);
void FillChIPQualInfor(const int snvpos,const char readnt,const char readbq,map<int,PosReadsInfor> &pos2Readsinfo);
void FillControlQualInfor(const int snvpos,const char readnt,const char readbq,map<int,PosReadsInfor> &pos2Readsinfo);

void FillContigInfor(const int snvpos,const char snvref,const vector<char> contigNTs,map<int,PosReadsInfor> &pos2Readsinfo);

double ChangeQualToProb(const char &qual);
bool ConsistentWithContig(const vector<int> &readscoor,const string &readsseq,vector<vector<int> > &contigcoor,const vector<string> &contigseq);

void GetReadSeqCoor(const string seq,const string bq,const int startpos,const string cigar,string &revisedseq,string &revisedbq,vector<int> &revisedcoor);

void OutputVcfResultHasInput(const string outputfile,const string regionchr,map<int,PosReadsInfor> &pos2Readsinfo);
void OutputVcfResultHasInput_header(const string outputfile,const int argc,char *argv[]);

#endif
