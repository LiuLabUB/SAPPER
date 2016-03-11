#include <string>
#include <vector>
#include <cstdlib>
#include <cstring>
#include <unistd.h>
#include "io.hpp"

using namespace std;

int main_call(int argc, char *argv[])
{
  int c;
  string Peakbedfile, Bamfile, InputBamfile, OutputVcffile, fermi_location;
  double Fermi_overlap_minpercent=0.5;//fermi min match should be >= read_length*Fermi_overlap_percent
  double top2nt_minpercent=0.8;//after contig filtering (only keep the reads consistent with any contig), top1nt_readsNo+top2nt_readsNo should be >= all_readsNo*top2nt_min_percent

  while ((c = getopt(argc, argv, "b:t:c:o:p:q:")) >= 0) {
    switch (c) {
    case 'b': Peakbedfile.assign(optarg); break;
    case 't': Bamfile.assign(optarg); break;
    case 'c': InputBamfile.assign(optarg); break;
    case 'o': OutputVcffile.assign(optarg); break;
    case 'p': Fermi_overlap_minpercent=atof(optarg); break;
    case 'q': top2nt_minpercent=atof(optarg); break;
    }
  }
  if ( Peakbedfile=="" or Bamfile=="" or InputBamfile=="" or OutputVcffile=="" ) {
    cerr<<"Program: SAPPER call -- Call all possible SNVs from ChIP-Seq\n";
    cerr<<"Version: 0.1\n";
    cerr<<"Contacts: Liqing Tian <liqingti@buffalo.edu> & Tao Liu <tliu4@buffalo.edu>\n";
    cerr<<"Usage: SAPPER call <-b peaks.bed> <-t peaksChIP.bam> <-c peaksControl.bam> <-o output.vcf> [-p FermiOverlapMinPercent] [-q top2ntMinPercent]\n\n";
    cerr<<"Required arguments: <-b peaks.bed>                  sorted bed file for peak regions\n";
    cerr<<"                    <-t peaksChIP.bam>              sorted bam file for peak regions\n";
    cerr<<"                    <-c peaksControl.bam>           sorted control bam file for peak regions\n";
    cerr<<"                    <-o output.vcf>                 output vcf file\n\n";
    cerr<<"Options:\n"
    	<<"                    [-p FermiOverlapMinRatio]       The minimal ratio of a read to be overlapped during local assembly.\n"
    	<<"                                                    (Default:0.5 which means half of the read). Must be a float between 0 and 1\n"
    	<<"                    [-q top2ntMinRatio]             The reads for the top 2 most frequent nucleotides at a loci shouldn't be too\n"
	<<"                                                    few comparing to total reads mapped there. The minimum ratio is set here.\n"
	<<"                                                    (Default:0.8 which means at least 80% of reads contain the top 2 nucleostides).\n"
    	<<"                                                    Must be a float between 0.5 and 1\n\n";
    cerr<<"Tips to prepare your input files from ChIP-Seq IP and CTRL BAM files:\n*Note: You need to modify the following sample command lines.*\n\n";
    cerr<<"1. Clean the BAM files:\n";
    cerr<<"    $ samtools view -q 30 -F 4 -F 256 -F 2048 -b IP.bam -o IP_clean.bam\n";
    cerr<<"    $ samtools view -q 30 -F 4 -F 256 -F 2048 -b CTRL.bam -o CTRL_clean.bam\n";
    cerr<<"2. Sort the BAM file:\n";
    cerr<<"    $ samtools sort  IP_clean.bam  IP_clean_sorted\n";
    cerr<<"    $ samtools sort  CTRL_clean.bam  CTRL_clean_sorted\n";
    cerr<<"3. Peak calling (example is for paired-end data):\n";
    cerr<<"    $ macs2 callpeak -f BAMPE -t IP_clean_sort.bam -c CTRL_clean_sort.bam -n MyFactor\n";
    cerr<<"4. Sort peak file:\n";
    cerr<<"    $ sort -k1,1 -k2,2n MyFactor_peaks.narrowPeak > MyFactor_peaks.sorted.bed\n";
    cerr<<"5. Extract reads in peak regions:\n";
    cerr<<"    $ samtools view -b IP_clean_sorted.bam -L MyFactor_peaks.sorted.bed -o IP_peaks.bam\n";
    cerr<<"    $ samtools view -b CTRL_clean_sorted.bam -L MyFactor_peaks.sorted.bed -o CTRL_peaks.bam\n\n";
    cerr<<"To run SAPPER:\n\n";
    cerr<<"    $ SAPPER call -b MyFactor_peaks.sorted.bed -t IP_peaks.bam -c CTRL_peaks.bam -o MyFactor.vcf\n";
    return 1;
  }
  
  fermi_location = "";

  if(Fermi_overlap_minpercent<0 || Fermi_overlap_minpercent>1) {cerr<<"Wrong FermiOverlapMinPercent, which must be must be an float between 0 and 1"<<endl;return 1;}
  if(top2nt_minpercent<0.5 || top2nt_minpercent>1) {cerr<<"Wrong top2ntMinPercent, which must be must be an float between 0.5 and 1"<<endl;return 1;}

  // make temp dir 
  char ctemplate[] = "SAPPER.XXXXXX";
  char * ctmpfilefolder = mkdtemp( ctemplate );
  const string tmpfilefolder( ctmpfilefolder );

  //read peak region bed file
  vector<BedRegion> peakbedregion_set;//sorted by bam header chr order
  cout<<"* Start to read regions defined in BED file..."<<endl;
  ReadPeakBedFile(Peakbedfile,Bamfile,peakbedregion_set);
  cout<<" Finish reading BED file."<<endl;
  //read Input bam file
  vector<vector<BamInfor> > AllPeakInputBamInfor;
  cout<<"* Start to read information within given regions from input BAM file..."<<endl;
  ReadInputBamfile(peakbedregion_set,InputBamfile,AllPeakInputBamInfor);
  cout<<" Finish reading input BAM file."<<endl;

  //get the read length from the first 100 reads of bam file
  int ReadLength=GetReadLengthFromBamFile(Bamfile);

  //read bam file and calculate
  OutputVcfResultHasInput_header(OutputVcffile,argc,argv);
  cout<<"* Start to read information within given regions from ChIP BAM file..."<<endl;
  cout<<" Assemble Unitigs in each given region, call SNVs and write to VCF file..."<<endl;
  ReadBamfile(top2nt_minpercent,Fermi_overlap_minpercent,ReadLength,peakbedregion_set,Bamfile,AllPeakInputBamInfor,fermi_location,tmpfilefolder,OutputVcffile);
  cout<<"* All done!"<<endl;

  // remove temp dir
  remove( ctmpfilefolder );
  return 0;
}
