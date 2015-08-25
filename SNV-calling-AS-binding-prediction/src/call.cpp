#include <string>
#include <vector>
#include <cstdlib>
#include <cstring>
#include <unistd.h>
#include "io.hpp"

using namespace std;

string lookupFermi()
{
  vector<string> path_array;

  char *dup = strdup(getenv("PATH"));
  char *s = dup;
  char *p = NULL;

  do {
    p = strchr(s, ':');
    if (p != NULL) {
      p[0] = 0;
    }
    path_array.push_back(s);
    s = p+1;
  } while (p != NULL);

  for (int i = 0; i < path_array.size(); ++i)
    {
      string temp = path_array[i];

      temp = temp + "/" + "SNVAS_fermi";
      // let execv determine if it is executable
      // or you can do that here if required
      if ( access(temp.c_str(), F_OK) == 0)
	{
	  cout << "Executable fermi binary file found at " << temp << "\n";
	  return temp;
	}
    }

  cout << "Executable fermi binary file can't be found in PATH. Please make sure fermi is installed in the system!\n";
  return("");
}

int main_call(int argc, char *argv[])
{
  int c;
  string Peakbedfile, Bamfile, InputBamfile, OutputVcffile, fermi_location;

  while ((c = getopt(argc, argv, "b:t:c:o:")) >= 0) {
    switch (c) {
    case 'b': Peakbedfile.assign(optarg); break;
    case 't': Bamfile.assign(optarg); break;
    case 'c': InputBamfile.assign(optarg); break;
    case 'o': OutputVcffile.assign(optarg); break;
    }
  }
  if ( Peakbedfile=="" or Bamfile=="" or InputBamfile=="" or OutputVcffile=="" ) {
    cerr<<"Program: SNVAS call -- Call all possible SNVs from ChIP-Seq\n";
    cerr<<"Version: 0.1\n";
    cerr<<"Contacts: Liqing Tian <liqingti@buffalo.edu> & Tao Liu <tliu4@buffalo.edu>\n";
    cerr<<"Usage: SNVAS call <-b peaks.bed> <-t peaksChIP.bam> <-c peaksControl.bam> <-o output.vcf>\n\n";
    cerr<<"Required arguments: <-b peaks.bed>            sorted bed file of peak regions\n";
    cerr<<"                    <-t peaksChIP.bam>            sorted bam file of peak regions\n";
    cerr<<"                    <-c peaksControl.bam>     sorted control bam file of peak regions\n";
    cerr<<"                    <-o output.vcf>           output vcf file\n\n";
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
    cerr<<"To run SNVAS:\n\n";
    cerr<<"    $ SNVAS call -b MyFactor_peaks.sorted.bed -t IP_peaks.bam -c CTRL_peaks.bam -o MyFactor.vcf\n";
    return 1;
  }
  
  fermi_location = lookupFermi();
  if ( fermi_location == "" ) {
    cerr<<"The modified fermi 'SNVAS_fermi' can't be found in the system.\n";
    return 1;
  }

  // make temp dir 
  char ctemplate[] = "SNVAS.XXXXXX";
  char * ctmpfilefolder = mkdtemp( ctemplate );
  const string tmpfilefolder( ctmpfilefolder );

  //read peak region bed file
  vector<BedRegion> peakbedregion_set;//sorted by bam header chr order
  ReadPeakBedFile(Peakbedfile,Bamfile,peakbedregion_set);
  cout<<"finish read bed"<<endl;
  //read Input bam file
  vector<vector<BamInfor> > AllPeakInputBamInfor;
  ReadInputBamfile(peakbedregion_set,InputBamfile,AllPeakInputBamInfor);
  cout<<"finish read input bam"<<endl;

  //get the read length from the first 100 reads of bam file
  int ReadLength=GetReadLengthFromBamFile(Bamfile);

  //read bam file and calculate
  OutputVcfResultHasInput_header(OutputVcffile,argv);
  ReadBamfile(ReadLength,peakbedregion_set,Bamfile,AllPeakInputBamInfor,fermi_location,tmpfilefolder,OutputVcffile);
  cout<<"finish all"<<endl;

  // remove temp dir
  remove( ctmpfilefolder );
  return 0;
}
