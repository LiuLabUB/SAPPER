#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <cstdlib>
#include <sstream>
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

      temp = temp + "/" + "fermi";
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

int main(int argc,char *argv[])
{
  if ( argc < 6 )
    {
      cout<<"Program: SNVAS (call SNV and allele-specific events from ChIP-seq data)\n";
      cout<<"Version: 0.1\n";
      cout<<"Contacts: Liqing Tian <liqingti@buffalo.edu> & Tao Liu <tliu4@buffalo.edu>\n";
      cout<<"Usage: SNVAS <peaks.bed> <peaks.bam> <controlpeaks.bam> {PE,SE} <output.vcf>\n\n";
      cout<<"Options: <peaks.bed>            sorted bed file of peak regions\n";
      cout<<"         <peaks.bam>            sorted bam file of peak regions\n";
      cout<<"         <controlpeaks.bam>     sorted control bam file of peak regions\n";
      cout<<"         {PE,SE}                format of bam file, PE (paired-end) or SE (single-end)\n";
      cout<<"         <output.vcf>           output vcf file\n\n";
      cout<<"Tips to prepare your input files from ChIP-Seq IP and CTRL BAM files:\n*Note: You need to modify the following sample command lines.*\n\n";
      cout<<"1. Clean the BAM files:\n";
      cout<<"    $ samtools view -q 30 -F 4 -F 256 -F 2048 -b IP.bam -o IP_clean.bam\n";
      cout<<"    $ samtools view -q 30 -F 4 -F 256 -F 2048 -b CTRL.bam -o CTRL_clean.bam\n";
      cout<<"2. Sort the BAM file:\n";
      cout<<"    $ samtools sort  IP_clean.bam  IP_clean_sorted\n";
      cout<<"    $ samtools sort  CTRL_clean.bam  CTRL_clean_sorted\n";
      cout<<"3. Peak calling (example is for paired-end data):\n";
      cout<<"    $ macs2 callpeak -f BAMPE -t IP_clean_sort.bam -c CTRL_clean_sort.bam -n MyFactor\n";
      cout<<"4. Sort peak file:\n";
      cout<<"    $ sort -k1,1 -k2,2n MyFactor_peaks.narrowPeak > MyFactor_peaks.sorted.bed\n";
      cout<<"5. Extract reads in peak regions:\n";
      cout<<"    $ samtools view -b IP_clean_sorted.bam -L MyFactor_peaks.sorted.bed -o IP_peaks.bam\n";
      cout<<"    $ samtools view -b CTRL_clean_sorted.bam -L MyFactor_peaks.sorted.bed -o CTRL_peaks.bam\n\n";
      cout<<"To run SNVAS:\n\n";
      cout<<"    $ SNVAS  MyFactor_peaks.sorted.bed IP_peaks.bam CTRL_peaks.bam PE MyFactor.vcf\n";
      return(1);
    }
  
  const string Peakbedfile(argv[1]);
  const string Bamfile(argv[2]);
  const string InputBamfile(argv[3]);
  const string PEorSE(argv[4]);
  const string OutputVcffile(argv[5]);
  const string fermi_location(lookupFermi());

  if (fermi_location == "")
    return(0);

  // make temp dir 
  char ctemplate[] = "SNVAS.XXXXXX";
  char * ctmpfilefolder = mkdtemp( ctemplate );
  const string tmpfilefolder( ctmpfilefolder );

  // check PE/SE mode
  if(PEorSE!="PE" && PEorSE!="SE") {cout<<"wrong mode! Choose from {PE,SE}"<<endl;exit(0);}
      
  //read peak region bed file
  vector<BedRegion> peakbedregion_set;//sorted by bam header chr order
  ReadPeakBedFile(Peakbedfile,Bamfile,peakbedregion_set);
  cout<<"finish read bed"<<endl;
  //read Input bam file
  vector<vector<BamInfor> > AllPeakInputBamInfor;
  ReadInputBamfile(peakbedregion_set,InputBamfile,AllPeakInputBamInfor);
  cout<<"finish read input bam"<<endl;
  //read bam file and calculate
  OutputVcfResultHasInput_header(OutputVcffile);
  ReadBamfile(PEorSE,peakbedregion_set,Bamfile,AllPeakInputBamInfor,fermi_location,tmpfilefolder,OutputVcffile);
  cout<<"finish all"<<endl;

  // remove temp dir
  remove( ctmpfilefolder );
}
