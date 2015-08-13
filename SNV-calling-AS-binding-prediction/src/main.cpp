#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <cstdlib>
#include <sstream>
#include "io.hpp"

using namespace std;


std::string lookupFermi()
{
  std::vector<std::string> path_array;

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
      std::string temp = path_array[i];

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

  if(argc<6)
    {
      cout<<"\nProgram: SNVAS (call SNV and allele-specific events from ChIP-seq data)\n";
      cout<<"Version: 0.1\n";
      cout<<"Usage: SNVAS <peaks.bed> <peaks.bam> <controlpeaks.bam> {PE,SE} <output.vcf> <fermi> <tmp>\n\n";
      cout<<"Options: <peaks.bed>            sorted bed file of peak regions\n";
      cout<<"         <peaks.bam>            sorted bam file of peak regions\n";
      cout<<"         <controlpeaks.bam>     sorted control bam file of peak regions\n";
      cout<<"         {PE,SE}                format of bam file, PE (paired-end) or SE (single-end)\n";
      cout<<"         <output.vcf>           output vcf file\n";
    }
  else
    {
      const string Peakbedfile(argv[1]);
      const string Bamfile(argv[2]);
      const string InputBamfile(argv[3]);
      const string PEorSE(argv[4]);
      const string OutputVcffile(argv[5]);
      const string fermi_location = lookupFermi();

      //if (fermi_location == "")
      //return(0);

      // make temp dir 
      const string tmpfilefolder( mkdtemp( "/tmp/SNVAS.XXXXXX" ) );

      if(PEorSE!="PE" && PEorSE!="SE") {cout<<"wrong format {PE,SE}"<<endl;exit(0);}
      
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
    }
}
