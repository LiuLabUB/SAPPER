
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <cstdlib>
#include <sstream>
#include "io.hpp"

using namespace std;



int main(int argc,char *argv[])
{
	if(argc<8)
	{
		cout<<"\nProgram: SNVAS (call SNV and allele-specific events from ChIP-seq data)\n";
		cout<<"Version: 0.1\n";
		cout<<"Usage: SNVAS <peaks.bed> <peaks.bam> <controlpeaks.bam> {PE,SE} <output.vcf> <fermi> <tmp>\n\n";
		cout<<"Options: <peaks.bed>            sorted bed file of peak regions\n";
		cout<<"         <peaks.bam>            sorted bam file of peak regions\n";
		cout<<"         <controlpeaks.bam>     sorted control bam file of peak regions\n";
		cout<<"         {PE,SE}                format of bam file, PE (paired-end) or SE (single-end)\n";
		cout<<"         <output.vcf>           output vcf file\n";
		cout<<"         <fermi>                path of fermi executable file\n";
		cout<<"         <tmp>                  temparary file folder\n\n";
	}
	else
	{
		const string Peakbedfile(argv[1]);
		const string Bamfile(argv[2]);
		const string InputBamfile(argv[3]);
		const string PEorSE(argv[4]);
		const string OutputVcffile(argv[5]);
		const string fermi_location(argv[6]);
		const string tmpfilefolder(argv[7]);

		if(PEorSE!="PE" && PEorSE!="SE") {cout<<"wrong format {PE,SE}"<<endl;exit(0);}

		/*string sbuf,fermi_location,tmpfilefolder;
		ifstream is("config.txt");
		if(!is) {cout<<"can not open bedtools location file: config.txt"<<endl;exit(0);}
		is>>sbuf>>fermi_location;
		if(sbuf!="fermi:") {cout<<"wrong format of config.txt"<<endl;exit(0);}
		if(fermi_location.empty()) {cout<<"wrong format of config.txt"<<endl;exit(0);}
		sbuf.clear();
		is>>sbuf>>tmpfilefolder;
		if(sbuf!="tmp_filefolder:") {cout<<"wrong format of config.txt"<<endl;exit(0);}
		if(tmpfilefolder.empty()) {cout<<"wrong format of config.txt"<<endl;exit(0);}
		is.close();*/

		mkdir(tmpfilefolder.c_str(), 0755);

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
cout<<"finishall"<<endl;
	}
}
