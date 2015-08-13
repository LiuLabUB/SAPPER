
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
	if(argc<15)
	{
		cout<<"\nProgram: SNVAS (call SNV and allele-specific events from ChIP-seq data)\n";
		cout<<"Version: 0.1\n";
		cout<<"Usage: SNVAS -bed <peaks.bed> -bam <peaks.bam> -controlbam <controlpeaks.bam> -f {PE,SE} -o <output.vcf> -fermi <fermi> -tmp <tmp>\n\n";
		cout<<"Options: -bed            sorted bed file of peak regions\n";
		cout<<"         -bam            sorted bam file of peak regions\n";
		cout<<"         -controlbam     sorted control bam file of peak regions\n";
		cout<<"         -f              {PE,SE}, format of bam file, PE (paired-end) or SE (single-end)\n";
		cout<<"         -o              output vcf file\n";
		cout<<"         -fermi          path of fermi executable file\n";
		cout<<"         -tmp            temparary file folder\n\n";
	}
	else
	{
		const string Peakbedfile(argv[2]);
		const string Bamfile(argv[4]);
		const string InputBamfile(argv[6]);
		const string PEorSE(argv[8]);
		const string OutputVcffile(argv[10]);
		const string fermi_location(argv[12]);
		const string tmpfilefolder(argv[14]);

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
