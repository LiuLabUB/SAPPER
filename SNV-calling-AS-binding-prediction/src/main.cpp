#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include "io.hpp"

using namespace std;

int main_call(int argc, char *argv[]);
int main_filter(int argc, char *argv[]);
int main_stat(int argc, char *argv[]);

int main(int argc,char *argv[])
{
  int ret=0;
  
  if (argc == 1) {
    cerr<<"Program: SNVAS (call SNV and allele-specific events from ChIP-seq data)\n";
    cerr<<"Version: 0.1\n";
    cerr<<"Contacts: Liqing Tian <liqingti@buffalo.edu> & Tao Liu <tliu4@buffalo.edu>\n\n";
    cerr<<"Usage:   SNVAS <command> [arguments]\n\n";
    cerr<<"Command: call     Call all possible SNVs\n";
    cerr<<"         filter   Filter the VCF file from 'call' command with cutoffs\n";
    cerr<<"         stat     Calculate SNV density and Ts/Tv for different GQ cutoffs\n";
    cerr<<"\n";
    return 1;
  }

  if (strcmp(argv[1], "call") == 0) ret = main_call(argc-1, argv+1);
  else if (strcmp(argv[1], "filter") == 0) ret = main_filter(argc-1, argv+1);
  else if (strcmp(argv[1], "stat") == 0) ret = main_stat(argc-1, argv+1);
  else {
    cerr<<"Unrecognized command.\n";
    return -1;
  }
  return ret;
}

