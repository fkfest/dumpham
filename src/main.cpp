#include <string>
#include <cstring> 
#include <iostream>
#include <fstream>
#include <vector>
#include "utilities.h"
#include "globals.h"
#include "finput.h"
#include "hdump.h"

int main(int argc, char **argv)
{
  // handle input and output
  std::string inputfile, outputfile,
    exePath = exepath();
  int iarg=1;
  bool fcidump = false;
  while ( iarg<argc && argv[iarg][0]=='-') {// handle options
    if (strcmp(argv[iarg],"-h")==0 || strcmp(argv[iarg],"--help")==0) {
      xout << "dumpham <input-file> [<output-file>]" << std::endl;
      // print README file if exists
      std::ifstream readme;
      readme.open((exePath+"README.md").c_str());
      if (readme.is_open()) {
        std::string line;
        while (readme.good()) {
          getline (readme,line);
          xout << line << std::endl;
        }
      }
      return 0;
    } else if (strcmp(argv[iarg],"-v")==0 || strcmp(argv[iarg],"--verbose")==0) {
      if ( iarg == argc-1 || !str2num<int>(Input::verbose,argv[iarg+1],std::dec)){
        Input::verbose = 1;
      } else {
        ++iarg;
      }
    } else if (strcmp(argv[iarg],"-d")==0 || strcmp(argv[iarg],"--dump")==0) {
      // the input file is an FCIDUMP file
      fcidump = true;
    }
    ++iarg;
  }
  if (iarg >= argc) error("Please provide an input file!");
  inputfile=argv[iarg];
  if (argc>iarg+1)
    outputfile=argv[iarg+1];
  else 
    outputfile = FileName(inputfile,true)+"_NEW.FCIDUMP";
  
  if (fcidump) {
    // de-symmetrize FCIDUMP
    Hdump dump(inputfile);
    dump.store(outputfile);
  } else {
      
    // read input
    Finput finput(exePath);
    std::ifstream fin;
    fin.open(inputfile.c_str());
    // save input file
    std::vector< std::string > inp;
    if (fin.is_open()) {
      std::string line;
      while (fin.good()) {
        std::getline (fin,line);
        _xout1(line << std::endl);
        inp.push_back(line);
      }
    } else
      error("Bad input file!");
    fin.close();
    std::ofstream fout;
    fout.open(outputfile.c_str());

    if ( inp.size() == 0 ){
      say("Empty input file!");
      return 1;
    }

  }
  return 0;
}
