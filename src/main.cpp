#include <string>
#include <cstring> 
#include <iostream>
#include <fstream>
#include <vector>
#include "utilities.h"
#include "globals.h"
#include "finput.h"
#include "hdump.h"
#include "odump.h"

int main(int argc, char **argv)
{
  // handle input and output
  int iarg=1;
  std::vector<std::string> options; 
  while ( iarg < argc && argv[iarg][0]=='-') {
    // get options
    if ( strlen(argv[iarg]) > 2 && argv[iarg][1]=='-' ){
      // long option "--word"
      options.push_back(std::string(&argv[iarg][1]));
    } else {
      // short options "-wo"
      for ( uint i = 1; i < strlen(argv[iarg]); ++i ){
        options.push_back(std::string(1,argv[iarg][i]));
      }
    }
    ++iarg;
  }  
  std::string inputfile, outputfile, orboutputfile,
    exePath = exepath();
  bool fcidump = false;
  bool orbdump = false;
  // handle options  
  for ( uint iopt = 0; iopt < options.size(); ++iopt ) {
    const std::string & opt = options[iopt];
    if ( opt == "h" || opt == "-help" ) {
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
    } else if ( opt == "v" || opt == "-verbose" ) {
      if ( iarg == argc-1 || !str2num<int>(Input::verbose,argv[iarg+1],std::dec)){
        Input::verbose = 1;
      } else {
        ++iarg;
      }
    } else if ( opt == "d" || opt == "-dump" ) {
      // the input file is an FCIDUMP file
      fcidump = true;
    } else if ( opt == "o" || opt == "-orbs" ) {
      // generate the corresponding orbital file as unity
      orbdump = true;
    } else {
      error("Unknown paratemer -"+opt);
    }
  }
  if (iarg >= argc) error("Please provide an input file!");
  inputfile=argv[iarg];
  ++iarg;
  if ( argc > iarg ) {
    outputfile=argv[iarg];
    ++iarg;
  } else 
    outputfile = FileName(inputfile,true)+"_NEW.FCIDUMP";
  if ( orbdump && argc > iarg ) {
    orboutputfile=argv[iarg];
    ++iarg;
  } else if ( orbdump ) 
    orboutputfile = FileName(inputfile,true)+"_NEW.ORBDUMP";
 
  // read input
  Finput finput(exePath);
  if ( Input::iPars["output"]["fcinamtoupper"] > 0 )
    outputfile = uppercase(outputfile);
  if ( orbdump && Input::iPars["output"]["orbnamtolower"] > 0 )
    orboutputfile = lowercase(orboutputfile);
  if (fcidump) {
    // de-symmetrize FCIDUMP
    Hdump dump(inputfile);
    dump.store(outputfile);
    Odump odump(dump.norb());
    odump.store(orboutputfile);
  } else {
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
//     std::ofstream fout;
//     fout.open(outputfile.c_str());

    if ( inp.size() == 0 ){
      say("Empty input file!");
      return 1;
    }
    for ( uint il = 0; il < inp.size(); ++il ){
      if ( finput.addline(inp[il]) ){
      }
    }
  }
  return 0;
}
