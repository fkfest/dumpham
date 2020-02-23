#include <string>
#include <cstring> 
#include <iostream>
#include <fstream>
#include <vector>
#include "argpars.h"
#include "utilities.h"
#include "globals.h"
#include "finput.h"
#include "hdump.h"
#include "odump.h"

using namespace HamDump;
using namespace ArgParser;

int main(int argc, char **argv)
{
  ArgPars args(argc,argv);
  std::string arg;
  std::string inputfile, outputfile, orboutputfile,
    exePath = exepath();
  bool fcidump = false;
  bool fcidump_out = true;
  bool orbdump = false;
  bool use_pgs = false;
  // handle options  
  while ( args.nextoption() ) {
    if ( args.check(ArgOpt("Verbosity level","v","-verbose" ))) {
      if ( args.optarg(arg) && str2num<int>(Input::verbose,arg,std::dec)){
        args.markasoption();
      } else {
        Input::verbose = 1;
      }
//       xout << "Verbosity: " << Input::verbose << std::endl;
    } else if ( args.check(ArgOpt("the input file is an FCIDUMP file","d","-dump")) ) {
      fcidump = true;
    } else if ( args.check(ArgOpt("don't write a new FCIDUMP file (for -d option)","n","-nodump")) ) {
      fcidump_out = false;
    } else if ( args.check(ArgOpt("generate the corresponding orbital file","o","-orbs")) ) {
      orbdump = true;
    } else if ( args.check(ArgOpt("generate files with point-group symmetry","s","-sym")) ) {
      use_pgs = true;
    } else if ( args.check(ArgOpt("print this help","h","-help")) ) {
      args.printhelp(xout,"dumpham [OPTIONS] <input-file> [<output-file>]",
                     "Dump various model Hamiltonians as FCIDUMP files");
//       // print README file if exists
//       std::ifstream readme;
//       readme.open((exePath+"README.md").c_str());
//       if (readme.is_open()) {
//         std::string line;
//         while (readme.good()) {
//           getline (readme,line);
//           xout << line << std::endl;
//         }
//       }
      return 0;
    } else {
      error("Unknown parameter -"+args.get_current_option());
    }
  }
  
  if ( !args.nextremaining(arg) ) error("Please provide an input file or -h option!");
  inputfile=arg;
  if ( args.nextremaining(arg) ) {
    outputfile=arg;
  } else 
    outputfile = FileName(inputfile,true)+"_NEW.FCIDUMP";
  if ( orbdump && args.nextremaining(arg) ) {
    orboutputfile=arg;
  } else if ( orbdump ) 
    orboutputfile = FileName(inputfile,true)+"_NEW.ORBDUMP";
  
  // read input
  Finput finput(exePath);
  if ( Input::iPars["output"]["fcinamtoupper"] > 0 )
    outputfile = uppercase(outputfile);
  if ( orbdump && Input::iPars["output"]["orbnamtolower"] > 0 )
    orboutputfile = lowercase(orboutputfile);
  if (fcidump) {
    if ( use_pgs ) {
      Input::iPars["ham"]["nosym"] = 0;
    } else { 
      // de-symmetrize FCIDUMP
      Input::iPars["ham"]["nosym"] = 1;
    }
    Hdump dump(inputfile);
    if ( fcidump_out ) {
      dump.read_dump();
      dump.store(outputfile);
    }
    if ( orbdump ) {
      Odump odump(dump.pgs());
      odump.store(orboutputfile);
    }
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
