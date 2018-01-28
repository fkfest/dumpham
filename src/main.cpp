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

int main(int argc, char **argv)
{
  ArgPars args(argc,argv);
  std::string opt, arg;
  std::string inputfile, outputfile, orboutputfile,
    exePath = exepath();
  bool fcidump = false;
  bool orbdump = false;
  // handle options  
  while ( args.nextoption(opt) ) {
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
      if ( args.optarg(arg) && str2num<int>(Input::verbose,arg,std::dec)){
        args.markasoption();
      } else {
        Input::verbose = 1;
      }
//       xout << "Verbosity: " << Input::verbose << std::endl;
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
  
  if ( !args.nextremaining(arg) ) error("Please provide an input file!");
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
    // de-symmetrize FCIDUMP
    Hdump dump(inputfile);
    dump.store(outputfile);
    if ( orbdump ) {
      Odump odump(dump.norb());
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
