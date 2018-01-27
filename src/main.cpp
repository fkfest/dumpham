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
  // options
  std::vector<std::string> options; 
  // position of the option in the argv list
  std::vector<int> argpos;
  for ( int iarg = 1; iarg < argc; ++iarg ) {
    if ( argv[iarg][0]=='-' ) {
      // get options
      if ( strlen(argv[iarg]) > 2 && argv[iarg][1]=='-' ){
        // long option "--word"
        options.push_back(std::string(&argv[iarg][1]));
        argpos.push_back(iarg);
      } else {
        // short options "-wo"
        for ( uint i = 1; i < strlen(argv[iarg]); ++i ){
          options.push_back(std::string(1,argv[iarg][i]));
          argpos.push_back(iarg);
        }
      }
    }
  }  
  std::string inputfile, outputfile, orboutputfile,
    exePath = exepath();
  bool fcidump = false;
  bool orbdump = false;
  // handle options  
  for ( uint iopt = 0; iopt < options.size(); ++iopt ) {
    const std::string & opt = options[iopt];
    int iarg = argpos[iopt];
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
      ++iarg;
      if ( iarg < argc && str2num<int>(Input::verbose,argv[iarg],std::dec)){
        argpos.push_back(iarg);
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
  // not-a-option arguments
  std::vector<std::string> args;
  for ( int iarg = 1; iarg < argc; ++iarg ) {
    if ( std::find(argpos.begin(),argpos.end(),iarg) == argpos.end() ){
      args.push_back(argv[iarg]);
    }
  }
  if ( args.size() == 0 ) error("Please provide an input file!");
  uint iarg = 0;
  inputfile=args[iarg]; ++iarg;
  if ( iarg < args.size() ) {
    outputfile=args[iarg]; ++iarg;
  } else 
    outputfile = FileName(inputfile,true)+"_NEW.FCIDUMP";
  if ( orbdump && iarg < args.size() ) {
    orboutputfile=args[iarg]; ++iarg;
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
