#ifndef GLOBAL_H 
#define GLOBAL_H

#include <iostream>
#include <sstream>
#include <map>
#include <list>
#include <cmath>
#include <assert.h>
#include "hdtypes.h"

#ifndef NDEBUG
#define _DEBUG
#endif

typedef long unsigned int lui;
typedef std::list< std::string > TParArray;
typedef std::map< std::string, std::string > TsPar;
typedef std::map< std::string, int > TiPar;
typedef std::map< std::string, double > TfPar;
typedef std::map< std::string, TParArray > TaPar;

typedef std::map< std::string, TsPar > TsParSet;
typedef std::map< std::string, TiPar > TiParSet;
typedef std::map< std::string, TfPar > TfParSet;
typedef std::map< std::string, TaPar > TaParSet;


namespace Numbers
{
  // "zero"
  static const double verysmall=1.e-16;
  static const double small=1.e-10;
  static const int big=1000;
}

class Return{
public:
  enum Vals{
    Done = 0,
    Delete = 1,
    Change_sign = 2,
    Repeat = 4
  };
  Return(Vals val = Done) : _val(val){};
  Return & operator +=(const Return& ret);
  char _val;
};

// global variables and functions from input
namespace Input
{
  // verbosity
  extern int verbose;
  // input-parameters
  //string parameters
  extern TsParSet sPars;
  //integer parameters
  extern TiParSet iPars;
  //float parameters
  extern TfParSet fPars;
  //array of string parameters
  extern TaParSet aPars;
}

#define xout std::cout
#define _xout(level,x) \
    { \
       if ( Input::verbose >= (level) ) \
       { xout << x; } \
    }
// output unless Input::verbose is below 0
#define _xout0(x) _xout(0,x)
// output unless Input::verbose is below 1
#define _xout1(x) _xout(1,x)
// output unless Input::verbose is below 2
#define _xout2(x) _xout(2,x)
// output unless Input::verbose is below 3
#define _xout3(x) _xout(3,x)

// print timing
#define _CPUtiming(what,start,end) xout << std::fixed << std::setprecision(2) << "CPU time " << what << 1000.0*(end-start)/CLOCKS_PER_SEC << " ms\n";

#endif
