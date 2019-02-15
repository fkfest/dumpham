#ifndef HDTYPES_H 
#define HDTYPES_H
#include <vector>
namespace HamDump {
typedef unsigned int uint;
typedef std::vector<int> FDPar;
typedef uint Irrep;

enum Spin{
  alpha = 0,
  beta = 1
};

} //namespace HamDump

#endif