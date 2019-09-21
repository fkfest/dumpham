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
struct SpinOrb {
  SpinOrb(){};
  SpinOrb(uint orb_, Spin spin_ ) : orb(orb_), spin(spin_){};
  // spatial orbital
  uint orb = 0;
  // spin
  Spin spin = alpha;
};

} //namespace HamDump

#endif
