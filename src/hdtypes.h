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

// various orbspaces, don't overlap with Spin type 
enum OrbSpace{
  OS_occupied = 0x02,
  OS_empty =    0x04,
  OS_core =     0x08,
  OS_closed =   0x10,
  OS_open =     0x20,
  OS_active =   0x40
};

struct SpinOrb {
  SpinOrb(){};
  SpinOrb(uint orb_, Spin spin_ ) : orb(orb_), spin(spin_){};
  // spatial orbital
  uint orb = 0;
  // spin
  Spin spin = alpha;
};

namespace Numbers
{
  // "zero"
  static const double verysmall=1.e-16;
  static const double small=1.e-10;
  static const int big=1000;
}


} //namespace HamDump

#endif
