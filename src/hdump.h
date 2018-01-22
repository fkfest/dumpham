#ifndef Hdump_H
#define Hdump_H
#include <string>
#include <vector>
#include "globals.h"
#include "utilities.h"
#include "FCIdump.h"
/*!
    Hamiltonian dump without point-group symmetry
*/
class Hdump {
public:  
  Hdump():_escal(0){};
  // construct from FCIdump
  Hdump(std::string fcidump);
  
  void store(std::string fcidump);
private:
  // Two-electron integrals
  Integrals _twoel;
  // One-electron integrals
  Integrals _oneel;
  // Scalar
  double _escal;
  FCIdump _dump;
  
};




#endif