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
  Hdump():_escal(0),_norb(0),_nelec(0),_ms2(0){};
  // construct from FCIdump
  Hdump(std::string fcidump);
  
  void store(std::string fcidump);
  uint norb() const { return _norb; }
  uint nelec() const { return _nelec; }
  uint ms2() const { return _ms2; }
  // number of closed-shell orbitals
  uint nclosed() const;
  // number of open-shell orbitals
  uint nopen() const { return _ms2; }
  // number of occupied orbitals
  uint nocc() const { return nclosed()+nopen(); }
  
private:
  // Two-electron integrals
  Integrals _twoel;
  // One-electron integrals
  Integrals _oneel;
  // Scalar
  double _escal;
  FCIdump _dump;
  uint _norb; 
  uint _nelec;
  uint _ms2;
};


#endif