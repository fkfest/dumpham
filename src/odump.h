#ifndef Odump_H
#define Odump_H
#include <string>
#include <vector>
#include <fstream>
#include <iomanip>
#include "globals.h"
#include "utilities.h"
/*!
    Orbitals dump without point-group symmetry
*/
class Odump {
public:  
  Odump() : _nAO(0),_norb(0) {};
  // construct unity matrix for norb orbitals
  Odump(uint norb);
  // zero orbitals
  void zero() { _orbs.assign(_nAO*_norb,0.0);};
  // element (i,j) (i and j zero based)
  double & operator()(const uint i, const uint j) {
      assert(_orbs.size() == _nAO*_norb);
      assert(i < _nAO && j < _norb);
      return _orbs[i+j*_nAO];
    }
  // store orbitals in file orbdump
  void store(std::string orbdump);
private:
  // Two-electron integrals
  Integrals _orbs;
  // number of AO and MO orbitals
  uint _nAO, _norb;
};




#endif