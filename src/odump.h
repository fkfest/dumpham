#ifndef Odump_H
#define Odump_H
#include <string>
#include <vector>
#include <fstream>
#include <iomanip>
#include "globals.h"
#include "utilities.h"
#include "inpline.h"

/*!
 *  Occupation vector, i.e., list of orbital indices corresponding to {doubly occupied orbitals, singly occupied orbitals} 
 *  The orbital indices are zero based!
 */
class Occupation : public std::vector<int> {
public:
  Occupation() : std::vector<int>() {}
  // occupy first norb orbitals
  Occupation(uint norb) { for (uint i = 0; i < norb; ++i ) push_back(i);}
  // occupation from a list of occupied spin orbitals
  // ibase - base for the indices in occs
  Occupation(const std::vector<int>& occs, int ibase = 1);
};

std::ostream & operator << (std::ostream & o, Occupation const & occ);

/*!
    Orbitals dump without point-group symmetry
*/
class Odump {
public:  
  Odump() : _nAO(0),_norb(0) {}
  // construct a unity matrix for norb orbitals
  // orbitals can be swapped according to the occupation vector
  Odump(uint norb, const Occupation & occs = Occupation());
  // construct from an orbdump file (comma-separated)
  Odump(std::string orbdump, uint norb = 0);
  // zero orbitals
  void zero() { _orbs.assign(_nAO*_norb,0.0);}
  // element (i,j) (i and j zero based)
  double & operator()(const uint i, const uint j) {
      assert(_orbs.size() == _nAO*_norb);
      assert(i < _nAO && j < _norb);
      return _orbs[i+j*_nAO];
    }
  const double & operator()(const uint i, const uint j) const {
      assert(_orbs.size() == _nAO*_norb);
      assert(i < _nAO && j < _norb);
      return _orbs[i+j*_nAO];
    }
  // store orbitals in file orbdump (comma-separated)
  void store(std::string orbdump);
  // guess AO-occupation vector from orbital coefficients
  // if nocc=0 - print all orbitals
  Occupation guess_occupation(uint nocc = 0) const;
private:
  // Two-electron integrals
  Integrals _orbs;
  // number of AO and MO orbitals
  uint _nAO, _norb;
};




#endif
