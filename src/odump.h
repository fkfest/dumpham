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
  // return a list of occupied spin orbitals
  std::vector<int> spinocc(uint nclos) const;
};

std::ostream & operator << (std::ostream & o, Occupation const & occ);

/*!
    Orbitals dump without point-group symmetry
*/
class Odump {
public:  
  Odump() : _nbas(0),_norb(0) {}
  // construct a unity matrix for norb orbitals
  // orbitals can be swapped according to the occupation vector
  Odump(uint norb, const Occupation & occs = Occupation());
  // construct from an orbdump file (comma-separated)
  Odump(std::string orbdump, uint norb = 0);
  // zero orbitals
  void zero() { _orbs.assign(_nbas*_norb,0.0);}
  // element (i,j) (i and j zero based)
  double & operator()(const uint i, const uint j) {
      assert(_orbs.size() == _nbas*_norb);
      assert(i < _nbas && j < _norb);
      return _orbs[i+j*_nbas];
    }
  const double & operator()(const uint i, const uint j) const {
      assert(_orbs.size() == _nbas*_norb);
      assert(i < _nbas && j < _norb);
      return _orbs[i+j*_nbas];
    }
  // store orbitals in file orbdump (comma-separated)
  void store(std::string orbdump);
  // guess Basis-occupation vector from orbital coefficients
  // if nclos=nopen=0 - print all orbitals
  Occupation guess_occupation(uint nclos = 0, uint nopen = 0) const;
private:
  // Two-electron integrals
  IntegralsD _orbs;
  // number of basis and molecular orbitals
  uint _nbas, _norb;
};




#endif
