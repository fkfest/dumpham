#ifndef Odump_H
#define Odump_H
#include <string>
#include <vector>
#include <fstream>
#include <iomanip>
#include "pgsym.h"
#include "integs.h"
#include "globals.h"
#include "utilities.h"
#include "inpline.h"

/*!
 *  Occupation vector for an irrep, i.e., list of orbital indices corresponding to {doubly occupied orbitals, singly occupied orbitals} 
 *  The orbital indices are zero based!
 */
class Occupation4Irrep : public std::vector<int> {
public:
  Occupation4Irrep() : std::vector<int>() {}
  // occupy first norb orbitals
  Occupation4Irrep(uint norb) { for (uint i = 0; i < norb; ++i ) push_back(i);}
  // add to a list of occupied spin orbitals
  void spinocc(std::vector<int>& socc) const;
  // number of closed shell orbitals in this symmetry
  uint _nclos;
};

class Occupation : public std::vector<Occupation4Irrep> {
public:
  Occupation() : std::vector<Occupation4Irrep>(), p_pgs(0) {}
  Occupation(const PGSym& pgs) : std::vector<Occupation4Irrep>(pgs.nIrreps()), p_pgs(&pgs) {}
  // occupy first nocc orbitals in each symmetry
  Occupation(const PGSym& pgs, const FDPar& nclos, const FDPar& nocc);
  // occupation from a list of occupied spin orbitals
  // ibase - base for the indices in occs
  Occupation(const PGSym& pgs, const std::vector<int>& occs, int ibase = 1);
  // return a list of occupied spin orbitals
  std::vector<int> spinocc() const;
private:
  const PGSym * p_pgs;
};

std::ostream & operator << (std::ostream & o, Occupation const & occ);

/*!
    Orbitals dump without point-group symmetry
*/
class Odump {
public:  
  Odump() : p_pgs(0) {}
  // construct a unity matrix for each symmetry block of orbitals
  // orbitals can be swapped according to the occupation vector
  Odump(const PGSym& pgs, Occupation occs = Occupation());
  // construct from an orbdump file (comma-separated)
  Odump(const PGSym& pgs, std::string orbdump);
  // store orbitals in file orbdump (comma-separated)
  void store(std::string orbdump);
  // guess Basis-occupation vector from orbital coefficients
  // if nclos=nocc=empty- print all orbitals
  Occupation guess_occupation(const FDPar& nclos, const FDPar& nocc) const;
private:
  // Two-electron integrals
  Integ2ab _orbs;
  // point group symmetry
  const PGSym * p_pgs;
};




#endif
