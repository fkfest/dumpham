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
  // occupy first norb orbitals starting from startorb
  Occupation4Irrep(uint startorb, uint norb) { for (uint i = startorb; i < startorb+norb; ++i ) push_back(i);}
  // add to a list of occupied spin orbitals
  void spinocc(std::vector<int>& socc, int ibase = 0) const;
  // number of closed shell orbitals in this symmetry
  uint _nclos;
};
/*! Occupation vectors, i.e., list of orbital indices corresponding to {doubly occupied orbitals, singly occupied orbitals}
 *  The orbital indices are zero based
 */
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
  std::vector<int> spinocc(int ibase = 0) const;
private:
  const PGSym * p_pgs;
};

std::ostream & operator << (std::ostream & o, Occupation const & occ);

/*!
    Orbitals dump with point-group symmetry, C(AO,MO)
*/
class Odump {
public:  
  Odump() : p_pgs(0) {}
  // construct a unity matrix for each symmetry block of orbitals
  // orbitals can be swapped according to the occupation vector
  Odump(const PGSym& pgs, Occupation occs = Occupation());
  // construct from an orbdump file (comma-separated)
  Odump(const PGSym& pgs, const FDPar& ncore, std::string orbdump);
  // construct from an Integ2ab file. Can be with or without core orbitals (will be filled with zeros)
  Odump(const PGSym& pgs, const FDPar& ncore, const Integ2ab& orbs);
  // store orbitals in file orbdump (comma-separated)
  void store(std::string orbdump);
  // return orbital index without core. orb_with_core has to be a valence orbital!
  uint rmcore( uint orb_with_core ) const 
        { assert( _ncoreaccu.size() == p_pgs->nIrreps() );
          uint ncor = _ncoreaccu[p_pgs->irrep(orb_with_core)];
          assert( ncor <= orb_with_core );
          return orb_with_core - ncor; }
  // is the input orbital core?
  bool is_core(uint imo) const 
        { assert( _ncore.size() == p_pgs->nIrreps() );
          Irrep ir = p_pgs->irrep(imo);
          return ( int(imo) < p_pgs->_firstorb4irrep[ir]+_ncore[ir] );}
  // return value c(iao,imo). imo runs over all orbitals including core. 
  // the point group symmetry has to be taken care outside
  double get(uint iao, uint imo) const { return _orbs.get(iao,imo); }
  // return value c(iao,imo). imo runs over all orbitals including core. 
  // resolves the point group symmetry
  double get_with_pgs(uint iao, uint imo) const { return _orbs.get_with_pgs(iao,imo); }
  // guess Basis-occupation vector from orbital coefficients (without core orbitals!)
  // if nclos=nocc=empty- print all orbitals
  Occupation guess_occupation(const FDPar& nclos, const FDPar& nocc) const;
  // transform as c_{\mu p} X_{pq}
  // if frozcore true - trmat doesn't contain core orbitals
  void transform(const Integ2ab& trmat, bool frozcore = true);
private:
  // print value
  void printval(std::ofstream& outputStream, double val, uint j, uint maxlen, bool scientific);
  // Two-electron integrals
  Integ2ab _orbs;
  // point group symmetry
  const PGSym * p_pgs;
  // number of core orbitals in each symmetry (have to be accounted in p_pgs)
  FDPar _ncore, _ncoreaccu;
};




#endif
