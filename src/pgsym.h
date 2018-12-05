#ifndef PGSym_H
#define PGSym_H
#include <string>
#include <vector>
#include "globals.h"
#include "utilities.h"

/*!
 *  Point group symmetry, Irrep is 0 based. 
*/
class PGSym {
public:
  typedef std::vector<Irrep> IrrepVec;
  PGSym() : _nIrreps(0){};
  // construct from orbsym string
  PGSym(const FDPar& orbsym) {
    if ( orbsym.size() == 0 ) 
      error("No orbital symmetries given!");
    _nIrreps = orbsym.back();
    _norb4irrep.resize(_nIrreps,0);
    _firstorb4irrep.resize(_nIrreps,0);
    uint orb = 0;
    _foreach_cauto(FDPar,os,orbsym) {
      if ( *os < 1 ) 
        error("Orbital symmetry below 1!");
      _irrep4orb.push_back(*os-1);
      ++_norb4irrep[*os-1];
      ++orb;
      if ( _firstorb4irrep[*os-1] == 0 ) _firstorb4irrep[*os-1] = orb;
    }
    // make _firstorb4irrep zero based
    _foreach_auto(FDPar,fo,_firstorb4irrep)
      --(*fo);
  };
  // irrep of orbital orb. Orbitals are 0 based
  Irrep irrep( uint orb ) const { return _irrep4orb[orb]; }
  // orbitals are 0 based
  Irrep totIrrep( uint orb1, uint orb2 ) const { return product(_irrep4orb[orb1],_irrep4orb[orb2]); }
//   Irrep totIrrep( uint orb1, uint orb2, uint orb3 ) { return product(_irrep4orb[orb1], totIrrep(orb2,orb3)); }
  Irrep totIrrep( uint orb1, uint orb2, uint orb3, uint orb4 ) const 
                { return product(totIrrep(orb1,orb2), totIrrep(orb3,orb4)); }
  // product of two irreps
  Irrep product(Irrep i, Irrep j) const { return (i^j);}
  // return the original orbsym for fcidump
  FDPar orbsym() const { FDPar osym; _foreach_cauto(IrrepVec,ir,_irrep4orb) osym.push_back(*ir+1); return osym;  }
  // number of orbitals in each irrep
  FDPar norbs_in_irreps() const { return _norb4irrep; }
  // number of orbitals in irrep
  uint norbs( Irrep ir ) const { assert(ir < _norb4irrep.size()); return _norb4irrep[ir]; }
  // total number of orbitals
  uint ntotorbs() const { assert(_nIrreps > 0); return _norb4irrep[_nIrreps-1] + _firstorb4irrep[_nIrreps-1]; }
  uint nIrreps() const { return _nIrreps; }
  
  uint _nIrreps;
  // Irreps of each orbital, Irrep is 0 based!
  IrrepVec _irrep4orb;
  // sizes for each irrep, Irrep is 0 based!
  FDPar _norb4irrep;
  // start of each irrep, Irrep and orbs are 0 based!
  FDPar _firstorb4irrep;
};

#endif
