#ifndef PGSym_H
#define PGSym_H

#include <string>
#include <vector>
#include <algorithm>
#include <numeric>
#ifdef MOLPRO
#include "hdtypes.h"
#include "hdcommon.h"
#else
#include "globals.h"
#include "utilities.h"
#endif

namespace HamDump {
/*!
 *  Point group symmetry, Irrep is 0 based. 
*/
class PGSym {
public:
  typedef std::vector<Irrep> IrrepVec;
  PGSym(){};
  // construct from orbsym string
  // if listoforbs: orbsym is a list of orbitals, otherwise it's number of orbitals in each symmetry
  PGSym(const FDPar& orbsym, bool listoforbs = true) {
    if ( orbsym.size() == 0 ) 
      error("No orbital symmetries given!");
    if ( listoforbs ) {
      uint nIrreps = orbsym.back();
      roundof_nIrreps(nIrreps);
      _norb4irrep.resize(nIrreps,0);
      _firstorb4irrep.resize(nIrreps,0);
      uint orb = 0;
      for ( const auto& os: orbsym ) {
        if ( os < 1 ) 
          error("Orbital symmetry below 1!");
        _irrep4orb.push_back(os-1);
        ++_norb4irrep[os-1];
        ++orb;
        if ( _firstorb4irrep[os-1] == 0 ) _firstorb4irrep[os-1] = orb;
      }
      // make _firstorb4irrep zero based
      orb = 0;
      uint ir = 0;
      for ( auto& fo: _firstorb4irrep ) {
        if ( fo > 0 ) {
          --fo;
        } else {
          fo = orb;
        }
        orb = fo + _norb4irrep[ir];
        ++ir;
      }
    } else {
      _norb4irrep = orbsym;
      for ( Irrep ir = 0; ir < _norb4irrep.size(); ++ir ){
        _firstorb4irrep.push_back(_irrep4orb.size());
        for ( int imo = 0; imo < _norb4irrep[ir]; ++imo ) _irrep4orb.push_back(ir);
      }
    }
  };
  // irrep of orbital orb. Orbitals are 0 based
  Irrep irrep( uint orb ) const { return _irrep4orb[orb]; }
  // orbitals are 0 based
  Irrep totIrrep( uint orb1, uint orb2 ) const { return product(_irrep4orb[orb1],_irrep4orb[orb2]); }
//   Irrep totIrrep( uint orb1, uint orb2, uint orb3 ) { return product(_irrep4orb[orb1], totIrrep(orb2,orb3)); }
  Irrep totIrrep( uint orb1, uint orb2, uint orb3, uint orb4 ) const 
                { return product(totIrrep(orb1,orb2), totIrrep(orb3,orb4)); }
  Irrep totIrrep( uint orb1, uint orb2, uint orb3, uint orb4, uint orb5, uint orb6 ) const
                { return product(product(totIrrep(orb1,orb2), totIrrep(orb3,orb4)),totIrrep(orb5,orb6)); }
  // product of two irreps
  Irrep product(Irrep i, Irrep j) const { assert((i^j)< nIrreps()); return (i^j);}
  // return the original orbsym for fcidump
  FDPar orbsym() const { FDPar osym; for (const auto& ir:_irrep4orb) osym.push_back(ir+1); return osym;  }
  // number of orbitals in each irrep
  FDPar norbs_in_irreps() const { return _norb4irrep; }
  // number of orbitals in irrep
  uint norbs( Irrep ir ) const { assert(ir < _norb4irrep.size()); return _norb4irrep[ir]; }
  // first orbital in irrep
  uint beginorb( Irrep ir ) const { assert(ir < _firstorb4irrep.size()); return _firstorb4irrep[ir]; }
  // next after last orbital in irrep
  uint endorb( Irrep ir ) const { assert(ir < _norb4irrep.size()); return _norb4irrep[ir]+_firstorb4irrep[ir]; }
  // total number of orbitals
  uint ntotorbs() const { assert(nIrreps() > 0); return _norb4irrep.back() + _firstorb4irrep.back(); }
  uint nIrreps() const { return _norb4irrep.size(); }
  void roundof_nIrreps( uint& nIrreps ) const {
      // number of irreps can be only 1,2,4, or 8
      // if it's not, it means that the last irreps don't have any orbitals!
      if ( nIrreps > 4 )
        nIrreps = 8;
      else if ( nIrreps > 2 )
        nIrreps = 4;
  }
  // guess nIrreps for a given norbs string. The result is always between or equal nIrreps() and 8. 
  uint guess_nIrreps(const FDPar& norbs4irreps) const {
        uint nIrr = norbs4irreps.size();
        for ( ; nIrr > 0 && norbs4irreps[nIrr-1] == 0; --nIrr ) {}
        roundof_nIrreps(nIrr);
        return std::max(nIrr,nIrreps());
  }
  bool operator==(const PGSym& pgs) const { 
       // if _irrep4orb is equal than all other arrays have to be equal, too.
       assert(_irrep4orb != pgs._irrep4orb || (_norb4irrep == pgs._norb4irrep && _firstorb4irrep == pgs._firstorb4irrep));
       return (_irrep4orb == pgs._irrep4orb);
  }
  bool operator!=(const PGSym& pgs) const { return !(*this == pgs); } 
  
  // Irreps of each orbital, Irrep is 0 based!
  IrrepVec _irrep4orb;
  // sizes for each irrep, Irrep is 0 based!
  FDPar _norb4irrep;
  // start of each irrep, Irrep and orbs are 0 based!
  FDPar _firstorb4irrep;
};

} //namespace HamDump
#endif
