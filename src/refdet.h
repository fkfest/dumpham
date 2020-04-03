#ifndef RefDet_H
#define RefDet_H

#include <string>
#include <vector>
#include <array>
#include <algorithm>
#include <numeric>
#include <iostream>
#ifdef MOLPRO
#include "hdtypes.h"
#include "hdcommon.h"
#else
#include "globals.h"
#include "utilities.h"
#endif
#include "pgsym.h"

namespace HamDump {

// internal orbital order. 
// e.g. from sorting ORBSYM symmetry-wise and store the order
struct OrbOrder : public std::vector<std::size_t> {
  OrbOrder() : std::vector<std::size_t>() {}
  // order from list of orbitals
  OrbOrder(const FDPar& ORBSYM);
  // put occupied orbitals first for each symmetry
  OrbOrder(const std::vector<uint>& occorb, const uint* nocc, const PGSym& pgs);
  // unity order
  OrbOrder(uint norbs) : std::vector<std::size_t>(norbs) { std::iota(begin(),end(),0); }
  void reorder(FDPar& orbsym) {
    if (!reordered) return;
    assert(orbsym.size() == size());
    FDPar oldorbsym(orbsym);
    for ( uint i = 0; i < size(); ++i ) orbsym[(*this)[i]] = oldorbsym[i];
  }
  bool reordered = false;
};

std::ostream & operator << (std::ostream& o, const OrbOrder& oo);

// reference determinant information
// also contains data from the embedding (e.g., number of core orbitals)
// data should be modified by member functions only
struct RefDet {
  RefDet(){}
  RefDet(const PGSym& pgs);
  // set from CLOSED, OCC and CORE (note that OCC includes CLOSED, which includes CORE!)
  // if wcore = false: core is not included in OCC and CLOSED specifications
  RefDet(const PGSym& pgs, const FDPar& occ_, const FDPar& clos_, 
         const FDPar& core_ = FDPar(), bool wcore = false);
  RefDet(const PGSym& pgs, const std::vector<uint>& occorba, const uint* nocca, 
         const std::vector<uint>& occorbb, const uint* noccb);
  bool operator==(const RefDet& rd) const;
  bool operator!=(const RefDet& rd) const {return !(*this == rd);}
  
  // number of closed shell orbitals in each irrep including core
  FDPar nclos_wcore() const;
  // number of occupied orbitals in each irrep including core
  FDPar nocc_wcore() const;
  // number of core orbitals in each symmetry
  FDPar ncore() const { if ( core.size() > 0 ) return core;
                        else return FDPar(clos.size(),0); }
  // set core from an 8-integers array. Check consistancy if core already set
  template<typename TINT> 
  void set_ncore( TINT* pNCore );
  void set_ncore( const FDPar& core_ );
  // subtract or add core orbitals
  void add_or_subtract_core(FDPar& orb, const std::string& kind, bool add = true) const;
  // generate refso from refao and refbo and from occa/occb
  void gen_refso();
  
  void sanity_check(uint nelec, uint ms2) const;
  void print(int verbosity = -1) const;
  
  FDPar occ, clos, core;
  // number of occupied orbitals alpha and beta
  std::array<FDPar,2> nocc;
  // reference determinant for spinorbitals: occ + virt
  std::vector<SpinOrb> refso;
  // reference determinant for spatial orbitals (alpha and beta): occ + virt
  std::array<OrbOrder,2> ref;
  const PGSym * p_pgs = 0;
  PGSym pgs_wcore;
};

template <typename TINT>
void RefDet::set_ncore(TINT* pNCore)
{
  if ( core.empty() ) {
    // check if the new core is empty
    if (*pNCore < 0) return;
    int ntotcore = 0;
    for ( uint ir = 0; ir < 8; ++ir ) {
      assert(*(pNCore+ir) >= 0);
      ntotcore += *(pNCore+ir);
    }
    if (ntotcore == 0) return;
    // set core
    core.insert(core.begin(),pNCore,pNCore+8);
    uint nIrrepsCore = p_pgs->guess_nIrreps(core);
    FDPar norb4irs(p_pgs->norbs_in_irreps());
    norb4irs.resize(nIrrepsCore,0);
    core.resize(norb4irs.size(),0);
    xout << "Set core to: ";
    for ( auto io: core) xout << io << ",";
    xout << std::endl;
    for ( Irrep ir = 0; ir < norb4irs.size(); ++ir )
      norb4irs[ir] += core[ir];
    pgs_wcore = PGSym(norb4irs,false);
  } else {
    // check consistency 
    for ( uint ir = 0; ir < core.size(); ++ir ) {
      if ( core[ir] != *(pNCore+ir) ) {
        for ( auto io: core) xout << io << " "; 
        xout << std::endl;
        for ( uint i = 0; i < 8; ++i ) xout << *(pNCore+i) << " "; 
        xout << std::endl;
        error("Core already set and differs!");
      }
    }
  }
}

} //namespace HamDump
#endif
