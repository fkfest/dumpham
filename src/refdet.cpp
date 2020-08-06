#include "refdet.h"
#include <algorithm>

namespace HamDump {

OrbOrder::OrbOrder(const FDPar& ORBSYM) : std::vector<std::size_t>(ORBSYM.size()) 
{
  std::vector<std::size_t> oldorder(ORBSYM.size());
  std::iota(oldorder.begin(),oldorder.end(),0);
  uint nSwaps = InsertionSort(&ORBSYM[0],&oldorder[0],oldorder.size());
  for ( uint i = 0; i < size(); ++i ) (*this)[oldorder[i]]=i;
  if (nSwaps > 0) {
    reordered = true;
    xout << "The orbitals in the fcidump will be reordered according to the symmetries" << std::endl;
    xout << "ORBSYMORDER="; 
    for ( const auto& s: (*this)) xout << s << ","; 
    xout<<std::endl;
  }
}

OrbOrder::OrbOrder(const std::vector<uint>& occorb, const uint* nocc, const PGSym& pgs)
                  : std::vector<std::size_t>(pgs.ntotorbs())
{
  uint iocc = 0;
  uint iorb = 0;
  for (uint ir = 0; ir < pgs.nIrreps(); ++ir ) {
    uint ioff = pgs.beginorb(ir);
    std::vector<bool> occupied(pgs.norbs(ir),false);
    // occupied orbitals first
    for ( uint i = 0; i < nocc[ir]; ++i, ++iocc, ++iorb ) {
      (*this)[iorb] = occorb[iocc];
      assert ( ioff <= occorb[iocc] );
      occupied[occorb[iocc]-ioff] = true;
    }
    // virtual orbitals
    for ( uint i = 0; i < pgs.norbs(ir); ++i ) {
      if ( !occupied[i] ) {
        (*this)[iorb] = i+ioff;
        ++iorb;
      }
    }
  }
  for ( uint i = 0; i < size() && !reordered; ++i ) {
    reordered = ((*this)[i] != i );       
  }
}

std::ostream & operator << (std::ostream& o, const OrbOrder& oo) {
  std::size_t lastorb = oo.size()+100;
  bool print = false;
  for ( auto orb: oo ) {
    if ( orb != lastorb + 1 ) { 
      if (print) o << lastorb << " ";
      o << orb << " ";
      print = false;
    } else {
      if (!print) o << "- ";
      print = true;
    }
    lastorb = orb;
  }
  if (print ) o << lastorb;
//   for ( auto orb: oo ) o << " " << orb;
  if ( oo.reordered ) 
    o << "(REORDERED)";
  else
    o << "(UNITY)";
  return o;
}

RefDet::RefDet(const PGSym& pgs)
{
  p_pgs = &pgs;
  ref[alpha] = OrbOrder(p_pgs->ntotorbs());
  ref[beta] = ref[alpha];
  gen_refso();
}
RefDet::RefDet(const PGSym& pgs, const FDPar& occ_, const FDPar& clos_, 
         const FDPar& core_, bool wcore)
{
  p_pgs = &pgs;
  if ( (core_.size() > 0 && core_[0] > 0) || core_.size() > 1 ) {
    set_ncore(core_);
  }
  if ( (occ_.size() > 0 && occ_[0] > 0) || occ_.size() > 1 ) {
    occ = occ_;
    uint nIrrepsOcc = pgs.guess_nIrreps(occ);
    occ.resize(nIrrepsOcc,0);
    if ( wcore ) {
      // remove core
      add_or_subtract_core(occ,"OCC",false);
    }
    occ.resize(pgs.nIrreps(),0);
  }
  if ( (clos_.size() > 0 && clos_[0] > 0) || clos_.size() > 1 ) {
    clos = clos_;
    uint nIrrepsOcc = pgs.guess_nIrreps(clos);
    clos.resize(nIrrepsOcc,0);
    if ( wcore ) {
      // remove core
      add_or_subtract_core(clos,"CLOSED",false);
    }
    clos.resize(pgs.nIrreps(),0);
    if ( occ.empty() ) {
      // assume closed-shell case
      occ = clos;
    }
  } else if ( !occ.empty() ) {
    clos.resize(pgs.nIrreps(),0);
  }
  ref[alpha] = OrbOrder(p_pgs->ntotorbs());
  ref[beta] = ref[alpha];
  // all open-shell electrons are alpha here
  nocc[alpha] = occ;
  nocc[beta] = clos;
  
  gen_refso();
  print();
}
RefDet::RefDet(const PGSym& pgs, const std::vector<uint>& occorba, const uint* nocca, 
         const std::vector<uint>& occorbb, const uint* noccb)
{
  p_pgs = &pgs;
  auto ita = occorba.begin();
  auto itb = occorbb.begin();
  occ.resize(pgs.nIrreps());
  clos.resize(pgs.nIrreps());
  nocc[alpha].resize(pgs.nIrreps());
  nocc[beta].resize(pgs.nIrreps());
  for (uint ir = 0; ir < pgs.nIrreps(); ++ir ) {
      assert(ita+nocca[ir]-occorba.begin() <= int(occorba.size()));
      assert(itb+noccb[ir]-occorbb.begin() <= int(occorbb.size()));
      std::vector<uint> occo(nocca[ir]+noccb[ir]);
      std::vector<uint>::iterator last;
      last = std::set_union(ita,ita+nocca[ir],itb,itb+noccb[ir],occo.begin());
      occ[ir] = last - occo.begin();
      clos[ir] = nocca[ir]+noccb[ir]-occ[ir];
      nocc[alpha][ir] = nocca[ir];
      nocc[beta][ir] = noccb[ir];
      ita += nocca[ir];
      itb += noccb[ir];
  }
  assert( ita == occorba.end() && itb == occorbb.end() );
  ref[alpha] = OrbOrder(occorba, nocca, pgs);
  ref[beta] = OrbOrder(occorbb, noccb, pgs);
  
  gen_refso();
  print(1);
}

bool RefDet::operator==(const RefDet& rd) const 
{
  return occ == rd.occ && clos == rd.clos && core == rd.core &&
         nocc[alpha] == rd.nocc[alpha] && nocc[beta] == rd.nocc[beta] && 
         ref[alpha] == rd.ref[alpha] && ref[beta] == rd.ref[beta] &&
         pgs_wcore == rd.pgs_wcore;
}

void RefDet::set_ncore(const FDPar& core_)
{
  // core can have more irreps than the hamdump
  assert(core_.size() <= 8);
  if (!core_.empty()){
    FDPar tmpcore(core_);
    tmpcore.resize(8,0);
    set_ncore(&tmpcore[0]);
  }
}

FDPar RefDet::nclos_wcore() const
{
  FDPar clco(clos); 
  add_or_subtract_core(clco,"CLOSED",true);
  return clco; 
}
FDPar RefDet::nocc_wcore() const
{
  FDPar occo(occ); 
  add_or_subtract_core(occo,"OCC",true);
  return occo; 
}
void RefDet::add_or_subtract_core(FDPar& orb, const std::string& kind, bool add) const 
{
  for ( uint io = 0; io < std::min(orb.size(),core.size()); ++io ) {
    if ( add ) {
      orb[io] += core[io];
    } else if ( orb[io] < core[io] ) {
      xout << io << ") " << kind << ": " << orb[io] << " CORE: " << core[io] << std::endl;
      error("Core orbitals should be included in the "+kind+" specification!");
    } else {
      orb[io] -= core[io];
    }
  }
  for ( uint io = orb.size(); io < core.size(); ++io ) {
    if ( add ) {
      orb.push_back(core[io]);
    } else if ( core[io] != 0 ) {
      error("Core orbitals should be included in the "+kind+" specification!");
    }
  }
}

void RefDet::gen_refso()
{
  // order of spin orbitals
  int isoord; 
  #ifndef MOLPRO
  isoord = Input::iPars["orbs"]["soorder"];
  #else
  // spin-blocked order for molpro
  isoord = 1;
  #endif
  uint norb = p_pgs->ntotorbs(); 
  assert( ref[alpha].size() == norb && ref[beta].size() == norb );
  refso.clear();
  refso.reserve(2*norb);
  if ( nocc[alpha].size() == 0 || isoord == 0 ) {
    // default: alternating alpha/beta
    for ( uint i = 0; i < norb; ++i ) {
      refso.emplace_back(SpinOrb(ref[alpha][i],alpha));
      refso.emplace_back(SpinOrb(ref[beta][i],beta));
    }
  } else {
    // blocked spin order: alpha(occ) + beta(occ) + alpha(virt) + beta(virt)
    assert( nocc[alpha].size() == nocc[beta].size() );
    // orbital offset for the irrep
    uint oorb = 0;
    for ( uint ir = 0; ir < nocc[alpha].size(); ++ir ) {
      // alpha occ
      for ( int i = 0; i < nocc[alpha][ir]; ++i ) {
        refso.emplace_back(SpinOrb(ref[alpha][i+oorb],alpha));
      }
      // beta occ
      for ( int i = 0; i < nocc[beta][ir]; ++i ) {
        refso.emplace_back(SpinOrb(ref[beta][i+oorb],beta));
      }
      int norbs = p_pgs->norbs(ir);
      // alpha virt
      for ( int i = nocc[alpha][ir]; i < norbs; ++i ) {
        refso.emplace_back(SpinOrb(ref[alpha][i+oorb],alpha));
      }
      // beta virt
      for ( int i = nocc[beta][ir]; i < norbs; ++i ) {
        refso.emplace_back(SpinOrb(ref[beta][i+oorb],beta));
      }
      oorb += norbs;
    }
  }
}

void RefDet::sanity_check(uint nelec, uint ms2) const
{
  if ( occ.size() < clos.size() )
    error("Number of values in OCC is smaller than in CLOSED!");
  if ( core.size() > 0 || occ.size() > 0 || clos.size() > 0 ) {
    uint norb_core = 0, norb_closed = 0, norb_occ = 0, na = 0, nb = 0;
    for ( auto io: core ) {
      if ( io < 0 ) error("Negative number of orbitals in CORE!");
      norb_core += io;
    }
    auto iocc = occ.begin();
    for ( auto io: clos ) {
      if ( io < 0 ) error("Negative number of orbitals in CLOSED!");
      if ( *iocc < io ) error("Number of OCC orbitals in a symmetry larger than number of CLOSED!");
      ++iocc;
      norb_closed += io;
    }
    for ( auto io: occ ) {
      if ( io < 0 ) error("Negative number of orbitals in OCC!");
      norb_occ += io;
    }
    for ( auto io: nocc[alpha] ) {
      if ( io < 0 ) error("Negative number of orbitals in OCCA!");
      na += io;
    }
    for ( auto io: nocc[beta] ) {
      if ( io < 0 ) error("Negative number of orbitals in OCCB!");
      nb += io;
    }
    if ( norb_occ > p_pgs->ntotorbs() ) {
      error("Number of OCC orbitals exceeds the total number of orbitals!");
    }
    // number of electrons without core!
    if ( norb_closed + norb_occ != nelec ) {
      error("Mismatch in number of electrons and OCC/CLOSED orbitals!");
    }
    if ( na + nb != nelec ) {
      error("Mismatch in number of electrons and ALPHA/BETA orbitals!");
    }
    if ( na - nb != ms2 ) {
      error("MS2 is not equal #ALPHA-#BETA!");
    }
      
  }
  
}

void RefDet::print(int verbosity) const
{
  if (verbosity >= 0) {
    xout << "OCC: "; for (auto io:occ) xout << io << " "; xout << std::endl;
    xout << "CLOS: "; for (auto io:clos) xout << io << " "; xout << std::endl;
  }
  if (verbosity > 0) {
    xout << "NOCCA: "; for (auto io:nocc[alpha]) xout << io << " "; xout << std::endl;
    xout << "NOCCB: "; for (auto io:nocc[beta]) xout << io << " "; xout << std::endl;
  }
  if (verbosity > 1) {
    xout << "Alpha: " << ref[alpha] << std::endl;
    xout << "Beta: " << ref[beta] << std::endl;
  }
}

} //namespace HamDump 
