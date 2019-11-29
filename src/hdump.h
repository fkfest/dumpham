#ifndef Hdump_H
#define Hdump_H
#include <string>
#include <vector>
#include <memory>

#ifdef MOLPRO
#include "hdtypes.h"
#include "hdcommon.h"
#else
#include "globals.h"
#include "utilities.h"
#endif
#include "FCIdump.h"

#include "pgsym.h"
#include "integs.h"

namespace HamDump {
typedef std::vector<uint64_t> BlockIndices;
/*!
    Hamiltonian dump with point-group symmetry
*/
class Hdump {
public:  
  Hdump(){};
  // construct from FCIdump
  Hdump(std::string fcidump, bool verbose = true);
  // construct from PGSym 
  Hdump(const PGSym& pgs_, uint nelec_ = 0, uint ms2_ = 0, 
        uint sym_ = 0, bool uhf_ = false, bool simtra_ = false, bool dm_ = false) :
        _pgs(pgs_),_norb(pgs_.ntotorbs()),_nelec(nelec_), 
        _ms2(ms2_),_sym(sym_),_uhf(uhf_),_simtra(simtra_), _dm(dm_) { gen_spinorbsref(); };
  // construct as a union of two dumps properties (only _uhf and _simtra can differ!)
  Hdump(const Hdump& hd1, const Hdump& hd2);
  // copy info from hd, changing optionally _uhf and _simtra (-1: false, 0: not changed, 1: true) 
  Hdump(const Hdump& hd, int i_uhf = 0, int i_simtra = 0);
  enum onetype {
    aa = 0,
    bb = 1
  };
  enum twotype {
    aaaa = 0,
    bbbb = 1,
    aabb = 2
  };
  void read_dump();
  void store(std::string fcidump);
  void alloc_ints();
  // import integrals from hd. If add is true: add to the current integrals, otherwise the integrals are allocated 
  void import(const Hdump& hd, bool add = false);
  // add integrals from hd
  void add(const Hdump& hd) { import(hd,true); };
  uint norb() const { return _norb; }
  uint nelec() const { return _nelec; }
  uint ms2() const { return _ms2; }
  uint sym() const { return _sym; }
  // number of closed-shell orbitals
  uint nclostot() const;
  // number of open-shell orbitals
  uint nopentot() const { return _ms2; }
  // number of occupied orbitals
  uint nocctot() const { return nclostot()+nopentot(); }
  // file name (empty if not set)
  std::string fcidump_filename() const { return _dump.fileName(); }
  // point group symmetry
  const PGSym& pgs() const { return _pgs; }
  // point group symmetry with core orbitals
  const PGSym& pgs_wcore() const { return _pgs_wcore; }
  // number of closed shell orbitals in each irrep
  const FDPar& nclos() const { return _clos; }
  // number of occupied orbitals in each irrep
  const FDPar& nocc() const { return _occ; }
  // number of core orbitals in each symmetry (are not counted in nclos() or nocc()!)
  FDPar ncore() const { if ( _core.size() > 0 ) return _core;
                        else return FDPar(_clos.size(),0); }
  // number of closed shell orbitals in each irrep including core
  FDPar nclos_wcore() const;
  // number of occupied orbitals in each irrep including core
  FDPar nocc_wcore() const;
  // set CLOSED, OCC and CORE (note that OCC includes CLOSED, which includes CORE!)
  // if wcore = false: core is not included in OCC and CLOSED specifications
  void set_noccorbs(const FDPar& core, const FDPar& closed, const FDPar& occ, bool wcore = true);
  // set core from an 8-integers array. Check consistancy if core already set
  template<typename TINT> 
  void set_ncore( TINT* pNCore );
  // similarity transformed
  bool simtra() const { return _simtra; }
  // in spatial orbitals, PG symmetry is handled outside
  double oneel_spa(uint p, uint q) const {  
          return (_oneel[aa].get())->get(p,q);}
  // in spatial orbitals, PG symmetry is handled outside
  inline void set_oneel_spa(uint p, uint q, double val) {  
          (_oneel[aa].get())->set(p,q,val);}
  // in spatial orbitals, PG symmetry is handled outside
  inline double twoel_spa(uint p, uint q, uint r, uint s) const { 
          return (_twoel[aaaa].get())->get(p,q,r,s);}
  // in spatial orbitals, PG symmetry is handled outside
  inline void set_twoel_spa(uint p, uint q, uint r, uint s, double val) { 
          (_twoel[aaaa].get())->set(p,q,r,s,val);}
  // in spin orbitals, PG symmetry is handled outside
  inline double oneel_spi(uint p, uint q) const; 
  // in spin orbitals, PG symmetry is handled outside
  inline void set_oneel_spi(uint p, uint q, double val); 
  // in spin orbitals, PG symmetry is handled outside
  inline double twoel_spi(uint p, uint q, uint r, uint s) const; 
  // in spin orbitals, PG symmetry is handled outside
  inline void set_twoel_spi(uint p, uint q, uint r, uint s, double val); 
  // in spatial orbitals
  double oneel_spa_pgs(uint p, uint q) const {  
          return (_oneel[aa].get())->get_with_pgs(p,q);}
  // in spatial orbitals
  inline double twoel_spa_pgs(uint p, uint q, uint r, uint s) const { 
          return (_twoel[aaaa].get())->get_with_pgs(p,q,r,s);}
  // in spin orbitals
  inline double oneel_spi_pgs(uint p, uint q) const; 
  // in spin orbitals
  inline double twoel_spi_pgs(uint p, uint q, uint r, uint s) const; 
  // spin-projection, PG symmetry is handled outside
  inline double spinprojector(uint p, uint q, bool alphaspin) const; 
  // get block of integrals defined by start and end indices
  // type of integrals depends on start.size()
  // index order given by Ham[i] <=> Data[order[i]]
  inline void get_block(double * pData, const BlockIndices& start, const BlockIndices& end,
                        const BlockIndices& order, bool spinorb = false) const;
  // set block of integrals defined by start and end indices
  // type of integrals depends on start.size()
  // index order given by Ham[i] <=> Data[order[i]]
  inline void set_block(double * pData, const BlockIndices& start, const BlockIndices& end, 
                        const BlockIndices& order, bool spinorb = false);
  // get block of spin-projectors, see get_block
  inline void get_block_spinprojector(double * pData, const BlockIndices& start, const BlockIndices& end,
                        const BlockIndices& order, bool alphaspin) const;
  double escal() const {return _escal;}
  void set_escal(double escal) { _escal = escal; }
  Spin spin(uint p) const { return Spin(p%2); }
  // scale the integrals by a number
  void scale(double scal);
  bool allocated() const { return _oneel.size() > 0 || _twoel.size() > 0; };
  void correct_2DM();
  void calc_Fock(const Hdump& DMmats);
  
private:
  // on input: first value (i,j,k,l,value,type)
  template<typename T>
  void readrec(T * pInt, int& i, int& j, int& k, int& l, double& value, FCIdump::integralType& curtype );
  template<typename T>
  void readrec(T * pInt, int& i, int& j, double& value, FCIdump::integralType& curtype );
  // simply skip the current record
  void skiprec(int& i, int& j, int& k, int& l, double& value, FCIdump::integralType& curtype );
  template<typename I2, typename I4aa, typename I4ab>
  void store_with_symmetry( I2 * pI2, I4aa * pI4aa, I4ab * pI4ab ) const; 
  template<typename I2, typename I4aa, typename I4ab>
  void store_without_symmetry( I2 * pI2, I4aa * pI4aa, I4ab * pI4ab ) const; 
  template<typename T>
  void storerec2_sym(const T * pInt) const;
  template<typename T>
  void storerec4_sym(const T * pInt) const;
  template<typename T>
  void storerec2_nosym(const T * pInt) const;
  template<typename T>
  void storerec4_nosym(const T * pInt) const;
  // copy integrals with types SI2, SI4aa, SI4ab from hd to this integrals with types DI2, DI4aa, DI4ab 
  // if add = true: add integrals to the existing ones
  // if sym = true: symmetrize simtra to normal integrals
  template<typename DI2, typename DI4aa, typename DI4ab, typename SI2, typename SI4aa, typename SI4ab>
  void copy_ints( DI2 * pDI2, DI4aa * pDI4aa, DI4ab * pDI4ab,
                  SI2 * pSI2, SI4aa * pSI4aa, SI4ab * pSI4ab, const Hdump& hd, bool add = false, bool sym = false);
  // copy 4-index integrals from pSrc to pDest
  // if add = true: add integrals to the existing ones
  // if sym = true: symmetrize simtra to normal integrals
  template<typename T, typename U>
  void copy_int4( T * pDest, const U * pSrc, bool add = false, bool sym = false );
  // copy 2-index integrals from pSrc to pDest
  // if add = true: add integrals to the existing ones
  // if sym = true: symmetrize simtra to normal integrals
  template<typename T, typename U>
  void copy_int2( T * pDest, const U * pSrc, bool add = false, bool sym = false );
  // add antisymmetrized product of 1RDM pD1 and pD11 to 2RDM pD2
  template<typename T, typename U, typename V>
  void add1RDMto2RDM( T * pD2, const U * pD1, const V * pD11, double fact, bool exchange);
  void check_addressing_integrals() const;
  // check input file for the number of orbitals in each symmetry
  void check_input_norbs(FDPar& orb, const std::string& kind, bool verbose) const;
  // subtract or add core orbitals
  void add_or_subtract_core(FDPar& orb, const std::string& kind, bool add = true) const;
  // generate _spinorbs from _occ/_clos/_core
  void gen_spinorbsref();
  // check for nelec, ms2, norb, occ, clos, core...
  void sanity_check() const;
  // Two-electron integrals (one set for cs rhf, otherwise aa, bb, and ab)
  std::vector< std::unique_ptr<BaseTensors> > _twoel;
  // One-electron integrals (one set for cs rhf, otherwise alpha and beta)
  std::vector< std::unique_ptr<BaseTensors> > _oneel;
  // point group symmetry
  PGSym _pgs, _pgs_wcore;
  // Scalar
  double _escal = 0.0;
  FCIdump _dump;
  uint _norb = 0; 
  uint _nelec = 0;
  uint _ms2 = 0;
  // wf symmetry (zero based) 
  uint _sym = 0;
  FDPar _occ, _clos, _core;
  // dump file in uhf orbitals
  bool _uhf = false;
  // doesn't have the bra-ket symmetry
  bool _simtra = false;
  // density matrix instead of integrals
  bool _dm = false;
  // reference determinant for spinorbitals
  std::vector<SpinOrb> _spinorbs;
};

template <typename TINT>
void Hdump::set_ncore(TINT* pNCore)
{
  if ( _core.empty() ) {
    // check if the new core is empty
    if (*pNCore < 0) return;
    int ntotcore = 0;
    for ( uint ir = 0; ir < 8; ++ir ) {
      assert(*(pNCore+ir) >= 0);
      ntotcore += *(pNCore+ir);
    }
    if (ntotcore == 0) return;
    // set core
    _core.insert(_core.begin(),pNCore,pNCore+8);
    uint nIrrepsCore = _pgs.guess_nIrreps(_core);
    FDPar norb4irs(_pgs.norbs_in_irreps());
    norb4irs.resize(nIrrepsCore,0);
    _core.resize(norb4irs.size(),0);
    xout << "Set core to: ";
    for ( auto io: _core) xout << io << ",";
    xout << std::endl;
    for ( Irrep ir = 0; ir < norb4irs.size(); ++ir )
      norb4irs[ir] += _core[ir];
    _pgs_wcore = PGSym(norb4irs,false);
  } else {
    // check consistency 
    for ( uint ir = 0; ir < _core.size(); ++ir ) {
      if ( _core[ir] != *(pNCore+ir) ) {
        for ( auto io: _core) xout << io << " "; 
        xout << std::endl;
        for ( uint i = 0; i < 8; ++i ) xout << *(pNCore+i) << " "; 
        xout << std::endl;
        error("Core already set and differs!");
      }
    }
  }
}

inline double Hdump::oneel_spi(uint p, uint q) const {
  assert( p < _spinorbs.size() && q < _spinorbs.size() );
  const SpinOrb 
    & pso = _spinorbs[p],
    & qso = _spinorbs[q];    
  if ( pso.spin != qso.spin ) return 0.0;
  if (!_uhf || pso.spin == alpha)
     return (_oneel[aa].get())->get(pso.orb,qso.orb);
  return (_oneel[bb].get())->get(pso.orb,qso.orb);
}

inline double Hdump::twoel_spi(uint p, uint q, uint r, uint s) const {
  assert( p < _spinorbs.size() && q < _spinorbs.size() && r < _spinorbs.size() && s < _spinorbs.size() );
  const SpinOrb 
    & pso = _spinorbs[p],
    & qso = _spinorbs[q],    
    & rso = _spinorbs[r],
    & sso = _spinorbs[s];    
  if (pso.spin == qso.spin && rso.spin == sso.spin) {
    if (!_uhf || ( pso.spin == alpha && rso.spin == alpha )) 
      return (_twoel[aaaa].get())->get(pso.orb,qso.orb,rso.orb,sso.orb);
    if ( pso.spin == alpha ) 
      return (_twoel[aabb].get())->get(pso.orb,qso.orb,rso.orb,sso.orb);
    if ( rso.spin == alpha ) 
      return (_twoel[aabb].get())->get(rso.orb,sso.orb,pso.orb,qso.orb);
    return (_twoel[bbbb].get())->get(pso.orb,qso.orb,rso.orb,sso.orb);
  }
  if (_dm && pso.spin == sso.spin && rso.spin == qso.spin) {
    if ( pso.spin == alpha ) 
      return -(_twoel[aabb].get())->get(pso.orb,sso.orb,rso.orb,qso.orb);
    else 
      return -(_twoel[aabb].get())->get(rso.orb,qso.orb,pso.orb,sso.orb);
  }
  return 0.0;
}

inline void Hdump::set_oneel_spi(uint p, uint q, double val) {
  assert( p < _spinorbs.size() && q < _spinorbs.size() );
  const SpinOrb 
    & pso = _spinorbs[p],
    & qso = _spinorbs[q];    
  if ( pso.spin != qso.spin ) return;
  if (!_uhf || pso.spin == alpha)
    (_oneel[aa].get())->set(pso.orb,qso.orb,val);
  else
    (_oneel[bb].get())->set(pso.orb,qso.orb,val);
}
void Hdump::set_twoel_spi(uint p, uint q, uint r, uint s, double val)
{
  assert( p < _spinorbs.size() && q < _spinorbs.size() && r < _spinorbs.size() && s < _spinorbs.size() );
  const SpinOrb 
    & pso = _spinorbs[p],
    & qso = _spinorbs[q],    
    & rso = _spinorbs[r],
    & sso = _spinorbs[s];    
  if (pso.spin != qso.spin || rso.spin != sso.spin) return;
  if (!_uhf || ( pso.spin == alpha && rso.spin == alpha )) {
    (_twoel[aaaa].get())->set(pso.orb,qso.orb,rso.orb,sso.orb,val);
  } else if( pso.spin == alpha ) {
    (_twoel[aabb].get())->set(pso.orb,qso.orb,rso.orb,sso.orb,val);
  } else if ( rso.spin == alpha ) {
    (_twoel[aabb].get())->set(rso.orb,sso.orb,pso.orb,qso.orb,val);
  } else { 
    (_twoel[bbbb].get())->set(pso.orb,qso.orb,rso.orb,sso.orb,val);
  }
}

inline double Hdump::oneel_spi_pgs(uint p, uint q) const {
  assert( p < _spinorbs.size() && q < _spinorbs.size() );
  const SpinOrb 
    & pso = _spinorbs[p],
    & qso = _spinorbs[q];    
  if ( pso.spin != qso.spin ) return 0.0;
  if (!_uhf || pso.spin == alpha)
     return (_oneel[aa].get())->get_with_pgs(pso.orb,qso.orb);
  return (_oneel[bb].get())->get_with_pgs(pso.orb,qso.orb);
}

inline double Hdump::twoel_spi_pgs(uint p, uint q, uint r, uint s) const {
  assert( p < _spinorbs.size() && q < _spinorbs.size() && r < _spinorbs.size() && s < _spinorbs.size() );
  const SpinOrb 
    & pso = _spinorbs[p],
    & qso = _spinorbs[q],    
    & rso = _spinorbs[r],
    & sso = _spinorbs[s];    
  if (pso.spin == qso.spin && rso.spin == sso.spin) {
    if (!_uhf || ( pso.spin == alpha && rso.spin == alpha )) 
      return (_twoel[aaaa].get())->get_with_pgs(pso.orb,qso.orb,rso.orb,sso.orb);
    if ( pso.spin == alpha ) 
      return (_twoel[aabb].get())->get_with_pgs(pso.orb,qso.orb,rso.orb,sso.orb);
    if ( rso.spin == alpha ) 
      return (_twoel[aabb].get())->get_with_pgs(rso.orb,sso.orb,pso.orb,qso.orb);
    return (_twoel[bbbb].get())->get_with_pgs(pso.orb,qso.orb,rso.orb,sso.orb);
  }
  if (_dm && pso.spin == sso.spin && rso.spin == qso.spin) {
    if ( pso.spin == alpha ) 
      return -(_twoel[aabb].get())->get_with_pgs(pso.orb,sso.orb,rso.orb,qso.orb);
    else 
      return -(_twoel[aabb].get())->get_with_pgs(rso.orb,qso.orb,pso.orb,sso.orb);
  }
  return 0.0;  
}

inline void Hdump::get_block(double* pData, const BlockIndices& start, const BlockIndices& end, 
                             const BlockIndices& order, bool spinorb) const
{
  assert(start.size() == end.size() && start.size() == order.size());
  double * p_Data = pData;
  if (start.size() == 2) {
    // h_pq
    uint64_t indx[2];
    double (Hdump::*getHDElement)(uint,uint) const;
    if (spinorb)
      getHDElement = &Hdump::oneel_spi;
    else
      getHDElement = &Hdump::oneel_spa;
    for ( indx[1] = start[1]; indx[1] < end[1]; ++indx[1] ) {
      for ( indx[0] = start[0]; indx[0] < end[0]; ++indx[0] ) {
        *p_Data = (this->*getHDElement)(indx[order[0]],indx[order[1]]);
        ++p_Data;
      }
    }
  } else if (start.size() == 4) {
    // (pq|rs) 
    uint64_t indx[4];
    double (Hdump::*getHDElement)(uint,uint,uint,uint) const;
    if (spinorb)
      getHDElement = &Hdump::twoel_spi;
    else
      getHDElement = &Hdump::twoel_spa;
    for ( indx[3] = start[3]; indx[3] < end[3]; ++indx[3] ) {
      for ( indx[2] = start[2]; indx[2] < end[2]; ++indx[2] ) {
        for ( indx[1] = start[1]; indx[1] < end[1]; ++indx[1] ) {
          for ( indx[0] = start[0]; indx[0] < end[0]; ++indx[0] ) {
            *p_Data = (this->*getHDElement)(indx[order[0]],indx[order[1]],indx[order[2]],indx[order[3]]);
            ++p_Data;
          }
        }
      }
    }
  } else {
    error("Number of indices is neither 2 nor 4","Hdump::get_block");
  }
}
inline void Hdump::set_block(double* pData, const BlockIndices& start, const BlockIndices& end, 
                             const BlockIndices& order, bool spinorb)
{
  assert(start.size() == end.size() && start.size() == order.size());
  double * p_Data = pData;
  if (start.size() == 2) {
    // h_pq
    uint64_t indx[2];
    void (Hdump::*setHDElement)(uint,uint,double);
    if (spinorb)
      setHDElement = &Hdump::set_oneel_spi;
    else
      setHDElement = &Hdump::set_oneel_spa;
    for ( indx[1] = start[1]; indx[1] < end[1]; ++indx[1] ) {
      for ( indx[0] = start[0]; indx[0] < end[0]; ++indx[0] ) {
        (this->*setHDElement)(indx[order[0]],indx[order[1]],*p_Data);
        ++p_Data;
      }
    }
  } else if (start.size() == 4) {
    // (pq|rs) 
    uint64_t indx[4];
    void (Hdump::*setHDElement)(uint,uint,uint,uint,double);
    if (spinorb)
      setHDElement = &Hdump::set_twoel_spi;
    else
      setHDElement = &Hdump::set_twoel_spa;
    for ( indx[3] = start[3]; indx[3] < end[3]; ++indx[3] ) {
      for ( indx[2] = start[2]; indx[2] < end[2]; ++indx[2] ) {
        for ( indx[1] = start[1]; indx[1] < end[1]; ++indx[1] ) {
          for ( indx[0] = start[0]; indx[0] < end[0]; ++indx[0] ) {
            (this->*setHDElement)(indx[order[0]],indx[order[1]],indx[order[2]],indx[order[3]],*p_Data);
            ++p_Data;
          }
        }
      }
    }
  } else {
    error("Number of indices is neither 2 nor 4","Hdump::set_block");
  }
}

inline double Hdump::spinprojector(uint p, uint q, bool alphaspin) const {
  assert( p < _spinorbs.size() && q < _spinorbs.size() );
  if ( p != q ) return 0.0;
  bool isalpha = (_spinorbs[p].spin == alpha);
  if ( alphaspin == isalpha ) return 1.0;
  return 0.0;
}
inline void Hdump::get_block_spinprojector(double* pData, const BlockIndices& start, const BlockIndices& end, 
                             const BlockIndices& order, bool alphaspin) const
{
  assert(start.size() == end.size() && start.size() == order.size());
  if ( start.size() != 2 )
    error("Number of indices is not 2","Hdump::get_block_spinprojector");
  double * p_Data = pData;
  // P_pq
  uint64_t indx[2];
  for ( indx[1] = start[1]; indx[1] < end[1]; ++indx[1] ) {
    for ( indx[0] = start[0]; indx[0] < end[0]; ++indx[0] ) {
      *p_Data = spinprojector(indx[order[0]],indx[order[1]],alphaspin);
      ++p_Data;
    }
  }
}


} //namespace HamDump

#endif
