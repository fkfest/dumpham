#ifndef Hdump_H
#define Hdump_H
#include <string>
#include <vector>
#include <memory>

#ifdef MOLPRO
#include "hdtypes.h"
#include "hdcommon.h"
#include "global/FCIdump.h"
#else
#include "globals.h"
#include "utilities.h"
#include "FCIdump.h"
#endif

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
        uint sym_ = 0, bool uhf_ = false, bool simtra_ = false) :
        _pgs(pgs_),_norb(pgs_.ntotorbs()),_nelec(nelec_), 
        _ms2(ms2_),_sym(sym_),_uhf(uhf_),_simtra(simtra_) {};
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
  // number of core orbitals in each symmetry (are not count in nclos or nocc!)
  const FDPar& ncore() const { return _core; }
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
  double escal() const {return _escal;}
  void set_escal(double escal) { _escal = escal; }
  Spin spin(uint p) const { return Spin(p%2); }
  // scale the integrals by a number
  void scale(double scal);
  bool allocated() const { return _oneel.size() > 0 || _twoel.size() > 0; };
  
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
  void check_addressing_integrals() const;
  // check input file for the number of orbitals in each symmetry
  void check_input_norbs(FDPar& orb, const std::string& kind, bool verbose) const;
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
};

inline double Hdump::oneel_spi(uint p, uint q) const {
  Spin sp = spin(p), sq = spin(q);
  if ( sp != sq) return 0;
  uint i = p/2, j = q/2;
  if (!_uhf || sp == alpha)
     return static_cast<Integ2*>(_oneel[aa].get())->get(i,j);
  return static_cast<Integ2*>(_oneel[bb].get())->get(i,j);
}

inline double Hdump::twoel_spi(uint p, uint q, uint r, uint s) const {
  Spin sp = spin(p), sq = spin(q), sr = spin(r), ss = spin(s);
  if (sp != sq || sr != ss) return 0;
  uint i = p/2, j = q/2, k = r/2, l = s/2;
  if (!_uhf || ( sp == alpha && sr == alpha )) 
    return (_twoel[aaaa].get())->get(i,j,k,l);
  if ( sp == alpha ) 
    return (_twoel[aabb].get())->get(i,j,k,l);
  return (_twoel[bbbb].get())->get(i,j,k,l);
}

inline void Hdump::set_oneel_spi(uint p, uint q, double val) {
  Spin sp = spin(p), sq = spin(q);
  if ( sp != sq) return;
  uint i = p/2, j = q/2;
  if (!_uhf || sp == alpha)
    (_oneel[aa].get())->set(i,j,val);
  else
    (_oneel[bb].get())->set(i,j,val);
}
void Hdump::set_twoel_spi(uint p, uint q, uint r, uint s, double val)
{
  Spin sp = spin(p), sq = spin(q), sr = spin(r), ss = spin(s);
  if (sp != sq || sr != ss) return;
  uint i = p/2, j = q/2, k = r/2, l = s/2;
  if (!_uhf || ( sp == alpha && sr == alpha )) 
    (_twoel[aaaa].get())->set(i,j,k,l,val);
  else if( sp == alpha ) 
    (_twoel[aabb].get())->set(i,j,k,l,val);
  else 
    (_twoel[bbbb].get())->set(i,j,k,l,val);
}

inline double Hdump::oneel_spi_pgs(uint p, uint q) const {
  Spin sp = spin(p), sq = spin(q);
  if ( sp != sq) return 0;
  uint i = p/2, j = q/2;
  if (!_uhf || sp == alpha)
     return (_oneel[aa].get())->get_with_pgs(i,j);
  return (_oneel[bb].get())->get_with_pgs(i,j);
}

inline double Hdump::twoel_spi_pgs(uint p, uint q, uint r, uint s) const {
  Spin sp = spin(p), sq = spin(q), sr = spin(r), ss = spin(s);
  if (sp != sq || sr != ss) return 0;
  uint i = p/2, j = q/2, k = r/2, l = s/2;
  if (!_uhf || ( sp == alpha && sr == alpha )) 
    return (_twoel[aaaa].get())->get_with_pgs(i,j,k,l);
  if ( sp == alpha ) 
    return (_twoel[aabb].get())->get_with_pgs(i,j,k,l);
  return (_twoel[bbbb].get())->get_with_pgs(i,j,k,l);
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


} //namespace HamDump

#endif
