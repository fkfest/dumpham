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
  Hdump():_escal(0),_norb(0),_nelec(0),_ms2(0){};
  // construct from FCIdump
  Hdump(std::string fcidump, bool verbose = true);
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
  // in spatial orbitals, PG symmetry is handled outside
  double oneel_spa(uint p, uint q) const {  
          return static_cast<Integ2*>(_oneel[aa].get())->get(p,q);}
  // in spatial orbitals, PG symmetry is handled outside
  inline double twoel_spa(uint p, uint q, uint r, uint s) const { 
          return static_cast<Integ4*>(_twoel[aaaa].get())->get(p,q,r,s);}
  // in spin orbitals, PG symmetry is handled outside
  inline double oneel_spi(uint p, uint q) const; 
  // in spin orbitals, PG symmetry is handled outside
  inline double twoel_spi(uint p, uint q, uint r, uint s) const; 
  // in spatial orbitals
  double oneel_spa_pgs(uint p, uint q) const {  
          return static_cast<Integ2*>(_oneel[aa].get())->get_with_pgs(p,q);}
  // in spatial orbitals
  inline double twoel_spa_pgs(uint p, uint q, uint r, uint s) const { 
          return static_cast<Integ4*>(_twoel[aaaa].get())->get_with_pgs(p,q,r,s);}
  // in spin orbitals
  inline double oneel_spi_pgs(uint p, uint q) const; 
  // in spin orbitals
  inline double twoel_spi_pgs(uint p, uint q, uint r, uint s) const; 
  // get block of integrals defined by start and end indices
  // type of integrals depends on start.size()
  inline void get_block(double * pData, const BlockIndices& start, const BlockIndices& end, bool spinorb = false) const;
  double escal() const {return _escal;}
  Spin spin(uint p) const { return Spin(p%2); }
  
private:
  void store_with_symmetry() const;
  void store_without_symmetry() const;
  // on input: first value (i,j,k,l,value,type)
  template<typename T>
  void readrec(T * pInt, int& i, int& j, int& k, int& l, double& value, FCIdump::integralType& curtype );
  template<typename T>
  void readrec(T * pInt, int& i, int& j, double& value, FCIdump::integralType& curtype );
  void storerec_sym(const Integ4 * pInt) const;
  void storerec_sym(const Integ4ab * pInt) const;
  void storerec_sym(const Integ2 * pInt) const;
  void storerec_nosym(const Integ4 * pInt) const;
  void storerec_nosym(const Integ4ab * pInt) const;
  void storerec_nosym(const Integ2 * pInt) const;
  void check_addressing_integrals() const;
  // check input file for the number of orbitals in each symmetry
  void check_input_norbs(FDPar& orb, const std::string& kind, bool verbose) const;
  // Two-electron integrals (one set for cs rhf, otherwise aa, bb, and ab)
  std::vector< std::unique_ptr<BaseTensors> > _twoel;
  // One-electron integrals (one set for cs rhf, otherwise alpha and beta)
  std::vector< std::unique_ptr<BaseTensors> > _oneel;
  // Scalar
  double _escal;
  FCIdump _dump;
  uint _norb; 
  uint _nelec;
  uint _ms2;
  // wf symmetry (zero based) 
  uint _sym;
  // point group symmetry
  PGSym _pgs, _pgs_wcore;
  FDPar _occ, _clos, _core;
  // dump file in uhf orbitals
  bool _uhf;
  // doesn't have the bra-ket symmetry
  bool _simtra;
};

inline double Hdump::oneel_spi(uint p, uint q) const {
  Spin sp = spin(p), sq = spin(q);
  if( sp != sq) return 0;
  uint i = p/2, j = q/2;
  if(!_uhf || sp == alpha)
     return static_cast<Integ2*>(_oneel[aa].get())->get(i,j);
  return static_cast<Integ2*>(_oneel[bb].get())->get(i,j);
}

inline double Hdump::twoel_spi(uint p, uint q, uint r, uint s) const {
  Spin sp = spin(p), sq = spin(q), sr = spin(r), ss = spin(s);
  if(sp != sq || sr != ss) return 0;
  uint i = p/2, j = q/2, k = r/2, l = s/2;
  if(!_uhf || ( sp == alpha && sr == alpha )) 
    return static_cast<Integ4*>(_twoel[aaaa].get())->get(i,j,k,l);
  if( sp == alpha ) 
    return static_cast<Integ4*>(_twoel[aabb].get())->get(i,j,k,l);
  return static_cast<Integ4*>(_twoel[bbbb].get())->get(i,j,k,l);
}
inline double Hdump::oneel_spi_pgs(uint p, uint q) const {
  Spin sp = spin(p), sq = spin(q);
  if( sp != sq) return 0;
  uint i = p/2, j = q/2;
  if(!_uhf || sp == alpha)
     return static_cast<Integ2*>(_oneel[aa].get())->get_with_pgs(i,j);
  return static_cast<Integ2*>(_oneel[bb].get())->get_with_pgs(i,j);
}

inline double Hdump::twoel_spi_pgs(uint p, uint q, uint r, uint s) const {
  Spin sp = spin(p), sq = spin(q), sr = spin(r), ss = spin(s);
  if(sp != sq || sr != ss) return 0;
  uint i = p/2, j = q/2, k = r/2, l = s/2;
  if(!_uhf || ( sp == alpha && sr == alpha )) 
    return static_cast<Integ4*>(_twoel[aaaa].get())->get_with_pgs(i,j,k,l);
  if( sp == alpha ) 
    return static_cast<Integ4*>(_twoel[aabb].get())->get_with_pgs(i,j,k,l);
  return static_cast<Integ4*>(_twoel[bbbb].get())->get_with_pgs(i,j,k,l);
}

inline void Hdump::get_block(double* pData, const BlockIndices& start, const BlockIndices& end, bool spinorb) const
{
  assert(start.size() == end.size());
  double * p_Data = pData;
  if (start.size() == 2) {
    // h_pq
    double (Hdump::*getHDElement)(uint,uint) const;
    if (spinorb)
      getHDElement = &Hdump::oneel_spi;
    else
      getHDElement = &Hdump::oneel_spa;
    for ( uint64_t q = start[1]; q < end[1]; ++q ) {
      for ( uint64_t p = start[0]; p < end[0]; ++p ) {
        *p_Data = (this->*getHDElement)(p,q);
        ++p_Data;
      }
    }
  } else if (start.size() == 4) {
    // (pq|rs) as <pr | qs>
    double (Hdump::*getHDElement)(uint,uint,uint,uint) const;
    if (spinorb)
      getHDElement = &Hdump::twoel_spi;
    else
      getHDElement = &Hdump::twoel_spa;
    for ( uint64_t s = start[3]; s < end[3]; ++s ) {
      for ( uint64_t q = start[2]; q < end[2]; ++q ) {
        for ( uint64_t r = start[1]; r < end[1]; ++r ) {
          for ( uint64_t p = start[0]; p < end[0]; ++p ) {
            *p_Data = (this->*getHDElement)(p,q,r,s);
            ++p_Data;
          }
        }
      }
    }
  } else {
    error("Number of indices is neither 2 nor 4","Hdump::get_block");
  }
}


} //namespace HamDump

#endif
