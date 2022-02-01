#ifndef Hdump_H
#define Hdump_H
#include <string>
#include <vector>
#include <memory>

#ifdef MOLPRO
#include "hdtypes.h"
#include "hdcommon.h"
#include "global/fmt.h"
#else
#include "globals.h"
#include "utilities.h"
#endif
#include "FCIdump.h"

#include "pgsym.h"
#include "integs.h"
#include "refdet.h"
#include "periodic.h"

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
        _ms2(ms2_),_sym(sym_),_uhf(uhf_),_simtra(simtra_), _dm(dm_) { _rd = RefDet(_pgs); }
  Hdump(std::vector<uint> dims, int charge, int ms2, std::vector<int> pbcs, double Upar, const std::vector<double>& tpar);
  Hdump(const Periodic& pers, int charge, int ms2, double Upar, const std::vector<double>& tpar);
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

  static constexpr Spin spin4el[2][3] = { { alpha, beta, alpha }, { alpha, beta, beta} };
  // set occupation from alpha and beta orbital sets
  // adapt _clos and _occ accordingly
  void set_occupationAB(const std::vector<uint>& occorba, const uint* nocca,
                        const std::vector<uint>& occorbb, const uint* noccb);
  void read_dump();
  // for 3 body integrals
  void read_3body_dump();
  void read_3body_dump_nosym();
  void writeIntegral_3body(int i, int j, int k, int l, int m, int n, double value, std::ofstream& outputStream) const;
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
  bool uhf() const { return _uhf; }
  // file name (empty if not set)
  std::string fcidump_filename() const { return _dump.fileName(); }
  // point group symmetry
  const PGSym& pgs() const { return _pgs; }
  // point group symmetry with core orbitals
  const PGSym& pgs_wcore() const { return _rd.pgs_wcore; }
  // internal order of orbital indices (idx_in_hdump = osord[idx_in_fcidump])
  const OrbOrder& orborder() const { return _osord; }
  // number of closed shell orbitals in each irrep
  const FDPar& nclos() const { return _rd.clos; }
  // number of occupied orbitals in each irrep
  const FDPar& nocc() const { return _rd.occ; }
  // number of core orbitals in each symmetry (are not counted in nclos() or nocc()!)
  FDPar ncore() const { return _rd.ncore(); }
  // number of closed shell orbitals in each irrep including core
  FDPar nclos_wcore() const { return _rd.nclos_wcore(); }
  // number of occupied orbitals in each irrep including core
  FDPar nocc_wcore() const { return _rd.nocc_wcore(); }
  // set CLOSED, OCC and CORE (note that OCC includes CLOSED, which includes CORE!)
  // if wcore = false: core is not included in OCC and CLOSED specifications
  void set_noccorbs(const FDPar& core, const FDPar& closed, const FDPar& occ, bool wcore = true)
       { _rd = RefDet(_pgs,occ,closed,core,wcore); }
  // set core from an 8-integers array. Check consistancy if core already set
  template<typename TINT>
  void set_ncore( TINT* pNCore );
  // similarity transformed
  bool simtra() const { return _simtra; }
  bool threebody_nosym() const { return _3body_nosym; }
  // in spatial orbitals (hdump order), PG symmetry is handled outside
  double oneel_spa(uint p, uint q, onetype xx = aa) const {
          if (!( (_uhf?xx:aa) < _oneel.size() )) {
            xout << (_uhf?xx:aa) << " xx: " << int(xx) << " aa " << int(aa) << " onelsize: " << _oneel.size() << std::endl;
          }
          assert( (_uhf?xx:aa) < _oneel.size() );
          const OrbOrder& oo = _rd.ref[spin4el[0][xx]];
          return (_oneel[_uhf?xx:aa].get())->get(oo[p],oo[q]);}
  // in spatial orbitals (hdump order), PG symmetry is handled outside
  inline void set_oneel_spa(uint p, uint q, double val, onetype xx = aa) {
          assert( xx < _oneel.size() );
          const OrbOrder& oo = _rd.ref[spin4el[0][xx]];
          (_oneel[xx].get())->set(oo[p],oo[q],val);}
  // in spatial orbitals (hdump order), PG symmetry is handled outside
  inline double twoel_spa(uint p, uint q, uint r, uint s, twotype xxxx = aaaa) const {
          assert( (_uhf?xxxx:aaaa) < _twoel.size() );
          const OrbOrder& o1 = _rd.ref[spin4el[0][xxxx]];
          const OrbOrder& o2 = _rd.ref[spin4el[1][xxxx]];
          return (_twoel[_uhf?xxxx:aaaa].get())->get(o1[p],o1[q],o2[r],o2[s]);}
  // in spatial orbitals (hdump order), PG symmetry is handled outside
  inline void set_twoel_spa(uint p, uint q, uint r, uint s, double val, twotype xxxx = aaaa) {
          assert( xxxx < _twoel.size() );
          const OrbOrder& o1 = _rd.ref[spin4el[0][xxxx]];
          const OrbOrder& o2 = _rd.ref[spin4el[1][xxxx]];
          (_twoel[xxxx].get())->set(o1[p],o1[q],o2[r],o2[s],val);}
  // in spatial orbitals (hdump order), PG symmetry is handled outside
  inline double threeel_spa(uint p, uint q, uint r, uint s, uint t, uint u, twotype xxxx = aaaa) const {
          assert(!_uhf);
          const OrbOrder& o1 = _rd.ref[alpha];
          const OrbOrder& o2 = _rd.ref[alpha];
          const OrbOrder& o3 = _rd.ref[alpha];
          return (_threeel[_uhf?xxxx:aaaa].get())->get(o1[p],o1[q],o2[r],o2[s],o3[t],o3[u]);}
  // in spin orbitals, PG symmetry is handled outside
  inline double oneel_spi(uint p, uint q) const;
  // in spin orbitals, PG symmetry is handled outside
  inline void set_oneel_spi(uint p, uint q, double val);
  // in spin orbitals, PG symmetry is handled outside
  inline double twoel_spi(uint p, uint q, uint r, uint s) const;
  // in spin orbitals, PG symmetry is handled outside
  inline void set_twoel_spi(uint p, uint q, uint r, uint s, double val);
  // in spin orbitals, PG symmetry is handled outside
  inline double threeel_spi(uint p, uint q, uint r, uint s, uint t, uint u) const;
  // in spatial orbitals (hdump order)
  double oneel_spa_pgs(uint p, uint q, onetype xx = aa) const {
          assert( (_uhf?xx:aa) < _oneel.size() );
          const OrbOrder& oo = _rd.ref[spin4el[0][xx]];
          return (_oneel[_uhf?xx:aa].get())->get_with_pgs(oo[p],oo[q]);}
  // in spatial orbitals (hdump order)
  inline double twoel_spa_pgs(uint p, uint q, uint r, uint s, twotype xxxx = aaaa) const {
          assert( (_uhf?xxxx:aaaa) < _twoel.size() );
          const OrbOrder& o1 = _rd.ref[spin4el[0][xxxx]];
          const OrbOrder& o2 = _rd.ref[spin4el[1][xxxx]];
          return (_twoel[_uhf?xxxx:aaaa].get())->get_with_pgs(o1[p],o1[q],o2[r],o2[s]);}
  // in spin orbitals
  inline double oneel_spi_pgs(uint p, uint q) const;
  // in spin orbitals
  inline double twoel_spi_pgs(uint p, uint q, uint r, uint s) const;
  // spin-projection, PG symmetry is handled outside
  inline double spinprojector(uint p, uint q, Spin spintype) const;
  // get block of integrals defined by start and end indices
  // type of integrals depends on start.size()
  // index order given by Ham[i] <=> Data[order[i]]
  // spintype from onetype or twotype
  inline void get_block(double * pData, const BlockIndices& start, const BlockIndices& end,
                        const BlockIndices& order, bool spinorb = false, uint spintype = 0) const;
  // set block of integrals defined by start and end indices
  // type of integrals depends on start.size()
  // index order given by Ham[i] <=> Data[order[i]]
  // spintype from onetype or twotype
  inline void set_block(double * pData, const BlockIndices& start, const BlockIndices& end,
                        const BlockIndices& order, bool spinorb = false, uint spintype = 0);
  // get block of spin-projectors, see get_block
  inline void get_block_spinprojector(double * pData, const BlockIndices& start, const BlockIndices& end,
                        const BlockIndices& order, uint spintype) const;
  double escal() const {return _escal;}
  void set_escal(double escal) { _escal = escal; }
  Spin spin(uint p) const { return Spin(p%2); }
  // scale the integrals by a number
  void scale(double scal);
  bool allocated() const { return _oneel.size() > 0 || _twoel.size() > 0; };
  void correct_2DM();
  void calc_Fock(const Hdump& DMmats);
  // brings 2RDM with only particle-exchange symmetry in compressed Molcas format
  void molcas_2DM(const std::string& filename, bool plus = true);
  //calculates total spin expectation value
  void calc_Spin2(const Integ2ab& Overlap, bool alt = false);
  //store 1RDM in ASCII format
  void store1RDM(std::string filepname, std::string filename, bool uhf);
  //gen integrals for Hubbard model
  void gen_hubbard(const std::vector<uint>& dims, const std::vector<int>& pbcs, double Upar, const std::vector<double>& tpar);
  void gen_hubbard(const Periodic& pers, double Upar, const std::vector<double>& tpar);
  // add scal*S^2 to the Hamiltonian
  void addS2(double scal = 1.0);

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
  template<typename I6>
  void store_with_symmetry( I6 * pI6 ) const;
  template<typename I6>
  void store_without_symmetry( I6 * pI6 ) const;
  template<typename I2, typename I4aa, typename I4ab>
  void store_without_symmetry( I2 * pI2, I4aa * pI4aa, I4ab * pI4ab ) const;
  template<typename T>
  void storerec2_sym(const T * pInt) const;
  template<typename T>
  void storerec4_sym(const T * pInt) const;
  template<typename T>
  void storerec6_sym(const T * pInt) const;
  template<typename T>
  void storerec6_nosym(const T * pInt) const;
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
  template<typename T, typename U>
  void copy_int6( T * pDest, const U * pSrc, bool add );
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
  // scale integrals
  template<typename SI2, typename SI4aa, typename SI4ab>
  void scale_ints( SI2 * pSI2, SI4aa * pSI4aa, SI4ab * pSI4ab, double scal);
  // scale 4-index integrals
  template<typename T>
  void scale_int4( T * pInt, double scal);
  // scale 2-index integrals
  template<typename T>
  void scale_int2( T * pInt, double scal);
  void check_addressing_integrals() const;
  // check input file for the number of orbitals in each symmetry
  void check_input_norbs(FDPar& orb, const std::string& kind, bool verbose) const;
  // calculates spin-summed 2RDM value from spatial orbital indices
  double spinsum_2DM(uint p, uint q, uint r, uint s);
  // calculates spin-summed 1RDM value from spatial orbital indices
  double spinsum_1DM(uint p, uint q, Irrep ir);
  // check for nelec, ms2, norb, occ, clos, core...
  void sanity_check() const;
  // Two-electron integrals (one set for rhf, otherwise aa, bb, and ab)
  std::vector< std::unique_ptr<BaseTensors> > _twoel;
  // One-electron integrals (one set for rhf, otherwise alpha and beta)
  std::vector< std::unique_ptr<BaseTensors> > _oneel;
  // Three-electron integrals (one set for rhf, otherwise aaa, aab, abb, bbb)
  std::vector< std::unique_ptr<BaseTensors> > _threeel;
  std::vector< std::unique_ptr<BaseTensors_nosym> > _threeel_nosym;
  // point group symmetry
  PGSym _pgs;
  // Scalar
  double _escal = 0.0;
  FCIdump _dump;
  uint _norb = 0;
  uint _nelec = 0;
  uint _ms2 = 0;
  // wf symmetry (zero based)
  uint _sym = 0;
  // order of orbitals to bring them to symmetry-sorted order
  // used for reading in only
  OrbOrder _osord;
  // dump file in uhf orbitals
  bool _uhf = false;
  // doesn't have the bra-ket symmetry
  bool _simtra = false;
  // three-body integrals
  bool _3body = false;
  // three-body integrals with no symmetry at all
  bool _3body_nosym = false;
  std::string _3body_file;
  // density matrix instead of integrals
  bool _dm = false;
  // reference determinant
  RefDet _rd;
};

template <typename TINT>
void Hdump::set_ncore(TINT* pNCore)
{
  _rd.set_ncore(pNCore);
}

inline double Hdump::oneel_spi(uint p, uint q) const {
  assert( p < _rd.refso.size() && q < _rd.refso.size() );
  const SpinOrb
    & pso = _rd.refso[p],
    & qso = _rd.refso[q];
  if ( pso.spin != qso.spin ) return 0.0;
  if (!_uhf || pso.spin == alpha)
     return (_oneel[aa].get())->get(pso.orb,qso.orb);
  return (_oneel[bb].get())->get(pso.orb,qso.orb);
}

inline double Hdump::twoel_spi(uint p, uint q, uint r, uint s) const {
  assert( p < _rd.refso.size() && q < _rd.refso.size() && r < _rd.refso.size() && s < _rd.refso.size() );
  const SpinOrb
    & pso = _rd.refso[p],
    & qso = _rd.refso[q],
    & rso = _rd.refso[r],
    & sso = _rd.refso[s];
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

inline double Hdump::threeel_spi(uint p, uint q, uint r, uint s, uint t, uint u) const {
  assert( p < _rd.refso.size() && q < _rd.refso.size() && r < _rd.refso.size() && s < _rd.refso.size() && t < _rd.refso.size() && u < _rd.refso.size() );
  assert( !_dm && !_uhf );
  const SpinOrb
    & pso = _rd.refso[p],
    & qso = _rd.refso[q],
    & rso = _rd.refso[r],
    & sso = _rd.refso[s],
    & tso = _rd.refso[t],
    & uso = _rd.refso[u];
  if (pso.spin == qso.spin && rso.spin == sso.spin && tso.spin == uso.spin) {
    return (_threeel[aa].get())->get(pso.orb,qso.orb,rso.orb,sso.orb,tso.orb,uso.orb);
  }
  return 0.0;
}

inline void Hdump::set_oneel_spi(uint p, uint q, double val) {
  assert( p < _rd.refso.size() && q < _rd.refso.size() );
  const SpinOrb
    & pso = _rd.refso[p],
    & qso = _rd.refso[q];
  if ( pso.spin != qso.spin ) return;
  if (!_uhf || pso.spin == alpha)
    (_oneel[aa].get())->set(pso.orb,qso.orb,val);
  else
    (_oneel[bb].get())->set(pso.orb,qso.orb,val);
}
void Hdump::set_twoel_spi(uint p, uint q, uint r, uint s, double val)
{
  assert( p < _rd.refso.size() && q < _rd.refso.size() && r < _rd.refso.size() && s < _rd.refso.size() );
  const SpinOrb
    & pso = _rd.refso[p],
    & qso = _rd.refso[q],
    & rso = _rd.refso[r],
    & sso = _rd.refso[s];
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
  assert( p < _rd.refso.size() && q < _rd.refso.size() );
  const SpinOrb
    & pso = _rd.refso[p],
    & qso = _rd.refso[q];
  if ( pso.spin != qso.spin ) return 0.0;
  if (!_uhf || pso.spin == alpha)
     return (_oneel[aa].get())->get_with_pgs(pso.orb,qso.orb);
  return (_oneel[bb].get())->get_with_pgs(pso.orb,qso.orb);
}

inline double Hdump::twoel_spi_pgs(uint p, uint q, uint r, uint s) const {
  assert( p < _rd.refso.size() && q < _rd.refso.size() && r < _rd.refso.size() && s < _rd.refso.size() );
  const SpinOrb
    & pso = _rd.refso[p],
    & qso = _rd.refso[q],
    & rso = _rd.refso[r],
    & sso = _rd.refso[s];
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
                             const BlockIndices& order, bool spinorb, uint spintype) const
{
  assert(start.size() == end.size() && start.size() == order.size());
  double * p_Data = pData;
  if (start.size() == 2) {
    // h_pq
    uint64_t indx[2];
    for ( indx[1] = start[1]; indx[1] < end[1]; ++indx[1] ) {
      for ( indx[0] = start[0]; indx[0] < end[0]; ++indx[0] ) {
        if (spinorb)
          *p_Data = oneel_spi(indx[order[0]],indx[order[1]]);
        else
          *p_Data = oneel_spa(indx[order[0]],indx[order[1]],onetype(spintype));
        ++p_Data;
      }
    }
  } else if (start.size() == 4) {
    // (pq|rs)
    uint64_t indx[4];
    for ( indx[3] = start[3]; indx[3] < end[3]; ++indx[3] ) {
      for ( indx[2] = start[2]; indx[2] < end[2]; ++indx[2] ) {
        for ( indx[1] = start[1]; indx[1] < end[1]; ++indx[1] ) {
          for ( indx[0] = start[0]; indx[0] < end[0]; ++indx[0] ) {
            if (spinorb)
              *p_Data = twoel_spi(indx[order[0]],indx[order[1]],indx[order[2]],indx[order[3]]);
            else
              *p_Data = twoel_spa(indx[order[0]],indx[order[1]],indx[order[2]],indx[order[3]],twotype(spintype));
            ++p_Data;
          }
        }
      }
    }
  } else if (start.size() == 6) {
    // (pq|rs|tu)
    uint64_t indx[6];
    for ( indx[5] = start[5]; indx[5] < end[5]; ++indx[5] ) {
      for ( indx[4] = start[4]; indx[4] < end[4]; ++indx[4] ) {
        for ( indx[3] = start[3]; indx[3] < end[3]; ++indx[3] ) {
          for ( indx[2] = start[2]; indx[2] < end[2]; ++indx[2] ) {
            for ( indx[1] = start[1]; indx[1] < end[1]; ++indx[1] ) {
              for ( indx[0] = start[0]; indx[0] < end[0]; ++indx[0] ) {
                if (spinorb)
                  *p_Data = threeel_spi(indx[order[0]],indx[order[1]],indx[order[2]],indx[order[3]],indx[order[4]],indx[order[5]]);
                else
                  *p_Data = threeel_spa(indx[order[0]],indx[order[1]],indx[order[2]],indx[order[3]],indx[order[4]],indx[order[5]],aaaa);
                ++p_Data;
              }
            }
          }
        }
      }
    }
  }
  else {
    error("Number of indices is neither 2 nor 4 nor 6","Hdump::get_block");
  }
}
inline void Hdump::set_block(double* pData, const BlockIndices& start, const BlockIndices& end,
                             const BlockIndices& order, bool spinorb, uint spintype)
{
  assert(start.size() == end.size() && start.size() == order.size());
  // other spintypes defined only in _uhf case!
  assert(_uhf || spintype == 0);
  double * p_Data = pData;
  if (start.size() == 2) {
    // h_pq
    uint64_t indx[2];
    for ( indx[1] = start[1]; indx[1] < end[1]; ++indx[1] ) {
      for ( indx[0] = start[0]; indx[0] < end[0]; ++indx[0] ) {
        if (spinorb)
          set_oneel_spi(indx[order[0]],indx[order[1]],*p_Data);
        else
          set_oneel_spa(indx[order[0]],indx[order[1]],*p_Data,onetype(spintype));
        ++p_Data;
      }
    }
  } else if (start.size() == 4) {
    // (pq|rs)
    uint64_t indx[4];
    for ( indx[3] = start[3]; indx[3] < end[3]; ++indx[3] ) {
      for ( indx[2] = start[2]; indx[2] < end[2]; ++indx[2] ) {
        for ( indx[1] = start[1]; indx[1] < end[1]; ++indx[1] ) {
          for ( indx[0] = start[0]; indx[0] < end[0]; ++indx[0] ) {
            if (spinorb)
              set_twoel_spi(indx[order[0]],indx[order[1]],indx[order[2]],indx[order[3]],*p_Data);
            else
              set_twoel_spa(indx[order[0]],indx[order[1]],indx[order[2]],indx[order[3]],*p_Data,twotype(spintype));
            ++p_Data;
          }
        }
      }
    }
  } else {
    error("Number of indices is neither 2 nor 4","Hdump::set_block");
  }
}

inline double Hdump::spinprojector(uint p, uint q, Spin spintype) const {
  assert( p < _rd.refso.size() && q < _rd.refso.size() );
  if ( p != q ) return 0.0;
  if ( _rd.refso[p].spin == spintype ) return 1.0;
  return 0.0;
}
inline void Hdump::get_block_spinprojector(double* pData, const BlockIndices& start, const BlockIndices& end,
                             const BlockIndices& order, uint spintype) const
{
  assert(start.size() == end.size() && start.size() == order.size());
  if ( start.size() != 2 )
    error("Number of indices is not 2","Hdump::get_block_spinprojector");
  double * p_Data = pData;
  // P_pq
  uint64_t indx[2];
  for ( indx[1] = start[1]; indx[1] < end[1]; ++indx[1] ) {
    for ( indx[0] = start[0]; indx[0] < end[0]; ++indx[0] ) {
      *p_Data = spinprojector(indx[order[0]],indx[order[1]],Spin(spintype));
      ++p_Data;
    }
  }
}

inline double Hdump::spinsum_2DM(uint p, uint q, uint r, uint s) {
  assert( p < _norb && q < _norb && r < _norb && s < _norb );
  double twoRDM;
  if ( _pgs.nIrreps() > 1 ) error("spin summation of 2DM only implemented for nosym.");
  twoRDM = _twoel[aaaa].get()->get_with_pgs(p,q,r,s)
          +_twoel[aabb].get()->get_with_pgs(p,q,r,s)
          +_twoel[aabb].get()->get_with_pgs(r,s,p,q)
          +_twoel[bbbb].get()->get_with_pgs(p,q,r,s);
  return twoRDM;
}

inline double Hdump::spinsum_1DM(uint p, uint q, Irrep ir) {
  assert( p < _norb && q < _norb );
  double oneRDM;
  oneRDM = _oneel[aa].get()->get(p,q,ir)
          +_oneel[bb].get()->get(p,q,ir);
  return oneRDM;
}

/*!
    Overlap dump with point-group symmetry
*/
struct Overlapdump {
public:
  Overlapdump() {};
  // construct from file
  Overlapdump(std::string ovdump, const PGSym& pgs, const FDPar& ncore, bool verbose = true);
  // construct unity matrix (for restricted orbitals)
  Overlapdump(const PGSym& pgs);
  Integ2ab overlap;
};

/*!
 * Hubbard Site
 */
struct HubSite : public std::vector<int> {
  HubSite() : std::vector<int>() {};
  HubSite(const std::vector<uint>& dims, const std::vector<int>& pbcs) {
    _dims = dims;
    resize(dims.size(),0);
    _pbcs = pbcs;
    assert(_dims.size() == _pbcs.size());
  }
  void zero() {
    for (auto & i: *this) i = 0;
  }
  bool next() {
    assert(_dims.size() == size());

    for ( uint id = 0; id < size(); ++id ) {
      ++(*this)[id];
      if ((*this)[id] < int(_dims[id]) )
        return true;
      else {
        (*this)[id] = 0;
      }
    }
    return false;
  }
  // square of the distance
  uint dist2(const HubSite& hs) const {
    assert(hs.size() == size());
    assert(_dims.size() == size());
    assert(_dims.size() == _pbcs.size());
    uint dd = 0;
    for ( uint i = 0; i < size(); ++i ) {
      uint r = std::abs((*this)[i] - hs[i]);
      if ( _pbcs[i] && r > _dims[i] - r )
        r = _dims[i] - r;
      dd += r*r;
    }
    return dd;
  }
  uint id() const {
    assert(_dims.size() == size());
    uint id = 0;
    uint nn = 1;
    for ( uint i = 0; i < size(); ++i ) {
      id += (*this)[i]*nn;
      nn *= _dims[i];
    }
    return id;
  }
  // dimensions of the Hubbard model
  std::vector<uint> _dims;
  // periodic boundary condition
  std::vector<int> _pbcs;
};
std::ostream & operator << (std::ostream& o, const HubSite& hs);

} //namespace HamDump

#endif
