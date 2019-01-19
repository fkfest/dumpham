#ifndef Hdump_H
#define Hdump_H
#include <string>
#include <vector>
#include <memory>
#include "globals.h"
#include "utilities.h"
#include "FCIdump.h"
#include "pgsym.h"
#include "integs.h"
/*!
    Hamiltonian dump without point-group symmetry
*/
class Hdump {
public:  
  Hdump():_escal(0),_norb(0),_nelec(0),_ms2(0){};
  // construct from FCIdump
  Hdump(std::string fcidump);
  enum onetype {
    aa = 0,
    bb = 1
  };
  enum twotype {
    aaaa = 0,
    bbbb = 1,
    aabb = 2
  };
  void store(std::string fcidump);
  uint norb() const { return _norb; }
  uint nelec() const { return _nelec; }
  uint ms2() const { return _ms2; }
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
  // in spin orbitals
  inline double oneel(uint p, uint q) const; 
  // in spin orbitals
  inline double twoel(uint p, uint q, uint r, uint s) const; 
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
  void check_input_norbs(FDPar& orb, const std::string& kind ) const;
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
  // point group symmetry
  PGSym _pgs, _pgs_wcore;
  FDPar _occ, _clos, _core;
  // dump file in uhf orbitals
  bool _uhf;
  // doesn't have the bra-ket symmetry
  bool _simtra;
};

inline double Hdump::oneel(uint p, uint q) const {
  Spin sp = spin(p), sq = spin(q);
  if( sp != sq) return 0;
  uint i = p/2, j = q/2;
  if(!_uhf || sp == alpha)
     return static_cast<Integ2*>(_oneel[aa].get())->get_with_pgs(i,j);
  return static_cast<Integ2*>(_oneel[bb].get())->get_with_pgs(i,j);
}

inline double Hdump::twoel(uint p, uint q, uint r, uint s) const {
  Spin sp = spin(p), sq = spin(q), sr = spin(r), ss = spin(s);
  if(sp != sq || sr != ss) return 0;
  uint i = p/2, j = q/2, k = r/2, l = s/2;
  if(!_uhf || ( sp == alpha && sr == alpha )) 
    return static_cast<Integ4*>(_twoel[aaaa].get())->get_with_pgs(i,j,k,l);
  if( sp == alpha ) 
    return static_cast<Integ4*>(_twoel[aabb].get())->get_with_pgs(i,j,k,l);
  return static_cast<Integ4*>(_twoel[bbbb].get())->get_with_pgs(i,j,k,l);
}


#endif
