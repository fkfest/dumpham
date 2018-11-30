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
  uint nclosed() const;
  // number of open-shell orbitals
  uint nopen() const { return _ms2; }
  // number of occupied orbitals
  uint nocc() const { return nclosed()+nopen(); }
  
private:
  void store_with_symmetry();
  void store_without_symmetry();
  // on input: first value (i,j,k,l,value,type)
  template<typename T>
  void readrec(T * pInt, int& i, int& j, int& k, int& l, double& value, FCIdump::integralType& curtype );
  template<typename T>
  void readrec(T * pInt, int& i, int& j, double& value, FCIdump::integralType& curtype );
  void storerec_sym(Integ4 * pInt);
  void storerec_sym(Integ4ab * pInt);
  void storerec_sym(Integ2 * pInt);
  void storerec_nosym(Integ4 * pInt);
  void storerec_nosym(Integ4ab * pInt);
  void storerec_nosym(Integ2 * pInt);
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
  PGSym _pgs;
  FDPar _occ, _clos;
  // dump file in uhf orbitals
  bool _uhf;
  // doesn't have the bra-ket symmetry
  bool _simtra;
};


#endif
