#include <vector>
#include "integs.h"
#include "density_mat.h"
#include "hdump.h"
#include <iomanip>

// Fock matrix in spin orbitals
class Fock_matrices{
public:
 Fock_matrices() {};
 Fock_matrices(const Hdump& hdump, const DMdump& dmdump);
 int oneif(uint p, uint t) const { assert( p < _nsorb && t < _nsorb ); return p + t*_nsorb; }
 // with symmetry
 int oneif4ir(uint p, uint t, uint ns4ir) const { assert( p < ns4ir && t < ns4ir ); return p + t*ns4ir; }
 void diagonalize(const DMdump& dmdump);
 void store(const std::string& filename) const;
private:
  // diagonalize fock as FC=DCe. The eigenvalues/vectors will be sorted
  void diagonalize4irrep(std::vector<double>& eigvalues, std::vector<double>& eigvectors, std::vector<double>& seigvec,
                         const DMdump& dmdump, Irrep ir);
  uint _nsorb;
  // point group symmetry (in spatial orbitals!)
  const PGSym * p_pgs;
  // original hdump
  const Hdump * p_dump;
  // fock matrix for each symmetry
  std::vector< std::vector<double> >_FMAT;
};
