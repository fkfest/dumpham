#include <vector>
#include "integs.h"
#include "density_mat.h"
#include "hdump.h"
#include <iomanip>

using namespace HamDump;

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
  // return spin summed value from matrix mat for spatial orbitals i and j and norb spatial orbitals
  inline double spinsummedvalue(const std::vector<double> & mat, uint i, uint j, uint norb) const;
  uint _nsorb;
  // point group symmetry (in spatial orbitals!)
  const PGSym * p_pgs;
  // original hdump
  const Hdump * p_dump;
  // fock matrix for each symmetry
  std::vector< std::vector<double> >_FMAT;
};

inline double Fock_matrices::spinsummedvalue(const std::vector< double >& mat, uint i, uint j, uint norb) const
{
  return mat[oneif4ir(2*i,2*j,2*norb)]+mat[oneif4ir(2*i+1,2*j,2*norb)]+
         mat[oneif4ir(2*i,2*j+1,2*norb)]+mat[oneif4ir(2*i+1,2*j+1,2*norb)];
}