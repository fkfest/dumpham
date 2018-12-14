#include <vector>
#include "integs.h"
#include "density_mat.h"
#include "hdump.h"
#include <iomanip>

class Fock_matrices{
public:
 Fock_matrices() {};
 Fock_matrices(const Hdump& hdump, const DMdump& dmdump);
 int oneif(uint p, uint t) const { assert( p < _nsorb && t < _nsorb ); return p + t*_nsorb; }
 void diagonalize(const DMdump& dmdump);
 void store(const std::string& filename) const;
private:
  uint _nsorb;
  std::vector<double> _FMAT;
};
