#include <vector>
#include "integs.h"
#include "density_mat.h"
#include "hdump.h"
#include <iomanip>
class Fock_matrices{
public:
 Fock_matrices() {};
 Fock_matrices(const Hdump& hdump, const DMdump& dmdump);
 int oneif(int p, int t) const;
 
private:
  uint _nsorb;
  std::vector<double> _FMAT;
 
    
};
