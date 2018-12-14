#ifndef dichtematrix
#define dichtematrix
#include <string>
#include <vector>
#include "globals.h"
#include "utilities.h"
#include "odump.h"
#include <fstream>
#include <iomanip>
#include <typeinfo>
#include <sstream>

class DMdump
{
public:
  /*!
     * \brief Construct an empty FCIdump object
     */
  
  //DMdump for Hartree-Fock RDMs
  DMdump(uint norb, const Occupation& occ);
  //Dmdump for RDMs from fciqmc
  DMdump(const std::string filename, uint norb, uint nelec);

  bool nextdm(std::ifstream& inFile, int& i,int& j, int& k, int& l, double& S);
  int onei(int a, int b, int c, int d, int& sign) const;
  int oneid(int a, int c) const;
  void read_2rdm(std::string filename, Spin sa, Spin sb, Spin sc, Spin sd);
  void store_rdm() const;
  double value(uint p, uint q) const {uint indx = oneid(p,q); assert(indx<_RDM1.size()); return _RDM1[indx];}
  double value(uint p, uint q, uint r, uint s) const {int sign; uint indx = onei(p,q,r,s,sign); assert(indx<_RDM2.size()); return sign * _RDM2[indx];}
  std::vector<double> _RDM2;
  std::vector<double> _RDM1;
private:
  
  uint _nsorb;
};






#endif