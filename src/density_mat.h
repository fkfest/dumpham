#ifndef dichtematrix
#define dichtematrix
#include <string>
#include <vector>
#include "globals.h"
#include "utilities.h"
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
  DMdump(uint norb, uint nelec);
  //Dmdump for RDMs from fciqmc
  DMdump(const std::string filename, uint norb, uint nelec);

  bool nextdm(std::ifstream& inFile, int& i,int& j, int& k, int& l, double& S);
  int onei(int a, int b, int c, int d, int& sign) const;
  int oneid(int a, int c) const;
  void read_2rdm(std::string filename, Spin sa, Spin sb, Spin sc, Spin sd);
  void store_rdm() const;
  double value(uint p, uint q) const {return _RDM1[oneid(p,q)];}
  double value(uint p, uint q, uint r, uint s) const {int sign, indx = onei(p,q,r,s,sign); return sign * _RDM2[indx];}
  std::vector<double> _RDM2;
private:
  std::vector<double> _HFRDM1;
  std::vector<double> _RDM1;
  uint _nsorb;
};






#endif