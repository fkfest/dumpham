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
  DMdump() : _nsorb(0) {};
  //DMdump for Hartree-Fock RDMs
  DMdump(uint norb, const Occupation& occ);
  //Dmdump for RDMs from fciqmc
  DMdump(const std::string filename, uint norb, uint nelec);

  bool nextdm(std::ifstream& inFile, int& i,int& j, int& k, int& l, double& S);
  inline int onei(int a, int b, int c, int d, int& sign) const;
  inline int oneid(int a, int c) const;
  void read_2rdm(std::string filename, Spin sa, Spin sb, Spin sc, Spin sd);
  void store_rdm() const;
  double value(uint p, uint q) const {uint indx = oneid(p,q); assert(indx<_RDM1.size()); return _RDM1[indx];}
  double value(uint p, uint q, uint r, uint s) const {int sign; uint indx = onei(p,q,r,s,sign); assert(indx<_RDM2.size()); return sign * _RDM2[indx];}
  std::vector<double> _RDM2;
  std::vector<double> _RDM1;
private:
  
  uint _nsorb;
};


inline int DMdump::onei(int a, int b, int c, int d, int& sign) const
/* function to obtain one unique index out of four spin-orbital indices.
Density matrix D^{a b}_{c d} = <a^+ b^+ d c>
*/
{
  sign = 1;
  int ab, cd, abcd;
  if(a>b){
  ab = a * (a+1)/2 + b;
  }
  else{
  ab = b * (b+1)/2 + a;
  sign = -sign;
  }
  if(c>d){
  cd = c * (c+1)/2 + d; 
  }
  else{
  cd = d * (d+1)/2 + c;
  sign = -sign;
  }
  if(ab > cd){
  abcd = ab * (ab +1)/2 + cd;  
  }
  else{
  abcd = cd * (cd+1)/2 + ab;
  }
  return abcd;
}

inline int DMdump::oneid(int a, int c) const
{
  int ac;
  if(a>c){
  ac = a * (a+1)/2 + c;
  }
  else{
  ac = c * (c+1)/2 + a;
  }
  return ac;
}




#endif