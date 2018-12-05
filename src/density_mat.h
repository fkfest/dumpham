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
  DMdump(const std::string filename, uint norb, uint nelec);

  bool nextdm(int& i,int& j, int& k, int& l, double& S);
  int onei(int a, int b, int c, int d, int& sign);
  int oneid(int a, int c);
  mutable std::ifstream inFile;
  mutable std::ofstream outFile;
  std::string f_aaaa, f_abba, f_abab, f_bbbb, f_baab, f_baba;
  
};






#endif