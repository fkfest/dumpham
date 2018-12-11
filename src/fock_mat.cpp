#include "fock_mat.h"
//construct fock matrix
Fock_matrices::Fock_matrices(const Hdump& hdump, const DMdump& dmdump)
{
  _nsorb = hdump.norb() * 2; 
  uint maxid = _nsorb-1;
  uint nelem = oneif(maxid,maxid)+1;
  _FMAT.resize(nelem,0);
  std::ofstream outFile;

  for(uint p = 0; p < _nsorb; p++){
    for(uint t = 0; t < _nsorb; t++){
      uint pt = oneif(p,t);
      for(uint q = 0; q < _nsorb; q++){
	_FMAT[pt] += hdump.oneel(t,q) * dmdump.value(p,q);
	for(uint r = 0; r < _nsorb; r++){
	  for(uint s = 0; s < _nsorb; s++){
	    _FMAT[pt] += 0.5*(hdump.twoel(t,r,q,s) * dmdump.value(p,q,r,s)); 
	  }
	}
      }
    }
  }

  outFile.open("Fmat.dat");
  for(uint p = 0; p < _nsorb; p++){
   outFile << " " << std::endl;
   for(uint t = 0; t < _nsorb; t++){  
    uint pt = oneif(p,t);
    outFile << std::setw(15) << std::left << std::fixed << _FMAT[pt];
   }
  }
  outFile.close();
  
  double Tr = 0.0;
  for(uint p = 0; p < _nsorb; p++){
    uint pt = oneif(p,p);
    Tr = Tr + _FMAT[pt];
  }
  
  xout << "Trace:" << Tr+hdump.escal() << std::endl;

}

int Fock_matrices::oneif(int p, int t) const
{
  int pt;
  pt = p + t*_nsorb;
return pt;
}