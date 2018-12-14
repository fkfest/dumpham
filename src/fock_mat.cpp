#include "fock_mat.h"
//SUBROUTINE DGGEV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO)
  extern "C" {
  extern void dggev_(char*, char*, uint*, double*, uint*, const double*, uint*, double*, double*, double*, double*, uint*, double*, uint*, double*, int*, int*);
  }

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
	    // 1.0 when printing FMAT; 0.5 when calculating trace 
	    _FMAT[pt] += hdump.twoel(t,r,q,s) * dmdump.value(p,q,r,s);
// 	    _FMAT[pt] += 0.5*(hdump.twoel(t,r,q,s) * dmdump.value(p,q,r,s));
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

  
//Calculate eigenvalues using LAPACK DGGEV for generalized eigenvalue problem
  uint nmax = 52;
  std::vector<double> eigReal;
  eigReal.resize(_nsorb,1);
  std::vector<double> eigImag;
  eigImag.resize(_nsorb,1);
  std::vector<double> beta;
  beta.resize(_nsorb,0);
  std::vector<double> v;
  v.resize(_nsorb*nmax,0);
  int lwork=-1;
  double workdummy;
  int info = 0;
  double zero = 0;
  char Nchar = 'N';
  char Vchar = 'V';
//   //write the vectors in matrices
//   std::vector<std::vector<double> > FM;
//   FM.resize(_nsorb);
//   for(int i =0; i < _nsorb; i++)
//     FM<i>.resize(_nsorb);
//   
//   for(uint i=0; i < _nsorb;i++){
//    for(uint j=0; j<_nsorb;j++){
//     FM[i][j] = _FMAT[i*_nsorb+j]; 
//    }
//   }
  // calculate eigenvalues using the DGEEV subroutine
  dggev_(&Nchar, &Vchar, &_nsorb, _FMAT.data(), &_nsorb, dmdump._RDM1.data(), &_nsorb, eigReal.data(), eigImag.data(), beta.data(), &zero, &nmax, v.data(), &nmax, &workdummy, &lwork, &info);
  lwork = int(workdummy) + 64;
  std::vector<double> work;
  work.resize(lwork,0);
  xout << "lwork:" << lwork << " " << "workdummy:" << workdummy << std::endl;
  dggev_(&Nchar, &Vchar, &_nsorb, _FMAT.data(), &_nsorb, dmdump._RDM1.data(), &_nsorb, eigReal.data(), eigImag.data(), beta.data(), &zero, &nmax, v.data(), &nmax, work.data(), &lwork, &info);
  // check for errors
  if (info!=0){
    xout << "Error: dggev returned error code " << info << std::endl;
    exit(1);
  }
//   output eigenvalues to stdout
  xout << "--- Eigenvalues ---" << std::endl;
  for (uint i=0;i<_nsorb;i++){
    xout << i+1 << " " <<"( " << eigReal[i] << " , " << eigImag[i] << " )\n";
  }
}

int Fock_matrices::oneif(int p, int t) const
{
  int pt;
  pt = p + t*_nsorb;
return pt;
}