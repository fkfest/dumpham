#include "fock_mat.h"
#include <numeric>
//SUBROUTINE DGGEV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO)
  extern "C" {
  extern void dggev_(char*, char*, uint*, double*, uint*, const double*, uint*, double*, double*, double*, double*, uint*, double*, uint*, double*, int*, int*);
  extern void dgemm_(char*, char*, int*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);
  }

//construct fock matrix
Fock_matrices::Fock_matrices(const Hdump& hdump, const DMdump& dmdump)
{
  _nsorb = hdump.norb() * 2; 
  p_pgs = &(hdump.pgs());
  p_dump = &hdump;
  _FMAT.resize(p_pgs->nIrreps());
  // energy from the density matrices
  double Tr = 0.0;
  for ( Irrep ir = 0; ir < p_pgs->nIrreps(); ++ir ) {
    uint nsorb4ir = p_pgs->norbs(ir) * 2;
    if (nsorb4ir == 0) continue;
    uint ioff4ir = p_pgs->_firstorb4irrep[ir] * 2;
    uint maxid = nsorb4ir-1;
    uint nelem = oneif4ir(maxid,maxid,nsorb4ir)+1;
    std::vector<double> & fmat = _FMAT[ir];
    fmat.resize(nelem,0);
    std::vector<double> fmat2(nelem,0);

    for(uint p4ir = 0; p4ir < nsorb4ir; p4ir++){
      uint p = p4ir+ioff4ir;
      for(uint t4ir = 0; t4ir < nsorb4ir; t4ir++){
        uint t = t4ir+ioff4ir;
        uint pt = oneif4ir(p4ir,t4ir,nsorb4ir);
        for(uint q = 0; q < _nsorb; q++){
          fmat[pt] += hdump.oneel_spi_pgs(t,q) * dmdump.value(p,q);
        }
      }
    }
    std::vector<double> dm_qrs(nsorb4ir);
    for(uint q = 0; q < _nsorb; q++){
      for(uint r = 0; r < _nsorb; r++){
        for(uint s = 0; s < _nsorb; s++){
          for(uint t4ir = 0; t4ir < nsorb4ir; t4ir++){
            uint t = t4ir+ioff4ir;
            dm_qrs[t4ir] = dmdump.value(t,q,r,s);
          }
          for(uint p4ir = 0; p4ir < nsorb4ir; p4ir++){
            uint p = p4ir+ioff4ir;
            double twoel_prqs = hdump.twoel_spi_pgs(p,r,q,s); 
            for(uint t4ir = 0; t4ir < nsorb4ir; t4ir++){
              fmat2[oneif4ir(t4ir,p4ir,nsorb4ir)] += twoel_prqs * dm_qrs[t4ir];
            }
          }
        }
      }
    }
    // add two electron parts to fock
    for(uint pt = 0; pt < fmat.size(); pt++){
      fmat[pt] += fmat2[pt]; 
    }
    xout << nsorb4ir << std::endl;
    for(uint p = 0; p < nsorb4ir; p++){
      uint pt = oneif4ir(p,p,nsorb4ir);
      Tr = Tr + fmat[pt]-0.5*fmat2[pt];
    }
//    xout << " Fock " << ir << std::endl;
//    for(uint p4ir = 0; p4ir < nsorb4ir; p4ir++){
//      for(uint t4ir = 0; t4ir < nsorb4ir; t4ir++){
//        uint pt = oneif4ir(p4ir,t4ir,nsorb4ir);
//        xout << fmat[pt] << " ";
//      }
//      xout << std::endl;
//    }
//    xout << " RDM1 " << ir << std::endl;
//    for(uint p4ir = 0; p4ir < nsorb4ir; p4ir++){
//      uint p = p4ir+ioff4ir;
//      for(uint t4ir = 0; t4ir < nsorb4ir; t4ir++){
//          uint q = t4ir+ioff4ir;
//          xout << dmdump.value(p,q) << " ";
//      }
//      xout << std::endl;
//    }
  }
  xout << "Energy:" << Tr+hdump.escal() << std::endl;
}
void Fock_matrices::store(const std::string& filename) const
{
  std::ofstream outFile;
  outFile.open(filename.c_str());
  if ( (outFile.rdstate() & std::ifstream::failbit ) != 0 ) {
      outFile.close();
      error("Fock_matrices::store failed to open "+ filename);
  }
  
  for(uint p = 0; p < _nsorb; p++){
    Irrep pir = p_pgs->irrep(p/2); 
    uint ns4ir = p_pgs->norbs(pir) * 2;
    outFile << " " << std::endl;
    for(uint t = 0; t < _nsorb; t++){
      Irrep tir = p_pgs->irrep(t/2);
      double val = 0;
      if ( pir == tir ) {
        uint p4ir = p - 2 * p_pgs->_firstorb4irrep[pir];
        uint t4ir = t - 2 * p_pgs->_firstorb4irrep[tir];
        uint pt = oneif4ir(p4ir,t4ir,ns4ir);
        val =  _FMAT[pir][pt];
      }
      outFile << std::setw(15) << std::left << std::fixed << val;
    }
  }
  outFile.close();
}

void Fock_matrices::diagonalize(const DMdump& dmdump)
{
  std::string ektorb = Input::sPars["orbs"]["ektorb"];
  std::string dysorb = Input::sPars["orbs"]["dysorb"];
  std::string ips = Input::sPars["fm"]["ips"];
  std::string orbsin = Input::sPars["orbs"]["in"];
  bool
    calc_ektorb = !ektorb.empty(),
    calc_dysorb = !dysorb.empty(),
    store_ips = !ips.empty();
  Integ2ab ekttrma, dystrma, ipsmat;
  if ( calc_ektorb || calc_dysorb ) {
    if ( orbsin.empty() ) {
      error("Provide original orbitals in order to calculate EKT or Dyson orbs");
    }
    if ( calc_ektorb ) ekttrma = Integ2ab(*p_pgs);
    if ( calc_dysorb ) dystrma = Integ2ab(*p_pgs);
  }
  if ( store_ips ) ipsmat = Integ2ab(*p_pgs);
  for ( Irrep ir = 0; ir < p_pgs->nIrreps(); ++ir ) {
    uint norb4ir = p_pgs->norbs(ir),
         ioff4ir = p_pgs->_firstorb4irrep[ir];
    if (norb4ir == 0) continue;
    std::vector<double> eigval,eigvec,seigvec;
    diagonalize4irrep(eigval,eigvec,seigvec,dmdump,ir);
    for ( uint i4ir = 0; i4ir < norb4ir; ++i4ir ) {
      for ( uint j4ir = 0; j4ir < norb4ir; ++j4ir ) {
        if ( calc_ektorb ) {
          ekttrma.set(i4ir+ioff4ir,j4ir+ioff4ir,spinsummedvalue(eigvec,i4ir,j4ir,norb4ir));
        }
        if ( calc_dysorb ) {
          dystrma.set(i4ir+ioff4ir,j4ir+ioff4ir,spinsummedvalue(seigvec,i4ir,j4ir,norb4ir));
        }
      }
      if ( store_ips ) {
        ipsmat.set(i4ir+ioff4ir,i4ir+ioff4ir,eigval[2*i4ir]);
      }
    }
  }
  if ( calc_ektorb ) {
    Odump odump(p_dump->pgs_wcore(),p_dump->ncore(),orbsin);
    odump.transform(ekttrma);
    odump.store(ektorb);
  }
  if ( calc_dysorb ) {
    Odump odump(p_dump->pgs_wcore(),p_dump->ncore(),orbsin);
    odump.transform(dystrma);
    odump.store(dysorb);
  }
  if ( store_ips ) {
    Odump ipdump(p_dump->pgs_wcore(),p_dump->ncore(),ipsmat);
    ipdump.store(ips);
  }
}

#ifdef _LAPACK
// make the largest coeff of the vector positive
static void make_largest_positive(double * pBeg, uint nElem) {
  double maxelem = 0.0;
  bool change_sign = false;
  for ( double * pVec = pBeg; pVec != pBeg + nElem; ++pVec ) {
    if ( maxelem < std::abs(*pVec) ) {
      maxelem = std::abs(*pVec);
      change_sign = ( *pVec < 0.0);
    }
  }
  if ( change_sign ) {
    for ( double * pVec = pBeg; pVec != pBeg + nElem; ++pVec ) {
      *pVec = - *pVec;
    }
  }
}
static bool normalize(double * pBeg, uint nElem) {
  double norm = 0.0;
  for ( double * pVec = pBeg; pVec != pBeg+nElem; ++pVec ) {
    norm += *pVec * *pVec;
  }
  if (norm < Numbers::verysmall) return false;
  norm = 1.0/sqrt(norm);
  for ( double * pVec = pBeg; pVec != pBeg+nElem; ++pVec ) {
    *pVec *= norm;
  }
  return true;
}
#endif

void Fock_matrices::diagonalize4irrep(std::vector< double >& eigvalues, std::vector< double >& eigvectors, 
                                      std::vector<double>& seigvec, const DMdump& dmdump, Irrep ir)
{
#ifdef _LAPACK
//Calculate eigenvalues using LAPACK DGGEV for generalized eigenvalue problem
  uint nsorb4ir = p_pgs->norbs(ir) * 2;
  uint ioff4ir = p_pgs->_firstorb4irrep[ir] * 2;
  std::vector<double> eigReal(nsorb4ir,1);
  std::vector<double> eigImag(nsorb4ir,1);
  std::vector<double> beta(nsorb4ir,0);
  std::vector<double> v(nsorb4ir*nsorb4ir,0);
  int lwork=-1;
  double workdummy;
  int info = 0;
  char Nchar = 'N';
  char Vchar = 'V';
  
  std::vector<double> RDM1(nsorb4ir*nsorb4ir);
  for ( uint i4ir = 0; i4ir < nsorb4ir; ++i4ir ){
    uint i = i4ir + ioff4ir;
    for ( uint j4ir = 0; j4ir < nsorb4ir; ++j4ir ){
      uint j = j4ir + ioff4ir;
      RDM1[oneif4ir(i4ir,j4ir,nsorb4ir)] = dmdump.value(i,j);
    }
  }
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
  dggev_(&Nchar, &Vchar, &nsorb4ir, _FMAT[ir].data(), &nsorb4ir, RDM1.data(), &nsorb4ir, 
         eigReal.data(), eigImag.data(), beta.data(), 0, &nsorb4ir, v.data(), &nsorb4ir, &workdummy, &lwork, &info);
  lwork = int(workdummy) + 64;
  std::vector<double> work;
  work.resize(lwork,0);
//   xout << "lwork:" << lwork << " " << "workdummy:" << workdummy << std::endl;
  dggev_(&Nchar, &Vchar, &nsorb4ir, _FMAT[ir].data(), &nsorb4ir, RDM1.data(), &nsorb4ir, 
         eigReal.data(), eigImag.data(), beta.data(), 0, &nsorb4ir, v.data(), &nsorb4ir, work.data(), &lwork, &info);
  // check for errors
  if (info!=0){
    xout << "Error: dggev returned error code " << info << std::endl;
    exit(1);
  }
  std::vector< std::pair<double,uint> > eigvalreal;
//   output eigenvalues to stdout
  xout << "--- Eigenvalues in symmetry " << ir << " ---" << std::endl;
  for (uint i = 0; i < nsorb4ir; i++){
    double den = beta[i];
    if (std::abs(den) < Numbers::small) {
      xout << "Small denominator for value: ";
      den = 1.0;
    }
    xout << i+1 << " " <<"( " << eigReal[i]/den << " , " << eigImag[i]/den << " )\n";
    eigvalreal.push_back(std::make_pair(-eigReal[i]/den,i)); 
  }
  std::sort(eigvalreal.begin(),eigvalreal.end());
  xout << "Sorted eigenvalues:" << std::endl;
  for (uint i = 0; i < nsorb4ir; i++){
    xout << eigvalreal[i].first << " " << eigvalreal[i].second << std::endl;
  }
  // return in sorted order
  eigvalues.resize(nsorb4ir);
  eigvectors.resize(nsorb4ir*nsorb4ir);
  for (uint i = 0; i < nsorb4ir; i++){
    eigvalues[i] = eigvalreal[i].first;
    uint ivec = eigvalreal[i].second;
    std::copy(v.begin()+ivec*nsorb4ir,v.begin()+(ivec+1)*nsorb4ir,eigvectors.begin()+i*nsorb4ir);
    make_largest_positive(&(eigvectors[i*nsorb4ir]),nsorb4ir);
  }
  seigvec.resize(nsorb4ir*nsorb4ir);
  char Tchar = 'T';
  double One = 1.0, Zero = 0.0;
  int len = nsorb4ir;
  // restore RDM back
  for ( uint i4ir = 0; i4ir < nsorb4ir; ++i4ir ){
    uint i = i4ir + ioff4ir;
    for ( uint j4ir = 0; j4ir < nsorb4ir; ++j4ir ){
      uint j = j4ir + ioff4ir;
      RDM1[oneif4ir(i4ir,j4ir,nsorb4ir)] = dmdump.value(i,j);
    }
  }
  dgemm_(&Tchar,&Nchar,&len,&len,&len,&One,RDM1.data(),&len,eigvectors.data(),&len,&Zero,seigvec.data(),&len);
  // norm
  for (uint i = 0; i < nsorb4ir; i++){
    if (!normalize(&(seigvec[nsorb4ir*i]),nsorb4ir)){
      warning("Small norm of Dyson transform vector number " << i);
    }
  }
}
#else
 (void) eigvalues; (void) eigvectors; (void) seigvec; (void) dmdump; (void) ir;
  error("Cannot diagonalize: Compiled without LAPACK","Fock_matrices::diagonalize4irrep");
}
#endif
