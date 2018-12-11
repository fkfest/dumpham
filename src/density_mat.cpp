#include "density_mat.h"

DMdump::DMdump(uint norb, uint nelec)
{
//build HF-1RDM
uint nelem1d;
nelem1d = oneid(nelec,nelec);
_HFRDM1.resize(nelem1d,0);  

}

DMdump::DMdump(const std::string filename, uint norb, uint nelec)
{
 
  int sign;
  uint nelem;
  uint nelem_1d;
  _nsorb = 2 * norb;
  int maxid = _nsorb -1;
  nelem = onei(maxid,maxid,maxid,maxid,sign)+1;
  _RDM2.resize(nelem,0);
  xout << "nelem " << nelem << std::endl;
  nelem_1d = oneid(maxid,maxid)+1;
  xout << "nelem_1d " << nelem_1d << std::endl; 
  _RDM1.resize(nelem_1d,0);
  double Tr;
  int shell;
  shell = 0;
  std::string f_aaaa, f_abba, f_abab, f_bbbb, f_baab, f_baba;
  std::string str = "aaaa";
  f_aaaa = filename;
  f_abba = f_aaaa;
  f_abba.replace(f_abba.find(str),str.length(),"abba");
  f_abab = f_aaaa;
  f_abab.replace(f_abab.find(str),str.length(),"abab");
  if(shell == 1){
    //open shell
    f_bbbb = f_aaaa;
    f_bbbb.replace(f_bbbb.find(str),str.length(),"bbbb");
    f_baab = f_aaaa;
    f_baab.replace(f_baab.find(str),str.length(),"baab");
    f_baba = f_aaaa;
    f_baba.replace(f_baba.find(str),str.length(),"baba"); 
  } else {
    //for closed shell aaaa = bbbb, abab = baba, abba = baab
    f_bbbb = f_aaaa;
    f_baab = f_abba;
    f_baba = f_abab;
  }
  
  read_2rdm(f_aaaa,alpha,alpha,alpha,alpha);
  read_2rdm(f_abab,alpha,beta,alpha,beta);
  read_2rdm(f_abba,alpha,beta,beta,alpha);
  read_2rdm(f_bbbb,beta,beta,beta,beta);
  read_2rdm(f_baba,beta,alpha,beta,alpha);
  read_2rdm(f_baab,beta,alpha,alpha,beta);
    
  uint ac, abcd;
  //Construct RDM1
  for(uint a = 0; a < _nsorb; a++){
    for(uint c = 0; c <= a; c++){
      ac = oneid(a,c);
      assert(ac < nelem_1d);
      for(uint b = 0; b < _nsorb; b++){
        abcd = onei(a,b,c,b,sign);
        //std::cout << ac << std::endl;
        assert(abcd < nelem);
        _RDM1[ac] = _RDM1[ac] + sign*_RDM2[abcd];
      }
      _RDM1[ac] = _RDM1[ac] / (nelec-1);
    }
  }

  //calculate trace of 1RDM 
  Tr = 0.0;
  for(uint a = 0; a < _nsorb; a++){
    ac = oneid(a,a);
    Tr = Tr + _RDM1[ac];
  }
  //std::cout << typeid(Tr).name() << '\n';
  std::cout << std::setprecision(16) << std::fixed << Tr << std::endl; 

  store_rdm();
 
}
void DMdump::read_2rdm(std::string filename, Spin sa, Spin sb, Spin sc, Spin sd)
{
  std::ifstream inFile;
  inFile.open(filename);
  if (!inFile){
    std::cout << "Unable to open file";
    exit(1);   // call system to stop
  }
  int i,j,k,l,sign;
  uint a,b,c,d;
  uint abcd;
  double S;
  while (nextdm(inFile,i,j,k,l,S)){
    a = 2*(i-1)+sa;
    b = 2*(j-1)+sb;
    c = 2*(k-1)+sc;
    d = 2*(l-1)+sd;
    abcd = onei(a,b,c,d,sign);
    assert(abcd < _RDM2.size());
#ifdef _DEBUG
    if (std::abs(_RDM2[abcd]) > 1.e-10 && std::abs(_RDM2[abcd]-S * sign) > 1.e-10) 
      xout << _RDM2[abcd] << " " << S * sign << " " << a << " " << b << " " << c << " " << d << std::endl;
#endif
    _RDM2[abcd] = S * sign;
  }
}

void DMdump::store_rdm() const
{
  std::ofstream outFile;
  uint ac;
  std::string out1rdm = Input::sPars["dm"]["out1rdm"];
  std::string out2rdm = Input::sPars["dm"]["out2rdm"];
  std::string out1rdmmat = Input::sPars["dm"]["out1rdmmat"];
  if ( out2rdm != "" ) {
    //write RDM2 to a file
    outFile.open(out2rdm.c_str());
    if ( (outFile.rdstate() & std::ifstream::failbit ) != 0 ) {
      outFile.close();
      error("DMdump::store_rdm failed to open "+ out2rdm);
    }
    int sign;
    uint abcd;
    for(uint a=0; a<=_nsorb-1; a++){
      for(uint b=0; b<=_nsorb-1; b++){
        for(uint c=0; c<=_nsorb-1; c++){
          for(uint d=0; d<=_nsorb-1; d++){
            abcd = onei(a,b,c,d,sign);
            outFile << std::setw(15) << sign*_RDM2[abcd] 
                    << std::setw(4) << a << std::setw(4) << b << std::setw(4) 
                    << c << std::setw(4) << d << std::endl;
          }
        }
      }
    }
    outFile.close();
  }
  if ( out1rdmmat != "" ) {
    // write RDM1 as a matrix to a file
    outFile.open(out1rdmmat.c_str());
    if ( (outFile.rdstate() & std::ifstream::failbit ) != 0 ) {
      outFile.close();
      error("DMdump::store_rdm failed to open "+ out1rdmmat);
    }
    for(uint a=0; a<=_nsorb-1; a++){
      outFile << " " << std::endl;
      for(uint c=0; c <=_nsorb-1; c++){
        ac = oneid(a,c);
        outFile << std::setw(15) << std::left << std::fixed << _RDM1[ac];
      }
    }
    outFile.close(); 
  }

  if ( out1rdm != "" ) {
    // write RDM1 to a file
    outFile.open(out1rdm.c_str());
    if ( (outFile.rdstate() & std::ifstream::failbit ) != 0 ) {
      outFile.close();
      error("DMdump::store_rdm failed to open "+ out1rdm);
    }
    for(uint a=0; a<=_nsorb-1; a++){
      for(uint c=0; c <=_nsorb-1; c++){
        ac = oneid(a,c);
        outFile << std::setw(15) << std::left << std::fixed << _RDM1[ac] 
                << std::setw(4) << a << std::setw(4) << c << std::endl;
      }
    }
    outFile.close(); 
  }
}

bool DMdump::nextdm(std::ifstream& inFile, int& i, int& j, int& k, int& l, double& S)
 {
  if(inFile >> i >> j >> k >> l >> S){
   return true;
  }
  else{
  inFile.close();
  return false;  
  }
 }
 
int DMdump::onei(int a, int b, int c, int d, int& sign) const
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

int DMdump::oneid(int a, int c) const
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