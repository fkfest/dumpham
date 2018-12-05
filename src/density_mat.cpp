#include "density_mat.h"


DMdump::DMdump(const std::string filename, uint norb, uint nelec)
{
 
  int i, j, k, l;
  uint a, b, c, d, ac, abcd;
  int sign;
  uint nelem;
  uint nelem_1d;
  double S;
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
  std::vector<double> RDM2aaaa(nelem,0);
  std::vector<double> RDM2abba(nelem,0);
  std::vector<double> RDM2abab(nelem,0);
  std::vector<double> RDM2bbbb(nelem,0);
  std::vector<double> RDM2baab(nelem,0);
  std::vector<double> RDM2baba(nelem,0);

  std::string str = "aaaa";
  f_aaaa = filename;
  f_abba = f_aaaa;
  f_abba.replace(f_abba.find(str),str.length(),"abba");
  f_abab = f_aaaa;
  f_abab.replace(f_abab.find(str),str.length(),"abab");
  std::ifstream inFile;
  
  inFile.open(f_aaaa);
  if (!inFile){
    std::cout << "Unable to open file";
    exit(1);   // call system to stop
  }
  while (nextdm(inFile,i,j,k,l,S)){
  a = 2*(i-1);
  b = 2*(j-1);
  c = 2*(k-1);
  d = 2*(l-1);
  abcd = onei(a,b,c,d,sign);
  assert(abcd < nelem);
  RDM2aaaa[abcd] = S * sign;
  if (std::abs(_RDM2[abcd]) > 1.e-10) 
    xout << _RDM2[abcd] << " " << S * sign << " " << a << " " << b << " " << c << " " << d << std::endl;
  _RDM2[abcd] = S * sign;
  
  }
  
  inFile.open(f_abab);
  if (!inFile){
    std::cout << "Unable to open file";
    exit(1);   // call system to stop
  }
  while (nextdm(inFile,i,j,k,l,S)){
  a = 2*(i-1);
  b = 2*j-1;
  c = 2*(k-1);
  d = 2*l-1;
  abcd = onei(a,b,c,d,sign);
  assert(abcd < nelem);
  RDM2abab[abcd] = S * sign;
  if (std::abs(_RDM2[abcd]) > 1.e-10) 
    xout << _RDM2[abcd] << " " << S * sign << " " << a << " " << b << " " << c << " " << d << std::endl;
  _RDM2[abcd] = S * sign;
  }
  
  inFile.open(f_abba);
  if (!inFile){
    std::cout << "Unable to open file";
    exit(1);   // call system to stop
  }
  while (nextdm(inFile,i,j,k,l,S)){
  a = 2*(i-1);
  b = 2*j-1;
  c = 2*k-1;
  d = 2*(l-1);
  abcd = onei(a,b,c,d,sign);
  assert(abcd < nelem);
  RDM2abba[abcd] = S * sign;
  if (std::abs(_RDM2[abcd]) > 1.e-10) 
    xout << _RDM2[abcd] << " " << S * sign << " " << a << " " << b << " " << c << " " << d << std::endl;
  _RDM2[abcd] = S * sign;
  }
  
 
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
    
  inFile.open(f_bbbb);
  if (!inFile){
    std::cout << "Unable to open file";
    exit(1);   // call system to stop
  }
  while (nextdm(inFile,i,j,k,l,S)){
  a = 2*i-1;
  b = 2*j-1;
  c = 2*k-1;
  d = 2*l-1;
  abcd = onei(a,b,c,d,sign);
  assert(abcd < nelem);
  RDM2bbbb[abcd] = S * sign;
  if (std::abs(_RDM2[abcd]) > 1.e-10) 
    xout << _RDM2[abcd] << " " << S * sign << " " << a << " " << b << " " << c << " " << d << std::endl;
  _RDM2[abcd] = S * sign;
  }
  
  inFile.open(f_baba);
  if (!inFile){
    std::cout << "Unable to open file";
    exit(1);   // call system to stop
  }
  while (nextdm(inFile,i,j,k,l,S)){
  a = 2*i-1;
  b = 2*(j-1);
  c = 2*k-1;
  d = 2*(l-1);
  abcd = onei(a,b,c,d,sign);
  assert(abcd < nelem);
  RDM2baba[abcd] = S * sign;
  if (std::abs(_RDM2[abcd]) > 1.e-10) 
    xout << _RDM2[abcd] << " " << S * sign << " " << a << " " << b << " " << c << " " << d << std::endl;
  _RDM2[abcd] = S * sign;
  }
  
  inFile.open(f_baab);
  if (!inFile){
    std::cout << "Unable to open file";
    exit(1);   // call system to stop
  }
  while (nextdm(inFile,i,j,k,l,S)){
  a = 2*i-1;
  b = 2*(j-1);
  c = 2*(k-1);
  d = 2*l-1;
  abcd = onei(a,b,c,d,sign);
  assert(abcd < nelem);
  RDM2baab[abcd] = S * sign;
  if (std::abs(_RDM2[abcd]) > 1.e-10) 
    xout << _RDM2[abcd] << " " << S * sign << " " << a << " " << b << " " << c << " " << d << std::endl;
  _RDM2[abcd] = S * sign;
  }

  
//   //write RDM2aaaa to a file
//   outFile.open("2RDMaaaa.txt");
//   for(a=0; a<=_nsorb-1; a++){
//     for(b=0; b<=_nsorb-1; b++){
//       for(c=0; c<=_nsorb-1; c++){
// 	for(d=0; d<=_nsorb-1; d++){
// 	abcd = onei(a,b,c,d,sign);
// 	outFile << sign*RDM2aaaa[abcd] << std::endl;
// 	}
//       }
//     }
//   }
//   outFile.close();
//   
//   //write RDM2abba to a file
//   outFile.open("2RDMabba.txt");
//   for(a=0; a<=_nsorb-1; a++){
//     for(b=0; b<=_nsorb-1; b++){
//       for(c=0; c<=_nsorb-1; c++){
// 	for(d=0; d<=_nsorb-1; d++){
// 	abcd = onei(a,b,c,d,sign);
// 	outFile << sign*RDM2abba[abcd] << std::endl;
// 	}
//       }
//     }
//   }
//   outFile.close();
//   
//   //write RDM2abab to a file
//   outFile.open("2RDMabab.txt");
//   for(a=0; a<=_nsorb-1; a++){
//     for(b=0; b<=_nsorb-1; b++){
//       for(c=0; c<=_nsorb-1; c++){
// 	for(d=0; d<=_nsorb-1; d++){
// 	abcd = onei(a,b,c,d,sign);
// 	outFile << sign*RDM2abab[abcd] << std::endl;
// 	}
//       }
//     }
//   }
//   
//   outFile.close();
//   
//    //write RDM2bbbb to a file
//   outFile.open("2RDMbbbb.txt");
//   for(a=0; a<=_nsorb-1; a++){
//     for(b=0; b<=_nsorb-1; b++){
//       for(c=0; c<=_nsorb-1; c++){
// 	for(d=0; d<=_nsorb-1; d++){
// 	abcd = onei(a,b,c,d,sign);
// 	outFile << sign*RDM2bbbb[abcd] << std::endl;
// 	}
//       }
//     }
//   }
//   
//   outFile.close();
//   
//      //write RDM2baab to a file
//   outFile.open("2RDMbaab.txt");
//   for(a=0; a<=_nsorb-1; a++){
//     for(b=0; b<=_nsorb-1; b++){
//       for(c=0; c<=_nsorb-1; c++){
// 	for(d=0; d<=_nsorb-1; d++){
// 	abcd = onei(a,b,c,d,sign);
// 	outFile << sign*RDM2baab[abcd] << std::endl;
// 	}
//       }
//     }
//   }
//   
//   outFile.close();
//   
//     //write RDM2baba to a file
//   outFile.open("2RDMbaba.txt");
//   for(a=0; a<=_nsorb-1; a++){
//     for(b=0; b<=_nsorb-1; b++){
//       for(c=0; c<=_nsorb-1; c++){
// 	for(d=0; d<=_nsorb-1; d++){
// 	abcd = onei(a,b,c,d,sign);
// 	outFile << sign*RDM2baba[abcd] << std::endl;
// 	}
//       }
//     }
//   }
  
//   outFile.close(); 
  //Construct RDM1 for testing
  int sign1;
  uint abdc;
  for(a=0; a<_nsorb; a++){
    for(c=0; c<_nsorb; c++){
      ac = oneid(a,c);
      for(b=0; b<_nsorb; b++){
	d=b;
	abdc = onei(a,b,d,c,sign1);
 	abcd = onei(a,b,c,d,sign);
  	//std::cout << ac << std::endl;
	assert(abcd < nelem);
	assert(abdc < nelem);
	assert(ac < nelem_1d);
	_RDM1[ac] = _RDM1[ac] + sign*_RDM2[abcd]-sign1*_RDM2[abdc];
      }
      _RDM1[ac] = _RDM1[ac] / (2*(nelec-1));
    }
  }



//calculate trace of 1RDM 
Tr = 0.0;
for(a=0;a<_nsorb;a++){
 c=a;
 ac = oneid(a,c);
 Tr = Tr + _RDM1[ac];
}
//std::cout << typeid(Tr).name() << '\n';
std::cout << std::setprecision(16) << std::fixed << Tr << std::endl; 

  store_rdm();
 
}

void DMdump::store_rdm() const
{
  std::ofstream outFile;
  int sign;
  uint abcd,ac;
  //write RDM2 to a file
  outFile.open("2RDM.txt");
  for(uint a=0; a<=_nsorb-1; a++){
    for(uint b=0; b<=_nsorb-1; b++){
      for(uint c=0; c<=_nsorb-1; c++){
        for(uint d=0; d<=_nsorb-1; d++){
          abcd = onei(a,b,c,d,sign);
          outFile << std::setw(15) << sign*_RDM2[abcd] << std::setw(4) << a << std::setw(4) << b << std::setw(4) << c << std::setw(4) << d << std::endl;
        }
      }
    }
  }
  outFile.close();
//    write RDM1 to a file
  outFile.open("1RDM_xx.txt");
  for(uint a=0; a<=_nsorb-1; a++){
    outFile << " " << std::endl;
    for(uint c=0; c <=_nsorb-1; c++){
      ac = oneid(a,c);
      outFile << std::setw(15) << std::left << std::fixed << _RDM1[ac];
    }
  }
  outFile.close(); 

//    write RDM1 to a file
  outFile.open("1RDM.txt");
  for(uint a=0; a<=_nsorb-1; a++){
    for(uint c=0; c <=_nsorb-1; c++){
      ac = oneid(a,c);
      outFile << std::setw(15) << std::left << std::fixed << _RDM1[ac] << std::setw(4) << a << std::setw(4) << c << std::endl;
    }
  }
  outFile.close(); 

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