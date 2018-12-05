#include "density_mat.h"


DMdump::DMdump(const std::string filename, uint norb, uint nelec)
{
 
  int i, j, k, l;
  int a, b, c, d, ac, abcd;
  int sign;
  int nelem;
  int nelem_1d;
  double S;
  int nsorb;
  nsorb = 2 * norb;
  int maxid = nsorb -1;
  nelem = onei(maxid,maxid,maxid,maxid,sign)+1;
  std::vector<double> RDM2(nelem,0);
  xout << "nelem " << nelem << std::endl;
  nelem_1d = oneid(maxid,maxid)+1;
  xout << "nelem_1d " << nelem_1d << std::endl; 
  std::vector<double> RDM1(nelem_1d,0);
  double Tr;
  int shell;
  shell = 0;
  int test;
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
  
  inFile.open(f_aaaa);
  if (!inFile){
    std::cout << "Unable to open file";
    exit(1);   // call system to stop
  }
  while (nextdm(i,j,k,l,S)){
  a = 2*(i-1);
  b = 2*(j-1);
  c = 2*(k-1);
  d = 2*(l-1);
  abcd = onei(a,b,c,d,sign);
  assert(abcd < nelem);
  RDM2aaaa[abcd] = S * sign;
  if (std::abs(RDM2[abcd]) > 1.e-10) 
    xout << RDM2[abcd] << " " << S * sign << " " << a << " " << b << " " << c << " " << d << std::endl;
  RDM2[abcd] = S * sign;
  
  }
  
  inFile.open(f_abab);
  if (!inFile){
    std::cout << "Unable to open file";
    exit(1);   // call system to stop
  }
  while (nextdm(i,j,k,l,S)){
  a = 2*(i-1);
  b = 2*j-1;
  c = 2*(k-1);
  d = 2*l-1;
  abcd = onei(a,b,c,d,sign);
  assert(abcd < nelem);
  RDM2abab[abcd] = S * sign;
  if (std::abs(RDM2[abcd]) > 1.e-10) 
    xout << RDM2[abcd] << " " << S * sign << " " << a << " " << b << " " << c << " " << d << std::endl;
  RDM2[abcd] = S * sign;
  }
  
  inFile.open(f_abba);
  if (!inFile){
    std::cout << "Unable to open file";
    exit(1);   // call system to stop
  }
  while (nextdm(i,j,k,l,S)){
  a = 2*(i-1);
  b = 2*j-1;
  c = 2*k-1;
  d = 2*(l-1);
  abcd = onei(a,b,c,d,sign);
  assert(abcd < nelem);
  RDM2abba[abcd] = S * sign;
  if (std::abs(RDM2[abcd]) > 1.e-10) 
    xout << RDM2[abcd] << " " << S * sign << " " << a << " " << b << " " << c << " " << d << std::endl;
  RDM2[abcd] = S * sign;
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
  while (nextdm(i,j,k,l,S)){
  a = 2*i-1;
  b = 2*j-1;
  c = 2*k-1;
  d = 2*l-1;
  abcd = onei(a,b,c,d,sign);
  assert(abcd < nelem);
  RDM2bbbb[abcd] = S * sign;
  if (std::abs(RDM2[abcd]) > 1.e-10) 
    xout << RDM2[abcd] << " " << S * sign << " " << a << " " << b << " " << c << " " << d << std::endl;
  RDM2[abcd] = S * sign;
  }
  
  inFile.open(f_baba);
  if (!inFile){
    std::cout << "Unable to open file";
    exit(1);   // call system to stop
  }
  while (nextdm(i,j,k,l,S)){
  a = 2*i-1;
  b = 2*(j-1);
  c = 2*k-1;
  d = 2*(l-1);
  abcd = onei(a,b,c,d,sign);
  assert(abcd < nelem);
  RDM2baba[abcd] = S * sign;
  if (std::abs(RDM2[abcd]) > 1.e-10) 
    xout << RDM2[abcd] << " " << S * sign << " " << a << " " << b << " " << c << " " << d << std::endl;
  RDM2[abcd] = S * sign;
  }
  
  inFile.open(f_baab);
  if (!inFile){
    std::cout << "Unable to open file";
    exit(1);   // call system to stop
  }
  while (nextdm(i,j,k,l,S)){
  a = 2*i-1;
  b = 2*(j-1);
  c = 2*(k-1);
  d = 2*l-1;
  abcd = onei(a,b,c,d,sign);
  assert(abcd < nelem);
  RDM2baab[abcd] = S * sign;
  if (std::abs(RDM2[abcd]) > 1.e-10) 
    xout << RDM2[abcd] << " " << S * sign << " " << a << " " << b << " " << c << " " << d << std::endl;
  RDM2[abcd] = S * sign;
  }

  //write RDM2 to a file
  outFile.open("2RDM.txt");
  for(a=0; a<=nsorb-1; a++){
    for(b=0; b<=nsorb-1; b++){
      for(c=0; c<=nsorb-1; c++){
	for(d=0; d<=nsorb-1; d++){
	abcd = onei(a,b,c,d,sign);
	outFile << std::setw(15) << sign*RDM2[abcd] << std::setw(4) << a << std::setw(4) << b << std::setw(4) << c << std::setw(4) << d << std::endl;
	}
      }
    }
  }
  outFile.close();
  
//   //write RDM2aaaa to a file
//   outFile.open("2RDMaaaa.txt");
//   for(a=0; a<=nsorb-1; a++){
//     for(b=0; b<=nsorb-1; b++){
//       for(c=0; c<=nsorb-1; c++){
// 	for(d=0; d<=nsorb-1; d++){
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
//   for(a=0; a<=nsorb-1; a++){
//     for(b=0; b<=nsorb-1; b++){
//       for(c=0; c<=nsorb-1; c++){
// 	for(d=0; d<=nsorb-1; d++){
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
//   for(a=0; a<=nsorb-1; a++){
//     for(b=0; b<=nsorb-1; b++){
//       for(c=0; c<=nsorb-1; c++){
// 	for(d=0; d<=nsorb-1; d++){
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
//   for(a=0; a<=nsorb-1; a++){
//     for(b=0; b<=nsorb-1; b++){
//       for(c=0; c<=nsorb-1; c++){
// 	for(d=0; d<=nsorb-1; d++){
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
//   for(a=0; a<=nsorb-1; a++){
//     for(b=0; b<=nsorb-1; b++){
//       for(c=0; c<=nsorb-1; c++){
// 	for(d=0; d<=nsorb-1; d++){
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
//   for(a=0; a<=nsorb-1; a++){
//     for(b=0; b<=nsorb-1; b++){
//       for(c=0; c<=nsorb-1; c++){
// 	for(d=0; d<=nsorb-1; d++){
// 	abcd = onei(a,b,c,d,sign);
// 	outFile << sign*RDM2baba[abcd] << std::endl;
// 	}
//       }
//     }
//   }
  
//   outFile.close(); 
  //Construct RDM1 for testing
  int sign1;
  int abdc;
  for(a=0; a<nsorb; a++){
    for(c=0; c<nsorb; c++){
      ac = oneid(a,c);
      for(b=0; b<nsorb; b++){
	d=b;
	abdc = onei(a,b,d,c,sign1);
 	abcd = onei(a,b,c,d,sign);
  	//std::cout << ac << std::endl;
	assert(abcd < nelem);
	assert(abdc < nelem);
	assert(ac < nelem_1d);
	RDM1[ac] = RDM1[ac] + sign*RDM2[abcd]-sign1*RDM2[abdc];
      }
      RDM1[ac] = RDM1[ac] / (2*(nelec-1));
    }
  }

//    write RDM1 to a file
  outFile.open("1RDM_xx.txt");
  for(a=0; a<=nsorb-1; a++){
      outFile << " " << std::endl;
      for(c=0; c <=nsorb-1; c++){
	ac = oneid(a,c);
	outFile << std::setw(15) << std::left << std::fixed << RDM1[ac];
      }
  }
  outFile.close(); 

//    write RDM1 to a file
  outFile.open("1RDM.txt");
  for(a=0; a<=nsorb-1; a++){
      for(c=0; c <=nsorb-1; c++){
	ac = oneid(a,c);
	outFile << std::setw(15) << std::left << std::fixed << RDM1[ac] << std::setw(4) << a << std::setw(4) << c << std::endl;
      }
  }
  outFile.close(); 


//calculate trace of 1RDM 
Tr = 0.0;
for(a=0;a<nsorb;a++){
 c=a;
 ac = oneid(a,c);
 Tr = Tr + RDM1[ac];
}
//std::cout << typeid(Tr).name() << '\n';
std::cout << std::setprecision(16) << std::fixed << Tr << std::endl; 
 
}
 
bool DMdump::nextdm(int& i, int& j, int& k, int& l, double& S)
 {
  if(inFile >> i >> j >> k >> l >> S){
   return true;
  }
  else{
  inFile.close();
  return false;  
  }
 }
 
int DMdump::onei(int a, int b, int c, int d, int& sign)
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

int DMdump::oneid(int a, int c)
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