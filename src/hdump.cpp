#include "hdump.h"

Hdump::Hdump(std::string fcidump) : _dump(fcidump)
{
  xout << "\n Process file "<< fcidump <<std::endl;
  FDPar NELEC = _dump.parameter("NELEC");
  xout << "NELEC=" << NELEC[0] << std::endl;
  FDPar MS2 = _dump.parameter("MS2");
  xout << "MS2=" << MS2[0] << std::endl;
  FDPar NORB = _dump.parameter("NORB");
  xout << "NORB=" << NORB[0] << std::endl;
  FDPar IUHF = _dump.parameter("IUHF");
  xout << "IUHF=" << IUHF[0] << std::endl;
  FDPar ORBSYM = _dump.parameter("ORBSYM");
  xout << "ORBSYM="; 
  _foreach_cauto(FDPar,s,ORBSYM)
    xout << *s << ","; 
  xout<<std::endl;
  
  int i,j,k,l;
  double value;
  int nn = NORB[0];
  if ( nn < 0 ) {
    error("NORB < 0 in FCIDUMP!");
  }
  _norb = nn;
  _nelec = NELEC[0];
  _ms2 = MS2[0];
  int n1el = nn*(nn+1)/2;
  int n2el = n1el*(n1el+1)/2;
  xout << "n2el: " << n2el << std::endl;
  
  _twoel.resize(n2el,0.0);
  _oneel.resize(n1el,0.0);
  _escal = 0.0;
  
  FCIdump::integralType type;
  _dump.rewind();
  
  bool redunwarn = Input::iPars["ham"]["redunwarning"];
  
  while ((type = _dump.nextIntegral(i,j,k,l,value)) != FCIdump::endOfFile) {
//     if ( type == FCIdump::endOfRecord ) {
//         // store
//         for (i = 1; i <= nn; ++i )
//             for (j = 1; j <= i; ++j )
//                 for (k = 1; k <= nn; ++k )
//                     for (l = 1; l <= k; ++l ) {
//                         int ij = (i-1)*i/2 + j;
//                         int kl = (k-1)*k/2 + l;
//                         if ( kl <= ij ){
//                             int ijkl = (ij-1)*ij/2 + kl - 1;
//                             std::cout << twoel[ijkl] << ": " << i << " " << j << " " << k << " " << l << std::endl;
//                         }
//                     }
//     }
    if ( k != 0 && l !=0 ) {
      // Two-electron integrals
      int ij = (i-1)*i/2 + j;
      int kl = (k-1)*k/2 + l;
      if ( j <= i && l <= k && kl <= ij ) {
        int ijkl = (ij-1)*ij/2 + kl - 1;
        _twoel[ijkl] = value;
      } else if ( redunwarn ){
        warning("Redundant entry in " << fcidump << " : " << i << " " << j << " " << k << " " << l );
      }
    } else if ( i != 0 && j != 0 ) {
      // One-electron integrals
      if ( j <= i ) {
        int ij = (i-1)*i/2 + j - 1;
        _oneel[ij] = value;
      } else if ( redunwarn ){
        warning("Redundant entry in " << fcidump << " : " << i << " " << j );
      }
    } else if ( type == FCIdump::I0 ){
      _escal = value;
    }
  }
  //remove the point-group symmetry
  _foreach_auto(FDPar,s,ORBSYM)
    *s = 1;
  _dump.modifyParameter("ORBSYM",ORBSYM);
    
}

void Hdump::store(std::string fcidump)
{
  if (_dump.write(fcidump,FCIdump::FileFormatted,false))
    xout << "will be written to file " << fcidump << std::endl;
  else
    error("failure to write to file "+fcidump);
  
  for ( uint i = 1; i <= _norb; ++i)
    for ( uint j = 1; j <= i; ++j )
      for ( uint k = 1; k <= _norb; ++k )
        for ( uint l = 1; l <= k; ++l ) {
          uint ij = (i-1)*i/2 + j;
          uint kl = (k-1)*k/2 + l;
          if ( kl <= ij ){
            int ijkl = (ij-1)*ij/2 + kl - 1;
            _dump.writeIntegral(i,j,k,l,_twoel[ijkl]);
          }
        }
  for ( uint i = 1; i <= _norb; ++i)
    for ( uint j = 1; j <= i; ++j ) {
      uint ij = (i-1)*i/2 + j - 1;
      _dump.writeIntegral(i,j,0,0,_oneel[ij]);
    }
  _dump.writeIntegral(0,0,0,0,_escal); 
}

uint Hdump::nclosed() const
{
  uint nclos = (_nelec - _ms2)/2;
  if ( nclos*2 + _ms2 != _nelec ) {
    xout << "NELEC: " << _nelec << " MS2: " << _ms2 << std::endl;
    error("Mismatch in NELEC and MS2 in Hdump");
  }
  return nclos;
}
