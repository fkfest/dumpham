#include "odump.h"

Occupation::Occupation(const std::vector< int >& occs, int ibase)
{
  std::vector<int> singly;
  // For 1-based: spatial obital index = int( (spin orbital index + 1) / 2 )
  int iorb0 = -1;
  _foreach_cauto( std::vector< int >, iocc, occs ){
    int iorb1 = int((*iocc-ibase)/2);
    if ( iorb0 == iorb1 ) {
      // doubly occupied
      push_back(iorb0);
      iorb0 = -1;
    } else if ( iorb0 != -1 ) {
      // iorb0 is singly occupied
      singly.push_back(iorb0);
      iorb0 = iorb1;
    } else {
      iorb0 = iorb1;
    }
  }
  // add singly occupied orbitals
  _foreach_cauto( std::vector< int >, iocc, singly ){
    push_back(*iocc);
  }
}

std::ostream& operator<<(std::ostream& o, const Occupation& occ)
{
  _foreach_cauto(Occupation,iocc,occ){
    o << *iocc << ",";
  }
  return o;
}


Odump::Odump(uint norb, const Occupation& occs ) : _nAO(norb), _norb(norb)
{
  this->zero();
  uint iao;
  xout << "Spatial orbital occupation: " << occs << std::endl;
  std::vector<bool> occupied(norb,false);
  for ( uint imo = 0; imo < occs.size(); ++imo ) {
    iao = occs[imo];
    (*this)(iao,imo) = 1.0;
    occupied[iao] = true;
  }
  uint imo = occs.size(); 
  for ( uint iao = 0; iao < norb; ++iao ) {
    if (!occupied[iao]) {
      (*this)(iao,imo) = 1.0;
      ++imo;
    }
  }
  assert( imo == norb );
}

void Odump::store(std::string orbdump)
{
  std::ofstream outputStream;
  outputStream.open(orbdump.c_str());
  if ( (outputStream.rdstate() & std::ifstream::failbit ) != 0 ) {
    outputStream.close();
    error("Odump::store failed to open "+ orbdump);
  }
  int precision = Input::iPars["output"]["precisioncoef"];
  outputStream<<std::scientific<<std::setprecision(precision);
  int maxlen = Input::iPars["output"]["maxncoef"];
  for ( uint i = 0; i < _nAO; ++i ) {
    for ( uint j = 0; j < _norb; ++j ) {
      if ( j > 0 && maxlen > 0 && j % maxlen == 0 ) outputStream << std::endl;
      outputStream << (*this)(i,j) << ", ";
    }
    outputStream << std::endl;
  }
}
