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

Odump::Odump(uint norb, const Occupation& occs ) : _nAO(norb), _norb(norb)
{
  this->zero();
  uint iao;
  for ( uint imo = 0; imo < norb; ++imo ) {
    if ( imo < occs.size() ) 
      iao = occs[imo];
    else
      iao = imo;
    (*this)(iao,imo) = 1.0;
  }
}

void Odump::store(std::string orbdump)
{
  std::ofstream outputStream;
  outputStream.open(orbdump.c_str());
  if ( (outputStream.rdstate() & std::ifstream::failbit ) != 0 ) {
    outputStream.close();
    error("Odump::store failed to open "+ orbdump);
  }
  outputStream<<std::scientific<<std::setprecision(15);
  int maxlen = Input::iPars["output"]["maxncoef"];
  for ( uint j = 0; j < _norb; ++j ) {
    for ( uint i = 0; i < _nAO; ++i ) {
      if ( i > 0 && maxlen > 0 && i % maxlen == 0 ) outputStream << std::endl;
      outputStream << (*this)(i,j) << ", ";
    }
    outputStream << std::endl;
  }
}
