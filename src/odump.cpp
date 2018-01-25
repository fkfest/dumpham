#include "odump.h"

Odump::Odump(uint norb) : _nAO(norb), _norb(norb)
{
  this->zero();
  for ( uint i = 0; i < norb; ++i )
    (*this)(i,i) = 1.0;
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
