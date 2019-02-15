#include "integs.h"

namespace HamDump {
Integ2::Integ2(const PGSym& pgs): BaseTensors(pgs,2)
{
  FDPar norb4ir = p_pgs->norbs_in_irreps();
  uint nIrreps = p_pgs->nIrreps();
  _blocks.resize(nIrreps*nIrreps);
  Irrep ir = 0;
  BlkIdx nint = 0;
  BlkIdx BlkLen;
  for ( const auto& ni: norb4ir) {
    // triangular symmetry
    BlkLen = ni * (ni + 1)/2;
    _blocks[ir+ir*nIrreps] = nint;
    ++ir;
    nint += BlkLen;
  }
  _data.resize(nint,0.0);
}

Integ2ab::Integ2ab(const PGSym& pgs): BaseTensors(pgs,2)
{
  FDPar norb4ir = p_pgs->norbs_in_irreps();
  uint nIrreps = p_pgs->nIrreps();
  _blocks.resize(nIrreps*nIrreps);
  Irrep ir = 0;
  BlkIdx nint = 0;
  BlkIdx BlkLen;
  for ( const auto& ni: norb4ir) {
    BlkLen = ni * ni;
    _blocks[ir+ir*nIrreps] = nint;
    ++ir;
    nint += BlkLen;
  }
  _data.resize(nint,0.0);
}


Integ4::Integ4(const PGSym& pgs): BaseTensors(pgs,4)
{
  FDPar norb4ir = p_pgs->norbs_in_irreps();
  uint nIrreps = p_pgs->nIrreps();
  _blocks.resize(nIrreps*nIrreps*nIrreps*nIrreps);
  BlkIdx nint = 0;
  for ( Irrep isym = 0; isym < p_pgs->nIrreps(); ++isym ) {
    std::vector<BlkIdx> len_of2, i_indx, j_indx;
    for ( Irrep iri = 0; iri < p_pgs->nIrreps(); ++iri ) {
      for ( Irrep irj = 0; irj <= iri; ++irj ) {
        if ( p_pgs->product(iri,irj) != isym ) continue;
        if ( iri == irj ) {
          len_of2.push_back(norb4ir[iri]*(norb4ir[iri] + 1)/2);
        } else {
          len_of2.push_back(norb4ir[iri]*norb4ir[irj]);
        }
        i_indx.push_back(iri);
        j_indx.push_back(irj);
      }
    }
    for ( uint iblk = 0; iblk < len_of2.size(); ++iblk ) {
      for ( uint jblk = 0; jblk < len_of2.size(); ++jblk ) {
        // irreps
        uint
          p = i_indx[iblk],
          q = j_indx[iblk],
          r = i_indx[jblk],
          s = j_indx[jblk];
        if ( (p+1)*p/2 + q < (r+1)*r/2 + s ) continue;
        _blocks[p+nIrreps*(q+nIrreps*(r+nIrreps*s))] = nint;
        if ( iblk == jblk ) {
          nint += len_of2[iblk]*(len_of2[iblk]+1)/2;
        } else {
          nint += len_of2[iblk]*len_of2[jblk];
        }
      }
    }
  }
  _data.resize(nint,0.0);
}

Integ4ab::Integ4ab(const PGSym& pgs): BaseTensors(pgs,4)
{
  FDPar norb4ir = p_pgs->norbs_in_irreps();
  uint nIrreps = p_pgs->nIrreps();
  _blocks.resize(nIrreps*nIrreps*nIrreps*nIrreps);
  BlkIdx nint = 0;
  for ( Irrep isym = 0; isym < p_pgs->nIrreps(); ++isym ) {
    std::vector<BlkIdx> len_of2, i_indx, j_indx;
    for ( Irrep iri = 0; iri < p_pgs->nIrreps(); ++iri ) {
      for ( Irrep irj = 0; irj <= iri; ++irj ) {
        if ( p_pgs->product(iri,irj) != isym ) continue;
        if ( iri == irj ) {
          len_of2.push_back(norb4ir[iri]*(norb4ir[iri] + 1)/2);
        } else {
          len_of2.push_back(norb4ir[iri]*norb4ir[irj]);
        }
        i_indx.push_back(iri);
        j_indx.push_back(irj);
      }
    }
    for ( uint iblk = 0; iblk < len_of2.size(); ++iblk ) {
      for ( uint jblk = 0; jblk < len_of2.size(); ++jblk ) {
        // irreps
        uint
          p = i_indx[iblk],
          q = j_indx[iblk],
          r = i_indx[jblk],
          s = j_indx[jblk];
        _blocks[p+nIrreps*(q+nIrreps*(r+nIrreps*s))] = nint;
        nint += len_of2[iblk]*len_of2[jblk];
      }
    }
  }
  _data.resize(nint,0.0);
}

} //namespace HamDump
