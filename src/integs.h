#ifndef Integs_H
#define Integs_H
#include <vector>
#include "pgsym.h"

#ifdef _DEBUG
#define REDUNWAR_ ,bool redunwar = false
#define REDUNWAR__ ,bool redunwar
#define REDUNWAR ,redunwar
#define WARNRED2(p,q) if(redunwar) warning("Redundant entry in FCIDUMP: " << p << " " << q );
#define WARNRED4(p,q,r,s) if(redunwar) warning("Redundant entry in FCIDUMP: " << p << " " << q << " " << r << " " << s );
#else
#define REDUNWAR_
#define REDUNWAR__
#define REDUNWAR
#define WARNRED2(p,q)
#define WARNRED4(p,q,r,s)
#endif

typedef std::vector<double> DData; 
typedef uint64_t BlkIdx;

/*! 
 * Base class for tensors with permutational and group symmetry
 */
class BaseTensors {
public:
  BaseTensors(uint nidx = 0) : _nidx(nidx) {};
  BaseTensors(PGSym pgs, uint nidx = 0) : _pgs(pgs), _nidx(nidx) {};
  virtual BlkIdx nelem() const { return _data.size(); };
protected:
  PGSym _pgs;
  // number of indices
  uint _nidx;
  // pointers to data for each block (according to irreps)
  // addressed as irrep(p)+irrep(q)*nIrrep+...
  // but actually only symmetry unique blocks will be set! 
  std::vector<BlkIdx> _blocks;
  DData _data;
};

/*! 
 * Class for (pq) with permutational (triangular) and group symmetry
 * (pq)=(qp)
 */
class Integ2 : public BaseTensors {
public:
  Integ2() : BaseTensors(2) {};
  Integ2(PGSym pgs);
  // return (pq) value
  double & value( uint p, uint q REDUNWAR_) { return _data[index(p,q REDUNWAR)];};
  inline BlkIdx index( uint p, uint q REDUNWAR_) const;
};

/*! 
 * Class for (pq|rs) with permutational (triangular) and group symmetry
 * (pq|rs)=(qp|rs)=(pq|sr)=(qp|sr)=(rs|pq)=(sr|pq)=(rs|qp)=(sr|qp)
 */
class Integ4 : public BaseTensors {
public:
  Integ4() : BaseTensors(4) {};
  Integ4(PGSym pgs);
  // return (pq|rs) value
  double & value( uint p, uint q, uint r, uint s REDUNWAR_) { return _data[index(p,q,r,s REDUNWAR)];};
  inline BlkIdx index( uint p, uint q, uint r, uint s REDUNWAR_) const;
};

/*! 
 * Class for (pq|rs) with incomplete permutational (triangular) and group symmetry
 * (pq|rs)=(qp|rs)=(pq|sr)=(qp|sr) 
 * used for (aa|bb) integrals
 */
class Integ4ab : public BaseTensors {
public:
  Integ4ab() : BaseTensors(4) {};
  Integ4ab(PGSym pgs);
  // return (pq|rs) value
  double & value( uint p, uint q, uint r, uint s REDUNWAR_) { return _data[index(p,q,r,s REDUNWAR)];};
  inline BlkIdx index( uint p, uint q, uint r, uint s REDUNWAR_) const;
};

// inline functions
BlkIdx Integ2::index(uint p, uint q REDUNWAR__) const
{
  if ( p < q ) {
    WARNRED2(p,q)
    std::swap(p,q);
  }
  Irrep pir = _pgs.irrep(p),
        qir = _pgs.irrep(q);
  BlkIdx blk_idx = _blocks[pir+qir*_pgs.nIrreps()];
  // indices relative to the block
  uint pb = p - _pgs._firstorb4irrep[pir],
       qb = q - _pgs._firstorb4irrep[qir];
  // triangular index
  return pb*(pb+1)/2+qb+blk_idx;
}

inline BlkIdx Integ4::index(uint p, uint q, uint r, uint s REDUNWAR__) const
{
  if ( p < q ) {
    WARNRED4(p,q,r,s)
    std::swap(p,q);
  }
  if ( r < s ) {
    WARNRED4(p,q,r,s)
    std::swap(r,s);
  }
  if ( p*(p+1)/2+q < r*(r+1)/2+s ){
    WARNRED4(p,q,r,s)
    std::swap(p,r);
    std::swap(q,s);
  }
  Irrep pir = _pgs.irrep(p),
        qir = _pgs.irrep(q),
        rir = _pgs.irrep(r),
        sir = _pgs.irrep(s);
  uint nIrreps = _pgs.nIrreps();
  BlkIdx blk_idx = _blocks[pir+nIrreps*(qir+nIrreps*(rir+nIrreps*sir))];
  // indices relative to the block
  BlkIdx pb = p - _pgs._firstorb4irrep[pir],
         qb = q - _pgs._firstorb4irrep[qir],  
         rb = r - _pgs._firstorb4irrep[rir],  
         sb = s - _pgs._firstorb4irrep[sir];  
  BlkIdx pqb, rsb, lenrs;
  if ( pir == qir ) {
    // triangular index
    pqb = pb*(pb+1)/2+qb,
    rsb = rb*(rb+1)/2+sb;
    lenrs = _pgs._norb4irrep[rir]*(_pgs._norb4irrep[rir]+1)/2;
  } else {
    pqb = pb*_pgs._norb4irrep[qir]+qb,
    rsb = rb*_pgs._norb4irrep[sir]+sb;
    lenrs = _pgs._norb4irrep[rir]*_pgs._norb4irrep[sir];
  }
  if ( pir == rir ) {
    // triangular index
    return pqb*(pqb+1)/2 + rsb + blk_idx;
  } else {
    return pqb*lenrs + rsb + blk_idx;
  }
}

BlkIdx Integ4ab::index(uint p, uint q, uint r, uint s REDUNWAR__) const
{
  if ( p < q ) {
    WARNRED4(p,q,r,s)
    std::swap(p,q);
  }
  if ( r < s ) {
    WARNRED4(p,q,r,s)
    std::swap(r,s);
  }
  Irrep pir = _pgs.irrep(p),
        qir = _pgs.irrep(q),
        rir = _pgs.irrep(r),
        sir = _pgs.irrep(s);
  uint nIrreps = _pgs.nIrreps();
  BlkIdx blk_idx = _blocks[pir+nIrreps*(qir+nIrreps*(rir+nIrreps*sir))];
  // indices relative to the block
  BlkIdx pb = p - _pgs._firstorb4irrep[pir],
         qb = q - _pgs._firstorb4irrep[qir],  
         rb = r - _pgs._firstorb4irrep[rir],  
         sb = s - _pgs._firstorb4irrep[sir];  
  BlkIdx pqb, rsb, lenrs;
  if ( pir == qir ) {
    // triangular index
    pqb = pb*(pb+1)/2+qb,
    rsb = rb*(rb+1)/2+sb;
    lenrs = _pgs._norb4irrep[rir]*(_pgs._norb4irrep[rir]+1)/2;
  } else {
    pqb = pb*_pgs._norb4irrep[qir]+qb,
    rsb = rb*_pgs._norb4irrep[sir]+sb;
    lenrs = _pgs._norb4irrep[rir]*_pgs._norb4irrep[sir];
  }
  return pqb*lenrs + rsb + blk_idx;
}

#undef REDUNWAR_
#undef REDUNWAR__
#undef REDUNWAR
#undef WARNRED2
#undef WARNRED4

#endif

