#ifndef Integs_H
#define Integs_H
#include <vector>
#include "pgsym.h"
#ifdef MOLPRO
#include "hdtypes.h"
#include "hdcommon.h"
#else
#include "globals.h"
#include "utilities.h"
#endif

namespace HamDump {
#ifdef _DEBUG
#define REDUNWAR_ ,bool redunwar = false
#define REDUNWAR__ ,bool redunwar
#define REDUNWAR ,redunwar
#define USERW (void)redunwar
#define WARNRED2(p,q) if(redunwar) warning("Redundant entry in FCIDUMP: "+ std::to_string(p) + " " + std::to_string(q) );
#define WARNRED4(p,q,r,s) if(redunwar) warning("Redundant entry in FCIDUMP: "+std::to_string(p)+" "+std::to_string(q)+" "+std::to_string(r)+" "+std::to_string(s) );
#else
#define REDUNWAR_
#define REDUNWAR__
#define REDUNWAR
#define USERW 
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
  BaseTensors(uint nidx = 0) : _nidx(nidx) {}
  BaseTensors(const PGSym& pgs, uint nidx = 0) : p_pgs(&pgs), _nidx(nidx) {}
  BlkIdx nelem() const { return _data.size(); }
  // set value using the tuple index. Note that in most of the tensors(p,q,..) p is slow running!
  void set( BlkIdx idx, double val ) { assert(idx < _data.size()); _data[idx] = val; }
  // set (pq) value
  void set( uint p, uint q, double val REDUNWAR_) { _data[index(p,q REDUNWAR)] = val; }
  // set (pq|rs) value
  void set( uint p, uint q, uint r, uint s, double val REDUNWAR_) { _data[index(p,q,r,s REDUNWAR)] = val; }
  // get (pq) value
  double get( uint p, uint q REDUNWAR_) const { return _data[index(p,q REDUNWAR)]; }
  // get (pq) value with p,q: indices in irrep ir
  double get( uint p, uint q, Irrep ir REDUNWAR_) const { return _data[index(p,q,ir REDUNWAR)]; }
  // get (pq|rs) value
  double get( uint p, uint q, uint r, uint s REDUNWAR_) const { return _data[index(p,q,r,s REDUNWAR)]; }
  // return (pq) value with point-group symmetry handling
  double get_with_pgs( uint p, uint q ) const {
    if ( p_pgs->totIrrep(p,q) == 0 ) return get(p,q);
    else return 0.0;
  }
  // return (pq|rs) value with point-group symmetry handling
  double get_with_pgs( uint p, uint q, uint r, uint s) const {
    if ( p_pgs->totIrrep(p,q,r,s) == 0 ) return get(p,q,r,s);
    else return 0.0;
  }
  virtual inline BlkIdx index( uint p, uint q REDUNWAR_) const 
          {error("Incompatible index call!","BaseTensors");(void)p;(void)q;USERW;return 0;};
  virtual inline BlkIdx index( uint p, uint q, Irrep ir REDUNWAR_) const 
          {error("Incompatible index call!","BaseTensors");(void)p;(void)q;(void)ir;USERW;return 0;};
  virtual inline BlkIdx index( uint p, uint q, uint r, uint s REDUNWAR_) const 
          { error("Incompatible index call!","BaseTensors");(void)p;(void)q;(void)r;(void)s;USERW;return 0;}
  const PGSym * pgs() const {return p_pgs;}
protected:
  const PGSym * p_pgs;
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
  Integ2(const PGSym& pgs);
  // index of p,q value
  inline BlkIdx index( uint p, uint q REDUNWAR_) const; 
  // index of p,q value with p,q indices in irrep ir
  inline BlkIdx index( uint p, uint q, Irrep ir REDUNWAR_) const; 
};

/*! 
 * Class for (pq) with group symmetry
 */
class Integ2ab : public BaseTensors {
public:
  Integ2ab() : BaseTensors(2) {};
  Integ2ab(const PGSym& pgs);
  // index of p,q value
  inline BlkIdx index( uint p, uint q REDUNWAR_) const; 
  // index of p,q value with p,q indices in irrep ir
  inline BlkIdx index( uint p, uint q, Irrep ir REDUNWAR_) const; 
};

typedef Integ2ab Integ2st;

/*! 
 * Class for (pq|rs) with permutational (triangular) and group symmetry
 * (pq|rs)=(qp|rs)=(pq|sr)=(qp|sr)=(rs|pq)=(sr|pq)=(rs|qp)=(sr|qp)
 */
class Integ4 : public BaseTensors {
public:
  Integ4() : BaseTensors(4) {};
  Integ4(const PGSym& pgs);
  // index of p,q,r,s value
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
  Integ4ab(const PGSym& pgs);
  // index of p,q,r,s value
  inline BlkIdx index( uint p, uint q, uint r, uint s REDUNWAR_) const;
};

/*! 
 * Class for (pq|rs) with incomplete permutational (particle) and group symmetry
 * (pq|rs)=(rs|pq)
 * used for similarity transformed integrals
 */
class Integ4st : public BaseTensors {
public:
  Integ4st() : BaseTensors(4) {};
  Integ4st(const PGSym& pgs);
  // index of p,q,r,s value
  inline BlkIdx index( uint p, uint q, uint r, uint s REDUNWAR_) const;
};

/*! 
 * Class for (pq|rs) with group symmetry
 * (pq|rs) without permutational symmetry
 * used for similarity transformed (aa|bb) integrals
 */
class Integ4stab : public BaseTensors {
public:
  Integ4stab() : BaseTensors(4) {};
  Integ4stab(const PGSym& pgs);
  // index of p,q,r,s value
  inline BlkIdx index( uint p, uint q, uint r, uint s REDUNWAR_) const;
};

// inline functions
inline BlkIdx Integ2::index(uint p, uint q REDUNWAR__) const
{
  if ( p < q ) {
    WARNRED2(p,q)
    std::swap(p,q);
  }
  Irrep pir = p_pgs->irrep(p),
        qir = p_pgs->irrep(q);
  BlkIdx blk_idx = _blocks[pir+qir*p_pgs->nIrreps()];
  // indices relative to the block
  uint pb = p - p_pgs->_firstorb4irrep[pir],
       qb = q - p_pgs->_firstorb4irrep[qir];
  // triangular index
  return pb*(pb+1)/2+qb+blk_idx;
}
inline BlkIdx Integ2::index(uint p, uint q, Irrep ir REDUNWAR__) const
{
  if ( p < q ) {
    WARNRED2(p,q)
    std::swap(p,q);
  }
  // triangular index
  return p*(p+1)/2+q+_blocks[ir+ir*p_pgs->nIrreps()];
}

inline BlkIdx Integ2ab::index(uint p, uint q REDUNWAR__) const
{
  USERW;
  Irrep pir = p_pgs->irrep(p),
        qir = p_pgs->irrep(q);
  BlkIdx blk_idx = _blocks[pir+qir*p_pgs->nIrreps()];
  // indices relative to the block
  uint pb = p - p_pgs->_firstorb4irrep[pir],
       qb = q - p_pgs->_firstorb4irrep[qir];
  return pb*p_pgs->_norb4irrep[qir]+qb+blk_idx;
}
inline BlkIdx Integ2ab::index(uint p, uint q, Irrep ir REDUNWAR__) const
{
  if ( p < q ) {
    WARNRED2(p,q)
    std::swap(p,q);
  }
  // triangular index
  return p*p_pgs->_norb4irrep[ir]+q+_blocks[ir+ir*p_pgs->nIrreps()];
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
  Irrep pir = p_pgs->irrep(p),
        qir = p_pgs->irrep(q),
        rir = p_pgs->irrep(r),
        sir = p_pgs->irrep(s);
  uint nIrreps = p_pgs->nIrreps();
  BlkIdx blk_idx = _blocks[pir+nIrreps*(qir+nIrreps*(rir+nIrreps*sir))];
  // indices relative to the block
  BlkIdx pb = p - p_pgs->_firstorb4irrep[pir],
         qb = q - p_pgs->_firstorb4irrep[qir],  
         rb = r - p_pgs->_firstorb4irrep[rir],  
         sb = s - p_pgs->_firstorb4irrep[sir];  
  BlkIdx pqb, rsb, lenrs;
  if ( pir == qir ) {
    // triangular index
    pqb = pb*(pb+1)/2+qb,
    rsb = rb*(rb+1)/2+sb;
    lenrs = p_pgs->_norb4irrep[rir]*(p_pgs->_norb4irrep[rir]+1)/2;
  } else {
    pqb = pb*p_pgs->_norb4irrep[qir]+qb,
    rsb = rb*p_pgs->_norb4irrep[sir]+sb;
    lenrs = p_pgs->_norb4irrep[rir]*p_pgs->_norb4irrep[sir];
  }
  if ( pir == rir ) {
    // triangular index
    return pqb*(pqb+1)/2 + rsb + blk_idx;
  } else {
    return pqb*lenrs + rsb + blk_idx;
  }
}

inline BlkIdx Integ4ab::index(uint p, uint q, uint r, uint s REDUNWAR__) const
{
  if ( p < q ) {
    WARNRED4(p,q,r,s)
    std::swap(p,q);
  }
  if ( r < s ) {
    WARNRED4(p,q,r,s)
    std::swap(r,s);
  }
  Irrep pir = p_pgs->irrep(p),
        qir = p_pgs->irrep(q),
        rir = p_pgs->irrep(r),
        sir = p_pgs->irrep(s);
  uint nIrreps = p_pgs->nIrreps();
  BlkIdx blk_idx = _blocks[pir+nIrreps*(qir+nIrreps*(rir+nIrreps*sir))];
  // indices relative to the block
  BlkIdx pb = p - p_pgs->_firstorb4irrep[pir],
         qb = q - p_pgs->_firstorb4irrep[qir],  
         rb = r - p_pgs->_firstorb4irrep[rir],  
         sb = s - p_pgs->_firstorb4irrep[sir];  
  BlkIdx pqb, rsb, lenrs;
  if ( pir == qir ) {
    // triangular index
    pqb = pb*(pb+1)/2+qb,
    rsb = rb*(rb+1)/2+sb;
    lenrs = p_pgs->_norb4irrep[rir]*(p_pgs->_norb4irrep[rir]+1)/2;
  } else {
    pqb = pb*p_pgs->_norb4irrep[qir]+qb,
    rsb = rb*p_pgs->_norb4irrep[sir]+sb;
    lenrs = p_pgs->_norb4irrep[rir]*p_pgs->_norb4irrep[sir];
  }
  return pqb*lenrs + rsb + blk_idx;
}

inline BlkIdx Integ4st::index(uint p, uint q, uint r, uint s REDUNWAR__) const
{
  if ( p < r ){
    WARNRED4(p,q,r,s)
    std::swap(p,r);
    std::swap(q,s);
  } else if ( p == r && q < s) {
    WARNRED4(p,q,r,s)
    std::swap(q,s);
  }
  Irrep pir = p_pgs->irrep(p),
        qir = p_pgs->irrep(q),
        rir = p_pgs->irrep(r),
        sir = p_pgs->irrep(s);
  uint nIrreps = p_pgs->nIrreps();
  BlkIdx blk_idx = _blocks[pir+nIrreps*(qir+nIrreps*(rir+nIrreps*sir))];
  // indices relative to the block
  BlkIdx pb = p - p_pgs->_firstorb4irrep[pir],
         qb = q - p_pgs->_firstorb4irrep[qir],  
         rb = r - p_pgs->_firstorb4irrep[rir],  
         sb = s - p_pgs->_firstorb4irrep[sir];  
  BlkIdx pqb, rsb, lenrs;
  pqb = pb*p_pgs->_norb4irrep[qir]+qb,
  rsb = rb*p_pgs->_norb4irrep[sir]+sb;
  if ( pir == rir ) {
    // triangular index
    return pqb*(pqb+1)/2 + rsb + blk_idx;
  } else {
    lenrs = p_pgs->_norb4irrep[rir]*p_pgs->_norb4irrep[sir];
    return pqb*lenrs + rsb + blk_idx;
  }
}

inline BlkIdx Integ4stab::index(uint p, uint q, uint r, uint s REDUNWAR__) const
{
  USERW;
  Irrep pir = p_pgs->irrep(p),
        qir = p_pgs->irrep(q),
        rir = p_pgs->irrep(r),
        sir = p_pgs->irrep(s);
  uint nIrreps = p_pgs->nIrreps();
  BlkIdx blk_idx = _blocks[pir+nIrreps*(qir+nIrreps*(rir+nIrreps*sir))];
  // indices relative to the block
  BlkIdx pb = p - p_pgs->_firstorb4irrep[pir],
         qb = q - p_pgs->_firstorb4irrep[qir],  
         rb = r - p_pgs->_firstorb4irrep[rir],  
         sb = s - p_pgs->_firstorb4irrep[sir];  
  BlkIdx pqb, rsb, lenrs;
  pqb = pb*p_pgs->_norb4irrep[qir]+qb,
  rsb = rb*p_pgs->_norb4irrep[sir]+sb;
  lenrs = p_pgs->_norb4irrep[rir]*p_pgs->_norb4irrep[sir];
  return pqb*lenrs + rsb + blk_idx;
}

#undef REDUNWAR_
#undef REDUNWAR__
#undef REDUNWAR
#undef USERW
#undef WARNRED2
#undef WARNRED4
} // namespace HamDump
#endif

