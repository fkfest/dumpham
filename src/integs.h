#ifndef Integs_H
#define Integs_H
#include <vector>
#include "pgsym.h"
#ifdef MOLPRO
#include "hdtypes.h"
#include "hdcommon.h"
#include "memory.h"
#else
#include "globals.h"
#include "utilities.h"
#endif

namespace HamDump {
#ifdef _DEBUG
#define WARNRED2(p,q) if(redunwarn) warning("Redundant entry in FCIDUMP: "+ std::to_string(p) + " " + std::to_string(q) );
#define WARNRED4(p,q,r,s) if(redunwarn) warning("Redundant entry in FCIDUMP: "+std::to_string(p)+" "+std::to_string(q)+" "+std::to_string(r)+" "+std::to_string(s) );
#else
#define WARNRED2(p,q)
#define WARNRED4(p,q,r,s)
#endif

#ifdef MOLPRO
typedef memory::vector<double> DData;
#else
typedef std::vector<double> DData; 
#endif
typedef uint64_t BlkIdx;

/*! 
 * Base class for tensors with permutational and group symmetry
 */
class BaseTensors {
public:
  BaseTensors(uint nidx = 0) : _nidx(nidx) {}
  BaseTensors(const PGSym& pgs, uint nidx = 0) : p_pgs(&pgs), _nidx(nidx) {}
  virtual ~BaseTensors() = default;
  BlkIdx nelem() const { return _data.size(); }
  // set value using the tuple index. Note that in most of the tensors(p,q,..) p is slow running!
  void set( BlkIdx idx, double val ) { assert(idx < _data.size()); _data[idx] = val; }
  // set (pq) value
  void set( uint p, uint q, double val ) { _data[index(p,q)] = val; }
  // set (pq) value with p,q: indices in irrep ir
  void set( uint p, uint q, Irrep ir, double val ) { _data[index(p,q,ir)] = val; }
  // set (pq|rs) value
  void set( uint p, uint q, uint r, uint s, double val ) { _data[index(p,q,r,s)] = val; }
  // set value using the tuple index. Note that in most of the tensors(p,q,..) p is slow running!
  double get( BlkIdx idx ) { assert(idx < _data.size()); return _data[idx]; }
  // get (pq) value
  double get( uint p, uint q ) const { return _data[index(p,q)]; }
  // get (pq) value with p,q: indices in irrep ir
  double get( uint p, uint q, Irrep ir ) const { return _data[index(p,q,ir)]; }
  // get (pq|rs) value
  double get( uint p, uint q, uint r, uint s ) const { return _data[index(p,q,r,s)]; }
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
  virtual inline BlkIdx index( uint p, uint q ) const 
          {error("Incompatible index call!","BaseTensors");(void)p;(void)q;return 0;};
  virtual inline BlkIdx index( uint p, uint q, Irrep ir ) const 
          {error("Incompatible index call!","BaseTensors");(void)p;(void)q;(void)ir;return 0;};
  virtual inline BlkIdx index( uint p, uint q, uint r, uint s ) const 
          { error("Incompatible index call!","BaseTensors");(void)p;(void)q;(void)r;(void)s;return 0;}
  virtual inline bool next_indices( uint& p, uint& q ) const 
          {error("Incompatible next_indices call!","BaseTensors");(void)p;(void)q;return false;};
  virtual inline bool next_indices_nosym( uint& p, uint& q ) const 
          {error("Incompatible next_indices call!","BaseTensors");(void)p;(void)q;return false;};
  virtual inline bool next_indices( uint& p, uint& q, uint& r, uint& s, Irrep& ir ) const 
          {error("Incompatible next_indices call!","BaseTensors");(void)p;(void)q;(void)r;(void)s;(void)ir;return false;};
  virtual inline bool next_indices_nosym( uint& p, uint& q, uint& r, uint& s ) const 
          {error("Incompatible next_indices call!","BaseTensors");(void)p;(void)q;(void)r;(void)s;return false;};
          
  const PGSym * pgs() const {return p_pgs;}
#ifdef _DEBUG
  bool redunwarn = false;
#endif
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
  inline BlkIdx index( uint p, uint q ) const; 
  // index of p,q value with p,q indices in irrep ir
  inline BlkIdx index( uint p, uint q, Irrep ir ) const; 
  // iterate to next index
  inline bool next_indices( uint& p, uint& q ) const;
  // iterate to next index without symmetry
  inline bool next_indices_nosym( uint& p, uint& q ) const;
};

/*! 
 * Class for (pq) with group symmetry
 */
class Integ2ab : public BaseTensors {
public:
  Integ2ab() : BaseTensors(2) {};
  Integ2ab(const PGSym& pgs);
  // index of p,q value
  inline BlkIdx index( uint p, uint q ) const; 
  // index of p,q value with p,q indices in irrep ir
  inline BlkIdx index( uint p, uint q, Irrep ir ) const; 
  // iterate to next index
  inline bool next_indices( uint& p, uint& q ) const;
  // iterate to next index without symmetry
  inline bool next_indices_nosym( uint& p, uint& q ) const;
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
  inline BlkIdx index( uint p, uint q, uint r, uint s ) const;
  // iterate to next index
  inline bool next_indices( uint& p, uint& q, uint& r, uint& s, Irrep& ir ) const;
  // iterate to next index without symmetry
  inline bool next_indices_nosym( uint& p, uint& q, uint& r, uint& s ) const;
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
  inline BlkIdx index( uint p, uint q, uint r, uint s ) const;
  // iterate to next index
  inline bool next_indices( uint& p, uint& q, uint& r, uint& s, Irrep& ir ) const;
  // iterate to next index without symmetry
  inline bool next_indices_nosym( uint& p, uint& q, uint& r, uint& s ) const;
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
  inline BlkIdx index( uint p, uint q, uint r, uint s ) const;
  // iterate to next index
  inline bool next_indices( uint& p, uint& q, uint& r, uint& s, Irrep& ir ) const;
  // iterate to next index without symmetry
  inline bool next_indices_nosym( uint& p, uint& q, uint& r, uint& s ) const;
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
  inline BlkIdx index( uint p, uint q, uint r, uint s ) const;
  // iterate to next index
  inline bool next_indices( uint& p, uint& q, uint& r, uint& s, Irrep& ir ) const;
  // iterate to next index without symmetry
  inline bool next_indices_nosym( uint& p, uint& q, uint& r, uint& s ) const;
};

// inline functions
inline BlkIdx Integ2::index(uint p, uint q ) const
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
inline BlkIdx Integ2::index(uint p, uint q, Irrep ir ) const
{
  if ( p < q ) {
    WARNRED2(p,q)
    std::swap(p,q);
  }
  // triangular index
  return p*(p+1)/2+q+_blocks[ir+ir*p_pgs->nIrreps()];
}
inline bool Integ2::next_indices(uint& p, uint& q) const
{
  ++q;
  do {
    for ( ; q <= p; ++q ) {
      if ( p_pgs->totIrrep(p,q) == 0 ) return true;
    } q = 0;
  ++p; } while ( p < p_pgs->ntotorbs() );
  return false;
}
bool Integ2::next_indices_nosym(uint& p, uint& q) const
{
  ++q;
  do {
    for ( ; q <= p; ++q ) {
        return true;
    } q = 0;
  ++p; } while ( p < p_pgs->ntotorbs());
  return false;
}

inline BlkIdx Integ2ab::index(uint p, uint q ) const
{
  Irrep pir = p_pgs->irrep(p),
        qir = p_pgs->irrep(q);
  BlkIdx blk_idx = _blocks[pir+qir*p_pgs->nIrreps()];
  // indices relative to the block
  uint pb = p - p_pgs->_firstorb4irrep[pir],
       qb = q - p_pgs->_firstorb4irrep[qir];
  return pb*p_pgs->_norb4irrep[qir]+qb+blk_idx;
}
inline BlkIdx Integ2ab::index(uint p, uint q, Irrep ir ) const
{
  return p*p_pgs->_norb4irrep[ir]+q+_blocks[ir+ir*p_pgs->nIrreps()];
}
inline bool Integ2ab::next_indices(uint& p, uint& q) const
{
  ++q;
  do {
    for ( ; q < p_pgs->ntotorbs(); ++q ) {
        if ( p_pgs->totIrrep(p,q) == 0 ) return true;
    } q = 0;
  ++p; } while ( p < p_pgs->ntotorbs() );
  return false;
}
bool Integ2ab::next_indices_nosym(uint& p, uint& q) const
{
  ++q;
  do {
    for ( ; q < p_pgs->ntotorbs(); ++q ) {
        return true;
    } q = 0;
  ++p; } while ( p < p_pgs->ntotorbs());
  return false;
}

inline BlkIdx Integ4::index(uint p, uint q, uint r, uint s ) const
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
inline bool Integ4::next_indices(uint& p, uint& q, uint& r, uint& s, Irrep& ir) const
{
  ++s;
  do { // ir
    do { // p
      do { // q
        if ( p_pgs->totIrrep(p,q) == ir ) do { // r
          for ( ; s <= r; ++s ) {
            if ( p_pgs->totIrrep(r,s) == ir && (r+1)*r/2 + s <= (p+1)*p/2 + q ) return true;
          } s = 0;
        ++r; } while ( r <= p ); r = 0; 
      ++q; } while ( q <= p ); q = 0;
    ++p; } while ( p < p_pgs->ntotorbs() ); p = 0;
  ++ir; } while ( ir < p_pgs->nIrreps() );
  return false;
}
inline bool Integ4::next_indices_nosym(uint& p, uint& q, uint& r, uint& s) const
{
  ++s;
  do { // p
    do { // q
      do { // r
        for ( ; s <= r; ++s ) {
          if ( (r+1)*r/2 + s <= (p+1)*p/2 + q ) return true;
        } s = 0;
      ++r; } while ( r <= p ); r = 0;
    ++q; } while ( q <= p ); q = 0;
  ++p; } while ( p < p_pgs->ntotorbs());
  return false;
}

inline BlkIdx Integ4ab::index(uint p, uint q, uint r, uint s ) const
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
inline bool Integ4ab::next_indices(uint& p, uint& q, uint& r, uint& s, Irrep& ir) const
{
  ++s;
  do { // ir
    do { // p
      do { // q
        if ( p_pgs->totIrrep(p,q) == ir ) do { // r
          for ( ; s <= r; ++s ) {
            if ( p_pgs->totIrrep(r,s) == ir ) return true;
          } s = 0;
        ++r; } while ( r < p_pgs->ntotorbs() ); r = 0; 
      ++q; } while ( q <= p ); q = 0;
    ++p; } while ( p < p_pgs->ntotorbs() ); p = 0;
  ++ir; } while ( ir < p_pgs->nIrreps() );
  return false;
}
inline bool Integ4ab::next_indices_nosym(uint& p, uint& q, uint& r, uint& s) const
{
  ++s;
  do { // p
    do { // q
      do { // r
        if( s <= r ) return true;
        s = 0;
      ++r; } while ( r < p_pgs->ntotorbs() ); r = 0;
    ++q; } while ( q <= p ); q = 0;
  ++p; } while ( p < p_pgs->ntotorbs() );
  
  return false;
}

inline BlkIdx Integ4st::index(uint p, uint q, uint r, uint s ) const
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
inline bool Integ4st::next_indices(uint& p, uint& q, uint& r, uint& s, Irrep& ir) const
{
  ++s;
  do { // ir
    do { // p
      do { // q
        if ( p_pgs->totIrrep(p,q) == ir ) do { // r 
          if ( p == r ) {
            for ( ; s <= q; ++s ) {
              if ( p_pgs->totIrrep(r,s) == ir ) return true;
            } 
          } else {
            for ( ; s < p_pgs->ntotorbs(); ++s ) {
              if ( p_pgs->totIrrep(r,s) == ir ) return true;
            } 
          } 
          s = 0;
        ++r; } while ( r <= p ); r = 0; 
      ++q; } while ( q < p_pgs->ntotorbs() ); q = 0;
    ++p; } while ( p < p_pgs->ntotorbs() ); p = 0;
  ++ir; } while ( ir < p_pgs->nIrreps() );
  
  return false;
}
inline bool Integ4st::next_indices_nosym(uint& p, uint& q, uint& r, uint& s) const
{
  ++s;
  do { // p
    do { // q
      do { // r
        if ( p == r ) {
          if ( s <= q ) return true;
        } else {
          if ( s < p_pgs->ntotorbs() ) return true;
        }
        s = 0;
      ++r; } while ( r <= p ); r = 0;
    ++q; } while ( q < p_pgs->ntotorbs() ); q = 0;
  ++p; } while ( p < p_pgs->ntotorbs() ); 
  
  return false;
}

inline BlkIdx Integ4stab::index(uint p, uint q, uint r, uint s ) const
{
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
inline bool Integ4stab::next_indices(uint& p, uint& q, uint& r, uint& s, Irrep& ir) const
{
  ++s;
  do { // ir
    do { // p
      do { // q
        if ( p_pgs->totIrrep(p,q) == ir ) do { // r
          for ( ; s < p_pgs->ntotorbs(); ++s ) {
            if ( p_pgs->totIrrep(r,s) == ir ) return true;
          } s = 0;
        ++r; } while ( r < p_pgs->ntotorbs() ); r = 0; 
      ++q; } while ( q < p_pgs->ntotorbs() ); q = 0;
    ++p; } while ( p < p_pgs->ntotorbs() ); p = 0;
  ++ir; } while ( ir < p_pgs->nIrreps() );
  
  return false;
}
inline bool Integ4stab::next_indices_nosym(uint& p, uint& q, uint& r, uint& s) const
{
  ++s;
  do { // p
    do { // q
      do { // r
        if ( s < p_pgs->ntotorbs() ) return true;
        s = 0;
      ++r; } while ( r < p_pgs->ntotorbs() ); r = 0;
    ++q; } while ( q < p_pgs->ntotorbs() ); q = 0;
  ++p; } while ( p < p_pgs->ntotorbs() );
  
  return false;
}

#undef WARNRED2
#undef WARNRED4
} // namespace HamDump
#endif

