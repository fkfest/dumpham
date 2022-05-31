#ifndef Integs_H
#define Integs_H
#include <vector>
#include "pgsym.h"
#ifdef MOLPRO
#include "hdtypes.h"
#include "hdcommon.h"
#include "molpro/memory.h"
#else
#include "globals.h"
#include "utilities.h"
#endif

namespace HamDump {
#ifdef _DEBUG
#define WARNRED2(p,q) if(redunwarn) warning("Redundant entry in FCIDUMP: "+ std::to_string(p) + " " + std::to_string(q) );
#define WARNRED4(p,q,r,s) if(redunwarn) warning("Redundant entry in FCIDUMP: "+std::to_string(p)+" "+std::to_string(q)+" "+std::to_string(r)+" "+std::to_string(s) );
#define WARNRED6(p,q,r,s,t,u) if(redunwarn) warning("Redundant entry in FCIDUMP: "+std::to_string(p)+" "+std::to_string(q)+" "+std::to_string(r)+" "+std::to_string(s)+" "+std::to_string(t)+" "+std::to_string(u) );
#else
#define WARNRED2(p,q)
#define WARNRED4(p,q,r,s)
#define WARNRED6(p,q,r,s,t,u)
#endif

#ifdef MOLPRO
typedef molpro::vector<double> DData;
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
  // set (pq|rs|tu) value
  void set( uint p, uint q, uint r, uint s, uint t, uint u, double val ) { _data[index(p,q,r,s,t,u)] = val; }
  // set value using the tuple index. Note that in most of the tensors(p,q,..) p is slow running!
  double get( BlkIdx idx ) { assert(idx < _data.size()); return _data[idx]; }
  // get (pq) value
  double get( uint p, uint q ) const { return _data[index(p,q)]; }
  // get (pq) value with p,q: indices in irrep ir
  double get( uint p, uint q, Irrep ir ) const { return _data[index(p,q,ir)]; }
  // get (pq|rs) value
  double get( uint p, uint q, uint r, uint s ) const { return _data[index(p,q,r,s)]; }
  // get (pq|rs|tu) value
  double get( uint p, uint q, uint r, uint s, uint t, uint u ) const { return _data[index(p,q,r,s,t,u)]; }
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
  // return (pq|rs|tu) value with point-group symmetry handling
  double get_with_pgs( uint p, uint q, uint r, uint s, uint t, uint u) const {
    if ( p_pgs->totIrrep(p,q,r,s,t,u) == 0 ) return get(p,q,r,s,t,u);
    else return 0.0;
  }
  virtual inline BlkIdx index( uint p, uint q ) const
          {error("Incompatible index call!","BaseTensors");(void)p;(void)q;return 0;};
  virtual inline BlkIdx index( uint p, uint q, Irrep ir ) const
          {error("Incompatible index call!","BaseTensors");(void)p;(void)q;(void)ir;return 0;};
  virtual inline BlkIdx index( uint p, uint q, uint r, uint s ) const
          { error("Incompatible index call!","BaseTensors");(void)p;(void)q;(void)r;(void)s;return 0;}
  virtual inline BlkIdx index( uint p, uint q, uint r, uint s, uint t, uint u ) const
          { error("Incompatible index call!","BaseTensors");(void)p;(void)q;(void)r;(void)s;(void)t;(void)u;return 0;}
  virtual inline bool next_indices( uint& p, uint& q ) const
          {error("Incompatible next_indices call!","BaseTensors");(void)p;(void)q;return false;};
  virtual inline bool next_indices_nosym( uint& p, uint& q ) const
          {error("Incompatible next_indices call!","BaseTensors");(void)p;(void)q;return false;};
  virtual inline bool next_indices( uint& p, uint& q, uint& r, uint& s, Irrep& ir ) const
          {error("Incompatible next_indices call!","BaseTensors");(void)p;(void)q;(void)r;(void)s;(void)ir;return false;};
  virtual inline bool next_indices_nosym( uint& p, uint& q, uint& r, uint& s ) const
          {error("Incompatible next_indices call!","BaseTensors");(void)p;(void)q;(void)r;(void)s;return false;};
  virtual inline bool next_indices( uint& p, uint& q, uint& r, uint& s, uint& t, uint& u, Irrep& ir, Irrep& ir2 ) const
          {error("Incompatible next_indices call!","BaseTensors");(void)p;(void)q;(void)r;(void)s;(void)t;(void)u;(void)ir;(void)ir2;return false;};
  virtual inline bool next_indices_nosym( uint& p, uint& q, uint& r, uint& s, uint& t, uint& u) const
          {error("Incompatible next_indices call!","BaseTensors");(void)p;(void)q;(void)r;(void)s;(void)t;(void)u;return false;};

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

/*!
 * Class for (pq|rs|tu) with permutational (triangular) and group symmetry
 * (pq|rs|tu)=(qp|rs|tu)=(pq|sr|tu)=(pq|rs|ut)=(qp|sr|tu)=....=(rs|pq|tu)=...
 */
class Integ6 : public BaseTensors {
public:
  Integ6() : BaseTensors(6) {};
  Integ6(const PGSym& pgs);
  // index of p,q,r,s,t,u value
  inline BlkIdx index( uint p, uint q, uint r, uint s, uint t, uint u ) const;
  // iterate to next index
  inline bool next_indices( uint& p, uint& q, uint& r, uint& s, uint& t, uint& u,
                            Irrep& irpq, Irrep& irrs ) const;
  // iterate to next index without symmetry
  inline bool next_indices_nosym( uint& p, uint& q, uint& r, uint& s, uint& t, uint& u ) const;
};

/*!
 * Class for (pq|rs|tu) with permutational (triangular) symmetry for two last electrons and group symmetry
 * (pq|rs|tu)=(qp|rs|tu)=(pq|sr|tu)=(pq|rs|ut)=(qp|sr|tu)=(pq|tu|rs)=...
 */
class Integ6abb : public BaseTensors {
public:
  Integ6abb() : BaseTensors(6) {};
  Integ6abb(const PGSym& pgs);
  // index of p,q,r,s,t,u value
  inline BlkIdx index( uint p, uint q, uint r, uint s, uint t, uint u ) const;
  // iterate to next index
  inline bool next_indices( uint& p, uint& q, uint& r, uint& s, uint& t, uint& u,
                            Irrep& irpq, Irrep& irrs ) const;
  // iterate to next index without symmetry
  inline bool next_indices_nosym( uint& p, uint& q, uint& r, uint& s, uint& t, uint& u ) const;
};

/*!
 * Class for (pq|rs|tu) with incomplete permutational (particle) and group symmetry
 * (pq|rs|tu)=(pq|tu|rs)=(rs|pq|tu)=(rs|tu|pq)=(tu|pq|rs)=(tu|rs|pq)
 */
class Integ6st : public BaseTensors {
public:
  Integ6st() : BaseTensors(6) {};
  Integ6st(const PGSym& pgs);
  // index of p,q,r,s,t,u value
  inline BlkIdx index( uint p, uint q, uint r, uint s, uint t, uint u ) const;
  // iterate to next index
  inline bool next_indices( uint& p, uint& q, uint& r, uint& s, uint& t, uint& u,
                            Irrep& irpq, Irrep& irrs ) const;
  // iterate to next index without symmetry
  inline bool next_indices_nosym( uint& p, uint& q, uint& r, uint& s, uint& t, uint& u ) const;
};

/*!
 * Class for (pq|rs|tu) with permutational (particle) symmetry for two last electrons and group symmetry
 * (pq|rs|tu)=(pq|tu|rs)
 */
class Integ6stabb : public BaseTensors {
public:
  Integ6stabb() : BaseTensors(6) {};
  Integ6stabb(const PGSym& pgs);
  // index of p,q,r,s,t,u value
  inline BlkIdx index( uint p, uint q, uint r, uint s, uint t, uint u ) const;
  // iterate to next index
  inline bool next_indices( uint& p, uint& q, uint& r, uint& s, uint& t, uint& u,
                            Irrep& irpq, Irrep& irrs ) const;
  // iterate to next index without symmetry
  inline bool next_indices_nosym( uint& p, uint& q, uint& r, uint& s, uint& t, uint& u ) const;
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

inline BlkIdx Integ6::index(uint p, uint q, uint r, uint s, uint t, uint u ) const
{
  if ( p < q ) {
    WARNRED6(p,q,r,s,t,u)
    std::swap(p,q);
  }
  if ( r < s ) {
    WARNRED6(p,q,r,s,t,u)
    std::swap(r,s);
  }
  if ( t < u ) {
    WARNRED6(p,q,r,s,t,u)
    std::swap(t,u);
  }
  if ( p*(p+1)/2+q < r*(r+1)/2+s ){
    WARNRED6(p,q,r,s,t,u)
    std::swap(p,r);
    std::swap(q,s);
  }
  if ( r*(r+1)/2+s < t*(t+1)/2+u ){
    WARNRED6(p,q,r,s,t,u)
    std::swap(r,t);
    std::swap(s,u);
    if ( p*(p+1)/2+q < r*(r+1)/2+s ){
      std::swap(p,r);
      std::swap(q,s);
    }
  }

  Irrep pir = p_pgs->irrep(p),
        qir = p_pgs->irrep(q),
        rir = p_pgs->irrep(r),
        sir = p_pgs->irrep(s),
        tir = p_pgs->irrep(t),
        uir = p_pgs->irrep(u);
  uint nIrreps = p_pgs->nIrreps();
  BlkIdx blk_idx = _blocks[pir+nIrreps*(qir+nIrreps*(rir+nIrreps*(sir+nIrreps*(tir+nIrreps*uir))))];
  // indices relative to the block
  BlkIdx pb = p - p_pgs->_firstorb4irrep[pir],
         qb = q - p_pgs->_firstorb4irrep[qir],
         rb = r - p_pgs->_firstorb4irrep[rir],
         sb = s - p_pgs->_firstorb4irrep[sir],
         tb = t - p_pgs->_firstorb4irrep[tir],
         ub = u - p_pgs->_firstorb4irrep[uir];
  BlkIdx pqb, rsb, tub, lenrs, lentu;
  if ( pir == qir ) {
    // triangular index
    pqb = pb*(pb+1)/2+qb;
  } else {
    pqb = pb*p_pgs->_norb4irrep[qir]+qb;
  }
  if ( rir == sir ) {
    // triangular index
    rsb = rb*(rb+1)/2+sb;
    lenrs = p_pgs->_norb4irrep[rir]*(p_pgs->_norb4irrep[rir]+1)/2;
  } else {
    rsb = rb*p_pgs->_norb4irrep[sir]+sb;
    lenrs = p_pgs->_norb4irrep[rir]*p_pgs->_norb4irrep[sir];
  }
  if ( tir == uir ) {
    // triangular index
    tub = tb*(tb+1)/2+ub;
    lentu = p_pgs->_norb4irrep[tir]*(p_pgs->_norb4irrep[tir]+1)/2;
  } else {
    tub = tb*p_pgs->_norb4irrep[uir]+ub;
    lentu = p_pgs->_norb4irrep[tir]*p_pgs->_norb4irrep[uir];
  }
  if ( pir == rir ) {
    // triangular index
    if ( pir == tir ) {
      return pqb*(pqb+1)*(pqb+2)/6 + rsb*(rsb+1)/2 + tub + blk_idx;
    } else {
      return (pqb*(pqb+1)/2 + rsb)*lentu + tub + blk_idx;
    }
  } else if ( rir == tir ) {
      return pqb*(lenrs*(lenrs-1)/2+lentu) + rsb*(rsb+1)/2 + tub + blk_idx;
  } else {
    return (pqb*lenrs + rsb)*lentu + tub + blk_idx;
  }
}
inline bool Integ6::next_indices(uint& p, uint& q, uint& r, uint& s, uint& t, uint& u,
                                 Irrep& irpq, Irrep& irrs) const
{
  ++u;
  do { // irpq
    do { // p
      do { // q
        if ( p_pgs->totIrrep(p,q) == irpq ) do { // irrs
          do { // r
            do { // s
              if ( p_pgs->totIrrep(r,s) == irrs && (r+1)*r/2 + s <= (p+1)*p/2 + q ) do { // t
                for ( ; u <= t; ++u ) {
                  if ( p_pgs->totIrrep(t,u) == p_pgs->product(irpq,irrs) &&
                       (t+1)*t/2 + u <= (r+1)*r/2 + s ) return true;
                } u = 0;
              ++t; } while ( t <= r ); t = 0;
            ++s; } while ( s <= r ); s = 0;
          ++r; } while ( r <= p ); r = 0;
        ++irrs; } while ( irrs < p_pgs->nIrreps() ); irrs = 0;
      ++q; } while ( q <= p ); q = 0;
    ++p; } while ( p < p_pgs->ntotorbs() ); p = 0;
  ++irpq; } while ( irpq < p_pgs->nIrreps() );
  return false;
}

inline bool Integ6::next_indices_nosym(uint& p, uint& q, uint& r, uint& s, uint& t, uint& u) const
{
  ++u;
  do { // p
    do { // q
      do { // r
        do { // s
          if ( (r+1)*r/2 + s <= (p+1)*p/2 + q ) do { // t
            for ( ; u <= t; ++u ) {
              if ( (t+1)*t/2 + u <= (r+1)*r/2 + s ) return true;
            } u = 0;
          ++t; } while ( t <= r ); t = 0;
        ++s; } while ( s <= r ); s = 0;
      ++r; } while ( r <= p ); r = 0;
    ++q; } while ( q <= p ); q = 0;
  ++p; } while ( p < p_pgs->ntotorbs());
  return false;
}

inline BlkIdx Integ6abb::index(uint p, uint q, uint r, uint s, uint t, uint u ) const
{
  if ( p < q ) {
    WARNRED6(p,q,r,s,t,u)
    std::swap(p,q);
  }
  if ( r < s ) {
    WARNRED6(p,q,r,s,t,u)
    std::swap(r,s);
  }
  if ( t < u ) {
    WARNRED6(p,q,r,s,t,u)
    std::swap(t,u);
  }
  if ( r*(r+1)/2+s < t*(t+1)/2+u ){
    WARNRED6(p,q,r,s,t,u)
    std::swap(r,t);
    std::swap(s,u);
  }

  Irrep pir = p_pgs->irrep(p),
        qir = p_pgs->irrep(q),
        rir = p_pgs->irrep(r),
        sir = p_pgs->irrep(s),
        tir = p_pgs->irrep(t),
        uir = p_pgs->irrep(u);
  uint nIrreps = p_pgs->nIrreps();
  BlkIdx blk_idx = _blocks[pir+nIrreps*(qir+nIrreps*(rir+nIrreps*(sir+nIrreps*(tir+nIrreps*uir))))];
  // indices relative to the block
  BlkIdx pb = p - p_pgs->_firstorb4irrep[pir],
         qb = q - p_pgs->_firstorb4irrep[qir],
         rb = r - p_pgs->_firstorb4irrep[rir],
         sb = s - p_pgs->_firstorb4irrep[sir],
         tb = t - p_pgs->_firstorb4irrep[tir],
         ub = u - p_pgs->_firstorb4irrep[uir];
  BlkIdx pqb, rsb, tub, lenrs, lentu;
  if ( pir == qir ) {
    // triangular index
    pqb = pb*(pb+1)/2+qb;
  } else {
    pqb = pb*p_pgs->_norb4irrep[qir]+qb;
  }
  if ( rir == sir ) {
    // triangular index
    rsb = rb*(rb+1)/2+sb;
    lenrs = p_pgs->_norb4irrep[rir]*(p_pgs->_norb4irrep[rir]+1)/2;
  } else {
    rsb = rb*p_pgs->_norb4irrep[sir]+sb;
    lenrs = p_pgs->_norb4irrep[rir]*p_pgs->_norb4irrep[sir];
  }
  if ( tir == uir ) {
    // triangular index
    tub = tb*(tb+1)/2+ub;
    lentu = p_pgs->_norb4irrep[tir]*(p_pgs->_norb4irrep[tir]+1)/2;
  } else {
    tub = tb*p_pgs->_norb4irrep[uir]+ub;
    lentu = p_pgs->_norb4irrep[tir]*p_pgs->_norb4irrep[uir];
  }
  if ( rir == tir ) {
      return pqb*(lenrs*(lenrs-1)/2+lentu) + rsb*(rsb+1)/2 + tub + blk_idx;
  } else {
    return (pqb*lenrs + rsb)*lentu + tub + blk_idx;
  }
}
inline bool Integ6abb::next_indices(uint& p, uint& q, uint& r, uint& s, uint& t, uint& u,
                                 Irrep& irpq, Irrep& irrs) const
{
  error("Not tested");
  ++u;
  do { // irpq
    do { // p
      do { // q
        if ( p_pgs->totIrrep(p,q) == irpq ) do { // irrs
          do { // r
            do { // s
              if ( p_pgs->totIrrep(r,s) == irrs ) do { // t
                for ( ; u <= t; ++u ) {
                  if ( p_pgs->totIrrep(t,u) == p_pgs->product(irpq,irrs) &&
                       (t+1)*t/2 + u <= (r+1)*r/2 + s ) return true;
                } u = 0;
              ++t; } while ( t <= r ); t = 0;
            ++s; } while ( s <= r ); s = 0;
          ++r; } while ( r <  p_pgs->ntotorbs() ); r = 0;
        ++irrs; } while ( irrs < p_pgs->nIrreps() ); irrs = 0;
      ++q; } while ( q <= p ); q = 0;
    ++p; } while ( p < p_pgs->ntotorbs() ); p = 0;
  ++irpq; } while ( irpq < p_pgs->nIrreps() );
  return false;
}

inline bool Integ6abb::next_indices_nosym(uint& p, uint& q, uint& r, uint& s, uint& t, uint& u) const
{
  error("Not tested");
  ++u;
  do { // p
    do { // q
      do { // r
        do { // s
          do { // t
            for ( ; u <= t; ++u ) {
              if ( (t+1)*t/2 + u <= (r+1)*r/2 + s ) return true;
            } u = 0;
          ++t; } while ( t <= r ); t = 0;
        ++s; } while ( s <= r ); s = 0;
      ++r; } while ( r < p_pgs->ntotorbs()); r = 0;
    ++q; } while ( q <= p ); q = 0;
  ++p; } while ( p < p_pgs->ntotorbs());
  return false;
}

inline BlkIdx Integ6st::index(uint p, uint q, uint r, uint s, uint t, uint u ) const
{
  if ( p < r ){
    WARNRED6(p,q,r,s,t,u)
    std::swap(p,r);
    std::swap(q,s);
  } else if ( p == r && q < s ){
    WARNRED6(p,q,r,s,t,u)
    std::swap(q,s);
  }
  if ( r < t || ( r == t && s < u ) ){
    WARNRED6(p,q,r,s,t,u)
    std::swap(r,t);
    std::swap(s,u);
    if ( p < r || ( p == r && q < s ) ){
      std::swap(p,r);
      std::swap(q,s);
    }
  }

  Irrep pir = p_pgs->irrep(p),
        qir = p_pgs->irrep(q),
        rir = p_pgs->irrep(r),
        sir = p_pgs->irrep(s),
        tir = p_pgs->irrep(t),
        uir = p_pgs->irrep(u);
  uint nIrreps = p_pgs->nIrreps();
  BlkIdx blk_idx = _blocks[pir+nIrreps*(qir+nIrreps*(rir+nIrreps*(sir+nIrreps*(tir+nIrreps*uir))))];
  // indices relative to the block
  BlkIdx pb = p - p_pgs->_firstorb4irrep[pir],
         qb = q - p_pgs->_firstorb4irrep[qir],
         rb = r - p_pgs->_firstorb4irrep[rir],
         sb = s - p_pgs->_firstorb4irrep[sir],
         tb = t - p_pgs->_firstorb4irrep[tir],
         ub = u - p_pgs->_firstorb4irrep[uir];
  BlkIdx pqb, rsb, tub, lenrs, lentu;
  pqb = pb*p_pgs->_norb4irrep[qir]+qb;
  rsb = rb*p_pgs->_norb4irrep[sir]+sb;
  lenrs = p_pgs->_norb4irrep[rir]*p_pgs->_norb4irrep[sir];
  tub = tb*p_pgs->_norb4irrep[uir]+ub;
  lentu = p_pgs->_norb4irrep[tir]*p_pgs->_norb4irrep[uir];
  if ( pir == rir ) {
    // triangular index
    if ( pir == tir ) {
      return pqb*(pqb+1)*(pqb+2)/6 + rsb*(rsb+1)/2 + tub + blk_idx;
    } else {
      return (pqb*(pqb+1)/2 + rsb)*lentu + tub + blk_idx;
    }
  } else if ( rir == tir ) {
    return pqb*(lenrs*(lenrs-1)/2+lentu) + rsb*(rsb+1)/2 + tub + blk_idx;
  } else {
    return (pqb*lenrs + rsb)*lentu + tub + blk_idx;
  }
}
inline bool Integ6st::next_indices(uint& p, uint& q, uint& r, uint& s, uint& t, uint& u,
                                 Irrep& irpq, Irrep& irrs) const
{
  ++u;
  do { // irpq
    do { // p
      do { // q
        if ( p_pgs->totIrrep(p,q) == irpq ) do { // irrs
          do { // r
            do { // s
              if ( p_pgs->totIrrep(r,s) == irrs ) do { // t
                if ( t == r ) {
                  for ( ; u <= s; ++u ) {
                    if ( p_pgs->totIrrep(t,u) == p_pgs->product(irpq,irrs) ) return true;
                  }
                } else {
                  for ( ; u < p_pgs->ntotorbs(); ++u ) {
                    if ( p_pgs->totIrrep(t,u) == p_pgs->product(irpq,irrs) ) return true;
                  }
                } u = 0;
              ++t; } while ( t <= r ); t = 0;
            ++s; } while (  s <= q || ( r != p && s < p_pgs->ntotorbs() ) ); s = 0;
          ++r; } while ( r <= p ); r = 0;
        ++irrs; } while ( irrs < p_pgs->nIrreps() ); irrs = 0;
      ++q; } while ( q < p_pgs->ntotorbs() ); q = 0;
    ++p; } while ( p < p_pgs->ntotorbs() ); p = 0;
  ++irpq; } while ( irpq < p_pgs->nIrreps() );
  return false;
}

inline bool Integ6st::next_indices_nosym(uint& p, uint& q, uint& r, uint& s, uint& t, uint& u) const
{
  ++u;
  do { // p
    do { // q
      do { // r
        do { // s
          do { // t
            if ( t == r ) {
              if ( u <= s ) return true;
            } else {
              if ( u <= p_pgs->ntotorbs() ) return true;
            }
            u = 0;
          ++t; } while ( t <= r ); t = 0;
        ++s; } while ( s <= q || ( r != p && s < p_pgs->ntotorbs() ) ); s = 0;
      ++r; } while ( r <= p ); r = 0;
    ++q; } while ( q < p_pgs->ntotorbs()); q = 0;
  ++p; } while ( p < p_pgs->ntotorbs());
  return false;
}

inline BlkIdx Integ6stabb::index(uint p, uint q, uint r, uint s, uint t, uint u ) const
{
  if ( r < t || ( r == t && s < u ) ){
    WARNRED6(p,q,r,s,t,u)
    std::swap(r,t);
    std::swap(s,u);
  }

  Irrep pir = p_pgs->irrep(p),
        qir = p_pgs->irrep(q),
        rir = p_pgs->irrep(r),
        sir = p_pgs->irrep(s),
        tir = p_pgs->irrep(t),
        uir = p_pgs->irrep(u);
  uint nIrreps = p_pgs->nIrreps();
  BlkIdx blk_idx = _blocks[pir+nIrreps*(qir+nIrreps*(rir+nIrreps*(sir+nIrreps*(tir+nIrreps*uir))))];
  // indices relative to the block
  BlkIdx pb = p - p_pgs->_firstorb4irrep[pir],
         qb = q - p_pgs->_firstorb4irrep[qir],
         rb = r - p_pgs->_firstorb4irrep[rir],
         sb = s - p_pgs->_firstorb4irrep[sir],
         tb = t - p_pgs->_firstorb4irrep[tir],
         ub = u - p_pgs->_firstorb4irrep[uir];
  BlkIdx pqb, rsb, tub, lenrs, lentu;
  pqb = pb*p_pgs->_norb4irrep[qir]+qb;
  rsb = rb*p_pgs->_norb4irrep[sir]+sb;
  lenrs = p_pgs->_norb4irrep[rir]*p_pgs->_norb4irrep[sir];
  tub = tb*p_pgs->_norb4irrep[uir]+ub;
  lentu = p_pgs->_norb4irrep[tir]*p_pgs->_norb4irrep[uir];
  if ( rir == tir ) {
    return pqb*(lenrs*(lenrs-1)/2+lentu) + rsb*(rsb+1)/2 + tub + blk_idx;
  } else {
    return (pqb*lenrs + rsb)*lentu + tub + blk_idx;
  }
}
inline bool Integ6stabb::next_indices(uint& p, uint& q, uint& r, uint& s, uint& t, uint& u,
                                 Irrep& irpq, Irrep& irrs) const
{
  error("Not tested");
  ++u;
  do { // irpq
    do { // p
      do { // q
        if ( p_pgs->totIrrep(p,q) == irpq ) do { // irrs
          do { // r
            do { // s
              if ( p_pgs->totIrrep(r,s) == irrs ) do { // t
                if ( t == r ) {
                  for ( ; u <= s; ++u ) {
                    if ( p_pgs->totIrrep(t,u) == p_pgs->product(irpq,irrs) ) return true;
                  }
                } else {
                  for ( ; u < p_pgs->ntotorbs(); ++u ) {
                    if ( p_pgs->totIrrep(t,u) == p_pgs->product(irpq,irrs) ) return true;
                  }
                } u = 0;
              ++t; } while ( t <= r ); t = 0;
            ++s; } while ( s < p_pgs->ntotorbs() ); s = 0;
          ++r; } while ( r <  p_pgs->ntotorbs() ); r = 0;
        ++irrs; } while ( irrs < p_pgs->nIrreps() ); irrs = 0;
      ++q; } while ( q < p_pgs->ntotorbs() ); q = 0;
    ++p; } while ( p < p_pgs->ntotorbs() ); p = 0;
  ++irpq; } while ( irpq < p_pgs->nIrreps() );
  return false;
}

inline bool Integ6stabb::next_indices_nosym(uint& p, uint& q, uint& r, uint& s, uint& t, uint& u) const
{
  error("Not tested");
  ++u;
  do { // p
    do { // q
      do { // r
        do { // s
          do { // t
            if ( t == r ) {
              if ( u <= s ) return true;
            } else {
              if ( u <= p_pgs->ntotorbs() ) return true;
            }
            u = 0;
          ++t; } while ( t <= r ); t = 0;
        ++s; } while ( s < p_pgs->ntotorbs()); s = 0;
      ++r; } while ( r < p_pgs->ntotorbs()); r = 0;
    ++q; } while ( q < p_pgs->ntotorbs()); q = 0;
  ++p; } while ( p < p_pgs->ntotorbs());
  return false;
}

#undef WARNRED2
#undef WARNRED4
#undef WARNRED6
} // namespace HamDump
#endif
