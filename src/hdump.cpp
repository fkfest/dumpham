#include "hdump.h"

namespace HamDump {

Hdump::Hdump(std::string fcidump, bool verbose) : _dump(fcidump)
{
  if (verbose) xout << "\n Process file "<< fcidump <<std::endl;
  FDPar NELEC = _dump.parameter("NELEC");
  if (verbose) xout << "NELEC=" << NELEC[0] << std::endl;
  FDPar MS2 = _dump.parameter("MS2");
  if (verbose) xout << "MS2=" << MS2[0] << std::endl;
  FDPar NORB = _dump.parameter("NORB");
  if (verbose) xout << "NORB=" << NORB[0] << std::endl;
  FDPar ISYM = _dump.parameter("ISYM");
  if (verbose) xout << "ISYM=" << ISYM[0] << std::endl;
  FDPar IUHF = _dump.parameter("IUHF");
  if (verbose) xout << "IUHF=" << IUHF[0] << std::endl;
  FDPar ORBSYM = _dump.parameter("ORBSYM");
  if (verbose) {
    xout << "ORBSYM="; 
    for ( const auto& s: ORBSYM)
      xout << s << ","; 
    xout<<std::endl;
  }
  FDPar OCC = _dump.parameter("OCC");
  check_input_norbs(OCC,"occ",verbose);
  FDPar CLOSED = _dump.parameter("CLOSED");
  check_input_norbs(CLOSED,"closed",verbose);
  FDPar CORE = _dump.parameter("CORE");
  check_input_norbs(CORE,"core",verbose);
  
  FDPar ST = _dump.parameter("ST");
  
  int nn = NORB[0];
  if ( nn < 0 ) {
    error("NORB < 0 in FCIDUMP!");
  }
  _norb = nn;
  _nelec = NELEC[0];
  _ms2 = MS2[0];
  _sym = ISYM[0];
  if ( _sym > 0 ) --_sym; // make zero-based
  _pgs = PGSym(ORBSYM);
  _escal = 0.0;
 
  assert( CORE.size() > 0 && OCC.size() > 0 && CLOSED.size() > 0 );
  if ( CORE.size() > 1 || CORE[0] > 0 ) {
    _core = CORE;
    FDPar norb4irs(_pgs.norbs_in_irreps());
    _core.resize(norb4irs.size(),0);
    for ( Irrep ir = 0; ir < norb4irs.size(); ++ir )
      norb4irs[ir] += _core[ir];
    _pgs_wcore = PGSym(norb4irs,false);
  }
  if ( OCC.size() > 1 || OCC[0] > 0 ) {
    _occ = OCC;
    _occ.resize(_pgs.nIrreps(),0);
  }
  if ( CLOSED.size() > 1 || CLOSED[0] > 0 ) {
    _clos = CLOSED;
    _clos.resize(_pgs.nIrreps(),0);
  }
  _uhf = bool(IUHF[0]);
  _simtra = bool(ST[0]);
}

Hdump::Hdump(const Hdump& hd1, const Hdump& hd2)
{
#define CHECKSET(xx,mes) if ( hd1.xx != hd2.xx ) error("Mismatch in "+ std::string(mes)); xx=hd1.xx;
  CHECKSET(_norb,"number of orbitals");
  CHECKSET(_nelec,"number of electrons");
  CHECKSET(_ms2,"spin state");
  CHECKSET(_sym,"total symmetry");
  CHECKSET(_pgs,"orbital space");
  _escal = 0.0;
  CHECKSET(_core,"core orbitals");
  CHECKSET(_pgs_wcore,"orbital space with core orbitals");
  CHECKSET(_occ,"occupied orbitals");
  CHECKSET(_clos,"closed-shell orbitals");
#undef CHECKSET
  _uhf = hd1._uhf || hd2._uhf;
  _simtra = hd1._simtra || hd2._simtra;
}

Hdump::Hdump(const Hdump& hd, int i_uhf, int i_simtra)
{
#define SETVAL(xx) xx=hd.xx;
  SETVAL(_norb);
  SETVAL(_nelec);
  SETVAL(_ms2);
  SETVAL(_sym);
  SETVAL(_pgs);
  _escal = 0.0;
  SETVAL(_core);
  SETVAL(_pgs_wcore);
  SETVAL(_occ);
  SETVAL(_clos);
  SETVAL(_uhf);
  SETVAL(_simtra);
#undef SETVAL
  if ( i_uhf < 0 ) _uhf = false;
  if ( i_uhf > 0 ) _uhf = true;
  if ( i_simtra < 0 ) _simtra = false;
  if ( i_simtra > 0 ) _simtra = true;
}

void Hdump::alloc_ints()
{
  assert(!allocated());
  if (_simtra) {
    // similarity transformed integrals
    _oneel.emplace_back(new Integ2st(_pgs));
    _twoel.emplace_back(new Integ4st(_pgs));
    if ( _uhf ) {
      _oneel.emplace_back(new Integ2st(_pgs)); // h_bb
      _twoel.emplace_back(new Integ4st(_pgs)); // (bb|bb)
      _twoel.emplace_back(new Integ4stab(_pgs)); // (aa|bb)
    }
  } else {
    // normal case
    _oneel.emplace_back(new Integ2(_pgs));
    _twoel.emplace_back(new Integ4(_pgs));
    if ( _uhf ) {
      _oneel.emplace_back(new Integ2(_pgs)); // h_bb
      _twoel.emplace_back(new Integ4(_pgs)); // (bb|bb)
      _twoel.emplace_back(new Integ4ab(_pgs)); // (aa|bb)
    }
  }
  
#ifdef _DEBUG
  if (_simtra) {
    assert( dynamic_cast<Integ4st*>(_twoel[aaaa].get()) );
    assert( dynamic_cast<Integ2st*>(_oneel[aa].get()) );
    if ( _uhf ) {
      assert( dynamic_cast<Integ2st*>(_oneel[bb].get()) );
      assert( dynamic_cast<Integ4st*>(_twoel[bbbb].get()) );
      assert( dynamic_cast<Integ4stab*>(_twoel[aabb].get()) );
    }
  } else {
    assert( dynamic_cast<Integ4*>(_twoel[aaaa].get()) );
    assert( dynamic_cast<Integ2*>(_oneel[aa].get()) );
    if ( _uhf ) {
      assert( dynamic_cast<Integ2*>(_oneel[bb].get()) );
      assert( dynamic_cast<Integ4*>(_twoel[bbbb].get()) );
      assert( dynamic_cast<Integ4ab*>(_twoel[aabb].get()) );
    }
  }
//   check_addressing_integrals();
#endif
}

void Hdump::read_dump()
{
  alloc_ints();
  int i,j,k,l;
  double value;
  
  FCIdump::integralType type;
  _dump.rewind();
  // read integrals
  type = _dump.nextIntegral(i,j,k,l,value);
  do {
    switch(type){
      case FCIdump::I2aa:
        xout << "I2aa" << std::endl;
        if (_simtra)
          readrec( static_cast<Integ4st*>(_twoel[aaaa].get()),i,j,k,l,value,type);
        else
          readrec( static_cast<Integ4*>(_twoel[aaaa].get()),i,j,k,l,value,type);
        break;
      case FCIdump::I2bb:
        xout << "I2bb" << std::endl;
        if (_simtra)
          readrec( static_cast<Integ4st*>(_twoel[bbbb].get()),i,j,k,l,value,type);
        else
          readrec( static_cast<Integ4*>(_twoel[bbbb].get()),i,j,k,l,value,type);
        break;
      case FCIdump::I2ab:
        xout << "I2ab" << std::endl;
        if (_simtra)
          readrec( static_cast<Integ4stab*>(_twoel[aabb].get()),i,j,k,l,value,type);
        else
          readrec( static_cast<Integ4ab*>(_twoel[aabb].get()),i,j,k,l,value,type);
        break;
      case FCIdump::I1a:
        xout << "I1a" << std::endl;
        if (_simtra)
          readrec( static_cast<Integ2st*>(_oneel[aa].get()),i,j,value,type);
        else
          readrec( static_cast<Integ2*>(_oneel[aa].get()),i,j,value,type);
        break;
      case FCIdump::I1b:
        xout << "I1b" << std::endl;
        if (_simtra)
          readrec( static_cast<Integ2st*>(_oneel[bb].get()),i,j,value,type);
        else
          readrec( static_cast<Integ2*>(_oneel[bb].get()),i,j,value,type);
        break;
      case FCIdump::I0:
        xout << "I0" << std::endl;
        _escal = value;
        type = _dump.nextIntegral(i,j,k,l,value);
        break;
      case FCIdump::endOfRecord:
        type = _dump.nextIntegral(i,j,k,l,value);
        break;
      default:
        if ( type != FCIdump::endOfFile ) error("Unknown integral type","hdump.cpp");
    }
      
  } while (type != FCIdump::endOfFile);
}
template<typename T>
void Hdump::readrec(T* pInt, int& i, int& j, int& k, int& l, double& value, FCIdump::integralType& curtype)
{
  #ifdef _DEBUG
    #ifndef MOLPRO
    bool redunwarn = Input::iPars["ham"]["redunwarning"];
    #else
    bool redunwarn = true;
    #endif
    #define REDUNWAR ,redunwarn
  #else
    #define REDUNWAR
  #endif
  FCIdump::integralType type;
  do {
    pInt->set(i-1,j-1,k-1,l-1, value REDUNWAR);
  } while ( (type = _dump.nextIntegral(i,j,k,l,value)) == curtype ); 
  curtype = type;
  #undef REDUNWAR
}
template<typename T>
void Hdump::readrec(T* pInt, int& i, int& j, double& value, FCIdump::integralType& curtype)
{
  #ifdef _DEBUG
    #ifndef MOLPRO
    bool redunwarn = Input::iPars["ham"]["redunwarning"];
    #else
    bool redunwarn = true;
    #endif
    #define REDUNWAR ,redunwarn
  #else
    #define REDUNWAR
  #endif
  FCIdump::integralType type;
  int k,l;
  do {
    pInt->set(i-1,j-1, value REDUNWAR);
  } while ( (type = _dump.nextIntegral(i,j,k,l,value)) == curtype ); 
  curtype = type;
  #undef REDUNWAR
}
void Hdump::skiprec(int& i, int& j, int& k, int& l, double& value, FCIdump::integralType& curtype)
{
  FCIdump::integralType type;
  do {
  } while ( (type = _dump.nextIntegral(i,j,k,l,value)) == curtype ); 
  curtype = type;
}

void Hdump::check_input_norbs(FDPar& orb, const std::string& kind, bool verbose) const
{
#ifndef MOLPRO
  const TParArray& inporb = Input::aPars["orbs"][kind];
  if ( inporb.size() > 0 ) {
    orb.clear();
    apars2nums<int>(orb,inporb,std::dec);
  }
#endif
  if ( verbose && ( orb.size() > 1 || orb[0] > 0 ) ) {
    xout << kind << "=";
    for ( const auto& s: orb)
      xout << s << ","; 
    xout<<std::endl;
  }
}
void Hdump::store(std::string fcidump)
{
  bool newfile = _dump.fileName().empty();
  bool nosym = false;
#ifndef MOLPRO
  nosym = Input::iPars["ham"]["nosym"];
  FDPar ORBSYM, ORBSYM_SAV;
#endif
  if ( newfile ) {
    // create a new FCIDUMP 
    if ( _simtra ) _dump.addParameter("ST",std::vector<int>(1,1));
    if ( _uhf ) _dump.addParameter("IUHF",std::vector<int>(1,1));
    _dump.addParameter("ISYM",std::vector<int>(1,_sym+1));
    _dump.addParameter("ORBSYM",_pgs.orbsym());
    _dump.addParameter("MS2",std::vector<int>(1,_ms2));
    _dump.addParameter("NELEC",std::vector<int>(1,_nelec));
    _dump.addParameter("NORB",std::vector<int>(1,_norb));
  } else {
#ifdef MOLPRO
    error("Cannot modify FCIDUMP "+_dump.fileName());
#else
    if ( nosym ) {
      //remove the point-group symmetry
      ORBSYM_SAV = _dump.parameter("ORBSYM");
      ORBSYM.resize(ORBSYM_SAV.size());
      for ( auto& s: ORBSYM)
        s = 1;
      _dump.modifyParameter("ORBSYM",ORBSYM);
    }
#endif
  }
  if (_dump.write(fcidump,FCIdump::FileFormatted,false))
    xout << "will be written to file " << fcidump << std::endl;
  else
    error("failure to write to file "+fcidump);
  
  if (_simtra) {
    // similarity transformed
    Integ4st * pI4aa = 0;
    Integ4stab * pI4ab = 0;
    Integ2st * pI2 = 0;
    if (nosym) {
      store_without_symmetry(pI2,pI4aa,pI4ab);
    } else {
      store_with_symmetry(pI2,pI4aa,pI4ab);
    }
  } else {
    // normal case
    Integ4 * pI4aa = 0;
    Integ4ab * pI4ab = 0;
    Integ2 * pI2 = 0;
    if (nosym) {
      store_without_symmetry(pI2,pI4aa,pI4ab);
    } else {
      store_with_symmetry(pI2,pI4aa,pI4ab);
    }
  }
#ifndef MOLPRO
  _dump.close_outputfile();
  if ( nosym && !newfile ) {
    //restore the point-group symmetry
    _dump.modifyParameter("ORBSYM",ORBSYM_SAV);
  }
#endif
}

template<typename I2, typename I4aa, typename I4ab>
void Hdump::store_with_symmetry(I2* pI2, I4aa* pI4aa, I4ab* pI4ab) const
{
  pI4aa = dynamic_cast<I4aa*>(_twoel[aaaa].get()); assert(pI4aa);
  storerec4_sym(pI4aa);
  if (_uhf) {
    _dump.writeIntegral(0,0,0,0,0.0);
    pI4aa = dynamic_cast<I4aa*>(_twoel[bbbb].get()); assert(pI4aa);
    storerec4_sym(pI4aa);
    _dump.writeIntegral(0,0,0,0,0.0);
    pI4ab = dynamic_cast<I4ab*>(_twoel[aabb].get()); assert(pI4ab);
    storerec4_sym(pI4ab);
    _dump.writeIntegral(0,0,0,0,0.0);
  }
  pI2 = dynamic_cast<I2*>(_oneel[aa].get()); assert(pI2);
  storerec2_sym(pI2);
  if (_uhf) {
    _dump.writeIntegral(0,0,0,0,0.0);
    pI2 = dynamic_cast<I2*>(_oneel[bb].get()); assert(pI2);
    storerec2_sym(pI2);
    _dump.writeIntegral(0,0,0,0,0.0);
  }
  _dump.writeIntegral(0,0,0,0,_escal); 
  
}

template<typename I2, typename I4aa, typename I4ab>
void Hdump::store_without_symmetry(I2* pI2, I4aa* pI4aa, I4ab* pI4ab) const
{
  pI4aa = dynamic_cast<I4aa*>(_twoel[aaaa].get()); assert(pI4aa);
  storerec4_nosym(pI4aa);
  if (_uhf) {
    _dump.writeIntegral(0,0,0,0,0.0);
    pI4aa = dynamic_cast<I4aa*>(_twoel[bbbb].get()); assert(pI4aa);
    storerec4_nosym(pI4aa);
    _dump.writeIntegral(0,0,0,0,0.0);
    pI4ab = dynamic_cast<I4ab*>(_twoel[aabb].get()); assert(pI4ab);
    storerec4_nosym(pI4ab);
    _dump.writeIntegral(0,0,0,0,0.0);
  }
  pI2 = dynamic_cast<I2*>(_oneel[aa].get()); assert(pI2);
  storerec2_nosym(pI2);
  if (_uhf) {
    _dump.writeIntegral(0,0,0,0,0.0);
    pI2 = dynamic_cast<I2*>(_oneel[bb].get()); assert(pI2);
    storerec2_nosym(pI2);
    _dump.writeIntegral(0,0,0,0,0.0);
  }
  _dump.writeIntegral(0,0,0,0,_escal); 
  
}
template<typename T>
void Hdump::storerec4_sym(const T * pInt) const
{
  uint i = 0, j = 0, k = 0, l = 0;
  Irrep isym = 0;
  do {
    _dump.writeIntegral(i+1,j+1,k+1,l+1,pInt->get(i,j,k,l));
  } while (pInt->next_indices(i,j,k,l,isym));
}
template<typename T>
void Hdump::storerec2_sym(const T * pInt) const
{
  uint i = 0, j = 0;
  do {
    _dump.writeIntegral(i+1,j+1,0,0,pInt->get(i,j));
  } while (pInt->next_indices(i,j));
}
template<typename T>
void Hdump::storerec4_nosym(const T * pInt) const
{
  uint i = 0, j = 0, k = 0, l = 0;
  do {
    _dump.writeIntegral(i+1,j+1,k+1,l+1,pInt->get_with_pgs(i,j,k,l));
  } while (pInt->next_indices_nosym(i,j,k,l));
}
template<typename T>
void Hdump::storerec2_nosym(const T * pInt) const
{
  uint i = 0, j = 0;
  do {
    _dump.writeIntegral(i+1,j+1,0,0,pInt->get_with_pgs(i,j));
  } while (pInt->next_indices_nosym(i,j));
}

template<typename T, typename U>
void Hdump::copy_int4(T* pDest, const U* pSrc, bool add, bool sym)
{
  uint i = 0, j = 0, k = 0, l = 0;
  Irrep isym = 0;
  if ( add ) {
    if ( sym ) { 
      // add, symmetrization
      do {
        BlkIdx indxD = pDest->index(i,j,k,l);
        pDest->set(indxD, 0.5*(pSrc->get(i,j,k,l)+pSrc->get(j,i,l,k))+pDest->get(indxD));
      } while (pDest->next_indices(i,j,k,l,isym));
    } else { 
      // add, no symmetrization
      do {
        BlkIdx indxD = pDest->index(i,j,k,l);
        pDest->set(indxD, pSrc->get(i,j,k,l)+pDest->get(indxD));
      } while (pDest->next_indices(i,j,k,l,isym));
    }
  } else {
    if ( sym ) { 
      // don't add, symmetrization
      do {
        BlkIdx indxD = pDest->index(i,j,k,l);
        pDest->set(indxD, 0.5*(pSrc->get(i,j,k,l)+pSrc->get(j,i,l,k)));
      } while (pDest->next_indices(i,j,k,l,isym));
    } else {
      // don't add, no symmetrization
      do {
        BlkIdx indxD = pDest->index(i,j,k,l);
        pDest->set(indxD, pSrc->get(i,j,k,l));
      } while (pDest->next_indices(i,j,k,l,isym));
    }
  }
}
template<typename T, typename U>
void Hdump::copy_int2(T* pDest, const U* pSrc, bool add, bool sym)
{
  uint i = 0, j = 0;
  if ( add ) {
    if ( sym ) { 
      // add, symmetrization
      do {
        BlkIdx indxD = pDest->index(i,j);
        pDest->set(indxD, 0.5*(pSrc->get(i,j)+pSrc->get(j,i))+pDest->get(indxD));
      } while (pDest->next_indices(i,j));
    } else { 
      // add, no symmetrization
      do {
        BlkIdx indxD = pDest->index(i,j);
        pDest->set(indxD, pSrc->get(i,j)+pDest->get(indxD));
      } while (pDest->next_indices(i,j));
    }
  } else {
    if ( sym ) { 
      // don't add, symmetrization
      do {
        BlkIdx indxD = pDest->index(i,j);
        pDest->set(indxD, 0.5*(pSrc->get(i,j)+pSrc->get(j,i)));
      } while (pDest->next_indices(i,j));
    } else {
      // don't add, no symmetrization
      do {
        BlkIdx indxD = pDest->index(i,j);
        pDest->set(indxD, pSrc->get(i,j));
      } while (pDest->next_indices(i,j));
    }
  }
}
template<typename DI2, typename DI4aa, typename DI4ab, typename SI2, typename SI4aa, typename SI4ab>
void Hdump::copy_ints(DI2* pDI2, DI4aa* pDI4aa, DI4ab* pDI4ab, 
                      SI2* pSI2, SI4aa* pSI4aa, SI4ab* pSI4ab, const Hdump& hd, bool add, bool sym)
{
  pDI4aa = dynamic_cast<DI4aa*>(_twoel[aaaa].get()); assert(pDI4aa);
  pSI4aa = dynamic_cast<SI4aa*>(hd._twoel[aaaa].get()); assert(pSI4aa);
  copy_int4(pDI4aa,pSI4aa,add,sym);
  pDI2 = dynamic_cast<DI2*>(_oneel[aa].get()); assert(pDI2);
  pSI2 = dynamic_cast<SI2*>(hd._oneel[aa].get()); assert(pSI2);
  copy_int2(pDI2,pSI2,add,sym);
  if ( _uhf ) {
    pDI4aa = dynamic_cast<DI4aa*>(_twoel[bbbb].get()); assert(pDI4aa);
    if ( hd._uhf ) {
      pSI4aa = dynamic_cast<SI4aa*>(hd._twoel[bbbb].get()); assert(pSI4aa);
    }
    copy_int4(pDI4aa,pSI4aa,add,sym);
    pDI4ab = dynamic_cast<DI4ab*>(_twoel[aabb].get()); assert(pDI4ab);
    if ( hd._uhf ) {
      pSI4ab = dynamic_cast<SI4ab*>(hd._twoel[aabb].get()); assert(pSI4ab);
      copy_int4(pDI4ab,pSI4ab,add,sym);
    } else {
      copy_int4(pDI4ab,pSI4aa,add,sym);
    }
    pDI2 = dynamic_cast<DI2*>(_oneel[bb].get()); assert(pDI2);
    if ( hd._uhf ) {
      pSI2 = dynamic_cast<SI2*>(hd._oneel[bb].get()); assert(pSI2);
    }
    copy_int2(pDI2,pSI2,add,sym);
  }
  if (add) {
    _escal += hd._escal;
  } else {
    _escal = hd._escal;
  }
}

void Hdump::import(const Hdump& hd, bool add)
{
  if (_pgs != hd._pgs) error("Different orbital spaces in two hamiltonians");
  if (!add) {
    alloc_ints();
  } else {
    assert(allocated());
  }
  if ( _simtra ) {
    // expand the permutation symmetry
    Integ4st * pDI4aa = 0;
    Integ4stab * pDI4ab = 0;
    Integ2st * pDI2 = 0;
    if ( hd._simtra ) {
      Integ4st * pSI4aa = 0;
      Integ4stab * pSI4ab = 0;
      Integ2st * pSI2 = 0;
      copy_ints(pDI2,pDI4aa,pDI4ab,pSI2,pSI4aa,pSI4ab,hd,add);
    } else {
      Integ4 * pSI4aa = 0;
      Integ4ab * pSI4ab = 0;
      Integ2 * pSI2 = 0;
      copy_ints(pDI2,pDI4aa,pDI4ab,pSI2,pSI4aa,pSI4ab,hd,add);
    }
  } else {
    // create the permutation symmetry
    Integ4 * pDI4aa = 0;
    Integ4ab * pDI4ab = 0;
    Integ2 * pDI2 = 0;
    if ( !hd._simtra ) {
      Integ4 * pSI4aa = 0;
      Integ4ab * pSI4ab = 0;
      Integ2 * pSI2 = 0;
      copy_ints(pDI2,pDI4aa,pDI4ab,pSI2,pSI4aa,pSI4ab,hd,add);
    } else {
      Integ4st * pSI4aa = 0;
      Integ4stab * pSI4ab = 0;
      Integ2st * pSI2 = 0;
      copy_ints(pDI2,pDI4aa,pDI4ab,pSI2,pSI4aa,pSI4ab,hd,add,true);
    }
  }
}

uint Hdump::nclostot() const
{
  uint nclos = (_nelec - _ms2)/2;
  if ( nclos*2 + _ms2 != _nelec ) {
    xout << "NELEC: " << _nelec << " MS2: " << _ms2 << std::endl;
    error("Mismatch in NELEC and MS2 in Hdump");
  }
  return nclos;
}
FDPar Hdump::nclos_wcore() const
{
  FDPar clco(_clos); 
  for ( uint io = 0; io < std::max(clco.size(),_core.size()); ++io )
    clco[io] += _core[io];
  for ( uint io = clco.size(); io < _core.size(); ++io )
    clco.push_back(_core[io]);
  return clco; 
}
FDPar Hdump::nocc_wcore() const
{
  FDPar occo(_occ); 
  for ( uint io = 0; io < std::max(occo.size(),_core.size()); ++io )
    occo[io] += _core[io];
  for ( uint io = occo.size(); io < _core.size(); ++io )
    occo.push_back(_core[io]);
  return occo; 
}

void Hdump::scale(double scal)
{
  (void) scal;
  error("Not implemented!");
}


void Hdump::check_addressing_integrals() const
{
  if (_simtra) return;
  xout << "test one el" << std::endl;
  int oldidx=-1;
  for (uint i = 1; i <= _norb; ++i )
    for (uint j = 1; j <= i; ++j ) {
      if ( _pgs.totIrrep(i-1,j-1) != 0 ) continue;
      BlkIdx idx = (_oneel[aa].get())->index(i-1,j-1);
      xout << i << " " << j << " " << idx << std::endl;
      if ( idx != BlkIdx(oldidx+1)) error("Indices are not consecutive in 1-el operator","check_addressing_integrals");
      oldidx = idx;
    }
  xout << "test two el" << std::endl;
  uint idxcan = 0;
  for ( Irrep isym = 0; isym < _pgs.nIrreps(); ++isym ) {
    for ( uint i = 1; i <= _norb; ++i) {
      for ( uint j = 1; j <= i; ++j ) {
        if ( _pgs.totIrrep(i-1,j-1) != isym ) continue;
        for ( uint k = 1; k <= i; ++k ) {
          for ( uint l = 1; l <= k; ++l ) {
            if ( _pgs.totIrrep(k-1,l-1) != isym ) continue;
            uint ij = (i-1)*i/2 + j;
            uint kl = (k-1)*k/2 + l;
            if ( kl <= ij ){
              BlkIdx idx = (_twoel[aaaa].get())->index(i-1,j-1,k-1,l-1);
              xout << i << " " << j << " " << k << " " << l << " " << idx << "  " << idxcan << std::endl;
              ++idxcan;
            }
          }
        }
      }
    }
  }
   xout << "n1el: " << _oneel[aa]->nelem() << std::endl;
   xout << "n2el: " << _twoel[aaaa]->nelem() << std::endl;
}

} //namespace HamDump
