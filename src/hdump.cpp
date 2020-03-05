#include "hdump.h"
#include <iomanip> 
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
  _osord = OrbOrder(ORBSYM);
  _osord.reorder(ORBSYM);
  FDPar OCC = _dump.parameter("OCC");
  check_input_norbs(OCC,"occ",verbose);
  FDPar CLOSED = _dump.parameter("CLOSED");
  check_input_norbs(CLOSED,"closed",verbose);
  FDPar CORE = _dump.parameter("CORE");
  check_input_norbs(CORE,"core",verbose);
  
  FDPar ST = _dump.parameter("ST");
  FDPar DM = _dump.parameter("DM");
  
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
    CORE.resize(8,0);
    set_ncore(&CORE[0]);
  }
  if ( OCC.size() > 1 || OCC[0] > 0 ) {
    _occ = OCC;
    uint nIrrepsOcc = _pgs.guess_nIrreps(_occ);
    _occ.resize(nIrrepsOcc,0);
    add_or_subtract_core(_occ,"OCC",false);
    _occ.resize(_pgs.nIrreps(),0);
  }
  if ( CLOSED.size() > 1 || CLOSED[0] > 0 ) {
    _clos = CLOSED;
    uint nIrrepsClos = _pgs.guess_nIrreps(_clos);
    _clos.resize(nIrrepsClos,0);
    add_or_subtract_core(_clos,"CLOSED",false);
    _clos.resize(_pgs.nIrreps(),0);
    if ( _occ.empty() ) {
      // assume closed-shell case
      _occ = _clos;
    }
  }
  _uhf = bool(IUHF[0]);
  _simtra = bool(ST[0]);
  _dm = bool(DM[0]);
  gen_spinorbsref();
  sanity_check();
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
  CHECKSET(_dm,"density matrix flag");
#undef CHECKSET
  _uhf = hd1._uhf || hd2._uhf;
  _simtra = hd1._simtra || hd2._simtra;
  gen_spinorbsref();
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
  SETVAL(_dm);
#undef SETVAL
  if ( i_uhf < 0 ) _uhf = false;
  if ( i_uhf > 0 ) _uhf = true;
  if ( i_simtra < 0 ) _simtra = false;
  if ( i_simtra > 0 ) _simtra = true;
  gen_spinorbsref();
}

void Hdump::sanity_check() const
{
  if ( _occ.size() < _clos.size() )
    error("Number of values in OCC is smaller than in CLOSED!");
  if ( _core.size() > 0 || _occ.size() > 0 || _clos.size() > 0 ) {
    uint norb_core = 0, norb_closed = 0, norb_occ = 0;
    for ( auto io: _core ) {
      if ( io < 0 ) error("Negative number of orbitals in CORE!");
      norb_core += io;
    }
    auto iocc = _occ.begin();
    for ( auto io: _clos ) {
      if ( io < 0 ) error("Negative number of orbitals in CLOSED!");
      if ( *iocc < io ) error("Number of OCC orbitals in a symmetry larger than number of CLOSED!");
      ++iocc;
      norb_closed += io;
    }
    for ( auto io: _occ ) {
      if ( io < 0 ) error("Negative number of orbitals in OCC!");
      norb_occ += io;
    }
    xout << "norb_closed: " << norb_closed << " norb_occ: " << norb_occ << std::endl;
    if ( norb_occ > _norb ) {
      error("Number of OCC orbitals exceeds the total number of orbitals!");
    }
    // number of electrons without core!
    if ( norb_closed + norb_occ != _nelec ) {
      error("Mismatch in number of electrons and OCC/CLOSED orbitals!");
    }
    if ( norb_occ - norb_closed != _ms2 ) {
      error("MS2 is not equal #OCC-#CLOSED!");
    }
      
  }
}

void Hdump::alloc_ints()
{
  assert(!allocated());
  assert(int(aaaa) == int(aa) && int(aa) == int(alpha) && int(bbbb) == int(bb) && int(bb) == int(beta));
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
    pInt->redunwarn = Input::iPars["ham"]["redunwarning"];
    #else
    pInt->redunwarn = true;
    #endif
  #endif
  FCIdump::integralType type;
  do {
    pInt->set(_osord[i-1],_osord[j-1],_osord[k-1],_osord[l-1], value);
  } while ( (type = _dump.nextIntegral(i,j,k,l,value)) == curtype ); 
  curtype = type;
  #ifdef _DEBUG
    pInt->redunwarn = false;
  #endif
}
template<typename T>
void Hdump::readrec(T* pInt, int& i, int& j, double& value, FCIdump::integralType& curtype)
{
  #ifdef _DEBUG
    #ifndef MOLPRO
    pInt->redunwarn = Input::iPars["ham"]["redunwarning"];
    #else
    pInt->redunwarn = true;
    #endif
  #endif
  FCIdump::integralType type;
  int k,l;
  do {
    pInt->set(_osord[i-1],_osord[j-1], value);
  } while ( (type = _dump.nextIntegral(i,j,k,l,value)) == curtype ); 
  curtype = type;
  #ifdef _DEBUG
    pInt->redunwarn = false;
  #endif
}
void Hdump::skiprec(int& i, int& j, int& k, int& l, double& value, FCIdump::integralType& curtype)
{
  FCIdump::integralType type;
  do {
  } while ( (type = _dump.nextIntegral(i,j,k,l,value)) == curtype ); 
  curtype = type;
}

void Hdump::gen_spinorbsref()
{
  // order of spin orbitals
  int isoord; 
  #ifndef MOLPRO
  isoord = Input::iPars["orbs"]["soorder"];
  #else
  // spin-blocked order for molpro
  isoord = 2;
  #endif
  _spinorbs.clear();
  _spinorbs.reserve(_norb*2);
  if ( _clos.size() == 0 || isoord == 0 ) {
    // default: alternating alpha/beta
    for ( uint i = 0; i < _norb; ++i ) {
      _spinorbs.emplace_back(SpinOrb(i,alpha));
      _spinorbs.emplace_back(SpinOrb(i,beta));
    }
  } else if ( isoord == 1 ) {
    // alpha/beta-alpha(open-shell) + beta(corresponding open-shell) + alpha/beta
    assert( _clos.size() == _occ.size() );
    // orbital offset for the irrep
    uint oorb = 0;
    for ( uint ir = 0; ir < _clos.size(); ++ir ) {
      for ( int i = 0; i < _clos[ir]; ++i ) {
        _spinorbs.emplace_back(SpinOrb(i+oorb,alpha));
        _spinorbs.emplace_back(SpinOrb(i+oorb,beta));
      }
      for ( int i = _clos[ir]; i < _occ[ir]; ++i ) {
        _spinorbs.emplace_back(SpinOrb(i+oorb,alpha));
      }
      for ( int i = _clos[ir]; i < _occ[ir]; ++i ) {
        _spinorbs.emplace_back(SpinOrb(i+oorb,beta));
      }
      int norbs = _pgs.norbs(ir);
      for ( int i = _occ[ir]; i < norbs; ++i ) {
        _spinorbs.emplace_back(SpinOrb(i+oorb,alpha));
        _spinorbs.emplace_back(SpinOrb(i+oorb,beta));
      }
      oorb += norbs;
    }
  } else {
    // blocked spin order: alpha(occ) + beta(occ) + alpha(virt) + beta(virt)
    assert( _clos.size() == _occ.size() );
    // orbital offset for the irrep
    uint oorb = 0;
    for ( uint ir = 0; ir < _clos.size(); ++ir ) {
      // alpha occ
      for ( int i = 0; i < _occ[ir]; ++i ) {
        _spinorbs.emplace_back(SpinOrb(i+oorb,alpha));
      }
      // beta occ
      for ( int i = 0; i < _clos[ir]; ++i ) {
        _spinorbs.emplace_back(SpinOrb(i+oorb,beta));
      }
      int norbs = _pgs.norbs(ir);
      // alpha virt
      for ( int i = _occ[ir]; i < norbs; ++i ) {
        _spinorbs.emplace_back(SpinOrb(i+oorb,alpha));
      }
      // beta virt
      for ( int i = _clos[ir]; i < norbs; ++i ) {
        _spinorbs.emplace_back(SpinOrb(i+oorb,beta));
      }
      oorb += norbs;
    }
  }
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
    if (_occ.size() > 0) _dump.addParameter("OCC",nocc_wcore());
    if (_clos.size() > 0) _dump.addParameter("CLOSED",nclos_wcore());
    if (_core.size() > 0) _dump.addParameter("CORE",_core);
    if ( _dm ) _dump.addParameter("DM",std::vector<int>(1,1));
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
    } else if ( _osord.reordered ) {
      //reorder orbitals according to the point-group symmetry
      ORBSYM_SAV = _dump.parameter("ORBSYM");
      ORBSYM = ORBSYM_SAV;
      _osord.reorder(ORBSYM);
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
  if ( ( nosym || _osord.reordered ) && !newfile ) {
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
        pDest->set(indxD, 0.25*(pSrc->get(i,j,k,l)+pSrc->get(j,i,l,k)+pSrc->get(j,i,k,l)+pSrc->get(i,j,l,k))+pDest->get(indxD));
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
        pDest->set(indxD, 0.25*(pSrc->get(i,j,k,l)+pSrc->get(j,i,l,k)+pSrc->get(j,i,k,l)+pSrc->get(i,j,l,k)));
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
void Hdump::add_or_subtract_core(FDPar& orb, const std::string& kind, bool add) const 
{
  for ( uint io = 0; io < std::min(orb.size(),_core.size()); ++io ) {
    if ( add ) {
      orb[io] += _core[io];
    } else if ( orb[io] < _core[io] ) {
      xout << io << ") " << kind << ": " << orb[io] << " CORE: " << _core[io] << std::endl;
      error("Core orbitals should be included in the "+kind+" specification!");
    } else {
      orb[io] -= _core[io];
    }
  }
  for ( uint io = orb.size(); io < _core.size(); ++io ) {
    if ( add ) {
      orb.push_back(_core[io]);
    } else if ( _core[io] != 0 ) {
      error("Core orbitals should be included in the "+kind+" specification!");
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
  add_or_subtract_core(clco,"CLOSED",true);
  return clco; 
}
FDPar Hdump::nocc_wcore() const
{
  FDPar occo(_occ); 
  add_or_subtract_core(occo,"OCC",true);
  return occo; 
}
void Hdump::set_noccorbs(const FDPar& core, const FDPar& closed, const FDPar& occ, bool wcore)
{
  assert(closed.size() <= _pgs.nIrreps());
  assert(occ.size() <= _pgs.nIrreps());
  // core can have more irreps than the hamdump
  assert(core.size() <= 8);
  _core.clear();
  if (!core.empty()){
    FDPar tmpcore(core);
    tmpcore.resize(8,0);
    set_ncore(&tmpcore[0]);
  }
  _clos = closed;
  _occ = occ;
  if ( wcore ) {
    // remove core
    add_or_subtract_core(_clos,"CLOSED",false);
    add_or_subtract_core(_occ,"OCC",false);
  }
  gen_spinorbsref();
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

template<typename T, typename U, typename V>
void Hdump::add1RDMto2RDM( T * pD2, const U * pD1, const V * pD11, double fact, bool exchange)
{
  Irrep isym = 0;
  // spatial orbitals
  uint p = 0, q = 0, r = 0, s = 0;
  double val;
  do { 
    BlkIdx indx = pD2->index(p,q,r,s);
    val = pD2->get(indx);
    if ( _pgs.totIrrep(p, q) == 0 ) {
      val += fact * pD1->get(p,q) * pD11->get(r,s);
    }
    if ( exchange && _pgs.totIrrep(p, s) == 0 ) {
      val -= fact * pD1->get(p,s) * pD11->get(r,q);
    }
    pD2->set(indx,val);
  } while (pD2->next_indices(p,q,r,s,isym));
  
}

void Hdump::correct_2DM()
// adds to the 2DM the (**|ll) parts etc. which where omitted in ITF StUccsd.itfaa
{
  assert(_simtra);
  Integ2st * pDa = dynamic_cast<Integ2st*>(_oneel[aa].get()); assert(pDa);
  Integ2st * pDb = dynamic_cast<Integ2st*>(_oneel[bb].get()); assert(pDb);
  Integ4st * pDaa = dynamic_cast<Integ4st*>(_twoel[aaaa].get()); assert(pDaa);
  Integ4st * pDbb = dynamic_cast<Integ4st*>(_twoel[bbbb].get()); assert(pDbb);
  Integ4stab * pDab = dynamic_cast<Integ4stab*>(_twoel[aabb].get()); assert(pDab);
  add1RDMto2RDM(pDaa,pDa,pDa,-1.0,true);
  add1RDMto2RDM(pDbb,pDb,pDb,-1.0,true);
  add1RDMto2RDM(pDab,pDa,pDb,-1.0,false);
  
  double val;
  for( uint ir = 0; ir < _pgs.nIrreps(); ir++ ){ 
    uint sporb = _pgs._firstorb4irrep[ir] * 2;
    for( int iorb = 0; iorb < _clos[ir]; iorb++ ){
      val = oneel_spi(sporb, sporb);
      val += 1.0;
      set_oneel_spi(sporb, sporb, val);
      sporb += 1;
      val = oneel_spi(sporb, sporb);
      val += 1.0;
      set_oneel_spi(sporb, sporb, val);
      sporb += 1;
    }
    for( int iorb = _clos[ir]; iorb < _occ[ir]; iorb++ ){
        val = oneel_spi(sporb, sporb);
        val += 1;
        set_oneel_spi(sporb, sporb, val);
        sporb += 1;
    }
  }
  assert(_simtra);
  add1RDMto2RDM(pDaa,pDa,pDa,1.0,true);
  add1RDMto2RDM(pDbb,pDb,pDb,1.0,true);
  add1RDMto2RDM(pDab,pDa,pDb,1.0,false);
}

void Hdump::calc_Fock(const Hdump& DMmats)
{
  std::vector< std::vector<double> >_FMAT; // fock matrix for each symmetry
  uint _nsorb = 2*_norb;
  _FMAT.resize(_pgs.nIrreps());
  // energy from the density matrices
  double Tr = 0.0;
  for ( uint ir = 0; ir < _pgs.nIrreps(); ++ir ) {
    uint nsorb4ir = _pgs.norbs(ir) * 2;
    uint isoff4ir = _pgs._firstorb4irrep[ir] * 2;
    uint maxid = nsorb4ir-1;
    uint nelem = maxid+maxid*nsorb4ir+1;
    std::vector<double> & fmat = _FMAT[ir];
    fmat.resize(nelem,0);
    std::vector<double> fmat2(nelem,0);
     
    for(uint p4ir = 0; p4ir < nsorb4ir; p4ir++){
      uint p = p4ir+isoff4ir;
      for(uint t4ir = 0; t4ir < nsorb4ir; t4ir++){
        uint t = t4ir+isoff4ir;
        uint pt = p4ir+t4ir*nsorb4ir;
        for(uint q4ir = 0; q4ir < nsorb4ir; q4ir++){
          uint q = q4ir+isoff4ir;
          fmat[pt] += oneel_spi(t,q) * DMmats.oneel_spi(p,q); 
        }
      }
    }
    std::vector<double> dm_qrs(nsorb4ir);
    
    for(uint q = 0; q < _nsorb; q++){
      for(uint r = 0; r < _nsorb; r++){
        for(uint s = 0; s < _nsorb; s++){
          for(uint t4ir = 0; t4ir < nsorb4ir; t4ir++){
            uint t = t4ir+isoff4ir;
            //ATTENTION hermitian symmetrisation only for comparison. TS
//             dm_qrs[t4ir] = (DMmats.twoel_spi_pgs(t,r,q,s) + DMmats.twoel_spi_pgs(r,t,s,q))/2.0;
            dm_qrs[t4ir] = DMmats.twoel_spi_pgs(t,r,q,s);
          }
          for(uint p4ir = 0; p4ir < nsorb4ir; p4ir++){
            uint p = p4ir+isoff4ir;
            double twoel_prqs = twoel_spi_pgs(p,r,q,s); 
            for(uint t4ir = 0; t4ir < nsorb4ir; t4ir++){
              fmat2[t4ir+p4ir*nsorb4ir] += twoel_prqs * dm_qrs[t4ir];
            }
          }
        }
      }
    }
//     // add two electron parts to fock
    for(uint pt = 0; pt < fmat.size(); pt++){
      fmat[pt] += fmat2[pt]; 
    }
    
    for(uint p = 0; p < nsorb4ir; p++){
      uint pt = p+p*nsorb4ir;
      Tr = Tr + fmat[pt]-0.5*fmat2[pt];
    }
  }
  xout << "Energy: " << std::setprecision (15) << Tr+escal() << std::endl;
}

void Hdump::molcas_2DM(bool plus){
    if( _pgs.nIrreps() > 1 ) error( "Molcas format for 2DM only implemented for no symmetry." );
    std::ofstream cDMfile;
    cDMfile.open("/home/schraivogel/cDM_dummy.dat");
    int cDMsize = _norb * _norb * _norb * _norb; //TODO cDMsize will definitively be smaller than that..
    std::vector<double> cDM(cDMsize); //compressed density matrix
    cDMfile << std::setprecision(10);
    double symdiff = 0;
    int tr, qs, trqs;
    for(uint t = 0; t < _norb; t++){
      for(uint r = 0; r < _norb; r++){
        for(uint q = 0; q < _norb; q++){
          for(uint s = 0; s < _norb; s++){
            if( t >= r ){
              tr = t*(t+1)/2 + r;
              if( q > s ){
                qs = q*(q+1)/2 +s;
                  if( tr >= qs ){
                    trqs = tr*(tr+1)/2 + qs;
                    if ( plus )
                      cDM[trqs] = (  ( spinsum_2DM(t,r,q,s) + spinsum_2DM(r,t,s,q) )/2.0 
                                    +( spinsum_2DM(t,r,s,q) + spinsum_2DM(r,t,q,s) )/2.0
                                  )/2.0;
                    else
                      cDM[trqs] = (  ( spinsum_2DM(t,r,q,s) + spinsum_2DM(r,t,s,q) )/2.0
                                    -( spinsum_2DM(t,r,s,q) + spinsum_2DM(r,t,q,s) )/2.0
                                  )/2.0;
                    symdiff += std::abs(twoel_spi_pgs(t,r,q,s) - twoel_spi_pgs(r,t,s,q));
                    if ( std::abs(cDM[trqs]) > 1.e-24 )
                      cDMfile << std::setw(20) << cDM[trqs] << " ";
                    else
                      cDMfile << std::setw(20)  << "0.0" << " ";
                  }
              }
              else if( q == s ){
                qs = q*(q+1)/2 +s;
                if( tr >= qs ){
                  trqs = tr*(tr+1)/2 + qs; 
                  symdiff += std::abs(twoel_spi_pgs(t,r,q,s) - twoel_spi_pgs(r,t,s,q));
                    cDM[trqs] = (spinsum_2DM(t,r,q,s)+spinsum_2DM(r,t,s,q))/4.0;
                  if ( std::abs(cDM[trqs]) > 1.e-24 )
                    cDMfile << std::setw(20) << cDM[trqs] << " ";
                  else
                    cDMfile << std::setw(20) << "0.0" << " ";
                }
              }
            }
          }
      }
        cDMfile << std::endl;
    }
    }
    cDMfile.close();
    //delete empty lines in cDM file
    std::ifstream inFile;
    inFile.open("/home/schraivogel/cDM_dummy.dat");
    std::ofstream outFile;
    outFile.open("/home/schraivogel/cDM.dat");
    std::string line;
    while( std::getline( inFile, line) ){
      if( ! line.empty() ){
        outFile << line << std::endl;
      }
    }
    inFile.close();
    outFile.close();
    
    xout << "A compressed 2RDM file called cDM.dat has been written to Thomas' home directory." << std::endl;
    xout << "Hermitian symmetry badness of 2RDM: " << symdiff << std::endl;
    remove( "/home/schraivogel/cDM_dummy.dat" );
}

std::size_t skip(const std::string& str, std::size_t ipos, const std::string& what) 
{
  while (ipos<str.size()&& what.find(str[ipos])!=std::string::npos )
    ++ipos;
  return ipos;
}
std::vector<std::string> split(const std::string& str, const std::string& listsep)
{
  std::vector<std::string> res;
  std::size_t 
    ipos = skip(str,0,listsep),
    ipend = str.find_first_of(listsep,ipos+1);
  
  while( ipend != ipos ){
    res.push_back(str.substr(ipos,ipend-ipos));
    ipos = ipend;
    if ( ipos < str.size() ) ++ipos;
    ipos = skip(str,ipos,listsep);
    ipend =str.find_first_of(listsep,ipos+1);
    if (ipend == std::string::npos ) ipend = str.size();
  }
  return res;
}
template<typename T>
void str2numvec(std::vector<T>& numvec, const std::vector<std::string>& strvec, uint startpos=0) {
#ifdef MOLPRO
  using namespace fmt;
#endif
  for ( uint pos = startpos; pos < strvec.size(); ++pos ) {
    T val;
    if ( str2num<T>(val,strvec[pos],std::dec) ) {
      numvec.push_back(val);
    } else
      error("Not a number: "+strvec[pos]);
  }
}

Overlapdump::Overlapdump(std::string ovdump, const PGSym& pgs, 
                         const FDPar& ncore, bool verbose)
{
  overlap = Integ2ab(pgs);
  uint nIrreps = pgs.nIrreps();
  std::ifstream oin;
    oin.open(ovdump.c_str());
  if ( !oin.is_open() ) {
    error("Error opening "+ ovdump);
  }
  std::string line;
  oin >> std::ws;
  std::getline(oin,line);
  if ( line.compare(0,3,"BAS") != 0 ) 
    error("No basis information in file "+ovdump);
  std::vector<int> nbasis;
  std::vector<std::string> strvec = split(line," ,;:");
  str2numvec(nbasis,strvec,1);
  if (verbose) {
    xout << "Basis dimensions of the overlap: ";
    for(auto bb: nbasis ) xout << bb << " ";
    xout << std::endl; 
  }
  // TODO sanity check for nbasis 
  
  
  std::vector<double> numbers;
  std::vector<double> pnumbers;
  uint p;
  for( uint isym = 0; isym < nIrreps; isym++ ){
    numbers.reserve(nbasis[isym]);
    pnumbers.reserve(nbasis[isym]);
    for( int row = 0; row < nbasis[isym]; row++ ){
      while (oin.good()) {
        std::getline(oin,line);
        if ( line.empty() || line[0] == '#' ||  
            line.compare("BEGIN_DATA,") == 0 || line.compare("END_DATA,") == 0  ) continue;
        strvec = split(line," ,;:");
        numbers.clear();
        str2numvec(numbers,strvec);
        pnumbers.insert(std::end(pnumbers), std::begin(numbers), std::end(numbers));
        if( pnumbers.size() > unsigned(nbasis[isym]) ) {
            xout << "Unexpected format of " << ovdump << std::endl;
            error("Basis dimensions do not correspond to the number of entries in a row");
        }
        if( pnumbers.size() == unsigned(nbasis[isym]) ) break;
      }   
      if( row >= ncore[isym] && unsigned(row) < ncore[isym] + pgs.norbs(isym) ){
        p = row - ncore[isym];
        for( uint q = 0; q < pgs.norbs(isym); q++ ){
          overlap.set(p,q,isym,pnumbers[q+ncore[isym]]);
        }
      }
      pnumbers.clear();
    }
//     xout << "Overlap, irrep=" << isym << std::endl;
//     for ( uint p = 0; p < pgs.norbs(isym); ++p ) {
//       for ( uint q = 0; q < pgs.norbs(isym); ++q ) {
//         xout << fmt::ff(overlap.get(p,q,isym),15,8);
//       }
//       xout << std::endl;
//     }
      
  }
}
void Hdump::calc_Spin2(const Integ2ab& Overlap){
  double spin = 0;
  double spinHF = 0;
  double ms = 0;
  for( uint p = 0; p < _norb; p++ ){
    ms += _oneel[aa].get()->get_with_pgs(p,p);
    ms -= _oneel[bb].get()->get_with_pgs(p,p);
    for( uint q = 0; q < _norb; q++ ){
      for( uint s = 0; s < _norb; s++ ){
        spin += Overlap.get_with_pgs(q,p) * Overlap.get_with_pgs(q,s) * _oneel[bb].get()->get_with_pgs(p,s);
        for( uint r = 0; r < _norb; r++ ){
          spin  -= Overlap.get_with_pgs(q,p) * Overlap.get_with_pgs(r,s) * _twoel[aabb].get()->get_with_pgs(r,q,p,s);
        }
      }
    }
  }
  spin += 0.5*ms * ( 0.5*ms + 1 );
  if( _pgs.nIrreps() > 1 ) xout << "ATTENTION: <0|S**2|0> only implemented for nosym." << std::endl;
  for( int ibeta = 0; ibeta < _clos[0]; ibeta++ ){
    for( uint aalpha = _occ[0]; aalpha < _norb; aalpha++ ){
      spinHF += Overlap.get_with_pgs(aalpha,ibeta) * Overlap.get_with_pgs(aalpha,ibeta);
    }
  }
  spinHF += 0.5*ms * ( 0.5*ms + 1 );
    
  xout << "<0|S**2|0>: " << spinHF << std::endl;
  xout << "Expectation value of S**2: " << spin << std::endl;
}
} //namespace HamDump
