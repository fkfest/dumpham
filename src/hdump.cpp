#include "hdump.h"

Hdump::Hdump(std::string fcidump) : _dump(fcidump)
{
  xout << "\n Process file "<< fcidump <<std::endl;
  FDPar NELEC = _dump.parameter("NELEC");
  xout << "NELEC=" << NELEC[0] << std::endl;
  FDPar MS2 = _dump.parameter("MS2");
  xout << "MS2=" << MS2[0] << std::endl;
  FDPar NORB = _dump.parameter("NORB");
  xout << "NORB=" << NORB[0] << std::endl;
  FDPar ISYM = _dump.parameter("ISYM");
  xout << "ISYM=" << ISYM[0] << std::endl;
  FDPar IUHF = _dump.parameter("IUHF");
  xout << "IUHF=" << IUHF[0] << std::endl;
  FDPar ORBSYM = _dump.parameter("ORBSYM");
  xout << "ORBSYM="; 
  _foreach_cauto(FDPar,s,ORBSYM)
    xout << *s << ","; 
  xout<<std::endl;
  FDPar OCC = _dump.parameter("OCC");
  xout << "OCC=";
  _foreach_cauto(FDPar,s,OCC)
    xout << *s << ","; 
  xout<<std::endl;
  FDPar CLOSED = _dump.parameter("CLOSED");
  xout << "CLOSED=";
  _foreach_cauto(FDPar,s,CLOSED)
    xout << *s << ","; 
  xout<<std::endl;
  FDPar ST = _dump.parameter("ST");
  
  int i,j,k,l;
  double value;
  int nn = NORB[0];
  if ( nn < 0 ) {
    error("NORB < 0 in FCIDUMP!");
  }
  _norb = nn;
  _nelec = NELEC[0];
  _ms2 = MS2[0];
  _pgs = PGSym(ORBSYM);
  if ( OCC.size() > 1 )
    _occ = OCC;
  if ( CLOSED.size() > 1 )
    _clos = CLOSED;
  _uhf = bool(IUHF[0]);
  _simtra = bool(ST[0]);
  
  _oneel.emplace_back(new Integ2(_pgs));
  _twoel.emplace_back(new Integ4(_pgs));
  if ( _uhf ) {
    _oneel.emplace_back(new Integ2(_pgs)); // h_bb
    _twoel.emplace_back(new Integ4(_pgs)); // (bb|bb)
    _twoel.emplace_back(new Integ4ab(_pgs)); // (aa|bb)
  }
  
#ifdef _DEBUG
//   check_addressing_integrals();
#endif
  _escal = 0.0;
  
  FCIdump::integralType type;
  _dump.rewind();
  
  type = _dump.nextIntegral(i,j,k,l,value);
  do {
    switch(type){
      case FCIdump::I2aa:
        readrec( static_cast<Integ4*> (_twoel[aaaa].get()),i,j,k,l,value,type);
        break;
      case FCIdump::I2bb:
        readrec( static_cast<Integ4*> (_twoel[bbbb].get()),i,j,k,l,value,type);
        break;
      case FCIdump::I2ab:
        readrec( static_cast<Integ4ab*> (_twoel[aabb].get()),i,j,k,l,value,type);
        break;
      case FCIdump::I1a:
        readrec( static_cast<Integ2*> (_oneel[aa].get()),i,j,value,type);
        break;
      case FCIdump::I1b:
        readrec( static_cast<Integ2*> (_oneel[bb].get()),i,j,value,type);
        break;
      case FCIdump::I0:
        _escal = value;
        [[fallthrough]];
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
    bool redunwarn = Input::iPars["ham"]["redunwarning"];
    #define REDUNWAR ,redunwarn
  #else
    #define REDUNWAR
  #endif
  FCIdump::integralType type;
  do {
    pInt->value(i-1,j-1,k-1,l-1 REDUNWAR) = value;
  } while ( (type = _dump.nextIntegral(i,j,k,l,value)) == curtype ); 
  curtype = type;
  #undef REDUNWAR
}
template<typename T>
void Hdump::readrec(T* pInt, int& i, int& j, double& value, FCIdump::integralType& curtype)
{
  #ifdef _DEBUG
    bool redunwarn = Input::iPars["ham"]["redunwarning"];
    #define REDUNWAR ,redunwarn
  #else
    #define REDUNWAR
  #endif
  FCIdump::integralType type;
  int k,l;
  do {
    pInt->value(i-1,j-1 REDUNWAR) = value;
  } while ( (type = _dump.nextIntegral(i,j,k,l,value)) == curtype ); 
  curtype = type;
  #undef REDUNWAR
}

void Hdump::store(std::string fcidump)
{
  bool nosym = Input::iPars["ham"]["nosym"];
  FDPar ORBSYM, ORBSYM_SAV ;
  if ( nosym ) {
    //remove the point-group symmetry
    ORBSYM_SAV = _dump.parameter("ORBSYM");
    ORBSYM.resize(ORBSYM_SAV.size());
    _foreach_auto(FDPar,s,ORBSYM)
      *s = 1;
    _dump.modifyParameter("ORBSYM",ORBSYM);
  }
  if (_dump.write(fcidump,FCIdump::FileFormatted,false))
    xout << "will be written to file " << fcidump << std::endl;
  else
    error("failure to write to file "+fcidump);
  
  if (nosym) {
    store_without_symmetry();
  } else {
    store_with_symmetry();
  }
  _dump.close_outputfile();
  if ( nosym ) {
    //restore the point-group symmetry
    _dump.modifyParameter("ORBSYM",ORBSYM_SAV);
  }
}

void Hdump::store_with_symmetry()
{
  storerec_sym(static_cast<Integ4*> (_twoel[aaaa].get()));
  if (_uhf) {
    _dump.writeIntegral(0,0,0,0,0.0);
    storerec_sym(static_cast<Integ4*> (_twoel[bbbb].get()));
    _dump.writeIntegral(0,0,0,0,0.0);
    storerec_sym(static_cast<Integ4ab*> (_twoel[aabb].get()));
    _dump.writeIntegral(0,0,0,0,0.0);
  }
  storerec_sym(static_cast<Integ2*> (_oneel[aa].get()));
  if (_uhf) {
    _dump.writeIntegral(0,0,0,0,0.0);
    storerec_sym(static_cast<Integ2*> (_oneel[bb].get()));
    _dump.writeIntegral(0,0,0,0,0.0);
  }
  _dump.writeIntegral(0,0,0,0,_escal); 
}

void Hdump::store_without_symmetry()
{
  storerec_nosym(static_cast<Integ4*> (_twoel[aaaa].get()));
  if (_uhf) {
    _dump.writeIntegral(0,0,0,0,0.0);
    storerec_nosym(static_cast<Integ4*> (_twoel[bbbb].get()));
    _dump.writeIntegral(0,0,0,0,0.0);
    storerec_nosym(static_cast<Integ4ab*> (_twoel[aabb].get()));
    _dump.writeIntegral(0,0,0,0,0.0);
  }
  storerec_nosym(static_cast<Integ2*> (_oneel[aa].get()));
  if (_uhf) {
    _dump.writeIntegral(0,0,0,0,0.0);
    storerec_nosym(static_cast<Integ2*> (_oneel[bb].get()));
    _dump.writeIntegral(0,0,0,0,0.0);
  }
  _dump.writeIntegral(0,0,0,0,_escal); 
}

void Hdump::storerec_sym(Integ4* pInt)
{
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
              _dump.writeIntegral(i,j,k,l,pInt->value(i-1,j-1,k-1,l-1));
            }
          }
        }
      }
    }
  }
}
void Hdump::storerec_sym(Integ4ab* pInt)
{
  for ( Irrep isym = 0; isym < _pgs.nIrreps(); ++isym ) {
    for ( uint i = 1; i <= _norb; ++i) {
      for ( uint j = 1; j <= i; ++j ) {
        if ( _pgs.totIrrep(i-1,j-1) != isym ) continue;
        for ( uint k = 1; k <= i; ++k ) {
          for ( uint l = 1; l <= k; ++l ) {
            if ( _pgs.totIrrep(k-1,l-1) != isym ) continue;
              _dump.writeIntegral(i,j,k,l,pInt->value(i-1,j-1,k-1,l-1));
          }
        }
      }
    }
  }
}
void Hdump::storerec_sym(Integ2* pInt)
{
  for ( uint i = 1; i <= _norb; ++i) {
    for ( uint j = 1; j <= i; ++j ) {
        if ( _pgs.totIrrep(i-1,j-1) != 0 ) continue;
        _dump.writeIntegral(i,j,0,0,pInt->value(i-1,j-1));
    }
  }
}

void Hdump::storerec_nosym(Integ4* pInt)
{
  for ( uint i = 1; i <= _norb; ++i)
    for ( uint j = 1; j <= i; ++j )
      for ( uint k = 1; k <= _norb; ++k )
        for ( uint l = 1; l <= k; ++l ) {
          uint ij = (i-1)*i/2 + j;
          uint kl = (k-1)*k/2 + l;
          if ( kl <= ij ){
//             int ijkl = (ij-1)*ij/2 + kl - 1;
            if ( _pgs.totIrrep(i-1,j-1,k-1,l-1) == 0 ) {
              _dump.writeIntegral(i,j,k,l,pInt->value(i-1,j-1,k-1,l-1));
            } else {
              _dump.writeIntegral(i,j,k,l,0.0);
            }
          }
        }
}
void Hdump::storerec_nosym(Integ4ab* pInt)
{
  for ( uint i = 1; i <= _norb; ++i)
    for ( uint j = 1; j <= i; ++j )
      for ( uint k = 1; k <= _norb; ++k )
        for ( uint l = 1; l <= k; ++l ) {
          if ( _pgs.totIrrep(i-1,j-1,k-1,l-1) == 0 ) {
            _dump.writeIntegral(i,j,k,l,pInt->value(i-1,j-1,k-1,l-1));
          } else {
            _dump.writeIntegral(i,j,k,l,0.0);
          }
        }
}
void Hdump::storerec_nosym(Integ2* pInt)
{
  for ( uint i = 1; i <= _norb; ++i)
    for ( uint j = 1; j <= i; ++j ) {
      if ( _pgs.totIrrep(i-1,j-1) == 0 ) {
        _dump.writeIntegral(i,j,0,0,pInt->value(i-1,j-1));
      } else {
        _dump.writeIntegral(i,j,0,0,0.0);
      }
    }
}


uint Hdump::nclosed() const
{
  uint nclos = (_nelec - _ms2)/2;
  if ( nclos*2 + _ms2 != _nelec ) {
    xout << "NELEC: " << _nelec << " MS2: " << _ms2 << std::endl;
    error("Mismatch in NELEC and MS2 in Hdump");
  }
  return nclos;
}

void Hdump::check_addressing_integrals() const
{
  xout << "test one el" << std::endl;
  int oldidx=-1;
  for (uint i = 1; i <= _norb; ++i )
    for (uint j = 1; j <= i; ++j ) {
      if ( _pgs.totIrrep(i-1,j-1) != 0 ) continue;
      BlkIdx idx = (static_cast<Integ2*>(_oneel[aa].get()))->index(i-1,j-1);
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
              BlkIdx idx = (static_cast<Integ4*>(_twoel[aaaa].get()))->index(i-1,j-1,k-1,l-1);
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

