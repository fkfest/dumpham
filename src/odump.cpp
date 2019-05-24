#include "odump.h"

void Occupation4Irrep::spinocc(std::vector<int>& socc, int ibase) const
{
  for ( uint imo = 0; imo < this->size(); ++imo ){
    socc.push_back((*this)[imo]*2+ibase);
    if ( imo < _nclos )
      socc.push_back((*this)[imo]*2+1+ibase);
  }
}

Occupation::Occupation(const PGSym& pgs, const FDPar& nclos, const FDPar& nocc)
    : p_pgs(&pgs)
{
  assert(nocc.size() >= p_pgs->nIrreps() && nclos.size() >= p_pgs->nIrreps());
  for ( uint ir = 0; ir < p_pgs->nIrreps(); ++ir ) {
    Occupation4Irrep o4ir(p_pgs->_firstorb4irrep[ir], nocc[ir]);
    o4ir._nclos = nclos[ir];
    push_back(o4ir);
  }
}


Occupation::Occupation(const PGSym& pgs, const std::vector< int >& occs, int ibase)
    : p_pgs(&pgs)
{
  resize(p_pgs->nIrreps());
  std::vector< std::vector<int> > singly(p_pgs->nIrreps());
  // For 1-based: spatial obital index = int( (spin orbital index + 1) / 2 )
  int iorb0 = -1;
  for ( const auto& iocc: occs ){
    int iorb1 = int((iocc-ibase)/2);
    if ( iorb0 == iorb1 ) {
      // doubly occupied
      (*this)[p_pgs->irrep(iorb0)].push_back(iorb0);
      iorb0 = -1;
    } else if ( iorb0 != -1 ) {
      // iorb0 is singly occupied
      singly[p_pgs->irrep(iorb0)].push_back(iorb0);
      iorb0 = iorb1;
    } else {
      iorb0 = iorb1;
    }
  }
  if ( iorb0 != -1 ) {
    // iorb0 is singly occupied
    singly[p_pgs->irrep(iorb0)].push_back(iorb0);
  }
  for ( uint ir = 0; ir < p_pgs->nIrreps(); ++ir ) {
    (*this)[ir]._nclos = (*this)[ir].size();
    // add singly occupied orbitals
    for ( const auto& iocc: singly[ir] ){
      (*this)[ir].push_back(iocc);
    }
  }
}

std::vector<int> Occupation::spinocc(int ibase) const
{
  std::vector<int> socc;
  for ( const auto& oc: *this)
    oc.spinocc(socc,ibase);
  std::sort(socc.begin(),socc.end());
  return socc;
}


std::ostream& operator<<(std::ostream& o, const Occupation& occ)
{
  for ( const auto& iocc: occ ){
    for ( uint i = 0; i < iocc.size(); ++i ) {
      if ( i < iocc._nclos )
        o << iocc[i] << ",";
      else
        o << iocc[i] << "(o),";
    }
  }
  return o;
}


Odump::Odump(const PGSym& pgs, Occupation occs ) : p_pgs(&pgs)
{
  _orbs = Integ2ab(pgs);
  _ncore.resize(p_pgs->nIrreps(),0);
  _ncoreaccu.resize(p_pgs->nIrreps(),0);
  if ( occs.empty() ) occs.resize(p_pgs->nIrreps()); // trivial occupation: all orbitals "empty"
  uint iao;
  xout << "Spatial orbital occupation: " << occs << std::endl;
  assert( occs.size() == p_pgs->nIrreps() );
  for ( uint ir = 0; ir < p_pgs->nIrreps(); ++ir ) {
    std::vector<bool> occupied(p_pgs->norbs(ir),false);
    for ( uint imo4ir = 0; imo4ir < occs[ir].size(); ++imo4ir ) {
      uint imo = imo4ir + p_pgs->_firstorb4irrep[ir];
      iao = occs[ir][imo4ir];
      _orbs.set(iao,imo,1.0);
      assert( iao >= uint(p_pgs->_firstorb4irrep[ir]) && iao < uint(p_pgs->_firstorb4irrep[ir]+p_pgs->norbs(ir)));
      occupied[iao-p_pgs->_firstorb4irrep[ir]] = true;
    }
    uint imo = occs[ir].size()+p_pgs->_firstorb4irrep[ir]; 
    for ( uint iao4ir = 0; iao4ir < p_pgs->norbs(ir); ++iao4ir ) {
      if (!occupied[iao4ir]) {
        iao = iao4ir + p_pgs->_firstorb4irrep[ir];
        _orbs.set(iao,imo,1.0);
        ++imo;
      }
    }
    assert( imo == p_pgs->norbs(ir) + p_pgs->_firstorb4irrep[ir] );
    _nclos.push_back(occs[ir]._nclos);
    _nocc.push_back(occs[ir].size());
  }
}
Odump::Odump(const PGSym& pgs, const FDPar& ncore, std::string orbdump) : p_pgs(&pgs)
{
  
  _orbs = Integ2ab(pgs);
  _ncore = ncore;
  assert( _ncore.size() == p_pgs->nIrreps() );
  _ncoreaccu.push_back(_ncore[0]);
  for ( Irrep ir = 1; ir < _ncore.size(); ++ir )
    _ncoreaccu.push_back(_ncoreaccu[ir-1]+_ncore[ir]);
//   xout << "orbs nelem: " << _orbs.nelem() << std::endl;
  std::ifstream oin;
  oin.open(orbdump.c_str());
  if ( !oin.is_open() ) {
    error("Error opening "+ orbdump);
  }
  std::string line;
  BlkIdx idx = 0;
  while (oin.good()) {
    std::getline(oin,line);
    if ( line.empty() || line[0] == '#' || 
         line.compare("BEGIN_DATA,") == 0 || line.compare("END_DATA,") == 0  ) continue;
    TParArray coefs_line = IL::parray(line);
    for ( const auto& ic: coefs_line ){
      double val;
      if ( str2num<double>(val,ic,std::dec) ) {
        _orbs.set(idx,val);
        ++idx;
      } else
        error("Not a number: "+ic);
    }
  }
  if ( idx != _orbs.nelem() )
    error("Number of orbital coefficients not consistens with the number of orbitals");
}
Odump::Odump(const PGSym& pgs, const FDPar& ncore, const Integ2ab& orbs) : p_pgs(&pgs)
{
  _ncore = ncore;
  assert( _ncore.size() == p_pgs->nIrreps() );
  _ncoreaccu.push_back(_ncore[0]);
  for ( Irrep ir = 1; ir < _ncore.size(); ++ir )
    _ncoreaccu.push_back(_ncoreaccu[ir-1]+_ncore[ir]);
  if ( orbs.pgs() != p_pgs ){
    assert( orbs.pgs()->nIrreps() == p_pgs->nIrreps() );
    _orbs = Integ2ab(pgs);
    for ( Irrep ir = 0; ir < p_pgs->nIrreps(); ++ir ){
      uint norb4ir_in = orbs.pgs()->norbs(ir);
      uint ioff4ir_in = orbs.pgs()->_firstorb4irrep[ir];
      uint ioff4ir = p_pgs->_firstorb4irrep[ir]+_ncore[ir];
      for ( uint i = 0; i < norb4ir_in; ++i ){
        for ( uint j = 0; j < norb4ir_in; ++j ){
          _orbs.set(i+ioff4ir,j+ioff4ir,orbs.get(i+ioff4ir_in,j+ioff4ir_in));
        }
      }
    }
  } else {
    _orbs = orbs;
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
  xout << "will be written to file " << orbdump << std::endl;
  bool scientific = Input::iPars["output"]["scientificoef"];
  int precision = Input::iPars["output"]["precisioncoef"];
  if ( scientific ) 
    outputStream << std::scientific;
  else
    outputStream << std::fixed;
  outputStream << std::setprecision(precision);
  int maxlen = Input::iPars["output"]["maxncoef"];
  bool nosym = Input::iPars["ham"]["nosym"];
  if (nosym) {
    if ( p_pgs->nIrreps() <= 1 ) {
      for ( uint i = 0; i < p_pgs->ntotorbs(); ++i ) {
        for ( uint j = 0; j < p_pgs->ntotorbs(); ++j ) {
          double val = _orbs.get_with_pgs(i,j);
          printval(outputStream,val,j,maxlen,scientific);
        } 
        outputStream << std::endl;
      }
    } else {
      if ( _nclos.size() != p_pgs->nIrreps() || _nocc.size() != p_pgs->nIrreps() ) 
        error("Cannot de-symmetrize orbitals without occupation information! Something to be implemented...");
      for ( uint i = 0; i < p_pgs->ntotorbs(); ++i ) {
        // closed-shell
        uint norb = 0;
        for ( Irrep ir = 0; ir < p_pgs->nIrreps(); ++ir ) {
          for ( uint j = p_pgs->_firstorb4irrep[ir]; j < uint(p_pgs->_firstorb4irrep[ir] + _nclos[ir]); ++j ) {
            double val = _orbs.get_with_pgs(i,j);
            printval(outputStream,val,norb,maxlen,scientific);
            ++norb;
          }
        }
        // open-shell
        for ( Irrep ir = 0; ir < p_pgs->nIrreps(); ++ir ) {
          for ( uint j = p_pgs->_firstorb4irrep[ir]+_nclos[ir]; j < uint(p_pgs->_firstorb4irrep[ir] + _nocc[ir]); ++j ) {
            double val = _orbs.get_with_pgs(i,j);
            printval(outputStream,val,norb,maxlen,scientific);
            ++norb;
          }
        }
        // virtual
        for ( Irrep ir = 0; ir < p_pgs->nIrreps(); ++ir ) {
          for ( uint j = p_pgs->_firstorb4irrep[ir]+_nocc[ir]; j < uint(p_pgs->_firstorb4irrep[ir] + p_pgs->_norb4irrep[ir]); ++j ) {
            double val = _orbs.get_with_pgs(i,j);
            printval(outputStream,val,norb,maxlen,scientific);
            ++norb;
          }
        }
        outputStream << std::endl;
      }
    }
  } else {
    for ( Irrep ir = 0; ir < p_pgs->nIrreps(); ++ir ) {
      for ( uint i = p_pgs->_firstorb4irrep[ir]; i < uint(p_pgs->_firstorb4irrep[ir] + p_pgs->_norb4irrep[ir]); ++i ) {
        for ( uint j = p_pgs->_firstorb4irrep[ir]; j < uint(p_pgs->_firstorb4irrep[ir] + p_pgs->_norb4irrep[ir]); ++j ) {
          double val = _orbs.get(i,j);
          printval(outputStream,val,j-p_pgs->_firstorb4irrep[ir],maxlen,scientific);
        } 
        outputStream << std::endl;
      }
    }
  }
}
void Odump::printval(std::ofstream& outputStream, double val, uint j, uint maxlen, bool scientific)
{
  if ( j > 0 && maxlen > 0 && j % maxlen == 0 ) outputStream << std::endl;
  if ( !scientific && val <= 0.0 && val >= -1.e-20 ) val = 1.e-20;
  if ( !scientific && val > -1.e-20 ) outputStream << " ";
  outputStream << val << ", ";
}

Occupation Odump::guess_occupation(const FDPar& nclos, const FDPar& nocc) const
{
  double throcc = Input::fPars["orbs"]["occuwarning"];
  bool sortocc = Input::iPars["orbs"]["occusort"];
  Occupation occ(*p_pgs);
  FDPar n_occ(nocc);
  if ( n_occ.empty()) {
    n_occ = p_pgs->norbs_in_irreps();
  }
  xout << n_occ.size() << " " << p_pgs->nIrreps() << " " << _ncore.size() << std::endl;
  assert( n_occ.size() == p_pgs->nIrreps());
  assert( _ncore.size() == n_occ.size() );
  for ( Irrep ir = 0; ir < n_occ.size(); ++ir ) {
    for ( uint j = uint(_ncore[ir]); j < uint(_ncore[ir]+n_occ[ir]); ++j ) {
      uint iao = 0;
      double maxcoef = 0.0;
      for ( uint i = 0; i < p_pgs->norbs(ir); ++i ) {
        double val = std::abs(_orbs.get(i,j));
        if ( val > maxcoef ) {
          iao = i;
          maxcoef = val;
        }
      }
      occ[ir].push_back(iao-_ncore[ir]);
      if ( maxcoef < throcc ) {
        warning("Small max occupation for the basis function "+num2str(iao+1,std::dec)+": "+num2str(maxcoef,std::dec));
      }
    }
    if ( sortocc && nclos.size() >= occ.size() ){
      // sort closed-shell orbitals
      std::sort(occ[ir].begin(),occ[ir].begin()+nclos[ir]);
      // sort open-shell orbitals
      std::sort(occ[ir].begin()+nclos[ir],occ[ir].end());
    }
  }
  return occ;
}

void Odump::transform(const Integ2ab& trmat, bool frozcore)
{
  assert( trmat.pgs()->nIrreps() == p_pgs->nIrreps() );
  for ( Irrep ir = 0; ir < p_pgs->nIrreps(); ++ir ) {
    uint firstao = p_pgs->_firstorb4irrep[ir],
         lenao = p_pgs->_norb4irrep[ir];
    uint firstmo_orbs = p_pgs->_firstorb4irrep[ir],
         lenmo_trmat = trmat.pgs()->_norb4irrep[ir];
    uint icor = 0;
    if (frozcore) icor = _ncore[ir];
    std::vector<double> coefs(lenmo_trmat);
    for ( uint iao = 0; iao < lenao; ++iao ) {
      for ( uint imo = 0; imo < lenmo_trmat; ++imo ) {
        coefs[imo] = 0.0;
        for ( uint kmo = 0; kmo < lenmo_trmat; ++kmo ) {
          coefs[imo] += _orbs.get(iao,kmo+icor,ir) * trmat.get(kmo,imo,ir);
        }
      }
      for ( uint imo = 0; imo < lenmo_trmat; ++imo ) {
        _orbs.set(firstao+iao,imo+icor+firstmo_orbs,coefs[imo]);
      }
    }
  }
}
