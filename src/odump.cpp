#include "odump.h"

void Occupation4Irrep::spinocc(std::vector<int>& socc) const
{
  for ( uint imo = 0; imo < this->size(); ++imo ){
    socc.push_back((*this)[imo]*2+1);
    if ( imo < _nclos )
      socc.push_back((*this)[imo]*2+2);
  }
}

Occupation::Occupation(const PGSym& pgs, const FDPar& nclos, const FDPar& nocc)
    : p_pgs(&pgs)
{
  assert(nocc.size() >= p_pgs->nIrreps() && nclos.size() >= p_pgs->nIrreps());
  for ( uint ir = 0; ir < p_pgs->nIrreps(); ++ir ) {
    Occupation4Irrep o4ir(nocc[ir]);
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
  _foreach_cauto( std::vector< int >, iocc, occs ){
    int iorb1 = int((*iocc-ibase)/2);
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
    _foreach_cauto( std::vector< int >, iocc, singly[ir] ){
      (*this)[ir].push_back(*iocc);
    }
  }
}

std::vector<int> Occupation::spinocc() const
{
  std::vector<int> socc;
  _foreach_cauto(Occupation, oc, *this)
    oc->spinocc(socc);
  std::sort(socc.begin(),socc.end());
  return socc;
}


std::ostream& operator<<(std::ostream& o, const Occupation& occ)
{
  _foreach_cauto(Occupation,iocc,occ){
    _foreach_cauto(Occupation4Irrep, ioir, *iocc)
      o << *ioir << ",";
  }
  return o;
}


Odump::Odump(const PGSym& pgs, Occupation occs ) : p_pgs(&pgs)
{
  _orbs = Integ2ab(pgs);
  if ( occs.empty() ) occs.resize(p_pgs->nIrreps());
  uint iao;
  xout << "Spatial orbital occupation: " << occs << std::endl;
  assert( occs.size() == p_pgs->nIrreps() );
  for ( uint ir = 0; ir < p_pgs->nIrreps(); ++ir ) {
    std::vector<bool> occupied(p_pgs->norbs(ir),false);
    for ( uint imo = 0; imo < occs[ir].size(); ++imo ) {
      iao = occs[ir][imo];
      _orbs.set(iao,imo,1.0);
      occupied[iao] = true;
    }
    uint imo = occs[ir].size(); 
    for ( uint iao = 0; iao < p_pgs->norbs(ir); ++iao ) {
      if (!occupied[iao]) {
        _orbs.set(iao,imo,1.0);
        ++imo;
      }
    }
    assert( imo == p_pgs->norbs(ir) );
  }
}
Odump::Odump(const PGSym& pgs, std::string orbdump) : p_pgs(&pgs)
{
  _orbs = Integ2ab(pgs);
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
    _foreach_cauto(TParArray,ic,coefs_line){
      double val;
      if ( str2num<double>(val,*ic,std::dec) ) {
        _orbs.set(idx,val);
        ++idx;
      } else
        error("Not a number: "+*ic);
    }
  }
  if ( idx != _orbs.nelem() )
    error("Number of orbital coefficients not consistens with the number of orbitals");
  // transpose to get the AO index fastest
  for ( Irrep ir = 0; ir < p_pgs->nIrreps(); ++ir ) {
    for ( uint j = p_pgs->_firstorb4irrep[ir]; j < uint(p_pgs->_firstorb4irrep[ir] + p_pgs->_norb4irrep[ir]); ++j ) {
      for ( uint i = p_pgs->_firstorb4irrep[ir]; i < j; ++i ) {
        double tmp = _orbs.get(i,j);
        _orbs.set(i,j,_orbs.get(j,i));
        _orbs.set(j,i,tmp);
      }
    }
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
    for ( uint i = 0; i < p_pgs->ntotorbs(); ++i ) {
      for ( uint j = 0; j < p_pgs->ntotorbs(); ++j ) {
        if ( j > 0 && maxlen > 0 && j % maxlen == 0 ) outputStream << std::endl;
        double val = _orbs.get_with_pgs(i,j);
        if ( !scientific && val <= 0.0 && val >= -1.e-20 ) val = 1.e-20;
        if ( !scientific && val > -1.e-20 ) outputStream << " ";
        outputStream << val << ", ";
      } 
      outputStream << std::endl;
    }
  } else {
    for ( Irrep ir = 0; ir < p_pgs->nIrreps(); ++ir ) {
      for ( uint i = p_pgs->_firstorb4irrep[ir]; i < uint(p_pgs->_firstorb4irrep[ir] + p_pgs->_norb4irrep[ir]); ++i ) {
        for ( uint j = p_pgs->_firstorb4irrep[ir]; j < uint(p_pgs->_firstorb4irrep[ir] + p_pgs->_norb4irrep[ir]); ++j ) {
          if ( j > 0 && maxlen > 0 && j % maxlen == 0 ) outputStream << std::endl;
          double val = _orbs.get(i,j);
          if ( !scientific && val <= 0.0 && val >= -1.e-20 ) val = 1.e-20;
          if ( !scientific && val > -1.e-20 ) outputStream << " ";
          outputStream << val << ", ";
        } 
        outputStream << std::endl;
      }
    }
  }
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
  assert( n_occ.size() == p_pgs->nIrreps());
  for ( Irrep ir = 0; ir < n_occ.size(); ++ir ) {
    for ( uint j = 0; j < uint(n_occ[ir]); ++j ) {
      uint iao = 0;
      double maxcoef = 0.0;
      for ( uint i = 0; i < p_pgs->norbs(ir); ++i ) {
        double val = std::abs(_orbs.get(i,j));
        if ( val > maxcoef ) {
          iao = i;
          maxcoef = val;
        }
      }
      occ[ir].push_back(iao);
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
