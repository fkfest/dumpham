#include "finput.h"

Finput::Finput(bool ham) :
_ham(ham){}

Finput::Finput(std::string paramspath, const std::vector<std::string>& cmd_inps) :
_ham(false)
{
  InitInpars(paramspath);
  // set input params from cmd_inps
  for (auto ps: cmd_inps) {
    IL::changePars(ps, 0);
  }
}

void Finput::InitInpars(std::string paramspath)
{
  std::string finp_file(paramspath+"params.reg");
  std::ifstream finp;
  finp.open(finp_file.c_str());
  if (finp.is_open()) {
    std::string
      line, type, set, name;
    while (finp.good()) {
      std::getline (finp,line);
      if ( !line.empty()){
        _xout3(line << std::endl);
        set = IL::key(line,"set");
        if ( set.size() == 0 ) continue;
        type = IL::key(line,"type");
        name = IL::key(line,"name");
        if ( type.size() == 0 )
          error("no type is given");
        else if ( type[0] == 's' ){
          Input::sPars[set][name] = IL::key(line,"value");
        } else if ( type[0] == 'i' ){
          int x;
          if (!str2num<int>(x,IL::key(line,"value"),std::dec))
            error("integer parameter is not integer :"+line);
          Input::iPars[set][name] = x;
        } else if ( type[0] == 'f' ){
          double x;
          if (!str2num<double>(x,IL::key(line,"value"),std::dec))
            error("float parameter is not float :"+line);
          Input::fPars[set][name] = x;
        } else if ( type[0] == 'a' ){
          Input::aPars[set][name] = IL::parray(IL::key(line,"value"));
        } else
          error("unknown type in params.reg");
//       finput+=line;
      }
    }
  }
  else
    error(finp_file+": Bad input-parameters file!");
  finp.close();
}

bool Finput::addline(const std::string& line)
{
  const TParArray& bham = Input::aPars["syntax"]["bham"];
  const TParArray& eham = Input::aPars["syntax"]["eham"];
  const TParArray& newcs = Input::aPars["syntax"]["newcommand"];
  const TParArray& newops = Input::aPars["syntax"]["newoperator"];
  const TParArray& comments = Input::aPars["syntax"]["comment"];
  short iprint = Input::iPars["output"]["level"];
  if ( iprint > 0 && _ham )
    _inham.push_back(line);
  if ( iprint > 1 && !_ham )
    _inlines.push_back(line);
  bool newham = false;
  lui ipos=0, ipend;
  // and skip " " on begin
  ipos = IL::skip(line,ipos," ");
  // line without front-spaces and comments
  std::string linesp;
  // remove comments
  for (unsigned long int i=ipos; i<line.size(); i++) {
    if(InSet(line.substr(i,1),comments)) // comment
      break;
    linesp += line[i];
  }
  // remove spaces at the end
  ipend = IL::skipr(linesp,linesp.size()," ");
  linesp = linesp.substr(0,ipend);
  ipos = 0;
  ipend = IL::nextwordpos(linesp,ipos,false);
  if (ipos == ipend){ // empty line
  } else if (InSet(linesp, bham)) {// begin hamiltonian specification
    _input="";
    if ( iprint > 1 && !_ham )
      _inlines.pop_back();
    _ham=true;
    //analyzenewops();
  } else if (InSet(linesp, eham)) {// end hamiltonian specification
    if ( iprint > 0 && _ham )
      _inham.pop_back();
    _ham=false;
    analyzeline();
    _input="";
    process_dump(*_dump);
    newham = true;
  } else if (InSet(linesp.substr(ipos,ipend-ipos), newcs)) {// newcommand
    ipos = IL::addnewcom(linesp,ipend);
  } else if (InSet(linesp.substr(ipos,ipend-ipos), newops)) {// newoperator
    ipos = IL::addnewcom(linesp,ipend,"newoperator");
  } else if (!_ham){
    IL::changePars(linesp,ipos);
  } else {
    _input += linesp;
    _input += " "; // add space for separation
  }
  return newham;
}
std::string Finput::input() const
{ return _input; }

bool Finput::analyzeline()
{
//   xout << "analyzeline " << _input << std::endl;
  char ch;
  lui i = 0;
  while (i < _input.size()) {
    ch=_input[i];
    if ( ch == '\\' ) { // command
      ++i;
      i = analyzecommand(i) - 1;
    } else if ( ch == '+' ) { // add next
      _add = true;
    } else if ( ch == '-' ) {
      _add = true;
      _scale = -1.0;
    } else if ( isdigit(ch) ) { // number
      lui ipos = IL::nextwordpos(_input,i);
      double scal;
      if ( str2num<double>(scal,_input.substr(i,ipos-i),std::dec) ) {
        _scale *= scal;
      } else {
        error("Could not transform to a number: "+_input.substr(i,ipos-i));
      }
    } else if ( ch == ' ' ) { // skip
    } else {
      lui ipos1 = IL::nextwordpos(_input,i);
      lui iend = IL::skipr(_input,ipos1," ,");
      analyzeham(_input.substr(i,iend-i));
      i = ipos1;
    }
    ++i;
  }
  return true;
}

lui Finput::analyzecommand(lui ipos)
{
  TsPar& commands = Input::sPars["command"];
  lui ipos1 = IL::nextwordpos(_input,ipos);
  std::string str = _input.substr(ipos,ipos1-ipos);
  ipos = ipos1;
  if ( str == commands["operator"] ) { // operators
  } else if ( str == commands["fcidump"] ) {
    ipos1 = IL::nextwordpos(_input,ipos);
    ipos = IL::skip(_input,ipos,"{} ");
    lui ipos2 = IL::skipr(_input,ipos1,"{} ");
    std::string filename = _input.substr(ipos,ipos2-ipos);
    analyzeham(filename);
  } else if ( str == commands["outfcidump"] ) {
    ipos1 = IL::nextwordpos(_input,ipos);
    ipos = IL::skip(_input,ipos,"{} ");
    lui ipos2 = IL::skipr(_input,ipos1,"{} ");
    Input::sPars["ham"]["out"] = _input.substr(ipos,ipos2-ipos);
  } else if ( str == commands["geom"] ) {
    ipos1 = IL::nextwordpos(_input,ipos);
    ipos = IL::skip(_input,ipos,"{} :");
    lui ipos2 = IL::skipr(_input,ipos1,"{} ");
    std::string geomdef = "geom,"+_input.substr(ipos,ipos2-ipos);
    IL::changePars(geomdef, 0);
  } else if ( str == commands["hubbard"] ) {
    ipos1 = IL::nextwordpos(_input,ipos);
    ipos = IL::skip(_input,ipos,"{} :");
    lui ipos2 = IL::skipr(_input,ipos1,"{} ");
    std::string hubdef = "hubbard,"+_input.substr(ipos,ipos2-ipos);
    IL::changePars(hubdef, 0);
    analyzehubham();
  } else if ( str == commands["heisenberg"] ) {
    ipos1 = IL::nextwordpos(_input,ipos);
    ipos = IL::skip(_input,ipos,"{} :");
    lui ipos2 = IL::skipr(_input,ipos1,"{} ");
    std::string heidef = "heisenberg,"+_input.substr(ipos,ipos2-ipos);
    IL::changePars(heidef, 0);
    analyzeheisham();
  } else if ( str == commands["ppp"] ) {
    ipos1 = IL::nextwordpos(_input,ipos);
    ipos = IL::skip(_input,ipos,"{} :");
    lui ipos2 = IL::skipr(_input,ipos1,"{} ");
    std::string pppdef = "ppp,"+_input.substr(ipos,ipos2-ipos);
    IL::changePars(pppdef, 0);
    analyzepppham();
  }
  return ipos1;
}

/*
uint string2orbs(std::vector<uint>& orbs, std::vector<uint>& norbs, std::string orbstr, 
                 const std::vector<uint>& ntotorbs4sym)
{
  IL::delbrack(orbstr);
  bool opposite = false;
  std::size_t indhat = orbstr.find('^');
  if (indhat != std::string::npos ) {
    // reverse meaning of the list
    opposite = true;
    if ( indhat != 0 ) error("USE ^ AS THE FIRST CHARACTER in OCCA/OCCB"); 
    orbstr=orbstr.substr(1);
  }
  if ( orbstr.find_first_not_of("1234567890.+- ") != std::string::npos ) {
    error("USE ONLY ^1234567890.+- in OCCA/OCCB");
  }
  uint norbtot = 0;
  std::vector<uint> orboff;
  for ( auto no4s: ntotorbs4sym ) {
    orboff.push_back(norbtot);
    norbtot += no4s;
  }
  std::vector<uint> orbsmask(norbtot,0);
  
  
}
*/

bool Finput::analyzeham(const std::string& inputfile)
{
//   xout << "analyzeham " << _input << std::endl;
  int
    ists = Input::iPars["ham"]["simtrasym"];
  if ( inputfile == "" ) {
    error("Empty hamiltonian specification!");
  }
  xout << "inputfile: *" << inputfile << "*" << std::endl;

  if ( !_dump ) {
    _dump = std::unique_ptr<Hdump>(new Hdump(inputfile));
    _dump->read_dump();
    if (scale()) _dump->scale(_scale);
    if ( ists != 0 ) { // check similarity transformation flag
      bool simtra = (ists > 0);
      if ( simtra != _dump->simtra() ) {
        // the output should have different sim.tra. symmetry
        auto newdump = std::unique_ptr<Hdump>(new Hdump(*_dump,0,ists));
        newdump->import(*_dump);
        _dump = std::move(newdump);
      }
    }
  } else {
    if ( !_add ) error ("The Hamiltonian " + _dump->fcidump_filename() + " will be overwritten by "+inputfile);
    Hdump dump(inputfile);
    // check whether we will have to enlarge _dump
    if ( ists ==0 && !_dump->simtra() && dump.simtra() ) {
      // transform old dump into similarity transformed version first!
      auto newdump = std::unique_ptr<Hdump>(new Hdump(*_dump,dump));
      newdump->import(*_dump);
      _dump = std::move(newdump);
    }
    dump.read_dump();
    if (scale()) dump.scale(_scale);
    // add to the old dump
    _dump->add(dump);
  }
  
  const std::string& occa = Input::sPars["orbs"]["occa"];
  const std::string& occb = Input::sPars["orbs"]["occb"];
  if ( !occa.empty() || !occb.empty() ) {
    // set refdet using occa and occb
    error("Not implemented yet!"); 
  }
  
  double addS2 = Input::fPars["pert"]["S2"];
  if ( addS2 < -Numbers::small || addS2 > Numbers::small ) {
    // add S^2 
    _dump->addS2(addS2);
  }
  _add = false;
  _scale = 1.0;
  return true;
}

Periodic Finput::analyzegeom()
{
  const TParArray& dim = Input::aPars["geom"]["dimension"];
  TParArray pbc = Input::aPars["geom"]["pbc"];
  const TParArray& ucell_str = Input::aPars["geom"]["ucell"];
  const TParArray& lat_str = Input::aPars["geom"]["lat"];

  // dimensions
  std::vector<uint> dims;
  apars2nums<uint>(dims,dim,std::dec,"Dimensions in Hubbard are not integer");;
  // PBC
  if ( pbc.size() == 0 ) pbc.push_back("0");
  if ( pbc.size() != dim.size() && pbc.size() != 1 )
    error("Define PBC-type for each dimension!");
  pbc.resize(dim.size(), pbc.front());
  std::vector<int> pbcs;
  apars2nums<int>(pbcs,pbc,std::dec,"PBC is not integer");

  // unit cell
  UCell ucell;
  for (auto site: ucell_str) {
    TParArray coord_str = IL::parray(site);
    Coords coords;
    apars2nums<double>(coords,coord_str,std::dec,"coord in ucell not float");
    ucell.add(coords);
  }
  if ( ucell.nsites() == 0 ) ucell.set_default(dims.size());
//  xout << "Unit cell: " << ucell << std::endl;
  // lattice vectors
  Lattice lat;
  for (auto lv: lat_str) {
    TParArray coord_str = IL::parray(lv);
    Coords coords;
    apars2nums<double>(coords,coord_str,std::dec,"coord in lat not float");
    lat.add(coords);
  }
  if ( lat.ndim() == 0 ) lat.set_default(dims.size());
//  xout << "Lattice: " << lat << std::endl;

  return Periodic(ucell,lat,dims,pbcs);
}

bool Finput::analyzehubham()
{
  int charge = Input::iPars["hubbard"]["charge"];
  int ms2 = Input::iPars["hubbard"]["ms2"];
  double Upar = Input::fPars["hubbard"]["U"];
  const TParArray& tparsarray = Input::aPars["hubbard"]["t"];

  // hopping
  std::vector<double> tpars;
  apars2nums<double>(tpars,tparsarray,std::dec,"t parameter is not float");

  Periodic persym = analyzegeom();
  _dump = std::unique_ptr<Hdump>(new Hdump(persym,charge,ms2,Upar,tpars));
  return true;
}

bool Finput::analyzeheisham()
{
  int ms2 = Input::iPars["heisenberg"]["ms2"];
  int norbs = Input::iPars["heisenberg"]["norbs"];
  const TParArray& jparsarray = Input::aPars["heisenberg"]["j"];
  const TParArray& kparsarray = Input::aPars["heisenberg"]["k"];


  // J
  std::vector<double> jpars, kpars;
  apars2nums<double>(jpars,jparsarray,std::dec,"j parameter is not float");
  apars2nums<double>(kpars,kparsarray,std::dec,"k parameter is not float");


  Periodic persym = analyzegeom();
  _dump = std::unique_ptr<Hdump>(new Hdump(persym,ms2,norbs,jpars,kpars));
  return true;
}

bool Finput::analyzepppham()
{
  int charge = Input::iPars["ppp"]["charge"];
  int ms2 = Input::iPars["ppp"]["ms2"];
  double Upar = Input::fPars["ppp"]["U"];
  double apar = Input::fPars["ppp"]["a"];
  const TParArray& tparsarray = Input::aPars["ppp"]["t"];
  double tdecay = Input::fPars["ppp"]["tdecay"];
  double udecay = Input::fPars["ppp"]["udecay"];
  double frozen = Input::fPars["ppp"]["frozen"];

  // hopping
  std::vector<double> tpars;
  apars2nums<double>(tpars,tparsarray,std::dec,"t parameter is not float");

  Periodic persym = analyzegeom();
  _dump = std::unique_ptr<Hdump>(new Hdump(persym,charge,ms2,Upar,apar,tpars,tdecay,udecay,frozen));
  return true;
}

void Finput::process_dump(Hdump& dump)
{
  std::string outputfile = Input::sPars["ham"]["out"];
  if ( outputfile == "" ) {
    if ( FileName(dump.fcidump_filename(),true) == "" ) {
      outputfile = "FCIDUMP";
    } else {
      outputfile = FileName(dump.fcidump_filename(),true)+"_NEW.FCIDUMP";
    }
  }
  if ( Input::iPars["ham"]["store"] <= 0 ) {
    outputfile.clear();
  }
  if ( Input::iPars["output"]["fcinamtoupper"] > 0 )
    outputfile = uppercase(outputfile);

  if ( outputfile != "" )
    dump.store(outputfile);

  std::string  dmfile = Input::sPars["dm"]["in"];
  std::string  fmfile = Input::sPars["fm"]["outmat"];
  DMdump dmdump;
  if ( dmfile != "" ) {
    // read density matrices
    dmdump = DMdump(dmfile,dump.norb(),dump.nelec());
  } else if ( fmfile != "" ) {
    // Build one-determinant density matrices
    Occupation hfocc(dump.pgs(),dump.nclos(),dump.nocc());
    dmdump = DMdump(dump.norb(),hfocc);
  }
  if ( fmfile != "" ) {
    // Fock matrix
    xout << "Calculate Fock matrix" << std::endl;
    Fock_matrices fock(dump, dmdump);
    xout << "Fock calculated " << std::endl;
    fock.store(fmfile);
    xout << "Fock stored " << std::endl;
    fock.diagonalize(dmdump);
    xout << "Fock diagonlized " << std::endl;
  }
//   Input::iPars["ham"]["nosym"] = 0;
//   dump.store(outputfile+"sym");
  handle_orbdump(dump);
}

void Finput::handle_orbdump(const Hdump& dump)
{
  std::string orboutputfile = Input::sPars["orbs"]["out"];
  bool orbdump = ( orboutputfile != "" );
  if ( orbdump && Input::iPars["output"]["orbnamtolower"] > 0 )
    orboutputfile = lowercase(orboutputfile);
  const std::string& orbcoefs_input = Input::sPars["orbs"]["in"];
  if ( orbcoefs_input != "" ) {
    if ( dump.pgs_wcore().nIrreps() == 0 ) {
      error("Please provide number of core orbitals in each symmetry to read orbital coeffs","Finput::handle_orbdump");
    }
    Odump odump(dump.pgs_wcore(),dump.ncore(),orbcoefs_input);
    if ( Input::iPars["orbs"]["guessoccu"] > 0 ) {
      Occupation occguess = odump.guess_occupation(dump.nclos(),dump.nocc());
      xout << "Guessed occupation: " << occguess << std::endl;
      xout << "Guessed spin occupation: ";
      std::vector<int> socc = occguess.spinocc(1);
      for ( const auto& iso: socc ) {
        xout << iso << "  ";
      }
      xout << std::endl;
    }
    if ( orbdump ) {
      odump.store(orboutputfile);
    }
  } else if ( orbdump ) {
    const TParArray& occ = Input::aPars["orbs"]["occvec"];
    if ( occ.size() > 0 ) {
      xout << "Occupation: ";
      for ( const auto& iocc: occ )
        xout << iocc << ",";
      xout << std::endl;
    }
    std::vector<int> occ_spin;
    apars2nums<int>(occ_spin,occ,std::dec);
    Occupation occupation(dump.pgs(), dump.orborder(), occ_spin);
    Odump odump(dump.pgs(), occupation);
    odump.store(orboutputfile);
  }
}

std::ostream& operator<<(std::ostream& o, const Finput& inp)
{
  for (unsigned long int i=0; i<inp.inlines().size(); i++)
    o << inp.inlines().at(i);
  return o;
}

