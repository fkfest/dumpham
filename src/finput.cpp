#include "finput.h"

Finput::Finput(bool ham) : 
_ham(ham){}

Finput::Finput(std::string paramspath) :
_ham(false)
{
  InitInpars(paramspath);
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
  }
  return ipos1;
}

bool Finput::analyzeham(const std::string& inputfile)
{
//   xout << "analyzeham " << _input << std::endl;
  if ( inputfile == "" ) {
    error("Empty hamiltonian specification!");
  }
  xout << "inputfile: *" << inputfile << "*" << std::endl;
  
  if ( !_dump ) { 
    _dump = std::unique_ptr<Hdump>(new Hdump(inputfile));
    _dump->read_dump();
    if (scale()) _dump->scale(_scale);
  } else {
    if ( !_add ) error ("The Hamiltonian " + _dump->fcidump_filename() + " will be overwritten by "+inputfile);
    Hdump dump(inputfile);
    // check whether we will have to enlarge _dump
    if ( !_dump->simtra() && dump.simtra() ) {
      // transform old dump into similarity transformed version first!
      Hdump newdump(*_dump,dump);
      newdump.import(*_dump);
      error("transform "+_dump->fcidump_filename() + " to ST!");
    }
    dump.read_dump();
    if (scale()) dump.scale(_scale);
    // add to the old dump
    
  }
  _add = false;
  _scale = 1.0;
  
  return true;
}

void Finput::process_dump(Hdump& dump)
{
  std::string outputfile = Input::sPars["ham"]["out"];
  if ( outputfile == "" ) {
    outputfile = FileName(dump.fcidump_filename(),true)+"_NEW.FCIDUMP";
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
    Occupation occupation(dump.pgs(), occ_spin);
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

