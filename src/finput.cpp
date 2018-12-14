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
  std::string  dmfile = Input::sPars["dm"]["in"];
  // TODO: move somewhere else!
  // TODO: the input file should be in the \input command!
  lui
    ibeg = IL::skip(_input,0," ,"),
    iend = IL::skipr(_input,_input.size()," ,");
  std::string inputfile = _input.substr(ibeg,iend-ibeg);
  if ( inputfile == "" ) {
    error("Empty hamiltonian specification!");
  }
  xout << "inputfile: *" << inputfile << "*" << std::endl;
  std::string outputfile = Input::sPars["ham"]["out"];
  if ( outputfile == "" ) {
    outputfile = FileName(inputfile,true)+"_NEW.FCIDUMP";
  }
  if ( Input::iPars["ham"]["store"] <= 0 ) {
    outputfile.clear();
  }
  if ( Input::iPars["output"]["fcinamtoupper"] > 0 )
    outputfile = uppercase(outputfile);
  
  Hdump dump(inputfile);
  if ( outputfile != "" )
    dump.store(outputfile);
  if ( dmfile != "" ) {
    DMdump dmdump(dmfile,dump.norb(),dump.nelec());
//     Occupation hfocc(dump.pgs(),dump.nclos(),dump.nocc());
//     DMdump hfdm(dump.norb(),hfocc);
    Fock_matrices fock(dump, dmdump);
//     Fock_matrices fock(dump, hfdm); 
    fock.diagonalize(dmdump);
    
  }
//   Input::iPars["ham"]["nosym"] = 0;
//   dump.store(outputfile+"sym");
  handle_orbdump(dump);
  return true;
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
      _foreach_cauto ( std::vector<int>, iso, socc ) {
        xout << *iso << "  ";
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
      _foreach_cauto(TParArray,iocc,occ)
        xout << *iocc << ",";
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

