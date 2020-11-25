#ifndef Finput_H
#define Finput_H
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <typeinfo>
#include "utilities.h"
#include "globals.h"
#include "inpline.h"
#include "hdump.h"
#include "odump.h"
#include "density_mat.h"
#include "fock_mat.h"
#include "periodic.h"

/*!
    Input analyzer 
*/
class Finput {
public:
  // constructor
  Finput ( bool ham = false );
  // constructor + init input-parameters
  Finput( std::string paramspath, const std::vector<std::string>& cmd_inps );
  // add string
  bool addline( const std::string& line );
  // get input
  std::string input() const;
  // analyze input line
  bool analyzeline();
  // analyze hamiltonian
  bool analyzeham(const std::string& inputfile);
  // analyze Hubbard Hamiltonian
  bool analyzehabham();
  // clear all arrays
  void clear() {_inlines.clear(); _inham.clear(); _input.clear(); _ham = false;};
  // return input lines
  const std::vector<std::string> & inlines() const { return _inlines;};
  const std::vector<std::string> & inham() const { return _inham;};
  
private:
  // initialyse default input-parameters 
  void InitInpars(std::string paramspath);
  // analyze command from the input line after backslash at ipos-1
  lui analyzecommand(lui ipos);
  // process the dump
  void process_dump(Hdump& dump);
  // handle orbdump
  void handle_orbdump(const Hdump& dump);
  bool scale() const { return std::abs(_scale - 1.0) > Numbers::small; }
  // variables
  std::string _input;
  bool _ham;
  std::vector<std::string> _inlines;
  std::vector<std::string> _inham;
  bool _add = false;
  double _scale = 1.0;
  std::unique_ptr<Hdump> _dump;
};

std::ostream & operator << (std::ostream & o, Finput const & inp);

#endif

