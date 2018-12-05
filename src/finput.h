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

/*!
    Input analyzer 
*/
class Finput {
public:
  // constructor
  Finput ( bool ham = false );
  // constructor + init input-parameters
  Finput( std::string paramspath );
  // add string
  bool addline( const std::string& line );
  // get input
  std::string input() const;
  // analyze input line
  bool analyzeline();
  // clear all arrays
  void clear() {_inlines.clear(); _inham.clear(); _input.clear(); _ham = false;};
  // return input lines
  const std::vector<std::string> & inlines() const { return _inlines;};
  const std::vector<std::string> & inham() const { return _inham;};
  
private:
  // initialyse default input-parameters 
  void InitInpars(std::string paramspath);
  // variables
  std::string _input;
  bool _ham;
  std::vector<std::string> _inlines;
  std::vector<std::string> _inham;
};

std::ostream & operator << (std::ostream & o, Finput const & inp);

#endif

