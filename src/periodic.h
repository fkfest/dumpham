#ifndef Periodic_H
#define Periodic_H

#include <vector>
#ifdef MOLPRO
#include "hdtypes.h"
#include "hdcommon.h"
#else
#include "globals.h"
#include "utilities.h"
#endif
namespace HamDump {
// typedef std::vector<double> Coords;
struct Coords : public std::vector<double> {
  Coords(uint ndim = 0, double val = 0.0) : std::vector<double>(ndim,val) {};
  // dot product
  double dot(const Coords& other) const { assert(other.size() == size());
                                          double res = 0.0;
                                          for (uint i = 0; i < size(); ++i)
                                            res += (*this)[i]*other[i];
                                          return res; }
  // square of vector norm
  double norm2() const { double res = 0.0; for(auto a: *this) res += a*a; return res; }
  // square of distance
  double dist2(const Coords& other) const { assert(other.size() == size());
                                            double res = 0.0;
                                            for (uint i = 0; i < size(); ++i) {
                                              double diff = (*this)[i] - other[i];
                                              res += diff*diff;
                                            }
                                            return res; }
  // "orthogonal" vectors
  // 2D: orthogonal vector
  Coords orthogvec() const { assert(size() == 2);
                          Coords ort(2); ort[0] = (*this)[1]; ort[1] = -(*this)[0];
                          return ort; }
  // 3D: cross product
  Coords cross(const Coords& vec) const { assert(size() == 3 && vec.size() == 3); Coords ort(3);
                                          ort[0] = (*this)[1]*vec[2]-(*this)[2]*vec[1];
                                          ort[1] = (*this)[2]*vec[0]-(*this)[0]*vec[2];
                                          ort[2] = (*this)[0]*vec[1]-(*this)[1]*vec[0];
                                          return ort; }
  void flip_sign() { for(auto & a: *this) a = -a; }
  // align to the vector (i.e., make the dot product positive)
  void align(const Coords& vec) { assert( size() == vec.size() );
                                  if ( dot(vec) < 0.0 ) flip_sign(); }
  void print_xyz(std::ostream& o) const;
};
std::ostream & operator << (std::ostream& o, const Coords& co);

typedef std::vector<Coords> CoordSet;

/*
 * Periodic site
 */
struct PSite : public std::vector<int> {
  // site in unit cell
  uint site = 0;
  PSite(uint ndim = 0) : std::vector<int>(ndim,0) {};
  void zero(bool full = true) {
    for (auto & i: *this) i = 0;
    if (full) site = 0;
  }
};
std::ostream & operator << (std::ostream& o, const PSite& ps);
/*
 * Unit cell
 */
struct UCell : public CoordSet {
  UCell() : CoordSet() {};
  void add(const Coords& site) {assert(size()==0||back().size()==site.size());
                                push_back(site);};
  void set_default(uint ndims);
  uint nsites() const { return size(); };
  const Coords& coords(uint site) const { assert(site<size()); return (*this)[site]; };
  // go to next site
  bool next(PSite& ps, uint neighbours) const {
    ++ps.site;
    for ( ; ps.site < size(); ++ps.site ) {
      if ( _redundant[ps.site] == 0 || (_redundant[ps.site] & neighbours) == 0 )
        // non-redundant site
        return true;
    }
    ps.site = 0;
    return false;
  }
  void set_redundancy(const std::vector<uint>& reds) {
    assert( reds.size() == size() );
    _redundant = reds;
  }
  std::string sitename(uint site) const { (void) site; return "H"; }
  // site on boundary and redundant by applying the listed translations to other sites
  std::vector<uint> _redundant;
  // names of sites?
  // std::vector<std::string> _names;
};
std::ostream & operator << (std::ostream& o, const UCell& uc);

/*
 * Lattice
 */
struct Lattice : public CoordSet {
  Lattice() : CoordSet() {};
  void add(const Coords& vect) {assert(size()==0||back().size()==vect.size());
                                push_back(vect);};
  void set_default(uint ndims);
  uint ndim() const { return size(); };
  // lattice vector
  const Coords& vec(uint ivec) const { assert(ivec<size()); return (*this)[ivec]; };
  // vector orthogonal to all other vectors
  Coords orthogvec2others(uint lv) const;
  // names of sites?
  // std::vector<std::string> _names;
};
std::ostream & operator << (std::ostream& o, const Lattice& lat);

/*
 * Periodic system
 */
struct Periodic {
  Periodic() {};
  Periodic(const UCell& ucell, const Lattice& lat,
          const std::vector<uint>& ncells, const std::vector<int>& pbcs);
  uint ndim() const { return _ncells.size(); }
  uint nsites_in_supercell() const { return _nsites_in_supercell; }
  Coords site_coords(const PSite& isite) const;
  const std::vector<uint>& ncells() const { return _ncells; }
  const std::vector<int>& pbcs() const { return _pbcs; }
  bool periodic() const { for (auto p: _pbcs) {if ( p != 0 ) return true;} return false; }
  bool antiperiodic() const { for (auto p: _pbcs) {if ( p < 0 ) return true;} return false; }
  void check_boundaries();
  // minimal square of distance between sites (considering PBC)
  double dist2(const PSite& ps1, const PSite& ps2) const;
  uint neighbours(const PSite& ps) const {
    assert(ps.size() == _ncells.size());
    assert(ps.size() == _pbcs.size());
    uint neigh = 0;
    for ( uint id = 0; id < ps.size(); ++id ) {
      if (_pbcs[id] != 0 || ps[id]+1 < int(_ncells[id]))
        neigh |= 1 << id;
    }
    return neigh;
  }
  // go to next ps
  bool next(PSite& ps) const {
    if ( _ucell.next(ps,neighbours(ps)) )
      return true;
    for ( uint id = 0; id < ps.size(); ++id ) {
      ++ps[id];
      if (ps[id] < int(_ncells[id]) )
        return true;
      else {
        ps[id] = 0;
      }
    }
    return false;
  };
  // ndist shortest square distances
  std::vector<double> dist2neighbours(uint ndist) const;

  // save supercell geometry in a file
  void xyz_supercell(const std::string& xyzfilename) const;

private:
  // coordinates of atoms in the unit cell
  UCell _ucell;
  // lattice vectors
  Lattice _lat;
  // number of cells in each dimension
  std::vector<uint> _ncells;
  // periodic boundary condition in each dimension
  std::vector<int> _pbcs;
  // number of sites in supercell
  uint _nsites_in_supercell = 0;
};

} //namespace HamDump

#endif
