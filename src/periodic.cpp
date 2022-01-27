#include "periodic.h"

#include <fstream>

namespace HamDump {

void UCell::set_default(uint ndims)
{
  assert(size() == 0);
  add(Coords(ndims,0.0));
}

void Lattice::set_default(uint ndims)
{
  assert(size() == 0);
  Coords c(ndims,0.0);
  for ( uint idim = 0; idim < ndims; ++idim ) {
    c[idim] = 1.0;
    add(c);
    c[idim] = 0.0;
  }
}

Coords Lattice::orthogvec2others( uint lv ) const
{
  assert( size() > 0 && lv < size() );
  Coords ret;
  if ( size() == 1 ) {
    ret = Coords(1,1.0);
  } else if ( size() == 2 ) {
    ret = vec(1-lv).orthogvec();
  } else if ( size() == 3 ) {
    uint i = 1, j = 2;
    if ( lv == 1 ) {
      i = 0;
    } else if ( lv == 2 ){
      j = 0;
    }
    ret = vec(i).cross(vec(j));
  } else {
    error("ndim > 3 in Lattice::orthogvec");
  }
  ret.align(vec(lv));
  return ret;
}


Periodic::Periodic(const UCell& ucell, const Lattice& lat,
                 const std::vector<uint>& ncells, const std::vector<int>& pbcs)
{
  assert( ncells.size() == pbcs.size() );
  assert( lat.size() > 0 && lat.size() == lat[0].size() );
  assert( ucell.size() > 0 && ucell[0].size() == lat[0].size() );
  assert( ncells.size() == ucell[0].size() );
  _ucell = ucell;
  _lat = lat;
  _ncells = ncells;
  _pbcs = pbcs;
  check_boundaries();
#ifndef MOLPRO
  std::string xyzfile = Input::sPars["periodic"]["xyzout"];
  if ( !xyzfile.empty() ) {
    xyz_supercell(xyzfile);
  }
#endif
}

void Periodic::check_boundaries()
{
#ifdef MOLPRO
  double tol = 1.e-6;
#else
  double tol = Input::fPars["periodic"]["thrdist"];
#endif
  // find sites on left-hand boundaries
  std::vector<std::vector<uint>> bounds(_ucell.nsites());
  assert( _lat.ndim() == _ncells.size() );
  for ( uint ld = 0; ld < _lat.ndim(); ++ld ) {
    Coords ortvec = _lat.orthogvec2others(ld);
    double len_ld = ortvec.dot(_lat.vec(ld));
    for ( uint is = 0; is < _ucell.nsites(); ++is ) {
      double dotprod = ortvec.dot(_ucell.coords(is));
      if ( std::abs(dotprod) < tol ) {
        // site is on the boundary in l dimension
        bounds[is].push_back(ld);
      } else if ( dotprod < -tol || dotprod > len_ld + tol ) {
        // the point is outside!
        xout << _ucell.coords(is);
        error(" The point is outside the lattice box!");
      }
    }
  }
  PSite ps(_lat.size());
  std::vector<uint> redundant(_ucell.nsites(),0);
  for ( uint is = 0; is < _ucell.nsites(); ++is ) {
    ps.site = is;
    std::vector<uint> & bound = bounds[is];
    // if site is on the boundary
    // check whether the site is non-unique
    // got through all combinations of lattices that apply
    for ( uint nlat = 1; nlat <= bound.size(); ++nlat ) {
      auto last = bound.begin()+nlat;
      do {
        ps.zero(false);
        for ( auto it = bound.begin(); it != last; ++it )
          ps[*it] = 1;
        Coords transl = site_coords(ps);
        // check whether it coinsides with any sites
        for ( uint is1 = 0; is1 < _ucell.nsites(); ++is1 ) {
          if ( transl.dist2(_ucell.coords(is1)) < tol ) {
            for ( auto it = bound.begin(); it != last; ++it )
              redundant[is1] |= 1 << *it;
          }
        }
      } while ( next_combination(bound.begin(),last,bound.end()) );
    }
    ps.zero(false);
//     xout << "Site: " << is << site_coords(ps) << ": " << redundant[is] << std::endl;
  }
  _ucell.set_redundancy(redundant);
  // calc number of sites in the supercell
  _nsites_in_supercell = 0;
  PSite s1(ndim());
  do {
    ++_nsites_in_supercell;
  } while ( next(s1) );
}

Coords Periodic::site_coords ( const PSite& isite ) const
{
  Coords scoords(_ucell.coords(isite.site));
  assert( isite.size() == _lat.ndim() );
  for ( uint l = 0; l < isite.size(); ++l ) {
    assert( scoords.size() == _lat.vec(l).size() );
    for ( uint i = 0; i < scoords.size(); ++i ) {
      scoords[i] += _lat.vec(l)[i]*isite[l];
    }
  }
  return scoords;
}

double Periodic::dist2(const PSite& ps1, const PSite& ps2) const
{
  Coords
    cps1 = site_coords(ps1),
    cps2 = site_coords(ps2);
  double dd = cps1.dist2(cps2);
  double dd1;
  PSite ps;
  std::vector<uint> perdims;
  for ( uint idim = 0; idim < _lat.ndim(); ++idim ) {
    if ( _pbcs[idim] != 0 ) perdims.push_back(idim);
  }
  for ( uint nlat = 1; nlat <= perdims.size(); ++nlat ) {
    auto last = perdims.begin()+nlat;
    do {
      for ( uint signs = 0; signs < (1u << nlat); ++signs ) {
        ps = ps2;
        uint pos = 0;
        for ( auto it = perdims.begin(); it != last; ++it, ++pos ) {
          if ( signs & (1 << pos) )
            ps[*it] += _ncells[*it];
          else
            ps[*it] -= _ncells[*it];
        }
        dd1 = cps1.dist2(site_coords(ps));
        if ( dd1 < dd ) dd = dd1;
      }
    } while ( next_combination(perdims.begin(),last,perdims.end()) );
  }
  return dd;
}


std::vector<double> Periodic::dist2neighbours(uint ndist) const
{
#ifdef MOLPRO
  double tol = 1.e-6;
#else
  double tol = Input::fPars["periodic"]["thrdist"];
#endif
  std::vector<double> res(ndist,1.e10);
  PSite s1(ndim()),s2(ndim());
  do {
    do {
      double dd = dist2(s1,s2);
      if ( dd < tol) continue;
      for (double& d2n: res) {
        if ( std::abs(dd-d2n) < tol ) {
          break;
        }
        if ( dd < d2n ) {
          std::swap(d2n,dd);
        }
      }
    } while ( next(s2) );
    s2.zero();
  } while ( next(s1) );
  s1.zero(); s2.zero();
  xout << "dist2neighbours "; for ( double d2n: res ) xout << d2n << " "; xout << std::endl;
  return res;
}

void Periodic::xyz_supercell(const std::string& xyzfilename) const
{
  std::ofstream xyzfile;
  xyzfile.open(xyzfilename);
  xyzfile << _nsites_in_supercell << std::endl;
  PSite ps(ndim());
  uint scsite = 1;
  do {
    xyzfile  << std::endl << _ucell.sitename(ps.site) << scsite << " ";
    site_coords(ps).print_xyz(xyzfile);
    ++scsite;
  } while ( next(ps) );
  xyzfile.close();
  xout << "XYZ coordinates written to " << xyzfilename << std::endl;
}

void Coords::print_xyz(std::ostream& o) const
{
  assert( size() <= 3 );
  for ( auto c: *this ) o << " " << c;
  for ( uint i = size(); i < 3; ++i ) o << " 0.0";
}

std::ostream & operator << (std::ostream& o, const Coords& co) {
  o << "{";
  for ( uint i = 0; i < co.size(); ++i ) {
    o << co[i];
    if ( i < co.size()-1 ) o << ",";
  }
  o << "}";
  return o;
}

std::ostream & operator << (std::ostream& o, const PSite& ps) {
  o << "{" << ps.site <<";";
  for ( uint i = 0; i < ps.size(); ++i ) {
    o << ps[i];
    if ( i < ps.size()-1 ) o << ",";
  }
  o << "}";
  return o;
}

std::ostream & operator << (std::ostream& o, const UCell& uc) {
  o << "{";
  for ( uint i = 0; i < uc.nsites(); ++i ) {
    o << uc[i];
    if ( i < uc.size()-1 ) o << ",";
  }
  o << "}";
  return o;
}

std::ostream & operator << (std::ostream& o, const Lattice& lat) {
  o << "{";
  for ( uint i = 0; i < lat.ndim(); ++i ) {
    o << lat[i];
    if ( i < lat.size()-1 ) o << ",";
  }
  o << "}";
  return o;
}

} //namespace HamDump
