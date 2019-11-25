/* -------------
 * | disco cat |
 * -------------
 * Discrete ordinates cartesian transport
 *
 * Gavin Ridley
 * 22.212 F2019
 */
#include <algorithm>
#include <iostream>
#include <experimental/filesystem>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include "lsquadrature.h"
using namespace std; // I love cluttered namespaces!
namespace fs = experimental::filesystem;

constexpr float EPS = 1e-8f;
constexpr float PI4 = M_PI * 4.0;

// --- A geometry class for the pedagogical fuel assembly in 22.212 ---
// This has been stripped down from my previous random ray solver.
struct SquarePinGeom
{
  unsigned mesh_dimx;
  float mesh_dx;
  array<unsigned, 6> index_endpoints;

  // Prescribed fuel dimensions:
  static constexpr float pitch = 1.2;
  static constexpr float assembly_width = 3.0 * pitch;
  static constexpr float assembly_radius = assembly_width / 2.0;
  static constexpr float pin_width = pitch / 3.0;

  public:
    SquarePinGeom(unsigned mesh_dimx);
    bool inside_fuel(unsigned i, unsigned j);
    bool inside_fuel(unsigned indx);
};
SquarePinGeom::SquarePinGeom(unsigned mesh_dimx) :
  mesh_dimx(mesh_dimx),
  mesh_dx(assembly_width / (float)mesh_dimx)
{
  if (mesh_dimx % 9 != 0)
  {
    cerr << "Mesh cell count in one direction must be divisible by 9. Got "
         << mesh_dimx << endl;
    exit(1);
  }
  unsigned fuelwide = mesh_dimx / 9;
  for (unsigned i=0; i<3; ++i)
  {
    index_endpoints[2*i] = fuelwide + i * mesh_dimx / 3;
    index_endpoints[2*i+1] = index_endpoints[2*i] + fuelwide - 1;
  }
}
bool SquarePinGeom::inside_fuel(unsigned i, unsigned j)
{
  // Loop over all 6 pin edge corner indices
  bool i_inside = false;
  bool j_inside = false;
  bool i_centered = false;
  bool j_centered = false;
  for (unsigned indx=0; indx < 3; ++indx)
  {
    if (i <= index_endpoints[2*indx+1] and i >= index_endpoints[2*indx]) i_inside=true;
    if (j <= index_endpoints[2*indx+1] and j >= index_endpoints[2*indx]) j_inside=true;
  }
  if (i <= index_endpoints[3] and i >= index_endpoints[2]) i_centered=true;
  if (j <= index_endpoints[3] and j >= index_endpoints[2]) j_centered=true;

  return (i_inside and j_inside) and not (i_centered and j_centered);
}
bool SquarePinGeom::inside_fuel(unsigned indx)
{
  // Calculate discrete cartesian coordinates
  unsigned i = indx / mesh_dimx;
  unsigned j = indx % mesh_dimx;
  return inside_fuel(i, j);
}

// General calculation settings
struct RunSettings
{
  string xslibrary;
  unsigned mesh_dimx;
  unsigned ngroups;
  unsigned s_n;
  unsigned maxiter;
  string quadrature_type;

  public:
    MomType quadType();
    RunSettings(string inputfile);
};
RunSettings::RunSettings(string inputfile)
{
  ifstream instream(inputfile);
  string word;
  while (instream >> word)
  {
    if (word == "xslibrary")
      instream >> xslibrary;
    else if (word == "mesh_dimx")
      instream >> mesh_dimx;
    else if (word == "ngroups")
      instream >> ngroups;
    else if (word == "s_n")
      instream >> s_n;
    else if (word == "quadrature_type")
      instream >> quadrature_type;
    else if (word == "maxiter")
      instream >> maxiter;
  }
  instream.close();
}
MomType RunSettings::quadType()
{
  if (quadrature_type == "ODD")
    return ODD;
  else if (quadrature_type == "EVEN")
    return EVEN;
  else
  {
    cerr << "unrecognized quadrature_type" << endl;
    exit(1);
  }
}

// Holds all of the macroscopic cross sections needed for a steady-state flux solution
struct Material
{
  string name;
  unsigned ngroups;
  bool fissile;
  bool diffusion; // if used in homogenized diffusion

  vector<float> trans, abs, nuscat, chi, nufiss, diff;

  static const array<const string, 3> xs_types;
  static const array<const string, 2> fiss_xs_types;

  public:
    Material(unsigned ngroups, bool fissile = false);
    void setFissile();

    // useful for homogenization
    void zeroEntries();
};
const array<const string, 3> Material::xs_types = {"trans", "abs", "nuscat"} ;
const array<const string, 2> Material::fiss_xs_types = {"chi", "nufiss"};
Material::Material(unsigned thisngroups, bool thisfissile) :
  ngroups(thisngroups),
  fissile(thisfissile),
  trans(ngroups, 0.0f),
  abs(ngroups, 0.0f),
  nuscat(ngroups*ngroups, 0.0f),
  chi(fissile ? ngroups : 0, 0.0f),
  nufiss(fissile ? ngroups : 0, 0.0f),
  diff(diffusion ? ngroups : 0, 0.0f)
{
}
void Material::setFissile()
{
  fissile = true;
  chi.resize(ngroups);
  nufiss.resize(ngroups);
}
void Material::zeroEntries()
{
  fill(trans.begin(), trans.end(), 0.0f);
  fill(abs.begin(), abs.end(), 0.0f);
  fill(nuscat.begin(), nuscat.end(), 0.0f);
  fill(chi.begin(), chi.end(), 0.0f);
  fill(nufiss.begin(), nufiss.end(), 0.0f);
  fill(diff.begin(), diff.end(), 0.0f);
}

// OK, so, everything representible on a computer is finite, so IDK why
// I put finite here. I suppose it means that macroscopics are not calculated
// on the fly, in contrast to what you'd do in depletion where you store micro
// XS and calculate macro from material density
class FiniteMaterialSet
{
  unsigned nmaterials;
  vector<Material> materials;
  map<string, unsigned> material_map;

  void loadVector(vector<float>& to_vec, fs::path infile);
  static unsigned getMaterialCount(string libname);
  public:

    const Material& getMaterial(string name);

    FiniteMaterialSet(string xslib, unsigned ngroups);
};
void FiniteMaterialSet::loadVector(vector<float>& to_vec, fs::path infile)
{
  // Checks correct number XS loaded
  unsigned loadCount = 0;
  ifstream instream(infile, ifstream::in);
  if (not instream.good())
  {
    cerr << "cannot load " << infile << endl;
    exit(1);
  }
  float value;
  while (instream >> value)
  {
    to_vec[loadCount++] = value;
    if (loadCount > to_vec.size())
    {
      cerr << "Tried to load too many XS from material " << infile << endl;
      cerr << "too many groups or too few?" << endl;
      exit(1);
    }
  }
  if (loadCount != to_vec.size())
  {
    cerr << "too few xs values in " << infile << endl;
    exit(1);
  }
}
FiniteMaterialSet::FiniteMaterialSet(string xslib, unsigned ngroups) :
  nmaterials(getMaterialCount(xslib)),
  materials(nmaterials, Material(ngroups))
{
  if (nmaterials == 0)
  {
    cout << "zero materials were found in xslib named: " << xslib << endl;
    exit(1);
  }

  unsigned mat_indx = 0;
  fs::path p(xslib);
  for (const auto& entry : fs::directory_iterator(p))
  {
    // Load required XS
    string materialname = entry.path().filename();
    cout << "Processing material " << materialname << endl;
    Material& mat = materials[mat_indx];
    for (string xs_type : Material::xs_types)
    {
      loadVector(mat.trans, entry/"trans");
      loadVector(mat.abs, entry/"abs");
      loadVector(mat.nuscat, entry/"nuscat");
    }

    // Maybe load fissile XS
    vector<bool> fissile_xs_present(Material::fiss_xs_types.size(), false);
    unsigned fiss_i = 0;
    for (string fiss_xs : Material::fiss_xs_types)
      if (fs::exists(entry/fiss_xs)) fissile_xs_present[fiss_i++] = true;
    bool no_fiss = none_of(fissile_xs_present.begin(), fissile_xs_present.end(),
        [](bool x){ return x; });
    bool all_fiss = any_of(fissile_xs_present.begin(), fissile_xs_present.end(),
        [](bool x){ return x; });
    if (not no_fiss ^ all_fiss) { cerr << "Some, but not all fiss. XS found. " << endl; exit(1); }
    if (all_fiss)
    {
      mat.setFissile();
      for (string fiss_xs : Material::fiss_xs_types)
      {
        loadVector(mat.chi, entry/"chi");
        loadVector(mat.nufiss, entry/"nufiss");
      }
    }

    pair<string, unsigned> mat_dict_entry(materialname, mat_indx++);
    material_map.insert(mat_dict_entry);
  }

}
const Material& FiniteMaterialSet::getMaterial(string name) { return materials[material_map[name]]; }
unsigned FiniteMaterialSet::getMaterialCount(string libname)
{
  fs::path p(libname);
  string filename;
  unsigned nmaterials = 0;

  // check all required cross sections present in each material
  for (const auto& entry : fs::directory_iterator(p))
  {
    if (fs::is_directory(entry))
    {
      for (string xs_type : Material::xs_types)
        if (not fs::exists(entry/xs_type))
        {
          cerr << "Required cross section " << xs_type <<
            " not found in material " << entry << endl;
          exit(1);
        }
      ++nmaterials;
    }
    else
    {
      cerr << "Found non-directory file in xslib." << endl;
      cerr << "That shouldn't be there. Name: " << entry << endl;
      exit(1);
    }
  }
  return nmaterials;
}

// OK, here is the code for my Sn solver
class Solver2D
{
  RunSettings settings;
  SquarePinGeom geom;
  LSQuadrature quad;
  FiniteMaterialSet materialSet;

  unsigned noct; // angular unknowns per octant
  unsigned ngroups;
  unsigned mesh_dimx;
  unsigned stride; // stride in boundary array between different rays
  float dx; // mesh spacing
  vector<float> fluxes; // scalar flux
  vector<float> source; // isotropic cell source

  // boundary fluxes
  vector<float> left_fluxes_fwd; // eta > 0
  vector<float> left_fluxes_bwd; // eta < 0
  vector<float> bottom_fluxes_fwd; // mu > 0
  vector<float> bottom_fluxes_bwd; // mu < 0
  vector<float> top_fluxes_fwd; // mu > 0
  vector<float> top_fluxes_bwd; // mu < 0
  vector<float> right_fluxes_fwd; // eta > 0
  vector<float> right_fluxes_bwd; // eta < 0

  // Core S_n kernel
  void processCell(unsigned row, unsigned col,
                   Ray ray,
                   vector<float>& side_from,
                   vector<float>& vert_from);

  // stuff for the built-in diffusion solver
  vector<float> diffusion_fluxes;
  vector<Material> homogenized_materials;

  public:
    Solver2D(RunSettings settings_a);

    void zeroScalarFlux();

    // Sets the source in a group-cell index
    void setSource(unsigned indx, float src);

     // Single fixed source sweep, return true if converged
    void sweepSource();

    // Result saving methods
    void dumpFluxes(string fname);

    // Calculate scattering source and add to source
    void scatter();

    void normalizeFlux();

    // Calculate fission source, and add to source. Returns integral
    // fission source
    float fission(float k);

    // guess a flat flux
    void setFlatFlux();

    void zeroSource();

    // homogenize materials in each unit cell for diffusion solve
    void homogenizeCells();
};
Solver2D::Solver2D(RunSettings settings_a) :
  settings(settings_a),
  geom(settings.mesh_dimx),
  quad(settings.s_n, settings.quadType()),
  materialSet(settings.xslibrary, settings.ngroups),
  noct(quad.nOct()),
  ngroups(settings.ngroups),
  mesh_dimx(settings.mesh_dimx),
  stride(mesh_dimx * ngroups),
  dx(geom.assembly_width / mesh_dimx),
  fluxes(mesh_dimx * mesh_dimx * ngroups),
  source(fluxes.size()),
  left_fluxes_fwd(stride * noct),
  left_fluxes_bwd(stride * noct),
  bottom_fluxes_fwd(stride * noct),
  bottom_fluxes_bwd(stride * noct),
  top_fluxes_fwd(stride * noct),
  top_fluxes_bwd(stride * noct),
  right_fluxes_fwd(stride * noct),
  right_fluxes_bwd(stride * noct),
  diffusion_fluxes(9 * ngroups),
  homogenized_materials(9, ngroups)
{
  for (auto& mat: homogenized_materials) mat.setFissile();
}
void Solver2D::homogenizeCells()
{
  // mesh_dimx is guaranteed to be divisible by 9, which makes
  // this work nicely
  constexpr unsigned nrow = 3;
  unsigned cell_width = mesh_dimx / nrow;

  for (unsigned i=0; i<nrow; ++i)
    for (unsigned j=0; j<nrow; ++j)
    {
      unsigned indx = i*nrow + j;

      // Lower left fine cell indices of this unit cell:
      unsigned min_x = i * cell_width;
      unsigned min_y = j * cell_width;

      // ensure that material's data is zeroed out
      Material& mat = homogenized_materials[indx];
      mat.zeroEntries();

      // chi is always the same as the fuel
      mat.chi = materialSet.getMaterial("fuel").chi;

      // loop over fine mesh cells:
      vector<float> groupflux_integral(ngroups);

      for (unsigned ii=0; ii<cell_width; ++ii)
        for (unsigned jj=0; jj<cell_width; ++jj)
        {
          unsigned fine_i = min_x + jj;
          unsigned fine_j = min_y + ii;
          unsigned fine_indx = fine_i * mesh_dimx + fine_j;

          if (geom.inside_fuel(fine_i, fine_j))
            Material& thismat = materialSet.getMaterial("fuel");
          else
            Material& thismat = materialSet.getMaterial("mod");

          // loop over groups
          for (unsigned g=0; g<ngroups; ++g)
          {
            float groupflux = fluxes[ngroups*(i * mesh_dimx + j) + g];
            groupflux_integral[g] += groupflux;
            mat.diff[g] += groupflux / (3.0f * thismat.trans[g]);
            mat.trans[g] += groupflux * thismat.trans[g];
            mat.nufiss[g] += groupflux * thismat. // TODO
          }
        }

      // Now divide by flux integrals on each group constant
    }
}
void Solver2D::normalizeFlux()
{
  // calculate norm of flux:
  float norm = 0.0;
  for (auto x: fluxes) norm += abs(x);

  // divide volumetric fluxes by norm:
  for (auto& x: fluxes) x /= norm;

  // divide boundary fluxes by norm:
  // (Case study of why C++ is better than Fortran)
  vector<vector<float>*> reflist = {
    &left_fluxes_fwd,
    &left_fluxes_bwd,
    &bottom_fluxes_fwd,
    &bottom_fluxes_bwd,
    &top_fluxes_fwd,
    &top_fluxes_bwd,
    &right_fluxes_fwd,
    &right_fluxes_bwd };

  for (auto ref: reflist)
    for (auto& x: *ref)
      x /= norm;
}
void Solver2D::zeroSource() { for (auto& x: source) x = 0.0; }
void Solver2D::setFlatFlux() { for (auto& x: fluxes) x = 1.0; }
void Solver2D::scatter()
{
  for (unsigned fsr=0; fsr<mesh_dimx*mesh_dimx; ++fsr)
  {
    string mat_name;
    if (geom.inside_fuel(fsr))  
      mat_name = "fuel";
    else
      mat_name = "mod";
    const Material& mat = materialSet.getMaterial(mat_name);
    const vector<float>& scatmat = mat.nuscat;
    for (unsigned g=0; g<ngroups; ++g)
    {
      source[ngroups * fsr + g] = 0.0;
      for (unsigned gprime=0; gprime<ngroups; ++gprime)
      {
        source[ngroups * fsr + g] += 
          scatmat[g*ngroups + gprime] * fluxes[ngroups * fsr + gprime] / PI4;
      }
    }
  }
} 
float Solver2D::fission(float k)
{
  float fissionSource = 0.0;
  for (unsigned fsr=0; fsr<mesh_dimx*mesh_dimx; ++fsr)
  {
    string mat_name;
    if (geom.inside_fuel(fsr)) 
      mat_name = "fuel";
    else
    {
      mat_name = "mod";
      continue;
    }
    const Material& mat = materialSet.getMaterial(mat_name); 
    const vector<float>& nusigf = mat.nufiss;
    const vector<float>& chi = mat.chi;
    // NOTE could be done more efficiently
    for (unsigned g=0; g<ngroups; ++g)
    {
      for (unsigned gprime=0; gprime<ngroups; ++gprime)
      {
        float this_fiss = chi[g] * fluxes[ngroups * fsr + gprime] * nusigf[gprime] / k;
        source[ngroups * fsr + g] += this_fiss / PI4;
        fissionSource += this_fiss;
      }
    }
  }
  return fissionSource;
}
void Solver2D::dumpFluxes(string fname)
{
  ofstream f(fname, ofstream::out);
  for (auto flx : fluxes) f << flx << endl;
  f.close();
}
void inline Solver2D::processCell(unsigned row, unsigned col,                                       
                        Ray ray,
                        vector<float>& side_from,                             
                        vector<float>& vert_from)         
{                                                     
  unsigned space_indx = row*mesh_dimx + col;
  // Get material cross section                      
  string mat_name;              
  if (geom.inside_fuel(space_indx))                  
    mat_name = "fuel";                               
  else                                        
    mat_name = "mod";                   
  const Material& mat = materialSet.getMaterial(mat_name);
  const vector<float>& sigt = mat.trans;                            
  unsigned flux_indx;                   
  float flux; // intermediate group flux
  
  // Loop on groups goes innermost
  for (unsigned g=0; g<ngroups; ++g)
  {
    // Calc. cell center flux:      
    flux_indx = ngroups * space_indx + g;
    flux = (source[flux_indx] + 2.0 * ray.mu / dx * side_from[g] +
           2.0 * ray.eta / dx * vert_from[col * ngroups + g] ) / 
           (sigt[g] + 2.0 * ray.mu / dx + 2.0 * ray.eta / dx);
                                                                
    // Calc flux of vertically next edge                        
    vert_from[col*ngroups+g] = 2.0 * flux - vert_from[col*ngroups+g];
                                        
    // Calc flux on the side's next edge:                                     
    side_from[g] = 2.0 * flux - side_from[g];                                
                                                                              
    // Add angular flux to total flux:                                     
    // Factor of two from z symmetry                                       
    fluxes[space_indx*ngroups+g] += ray.wgt * flux * 2.0;                     
  }                                                                        
}
void Solver2D::zeroScalarFlux() { fill(fluxes.begin(), fluxes.end(), 0.0); }
void Solver2D::setSource(unsigned indx, float src)
{
  if (indx < source.size()) source[indx] = src;
  else cout << "warn: attempt to set source out of bounds" << endl;
}
void Solver2D::sweepSource()
{
  // Current direction cosines and ray weight:
  unsigned row, col, g; // cell indices
  unsigned ray_id; // enumeration of rays in octant
  Ray ray;

  // Zero flux out, then add in angular components back 
  // piece by piece during the sweep
  zeroScalarFlux();

  // Temporary fluxes
  vector<float> vertical_tmp(ngroups * mesh_dimx); // below when going up, above going down
  vector<float> side_tmp(ngroups); // left when going right, right when going left

  // quadrant 1
  ray_id = 0;
  while (quad.iterateOctant(ray))
  {
    /*   sweep this way              
     *
     *     /o/    
     *      |            
     *     / \
     */
    copy(bottom_fluxes_fwd.begin()+ray_id*stride,
         bottom_fluxes_fwd.begin()+(ray_id+1)*stride,
         vertical_tmp.begin()); // apply lower BC
    for (row=0; row<mesh_dimx; ++row)
    {
      // Copy left boundary flux to lateral temporary fluxes
      copy(left_fluxes_fwd.begin()+row*ngroups+ray_id*stride,
           left_fluxes_fwd.begin()+(row+1)*ngroups+ray_id*stride,
           side_tmp.begin());

      // Sweep rightward
      for (col=0; col<mesh_dimx; ++col)
        processCell(row, col, ray, side_tmp, vertical_tmp);

      // Save right edge flux to boundary flux
      for (g=0; g<ngroups; ++g)
        right_fluxes_fwd[ray_id*stride+row*ngroups+g] = side_tmp[g];
    }
    // Save top flux result to top boundary flux
    copy(vertical_tmp.begin(),
         vertical_tmp.begin()+stride,
         top_fluxes_fwd.begin()+ray_id*stride);

    /*   sweep this way              
     *
     *     \o\
     *      |            
     *     / \
     */
    copy(bottom_fluxes_bwd.begin()+ray_id*stride,
         bottom_fluxes_bwd.begin()+(ray_id+1)*stride,
         vertical_tmp.begin()); // apply lower BC
    for (row=0; row<mesh_dimx; ++row)
    {
      // Copy right boundary flux to lateral temporary fluxes
      copy(right_fluxes_fwd.begin()+row*ngroups+ray_id*stride,
           right_fluxes_fwd.begin()+(row+1)*ngroups+ray_id*stride,
           side_tmp.begin());

      // Sweep leftward 
      col = mesh_dimx;
      while (col --> 0)
        processCell(row, col, ray, side_tmp, vertical_tmp);

      // Save left edge flux to boundary flux
      for (g=0; g<ngroups; ++g)
        left_fluxes_fwd[ray_id*stride+row*ngroups+g] = side_tmp[g];
    }
    // Save top flux result to top boundary flux
    copy(vertical_tmp.begin(),
         vertical_tmp.begin()+stride,
         top_fluxes_bwd.begin()+ray_id*stride);

    /*   sweep this way              
     *
     *      o 
     *     \|\
     *     / \
     */
    copy(top_fluxes_fwd.begin()+ray_id*stride,
         top_fluxes_fwd.begin()+(ray_id+1)*stride,
         vertical_tmp.begin()); // apply upper BC
    row = mesh_dimx;
    while (row --> 0)
    {
      // Copy left boundary flux to lateral temporary fluxes
      copy(left_fluxes_bwd.begin()+row*ngroups+ray_id*stride,
           left_fluxes_bwd.begin()+(row+1)*ngroups+ray_id*stride,
           side_tmp.begin());

      // Sweep rightward
      for (col=0; col<mesh_dimx; ++col)
        processCell(row, col, ray, side_tmp, vertical_tmp);

      // Save right edge flux to boundary flux
      for (g=0; g<ngroups; ++g)
        right_fluxes_bwd[ray_id*stride+row*ngroups+g] = side_tmp[g];
    }
    // Save top flux result to top boundary flux
    copy(vertical_tmp.begin(),
         vertical_tmp.begin()+stride,
         bottom_fluxes_fwd.begin()+ray_id*stride);

    /*   sweep this way              
     *
     *      o 
     *     /|/
     *     / \
     */
    copy(top_fluxes_bwd.begin()+ray_id*stride,
         top_fluxes_bwd.begin()+(ray_id+1)*stride,
         vertical_tmp.begin()); // apply upper BC
    row = mesh_dimx;
    while (row --> 0)
    {
      // Copy left boundary flux to lateral temporary fluxes
      copy(right_fluxes_bwd.begin()+row*ngroups+ray_id*stride,
           right_fluxes_bwd.begin()+(row+1)*ngroups+ray_id*stride,
           side_tmp.begin());

      // Sweep rightward
      col = mesh_dimx;
      while (col --> 0)
        processCell(row, col, ray, side_tmp, vertical_tmp);

      // Save right edge flux to boundary flux
      for (g=0; g<ngroups; ++g)
        left_fluxes_bwd[ray_id*stride+row*ngroups+g] = side_tmp[g];
    }
    // Save top flux result to top boundary flux
    copy(vertical_tmp.begin(),
         vertical_tmp.begin()+stride,
         bottom_fluxes_bwd.begin()+ray_id*stride);

    ray_id++;
  }
}

int main(int argc, char * argv[])
{
  // Get command line argument to input file
  if (argc != 2)                                                   
  {                                      
    cerr << "Should pass a settings file as a command line argument" << endl;
    exit(1);
  }                 
  string filename = argv[1];
  RunSettings settings(filename);

  Solver2D solver(settings);

  // Do power iteration:
  float k = 1.0;
  float fissionSource=1.0;
  float oldFissionSource, fiss_quo;
  solver.setFlatFlux();
  for (unsigned n=0; n<settings.maxiter; ++n)
  {
    // recalculate source
    solver.zeroSource();
    solver.scatter();
    oldFissionSource = fissionSource;
    fissionSource = solver.fission(k);

    // calculate k
    fiss_quo = fissionSource / oldFissionSource;
    k *= fiss_quo;
    cout << "k = " << k << endl;

    // check convergence
    if (abs(fiss_quo-1.0) < EPS) break;

    solver.sweepSource();
  }

  solver.dumpFluxes("flux");
}
