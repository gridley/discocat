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

// Using Eigen for the diffusion solves to save some coding
// Also, dense b/c it's not even worth going sparse w/ a 9x9
#include "Eigen/Dense"

constexpr float EPS = 1e-12f;
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
  static constexpr float assembly_width = 3.0f * pitch;
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
    void setDiffusion();

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
void Material::setDiffusion()
{
  diffusion = true;
  diff.resize(ngroups);
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
// Print material nicely:
ostream& operator<<(ostream& os, const Material& mat)
{
  os << "Material definition" << endl;
  os << "-------------------" << endl;
  os << "Abs. xs:" << endl;
  for (auto x: mat.abs) os << x << " ";
  os << endl;
  os << "Nuscat matrix:" << endl;
  for (unsigned g=0; g<mat.ngroups; ++g)
  {
    for (unsigned gp=0; gp<mat.ngroups; ++gp)
      os << mat.nuscat[mat.ngroups*g+gp] << " ";
    os << endl;
  }
  os << endl;
  os << "chi:" << endl;
  for (auto x: mat.chi ) os << x << " ";
  os << endl;
  os << "nufiss:" << endl;
  for (auto x: mat.nufiss) os << x << " ";
  os << endl;
  os << "diff coeff:" << endl;
  for (auto x: mat.nufiss) os << x << " ";
  os << endl;
  return os;
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
  vector<float> diffusion_source;
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

    // guess a flat fission source
    void setFlatSource();

    void zeroSource();

    void printPeakingFactors();

    // DIFFUSION SHI
    void homogenizeCells();
    void normalizeDiffusionFlux();
    void diffusionScatter();
    float diffusionFission(float k);
    void solveCoarseDiffusion();
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
  diffusion_source(9 * ngroups),
  homogenized_materials(9, ngroups)
{
  for (auto& mat: homogenized_materials)
  {
    // Assume homogenized materials store both a diffusion coefficient
    // and fission data by default
    mat.setFissile();
    mat.setDiffusion();
  }
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

          string mat_name;
          if (geom.inside_fuel(fine_indx))  
            mat_name = "fuel";
          else
            mat_name = "mod";
          const Material& thismat = materialSet.getMaterial(mat_name);

          // loop over groups
          for (unsigned g=0; g<ngroups; ++g)
          {
            float groupflux = fluxes[ngroups*fine_indx + g];
            groupflux_integral[g] += groupflux;
            mat.diff[g] += groupflux / (3.0f * thismat.trans[g]);
            mat.abs[g] += groupflux * thismat.abs[g];
            if (thismat.fissile)
              mat.nufiss[g] += groupflux * thismat.nufiss[g];

            // For the scattering matrix, the columns are homogenized
            // according to the flux
            for (unsigned gprime=0; gprime<ngroups; ++gprime)
              mat.nuscat[gprime*ngroups+g] += groupflux * thismat.nuscat[gprime*ngroups+g];
          }
        }

      // Now divide by flux integrals on each group constant
      for (unsigned g=0; g<ngroups; ++g)
      {
        mat.diff[g] /= groupflux_integral[g];
        mat.abs[g] /= groupflux_integral[g];
        mat.nufiss[g] /= groupflux_integral[g];
        for (unsigned gprime=0; gprime<ngroups; ++gprime)
          mat.nuscat[gprime*ngroups+g] /= groupflux_integral[g];
      }
    }

  // DBG print out obtained homogenized materials
  // for (auto& mat: homogenized_materials)
  //   cout << mat << endl;
}
void Solver2D::solveCoarseDiffusion()
{
  // guess flux to fission group
  constexpr unsigned nx = 3;
  constexpr unsigned ny = 3;
  constexpr unsigned n = nx * ny;

  // area of each pin cell
  float cell_width = geom.assembly_width / 3.0f;
  float cell_area = cell_width * cell_width;

  // Create diffusion matrices for each group (row major)
  vector<Eigen::Matrix<float,n,n>> diff_matrices(ngroups);
  for (unsigned g=0; g<ngroups; ++g)
    diff_matrices[g] = Eigen::Matrix<float,n,n>::Zero();

  // Create each matrix
  for (unsigned g=0; g<ngroups; ++g)
  {
    Eigen::Matrix<float,n,n>& mat = diff_matrices[g];
    for (unsigned i=0; i<ny; ++i)
      for (unsigned j=0; j<nx; ++j)
      {
        // get neighbor indices
        unsigned me = i * nx + j;
        unsigned left = me-1;
        unsigned right = me+1;
        unsigned top = me + nx;
        unsigned bottom = me - nx;

        // handle boundary conditions
        if (j == 0) left += nx;
        if (j == nx-1) right -= nx;
        if (i == 0) bottom += n;
        if (i == ny-1) top -= n;

        // calculate effective diffusion coefficients
        float dthis, dtop, dleft, dright, dbottom;
        dthis = homogenized_materials[me].diff[g];
        dtop = homogenized_materials[top].diff[g];
        dleft = homogenized_materials[left].diff[g];
        dright = homogenized_materials[right].diff[g];
        dbottom = homogenized_materials[bottom].diff[g];

        dtop *= 2.0f * dthis / (dtop + dthis);
        dleft *= 2.0f * dthis / (dleft + dthis);
        dright *= 2.0f * dthis / (dright + dthis);
        dbottom *= 2.0f * dthis / (dbottom + dthis);

        // put diffusion entries to matrix
        mat(me,me) += dtop + dleft + dright + dbottom;
        mat(me,left) -= dleft;
        mat(me,right) -= dright;
        mat(me,top) -= dtop;
        mat(me,bottom) -= dbottom;

        // Diffusion stuff needs to be divided by cell area
        for (auto col: {me, left, right, top, bottom}) mat(me,col) /= cell_area;

        // add removal coefficient to diagonal
        mat(me,me) += homogenized_materials[me].abs[g];
      }
  }

  // Pre-calculate some pivoted QR decompositions
  vector<Eigen::ColPivHouseholderQR<Eigen::Matrix<float,n,n>>> diff_qrs;
  for (auto& m: diff_matrices)
    diff_qrs.emplace_back(m);

  // check for rank deficiency. Should work with this factorization anyways tho
  for (unsigned g=0; g<ngroups; ++g)
  {
    if (diff_qrs[g].rank() < n)
    {
      cout << "Warning: group " << g << " has a rank deficient matrix" << endl;
      cout << diff_qrs[g].rank() << endl << endl;
    }
  }

  // solve multigroup diffusion using power iteration
  fill(diffusion_fluxes.begin(), diffusion_fluxes.end(), 1.0f);
  float oldFissionSource, fissionSource, fiss_quo, k;
  k = 1.0f;
  for (unsigned n=0; n<settings.maxiter; ++n)
  {
    normalizeDiffusionFlux();

    fill(diffusion_source.begin(), diffusion_source.end(), 0.0f);
    diffusionScatter();
    oldFissionSource = diffusionFission(k);

    // Update fluxes, given the previously calculated source
    Eigen::Matrix<float, 9, 1> rhs, newflux;

    // Loop over groups, constructing RHS and solving for new flux
    for (unsigned g=0; g<ngroups; ++g)
    {
      for (unsigned i=0; i<9; ++i)
        rhs(i) = diffusion_source[i*ngroups+g];
      newflux = diff_qrs[g].solve(rhs);
      for (unsigned i=0; i<9; ++i)
        diffusion_fluxes[i] = newflux(i);
    }

    fill(diffusion_source.begin(), diffusion_source.end(), 0.0f);
    diffusionScatter();
    fissionSource = diffusionFission(k);

    // calculate k
    fiss_quo = fissionSource / oldFissionSource;
    k *= fiss_quo;
    cout << "\r" << "k = " << k << "   " << "(" << n+1 << "/" <<
      settings.maxiter << ")" << flush;
    if (abs(fiss_quo-1.0f) < EPS)
    {
      cout << endl << "k convergence detected. stopping" << endl;
      break;
    }
  }
  cout << endl;
  // for (auto x: diffusion_fluxes) cout << x << endl;
}
void Solver2D::normalizeDiffusionFlux()
{
  float norm = 0.0f;
  for (auto x: diffusion_fluxes) norm += abs(x);
  for (auto& x: diffusion_fluxes) x /= norm;
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
void Solver2D::zeroSource() { for (auto& x: source) x = 0.0f; }
void Solver2D::setFlatSource()
{
  // Set source to be only in top fission group
  for (unsigned i=0; i < mesh_dimx*mesh_dimx; ++i)
    for (unsigned g=0; g<ngroups; ++g)
      if (g==0)
        source[i * ngroups] = 1.0f;
}
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
void Solver2D::diffusionScatter()
{
  for (unsigned fsr=0; fsr<9; ++fsr)
  {
    const Material& mat = homogenized_materials[fsr];
    const vector<float>& scatmat = mat.nuscat;
    for (unsigned g=0; g<ngroups; ++g)
    {
      diffusion_source[ngroups * fsr + g] = 0.0;
      for (unsigned gprime=0; gprime<ngroups; ++gprime)
      {
        diffusion_source[ngroups * fsr + g] += 
          scatmat[g*ngroups + gprime] * diffusion_fluxes[ngroups * fsr + gprime];
      }
    }
  }
} 
float Solver2D::fission(float k)
{
  float fissionSource = 0.0f;
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
float Solver2D::diffusionFission(float k)
{
  float fissionSource = 0.0f;
  for (unsigned fsr=0; fsr<9; ++fsr)
  {
    const Material& mat = homogenized_materials[fsr];
    const vector<float>& nusigf = mat.nufiss;
    const vector<float>& chi = mat.chi;
    // NOTE could be done more efficiently
    for (unsigned g=0; g<ngroups; ++g)
    {
      for (unsigned gprime=0; gprime<ngroups; ++gprime)
      {
        float this_fiss = chi[g] * diffusion_fluxes[ngroups * fsr + gprime] * nusigf[gprime] / k;
        diffusion_source[ngroups * fsr + g] += this_fiss;
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
    flux = (dx*source[flux_indx] + 2.0f * ray.mu * side_from[g] +
           2.0f * ray.eta * vert_from[col * ngroups + g] ) / 
           (dx * sigt[g] + 2.0f * ray.mu + 2.0f * ray.eta);
                                                                
    // Calc flux of vertically next edge                        
    vert_from[col*ngroups+g] = 2.0f * flux - vert_from[col*ngroups+g];
                                        
    // Calc flux on the side's next edge:                                     
    side_from[g] = 2.0f * flux - side_from[g];                                
                                                                              
    // Add angular flux to total flux:                                     
    // Factor of two from z symmetry                                       
    fluxes[flux_indx] += ray.wgt * flux * 2.0f;                     
  }                                                                        
}
void Solver2D::printPeakingFactors()
{
  // The assumption here is that the nu value is
  // pretty much constant throughout the problem
  float cornerFission, edgeFission;
  cornerFission = 0.0f;
  edgeFission = 0.0f;
  unsigned cellwide = mesh_dimx / 3;

  // get corner pin fission rate
  for (unsigned i=0; i<cellwide; ++i)
    for (unsigned j=0; j<cellwide; ++j)
    {
      unsigned fsr = i*mesh_dimx + j;
      string mat_name;              
      if (geom.inside_fuel(fsr))                  
        mat_name = "fuel";                               
      else                                        
        mat_name = "mod";                   
      const Material& mat = materialSet.getMaterial(mat_name);
      if (not mat.fissile) continue;
      const vector<float>& fiss = mat.nufiss;

      for (unsigned g=0; g<ngroups; ++g)
        cornerFission += fluxes[fsr*ngroups+g]*fiss[g];
    }

  // edge pin fission rate
  for (unsigned i=0; i<cellwide; ++i)
    for (unsigned j=cellwide; j<2*cellwide; ++j)
    {
      unsigned fsr = i*mesh_dimx + j;
      string mat_name;              
      if (geom.inside_fuel(fsr))                  
        mat_name = "fuel";                               
      else                                        
        mat_name = "mod";                   
      const Material& mat = materialSet.getMaterial(mat_name);
      if (not mat.fissile) continue;
      const vector<float>& fiss = mat.nufiss;

      for (unsigned g=0; g<ngroups; ++g)
        edgeFission += fluxes[fsr*ngroups+g]*fiss[g];
    }

  // get peaking factors:
  float avg = (cornerFission + edgeFission)/2.0f;
  cout << "Corner pin peaking factor is " << cornerFission/avg << endl;
  cout << "Edge pin peaking factor is " << edgeFission/avg << endl;
}
void Solver2D::zeroScalarFlux() { fill(fluxes.begin(), fluxes.end(), 0.0f); }
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
     *  (definitely not a dab)
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
  solver.setFlatSource();
  solver.sweepSource();
  for (unsigned n=0; n<settings.maxiter; ++n)
  {
    solver.normalizeFlux();
    solver.zeroSource();
    solver.scatter();
    oldFissionSource = solver.fission(k);
    solver.sweepSource();

    solver.zeroSource();
    solver.scatter();
    fissionSource = solver.fission(k);

    // calculate k
    fiss_quo = fissionSource / oldFissionSource;
    k *= fiss_quo;
    cout << "\r" << "k = " << k << "   " << "(" << n+1 << "/" <<
      settings.maxiter << ")" << flush;
    if (abs(fiss_quo-1.0f) < EPS)
    {
      cout << endl << "k convergence detected. stopping" << endl;
      break;
    }
  }
  cout << endl;
  solver.printPeakingFactors();

  // Homogenize group constants for diffusion solve:
  solver.homogenizeCells();

  cout << "Diffusion solve:" << endl;
  solver.solveCoarseDiffusion();

  solver.dumpFluxes("flux");
}
