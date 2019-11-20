/* -------------
 * | disco cat |
 * -------------
 * Discrete ordinates cartesian transport
 *
 * Gavin Ridley
 * 22.212 F2019
 */
#include <iostream>
#include <experimental/filesystem>
#include <fstream>
#include <string>
#include <vector>
using namespace std; // I love cluttered namespaces!

// --- A geometry class for the pedagogical fuel assembly in 22.212 ---
// This has been stripped down from my previous random ray solver.
// (S_n is easier than random ray IMO, why are we doing this second?)
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

  public:
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
  }
  instream.close();
}

// Holds all of the macroscopic cross sections needed for a steady-state flux solution
struct Material
{
  string name;
  unsigned ngroups;
  bool fissile;

  vector<float> trans, abs, nuscat, chi, nufiss;

  static const array<const string, 3> xs_types;
  static const array<const string, 2> fiss_xs_types;

  public:
    Material(unsigned ngroups, bool fissile = false);
    void setFissile();
};
const array<const string, 3> Material::xs_types = {"trans", "abs", "nuscat"} ;
const array<const string, 2> Material::fiss_xs_types = {"chi", "nufiss"};
Material::Material(unsigned thisngroups, bool thisfissile) :
  ngroups(thisngroups),
  fissile(thisfissile),
  trans(ngroups),
  abs(ngroups),
  nuscat(ngroups*ngroups),
  chi(fissile ? ngroups : 0),
  nufiss(fissile ? ngroups : 0)
{
}
void Material::setFissile()
{
  fissile = true;
  chi.resize(ngroups);
  nufiss.resize(ngroups);
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

// OK, here is the unique code for my Sn solver. It's quite simple!
class Solver2D
{
  vector<float> fluxes;
}

int main()
{
  SquarePinGeom geom(9);
}
