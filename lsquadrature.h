#pragma once
// Holds a level set quadrature class, which may match
// either purely even moments or both odd and even moments.
// The ones that match both tend to have negative weights
// after n=12, which is numerically undesirable.
#include <array>
#include <iostream>
#include <map>
#include <utility>
#include <cmath>
#include <vector>
#include <algorithm>

enum MomType {EVEN, ODD};

// Stores data for a given angular direction
struct Ray
{
  float mu, eta, theta, wgt;
  Ray(float mu_a, float eta_a, float theta_a, float wgt_a);
  Ray();
};

class LSQuadrature
{
  // Quadrature type
  MomType quadtype;

  // These keep track of angles when iterating over octants
  unsigned i, j;

  // Number of distinct cosines:
  unsigned n;

  std::vector<float> mu; // available cosine
  std::vector<float> weights; // 1D quadrature weights
  std::vector<float> point_weights; // 3D quadrature weights

  // Map i,j,k indices into point weights.
  // This is a kinda hard combinatorial problem which is easily
  // solvable by brute force then storing the solution.
  std::map<std::array<unsigned, 3>, unsigned> indices2wgt;

  // Fills entries based on template values
  void fillEntries();

  // Generate point weight map from symmetry considerations
  void generateMap();

  // Change from weight normalization used by Lathrop and Carlson
  // into 4 pi normalization
  void normalize4Pi();

  public:
    LSQuadrature(unsigned na, MomType ma);

    // True until all points in the octant have been covered
    bool iterateOctant(float& xi1, float& xi2, float& xi3, float& wgt);
    bool iterateOctant(Ray& ray);

    // Calculate number of unknowns per octant
    unsigned nOct();
};
