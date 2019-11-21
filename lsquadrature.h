// Holds a level set quadrature class, which may match
// either purely even moments or both odd and even moments.
// The ones that match both tend to have negative weights
// after n=12, which is numerically undesirable.
#include <array>

enum MomType {EVEN, ODD};

template <unsigned n, MomType m>
class LSQuadrature
{
  array<float, n> mu;
  array<float, n> level_weights; // 1D quadrature weights
  array<float, n> point_weights; // 3D quadrature weights
  public:
    LSQuadrature();
}
