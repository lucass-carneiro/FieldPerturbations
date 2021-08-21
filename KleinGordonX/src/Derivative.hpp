#ifndef DERIVATIVE_HPP
#define DERIVATIVE_HPP

/*************************
 * CarpetX includes part *
 *************************/
#include <loop.hxx>
#include <vect.hxx>

namespace KleinGordonX {

using Arith::vect;
using Loop::PointDesc;
  
/* Stencil naming convention:
 * 1) Bias: Does the stencil has more points to the left (l), right (r) or is
 * centered?
 *
 * 2) Points: How many points does this stencil have?
 */
enum class Stencil { c3, c5, c7, c9 };

template <Stencil stencil> class Derivative {
private:
  const Loop::GF3D2<const CCTK_REAL> &gf;

public:
  Derivative(const Loop::GF3D2<const CCTK_REAL> &gridFunction)
      : gf(gridFunction) {}

  const Arith::vect<CCTK_REAL, 3> grad(const Loop::PointDesc &) const;
  const Arith::vect<CCTK_REAL, 6> lapl(const Loop::PointDesc &) const;
};

#include "c3.hpp"
  
} // namespace KleinGordonX

#endif // DERIVATIVE_HPP
