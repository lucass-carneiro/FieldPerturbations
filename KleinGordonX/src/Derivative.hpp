#ifndef DERIVATIVE_HPP
#define DERIVATIVE_HPP

/*************************
 * CarpetX includes part *
 *************************/
#include <loop.hxx>
#include <vect.hxx>

/*************************
 * This thorn's includes *
 *************************/
#include "Arithmetic.hpp"

namespace KleinGordonX {

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
  const Arith::vect<CCTK_REAL, 6> hess(const Loop::PointDesc &) const;
};

#include "c3.hpp"
#include "c5.hpp"
#include "c7.hpp"
#include "c9.hpp"

} // namespace KleinGordonX

#endif // DERIVATIVE_HPP
