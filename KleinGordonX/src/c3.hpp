#ifndef C3_HPP
#define C3_HPP

template <> class Derivative<Stencil::c3> {
private:
  const Loop::GF3D2<const CCTK_REAL> &gf;

public:
  Derivative(const Loop::GF3D2<const CCTK_REAL> &gridFunction)
      : gf(gridFunction) {}

  const Arith::vect<CCTK_REAL, 3> grad(const Loop::PointDesc &p) const {
    const vect<CCTK_REAL, 3> h = {p.dx, p.dy, p.dz};
    vect<CCTK_REAL, 3> vals = {0.0, 0.0, 0.0};

    for (size_t i = 0; i < 3; i++)
      vals[i] = (-gf(p.I - p.DI[i]) + gf(p.I + p.DI[i])) / (2.0 * h[i]);

    return vals;
  }

  const Arith::vect<CCTK_REAL, 6> hess(const Loop::PointDesc &p) const {
    const vect<CCTK_REAL, 6> h = {p.dx * p.dx, p.dy * p.dy, p.dz * p.dz,
                                  p.dx * p.dy, p.dx * p.dz, p.dy * p.dz};

    vect<CCTK_REAL, 6> vals = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    // Diagonal components
    for (size_t i = 0; i < 3; i++)
      vals[i] = (gf(p.I - p.DI[i]) - 2 * gf(p.I) + gf(p.I + p.DI[i])) / h[i];

    // Off diagonal components. TODO: Find a smarter way to do this
    // (use a loop?)
    vals[3] = (gf(p.I - p.DI[0] - p.DI[1]) - gf(p.I - p.DI[0] + p.DI[1]) -
               gf(p.I + p.DI[0] - p.DI[1]) + gf(p.I + p.DI[0] + p.DI[1])) /
              (4.0 * h[3]);

    vals[4] = (gf(p.I - p.DI[0] - p.DI[2]) - gf(p.I - p.DI[0] + p.DI[2]) -
               gf(p.I + p.DI[0] - p.DI[2]) + gf(p.I + p.DI[0] + p.DI[2])) /
              (4.0 * h[4]);

    vals[5] = (gf(p.I - p.DI[1] - p.DI[2]) - gf(p.I - p.DI[1] + p.DI[2]) -
               gf(p.I + p.DI[1] - p.DI[2]) + gf(p.I + p.DI[1] + p.DI[2])) /
              (4.0 * h[5]);

    return vals;
  }
};

#endif // C3_HPP
