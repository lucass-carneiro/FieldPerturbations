#ifndef C3_HPP
#define C3_HPP

template <> class Derivative<Stencil::c3> {
private:
  const Loop::GF3D2<const CCTK_REAL> &gf;

public:
  Derivative(const Loop::GF3D2<const CCTK_REAL> &gridFunction)
      : gf(gridFunction) {}

  const Arith::vect<CCTK_REAL, 3> grad(const Loop::PointDesc &p) const {
    const vect<CCTK_REAL, 3> vals = {
        (-gf(-p.DI[0] + p.I) + gf(p.DI[0] + p.I)) * (1.0 / (2 * p.dx)),
        (-gf(-p.DI[1] + p.I) + gf(p.DI[1] + p.I)) * (1.0 / (2 * p.dy)),
        (-gf(-p.DI[2] + p.I) + gf(p.DI[2] + p.I)) * (1.0 / (2 * p.dz))};
    return vals;
  }

  const Arith::vect<CCTK_REAL, 6> hess(const Loop::PointDesc &p) const {
    const vect<CCTK_REAL, 6> vals = {
        (-2 * gf(p.I) + gf(-p.DI[0] + p.I) + gf(p.DI[0] + p.I)) *
            (1.0 / (power<2>(p.dx))),
        (-2 * gf(p.I) + gf(-p.DI[1] + p.I) + gf(p.DI[1] + p.I)) *
            (1.0 / (power<2>(p.dy))),
        (-2 * gf(p.I) + gf(-p.DI[2] + p.I) + gf(p.DI[2] + p.I)) *
            (1.0 / (power<2>(p.dz))),
        (gf(-p.DI[0] - p.DI[1] + p.I) - gf(p.DI[0] - p.DI[1] + p.I) -
         gf(-p.DI[0] + p.DI[1] + p.I) + gf(p.DI[0] + p.DI[1] + p.I)) *
            (1.0 / (4 * p.dx * p.dy)),
        (gf(-p.DI[0] - p.DI[2] + p.I) - gf(p.DI[0] - p.DI[2] + p.I) -
         gf(-p.DI[0] + p.DI[2] + p.I) + gf(p.DI[0] + p.DI[2] + p.I)) *
            (1.0 / (4 * p.dx * p.dz)),
        (gf(-p.DI[1] - p.DI[2] + p.I) - gf(p.DI[1] - p.DI[2] + p.I) -
         gf(-p.DI[1] + p.DI[2] + p.I) + gf(p.DI[1] + p.DI[2] + p.I)) *
            (1.0 / (4 * p.dy * p.dz))};
    return vals;
  }
};

#endif // C3_HPP
