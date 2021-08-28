#ifndef C5_HPP
#define C5_HPP

template <> class Derivative<Stencil::c5> {
private:
  const Loop::GF3D2<const CCTK_REAL> &gf;

public:
  Derivative(const Loop::GF3D2<const CCTK_REAL> &gridFunction)
      : gf(gridFunction) {}

  const Arith::vect<CCTK_REAL, 3> grad(const Loop::PointDesc &p) const {
    const vect<CCTK_REAL, 3> vals = {
        (gf(-2 * p.DI[0] + p.I) - 8 * gf(-p.DI[0] + p.I) +
         8 * gf(p.DI[0] + p.I) - gf(2 * p.DI[0] + p.I)) *
            (1.0 / (12 * p.dx)),
        (gf(-2 * p.DI[1] + p.I) - 8 * gf(-p.DI[1] + p.I) +
         8 * gf(p.DI[1] + p.I) - gf(2 * p.DI[1] + p.I)) *
            (1.0 / (12 * p.dy)),
        (gf(-2 * p.DI[2] + p.I) - 8 * gf(-p.DI[2] + p.I) +
         8 * gf(p.DI[2] + p.I) - gf(2 * p.DI[2] + p.I)) *
            (1.0 / (12 * p.dz))};
    return vals;
  }

  const Arith::vect<CCTK_REAL, 6> hess(const Loop::PointDesc &p) const {
    const vect<CCTK_REAL, 6> vals = {
        (-30 * gf(p.I) - gf(-2 * p.DI[0] + p.I) + 16 * gf(-p.DI[0] + p.I) +
         16 * gf(p.DI[0] + p.I) - gf(2 * p.DI[0] + p.I)) *
            (1.0 / (12 * power<2>(p.dx))),
        (-30 * gf(p.I) - gf(-2 * p.DI[1] + p.I) + 16 * gf(-p.DI[1] + p.I) +
         16 * gf(p.DI[1] + p.I) - gf(2 * p.DI[1] + p.I)) *
            (1.0 / (12 * power<2>(p.dy))),
        (-30 * gf(p.I) - gf(-2 * p.DI[2] + p.I) + 16 * gf(-p.DI[2] + p.I) +
         16 * gf(p.DI[2] + p.I) - gf(2 * p.DI[2] + p.I)) *
            (1.0 / (12 * power<2>(p.dz))),
        (gf(-2 * p.DI[0] - 2 * p.DI[1] + p.I) -
         8 * gf(-p.DI[0] - 2 * p.DI[1] + p.I) +
         8 * gf(p.DI[0] - 2 * p.DI[1] + p.I) -
         gf(2 * p.DI[0] - 2 * p.DI[1] + p.I) -
         8 * gf(-2 * p.DI[0] - p.DI[1] + p.I) +
         64 * gf(-p.DI[0] - p.DI[1] + p.I) - 64 * gf(p.DI[0] - p.DI[1] + p.I) +
         8 * gf(2 * p.DI[0] - p.DI[1] + p.I) +
         8 * gf(-2 * p.DI[0] + p.DI[1] + p.I) -
         64 * gf(-p.DI[0] + p.DI[1] + p.I) + 64 * gf(p.DI[0] + p.DI[1] + p.I) -
         8 * gf(2 * p.DI[0] + p.DI[1] + p.I) -
         gf(-2 * p.DI[0] + 2 * p.DI[1] + p.I) +
         8 * gf(-p.DI[0] + 2 * p.DI[1] + p.I) -
         8 * gf(p.DI[0] + 2 * p.DI[1] + p.I) +
         gf(2 * p.DI[0] + 2 * p.DI[1] + p.I)) *
            (1.0 / (144 * p.dx * p.dy)),
        (gf(-2 * p.DI[0] - 2 * p.DI[2] + p.I) -
         8 * gf(-p.DI[0] - 2 * p.DI[2] + p.I) +
         8 * gf(p.DI[0] - 2 * p.DI[2] + p.I) -
         gf(2 * p.DI[0] - 2 * p.DI[2] + p.I) -
         8 * gf(-2 * p.DI[0] - p.DI[2] + p.I) +
         64 * gf(-p.DI[0] - p.DI[2] + p.I) - 64 * gf(p.DI[0] - p.DI[2] + p.I) +
         8 * gf(2 * p.DI[0] - p.DI[2] + p.I) +
         8 * gf(-2 * p.DI[0] + p.DI[2] + p.I) -
         64 * gf(-p.DI[0] + p.DI[2] + p.I) + 64 * gf(p.DI[0] + p.DI[2] + p.I) -
         8 * gf(2 * p.DI[0] + p.DI[2] + p.I) -
         gf(-2 * p.DI[0] + 2 * p.DI[2] + p.I) +
         8 * gf(-p.DI[0] + 2 * p.DI[2] + p.I) -
         8 * gf(p.DI[0] + 2 * p.DI[2] + p.I) +
         gf(2 * p.DI[0] + 2 * p.DI[2] + p.I)) *
            (1.0 / (144 * p.dx * p.dz)),
        (gf(-2 * p.DI[1] - 2 * p.DI[2] + p.I) -
         8 * gf(-p.DI[1] - 2 * p.DI[2] + p.I) +
         8 * gf(p.DI[1] - 2 * p.DI[2] + p.I) -
         gf(2 * p.DI[1] - 2 * p.DI[2] + p.I) -
         8 * gf(-2 * p.DI[1] - p.DI[2] + p.I) +
         64 * gf(-p.DI[1] - p.DI[2] + p.I) - 64 * gf(p.DI[1] - p.DI[2] + p.I) +
         8 * gf(2 * p.DI[1] - p.DI[2] + p.I) +
         8 * gf(-2 * p.DI[1] + p.DI[2] + p.I) -
         64 * gf(-p.DI[1] + p.DI[2] + p.I) + 64 * gf(p.DI[1] + p.DI[2] + p.I) -
         8 * gf(2 * p.DI[1] + p.DI[2] + p.I) -
         gf(-2 * p.DI[1] + 2 * p.DI[2] + p.I) +
         8 * gf(-p.DI[1] + 2 * p.DI[2] + p.I) -
         8 * gf(p.DI[1] + 2 * p.DI[2] + p.I) +
         gf(2 * p.DI[1] + 2 * p.DI[2] + p.I)) *
            (1.0 / (144 * p.dy * p.dz))};
    return vals;
  }
};

#endif // C5_HPP
