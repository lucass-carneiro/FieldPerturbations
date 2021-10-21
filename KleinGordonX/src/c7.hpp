#ifndef C7_HPP
#define C7_HPP

template <> class Derivative<Stencil::c7> {
private:
  const Loop::GF3D2<const CCTK_REAL> &gf;

public:
  Derivative(const Loop::GF3D2<const CCTK_REAL> &gridFunction) : gf(gridFunction) {}

  const Arith::vect<CCTK_REAL, 3> grad(const Loop::PointDesc &p) const {
    const Arith::vect<CCTK_REAL, 3> vals
        = {(-gf(-3 * p.DI[0] + p.I) + 9 * gf(-2 * p.DI[0] + p.I)
            - 9 * (5 * gf(-p.DI[0] + p.I) - 5 * gf(p.DI[0] + p.I) + gf(2 * p.DI[0] + p.I))
            + gf(3 * p.DI[0] + p.I))
               * (1.0 / (60 * p.dx)),
           (-gf(-3 * p.DI[1] + p.I) + 9 * gf(-2 * p.DI[1] + p.I)
            - 9 * (5 * gf(-p.DI[1] + p.I) - 5 * gf(p.DI[1] + p.I) + gf(2 * p.DI[1] + p.I))
            + gf(3 * p.DI[1] + p.I))
               * (1.0 / (60 * p.dy)),
           (-gf(-3 * p.DI[2] + p.I) + 9 * gf(-2 * p.DI[2] + p.I)
            - 9 * (5 * gf(-p.DI[2] + p.I) - 5 * gf(p.DI[2] + p.I) + gf(2 * p.DI[2] + p.I))
            + gf(3 * p.DI[2] + p.I))
               * (1.0 / (60 * p.dz))};
    return vals;
  }

  const Arith::vect<CCTK_REAL, 6> hess(const Loop::PointDesc &p) const {
    const Arith::vect<CCTK_REAL, 6> vals
        = {(-490 * gf(p.I) + 2 * gf(-3 * p.DI[0] + p.I) - 27 * gf(-2 * p.DI[0] + p.I)
            + 270 * gf(-p.DI[0] + p.I) + 270 * gf(p.DI[0] + p.I) - 27 * gf(2 * p.DI[0] + p.I)
            + 2 * gf(3 * p.DI[0] + p.I))
               * (1.0 / (180 * power<2>(p.dx))),
           (-490 * gf(p.I) + 2 * gf(-3 * p.DI[1] + p.I) - 27 * gf(-2 * p.DI[1] + p.I)
            + 270 * gf(-p.DI[1] + p.I) + 270 * gf(p.DI[1] + p.I) - 27 * gf(2 * p.DI[1] + p.I)
            + 2 * gf(3 * p.DI[1] + p.I))
               * (1.0 / (180 * power<2>(p.dy))),
           (-490 * gf(p.I) + 2 * gf(-3 * p.DI[2] + p.I) - 27 * gf(-2 * p.DI[2] + p.I)
            + 270 * gf(-p.DI[2] + p.I) + 270 * gf(p.DI[2] + p.I) - 27 * gf(2 * p.DI[2] + p.I)
            + 2 * gf(3 * p.DI[2] + p.I))
               * (1.0 / (180 * power<2>(p.dz))),
           (gf(-3 * p.DI[0] - 3 * p.DI[1] + p.I) - 9 * gf(-2 * p.DI[0] - 3 * p.DI[1] + p.I)
            + 45 * gf(-p.DI[0] - 3 * p.DI[1] + p.I) - 45 * gf(p.DI[0] - 3 * p.DI[1] + p.I)
            + 9 * gf(2 * p.DI[0] - 3 * p.DI[1] + p.I) - gf(3 * p.DI[0] - 3 * p.DI[1] + p.I)
            - 9 * gf(-3 * p.DI[0] - 2 * p.DI[1] + p.I) + 81 * gf(-2 * p.DI[0] - 2 * p.DI[1] + p.I)
            - 405 * gf(-p.DI[0] - 2 * p.DI[1] + p.I) + 405 * gf(p.DI[0] - 2 * p.DI[1] + p.I)
            - 81 * gf(2 * p.DI[0] - 2 * p.DI[1] + p.I) + 45 * gf(-3 * p.DI[0] - p.DI[1] + p.I)
            - 405 * gf(-2 * p.DI[0] - p.DI[1] + p.I) + 2025 * gf(-p.DI[0] - p.DI[1] + p.I)
            - 2025 * gf(p.DI[0] - p.DI[1] + p.I) + 405 * gf(2 * p.DI[0] - p.DI[1] + p.I)
            - 45 * gf(-3 * p.DI[0] + p.DI[1] + p.I) + 405 * gf(-2 * p.DI[0] + p.DI[1] + p.I)
            - 2025 * gf(-p.DI[0] + p.DI[1] + p.I) + 2025 * gf(p.DI[0] + p.DI[1] + p.I)
            - 405 * gf(2 * p.DI[0] + p.DI[1] + p.I) + 9 * gf(-3 * p.DI[0] + 2 * p.DI[1] + p.I)
            - 81 * gf(-2 * p.DI[0] + 2 * p.DI[1] + p.I) + 405 * gf(-p.DI[0] + 2 * p.DI[1] + p.I)
            - 405 * gf(p.DI[0] + 2 * p.DI[1] + p.I) + 81 * gf(2 * p.DI[0] + 2 * p.DI[1] + p.I)
            + 9
                  * (gf(3 * p.DI[0] - 2 * p.DI[1] + p.I) - 5 * gf(3 * p.DI[0] - p.DI[1] + p.I)
                     + 5 * gf(3 * p.DI[0] + p.DI[1] + p.I) - gf(3 * p.DI[0] + 2 * p.DI[1] + p.I))
            - gf(-3 * p.DI[0] + 3 * p.DI[1] + p.I) + 9 * gf(-2 * p.DI[0] + 3 * p.DI[1] + p.I)
            - 45 * gf(-p.DI[0] + 3 * p.DI[1] + p.I) + 45 * gf(p.DI[0] + 3 * p.DI[1] + p.I)
            - 9 * gf(2 * p.DI[0] + 3 * p.DI[1] + p.I) + gf(3 * p.DI[0] + 3 * p.DI[1] + p.I))
               * (1.0 / (3600 * p.dx * p.dy)),
           (gf(-3 * p.DI[0] - 3 * p.DI[2] + p.I) - 9 * gf(-2 * p.DI[0] - 3 * p.DI[2] + p.I)
            + 45 * gf(-p.DI[0] - 3 * p.DI[2] + p.I) - 45 * gf(p.DI[0] - 3 * p.DI[2] + p.I)
            + 9 * gf(2 * p.DI[0] - 3 * p.DI[2] + p.I) - gf(3 * p.DI[0] - 3 * p.DI[2] + p.I)
            - 9 * gf(-3 * p.DI[0] - 2 * p.DI[2] + p.I) + 81 * gf(-2 * p.DI[0] - 2 * p.DI[2] + p.I)
            - 405 * gf(-p.DI[0] - 2 * p.DI[2] + p.I) + 405 * gf(p.DI[0] - 2 * p.DI[2] + p.I)
            - 81 * gf(2 * p.DI[0] - 2 * p.DI[2] + p.I) + 45 * gf(-3 * p.DI[0] - p.DI[2] + p.I)
            - 405 * gf(-2 * p.DI[0] - p.DI[2] + p.I) + 2025 * gf(-p.DI[0] - p.DI[2] + p.I)
            - 2025 * gf(p.DI[0] - p.DI[2] + p.I) + 405 * gf(2 * p.DI[0] - p.DI[2] + p.I)
            - 45 * gf(-3 * p.DI[0] + p.DI[2] + p.I) + 405 * gf(-2 * p.DI[0] + p.DI[2] + p.I)
            - 2025 * gf(-p.DI[0] + p.DI[2] + p.I) + 2025 * gf(p.DI[0] + p.DI[2] + p.I)
            - 405 * gf(2 * p.DI[0] + p.DI[2] + p.I) + 9 * gf(-3 * p.DI[0] + 2 * p.DI[2] + p.I)
            - 81 * gf(-2 * p.DI[0] + 2 * p.DI[2] + p.I) + 405 * gf(-p.DI[0] + 2 * p.DI[2] + p.I)
            - 405 * gf(p.DI[0] + 2 * p.DI[2] + p.I) + 81 * gf(2 * p.DI[0] + 2 * p.DI[2] + p.I)
            + 9
                  * (gf(3 * p.DI[0] - 2 * p.DI[2] + p.I) - 5 * gf(3 * p.DI[0] - p.DI[2] + p.I)
                     + 5 * gf(3 * p.DI[0] + p.DI[2] + p.I) - gf(3 * p.DI[0] + 2 * p.DI[2] + p.I))
            - gf(-3 * p.DI[0] + 3 * p.DI[2] + p.I) + 9 * gf(-2 * p.DI[0] + 3 * p.DI[2] + p.I)
            - 45 * gf(-p.DI[0] + 3 * p.DI[2] + p.I) + 45 * gf(p.DI[0] + 3 * p.DI[2] + p.I)
            - 9 * gf(2 * p.DI[0] + 3 * p.DI[2] + p.I) + gf(3 * p.DI[0] + 3 * p.DI[2] + p.I))
               * (1.0 / (3600 * p.dx * p.dz)),
           (gf(-3 * p.DI[1] - 3 * p.DI[2] + p.I) - 9 * gf(-2 * p.DI[1] - 3 * p.DI[2] + p.I)
            + 45 * gf(-p.DI[1] - 3 * p.DI[2] + p.I) - 45 * gf(p.DI[1] - 3 * p.DI[2] + p.I)
            + 9 * gf(2 * p.DI[1] - 3 * p.DI[2] + p.I) - gf(3 * p.DI[1] - 3 * p.DI[2] + p.I)
            - 9 * gf(-3 * p.DI[1] - 2 * p.DI[2] + p.I) + 81 * gf(-2 * p.DI[1] - 2 * p.DI[2] + p.I)
            - 405 * gf(-p.DI[1] - 2 * p.DI[2] + p.I) + 405 * gf(p.DI[1] - 2 * p.DI[2] + p.I)
            - 81 * gf(2 * p.DI[1] - 2 * p.DI[2] + p.I) + 45 * gf(-3 * p.DI[1] - p.DI[2] + p.I)
            - 405 * gf(-2 * p.DI[1] - p.DI[2] + p.I) + 2025 * gf(-p.DI[1] - p.DI[2] + p.I)
            - 2025 * gf(p.DI[1] - p.DI[2] + p.I) + 405 * gf(2 * p.DI[1] - p.DI[2] + p.I)
            - 45 * gf(-3 * p.DI[1] + p.DI[2] + p.I) + 405 * gf(-2 * p.DI[1] + p.DI[2] + p.I)
            - 2025 * gf(-p.DI[1] + p.DI[2] + p.I) + 2025 * gf(p.DI[1] + p.DI[2] + p.I)
            - 405 * gf(2 * p.DI[1] + p.DI[2] + p.I) + 9 * gf(-3 * p.DI[1] + 2 * p.DI[2] + p.I)
            - 81 * gf(-2 * p.DI[1] + 2 * p.DI[2] + p.I) + 405 * gf(-p.DI[1] + 2 * p.DI[2] + p.I)
            - 405 * gf(p.DI[1] + 2 * p.DI[2] + p.I) + 81 * gf(2 * p.DI[1] + 2 * p.DI[2] + p.I)
            + 9
                  * (gf(3 * p.DI[1] - 2 * p.DI[2] + p.I) - 5 * gf(3 * p.DI[1] - p.DI[2] + p.I)
                     + 5 * gf(3 * p.DI[1] + p.DI[2] + p.I) - gf(3 * p.DI[1] + 2 * p.DI[2] + p.I))
            - gf(-3 * p.DI[1] + 3 * p.DI[2] + p.I) + 9 * gf(-2 * p.DI[1] + 3 * p.DI[2] + p.I)
            - 45 * gf(-p.DI[1] + 3 * p.DI[2] + p.I) + 45 * gf(p.DI[1] + 3 * p.DI[2] + p.I)
            - 9 * gf(2 * p.DI[1] + 3 * p.DI[2] + p.I) + gf(3 * p.DI[1] + 3 * p.DI[2] + p.I))
               * (1.0 / (3600 * p.dy * p.dz))};
    return vals;
  }
};

#endif // C7_HPP
