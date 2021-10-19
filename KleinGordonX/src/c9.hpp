#ifndef C9_HPP
#define C9_HPP

template <> class Derivative<Stencil::c9> {
private:
  const Loop::GF3D2<const CCTK_REAL> &gf;

public:
  Derivative(const Loop::GF3D2<const CCTK_REAL> &gridFunction) : gf(gridFunction) {}

  const Arith::vect<CCTK_REAL, 3> grad(const Loop::PointDesc &p) const {
    const Arith::vect<CCTK_REAL, 3> vals
        = {(3 * gf(-4 * p.DI[0] + p.I) - 32 * gf(-3 * p.DI[0] + p.I)
            + 168
                  * (gf(-2 * p.DI[0] + p.I) - 4 * gf(-p.DI[0] + p.I) + 4 * gf(p.DI[0] + p.I)
                     - gf(2 * p.DI[0] + p.I))
            + 32 * gf(3 * p.DI[0] + p.I) - 3 * gf(4 * p.DI[0] + p.I))
               * (1.0 / (840 * p.dx)),
           (3 * gf(-4 * p.DI[1] + p.I) - 32 * gf(-3 * p.DI[1] + p.I)
            + 168
                  * (gf(-2 * p.DI[1] + p.I) - 4 * gf(-p.DI[1] + p.I) + 4 * gf(p.DI[1] + p.I)
                     - gf(2 * p.DI[1] + p.I))
            + 32 * gf(3 * p.DI[1] + p.I) - 3 * gf(4 * p.DI[1] + p.I))
               * (1.0 / (840 * p.dy)),
           (3 * gf(-4 * p.DI[2] + p.I) - 32 * gf(-3 * p.DI[2] + p.I)
            + 168
                  * (gf(-2 * p.DI[2] + p.I) - 4 * gf(-p.DI[2] + p.I) + 4 * gf(p.DI[2] + p.I)
                     - gf(2 * p.DI[2] + p.I))
            + 32 * gf(3 * p.DI[2] + p.I) - 3 * gf(4 * p.DI[2] + p.I))
               * (1.0 / (840 * p.dz))};
    return vals;
  }

  const Arith::vect<CCTK_REAL, 6> hess(const Loop::PointDesc &p) const {
    const Arith::vect<CCTK_REAL, 6> vals = {
        (-14350 * gf(p.I) - 9 * gf(-4 * p.DI[0] + p.I) + 128 * gf(-3 * p.DI[0] + p.I)
         - 1008 * gf(-2 * p.DI[0] + p.I) + 8064 * gf(-p.DI[0] + p.I) + 8064 * gf(p.DI[0] + p.I)
         - 1008 * gf(2 * p.DI[0] + p.I) + 128 * gf(3 * p.DI[0] + p.I) - 9 * gf(4 * p.DI[0] + p.I))
            * (1.0 / (5040 * power<2>(p.dx))),
        (-14350 * gf(p.I) - 9 * gf(-4 * p.DI[1] + p.I) + 128 * gf(-3 * p.DI[1] + p.I)
         - 1008 * gf(-2 * p.DI[1] + p.I) + 8064 * gf(-p.DI[1] + p.I) + 8064 * gf(p.DI[1] + p.I)
         - 1008 * gf(2 * p.DI[1] + p.I) + 128 * gf(3 * p.DI[1] + p.I) - 9 * gf(4 * p.DI[1] + p.I))
            * (1.0 / (5040 * power<2>(p.dy))),
        (-14350 * gf(p.I) - 9 * gf(-4 * p.DI[2] + p.I) + 128 * gf(-3 * p.DI[2] + p.I)
         - 1008 * gf(-2 * p.DI[2] + p.I) + 8064 * gf(-p.DI[2] + p.I) + 8064 * gf(p.DI[2] + p.I)
         - 1008 * gf(2 * p.DI[2] + p.I) + 128 * gf(3 * p.DI[2] + p.I) - 9 * gf(4 * p.DI[2] + p.I))
            * (1.0 / (5040 * power<2>(p.dz))),
        (9 * gf(-4 * p.DI[0] - 4 * p.DI[1] + p.I) - 96 * gf(-3 * p.DI[0] - 4 * p.DI[1] + p.I)
         + 504 * gf(-2 * p.DI[0] - 4 * p.DI[1] + p.I) - 2016 * gf(-p.DI[0] - 4 * p.DI[1] + p.I)
         + 2016 * gf(p.DI[0] - 4 * p.DI[1] + p.I) - 504 * gf(2 * p.DI[0] - 4 * p.DI[1] + p.I)
         + 96 * gf(3 * p.DI[0] - 4 * p.DI[1] + p.I) - 9 * gf(4 * p.DI[0] - 4 * p.DI[1] + p.I)
         - 96 * gf(-4 * p.DI[0] - 3 * p.DI[1] + p.I) + 1024 * gf(-3 * p.DI[0] - 3 * p.DI[1] + p.I)
         - 5376 * gf(-2 * p.DI[0] - 3 * p.DI[1] + p.I) + 21504 * gf(-p.DI[0] - 3 * p.DI[1] + p.I)
         - 21504 * gf(p.DI[0] - 3 * p.DI[1] + p.I) + 5376 * gf(2 * p.DI[0] - 3 * p.DI[1] + p.I)
         - 1024 * gf(3 * p.DI[0] - 3 * p.DI[1] + p.I) + 96 * gf(4 * p.DI[0] - 3 * p.DI[1] + p.I)
         + 504 * gf(-4 * p.DI[0] - 2 * p.DI[1] + p.I) - 5376 * gf(-3 * p.DI[0] - 2 * p.DI[1] + p.I)
         + 28224 * gf(-2 * p.DI[0] - 2 * p.DI[1] + p.I) - 112896 * gf(-p.DI[0] - 2 * p.DI[1] + p.I)
         + 112896 * gf(p.DI[0] - 2 * p.DI[1] + p.I) - 28224 * gf(2 * p.DI[0] - 2 * p.DI[1] + p.I)
         + 5376 * gf(3 * p.DI[0] - 2 * p.DI[1] + p.I) - 2016 * gf(-4 * p.DI[0] - p.DI[1] + p.I)
         + 21504 * gf(-3 * p.DI[0] - p.DI[1] + p.I) - 112896 * gf(-2 * p.DI[0] - p.DI[1] + p.I)
         + 451584 * gf(-p.DI[0] - p.DI[1] + p.I) - 451584 * gf(p.DI[0] - p.DI[1] + p.I)
         + 112896 * gf(2 * p.DI[0] - p.DI[1] + p.I) - 21504 * gf(3 * p.DI[0] - p.DI[1] + p.I)
         + 2016 * gf(-4 * p.DI[0] + p.DI[1] + p.I) - 21504 * gf(-3 * p.DI[0] + p.DI[1] + p.I)
         + 112896 * gf(-2 * p.DI[0] + p.DI[1] + p.I) - 451584 * gf(-p.DI[0] + p.DI[1] + p.I)
         + 451584 * gf(p.DI[0] + p.DI[1] + p.I) - 112896 * gf(2 * p.DI[0] + p.DI[1] + p.I)
         + 21504 * gf(3 * p.DI[0] + p.DI[1] + p.I) - 504 * gf(-4 * p.DI[0] + 2 * p.DI[1] + p.I)
         + 5376 * gf(-3 * p.DI[0] + 2 * p.DI[1] + p.I)
         - 28224 * gf(-2 * p.DI[0] + 2 * p.DI[1] + p.I) + 112896 * gf(-p.DI[0] + 2 * p.DI[1] + p.I)
         - 112896 * gf(p.DI[0] + 2 * p.DI[1] + p.I) + 28224 * gf(2 * p.DI[0] + 2 * p.DI[1] + p.I)
         - 5376 * gf(3 * p.DI[0] + 2 * p.DI[1] + p.I)
         - 504 * (gf(4 * p.DI[0] - 2 * p.DI[1] + p.I) - 4 * gf(4 * p.DI[0] - p.DI[1] + p.I) + 4 * gf(4 * p.DI[0] + p.DI[1] + p.I) - gf(4 * p.DI[0] + 2 * p.DI[1] + p.I))
         + 96 * gf(-4 * p.DI[0] + 3 * p.DI[1] + p.I) - 1024 * gf(-3 * p.DI[0] + 3 * p.DI[1] + p.I)
         + 5376 * gf(-2 * p.DI[0] + 3 * p.DI[1] + p.I) - 21504 * gf(-p.DI[0] + 3 * p.DI[1] + p.I)
         + 21504 * gf(p.DI[0] + 3 * p.DI[1] + p.I) - 5376 * gf(2 * p.DI[0] + 3 * p.DI[1] + p.I)
         + 1024 * gf(3 * p.DI[0] + 3 * p.DI[1] + p.I) - 96 * gf(4 * p.DI[0] + 3 * p.DI[1] + p.I)
         - 9 * gf(-4 * p.DI[0] + 4 * p.DI[1] + p.I) + 96 * gf(-3 * p.DI[0] + 4 * p.DI[1] + p.I)
         - 504 * gf(-2 * p.DI[0] + 4 * p.DI[1] + p.I) + 2016 * gf(-p.DI[0] + 4 * p.DI[1] + p.I)
         - 2016 * gf(p.DI[0] + 4 * p.DI[1] + p.I) + 504 * gf(2 * p.DI[0] + 4 * p.DI[1] + p.I)
         - 96 * gf(3 * p.DI[0] + 4 * p.DI[1] + p.I) + 9 * gf(4 * p.DI[0] + 4 * p.DI[1] + p.I))
            * (1.0 / (705600 * p.dx * p.dy)),
        (9 * gf(-4 * p.DI[0] - 4 * p.DI[2] + p.I) - 96 * gf(-3 * p.DI[0] - 4 * p.DI[2] + p.I)
         + 504 * gf(-2 * p.DI[0] - 4 * p.DI[2] + p.I) - 2016 * gf(-p.DI[0] - 4 * p.DI[2] + p.I)
         + 2016 * gf(p.DI[0] - 4 * p.DI[2] + p.I) - 504 * gf(2 * p.DI[0] - 4 * p.DI[2] + p.I)
         + 96 * gf(3 * p.DI[0] - 4 * p.DI[2] + p.I) - 9 * gf(4 * p.DI[0] - 4 * p.DI[2] + p.I)
         - 96 * gf(-4 * p.DI[0] - 3 * p.DI[2] + p.I) + 1024 * gf(-3 * p.DI[0] - 3 * p.DI[2] + p.I)
         - 5376 * gf(-2 * p.DI[0] - 3 * p.DI[2] + p.I) + 21504 * gf(-p.DI[0] - 3 * p.DI[2] + p.I)
         - 21504 * gf(p.DI[0] - 3 * p.DI[2] + p.I) + 5376 * gf(2 * p.DI[0] - 3 * p.DI[2] + p.I)
         - 1024 * gf(3 * p.DI[0] - 3 * p.DI[2] + p.I) + 96 * gf(4 * p.DI[0] - 3 * p.DI[2] + p.I)
         + 504 * gf(-4 * p.DI[0] - 2 * p.DI[2] + p.I) - 5376 * gf(-3 * p.DI[0] - 2 * p.DI[2] + p.I)
         + 28224 * gf(-2 * p.DI[0] - 2 * p.DI[2] + p.I) - 112896 * gf(-p.DI[0] - 2 * p.DI[2] + p.I)
         + 112896 * gf(p.DI[0] - 2 * p.DI[2] + p.I) - 28224 * gf(2 * p.DI[0] - 2 * p.DI[2] + p.I)
         + 5376 * gf(3 * p.DI[0] - 2 * p.DI[2] + p.I) - 2016 * gf(-4 * p.DI[0] - p.DI[2] + p.I)
         + 21504 * gf(-3 * p.DI[0] - p.DI[2] + p.I) - 112896 * gf(-2 * p.DI[0] - p.DI[2] + p.I)
         + 451584 * gf(-p.DI[0] - p.DI[2] + p.I) - 451584 * gf(p.DI[0] - p.DI[2] + p.I)
         + 112896 * gf(2 * p.DI[0] - p.DI[2] + p.I) - 21504 * gf(3 * p.DI[0] - p.DI[2] + p.I)
         + 2016 * gf(-4 * p.DI[0] + p.DI[2] + p.I) - 21504 * gf(-3 * p.DI[0] + p.DI[2] + p.I)
         + 112896 * gf(-2 * p.DI[0] + p.DI[2] + p.I) - 451584 * gf(-p.DI[0] + p.DI[2] + p.I)
         + 451584 * gf(p.DI[0] + p.DI[2] + p.I) - 112896 * gf(2 * p.DI[0] + p.DI[2] + p.I)
         + 21504 * gf(3 * p.DI[0] + p.DI[2] + p.I) - 504 * gf(-4 * p.DI[0] + 2 * p.DI[2] + p.I)
         + 5376 * gf(-3 * p.DI[0] + 2 * p.DI[2] + p.I)
         - 28224 * gf(-2 * p.DI[0] + 2 * p.DI[2] + p.I) + 112896 * gf(-p.DI[0] + 2 * p.DI[2] + p.I)
         - 112896 * gf(p.DI[0] + 2 * p.DI[2] + p.I) + 28224 * gf(2 * p.DI[0] + 2 * p.DI[2] + p.I)
         - 5376 * gf(3 * p.DI[0] + 2 * p.DI[2] + p.I)
         - 504 * (gf(4 * p.DI[0] - 2 * p.DI[2] + p.I) - 4 * gf(4 * p.DI[0] - p.DI[2] + p.I) + 4 * gf(4 * p.DI[0] + p.DI[2] + p.I) - gf(4 * p.DI[0] + 2 * p.DI[2] + p.I))
         + 96 * gf(-4 * p.DI[0] + 3 * p.DI[2] + p.I) - 1024 * gf(-3 * p.DI[0] + 3 * p.DI[2] + p.I)
         + 5376 * gf(-2 * p.DI[0] + 3 * p.DI[2] + p.I) - 21504 * gf(-p.DI[0] + 3 * p.DI[2] + p.I)
         + 21504 * gf(p.DI[0] + 3 * p.DI[2] + p.I) - 5376 * gf(2 * p.DI[0] + 3 * p.DI[2] + p.I)
         + 1024 * gf(3 * p.DI[0] + 3 * p.DI[2] + p.I) - 96 * gf(4 * p.DI[0] + 3 * p.DI[2] + p.I)
         - 9 * gf(-4 * p.DI[0] + 4 * p.DI[2] + p.I) + 96 * gf(-3 * p.DI[0] + 4 * p.DI[2] + p.I)
         - 504 * gf(-2 * p.DI[0] + 4 * p.DI[2] + p.I) + 2016 * gf(-p.DI[0] + 4 * p.DI[2] + p.I)
         - 2016 * gf(p.DI[0] + 4 * p.DI[2] + p.I) + 504 * gf(2 * p.DI[0] + 4 * p.DI[2] + p.I)
         - 96 * gf(3 * p.DI[0] + 4 * p.DI[2] + p.I) + 9 * gf(4 * p.DI[0] + 4 * p.DI[2] + p.I))
            * (1.0 / (705600 * p.dx * p.dz)),
        (
            9 * gf(-4 * p.DI[1] - 4 * p.DI[2] + p.I) - 96 * gf(-3 * p.DI[1] - 4 * p.DI[2] + p.I)
            + 504 * gf(-2 * p.DI[1] - 4 * p.DI[2] + p.I) - 2016 * gf(-p.DI[1] - 4 * p.DI[2] + p.I)
            + 2016 * gf(p.DI[1] - 4 * p.DI[2] + p.I) - 504 * gf(2 * p.DI[1] - 4 * p.DI[2] + p.I)
            + 96 * gf(3 * p.DI[1] - 4 * p.DI[2] + p.I) - 9 * gf(4 * p.DI[1] - 4 * p.DI[2] + p.I)
            - 96 * gf(-4 * p.DI[1] - 3 * p.DI[2] + p.I)
            + 1024 * gf(-3 * p.DI[1] - 3 * p.DI[2] + p.I)
            - 5376 * gf(-2 * p.DI[1] - 3 * p.DI[2] + p.I) + 21504 * gf(-p.DI[1] - 3 * p.DI[2] + p.I)
            - 21504 * gf(p.DI[1] - 3 * p.DI[2] + p.I) + 5376 * gf(2 * p.DI[1] - 3 * p.DI[2] + p.I)
            - 1024 * gf(3 * p.DI[1] - 3 * p.DI[2] + p.I) + 96 * gf(4 * p.DI[1] - 3 * p.DI[2] + p.I)
            + 504 * gf(-4 * p.DI[1] - 2 * p.DI[2] + p.I)
            - 5376 * gf(-3 * p.DI[1] - 2 * p.DI[2] + p.I)
            + 28224 * gf(-2 * p.DI[1] - 2 * p.DI[2] + p.I)
            - 112896 * gf(-p.DI[1] - 2 * p.DI[2] + p.I) + 112896 * gf(p.DI[1] - 2 * p.DI[2] + p.I)
            - 28224 * gf(2 * p.DI[1] - 2 * p.DI[2] + p.I)
            + 5376 * gf(3 * p.DI[1] - 2 * p.DI[2] + p.I) - 2016 * gf(-4 * p.DI[1] - p.DI[2] + p.I)
            + 21504 * gf(-3 * p.DI[1] - p.DI[2] + p.I) - 112896 * gf(-2 * p.DI[1] - p.DI[2] + p.I)
            + 451584 * gf(-p.DI[1] - p.DI[2] + p.I) - 451584 * gf(p.DI[1] - p.DI[2] + p.I)
            + 112896 * gf(2 * p.DI[1] - p.DI[2] + p.I) - 21504 * gf(3 * p.DI[1] - p.DI[2] + p.I)
            + 2016 * gf(-4 * p.DI[1] + p.DI[2] + p.I) - 21504 * gf(-3 * p.DI[1] + p.DI[2] + p.I)
            + 112896 * gf(-2 * p.DI[1] + p.DI[2] + p.I) - 451584 * gf(-p.DI[1] + p.DI[2] + p.I)
            + 451584 * gf(p.DI[1] + p.DI[2] + p.I) - 112896 * gf(2 * p.DI[1] + p.DI[2] + p.I)
            + 21504 * gf(3 * p.DI[1] + p.DI[2] + p.I) - 504 * gf(-4 * p.DI[1] + 2 * p.DI[2] + p.I)
            + 5376 * gf(-3 * p.DI[1] + 2 * p.DI[2] + p.I)
            - 28224 * gf(-2 * p.DI[1] + 2 * p.DI[2] + p.I)
            + 112896 * gf(-p.DI[1] + 2 * p.DI[2] + p.I) - 112896 * gf(p.DI[1] + 2 * p.DI[2] + p.I)
            + 28224 * gf(2 * p.DI[1] + 2 * p.DI[2] + p.I)
            - 5376 * gf(3 * p.DI[1] + 2 * p.DI[2] + p.I)
            - 504 * (gf(4 * p.DI[1] - 2 * p.DI[2] + p.I) - 4 * gf(4 * p.DI[1] - p.DI[2] + p.I) + 4 * gf(4 * p.DI[1] + p.DI[2] + p.I) - gf(4 * p.DI[1] + 2 * p.DI[2] + p.I))
            + 96 * gf(-4 * p.DI[1] + 3 * p.DI[2] + p.I)
            - 1024 * gf(-3 * p.DI[1] + 3 * p.DI[2] + p.I)
            + 5376 * gf(-2 * p.DI[1] + 3 * p.DI[2] + p.I) - 21504 * gf(-p.DI[1] + 3 * p.DI[2] + p.I)
            + 21504 * gf(p.DI[1] + 3 * p.DI[2] + p.I) - 5376 * gf(2 * p.DI[1] + 3 * p.DI[2] + p.I)
            + 1024 * gf(3 * p.DI[1] + 3 * p.DI[2] + p.I) - 96 * gf(4 * p.DI[1] + 3 * p.DI[2] + p.I)
            - 9 * gf(-4 * p.DI[1] + 4 * p.DI[2] + p.I) + 96 * gf(-3 * p.DI[1] + 4 * p.DI[2] + p.I)
            - 504 * gf(-2 * p.DI[1] + 4 * p.DI[2] + p.I) + 2016 * gf(-p.DI[1] + 4 * p.DI[2] + p.I)
            - 2016 * gf(p.DI[1] + 4 * p.DI[2] + p.I) + 504 * gf(2 * p.DI[1] + 4 * p.DI[2] + p.I)
            - 96 * gf(3 * p.DI[1] + 4 * p.DI[2] + p.I) + 9 * gf(4 * p.DI[1] + 4 * p.DI[2] + p.I))
            * (1.0 / (705600 * p.dy * p.dz))};
    return vals;
  }
};

#endif // C9_HPP
