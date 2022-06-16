#ifndef FC_KLEIN_GORDON_DERIVATIVES
#define FC_KLEIN_GORDON_DERIVATIVES

//clang-format off
#include <cctk.h>
//clang-format on

#include <array>
#include <tuple>
#include <vector>

namespace fckg {

enum class derivative_direction : CCTK_INT { x = 0, y = 1, z = 2 };

template <derivative_direction dir> class sbp_coefficients {
public:
  sbp_coefficients(const cGH *restrict const cctkGH, CCTK_INT n)
      : nsize(n), imin(nsize), imax(nsize), q(nsize * nsize) {
    compute(cctkGH);
  }

  auto get_n() const -> CCTK_INT { return nsize; }
  auto get_imin() const -> const std::vector<CCTK_INT> & { return imin; }
  auto get_imax() const -> const std::vector<CCTK_INT> & { return imax; }
  auto get_q() const -> const std::vector<CCTK_REAL> & { return q; }

  template <typename grid_function> auto local_d(const grid_function &f, const CCTK_INT i,
                                                 const CCTK_INT j, const CCTK_INT k,
                                                 const CCTK_REAL h) const -> CCTK_REAL {

    CCTK_REAL dval = 0.0;

    if constexpr (dir == derivative_direction::x) {
      for (CCTK_INT a = imin[i]; a <= imax[i]; a++) {
        dval += q[a + i * nsize] * f(a, j, k);
      }
    } else if constexpr (dir == derivative_direction::y) {
      for (CCTK_INT a = imin[j]; a <= imax[j]; a++) {
        dval += q[a + j * nsize] * f(i, a, k);
      }
    } else if constexpr (dir == derivative_direction::z) {
      for (CCTK_INT a = imin[k]; a <= imax[k]; a++) {
        dval += q[a + k * nsize] * f(i, j, a);
      }
    }

    return (1.0 / h) * dval;
  }

private:
  const CCTK_INT nsize;
  std::vector<CCTK_INT> imin;
  std::vector<CCTK_INT> imax;
  std::vector<CCTK_REAL> q;

  void compute(const cGH *restrict const cctkGH) {
    Diff_coeff(cctkGH, static_cast<CCTK_INT>(dir), nsize, imin.data(), imax.data(), q.data(), -1);

    // Fortran -> C index conversion.
    for (auto &i : imin) {
      i -= 1;
    }

    for (auto &i : imax) {
      i -= 1;
    }
  }
};

template <derivative_direction dir, typename grid_function>
auto global_d(const sbp_coefficients<derivative_direction::x> &cx,
              const sbp_coefficients<derivative_direction::y> &cy,
              const sbp_coefficients<derivative_direction::z> &cz, const grid_function &f,
              const std::array<std::array<CCTK_REAL, 3>, 3> &J, const CCTK_INT i, const CCTK_INT j,
              const CCTK_INT k, const CCTK_REAL hx, const CCTK_REAL hy, const CCTK_REAL hz) {

  if constexpr (dir == derivative_direction::x) {
    return J[0][0] * cx.local_d(f, i, j, k, hx) + J[1][0] * cy.local_d(f, i, j, k, hy)
           + J[2][0] * cz.local_d(f, i, j, k, hz);
  } else if constexpr (dir == derivative_direction::y) {
    return J[0][1] * cx.local_d(f, i, j, k, hx) + J[1][1] * cy.local_d(f, i, j, k, hy)
           + J[2][1] * cz.local_d(f, i, j, k, hz);
  } else if constexpr (dir == derivative_direction::z) {
    return J[0][2] * cx.local_d(f, i, j, k, hx) + J[1][2] * cy.local_d(f, i, j, k, hy)
           + J[2][2] * cz.local_d(f, i, j, k, hz);
  }
}

} // namespace fckg

#endif // FC_KLEIN_GORDON_DERIVATIVES