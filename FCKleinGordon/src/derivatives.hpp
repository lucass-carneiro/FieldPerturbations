#ifndef FC_KLEIN_GORDON_INITIAL_DERIVATIVES_HPP
#define FC_KLEIN_GORDON_INITIAL_DERIVATIVES_HPP

#include <cctk.h>
#include <cstddef>

namespace fckg {

struct deriv_data {
  CCTK_INT i;
  CCTK_INT j;
  CCTK_INT k;
  CCTK_REAL dx;
  CCTK_REAL dy;
  CCTK_REAL dz;
};

template <typename cctkgh_t>
static inline auto I(cctkgh_t cctkGH, CCTK_INT i, CCTK_INT j, CCTK_INT k) noexcept -> CCTK_INT {
  return CCTK_GFINDEX3D(cctkGH, i, j, k);
}

template <std::size_t order, typename cctkgh_t, typename cctk_gf_t>
static inline auto local_Dx(cctkgh_t cctkGH, const deriv_data &d, const cctk_gf_t &f) -> CCTK_REAL {
  if constexpr (order == 4) {
    const auto den{1.0 / (12.0 * d.dx)};
    const auto num{-f[I(cctkGH, d.i + 2, d.j, d.k)] + 8 * f[I(cctkGH, d.i + 1, d.j, d.k)]
                   - 8 * f[I(cctkGH, d.i - 1, d.j, d.k)] + f[I(cctkGH, d.i - 2, d.j, d.k)]};
    return den * num;
  }
}

template <std::size_t order, typename cctkgh_t, typename cctk_gf_t>
static inline auto local_Dy(cctkgh_t cctkGH, const deriv_data &d, const cctk_gf_t &f) -> CCTK_REAL {
  if constexpr (order == 4) {
    const auto den{1.0 / (12.0 * d.dy)};
    const auto num{-f[I(cctkGH, d.i, d.j + 2, d.k)] + 8 * f[I(cctkGH, d.i, d.j + 1, d.k)]
                   - 8 * f[I(cctkGH, d.i, d.j - 1, d.k)] + f[I(cctkGH, d.i, d.j - 2, d.k)]};
    return den * num;
  }
}

template <std::size_t order, typename cctkgh_t, typename cctk_gf_t>
static inline auto local_Dz(cctkgh_t cctkGH, const deriv_data &d, const cctk_gf_t &f) -> CCTK_REAL {
  if constexpr (order == 4) {
    const auto den{1.0 / (12.0 * d.dz)};
    const auto num{-f[I(cctkGH, d.i, d.j, d.k + 2)] + 8 * f[I(cctkGH, d.i, d.j, d.k + 1)]
                   - 8 * f[I(cctkGH, d.i, d.j, d.k - 1)] + f[I(cctkGH, d.i, d.j, d.k - 2)]};
    return den * num;
  }
}

template <std::size_t order = 4, typename cctkgh_t, typename cctk_gf_t>
static inline auto global_Dx(cctkgh_t cctkGH, const deriv_data &d, const cctk_gf_t &f,
                             CCTK_REAL J11, CCTK_REAL J21, CCTK_REAL J31) -> CCTK_REAL {
  return J11 * local_Dx<order>(cctkGH, d, f) + J21 * local_Dy<order>(cctkGH, d, f)
         + J31 * local_Dz<order>(cctkGH, d, f);
}

template <std::size_t order = 4, typename cctkgh_t, typename cctk_gf_t>
static inline auto global_Dy(cctkgh_t cctkGH, const deriv_data &d, const cctk_gf_t &f,
                             CCTK_REAL J12, CCTK_REAL J22, CCTK_REAL J32) -> CCTK_REAL {
  return J12 * local_Dx<order>(cctkGH, d, f) + J22 * local_Dy<order>(cctkGH, d, f)
         + J32 * local_Dz<order>(cctkGH, d, f);
}

template <std::size_t order = 4, typename cctkgh_t, typename cctk_gf_t>
static inline auto global_Dz(cctkgh_t cctkGH, const deriv_data &d, const cctk_gf_t &f,
                             CCTK_REAL J13, CCTK_REAL J23, CCTK_REAL J33) -> CCTK_REAL {
  return J13 * local_Dx<order>(cctkGH, d, f) + J23 * local_Dy<order>(cctkGH, d, f)
         + J33 * local_Dz<order>(cctkGH, d, f);
}

} // namespace fckg

#endif // FC_KLEIN_GORDON_INITIAL_DERIVATIVES_HPP
