#ifndef FC_KLEIN_GORDON_INITIAL_CONDITIONS_HPP
#define FC_KLEIN_GORDON_INITIAL_CONDITIONS_HPP

#include <cmath>
#include <cstdlib>
#include <limits>

namespace fckg {

template <typename fp_type> inline bool isapprox(fp_type x, fp_type y, fp_type atol = 0.0) {
  using std::abs;
  using std::max;
  using std::sqrt;

  const fp_type rtol = atol > 0.0 ? 0.0 : sqrt(std::numeric_limits<fp_type>::epsilon());
  return abs(x - y) <= max(atol, rtol * max(abs(x), abs(y)));
}

template <typename T> inline auto base_gaussian(T r, T sigma) -> T {
  using std::exp;
  return exp(-(r * r) / (2.0 * sigma * sigma));
}

template <typename T> inline auto d_base_gaussian_dr(T r, T sigma) -> T {
  using std::exp;
  return -base_gaussian(r, sigma) * r / (sigma * sigma);
}

template <typename T> inline auto exact_gaussian_solution(T t, T r, T sigma) -> T {
  if (isapprox(r, 0.0)) {
    return 2.0 * base_gaussian(t, sigma) * t / (sigma * sigma);
  } else {
    return (base_gaussian(r - t, sigma) - base_gaussian(r + t, sigma)) / r;
  }
}

template <typename T> inline auto d_exact_gaussian_solution_dt(T t, T r, T sigma) -> T {
  if (isapprox(r, 0.0)) {
    return 2.0 * (base_gaussian(t, sigma) + t * d_base_gaussian_dr(t, sigma)) / (sigma * sigma);
  } else {
    return -(d_base_gaussian_dr(r - t, sigma) + d_base_gaussian_dr(r + t, sigma)) / r;
  }
}

template <typename T> inline auto d_exact_gaussian_solution_dr(T t, T r, T sigma) -> T {
  if (isapprox(r, 0.0)) {
    return T{0};
  } else {
    return (-base_gaussian(r - t, sigma) + base_gaussian(r + t, sigma)
            + r * (d_base_gaussian_dr(r - t, sigma) - d_base_gaussian_dr(r + t, sigma)))
           / (r * r);
  }
}

} // namespace fckg

#endif // FC_KLEIN_GORDON_INITIAL_CONDITIONS_HPP
