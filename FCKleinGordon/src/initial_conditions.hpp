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

template <typename T> inline auto smooth_transition(T sigma, T x, T y, T z) -> T {
  using std::exp;
  using std::sqrt;

  const T r = sqrt(x * x + y * y + z * z);

  if (r > 0 && sigma > r) {
    return 1 / (1 + exp(sigma * (1 / (sigma - r) - 1 / r)));
  } else if (r > 0) {
    return 0;
  } else {
    return 1;
  }
}

template <typename T> inline auto sech(T x) -> T {
  using std::cosh;
  return 1 / cosh(x);
}

template <typename T> inline auto dsmooth_transition_dx(T sigma, T x, T y, T z) -> T {
  using std::exp;
  using std::pow;
  using std::sqrt;

  const T r = sqrt(x * x + y * y + z * z);

  if (r > 0 && sigma > r) {
    return -(x * sigma * (pow(sigma, 2) + 2 * r * (-sigma + r))
             * pow(sech((sigma * (1 / r + 1 / (-sigma + r))) / 2.), 2))
           / (4. * pow(sigma - r, 2) * pow(r, 3));
  } else {
    return 0;
  }
}

template <typename T> inline auto dsmooth_transition_dy(T sigma, T x, T y, T z) -> T {
  using std::exp;
  using std::pow;
  using std::sqrt;

  const T r = sqrt(x * x + y * y + z * z);

  if (r > 0 && sigma > r) {
    return -(y * sigma * (pow(sigma, 2) + 2 * r * (-sigma + r))
             * pow(sech((sigma * (1 / r + 1 / (-sigma + r))) / 2.), 2))
           / (4. * pow(sigma - r, 2) * pow(r, 3));
  } else {
    return 0;
  }
}

template <typename T> inline auto dsmooth_transition_dz(T sigma, T x, T y, T z) -> T {
  using std::exp;
  using std::pow;
  using std::sqrt;

  const T r = sqrt(x * x + y * y + z * z);

  if (r > 0 && sigma > r) {
    return -(z * sigma * (pow(sigma, 2) + 2 * r * (-sigma + r))
             * pow(sech((sigma * (1 / r + 1 / (-sigma + r))) / 2.), 2))
           / (4. * pow(sigma - r, 2) * pow(r, 3));
  } else {
    return 0;
  }
}

template <typename T> inline auto multipole_sum(const T *c, T sigma, T x, T y, T z) -> T {
  using std::pow;
  using std::sqrt;

  return c[0]
         + (y * z * c[5] + 2 * pow(z, 2) * c[6] + x * (y * c[4] + z * c[7])
            + pow(x, 2) * (-c[6] + c[8]) - pow(y, 2) * (c[6] + c[8]))
               / (pow(x, 2) + pow(y, 2) + pow(z, 2) + smooth_transition(sigma, x, y, z))
         + (y * c[1] + z * c[2] + x * c[3])
               / (sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)) + smooth_transition(sigma, x, y, z));
}

template <typename T> inline auto d_multipole_sum_dx(const T *c, T sigma, T x, T y, T z) -> T {
  using std::pow;
  using std::sqrt;

  return (y * c[4] - 2 * x * c[6] + z * c[7] + 2 * x * c[8])
             / (pow(x, 2) + pow(y, 2) + pow(z, 2) + smooth_transition(sigma, x, y, z))
         + c[3] / (sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)) + smooth_transition(sigma, x, y, z))
         - ((y * z * c[5] + 2 * pow(z, 2) * c[6] + x * (y * c[4] + z * c[7])
             + pow(x, 2) * (-c[6] + c[8]) - pow(y, 2) * (c[6] + c[8]))
            * (2 * x + dsmooth_transition_dx(sigma, x, y, z)))
               / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) + smooth_transition(sigma, x, y, z), 2)
         - ((y * c[1] + z * c[2] + x * c[3])
            * (x / sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)) + dsmooth_transition_dx(sigma, x, y, z)))
               / pow(sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)) + smooth_transition(sigma, x, y, z),
                     2);
}

template <typename T> inline auto d_multipole_sum_dy(const T *c, T sigma, T x, T y, T z) -> T {
  using std::pow;
  using std::sqrt;

  return (x * c[4] + z * c[5] - 2 * y * (c[6] + c[8]))
             / (pow(x, 2) + pow(y, 2) + pow(z, 2) + smooth_transition(sigma, x, y, z))
         + c[1] / (sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)) + smooth_transition(sigma, x, y, z))
         - ((y * z * c[5] + 2 * pow(z, 2) * c[6] + x * (y * c[4] + z * c[7])
             + pow(x, 2) * (-c[6] + c[8]) - pow(y, 2) * (c[6] + c[8]))
            * (2 * y + dsmooth_transition_dy(sigma, x, y, z)))
               / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) + smooth_transition(sigma, x, y, z), 2)
         - ((y * c[1] + z * c[2] + x * c[3])
            * (y / sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)) + dsmooth_transition_dy(sigma, x, y, z)))
               / pow(sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)) + smooth_transition(sigma, x, y, z),
                     2);
}

template <typename T> inline auto d_multipole_sum_dz(const T *c, T sigma, T x, T y, T z) -> T {
  using std::pow;
  using std::sqrt;

  return (y * c[5] + 4 * z * c[6] + x * c[7])
             / (pow(x, 2) + pow(y, 2) + pow(z, 2) + smooth_transition(sigma, x, y, z))
         + c[2] / (sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)) + smooth_transition(sigma, x, y, z))
         - ((y * z * c[5] + 2 * pow(z, 2) * c[6] + x * (y * c[4] + z * c[7])
             + pow(x, 2) * (-c[6] + c[8]) - pow(y, 2) * (c[6] + c[8]))
            * (2 * z + dsmooth_transition_dz(sigma, x, y, z)))
               / pow(pow(x, 2) + pow(y, 2) + pow(z, 2) + smooth_transition(sigma, x, y, z), 2)
         - ((y * c[1] + z * c[2] + x * c[3])
            * (z / sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)) + dsmooth_transition_dz(sigma, x, y, z)))
               / pow(sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2)) + smooth_transition(sigma, x, y, z),
                     2);
}

} // namespace fckg

#endif // FC_KLEIN_GORDON_INITIAL_CONDITIONS_HPP
