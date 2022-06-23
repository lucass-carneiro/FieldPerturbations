#ifndef FC_KLEIN_GORDON_INITIAL_CONDITIONS_HPP
#define FC_KLEIN_GORDON_INITIAL_CONDITIONS_HPP

#include <cmath>
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

template <typename T> auto sign(T val) -> T {
  if (val < T{0}) {
    return T{-1};
  } else if (val > T{0}) {
    return T{1};
  } else {
    return T{0};
  }
}

template <typename T> inline auto multipole_sum(const T *c, T x, T y, T z) -> T {
  using std::pow;
  using std::sqrt;

  const T r2 = x * x + y * y + z * z;
  const T r = sqrt(r2);

  if (isapprox(r, 0.0)) {
    return c[0];
  }

  if (isapprox(x, 0.0) && isapprox(y, 0.0)) {
    return c[0] + c[6] + c[2] * sign(z);

  } else if (isapprox(x, 0.0) && !isapprox(y, 0.0)) {
    return c[0]
           + (2 * y * (sqrt(pow(y, 2) + pow(z, 2)) * c[1] + 3 * z * (2 * c[5] / 3))
              + 2 * z * (sqrt(pow(y, 2) + pow(z, 2)) * c[2] + z * c[6])
              - pow(y, 2) * (c[6] + 6 * (c[8] / 3)))
                 / (2. * (pow(y, 2) + pow(z, 2)));

  } else if (!isapprox(x, 0.0) && isapprox(y, 0.0)) {
    return c[0]
           + (2 * z * sqrt(pow(x, 2) + pow(z, 2)) * c[2]
              + 2 * x * sqrt(pow(x, 2) + pow(z, 2)) * c[3] - (pow(x, 2) - 2 * pow(z, 2)) * c[6]
              + 6 * x * (z * (2 * c[7] / 3) + x * (c[8] / 3)))
                 / (2. * (pow(x, 2) + pow(z, 2)));

  } else if (!isapprox(x, 0.0) && !isapprox(y, 0.0)) {
    return c[0] + (y * c[1]) / r + (z * c[2]) / r + (x * c[3]) / r - c[6] / 2. - 3 * (c[8] / 3)
           + (3
              * (4 * x * y * (c[4] / 3) + 2 * y * z * (2 * c[5] / 3) + pow(z, 2) * c[6]
                 + 2 * x * z * (2 * c[7] / 3) + 2 * (2 * pow(x, 2) + pow(z, 2)) * (c[8] / 3)))
                 / (2. * (r2));
  } else {
    return T{0};
  }
}

template <typename T> inline auto d_multipole_sum_dx(const T *c, T x, T y, T z) -> T {
  using std::pow;
  using std::sqrt;

  const T r2 = x * x + y * y + z * z;
  const T r = sqrt(r2);

  if (isapprox(r, 0.0)) {
    return c[0];
  }

  if (isapprox(x, 0.0) && isapprox(y, 0.0)) {
    return 0.0;

  } else if (isapprox(x, 0.0) && !isapprox(y, 0.0)) {
    return (sqrt(pow(y, 2) + pow(z, 2)) * c[3] + 6 * y * (c[4] / 3) + 3 * z * (2 * c[7] / 3))
           / (pow(y, 2) + pow(z, 2));

  } else if (!isapprox(x, 0.0) && isapprox(y, 0.0)) {
    return (z
            * (-(pow(x, 3) * c[2])
               + pow(x, 2) * (z * c[3] - 3 * sqrt(pow(x, 2) + pow(z, 2)) * (2 * c[7] / 3))
               + pow(z, 2) * (z * c[3] + 3 * sqrt(pow(x, 2) + pow(z, 2)) * (2 * c[7] / 3))
               - x * z * (z * c[2] + 3 * sqrt(pow(x, 2) + pow(z, 2)) * (c[6] - 2 * (c[8] / 3)))))
           / pow(pow(x, 2) + pow(z, 2), 2.5);

  } else if (!isapprox(x, 0.0) && !isapprox(y, 0.0)) {
    return (-3 * pow(x, 2) * (2 * y * (c[4] / 3) + z * (2 * c[7] / 3))
            + (pow(y, 2) + pow(z, 2)) * (r * c[3] + 6 * y * (c[4] / 3) + 3 * z * (2 * c[7] / 3))
            - x
                  * (y * r * c[1] + z * r * c[2] + 6 * y * z * (2 * c[5] / 3)
                     + 3 * pow(z, 2) * (c[6] - 2 * (c[8] / 3)) - 12 * pow(y, 2) * (c[8] / 3)))
           / pow(r2, 2);
  } else {
    return T{0};
  }
}

template <typename T> inline auto d_multipole_sum_dy(const T *c, T x, T y, T z) -> T {
  using std::pow;
  using std::sqrt;

  const T r2 = x * x + y * y + z * z;
  const T r = sqrt(r2);

  if (isapprox(r, 0.0)) {
    return c[0];
  }

  if (isapprox(x, 0.0) && isapprox(y, 0.0)) {
    return 0.0;

  } else if (isapprox(x, 0.0) && !isapprox(y, 0.0)) {
    return (z
            * (-(pow(y, 3) * c[2])
               + pow(y, 2) * (z * c[1] - 3 * sqrt(pow(y, 2) + pow(z, 2)) * (2 * c[5] / 3))
               + pow(z, 2) * (z * c[1] + 3 * sqrt(pow(y, 2) + pow(z, 2)) * (2 * c[5] / 3))
               - y * z * (z * c[2] + 3 * sqrt(pow(y, 2) + pow(z, 2)) * (c[6] + 2 * (c[8] / 3)))))
           / pow(pow(y, 2) + pow(z, 2), 2.5);

  } else if (!isapprox(x, 0.0) && isapprox(y, 0.0)) {
    return (sqrt(pow(x, 2) + pow(z, 2)) * c[1] + 6 * x * (c[4] / 3) + 3 * z * (2 * c[5] / 3))
           / (pow(x, 2) + pow(z, 2));

  } else if (!isapprox(x, 0.0) && !isapprox(y, 0.0)) {
    return -((-((pow(x, 2) + pow(z, 2)) * r * c[1]) + y * z * r * c[2] + x * y * r * c[3]
              - 6 * x * (pow(x, 2) - pow(y, 2) + pow(z, 2)) * (c[4] / 3)
              - 3 * z * (pow(x, 2) - pow(y, 2) + pow(z, 2)) * (2 * c[5] / 3)
              + 3 * y * pow(z, 2) * c[6] + 6 * x * y * z * (2 * c[7] / 3)
              + 6 * y * (2 * pow(x, 2) + pow(z, 2)) * (c[8] / 3))
             / pow(r2, 2));
  } else {
    return T{0};
  }
}

template <typename T> inline auto d_multipole_sum_dz(const T *c, T x, T y, T z) -> T {
  using std::pow;
  using std::sqrt;

  const T r2 = x * x + y * y + z * z;
  const T r = sqrt(r2);

  if (isapprox(r, 0.0)) {
    return c[0];
  }

  if (isapprox(x, 0.0) && isapprox(y, 0.0)) {
    return 0.0;

  } else if (isapprox(x, 0.0) && !isapprox(y, 0.0)) {
    return (y
            * (pow(y, 3) * c[2]
               + pow(y, 2) * (-(z * c[1]) + 3 * sqrt(pow(y, 2) + pow(z, 2)) * (2 * c[5] / 3))
               - pow(z, 2) * (z * c[1] + 3 * sqrt(pow(y, 2) + pow(z, 2)) * (2 * c[5] / 3))
               + y * z * (z * c[2] + 3 * sqrt(pow(y, 2) + pow(z, 2)) * (c[6] + 2 * (c[8] / 3)))))
           / pow(pow(y, 2) + pow(z, 2), 2.5);

  } else if (!isapprox(x, 0.0) && isapprox(y, 0.0)) {
    return (x
            * (pow(x, 3) * c[2]
               + pow(x, 2) * (-(z * c[3]) + 3 * sqrt(pow(x, 2) + pow(z, 2)) * (2 * c[7] / 3))
               - pow(z, 2) * (z * c[3] + 3 * sqrt(pow(x, 2) + pow(z, 2)) * (2 * c[7] / 3))
               + x * z * (z * c[2] + 3 * sqrt(pow(x, 2) + pow(z, 2)) * (c[6] - 2 * (c[8] / 3)))))
           / pow(pow(x, 2) + pow(z, 2), 2.5);

  } else if (!isapprox(x, 0.0) && !isapprox(y, 0.0)) {
    return -((y * z * r * c[1] + pow(z, 2) * r * c[2] - pow(r2, 1.5) * c[2] + x * z * r * c[3]
              + 12 * x * y * z * (c[4] / 3)
              - 3 * y * (pow(x, 2) + pow(y, 2) - pow(z, 2)) * (2 * c[5] / 3)
              - 3 * (pow(x, 2) + pow(y, 2)) * z * c[6]
              - 3 * x * (pow(x, 2) + pow(y, 2) - pow(z, 2)) * (2 * c[7] / 3)
              + 6 * (x - y) * (x + y) * z * (c[8] / 3))
             / pow(r2, 2));
  } else {
    return T{0};
  }
}

} // namespace fckg

#endif // FC_KLEIN_GORDON_INITIAL_CONDITIONS_HPP
