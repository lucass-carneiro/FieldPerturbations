#ifndef ARITHMETIC_HPP
#define ARITHMETIC_HPP

namespace KleinGordonX {

/****************************************************************
 * template <size_t N, typename T>                              *
 * constexpr inline T power(T x);                               *
 *                                                              *
 * This function computes the power of a number. If the number  *
 * is know at compile time, the compiler computes the result    *
 * of the exponentiation. If not, the power gets "unrolled" and *
 * should be equivalent to writing each term in the power       *
 * sequency by hand.                                            *
 *                                                              *
 * Input N: The power                                           *
 * Input x: The number.                                         *
 *                                                              *
 * Output: The Nth power of the  number x.                      *
 ****************************************************************/
template <size_t N> constexpr inline CCTK_REAL power(CCTK_REAL x) { return x * power<N - 1>(x); }

template <> constexpr inline CCTK_REAL power<0>(CCTK_REAL x) { return 1.0; }

/****************************************************************
 * inline char todir(CCTK_INT d)                                *
 *                                                              *
 * This function converts an integer into a direction character *
 * this is usefull in output messages.                          *
 *                                                              *
 * Input: An integer                                            *
 *                                                              *
 * Output: The direction                                        *
 ****************************************************************/
inline char todir(CCTK_INT d) {
  if (d == 0)
    return 'x';
  else if (d == 1)
    return 'y';
  else
    return 'z';
}

} // namespace KleinGordonX

#endif // ARITHMETIC_HPP
