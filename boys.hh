#ifndef BOYS_HH
#define BOYS_HH

/*!
 * \file boys.hh
 * \brief Functions to evaluate the Boys function
 */

#include "quillmath.hh"

/*!
 * \brief Evaluate the Boys function
 *
 * Evaluate the Boys function
 * \f{eqnarray*}
 *    F_m(t) &=& \int_0^1 \exp(-t s^2) s^{2m} ds\\
 *           &=& \frac{1}{2t^{m+1/2}} \Gamma(m+1/2) P(m+1/2, t)\\
 *           &=& \frac{1}{2t^{m+1/2}} \Gamma(m+1/2) [1-Q(m+1/2, t)]
 * \f}
 * for a single value of \a m and \a t. The third argument \a expmt is
 * \f$e^{-t}\f$. This function can be used when the exponential is already
 * known, and saves the need to recompute it.
 */
double Fm(int m, double t, double expmt);
#ifdef __SSE2__
__m128d Fm(int m, __m128d t, __m128d expmt);
#endif

namespace
{

/*!
 * \brief Evaluate the Boys function
 *
 * Evaluate the Boys function
 * \f{eqnarray*}
 *    F_m(t) &=& \int_0^1 \exp(-t s^2) s^{2m} ds\\
 *           &=& \frac{1}{2t^{m+1/2}} \Gamma(m+1/2) P(m+1/2, t)\\
 *           &=& \frac{1}{2t^{m+1/2}} \Gamma(m+1/2) [1-Q(m+1/2, t)]
 * \f}
 * for a single value of \a m and \a t.
 */
double Fm(int m, double t)
{
	return ::Fm(m, t, qexp(-t));
}
#ifdef __SSE2__
__m128d Fm(int m, __m128d t)
{
	return ::Fm(m, t, qexp(_mm_xor_pd(_mm_set1_pd(-0.0), t)));
}
#endif

} // namespace

#endif // BOYS_HH
