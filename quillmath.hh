#ifndef QUILLMATH_HH
#define QUILLMATH_HH

/*!
 * \file quillmath.hh
 * \brief Basic mathematical function used in Quill
 */

#ifdef __SSE2__
#include <emmintrin.h>
#endif

/*!
 * \brief Exponential function
 *
 * Return the exponential of \a x, \f$e^x\f$
 */
double qexp(double x);
#ifdef __SSE2__
__m128d qexp(__m128d x);
#endif

/*!
 * \brief Natural logarithm
 *
 * Return the natural logarithm of \a x.
 */
double qlog(double x);
#ifdef __SSE2__
__m128d qlog(__m128d x);
#endif

/*!
 * \brief Error function
 *
 * Return the error function \f$\erf(x) = \frac{2}{\sqrt\pi}\int_0^x e^-t^2 dt\f$
 */
double qerf(double x);
#ifdef __SSE2__
__m128d qerf(__m128d x);
#endif

/*!
 * \brief Log gamma function
 *
 * Return the natural logarithm of the absolute value of the Gamma function.
 */
double qlgamma(double x);
#ifdef __SSE2__
__m128d qlgamma(__m128d x);
#endif

#ifdef __SSE2__
namespace
{
	
inline __m128d qabs(__m128d x)
{
	return _mm_andnot_pd(_mm_set1_pd(-0.0), x);
}

inline __m128d select(__m128d mask, __m128d a, __m128d b)
{
	return _mm_or_pd(_mm_and_pd(mask, a), _mm_andnot_pd(mask, b));
}

} // namespace
#endif

#endif // QUILLMATH_HH