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
 * \brief Error function
 *
 * Return the error function \f$\erf(x) = \frac{2}{\sqrt\pi}\int_0^x e^-t^2 dt\f$
 */
double qerf(double x);

#endif // QUILLMATH_HH