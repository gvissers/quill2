#include <limits>
#include <cmath>
#include "quillmath.hh"

static const double MAXLOG =  7.09782712893383996843e2;
static const double MINLOG = -7.08396418532264106224e2;

static const double EXP_C1 = 6.93145751953125e-1; // C1+C2 = log(2)
static const double EXP_C2 = 1.42860682030941723212e-6;
static const double EXP_P[] = {
	1.26177193074810590878e-4,
	3.02994407707441961300e-2,
	9.99999999999999999910e-1,
};
static const double EXP_Q[] = {
	3.00198505138664455042e-6,
	2.52448340349684104192e-3,
	2.27265548208155028766e-1,
	2.00000000000000000009e0,
};

double qexp(double x)
{
	// Taken from Cephes, http://netlib.org/cephes
	if (std::isnan(x))
		return x;

	/* Express e**x = e**g 2**n
	*   = e**g e**( n loge(2) )
	*   = e**( g + n loge(2) )
	*/
	int n;
	if (x >= 0)
	{
		if (x > MAXLOG)
			return std::numeric_limits<double>::infinity();
		n = int(M_LOG2E*x + 0.5);
	}
	else
	{
		if (x < MINLOG)
			return 0;
		n = int(M_LOG2E*x - 0.5);
	}
	x -= n * EXP_C1;
	x -= n * EXP_C2;

	/* rational approximation for exponential
	* of the fractional part:
	* e**x = 1 + 2x P(x**2)/( Q(x**2) - P(x**2) )
	*/
	double xx = x * x;
	double px = x * (EXP_P[2] + xx*(EXP_P[1] + xx*EXP_P[0]));
	x = px / (EXP_Q[3] + xx*(EXP_Q[2] + xx*(EXP_Q[1] + xx*EXP_Q[0])) - px);
	x = 1 + 2*x;

	/* multiply by power of 2 */
	return scalbn(x, n);
}

#ifdef __SSE2__
__m128d qexp(__m128d x)
{
	// Taken from Cephes, http://netlib.org/cephes
	__m128d fx, xx, p, q;
	__m128i n;
	unsigned int round_mode = _MM_GET_ROUNDING_MODE();
	
	/* Express e**x = e**g 2**n
	*   = e**g e**( n loge(2) )
	*   = e**( g + n loge(2) )
	*/
	x = _mm_min_pd(x, _mm_set1_pd(MAXLOG));
	x = _mm_max_pd(x, _mm_set1_pd(MINLOG));

	fx = _mm_mul_pd(x, _mm_set1_pd(M_LOG2E));
	_MM_SET_ROUNDING_MODE(_MM_ROUND_NEAREST);
	n = _mm_cvtpd_epi32(fx);
	_MM_SET_ROUNDING_MODE(round_mode);
	fx = _mm_cvtepi32_pd(n);
	x = _mm_sub_pd(x, _mm_mul_pd(fx, _mm_set1_pd(EXP_C1)));
	x = _mm_sub_pd(x, _mm_mul_pd(fx, _mm_set1_pd(EXP_C2)));

	/* rational approximation for exponential
	* of the fractional part:
	* e**x = 1 + 2x P(x**2)/( Q(x**2) - P(x**2) )
	*/
	xx = _mm_mul_pd(x, x);
	p = _mm_mul_pd(xx, _mm_set1_pd(EXP_P[0]));
	p = _mm_mul_pd(xx, _mm_add_pd(p, _mm_set1_pd(EXP_P[1])));
	p = _mm_mul_pd(x, _mm_add_pd(p, _mm_set1_pd(EXP_P[2])));
	q = _mm_mul_pd(xx, _mm_set1_pd(EXP_Q[0]));
	q = _mm_mul_pd(xx, _mm_add_pd(q, _mm_set1_pd(EXP_Q[1])));
	q = _mm_mul_pd(xx, _mm_add_pd(q, _mm_set1_pd(EXP_Q[2])));
	q = _mm_add_pd(q, _mm_set1_pd(EXP_Q[3]));
	x = _mm_div_pd(p, _mm_sub_pd(q, p));
	x = _mm_add_pd(x, _mm_add_pd(x, _mm_set1_pd(1)));

	/* multiply by power of 2 */
	n = _mm_add_epi32(n, _mm_set_epi32(0, 0, 0x3ff, 0x3ff));
	n = _mm_shuffle_epi32(n, _MM_SHUFFLE(3, 1, 2, 0));
	n = _mm_sll_epi64(n, _mm_cvtsi32_si128(52));
	x = _mm_mul_pd(x, _mm_castsi128_pd(n));

	return x;
}
#endif