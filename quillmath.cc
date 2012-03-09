#include <limits>
#include <cmath>
#include "quillmath.hh"

// qexp(), qerf() based on code from Cephes, http://netlib.org/cephes

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

static const double ERFC_P[] = {
	2.46196981473530512524E-10,
	5.64189564831068821977E-1,
	7.46321056442269912687E0,
	4.86371970985681366614E1,
	1.96520832956077098242E2,
	5.26445194995477358631E2,
	9.34528527171957607540E2,
	1.02755188689515710272E3,
	5.57535335369399327526E2
};
static const double ERFC_Q[] = {
	/* 1.00000000000000000000E0,*/
	1.32281951154744992508E1,
	8.67072140885989742329E1,
	3.54937778887819891062E2,
	9.75708501743205489753E2,
	1.82390916687909736289E3,
	2.24633760818710981792E3,
	1.65666309194161350182E3,
	5.57535340817727675546E2
};
static const double ERFC_R[] = {
	5.64189583547755073984E-1,
	1.27536670759978104416E0,
	5.01905042251180477414E0,
	6.16021097993053585195E0,
	7.40974269950448939160E0,
	2.97886665372100240670E0
};
static const double ERFC_S[] = {
	/* 1.00000000000000000000E0,*/
	2.26052863220117276590E0,
	9.39603524938001434673E0,
	1.20489539808096656605E1,
	1.70814450747565897222E1,
	9.60896809063285878198E0,
	3.36907645100081516050E0
};
static const double ERF_T[] = {
	9.60497373987051638749e0,
	9.00260197203842689217e1,
	2.23200534594684319226e3,
	7.00332514112805075473e3,
	5.55923013010394962768e4
};
static const double ERF_U[] = {
	/* 1.00000000000000000000E0,*/
	3.35617141647503099647e1,
	5.21357949780152679795e2,
	4.59432382970980127987e3,
	2.26290000613890934246e4,
	4.92673942608635921086e4
};

static inline double evalPoly(double x, const double *C, int n)
{
	double p = *C;
	for (++C; n > 0; --n, ++C)
		p = x*p + *C;
	return p;
}

static inline double evalPoly1(double x, const double *C, int n)
{
	double p = x + *C;
	for (--n, ++C; n > 0; --n, ++C)
		p = x*p + *C;
	return p;
}

double qexp(double x)
{
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
	double px = x * evalPoly(xx, EXP_P, 2);
	x = px / (evalPoly(xx, EXP_Q, 3) - px);
	x = 1 + 2*x;

	/* multiply by power of 2 */
	return scalbn(x, n);
}

static double qerfc_xgt1(double x)
{
	double mxx = -x*x;
	if (mxx < MINLOG)
		return x < 0 ? 2 : 0;

	double ax = std::abs(x);
	double p, q;
	if (std::abs(x) < 8)
	{
		p = evalPoly(ax, ERFC_P, 8);
		q = evalPoly1(ax, ERFC_Q, 8);
	}
	else
	{
		p = evalPoly(ax, ERFC_R, 5);
		q = evalPoly1(ax, ERFC_S, 6);
	}

	double y = (qexp(mxx) * p) / q;
	return x < 0 ? 2-y : y;
}

double qerf(double x)
{
	if (std::abs(x) > 1)
		return 1 - qerfc_xgt1(x);
	double xx = x*x;
	return x * evalPoly(xx, ERF_T, 4) / evalPoly1(xx, ERF_U, 5);
}

#ifdef __SSE2__
static inline __m128d evalPoly(__m128d x, const double *C, int n)
{
	__m128d p = _mm_set1_pd(*C);
	for (++C; n > 0; --n, ++C)
		p = _mm_add_pd(_mm_mul_pd(x, p), _mm_set1_pd(*C));
	return p;
}

static inline __m128d evalPoly1(__m128d x, const double *C, int n)
{
	__m128d p = _mm_add_pd(x, _mm_set1_pd(*C));
	for (--n, ++C; n > 0; --n, ++C)
		p = _mm_add_pd(_mm_mul_pd(x, p), _mm_set1_pd(*C));
	return p;
}

__m128d qexp(__m128d x)
{
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
	p = _mm_mul_pd(x, evalPoly(xx, EXP_P, 2));
	q = evalPoly(xx, EXP_Q, 3);
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