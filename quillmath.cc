#include <limits>
#include <cmath>
#include "quillmath.hh"
#include "constants.hh"

// qexp(), qlog(), qerf() based on code from Cephes, http://netlib.org/cephes

static const double MAXLOG =  7.09782712893383996843e2;
static const double MINLOG = -7.08396418532264106224e2;
static const double MAXLGM = 2.556348e305;

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

static const double LOG_P[] = {
	1.01875663804580931796e-4,
	4.97494994976747001425e-1,
	4.70579119878881725854e0,
	1.44989225341610930846e1,
	1.79368678507819816313e1,
	7.70838733755885391666e0,
};
static const double LOG_Q[] = {
	/* 1.00000000000000000000E0, */
	1.12873587189167450590e1,
	4.52279145837532221105e1,
	8.29875266912776603211e1,
	7.11544750618563894466e1,
	2.31251620126765340583e1,
};
static const double LOG_R[] = {
	-7.89580278884799154124e-1,
	1.63866645699558079767e1,
	-6.41409952958715622951e1,
};
static const double LOG_S[] = {
	/* 1.00000000000000000000e0,*/
	-3.56722798256324312549e1,
	3.12093766372244180303e2,
	-7.69691943550460008604e2,
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

static const double LGAMMA_A[] = {
	8.11614167470508450300e-4,
	-5.95061904284301438324e-4,
	7.93650340457716943945e-4,
	-2.77777777730099687205e-3,
	8.33333333333331927722e-2
};
static const double LGAMMA_B[] = {
	-1.37825152569120859100e3,
	-3.88016315134637840924e4,
	-3.31612992738871184744e5,
	-1.16237097492762307383e6,
	-1.72173700820839662146e6,
	-8.53555664245765465627e5
};
static const double LGAMMA_C[] = {
	/* 1.00000000000000000000e0, */
	-3.51815701436523470549e2,
	-1.70642106651881159223e4,
	-2.20528590553854454839e5,
	-1.13933444367982507207e6,
	-2.53252307177582951285e6,
	-2.01889141433532773231e6
};

int qsigngam;

namespace {

inline double evalPoly(double x, const double *C, int n)
{
	double p = *C;
	for (++C; n > 0; --n, ++C)
		p = x*p + *C;
	return p;
}

inline double evalPoly1(double x, const double *C, int n)
{
	double p = x + *C;
	for (--n, ++C; n > 0; --n, ++C)
		p = x*p + *C;
	return p;
}

} // namespace

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

double qlog(double x)
{
	int e;
	double y, z;

	if (std::isnan(x))
		return x;
	if (std::isinf(x) == 1)
		return x;
	if (x == 0)
		return -std::numeric_limits<double>::infinity();
	if (x < 0)
		return std::numeric_limits<double>::quiet_NaN();

	/* separate mantissa from exponent */
	/* Note, frexp is used so that denormal numbers
	* will be handled properly.
	*/
	x = std::frexp(x, &e);

	/* logarithm using log(x) = z + z**3 P(z)/Q(z),
	* where z = 2(x-1)/x+1)
	*/
	if (e > 2 || e < -2)
	{
		if (x < M_SQRT1_2)
		{ /* 2( 2x-1 )/( 2x+1 ) */
			--e;
			z = x - 0.5;
			y = 0.5 * z + 0.5;
		}
		else
		{ /*  2 (x-1)/(x+1)   */
			z = x - 0.5;
			z -= 0.5;
			y = 0.5 * x + 0.5;
		}
		x = z / y;

		/* rational form */
		z = x * x;
		z = x * (z * evalPoly(z, LOG_R, 2) / evalPoly1(z, LOG_S, 3));
		y = e;
		z -= y * 2.121944400546905827679e-4;
		z += x;
		z += e * 0.693359375;
		return z;
	}

	/* logarithm using log(1+x) = x - .5x**2 + x**3 P(x)/Q(x) */
	if (x < M_SQRT1_2)
	{
		--e;
		x = 2*x - 1.0; /*  2x - 1  */
	}
	else
	{
		x = x - 1.0;
	}

	/* rational form */
	z = x * x;
	y = x * (z * evalPoly(x, LOG_P, 5) / evalPoly1(x, LOG_Q, 5));
	if (e)
		y -= e * 2.121944400546905827679e-4;
	y -= 0.5*z;   /*  y - 0.5 * z  */
	z = x + y;
	if (e)
		z += e * 0.693359375;
	return z;
}

namespace
{

double qerfc_xgt1(double x)
{
	double mxx = -x*x;
	if (mxx < MINLOG)
		return x < 0 ? 2 : 0;

	double p, q;
	double ax = std::abs(x);
	if (ax < 8)
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

} // namespace

double qerf(double x)
{
	if (std::abs(x) > 1)
		return 1 - qerfc_xgt1(x);
	double xx = x*x;
	return x * evalPoly(xx, ERF_T, 4) / evalPoly1(xx, ERF_U, 5);
}


double qlgamma(double x)
{
	double p, q, u, w, z;
	int i;

	qsigngam = 1;
	if (std::isnan(x))
		return x;
	if (!std::isfinite(x))
		return std::numeric_limits<double>::infinity();

	if (x < -34.0)
	{
		q = -x;
		w = qlgamma(q); /* note this modifies qsigngam! */
		p = floor(q);
		if (p == q)
			return std::numeric_limits<double>::infinity();	

		i = p;
		if ((i & 1) == 0)
			qsigngam = -1;
		else
			qsigngam = 1;
		z = q - p;
		if( z > 0.5 )
		{
			p += 1.0;
			z = p - q;
		}
		z = q * std::sin(M_PI * z);
		if (z == 0.0)
			return std::numeric_limits<double>::infinity();

		z = Constants::log_pi - qlog(z) - w;
		return z;
	}

	if (x < 13.0)
	{
		z = 1.0;
		p = 0.0;
		u = x;
		while (u >= 3.0)
		{
			p -= 1.0;
			u = x + p;
			z *= u;
		}
		while (u < 2.0)
		{
			if (u == 0.0)
				return std::numeric_limits<double>::infinity();
			z /= u;
			p += 1.0;
			u = x + p;
		}
		if (z < 0.0)
		{
			qsigngam = -1;
			z = -z;
		}
		else
		{
			qsigngam = 1;
		}
		if (u == 2.0)
			return qlog(z);

		p -= 2.0;
		x += p;
		p = x * evalPoly(x, LGAMMA_B, 5) / evalPoly1(x, LGAMMA_C, 6);
		return qlog(z) + p;
	}

	if (x > MAXLGM)
		return qsigngam * std::numeric_limits<double>::infinity();

	q = (x - 0.5) * qlog(x) - x + Constants::log_sqrt_2pi;
	if (x > 1.0e8)
		return q;

	p = 1.0 / (x * x);
	if (x >= 1000.0)
		q += ((   7.9365079365079365079365e-4 * p
			- 2.7777777777777777777778e-3) *p
			+ 0.0833333333333333333333) / x;
	else
		q += evalPoly(p, LGAMMA_A, 4) / x;
	return q;
}

#ifdef __SSE2__
namespace {

inline __m128d evalPoly(__m128d x, const double *C, int n)
{
	__m128d p = _mm_set1_pd(*C);
	for (++C; n > 0; --n, ++C)
		p = _mm_add_pd(_mm_mul_pd(x, p), _mm_set1_pd(*C));
	return p;
}

inline __m128d evalPoly1(__m128d x, const double *C, int n)
{
	__m128d p = _mm_add_pd(x, _mm_set1_pd(*C));
	for (--n, ++C; n > 0; --n, ++C)
		p = _mm_add_pd(_mm_mul_pd(x, p), _mm_set1_pd(*C));
	return p;
}

inline __m128d qabs(__m128d x)
{
	return _mm_andnot_pd(_mm_set1_pd(-0.0), x);
}

inline __m128d select(__m128d mask, __m128d a, __m128d b)
{
	return _mm_or_pd(_mm_and_pd(mask, a), _mm_andnot_pd(mask, b));
}

inline __m128i select(__m128i mask, __m128i a, __m128i b)
{
	return _mm_or_si128(_mm_and_si128(mask, a), _mm_andnot_si128(mask, b));
}

} // namespace

__m128d qfrexp(__m128d x, __m128i& exp)
{
	__m128i expmask = _mm_set1_epi64x(0x7ff0000000000000LL);
	__m128i expoff = _mm_set1_epi64x(0x3fe0000000000000LL);
	__m128i ix, denorm_mask;
	__m128d spec_mask, denorm_scale, fx;

	// extract exponent bits
	exp = _mm_and_si128(_mm_castpd_si128(x), expmask);

	// check if number is zero, inf, or NaN
	spec_mask = _mm_castsi128_pd(
		_mm_shuffle_epi32(_mm_cmpeq_epi32(exp, expmask),
		_MM_SHUFFLE(3, 3, 1, 1)));          // +-inf or NaN
	spec_mask = _mm_or_pd(spec_mask,
		_mm_cmpeq_pd(x, _mm_setzero_pd())); // zero

	// Scale denormalized numbers. This will also scale zero and adjust its
	// exponent, but that is fixed at the end where the exponent is set
	// to zero.
	denorm_mask = _mm_shuffle_epi32(_mm_cmpeq_epi32(exp, _mm_setzero_si128()),
		_MM_SHUFFLE(3, 3, 1, 1));
	denorm_scale = select(_mm_castsi128_pd(denorm_mask),
		_mm_set1_pd(0x1p54), _mm_set1_pd(1));
	fx = _mm_mul_pd(x, denorm_scale);
	
	// Extract exponent, subtract exponent offset, and adjust for denorm scaling
	ix = _mm_castpd_si128(fx);
	exp = _mm_and_si128(ix, expmask);
	exp = _mm_srl_epi64(exp, _mm_cvtsi32_si128(52));         // e0, 0, e1, 0
	exp = _mm_sub_epi64(exp,
		select(denorm_mask, _mm_set1_epi64x(0x3fe + 54), _mm_set1_epi64x(0x3fe)));

	// Set exponent in fraction to zero
	ix = _mm_andnot_si128(expmask, ix);
	ix = _mm_or_si128(expoff, ix);
	fx = _mm_castsi128_pd(ix);

	// Restore 0, inf and NaN, and set the corresponding exponents to zero
	exp = select(_mm_castpd_si128(spec_mask), _mm_setzero_si128(), exp);
	exp = _mm_shuffle_epi32(exp, _MM_SHUFFLE(3, 1, 2, 0));
	return select(spec_mask, x, fx);
}

__m128d qexp(__m128d x)
{
	__m128d maxmask, minmask, fx, xx, p, q;
	__m128i n;
	unsigned int round_mode = _MM_GET_ROUNDING_MODE();

	maxmask = _mm_cmple_pd(x, _mm_set1_pd(MAXLOG));
	minmask = _mm_cmpge_pd(x, _mm_set1_pd(MINLOG));

	/* Express e**x = e**g 2**n
	*   = e**g e**( n loge(2) )
	*   = e**( g + n loge(2) )
	*/
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

	// Set result for too high values of x to positive infinity
	x = select(maxmask, x, _mm_set1_pd(std::numeric_limits<double>::infinity()));
	// Set result for too low values of x to zero
	return _mm_and_pd(minmask, x);
}

namespace {

__m128d qlog_PQ(__m128d x, __m128i e)
{
	/* logarithm using log(1+x) = x - .5x**2 + x**3 P(x)/Q(x) */

	__m128d mask = _mm_cmplt_pd(x, _mm_set1_pd(M_SQRT1_2));
	x = select(mask, _mm_mul_pd(x, _mm_set1_pd(2)), x);
	x = _mm_sub_pd(x, _mm_set1_pd(1));
	__m128i emask = _mm_shuffle_epi32(_mm_castpd_si128(mask), _MM_SHUFFLE(2, 0, 3, 1));
	e = _mm_sub_epi32(e, select(emask, _mm_set1_epi32(1), _mm_setzero_si128()));
	__m128d ee = _mm_cvtepi32_pd(e);

	__m128d z = _mm_mul_pd(x, x);
	__m128d y = _mm_mul_pd(
		x,
		_mm_div_pd(
			_mm_mul_pd(z, evalPoly(x, LOG_P, 5)),
			evalPoly1(x, LOG_Q, 5)));
	y = _mm_sub_pd(y, _mm_mul_pd(ee, _mm_set1_pd(2.121944400546905827679e-4)));
	y = _mm_sub_pd(y, _mm_mul_pd(z, _mm_set1_pd(0.5)));
	z = _mm_add_pd(_mm_add_pd(x, y), _mm_mul_pd(ee, _mm_set1_pd(0.693359375)));
	return z;
}

__m128d qlog_RS(__m128d x, __m128i e)
{
	/* logarithm using log(x) = z + z**3 P(z)/Q(z),
	 * where z = 2(x-1)/x+1)
	 */
	__m128d half = _mm_set1_pd(0.5);
	__m128d mask = _mm_cmplt_pd(x, _mm_set1_pd(M_SQRT1_2));
	__m128d z = _mm_sub_pd(x, select(mask, half, _mm_set1_pd(1)));
	__m128d y = _mm_add_pd(_mm_mul_pd(select(mask, z, x), half), half);
	__m128i emask = _mm_shuffle_epi32(_mm_castpd_si128(mask), _MM_SHUFFLE(2, 0, 3, 1));
	e = _mm_sub_epi32(e, select(emask, _mm_set1_epi32(1), _mm_setzero_si128()));
	__m128d ee = _mm_cvtepi32_pd(e);
	x = _mm_div_pd(z, y);
	
	z = _mm_mul_pd(x, x);
	y = _mm_mul_pd(
		x,
		_mm_div_pd(
			_mm_mul_pd(z, evalPoly(z, LOG_R, 2)),
			evalPoly1(z, LOG_S, 3)));
	y = _mm_sub_pd(y, _mm_mul_pd(ee, _mm_set1_pd(2.121944400546905827679e-4)));
	z = _mm_add_pd(_mm_add_pd(x, y), _mm_mul_pd(ee, _mm_set1_pd(0.693359375)));
	return z;
}
	
} // namespace

__m128d qlog(__m128d x)
{
	double inf = std::numeric_limits<double>::infinity();
	__m128d specmask, negmask, zeromask, specval, lnx;
	__m128i ix, e, mask;

	// x == inf or x == NaN. This also matches -inf, but that will be
	// overwritten in the next test
	ix = _mm_castpd_si128(_mm_andnot_pd(_mm_set1_pd(-0.0), x));
	ix = _mm_xor_si128(_mm_set1_epi32(0xffffffff),
		_mm_cmplt_epi32(ix, _mm_castpd_si128(_mm_set1_pd(inf))));
	ix = _mm_shuffle_epi32(ix, _MM_SHUFFLE(3, 3, 1, 1));
	specmask = _mm_castsi128_pd(ix);
	specval = x;
	// x < 0
	negmask = _mm_cmplt_pd(x, _mm_setzero_pd());
	specmask = _mm_or_pd(specmask, negmask);
	specval = select(negmask, _mm_set1_pd(std::numeric_limits<double>::quiet_NaN()),
		specval);
	// x == 0
	zeromask = _mm_cmpeq_pd(x, _mm_setzero_pd());
	specmask = _mm_or_pd(specmask, zeromask);
	specval = select(zeromask, _mm_set1_pd(-inf), specval);

	/* separate mantissa from exponent
	 * Note, frexp is used so that denormal numbers
	 * will be handled properly.
	 */
	x = qfrexp(x, e);

	mask = _mm_or_si128(_mm_cmplt_epi32(e, _mm_set1_epi32(-2)),
		_mm_cmpgt_epi32(e, _mm_set1_epi32(2)));
	mask = _mm_shuffle_epi32(mask, _MM_SHUFFLE(1, 1, 0, 0));
	int mbits = _mm_movemask_epi8(mask);
	if (mbits == 0)
		lnx = qlog_PQ(x, e);
	if (mbits == 0xffff)
		lnx = qlog_RS(x, e);
	lnx = select(_mm_castsi128_pd(mask), qlog_RS(x, e), qlog_PQ(x, e));

	return select(specmask, specval, lnx);
}

static __m128d qerfc_xgt1(__m128d x)
{
	__m128d p, q;
	__m128d ax = qabs(x);
	__m128d mask = _mm_cmplt_pd(x, _mm_set1_pd(8));
	int mbits = _mm_movemask_pd(mask);
	if (mbits == 0)
	{
		p = evalPoly(ax, ERFC_R, 5);
		q = evalPoly1(ax, ERFC_S, 6);
	}
	else if (mbits == 3)
	{
		p = evalPoly(ax, ERFC_P, 8);
		q = evalPoly1(ax, ERFC_Q, 8);
	}
	else
	{
		__m128d p1 = evalPoly(ax, ERFC_R, 5);
		__m128d q1 = evalPoly1(ax, ERFC_S, 6);
		__m128d p2 = evalPoly(ax, ERFC_P, 8);
		__m128d q2 = evalPoly1(ax, ERFC_Q, 8);
		p = select(mask, p2, p1);
		q = select(mask, q2, q1);
	}

	__m128d mxx = _mm_xor_pd(_mm_set1_pd(-0.0), _mm_mul_pd(x, x));
	__m128d y = _mm_div_pd(_mm_mul_pd(qexp(mxx), p), q);
	__m128d ty = _mm_sub_pd(_mm_set1_pd(2), y);
	mask = _mm_cmpge_pd(x, _mm_setzero_pd());
	y = select(mask, y, ty);
	return y;
}

__m128d qerf(__m128d x)
{
	__m128d ax = qabs(x);
	__m128d one = _mm_set1_pd(1);
	__m128d mask = _mm_cmpgt_pd(ax, one);
	int mbits = _mm_movemask_pd(mask);
	if (mbits == 0)
	{
		__m128d xx = _mm_mul_pd(x, x);
		return _mm_div_pd(_mm_mul_pd(x, evalPoly(xx, ERF_T, 4)),
			evalPoly1(xx, ERF_U, 5));
	}
	else if (mbits == 3)
	{
		return _mm_sub_pd(one, qerfc_xgt1(x));
	}
	else
	{
		__m128d xx = _mm_mul_pd(x, x);
		__m128d el = _mm_div_pd(_mm_mul_pd(x, evalPoly(xx, ERF_T, 4)),
			evalPoly1(xx, ERF_U, 5));
		__m128d eg = _mm_sub_pd(one, qerfc_xgt1(x));
		return select(mask, eg, el);
	}
}
#endif