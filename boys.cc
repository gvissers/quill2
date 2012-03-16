#include <limits>
#include <cmath>
#include "boys.hh"
#include "exceptions.hh"

static const int maxiter = 200;
static const double fpmin = std::numeric_limits<double>::min();
static const double ifpmin = 1.0 / std::numeric_limits<double>::min();
static const double depsilon = std::numeric_limits<double>::epsilon();

namespace {

double lgamma_cache(int m)
{
	static double lgammas[100];
	static int used = 0;

	if (m >= used)
	{
		for (int i = used; i <= m; ++i)
			lgammas[i] = qlgamma(i+0.5);
		used = m+1;
	}
	return lgammas[m];
}
	
double Fm_P_series_pos_a_t(double a, double t, double expmt)
{
	double ap = a,
		del = 0.5 / a,
		sum = del;
	int i;

	int nlow = t > a ? t - a: 0;
	for (i = 0; i < nlow; i++)
	{
		ap += 1;
		del *= t / ap;
		sum += del;
	}
	for ( ; i < maxiter; i++)
	{
		ap += 1;
		del *= t / ap;
		sum += del;
		if (del < sum * depsilon)
			return sum * expmt;
	}

	throw NoConvergence();
}

double Fm_Q_cf(double a, double t, double expmt)
{
	double b = t + 1 - a,
		c = ifpmin,
		d = 1.0 / b,
		h = 0.5 * d;

	for (int i = 1; i <= maxiter; i++)
	{
		double an = -i * (i-a), del;
		b += 2;
		d = an*d + b;
		if (std::abs(d) < fpmin)
			d = fpmin;
		c = b + an / c;
		if (std::abs(c) < fpmin)
			c = fpmin;
		d = 1.0 / d;
		del = d * c;
		h *= del;
		if (std::abs(del - 1.0) < depsilon)
			return 0.5 * qexp(lgamma_cache(int(a)) - a*qlog(t)) - h * expmt;
	}

	throw NoConvergence();
}

} // namespace

double Fm(int m, double t, double expmt)
{
	// compute the integral
	//      F_m(t) = \int_0^1 \exp[-t s^2] s^{2m} ds
	//             = 1/(2*t^{m+1/2}) \Gamma(m+1/2) P(m+1/2, t)
	//             = 1/(2*t^{m+1/2}) \Gamma(m+1/2) [1-Q(m+1/2, t)]
	// by using either the series expansion of the incomplete gamma
	// function P(m+0.5,t) or the continued fraction expansion of
	// Q(m+0.5,t) = 1-P(m+0.5,t). Algorithms for the expansions are
	// taken from Press, Teukolsky, Vetterling, Flannery, Numerical Recipes
	// in C, the Art of Scientific Computing, second edition.
	double a = m + 0.5;
	if (std::abs(t) < depsilon)
	{
		return 0.5 / a;
	}
	else if (m == 0)
	{
		double st = std::sqrt(t);
		return qerf(st) / (M_2_SQRTPI * st);
	}
	else if (t < 20 || t < 0.5*a)
	{
		return Fm_P_series_pos_a_t(a, t, expmt);
	}
	else
	{
		return Fm_Q_cf(a, t, expmt);
	}
}

#ifdef __SSE2__
namespace {

__m128d Fm_P_series_pos_a_t(double a, __m128d t, __m128d expmt)
{
	__m128d one = _mm_set1_pd(1);
	__m128d ap = _mm_set1_pd(a);
	__m128d del = _mm_div_pd(_mm_set1_pd(0.5), ap);
	__m128d sum = del;
	__m128d eps = _mm_set1_pd(depsilon);
	__m128d mask;
	
	for (int i = 0; i < maxiter; i++)
	{
		ap = _mm_add_pd(ap, one);
		del = _mm_mul_pd(del, _mm_div_pd(t, ap));
		sum = _mm_add_pd(sum, del);
		mask = _mm_cmplt_pd(del, _mm_mul_pd(sum, eps));
		if (_mm_movemask_pd(mask) == 3)
			return _mm_mul_pd(sum, expmt);
	}

	throw NoConvergence();
}

__m128d Fm_Q_cf(double a, __m128d t, __m128d expmt)
{
	__m128d one = _mm_set1_pd(1);
	__m128d two = _mm_set1_pd(2);
	__m128d pfpmin = _mm_set1_pd(fpmin);
	__m128d eps = _mm_set1_pd(depsilon);
	__m128d b = _mm_add_pd(t, _mm_set1_pd(1-a));
	__m128d c = _mm_set1_pd(ifpmin);
	__m128d d = _mm_div_pd(one, b);
	__m128d h = _mm_mul_pd(_mm_set1_pd(0.5), d);
	__m128d an, del, mask;

 	for (int i = 1; i <= maxiter; i++)
 	{
		an = _mm_set1_pd(-i * (i-a));
		b = _mm_add_pd(b, two);
		d = _mm_add_pd(_mm_mul_pd(an, d), b);
		mask = _mm_cmplt_pd(qabs(d), pfpmin);
		d = select(mask, pfpmin, d);
		c = _mm_add_pd(b, _mm_div_pd(an, c));
		mask = _mm_cmplt_pd(qabs(c), pfpmin);
		c = select(mask, pfpmin, c);
		d = _mm_div_pd(one, d);
		del = _mm_mul_pd(d, c);
		h = _mm_mul_pd(h, del);
		mask = _mm_cmplt_pd(qabs(_mm_sub_pd(del, one)), eps);
		if (_mm_movemask_pd(mask) == 3)
		{
			__m128d res = _mm_set1_pd(lgamma_cache(int(a)));
			res = _mm_sub_pd(res, _mm_mul_pd(_mm_set1_pd(a), qlog(t)));
			res = _mm_mul_pd(_mm_set1_pd(0.5), qexp(res));
			res = _mm_sub_pd(res, _mm_mul_pd(h, expmt));
			return res;
		}
	}

	throw NoConvergence();
}

__m128d Fm_notzero(int m, __m128d t, __m128d expmt)
{
	if (m == 0)
	{
		__m128d st = _mm_sqrt_pd(t);
		return _mm_div_pd(qerf(st), _mm_mul_pd(st, _mm_set1_pd(M_2_SQRTPI)));
	}

	double a = m + 0.5;
	__m128d mask = _mm_cmpge_pd(t, _mm_set1_pd(std::max(20.0, 0.5*a)));
	int mbits = _mm_movemask_pd(mask);
	if (mbits == 0)
		return Fm_P_series_pos_a_t(a, t, expmt);
	if (mbits == 3)
		return Fm_Q_cf(a, t, expmt);

	__m128d tl = _mm_unpacklo_pd(t, t);
	__m128d el = _mm_unpacklo_pd(expmt, expmt);
	__m128d th = _mm_unpackhi_pd(t, t);
	__m128d eh = _mm_unpackhi_pd(expmt, expmt);
	__m128d fp, fq;
	if (mbits == 1)
	{
		fp = Fm_P_series_pos_a_t(a, th, eh);
		fq = Fm_Q_cf(a, tl, el);
	}
	else
	{
		fp = Fm_P_series_pos_a_t(a, tl, el);
		fq = Fm_Q_cf(a, th, eh);
	}
	return select(mask, fq, fp);
}

} // namespace

__m128d Fm(int m, __m128d t, __m128d expmt)
{
	__m128d mask = _mm_cmplt_pd(qabs(t), _mm_set1_pd(depsilon));
	int mbits = _mm_movemask_pd(mask);
	if (mbits == 0)
		return Fm_notzero(m, t, expmt);
	if (mbits == 3)
		return _mm_set1_pd(0.5/(m+0.5));

	t = mbits == 1 ? _mm_unpackhi_pd(t, t) : _mm_unpacklo_pd(t, t);
	return select(mask, _mm_set1_pd(0.5/(m+0.5)), Fm_notzero(m, t, expmt));
}
#endif // __SSE2__
