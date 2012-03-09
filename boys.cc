#include <limits>
#include <cmath>
#include "boys.hh"
#include "exceptions.hh"

static const int maxiter = 200;
static const double fpmin = std::numeric_limits<double>::min();
static const double ifpmin = 1.0 / std::numeric_limits<double>::min();
static const double depsilon = std::numeric_limits<double>::epsilon();

static double Fm_P_series_pos_a_t(double a, double t, double expmt)
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

static double Fm_Q_cf(double a, double t, double expmt)
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
			return 0.5 * qexp(std::lgamma(a) - a*std::log(t))
				- h * expmt;
	}

	throw NoConvergence();
}

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
