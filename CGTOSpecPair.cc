#include "CGTOSpecPair.hh"

static Eigen::ArrayXXd overlap_ps(double x, const Eigen::ArrayXXd& beta,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& exp_ared)
{
	return x * beta * exp_ared / (asum.square() * asum.sqrt());
}

static Eigen::ArrayXXd overlap_pp(double x, const Eigen::ArrayXXd& asum,
	const Eigen::ArrayXXd& ared, const Eigen::ArrayXXd& exp_ared)
{
	return (0.5 - x*x*ared) * exp_ared / (asum.square() * asum.sqrt());
}

static Eigen::ArrayXXd overlap_pp(double x, double y, const Eigen::ArrayXXd& asum,
	const Eigen::ArrayXXd& ared, const Eigen::ArrayXXd& exp_ared)
{
	return -x * y * ared * exp_ared / (asum.square() * asum.sqrt());
}

static Eigen::ArrayXXd overlap_ds(double x, const Eigen::ArrayXXd& beta,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& exp_ared)
{
	return (0.5*asum + x*x*beta.square()) * exp_ared
		/ (asum.cube() * asum.sqrt());
}

static Eigen::ArrayXXd overlap_ds(double x, double y, const Eigen::ArrayXXd& beta,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& exp_ared)
{
	return x * y * beta.square() * exp_ared / (asum.cube() * asum.sqrt());
}

template <>
double CGTOSpecPair<0, 0, 0, 0, 0, 0>::overlap() const
{
	return Constants::pi_sqrt_pi * f().weights().transpose()
		* (exp_ared() / (widthsSum() * widthsSum().sqrt())).matrix()
		* g().weights();
}

template <>
double CGTOSpecPair<1, 0, 0, 0, 0, 0>::overlap() const
{
	return Constants::pi_sqrt_pi * f().weights().transpose()
		* overlap_ps(r(0), widthsB(), widthsSum(), exp_ared()).matrix()
		* g().weights();
}

template <>
double CGTOSpecPair<0, 1, 0, 0, 0, 0>::overlap() const
{
	return Constants::pi_sqrt_pi * f().weights().transpose()
		* overlap_ps(r(1), widthsB(), widthsSum(), exp_ared()).matrix()
		* g().weights();
}

template <>
double CGTOSpecPair<0, 0, 1, 0, 0, 0>::overlap() const
{
	return Constants::pi_sqrt_pi * f().weights().transpose()
		* overlap_ps(r(2), widthsB(), widthsSum(), exp_ared()).matrix()
		* g().weights();
}

template <>
double CGTOSpecPair<0, 0, 0, 1, 0, 0>::overlap() const
{
	return Constants::pi_sqrt_pi * f().weights().transpose()
		* overlap_ps(-r(0), widthsA(), widthsSum(), exp_ared()).matrix()
		* g().weights();
}

template <>
double CGTOSpecPair<0, 0, 0, 0, 1, 0>::overlap() const
{
	return Constants::pi_sqrt_pi * f().weights().transpose()
		* overlap_ps(-r(1), widthsA(), widthsSum(), exp_ared()).matrix()
		* g().weights();
}

template <>
double CGTOSpecPair<0, 0, 0, 0, 0, 1>::overlap() const
{
	return Constants::pi_sqrt_pi * f().weights().transpose()
		* overlap_ps(-r(2), widthsA(), widthsSum(), exp_ared()).matrix()
		* g().weights();
}

template <>
double CGTOSpecPair<2, 0, 0, 0, 0, 0>::overlap() const
{
	return Constants::pi_sqrt_pi * f().weights().transpose()
		* overlap_ds(r(0), widthsB(), widthsSum(), exp_ared()).matrix()
		* g().weights();
}

template <>
double CGTOSpecPair<1, 1, 0, 0, 0, 0>::overlap() const
{
	return Constants::pi_sqrt_pi * f().weights().transpose()
		* overlap_ds(r(0), r(1), widthsB(), widthsSum(), exp_ared()).matrix()
		* g().weights();
}

template <>
double CGTOSpecPair<1, 0, 1, 0, 0, 0>::overlap() const
{
	return Constants::pi_sqrt_pi * f().weights().transpose()
		* overlap_ds(r(0), r(2), widthsB(), widthsSum(), exp_ared()).matrix()
		* g().weights();
}

template <>
double CGTOSpecPair<1, 0, 0, 1, 0, 0>::overlap() const
{
	return Constants::pi_sqrt_pi * f().weights().transpose()
		* overlap_pp(r(0), widthsSum(), widthsReduced(), exp_ared()).matrix()
		* g().weights();
}

template <>
double CGTOSpecPair<1, 0, 0, 0, 1, 0>::overlap() const
{
	return Constants::pi_sqrt_pi * f().weights().transpose()
		* overlap_pp(r(0), r(1), widthsSum(), widthsReduced(), exp_ared()).matrix()
		* g().weights();
}

template <>
double CGTOSpecPair<1, 0, 0, 0, 0, 1>::overlap() const
{
	return Constants::pi_sqrt_pi * f().weights().transpose()
		* overlap_pp(r(0), r(2), widthsSum(), widthsReduced(), exp_ared()).matrix()
		* g().weights();
}

template <>
double CGTOSpecPair<0, 2, 0, 0, 0, 0>::overlap() const
{
	return Constants::pi_sqrt_pi * f().weights().transpose()
		* overlap_ds(r(1), widthsB(), widthsSum(), exp_ared()).matrix()
		* g().weights();
}

template <>
double CGTOSpecPair<0, 1, 1, 0, 0, 0>::overlap() const
{
	return Constants::pi_sqrt_pi * f().weights().transpose()
		* overlap_ds(r(1), r(2), widthsB(), widthsSum(), exp_ared()).matrix()
		* g().weights();
}

template <>
double CGTOSpecPair<0, 1, 0, 1, 0, 0>::overlap() const
{
	return Constants::pi_sqrt_pi * f().weights().transpose()
		* overlap_pp(r(1), r(0), widthsSum(), widthsReduced(), exp_ared()).matrix()
		* g().weights();
}

template <>
double CGTOSpecPair<0, 1, 0, 0, 1, 0>::overlap() const
{
	return Constants::pi_sqrt_pi * f().weights().transpose()
		* overlap_pp(r(1), widthsSum(), widthsReduced(), exp_ared()).matrix()
		* g().weights();
}

template <>
double CGTOSpecPair<0, 1, 0, 0, 0, 1>::overlap() const
{
	return Constants::pi_sqrt_pi * f().weights().transpose()
		* overlap_pp(r(1), r(2), widthsSum(), widthsReduced(), exp_ared()).matrix()
		* g().weights();
}

template <>
double CGTOSpecPair<0, 0, 2, 0, 0, 0>::overlap() const
{
	return Constants::pi_sqrt_pi * f().weights().transpose()
		* overlap_ds(r(2), widthsB(), widthsSum(), exp_ared()).matrix()
		* g().weights();
}

template <>
double CGTOSpecPair<0, 0, 1, 1, 0, 0>::overlap() const
{
	return Constants::pi_sqrt_pi * f().weights().transpose()
		* overlap_pp(r(2), r(0), widthsSum(), widthsReduced(), exp_ared()).matrix()
		* g().weights();
}

template <>
double CGTOSpecPair<0, 0, 1, 0, 1, 0>::overlap() const
{
	return Constants::pi_sqrt_pi * f().weights().transpose()
		* overlap_pp(r(2), r(1), widthsSum(), widthsReduced(), exp_ared()).matrix()
		* g().weights();
}

template <>
double CGTOSpecPair<0, 0, 1, 0, 0, 1>::overlap() const
{
	return Constants::pi_sqrt_pi * f().weights().transpose()
		* overlap_pp(r(2), widthsSum(), widthsReduced(), exp_ared()).matrix()
		* g().weights();
}

template <>
double CGTOSpecPair<0, 0, 0, 2, 0, 0>::overlap() const
{
	return Constants::pi_sqrt_pi * f().weights().transpose()
		* overlap_ds(r(0), widthsA(), widthsSum(), exp_ared()).matrix()
		* g().weights();
}

template <>
double CGTOSpecPair<0, 0, 0, 1, 1, 0>::overlap() const
{
	return Constants::pi_sqrt_pi * f().weights().transpose()
		* overlap_ds(r(0), r(1), widthsA(), widthsSum(), exp_ared()).matrix()
		* g().weights();
}

template <>
double CGTOSpecPair<0, 0, 0, 1, 0, 1>::overlap() const
{
	return Constants::pi_sqrt_pi * f().weights().transpose()
		* overlap_ds(r(0), r(2), widthsA(), widthsSum(), exp_ared()).matrix()
		* g().weights();
}

template <>
double CGTOSpecPair<0, 0, 0, 0, 2, 0>::overlap() const
{
	return Constants::pi_sqrt_pi * f().weights().transpose()
		* overlap_ds(r(1), widthsA(), widthsSum(), exp_ared()).matrix()
		* g().weights();
}

template <>
double CGTOSpecPair<0, 0, 0, 0, 1, 1>::overlap() const
{
	return Constants::pi_sqrt_pi * f().weights().transpose()
		* overlap_ds(r(1), r(2), widthsA(), widthsSum(), exp_ared()).matrix()
		* g().weights();
}

template <>
double CGTOSpecPair<0, 0, 0, 0, 0, 2>::overlap() const
{
	return Constants::pi_sqrt_pi * f().weights().transpose()
		* overlap_ds(r(2), widthsA(), widthsSum(), exp_ared()).matrix()
		* g().weights();
}
