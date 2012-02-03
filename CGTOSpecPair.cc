#include "CGTOSpecPair.hh"
#include "boys.hh"

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

// S_00 = exp_ared / (asum * asum.sqrt());
// S_10 = x beta S_00 / asum
// T_10 = (x_p-x_a)T_00 + 2 ared S_10
//      = x beta T_00 / asum + 2 ared S_10
//      = x beta T_00 / asum + 2 x beta ared S_00 / asum
//      = [x beta ared(3 - 2 r^2 ared) / asum + 2 x beta ared / asum] S_00
//      = x beta ared(5 - 2 r^2 ared) S_00 / asum
static Eigen::ArrayXXd kinetic_ps(double x, double rsq,
	const Eigen::ArrayXXd& beta, const Eigen::ArrayXXd& asum,
	const Eigen::ArrayXXd& ared, const Eigen::ArrayXXd& exp_ared)
{
	return x * beta * ared * (5 - 2*rsq*ared) * exp_ared
		/ (asum.square() * asum.sqrt());
}

// d/dx G_0 = -2 alpha G_1
// <0| d^2/dx^2 |0>
// = -4 alpha beta S_11
// = -4 alpha beta ((x_p-x_a)S_01 + S_00/2asum)
// = -4 alpha beta ((x_p-x_a)(x_p-xb) + 1/2asum) S_00
// = -4 alpha beta (beta(x_b-x_a)alpha(x_a-x_b)/asum^2 + 1/2asum) S_00
// = xi (4 xi (xb-xa)^2 - 2) S_00
// T_00 = -1/2 <0| d^2/dx^2+d^2/dy^2 + d^2/dz^2 |0>
//      = xi (3 - 2 xi r^2) S_00
template <>
double CGTOSpecPair<0, 0, 0, 0, 0, 0>::kineticEnergy() const
{
	return Constants::pi_sqrt_pi * f().weights().transpose()
		* (widthsReduced() * (3 - 2*r().squaredNorm()*widthsReduced())
			* exp_ared() / (widthsSum() * widthsSum().sqrt())).matrix()
		* g().weights();
}

template <>
double CGTOSpecPair<1, 0, 0, 0, 0, 0>::kineticEnergy() const
{
	return Constants::pi_sqrt_pi * f().weights().transpose()
		* kinetic_ps(r(0), r().squaredNorm(), widthsB(), widthsSum(),
			widthsReduced(), exp_ared()).matrix()
		* g().weights();
}

template <>
double CGTOSpecPair<0, 1, 0, 0, 0, 0>::kineticEnergy() const
{
	return Constants::pi_sqrt_pi * f().weights().transpose()
		* kinetic_ps(r(1), r().squaredNorm(), widthsB(), widthsSum(),
			widthsReduced(), exp_ared()).matrix()
		* g().weights();
}

template <>
double CGTOSpecPair<0, 0, 1, 0, 0, 0>::kineticEnergy() const
{
	return Constants::pi_sqrt_pi * f().weights().transpose()
		* kinetic_ps(r(2), r().squaredNorm(), widthsB(), widthsSum(),
			widthsReduced(), exp_ared()).matrix()
		* g().weights();
}

template <>
double CGTOSpecPair<0, 0, 0, 1, 0, 0>::kineticEnergy() const
{
	return Constants::pi_sqrt_pi * f().weights().transpose()
		* kinetic_ps(-r(0), r().squaredNorm(), widthsA(), widthsSum(),
			widthsReduced(), exp_ared()).matrix()
		* g().weights();
}

template <>
double CGTOSpecPair<0, 0, 0, 0, 1, 0>::kineticEnergy() const
{
	return Constants::pi_sqrt_pi * f().weights().transpose()
		* kinetic_ps(-r(1), r().squaredNorm(), widthsA(), widthsSum(),
			widthsReduced(), exp_ared()).matrix()
		* g().weights();
}

template <>
double CGTOSpecPair<0, 0, 0, 0, 0, 1>::kineticEnergy() const
{
	return Constants::pi_sqrt_pi * f().weights().transpose()
		* kinetic_ps(-r(2), r().squaredNorm(), widthsA(), widthsSum(),
			widthsReduced(), exp_ared()).matrix()
		* g().weights();
}

template<>
double CGTOSpecPair<0, 0, 0, 0, 0, 0>::nuclearAttraction(
	const Eigen::MatrixXd& nuc_pos, const Eigen::VectorXd& nuc_charge) const
{
	int nr_nuc = nuc_pos.cols();
	int nr_rows = f().size(), nr_cols = g().size();
	Eigen::ArrayXXd Pi;
	Eigen::ArrayXXd U = Eigen::ArrayXXd::Zero(nr_rows, nr_cols*nr_nuc);

	for (int i = 0; i < 3; i++)
	{
		Pi = P(i);
		for (int iC = 0; iC < nr_nuc; ++iC)
			U.block(0, iC*nr_cols, nr_rows, nr_cols) += (Pi - nuc_pos(i, iC)).square();
	}
	U *= widthsSum().replicate(1, nr_nuc);

	Eigen::ArrayXXd Am = Fm(0, U);
	Eigen::ArrayXXd A = Eigen::ArrayXXd::Zero(nr_rows, nr_cols);
	for (int iC = 0; iC < nr_nuc; ++iC)
		A += -nuc_charge[iC] * Am.block(0, iC*nr_cols, nr_rows, nr_cols);

	return 2 * M_PI * f().weights().transpose()
		* (exp_ared() * A / widthsSum()).matrix()
		* g().weights();
}

