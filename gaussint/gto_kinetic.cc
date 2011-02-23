#include <vector>
#include <Exception.hh>
#include "gaussint/gto_kinetic.hh"
#include "gaussint/gto_one_elec.hh"

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
Eigen::ArrayXXd gto_kinetic_primitive_specialized<0, 0, 0, 0, 0, 0>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return ared * (3 - 2*r.squaredNorm()*ared) * exp_ared
		/ (asum * asum.sqrt());
}

// S_00 = exp_ared / (asum * asum.sqrt());
// S_10 = x beta S_00 / asum
// T_10 = (x_p-x_a)T_00 + 2 ared S_10
//      = x beta T_00 / asum + 2 ared S_10
//      = x beta T_00 / asum + 2 x beta ared S_00 / asum
//      = [x beta ared(3 - 2 r^2 ared) / asum + 2 x beta ared / asum] S_00
//      = x beta ared(5 - 2 r^2 ared) S_00 / asum 

template <>
Eigen::ArrayXXd gto_kinetic_primitive_specialized<1, 0, 0, 0, 0, 0>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd& beta,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return r.x() * beta * ared * (5 - 2*r.squaredNorm()*ared) * exp_ared
		/ (asum.square() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_kinetic_primitive_specialized<0, 1, 0, 0, 0, 0>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd& beta,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return r.y() * beta * ared * (5 - 2*r.squaredNorm()*ared) * exp_ared
		/ (asum.square() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_kinetic_primitive_specialized<0, 0, 1, 0, 0, 0>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd& beta,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return r.z() * beta * ared * (5 - 2*r.squaredNorm()*ared) * exp_ared
		/ (asum.square() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_kinetic_primitive_specialized<0, 0, 0, 1, 0, 0>(
	const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return -r.x() * alpha * ared * (5 - 2*r.squaredNorm()*ared) * exp_ared
		/ (asum.square() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_kinetic_primitive_specialized<0, 0, 0, 0, 1, 0>(
	const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return -r.y() * alpha * ared * (5 - 2*r.squaredNorm()*ared) * exp_ared
		/ (asum.square() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_kinetic_primitive_specialized<0, 0, 0, 0, 0, 1>(
	const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return -r.z() * alpha * ared * (5 - 2*r.squaredNorm()*ared) * exp_ared
		/ (asum.square() * asum.sqrt());
}

Eigen::ArrayXXd gto_kinetic_primitive_generic(
	const Eigen::Vector3i& ls1, const Eigen::Vector3i& ls2,
	const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd& beta,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	// We need to compute overlaps to do the recursion for the kinetic
	// energy, so compute both and only return the kinetic energy.
	Eigen::ArrayXXd Sp, Tp;
	gto_one_elec_primitive_generic(ls1, ls2, alpha, beta, asum, ared,
		exp_ared, r, Sp, Tp);
	return Tp;
}

double gto_kinetic_generic(const Eigen::Vector3i& ls1,
	const Eigen::VectorXd& weights1, const Eigen::VectorXd& widths1,
	const Eigen::Vector3d& pos1,
	const Eigen::Vector3i& ls2,
	const Eigen::VectorXd& weights2, const Eigen::VectorXd& widths2,
	const Eigen::Vector3d& pos2)
{
	// We need to compute overlaps to do the recursion for the kinetic
	// energy, so compute both and only return the kinetic energy.
	double S, T;
	gto_one_elec_generic(ls1, weights1, widths1, pos1, ls2, weights2, widths2, pos2, &S, &T);
	return T;
}
