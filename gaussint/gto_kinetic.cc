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
