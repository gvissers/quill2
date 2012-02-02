#include <vector>
#include <Exception.hh>
#include "gaussint/gto_overlap.hh"

template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 0, 0, 0, 0, 0>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d&)
{
	return exp_ared / (asum * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<1, 0, 0, 0, 0, 0>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd& beta,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return r.x() * beta * exp_ared / (asum.square() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 1, 0, 0, 0, 0>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd& beta,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return r.y() * beta * exp_ared / (asum.square() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 0, 1, 0, 0, 0>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd& beta,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return r.z() * beta * exp_ared / (asum.square() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 0, 0, 1, 0, 0>(
	const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return -r.x() * alpha * exp_ared / (asum.square() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 0, 0, 0, 1, 0>(
	const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return -r.y() * alpha * exp_ared / (asum.square() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 0, 0, 0, 0, 1>(
	const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return -r.z() * alpha * exp_ared / (asum.square() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<2, 0, 0, 0, 0, 0>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd& beta,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return (0.5*asum + r.x()*r.x()*beta.square()) * exp_ared
		/ (asum.cube() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<1, 1, 0, 0, 0, 0>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd& beta,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return r.x() * r.y() * beta.square() * exp_ared
		/ (asum.cube() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<1, 0, 1, 0, 0, 0>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd& beta,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return r.x() * r.z() * beta.square() * exp_ared
		/ (asum.cube() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<1, 0, 0, 1, 0, 0>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return (0.5 - r.x()*r.x()*ared) * exp_ared
		/ (asum.square() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<1, 0, 0, 0, 1, 0>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return -r.x() * r.y() * ared * exp_ared
		/ (asum.square() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<1, 0, 0, 0, 0, 1>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return -r.x() * r.z() * ared * exp_ared
		/ (asum.square() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 2, 0, 0, 0, 0>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd& beta,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return (0.5*asum + r.y()*r.y()*beta.square()) * exp_ared
		/ (asum.cube() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 1, 1, 0, 0, 0>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd& beta,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return r.y() * r.z() * beta.square() * exp_ared
		/ (asum.cube() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 1, 0, 1, 0, 0>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return -r.y() * r.x() * ared * exp_ared
		/ (asum.square() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 1, 0, 0, 1, 0>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return (0.5 - r.y()*r.y()*ared) * exp_ared
		/ (asum.square() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 1, 0, 0, 0, 1>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return -r.y() * r.z() * ared * exp_ared
		/ (asum.square() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 0, 2, 0, 0, 0>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd& beta,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return (0.5*asum + r.z()*r.z()*beta.square()) * exp_ared
		/ (asum.cube() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 0, 1, 1, 0, 0>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return -r.z() * r.x() * ared * exp_ared
		/ (asum.square() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 0, 1, 0, 1, 0>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return -r.z() * r.y() * ared * exp_ared
		/ (asum.square() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 0, 1, 0, 0, 1>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return (0.5 - r.z()*r.z()*ared) * exp_ared
		/ (asum.square() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 0, 0, 2, 0, 0>(
	const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return (0.5*asum + r.x()*r.x()*alpha.square()) * exp_ared
		/ (asum.cube() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 0, 0, 1, 1, 0>(
	const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return r.x() * r.y() * alpha.square() * exp_ared
		/ (asum.cube() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 0, 0, 1, 0, 1>(
	const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return r.x() * r.z() * alpha.square() * exp_ared
		/ (asum.cube() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 0, 0, 0, 2, 0>(
	const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return (0.5*asum + r.y()*r.y()*alpha.square()) * exp_ared
		/ (asum.cube() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 0, 0, 0, 1, 1>(
	const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return r.y() * r.z() * alpha.square() * exp_ared
		/ (asum.cube() * asum.sqrt());
}

template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 0, 0, 0, 0, 2>(
	const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return (0.5*asum + r.z()*r.z()*alpha.square()) * exp_ared
		/ (asum.cube() * asum.sqrt());
}
