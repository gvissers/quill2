#ifndef GTO_OVERLAP_HH
#define GTO_OVERLAP_HH

/*!
 * \file gto_overlap.hh
 * \brief Functions for computing the overlap between contracted GTOs
 */

#include <Eigen/Dense>
#include "constants.hh"

/*!
 * \brief Compute the overlap between primitives
 *
 * Compute the overlap between the primitive Gaussians in a contracted GTO,
 * minus a normalization factor.
 * \param ls1      Angular momentum quantum numbers of the first orbital
 * \param ls2      Angular momentum quantum numbers of the second orbital
 * \param alpha    Primitive widths in the first orbital, for all second primitives
 * \param beta     Primitive widths in the second orbital, for all first primitives
 * \param asum     Sum of \a alpha and \a beta
 * \param exp_ared \f$\exp(-\xi r^2) with \f$r\$ the distance between the orbital centers
 * \param r        Vector from first to second orbital center
 */
Eigen::ArrayXXd gto_overlap_primitive_generic(
	const Eigen::Vector3i& ls1, const Eigen::Vector3i& ls2,
	const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd& beta,
	const Eigen::ArrayXXd& asum,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r);

/*!
 * \brief Compute the overlap between primitives
 *
 * Specializations of this template compute the overlap between the primitive
 * Gaussians in a contracted GTO, minus a normalization factor.
 * \tparam lxA     Angular momentum in the \f$x\f$ direction of the first function
 * \tparam lyA     Angular momentum in the \f$y\f$ direction of the first function
 * \tparam lzA     Angular momentum in the \f$z\f$ direction of the first function
 * \tparam lxB     Angular momentum in the \f$x\f$ direction of the second function
 * \tparam lyB     Angular momentum in the \f$y\f$ direction of the second function
 * \tparam lzB     Angular momentum in the \f$z\f$ direction of the second function
 * \param alpha    Primitive widths in the first orbital, for all second primitives
 * \param beta     Primitive widths in the second orbital, for all first primitives
 * \param asum     Sum of \a alpha and \a beta
 * \param ared     "Reduced" widths \f$\xi = \alpha\beta / (\alpha+\beta)\f$
 * \param exp_ared \f$\exp(-\xi r^2) with \f$r\$ the distance between the orbital centers
 * \param r        Vector from first to second orbital center
 */
template <int lxA, int lyA, int lzA, int lxB, int lyB, int lzB>
Eigen::ArrayXXd gto_overlap_primitive_specialized(
	const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd& beta,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r)
{
	return gto_overlap_primitive_generic(
		Eigen::Vector3i(lxA, lyA, lzA), Eigen::Vector3i(lxB, lyB, lzB),
		alpha, beta, asum, exp_ared, r);
}
template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 0, 0, 0, 0, 0>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d&);
template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<1, 0, 0, 0, 0, 0>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd& beta,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r);
template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 1, 0, 0, 0, 0>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd& beta,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r);
template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 0, 1, 0, 0, 0>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd& beta,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r);
template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 0, 0, 1, 0, 0>(
	const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r);
template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 0, 0, 0, 1, 0>(
	const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r);
template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 0, 0, 0, 0, 1>(
	const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r);
template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<2, 0, 0, 0, 0, 0>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd& beta,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r);
template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<1, 1, 0, 0, 0, 0>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd& beta,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r);
template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<1, 0, 1, 0, 0, 0>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd& beta,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r);
template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<1, 0, 0, 1, 0, 0>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r);
template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<1, 0, 0, 0, 1, 0>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r);
template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<1, 0, 0, 0, 0, 1>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r);
template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 2, 0, 0, 0, 0>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd& beta,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r);
template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 1, 1, 0, 0, 0>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd& beta,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r);
template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 1, 0, 1, 0, 0>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r);
template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 1, 0, 0, 1, 0>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r);
template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 1, 0, 0, 0, 1>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r);
template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 0, 2, 0, 0, 0>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd& beta,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r);
template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 0, 1, 1, 0, 0>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r);
template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 0, 1, 0, 1, 0>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r);
template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 0, 1, 0, 0, 1>(
	const Eigen::ArrayXXd&, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r);
template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 0, 0, 2, 0, 0>(
	const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r);
template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 0, 0, 1, 1, 0>(
	const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r);
template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 0, 0, 1, 0, 1>(
	const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r);
template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 0, 0, 0, 2, 0>(
	const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r);
template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 0, 0, 0, 1, 1>(
	const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r);
template <>
Eigen::ArrayXXd gto_overlap_primitive_specialized<0, 0, 0, 0, 0, 2>(
	const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd&,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r);

/*!
 * \brief Compute the overlap between contracted GTOs
 *
 * Compute the overlap between two contracted Gaussian type orbitals. This
 * template function only works for combinations of angular momentum
 * quantum number for which gto_overlap_primitive_specialized() is implemented.
 * For a generic overlap function, see gto_overlap_generic().
 * \tparam lx1     Angular momentum in  the \f$x\f$ direction of the first function
 * \tparam ly1     Angular momentum in  the \f$y\f$ direction of the first function
 * \tparam lz1     Angular momentum in  the \f$z\f$ direction of the first function
 * \tparam lx2     Angular momentum in  the \f$x\f$ direction of the second function
 * \tparam ly2     Angular momentum in  the \f$y\f$ direction of the second function
 * \tparam lz2     Angular momentum in  the \f$z\f$ direction of the second function
 * \param weights1 Weights of the primitives in the first orbital
 * \param widths1  Widths of the primitives in the first orbital
 * \param pos1     Center of the first orbital
 * \param weights2 Weights of the primitives in the second orbital
 * \param widths2  Widths of the primitives in the second orbital
 * \param pos2     Center of the second orbital
 * \return The overlap between the two orbitals
 */
template <int lx1, int ly1, int lz1, int lx2, int ly2, int lz2>
double gto_overlap_specialized(const Eigen::VectorXd& weights1,
	const Eigen::VectorXd& widths1, const Eigen::Vector3d& pos1,
	const Eigen::VectorXd& weights2, const Eigen::VectorXd& widths2,
	const Eigen::Vector3d& pos2)
{
	int n1 = widths1.size(), n2 = widths2.size();
	Eigen::ArrayXXd alpha = widths1.replicate(1, n2);
	Eigen::ArrayXXd beta = widths2.transpose().replicate(n1, 1);
	Eigen::ArrayXXd asum = alpha + beta;
        Eigen::ArrayXXd ared = (widths1 * widths2.transpose()).array() / asum;
	Eigen::Vector3d r = pos2 - pos1;
	Eigen::ArrayXXd exp_ared = (-r.squaredNorm()*ared).exp();

	return Constants::pi_sqrt_pi * weights1.transpose()
		* gto_overlap_primitive_specialized<lx1, ly1, lz1, lx2, ly2, lz2>(
			alpha, beta, asum, ared, exp_ared, r).matrix()
		* weights2;
}

#endif // GTO_OVERLAP_HH
