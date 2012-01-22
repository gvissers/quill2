#ifndef GTO_NUC_ATTR_HH
#define GTO_NUC_ATTR_HH

#include "Eigen/Dense"

/*!
 * \brief Compute nuclear attraction integrals between primitive GTOs
 *
 * Compute the nuclear attraction integrals between two sets of
 * primitive Gaussians in a contracted GTO, without normalization.
 * \param lsA       Angular momentum quantum numbers of the first orbital
 * \param lsB       Angular momentum quantum numbers of the second orbital
 * \param alpha     Primitive widths in the first orbital, for all second primitives
 * \param beta      Primitive widths in the second orbital, for all first primitives
 * \param posA      The center of the first orbital
 * \param posB      The center of the second orbital
 * \param asum      Sum of \a alpha and \a beta
 * \param exp_ared  \f$\exp(-\xi r^2) with \f$r\$ the distance between the orbital centers
 * \param nuc_pos   The positions of the nuclei
 * \param nuc_charge The charges of the nuclei
 */
Eigen::ArrayXXd gto_nuc_attr_primitive_generic(
	const Eigen::Vector3i& lsA, const Eigen::Vector3i& lsB,
	const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd& beta,
	const Eigen::Vector3d& posA, const Eigen::Vector3d& posB,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& exp_ared,
	const Eigen::MatrixXd& nuc_pos, const Eigen::VectorXd& nuc_charge);

template <int lxA, int lyA, int lzA, int lxB, int lyB, int lzB>
Eigen::ArrayXXd gto_nuc_attr_primitive_specialized(
	const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd& beta,
	const Eigen::Vector3d& posA, const Eigen::Vector3d& posB,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& exp_ared,
	const Eigen::MatrixXd& nuc_pos, const Eigen::VectorXd& nuc_charge)
{
	return gto_nuc_attr_primitive_generic(
		Eigen::Vector3i(lxA, lyA, lzA), Eigen::Vector3i(lxB, lyB, lzB),
		alpha, beta, posA, posB, asum, exp_ared, nuc_pos, nuc_charge);
}
template <>
Eigen::ArrayXXd gto_nuc_attr_primitive_specialized<0, 0, 0, 0, 0, 0>(
	const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd& beta,
	const Eigen::Vector3d& posA, const Eigen::Vector3d& posB,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& exp_ared,
	const Eigen::MatrixXd& nuc_pos, const Eigen::VectorXd& nuc_charge);

/*!
 * \brief Compute nuclear attraction integrals between primitive GTOs
 *
 * Compute the nuclear attraction integrals between contracted Gaussian type 
 * orbitals.
 * \param ls1       Angular momentum quantum numbers of the first orbital
 * \param weights1  Weights of the primitives in the first orbital
 * \param widths1   Widths of the primitives in the first orbital
 * \param pos1      The center of the first orbital
 * \param ls2       Angular momentum quantum numbers of the second orbital
 * \param weights2  Weights of the primitives in the second orbital
 * \param widths2   Widths of the primitives in the second orbital
 * \param pos2      The center of the second orbital
 * \param nuc_pos   The positions of the nuclei
 * \param nuc_charge The charges of the nuclei
 */
double gto_nuc_attr_primitive_generic(const Eigen::Vector3i& ls1,
	const Eigen::VectorXd& weights1, const Eigen::VectorXd& widths1,
	const Eigen::Vector3d& pos1,
	const Eigen::Vector3i& ls2,
	const Eigen::VectorXd& weights2, const Eigen::VectorXd& widths2,
	const Eigen::Vector3d& pos2,
	const Eigen::MatrixXd& nuc_pos, const Eigen::VectorXd& nuc_charge);

#endif // GTO_NUC_ATTR_HH
