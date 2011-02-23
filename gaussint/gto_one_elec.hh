#ifndef GTO_ONE_ELEC_HH
#define GTO_ONE_ELEC_HH

#include "Eigen/Core"
#include "gaussint/gto_overlap.hh"
#include "gaussint/gto_kinetic.hh"

/*!
 * \brief Compute the one electron integrals between contracted GTOs
 *
 * Compute the overlap and kinetic energy integral between two contracted
 * Gaussian type orbitals. This template function only works for combinations
 * of angular momentum quantum number for which specialized_primitive_overlap()
 * and specialized_kinetic_overlap() are implemented. For a generic
 * one-electron function, see gto_one_elec_generic().
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
 * \param S        Place to store the overlap
 * \param T        Place to store the kinetic energy
 */
template <int lx1, int ly1, int lz1, int lx2, int ly2, int lz2>
void gto_one_elec_specialized(const Eigen::VectorXd& weights1,
        const Eigen::VectorXd& widths1, const Eigen::Vector3d& pos1,
        const Eigen::VectorXd& weights2, const Eigen::VectorXd& widths2,
        const Eigen::Vector3d& pos2,
	double *S, double *T)
{
        int n1 = widths1.size(), n2 = widths2.size();
        Eigen::ArrayXXd alpha = widths1.replicate(1, n2);
        Eigen::ArrayXXd beta = widths2.transpose().replicate(n1, 1);
        Eigen::ArrayXXd asum = alpha + beta;
        Eigen::ArrayXXd ared = (widths1 * widths2.transpose()).array() / asum;
        Eigen::Vector3d r = pos2 - pos1;
        Eigen::ArrayXXd exp_ared = (-r.squaredNorm()*ared).exp();

	*S = Constants::pi_sqrt_pi * weights1.transpose()
                * gto_overlap_primitive_specialized<lx1, ly1, lz1, lx2, ly2, lz2>(
                        alpha, beta, asum, ared, exp_ared, r).matrix()
                * weights2;
	*T = Constants::pi_sqrt_pi * weights1.transpose()
                * gto_kinetic_primitive_specialized<lx1, ly1, lz1, lx2, ly2, lz2>(
                        alpha, beta, asum, ared, exp_ared, r).matrix()
                * weights2;
}

/*!
 * \brief Compute one-electron integrals between primtive GTOs
 *
 * Compute the overlap and kinetic energy integrals between two sets of
 * primitive Gaussians in a contracted GTO, without normalization.
 * \param ls1      Angular momentum quantum numbers of the first orbital
 * \param ls2      Angular momentum quantum numbers of the second orbital
 * \param alpha    Primitive widths in the first orbital, for all second primitives
 * \param beta     Primitive widths in the second orbital, for all first primitives
 * \param asum     Sum of \a alpha and \a beta
 * \param ared     "Reduced" widths \f$\xi = \alpha\beta / (\alpha+\beta)\f$
 * \param exp_ared \f$\exp(-\xi r^2) with \f$r\$ the distance between the orbital centers
 * \param r        Vector from first to second orbital center
 */
void gto_one_elec_primitive_generic(
	const Eigen::Vector3i& ls1, const Eigen::Vector3i& ls2,
	const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd& beta,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r,
	Eigen::ArrayXXd& Sp, Eigen::ArrayXXd& Tp);

/*!
 * \brief Compute one-electron integrals between contracted GTOs
 *
 * Compute the overlap and kinetic energy integrals between two contracted
 * Gaussian type orbitals.
 * \param ls1      Angular momentum quantum numbers of the first orbital
 * \param weights1 Weights of the primitives in the first orbital
 * \param widths1  Widths of the primitives in the first orbital
 * \param pos1     Center of the first orbital
 * \param ls2      Angular momentum quantum numbers of the second orbital
 * \param weights2 Weights of the primitives in the second orbital
 * \param widths2  Widths of the primitives in the second orbital
 * \param pos2     Center of the second orbital
 * \param S        Place to store the overlap
 * \param T        Place to store the kinetic energy
 */
void gto_one_elec_generic(const Eigen::Vector3i& ls1,
	const Eigen::VectorXd& weights1, const Eigen::VectorXd& widths1,
	const Eigen::Vector3d& pos1,
	const Eigen::Vector3i& ls2,
	const Eigen::VectorXd& weights2, const Eigen::VectorXd& widths2,
	const Eigen::Vector3d& pos2,
	double *S, double *T);

#endif // GTO_ONE_ELEC_HH
