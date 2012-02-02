#ifndef GTO_ONE_ELEC_HH
#define GTO_ONE_ELEC_HH

#include "Eigen/Core"
#include "gaussint/gto_overlap.hh"
#include "gaussint/gto_kinetic.hh"

/*!
 * \brief Compute one-electron integrals between primitive GTOs
 *
 * Compute the overlap and kinetic energy integrals between two sets of
 * primitive Gaussians in a contracted GTO, without normalization.
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
void gto_one_elec_primitive_specialized(
	const Eigen::ArrayXXd& alpha, const Eigen::ArrayXXd& beta,
	const Eigen::ArrayXXd& asum, const Eigen::ArrayXXd& ared,
	const Eigen::ArrayXXd& exp_ared, const Eigen::Vector3d& r,
	Eigen::ArrayXXd& Sp, Eigen::ArrayXXd& Tp)
{
	Sp = gto_overlap_primitive_specialized<lxA, lyA, lzA, lxB, lyB, lzB>(
		alpha, beta, asum, ared, exp_ared, r);
	Tp = gto_kinetic_primitive_specialized<lxA, lyA, lzA, lxB, lyB, lzB>(
		alpha, beta, asum, ared, exp_ared, r);
}

/*!
 * \brief Compute the one electron integrals between contracted GTOs
 *
 * Compute the overlap and kinetic energy integral between two contracted
 * Gaussian type orbitals. This template function only works for combinations
 * of angular momentum quantum number for which specialized_primitive_overlap()
 * and specialized_kinetic_overlap() are implemented. For a generic
 * one-electron function, see gto_one_elec_generic().
 * \tparam lxA     Angular momentum in the \f$x\f$ direction of the first function
 * \tparam lyA     Angular momentum in the \f$y\f$ direction of the first function
 * \tparam lzA     Angular momentum in the \f$z\f$ direction of the first function
 * \tparam lxB     Angular momentum in the \f$x\f$ direction of the second function
 * \tparam lyB     Angular momentum in the \f$y\f$ direction of the second function
 * \tparam lzB     Angular momentum in the \f$z\f$ direction of the second function
 * \param weightsA Weights of the primitives in the first orbital
 * \param widthsA  Widths of the primitives in the first orbital
 * \param posA     Center of the first orbital
 * \param weightsB Weights of the primitives in the second orbital
 * \param widthsB  Widths of the primitives in the second orbital
 * \param posB     Center of the second orbital
 * \param S        Place to store the overlap
 * \param T        Place to store the kinetic energy
 */
template <int lxA, int lyA, int lzA, int lxB, int lyB, int lzB>
void gto_one_elec_specialized(const Eigen::VectorXd& weightsA,
        const Eigen::VectorXd& widthsA, const Eigen::Vector3d& posA,
        const Eigen::VectorXd& weightsB, const Eigen::VectorXd& widthsB,
        const Eigen::Vector3d& posB,
	double *S, double *T)
{
        int nA = widthsA.size(), nB = widthsB.size();
        Eigen::ArrayXXd alpha = widthsA.replicate(1, nB);
        Eigen::ArrayXXd beta = widthsB.transpose().replicate(nA, 1);
        Eigen::ArrayXXd asum = alpha + beta;
        Eigen::ArrayXXd ared = (widthsA * widthsB.transpose()).array() / asum;
        Eigen::Vector3d r = posB - posA;
        Eigen::ArrayXXd exp_ared = (-r.squaredNorm()*ared).exp();
	Eigen::ArrayXXd Sp, Tp;
	
	gto_one_elec_primitive_specialized<lxA, lyA, lzA, lxB, lyB, lzB>(
		alpha, beta, asum, ared, exp_ared, r, Sp, Tp);

	*S = Constants::pi_sqrt_pi * weightsA.transpose() * Sp.matrix() * weightsB;
	*T = Constants::pi_sqrt_pi * weightsA.transpose() * Tp.matrix() * weightsB;
}

#endif // GTO_ONE_ELEC_HH
