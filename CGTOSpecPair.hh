#ifndef CGTOSPECPAIR_HH
#define CGTOSPECPAIR_HH

/*!
 * \file CGTOSpecPair.hh
 * \brief Definition of the CGTOSpecPair class
 */

#include <Eigen/Dense>
#include "CGTOPair.hh"
#include "CGTOSpec.hh"
#include "gaussint/gto_one_elec.hh"
#include "constants.hh"

/*!
 * \brief Class for CGTO specializations
 *
 * Class CGTOSpecPair holds a pair of specializations of contracted Gaussian
 * type orbitals, for specific angular momentum quantum numbers.
 * \tparam lx1 Angular momentum in the \f$x\f$ direction for the first orbital
 * \tparam ly1 Angular momentum in the \f$y\f$ direction for the first orbital
 * \tparam lz1 Angular momentum in the \f$z\f$ direction for the first orbital
 * \tparam lx2 Angular momentum in the \f$x\f$ direction for the second orbital
 * \tparam ly2 Angular momentum in the \f$y\f$ direction for the seoncd orbital
 * \tparam lz2 Angular momentum in the \f$z\f$ direction for the second orbital
 */
template <int lx1, int ly1, int lz1, int lx2, int ly2, int lz2>
struct CGTOSpecPair: public CGTOPair
{
	//! Constructor
	CGTOSpecPair(const CGTOSpec<lx1, ly1, lz1>& f,
		const CGTOSpec<lx2, ly2, lz2>& g): CGTOPair(f, g) {}

	//! Return the first orbital in this pair
	const CGTOSpec<lx1, ly1, lz1>& f() const
	{
		return static_cast< const CGTOSpec<lx1, ly1, lz1>& >(AbstractBFPair::f());
	}
	//! Return the second orbital in this pair
	const CGTOSpec<lx1, ly1, lz1>& g() const
	{
		return static_cast< const CGTOSpec<lx1, ly1, lz1>& >(AbstractBFPair::g());
	}

	//! Compute the overlap between the two orbitals in this pair
	double overlap() const;
	//! Compute the kinetic energy integral between the two orbitals in this pair
	double kineticEnergy() const;
	/*!
	 * \brief Compute one-electron integrals
	 *
	 * Compute the overlap and kinetic energy integral between the two
	 * functions in this pair.
	 * \param S Place to store the overlap
	 * \param T Place to store the kinetic energy
	 */
	void oneElectron(double& S, double& T) const;

	/*!
	 * \brief Create a new CGTOSpecPair
	 *
	 * Create a new pairs of CGTO specializations. This pseudo-constructor
	 * provides a common interface that allows the Dispatcher to add pair
	 * creation functions for arbitrary pairs of basis function types.
	 * \param f The first orbital in the pair. Should be a CGTOSpec<lx1, ly1, lz1>.
	 * \param g The second orbital in the pair. Should be a CGTOSpec<lx2, ly2, lz2>.
	 */
	static AbstractBFPair *create(const AbstractBF& f, const AbstractBF& g)
	{
		try
		{
			const CGTOSpec<lx1, ly1, lz1>& ff = dynamic_cast< const CGTOSpec<lx1, ly1, lz1>& >(f);
			const CGTOSpec<lx2, ly2, lz2>& gg = dynamic_cast< const CGTOSpec<lx2, ly2, lz2>& >(g);
			return new CGTOSpecPair<lx1, ly1, lz1, lx2, ly2, lz2>(ff, gg);
		}
		catch (const std::bad_cast&)
		{
			throw Li::Exception("Invalid basis function type");
		}
	}
};

template <int lx1, int ly1, int lz1, int lx2, int ly2, int lz2>
double CGTOSpecPair<lx1, ly1, lz1, lx2, ly2, lz2>::overlap() const
{
	return Constants::pi_sqrt_pi * f().weights().transpose()
		* gto_overlap_primitive_specialized<lx1, ly1, lz1, lx2, ly2, lz2>(
			alpha(), beta(), asum(), ared(), exp_ared(),
			r()).matrix()
		* g().weights();
}

template <int lx1, int ly1, int lz1, int lx2, int ly2, int lz2>
double CGTOSpecPair<lx1, ly1, lz1, lx2, ly2, lz2>::kineticEnergy() const
{
	return Constants::pi_sqrt_pi * f().weights().transpose()
		* gto_kinetic_primitive_generic(f().ls(), g().ls(),
			alpha(), beta(), asum(), ared(), exp_ared(),
			r()).matrix()
		* g().weights();
}

template <int lx1, int ly1, int lz1, int lx2, int ly2, int lz2>
void CGTOSpecPair<lx1, ly1, lz1, lx2, ly2, lz2>::oneElectron(double &S, double &T) const
{
	Eigen::ArrayXXd Sp, Tp;
	gto_one_elec_primitive_generic(f().ls(), g().ls(),
		alpha(), beta(), asum(), ared(), exp_ared(), r(), Sp, Tp);
	S = Constants::pi_sqrt_pi * f().weights().transpose()
		* Sp.matrix() * g().weights();
	T = Constants::pi_sqrt_pi * f().weights().transpose()
		* Tp.matrix() * g().weights();
}

#endif // CGTOSPECPAIR_HH
