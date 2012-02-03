#ifndef CGTOSPECPAIR_HH
#define CGTOSPECPAIR_HH

/*!
 * \file CGTOSpecPair.hh
 * \brief Definition of the CGTOSpecPair class
 */

#include <Eigen/Dense>
#include "CGTOPair.hh"
#include "CGTOSpec.hh"
#include "gaussint/gto_nuc_attr.hh"
#include "constants.hh"

/*!
 * \brief Class for CGTO specializations
 *
 * Class CGTOSpecPair holds a pair of specializations of contracted Gaussian
 * type orbitals, for specific angular momentum quantum numbers.
 * \tparam lxA Angular momentum in the \f$x\f$ direction for the first orbital
 * \tparam lyA Angular momentum in the \f$y\f$ direction for the first orbital
 * \tparam lzA Angular momentum in the \f$z\f$ direction for the first orbital
 * \tparam lxB Angular momentum in the \f$x\f$ direction for the second orbital
 * \tparam lyB Angular momentum in the \f$y\f$ direction for the seoncd orbital
 * \tparam lzB Angular momentum in the \f$z\f$ direction for the second orbital
 */
template <int lxA, int lyA, int lzA, int lxB, int lyB, int lzB>
struct CGTOSpecPair: public CGTOPair
{
	//! Constructor
	CGTOSpecPair(const CGTOSpec<lxA, lyA, lzA>& f,
		const CGTOSpec<lxB, lyB, lzB>& g): CGTOPair(f, g) {}

	//! Return the first orbital in this pair
	const CGTOSpec<lxA, lyA, lzA>& f() const
	{
		return static_cast< const CGTOSpec<lxA, lyA, lzA>& >(AbstractBFPair::f());
	}
	//! Return the second orbital in this pair
	const CGTOSpec<lxB, lyB, lzB>& g() const
	{
		return static_cast< const CGTOSpec<lxB, lyB, lzB>& >(AbstractBFPair::g());
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
	 * Compute the nuclear attraction integrals, due to the nuclei
	 * with positions \a nuc_pos and charges \a nuc_charge, between
	 * the functions in this pair.
	 * \param nuc_pos    The positions of the nuclei
	 * \param nuc_charge The nuclear charges
	 */
	double nuclearAttraction(const Eigen::MatrixXd& nuc_pos,
		const Eigen::VectorXd& nuc_charge) const;

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
			const CGTOSpec<lxA, lyA, lzA>& ff = dynamic_cast< const CGTOSpec<lxA, lyA, lzA>& >(f);
			const CGTOSpec<lxB, lyB, lzB>& gg = dynamic_cast< const CGTOSpec<lxB, lyB, lzB>& >(g);
			return new CGTOSpecPair<lxA, lyA, lzA, lxB, lyB, lzB>(ff, gg);
		}
		catch (const std::bad_cast&)
		{
			throw Li::Exception("Invalid basis function type");
		}
	}
};

template <int lxA, int lyA, int lzA, int lxB, int lyB, int lzB>
double CGTOSpecPair<lxA, lyA, lzA, lxB, lyB, lzB>::overlap() const
{
	return CGTOPair::overlap();
}
template <>
double CGTOSpecPair<0, 0, 0, 0, 0, 0>::overlap() const;
template <>
double CGTOSpecPair<1, 0, 0, 0, 0, 0>::overlap() const;
template <>
double CGTOSpecPair<0, 1, 0, 0, 0, 0>::overlap() const;
template <>
double CGTOSpecPair<0, 0, 1, 0, 0, 0>::overlap() const;
template <>
double CGTOSpecPair<0, 0, 0, 1, 0, 0>::overlap() const;
template <>
double CGTOSpecPair<0, 0, 0, 0, 1, 0>::overlap() const;
template <>
double CGTOSpecPair<0, 0, 0, 0, 0, 1>::overlap() const;
template <>
double CGTOSpecPair<2, 0, 0, 0, 0, 0>::overlap() const;
template <>
double CGTOSpecPair<1, 1, 0, 0, 0, 0>::overlap() const;
template <>
double CGTOSpecPair<1, 0, 1, 0, 0, 0>::overlap() const;
template <>
double CGTOSpecPair<1, 0, 0, 1, 0, 0>::overlap() const;
template <>
double CGTOSpecPair<1, 0, 0, 0, 1, 0>::overlap() const;
template <>
double CGTOSpecPair<1, 0, 0, 0, 0, 1>::overlap() const;
template <>
double CGTOSpecPair<0, 2, 0, 0, 0, 0>::overlap() const;
template <>
double CGTOSpecPair<0, 1, 1, 0, 0, 0>::overlap() const;
template <>
double CGTOSpecPair<0, 1, 0, 1, 0, 0>::overlap() const;
template <>
double CGTOSpecPair<0, 1, 0, 0, 1, 0>::overlap() const;
template <>
double CGTOSpecPair<0, 1, 0, 0, 0, 1>::overlap() const;
template <>
double CGTOSpecPair<0, 0, 2, 0, 0, 0>::overlap() const;
template <>
double CGTOSpecPair<0, 0, 1, 1, 0, 0>::overlap() const;
template <>
double CGTOSpecPair<0, 0, 1, 0, 1, 0>::overlap() const;
template <>
double CGTOSpecPair<0, 0, 1, 0, 0, 1>::overlap() const;
template <>
double CGTOSpecPair<0, 0, 0, 2, 0, 0>::overlap() const;
template <>
double CGTOSpecPair<0, 0, 0, 1, 1, 0>::overlap() const;
template <>
double CGTOSpecPair<0, 0, 0, 1, 0, 1>::overlap() const;
template <>
double CGTOSpecPair<0, 0, 0, 0, 2, 0>::overlap() const;
template <>
double CGTOSpecPair<0, 0, 0, 0, 1, 1>::overlap() const;
template <>
double CGTOSpecPair<0, 0, 0, 0, 0, 2>::overlap() const;

template <int lxA, int lyA, int lzA, int lxB, int lyB, int lzB>
double CGTOSpecPair<lxA, lyA, lzA, lxB, lyB, lzB>::kineticEnergy() const
{
	return CGTOPair::kineticEnergy();
}
template<>
double CGTOSpecPair<0, 0, 0, 0, 0, 0>::kineticEnergy() const;
template<>
double CGTOSpecPair<1, 0, 0, 0, 0, 0>::kineticEnergy() const;
template<>
double CGTOSpecPair<0, 1, 0, 0, 0, 0>::kineticEnergy() const;
template<>
double CGTOSpecPair<0, 0, 1, 0, 0, 0>::kineticEnergy() const;
template<>
double CGTOSpecPair<0, 0, 0, 1, 0, 0>::kineticEnergy() const;
template<>
double CGTOSpecPair<0, 0, 0, 0, 1, 0>::kineticEnergy() const;
template<>
double CGTOSpecPair<0, 0, 0, 0, 0, 1>::kineticEnergy() const;

template <int lxA, int lyA, int lzA, int lxB, int lyB, int lzB>
void CGTOSpecPair<lxA, lyA, lzA, lxB, lyB, lzB>::oneElectron(double &S, double &T) const
{
	CGTOPair::oneElectron(S, T);
}

template <int lxA, int lyA, int lzA, int lxB, int lyB, int lzB>
double CGTOSpecPair<lxA, lyA, lzA, lxB, lyB, lzB>::nuclearAttraction(
	const Eigen::MatrixXd& nuc_pos, const Eigen::VectorXd& nuc_charge) const
{
	return 2 * M_PI * f().weights().transpose()
		* gto_nuc_attr_primitive_specialized<lxA, lyA, lzA, lxB, lyB, lzB>(
			widthsA(), widthsB(), f().center(), g().center(),
			widthsSum(), exp_ared(), nuc_pos, nuc_charge).matrix()
		* g().weights();
}

#endif // CGTOSPECPAIR_HH
