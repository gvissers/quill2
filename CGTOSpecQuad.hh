#ifndef CGTOSPECQUAD_HH
#define CGTOSPECQUAD_HH

/*!
 * \file CGTOSpecQuad.hh
 * \brief Definition of the CGTOSpecQuad class
 */

#include "CGTOQuad.hh"

/*!
 * \brief Specialized orbital quartets
 *
 * Class CGTOSpecQuad holds specialized versions of quartets od contracted
 * Gaussian type orbitals that are used to compute the electron repulsion
 * integrals. The electronRepulsion() method that computes these integrals can
 * be overridden with faster specializations than the generic method.
 * \tparam pos_sym The symmetry in positions of the four orbital centers
 * \tparam l1      The total angular momentum of the first pair in the quartet
 * \tparam l2      The total angular momentum of the second pair in the quartet
 */
template <enum CGTOShellQuad::PositionSymmetry pos_sym, int l1, int l2>
struct CGTOSpecQuad: public CGTOQuad
{
	//! Constructor
	CGTOSpecQuad(const CGTOPair& p, const CGTOPair& q): CGTOQuad(p, q) {}

	/*!
	 * Compute the electron repulsion integral
	 * \f$\langle ac|\frac{1}{r_{12}}|bd\rangle\f$, where orbitals \f$a\f$
	 * and \f$b\f$ are stored in the first orbital pair of this quartet, and
	 * \f$c\f$ and \f$d\f$ are stored in the second pair.
	 */
	double electronRepulsion() const
	{
		return CGTOQuad::electronRepulsion();
	}
};

template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_AAAA, 0, 0>::electronRepulsion() const;
template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_AAAA, 1, 0>::electronRepulsion() const;
template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_AAAA, 0, 1>::electronRepulsion() const;
template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_AAAA, 2, 0>::electronRepulsion() const;
template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_AAAA, 1, 1>::electronRepulsion() const;
template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_AAAA, 0, 2>::electronRepulsion() const;
template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_AAAA, 2, 1>::electronRepulsion() const;
template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_AAAA, 1, 2>::electronRepulsion() const;

template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_AACC, 0, 0>::electronRepulsion() const;
template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_AACC, 1, 0>::electronRepulsion() const;
template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_AACC, 0, 1>::electronRepulsion() const;
template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_AACC, 2, 0>::electronRepulsion() const;
template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_AACC, 1, 1>::electronRepulsion() const;
template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_AACC, 0, 2>::electronRepulsion() const;

template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_ABCC, 0, 0>::electronRepulsion() const;
template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_ABCC, 1, 0>::electronRepulsion() const;
template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_ABCC, 0, 1>::electronRepulsion() const;
template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_ABCC, 2, 0>::electronRepulsion() const;
template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_ABCC, 1, 1>::electronRepulsion() const;
template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_ABCC, 0, 2>::electronRepulsion() const;

template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_AACD, 0, 0>::electronRepulsion() const;
template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_AACD, 1, 0>::electronRepulsion() const;
template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_AACD, 0, 1>::electronRepulsion() const;
template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_AACD, 2, 0>::electronRepulsion() const;
template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_AACD, 1, 1>::electronRepulsion() const;
template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_AACD, 0, 2>::electronRepulsion() const;

template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_ABCD, 0, 0>::electronRepulsion() const;
template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_ABCD, 1, 0>::electronRepulsion() const;
template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_ABCD, 0, 1>::electronRepulsion() const;
template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_ABCD, 2, 0>::electronRepulsion() const;
template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_ABCD, 1, 1>::electronRepulsion() const;
template <>
double CGTOSpecQuad<CGTOShellQuad::POS_SYM_ABCD, 0, 2>::electronRepulsion() const;

/*!
 * \brief Create a new CGTOSpecQuad
 *
 * Create a new quartet of CGTOs. This pseudo-constructor provides a common
 * interface that allows the Dispatcher to add quartet creation functions for
 * arbitrary quartets of basis function types.
 * \tparam l1 The total angular momentum of the first orbital pair.
 * \tparam l2 The total angular momentum of the second orbital pair.
 * \param p   The first orbital pair in the quartet. Should be a CGTOPair.
 * \param q   The second orbital pair in the quartet. Should be a CGTOPair.
 */
template <int l1, int l2>
AbstractBFQuad* createSpecQuad(const AbstractBFPair& p, const AbstractBFPair& q)
{
	try
	{
		const CGTOPair& pp = dynamic_cast< const CGTOPair& >(p);
		const CGTOPair& qq = dynamic_cast< const CGTOPair& >(q);
		int idx = CGTOShellList::pairIndex(pp.ishellPair(), qq.ishellPair());
		switch (CGTOShellList::singleton().quad(idx).positionSymmetry())
		{
			case CGTOShellQuad::POS_SYM_AAAA:
				return new CGTOSpecQuad<CGTOShellQuad::POS_SYM_AAAA, l1, l2>(pp, qq);
			case CGTOShellQuad::POS_SYM_AACC:
				return new CGTOSpecQuad<CGTOShellQuad::POS_SYM_AACC, l1, l2>(pp, qq);
			case CGTOShellQuad::POS_SYM_AACD:
				return new CGTOSpecQuad<CGTOShellQuad::POS_SYM_AACD, l1, l2>(pp, qq);
			case CGTOShellQuad::POS_SYM_ABCC:
				return new CGTOSpecQuad<CGTOShellQuad::POS_SYM_ABCC, l1, l2>(pp, qq);
			default:
				return new CGTOSpecQuad<CGTOShellQuad::POS_SYM_ABCD, l1, l2>(pp, qq);
		}
	}
	catch (const std::bad_cast&)
	{
		throw Li::Exception("Invalid basis function pair type");
	}
}

#endif // CGTOSPECQUAD_HH