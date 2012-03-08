#ifndef CGTOSPECQUADCREATE_HH
#define CGTOSPECQUADCREATE_HH

#include "BFQuadPool.hh"

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
AbstractBFQuad* createCGTOSpecQuad(const AbstractBFPair& p,
	const AbstractBFPair& q, BFQuadPool& pool)
{
	try
	{
#ifdef DEBUG
		const CGTOPair& pp = dynamic_cast< const CGTOPair& >(p);
		const CGTOPair& qq = dynamic_cast< const CGTOPair& >(q);
#else
		const CGTOPair& pp = static_cast< const CGTOPair& >(p);
		const CGTOPair& qq = static_cast< const CGTOPair& >(q);
#endif
		int idx = CGTOShellList::pairIndex(pp.ishellPair(), qq.ishellPair());
		switch (CGTOShellList::singleton().quad(idx).positionSymmetry())
		{
			case CGTOShellQuad::POS_SYM_AAAA:
				return new(pool)
					CGTOSpecQuad<CGTOShellQuad::POS_SYM_AAAA, l1, l2>(pp, qq);
			case CGTOShellQuad::POS_SYM_AACC:
				return new(pool)
					CGTOSpecQuad<CGTOShellQuad::POS_SYM_AACC, l1, l2>(pp, qq);
			case CGTOShellQuad::POS_SYM_AACD:
				return new(pool)
					CGTOSpecQuad<CGTOShellQuad::POS_SYM_AACD, l1, l2>(pp, qq);
			case CGTOShellQuad::POS_SYM_ABCC:
				return new(pool)
					CGTOSpecQuad<CGTOShellQuad::POS_SYM_ABCC, l1, l2>(pp, qq);
			default:
				//return new(pool)
				//	CGTOSpecQuad<CGTOShellQuad::POS_SYM_ABCD, l1, l2>(pp, qq);
				return new(pool) CGTOQuad(pp, qq);
		}
	}
	catch (const std::bad_cast&)
	{
		throw Li::Exception("Invalid basis function pair type");
	}
}

#endif // CGTOSPECQUADCREATE_HH