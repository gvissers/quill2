#ifndef CGTOQUAD_HH
#define CGTOQUAD_HH

#include "AbstractBFQuad.hh"
#include "CGTOPair.hh"

#include <iostream>
class CGTOQuad: public AbstractBFQuad
{
public:
	/*!
	 * \brief Constructor
	 * 
	 * Create a new CGTOQuad from pairs \a p and \a q.
	 * \param p The first orbital pair in the quartet.
	 * \param q The second orbital pair in the quartet.
	 */
	CGTOQuad(const CGTOPair& p, const CGTOPair& q): AbstractBFQuad(p, q) {}
	
	//! Return the first orbital pair in the quartet
	const CGTOPair& p() const
	{
		return static_cast< const CGTOPair& >(AbstractBFQuad::p());
	}
	//! Return the second orbital in the pair
	const CGTOPair& q() const
	{
		return static_cast< const CGTOPair& >(AbstractBFQuad::q());
	}

	//! Reduced widths for the first pair, for all primitives in the second pair
	Eigen::ArrayXXd alpha() const
	{
		return Eigen::ArrayXd::Map(p().ared().data(), p().size())
			.replicate(1, q().size());
	}
	//! Reduced widths for the second pair, for all primitives in the first pair
	Eigen::ArrayXXd beta() const
	{
		return Eigen::ArrayXd::Map(q().ared().data(), q().size())
			.transpose().replicate(p().size(), 1);
	}
	//! Sums of reduced widths, for all primitives in the quartet
	Eigen::ArrayXXd asum() const
	{
		return alpha() + beta();
	}
	//! Overall reduced widths \f$\xi = \alpha\beta / (\alpha+\beta)\f$
	Eigen::ArrayXXd ared() const
	{
		return alpha()*beta() / asum();
	}
	/*!
	 * \brief Return the weighted average \a i coordinate, for each
	 *    combination of primitive weights.
	 */
	Eigen::ArrayXXd P(int i) const
	{
		return (alpha()*p().P(i) + beta()*q().P(i)) / asum();
	}

	/*!
	 * \brief Create a new CGTOQuad
	 *
	 * Create a new quartets of CGTOs. This pseudo-constructor provides a
	 * common interface that allows the Dispatcher to add quartet creation
	 * functions for arbitrary quartets of basis function types.
	 * \param p The first orbital pair in the quartet. Should be a CGTOPair.
	 * \param q The second orbital pair in the quartet. Should be a CGTOPair.
	 */
	static AbstractBFQuad *create(const AbstractBFPair& p,
		const AbstractBFPair& q)
	{
		try
		{
			const CGTOPair& pp = dynamic_cast< const CGTOPair& >(p);
			const CGTOPair& qq = dynamic_cast< const CGTOPair& >(q);
			return new CGTOQuad(pp, qq);
		}
		catch (const std::bad_cast&)
		{
			throw Li::Exception("Invalid basis function pair type");
		}
	}
};

#endif // CGTOQUAD_HH