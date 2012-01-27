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
	
	double electronRepulsion() const;

	/*!
	 * \brief Sums of primitive widths for the first pair, for all
	 *     primitives in the second pair
	 */
	Eigen::ArrayXXd widthsAB() const
	{
		return Eigen::ArrayXd::Map(p().widthsSum().data(), p().size())
			.replicate(1, q().size());
	}
	/*!
	 * \brief Sums of primitive widths for the second pair, for all
	 *    primitives in the first pair
	 */
	Eigen::ArrayXXd widthsCD() const
	{
		return Eigen::ArrayXd::Map(q().widthsSum().data(), q().size())
			.transpose().replicate(p().size(), 1);
	}
	//! Sums of primitve widths, for all combination of primitives GTOs
	Eigen::ArrayXXd widthsSum() const
	{
		return widthsAB() + widthsCD();
	}
	//! Reduced widths \f$\rho = \frac{\zeta\eta}{\zeta+\eta}\f$
	Eigen::ArrayXXd widthsReduced() const
	{
		return widthsAB() * widthsCD() / widthsSum();
	}
	
	Eigen::VectorXd weightsAB() const
	{
		return Eigen::VectorXd::Map(p().weights().data(), p().size());
	}
	Eigen::VectorXd weightsCD() const
	{
		return Eigen::VectorXd::Map(q().weights().data(), q().size());
	}

	Eigen::ArrayXXd P(int i) const
	{
		return Eigen::ArrayXd::Map(p().P(i).data(), p().size())
			.replicate(1, q().size());
	}
	Eigen::ArrayXXd Q(int i) const
	{
		return Eigen::ArrayXd::Map(q().P(i).data(), q().size())
			.transpose().replicate(p().size(), 1);
	}
	
	/*!
	 * \brief Return the weighted average \a i coordinate, for each
	 *    combination of primitive weights.
	 */
	Eigen::ArrayXXd W(int i) const
	{
		return (widthsAB()*P(i) + widthsCD()*Q(i)) / widthsSum();
	}
	Eigen::ArrayXXd KK() const
	{
		Eigen::ArrayXd K1 = Eigen::ArrayXd::Map(p().K().data(), p().size());
		Eigen::ArrayXd K2 = Eigen::ArrayXd::Map(q().K().data(), q().size());
		return K1.replicate(1, q().size())
			* K2.transpose().replicate(p().size(), 1);
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