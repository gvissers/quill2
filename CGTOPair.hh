#ifndef CGTOPAIR_HH
#define CGTOPAIR_HH

/*!
 * \file CGTOPair.hh
 * \brief Definition of the CGTOPair class
 */

#include <typeinfo>
#include <Eigen/Core>
#include "AbstractBFPair.hh"
#include "CGTO.hh"

/*!
 * \brief Class for pairs of contracted GTOs
 *
 * Class CGTOPair holds a pair of contracted Gaussian type orbitals, and
 * defines a set of operations on this pair, like computing the overlap
 * or kinetic energy integral.
 */
class CGTOPair: public AbstractBFPair
{
	public:
		//! Constructor
		CGTOPair(const CGTO& f, const CGTO& g): AbstractBFPair(f, g) {}

		//! Return the first orbital in the pair
		const CGTO& f() const
		{
			return static_cast< const CGTO& >(AbstractBFPair::f());
		}
		//! Return the second orbital in the pair
		const CGTO& g() const
		{
			return static_cast< const CGTO& >(AbstractBFPair::g());
		}

		//! Compute the overlap between the two orbitals in this pair
		virtual double overlap() const;
		//! Compute the kinetic energy integral between the two orbitals in this pair
		virtual double kineticEnergy() const;
		/*!
		 * \brief Compute one-electron integrals
		 *
		 * Compute the overlap and kinetic energy integral between the
		 * two functions in this pair.
		 * \param S Place to store the overlap
		 * \param T Place to store the kinetic energy
		 */
		virtual void oneElectron(double& S, double& T) const;

		/*!
		 * Compute the nuclear attraction integrals, due to the nuclei
		 * with positions \a nuc_pos and charges \a nuc_charge, between
		 * the functions in this pair.
		 * \param nuc_pos    The positions of the nuclei
		 * \param nuc_charge The nuclear charges
		 */
		virtual double nuclearAttraction(const Eigen::MatrixXd& nuc_pos,
			const Eigen::VectorXd& nuc_charge) const;

		/*!
		 * \brief Return the widths of the primitives in the first
		 *    contraction, for all primitives in the second contraction.
		 */
		Eigen::ArrayXXd alpha() const
		{
			return f().widths().replicate(1, g().size());
		}
		/*!
		 * \brief Return the widths of the primitives in the second
		 *    contraction, for all primitives in the first contraction.
		 */
		Eigen::ArrayXXd beta() const
		{
			return g().widths().transpose().replicate(f().size(), 1);
		}
		//! The sums of primitive widths, equivalent to alpha() + beta()
		Eigen::ArrayXXd asum() const
		{
			return alpha() + beta();
		}
		//! The "reduced" primitive widths \f$\xi = \alpha\beta / (\alpha+\beta)\f$
		Eigen::ArrayXXd ared() const
		{
			return alpha()*beta() / asum();
		}
		//! \f$\exp(-\xi r^2)\f$ with \f$r\f$ the distance between the centers of the two orbitals
		Eigen::ArrayXXd exp_ared() const
		{
			return (-r().squaredNorm()*ared()).exp();
		}
		//! Return the vector from the first to the second orbital center
		Eigen::Vector3d r() const
		{
			return g().center() - f().center();
		}

		/*!
		 * \brief Create a new CGTOPair
		 *
		 * Create a new pairs of CGTOs. This pseudo-constructor
		 * provides a common interface that allows the Dispatcher to
		 * add pair creation functions for arbitrary pairs of basis
		 * function types.
		 * \param f The first orbital in the pair. Should be a CGTO.
		 * \param g The second orbital in the pair. Should be a CGTO.
		 */
		static AbstractBFPair *create(const AbstractBF& f, const AbstractBF& g)
		{
			try
			{
				const CGTO& ff = dynamic_cast< const CGTO& >(f);
				const CGTO& gg = dynamic_cast< const CGTO& >(g);
				return new CGTOPair(ff, gg);
			}
			catch (const std::bad_cast&)
			{
				throw Li::Exception("Invalid basis function type");
			}
		}
};

#endif // CGTOPAIR_HH
