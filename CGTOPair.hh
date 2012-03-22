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
#include "MultiArray.hh"
#include "constants.hh"

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
	//! Unique class ID, used in looking up integral calculation functions
	static const size_t cid;

	//! Constructor
	CGTOPair(const CGTO& f, const CGTO& g):
		AbstractBFPair(cid, f, g),
		_shell_pair(CGTOShellList::singleton().pair(f.shell(), g.shell()))
		{}

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

	//! Return the pair of orbitals shells
	const CGTOShellPair& shellPair() const { return _shell_pair; }

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
	 * \brief Compute nuclear attraction integrals
	 * 
	 * Compute the nuclear attraction integrals, due to the nuclei with
	 * positions \a nuc_pos and charges \a nuc_charge, between the functions
	 * in this pair.
	 * \param nuc_pos    The positions of the nuclei
	 * \param nuc_charge The nuclear charges
	 */
	double nuclearAttraction(const Eigen::MatrixXd& nuc_pos,
		const Eigen::VectorXd& nuc_charge) const;

	/*!
	 * \brief Return the total number of primitive pairs in this
	 *    pair of contractions
	 */
	int size() const
	{
		return f().size() * g().size();
	}
	//! Return the angular momentum in the \a i direction for the first orbital
	int lA(int i) const
	{
		return f().l(i);
	}
	//! Return the angular momentum in the \a i direction for the second orbital
	int lB(int i) const
	{
		return g().l(i);
	}
	//! Return the sum of angular momenta in the \a i direction for the two orbitals
	int lsum(int i) const
	{
		return lA(i) + lB(i);
	}
	/*!
	 * \brief Return the widths of the primitives in the first contraction,
	 *    for all primitives in the second contraction.
	 */
	const Eigen::ArrayXXd& widthsA() const
	{
		//return f().widths().replicate(1, g().size());
		return _shell_pair.widthsA();
	}
	/*!
	 * \brief Return the widths of the primitives in the second
	 *    contraction, for all primitives in the first contraction.
	 */
	const Eigen::ArrayXXd& widthsB() const
	{
		//return g().widths().transpose().replicate(f().size(), 1);
		return _shell_pair.widthsB();
	}
	//! The sums of primitive widths, equivalent to widthsA() + widthsB()
	const Eigen::ArrayXXd& widthsSum() const
	{
		//return widthsA() + widthsB();
		return _shell_pair.widthsSum();
	}
	//! The "reduced" primitive widths \f$\xi = \alpha\beta / (\alpha+\beta)\f$
	const Eigen::ArrayXXd& widthsReduced() const
	{
		return _shell_pair.widthsReduced();
	}
	/*!
	 * \brief \f\frac{$\exp(-\xi r^2)}{\alpha+\beta}\f$ with \f$r\f$ the
	 *    distance between the centers of the two orbitals.
	 */
	const Eigen::ArrayXXd& gaussReduced() const
	{
		return _shell_pair.gaussReduced();
	}
	//! Return half the inverse widths sums, \f$\frac{1}{2(\alpha+\beta)}\f$
	const Eigen::ArrayXXd& hInvWidths() const
	{
		return _shell_pair.hInvWidths();
	}
	//! Return the first orbital center
	const Eigen::Vector3d& centerA() const
	{
		return f().center();
	}
	//! Return the \a i coordinate of the first orbital center
	double centerA(int i) const
	{
		return f().center(i);
	}
	//! Return the second orbital center
	const Eigen::Vector3d& centerB() const
	{
		return g().center();
	}
	//! Return the \a i coordinate of the second orbital center
	double centerB(int i) const
	{
		return g().center(i);
	}
	//! Return the vector from the first to the second orbital center
	Eigen::Vector3d r() const
	{
		return g().center() - f().center();
	}
	//! Return the distance in the \a i direction between the two centers
	double r(int i) const
	{
		return g().center(i) - f().center(i);
	}
	/*!
	 * \brief Return the weighted average \a i coordinate, for each
	 *    combination of primitive widths.
	 */
	const Eigen::ArrayXXd& P(int i) const
	{
		return _shell_pair.P(i);
	}
	Eigen::ArrayXXd K() const
	{
		return Constants::sqrt_2_pi_5_4 * gaussReduced();
	}

	/*!
	 * \brief Create a new CGTOPair
	 *
	 * Create a new pairs of CGTOs. This pseudo-constructor provides a
	 * common interface that allows the Dispatcher to add pair creation
	 * functions for arbitrary pairs of basis function types.
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

private:
	/*!
	 * Multiply the contributions to an integral \a C of each primitive pair
	 * with the weights of the primitives, and sum the results to compute
	 * the integral over the contraction.
	 * \param C the primitve integrals
	 */
	template <typename Derived>
	double mulWeights(const Eigen::ArrayBase<Derived>& C) const
	{
		return _shell_pair.mulWeights(C);
	}

	//! The shells for this pair of orbitals
	CGTOShellPair _shell_pair;

	//! Integrate the overlap matrix for this pair over dimension \a i.
	void overlapPrim1D(int i, Eigen::ArrayXXd& Sp) const;
	//! Integrate the overlap and kinetic energy matrices for this pair over dimension \a i.
	void oneElecPrim1D(int i, Eigen::ArrayXXd& Sp, Eigen::ArrayXXd& Tp) const;
	//! Integrate the coefficients for nuclear attraction integrals over dimension \a i.
	void nucAttrPrim1D(int i, const Eigen::ArrayXXd& theta,
		const Eigen::ArrayXXd& Pi, const Eigen::ArrayXXd& dPC,
		MultiArray& res) const;
};

#endif // CGTOPAIR_HH
