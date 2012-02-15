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
		_widths_A(f.widths().replicate(1, g.size())),
		_widths_B(g.widths().transpose().replicate(f.size(), 1)),
		_widths_sum(_widths_A + _widths_B),
		_widths_red(_widths_A * _widths_B / _widths_sum),
		_gauss_red((-(g.center()-f.center()).squaredNorm() * _widths_red).exp()
			/ _widths_sum),
		_weights(f.weights() * g.weights().transpose())
	{
		for (int i = 0; i < 3; ++i)
		{
			_P[i] = (_widths_A*f.center(i) + _widths_B*g.center(i))
				/ _widths_sum;
		}
	}

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
	 * Compute the overlap and kinetic energy integral between the two
	 * functions in this pair.
	 * \param S Place to store the overlap
	 * \param T Place to store the kinetic energy
	 */
	virtual void oneElectron(double& S, double& T) const;
	/*!
	 * \brief Compute nuclear attraction integrals
	 * 
	 * Compute the nuclear attraction integrals, due to the nuclei with
	 * positions \a nuc_pos and charges \a nuc_charge, between the functions
	 * in this pair.
	 * \param nuc_pos    The positions of the nuclei
	 * \param nuc_charge The nuclear charges
	 */
	virtual double nuclearAttraction(const Eigen::MatrixXd& nuc_pos,
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
		return _widths_A;
	}
	/*!
	 * \brief Return the widths of the primitives in the second
	 *    contraction, for all primitives in the first contraction.
	 */
	const Eigen::ArrayXXd& widthsB() const
	{
		//return g().widths().transpose().replicate(f().size(), 1);
		return _widths_B;
	}
	//! The sums of primitive widths, equivalent to alpha() + beta()
	const Eigen::ArrayXXd& widthsSum() const
	{
		//return widthsA() + widthsB();
		return _widths_sum;
	}
	//! The "reduced" primitive widths \f$\xi = \alpha\beta / (\alpha+\beta)\f$
	const Eigen::ArrayXXd& widthsReduced() const
	{
		return _widths_red;
	}
	/*!
	 * \brief \f\frac{$\exp(-\xi r^2)}{\alpha+\beta}\f$ with \f$r\f$ the
	 *    distance between the centers of the two orbitals.
	 */
	const Eigen::ArrayXXd& gaussReduced() const
	{
		return _gauss_red;
	}
	int positionIdA() const
	{
		return f().positionId();
	}
	int positionIdB() const
	{
		return g().positionId();
	}
	bool samePositionId() const
	{
		return positionIdA() == positionIdB();
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
	 *    combination of primitive weights.
	 */
	const Eigen::ArrayXXd& P(int i) const
	{
		return _P[i];
	}
	Eigen::ArrayXXd K() const
	{
		return Constants::sqrt_2_pi_5_4 * gaussReduced();
	}
	/*!
	 * \brief Return the products of the weights for each combination of
	 *    primitives in the contraction
	 */
	const Eigen::ArrayXXd& weights() const
	{
		return _weights;
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

protected:
	/*!
	 * \brief Constructor
	 *
	 * Create a new pair of contracted gaussian type orbitals \a f and \a g.
	 * This constructor is called by specialized child classes, that pass
	 * their own class ID \a cid.
	 */
	CGTOPair(size_t cid, const CGTO& f, const CGTO& g):
		AbstractBFPair(cid, f, g),
		_widths_A(f.widths().replicate(1, g.size())),
		_widths_B(g.widths().transpose().replicate(f.size(), 1)),
		_widths_sum(_widths_A + _widths_B),
		_widths_red(_widths_A * _widths_B / _widths_sum),
		_gauss_red((-(g.center()-f.center()).squaredNorm() * _widths_red).exp()
			/ _widths_sum),
		_weights(f.weights() * g.weights().transpose())
	{
		for (int i = 0; i < 3; ++i)
		{
			_P[i] = (_widths_A*f.center(i) + _widths_B*g.center(i))
				/ _widths_sum;
		}
	}


private:
	//! Primitives widths in first orbital, for all primitives in second orbital.
	Eigen::ArrayXXd _widths_A;
	//! Primitives widths in second orbital, for all primitives in first orbital.
	Eigen::ArrayXXd _widths_B;
	//! Sum of primitive widths, for all combinations of primitives
	Eigen::ArrayXXd _widths_sum;
	//! Reduced primitive widths, for all combinations of primitives
	Eigen::ArrayXXd _widths_red;
	//! \f$\frac{\exp(-\xi r^2)}{\alpha+\beta}\f$
	Eigen::ArrayXXd _gauss_red;
	//! Weighted average coordinates
	Eigen::ArrayXXd _P[3];
	//! Products of the weights, for all combinations of primitives
	Eigen::ArrayXXd _weights;

	//! Integrate the overlap matrix for this pair over dimension \a i.
	void overlapPrim1D(int i, Eigen::ArrayXXd& Sp) const;
	//! Integrate the overlap and kinetic energy matrices for this pair over dimension \a i.
	void oneElecPrim1D(int i, Eigen::ArrayXXd& Sp, Eigen::ArrayXXd& Tp) const;
	//! Integrate the coefficients for nuclear attraction integrals over dimension \a i.
	void nucAttrPrim1D(int i, const Eigen::ArrayXXd& theta,
		const Eigen::ArrayXXd& Pi, const Eigen::ArrayXXd& dPC,
		std::vector<Eigen::ArrayXXd>& res) const;
};

#endif // CGTOPAIR_HH
