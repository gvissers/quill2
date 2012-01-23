#ifndef ABSTRACTBFPAIR_HH
#define ABSTRACTBFPAIR_HH

/*!
 * \file AbstractBFPair.hh
 * \brief Definition of the AbstractBDPair class
 */

#include "AbstractBF.hh"

/*!
 * \brief Class for a pair of basis functions
 *
 * Class AbstractBFPair defines the interface for a tuple of two basis
 * functions. It specifies the functions that can be computed on these
 * pairs, inclusing overlap and kinetic energy integrals. Concrete child
 * classes deriving from AbstractBFPair are responsible for providing
 * implementations for these integrals for the basis function types they
 * are holding.
 */
class AbstractBFPair
{
public:
	//! Unique class ID, used in looking up integral calculation functions
	size_t cid;

	//! Constructor
	AbstractBFPair(size_t cid, const AbstractBF& f, const AbstractBF& g):
		cid(cid), _f(f), _g(g) {}

	//! Return the first function in the pair
	const AbstractBF& f() const { return _f; }
	//! Return the second function in the pair
	const AbstractBF& g() const { return _g; }

	/*!
	 * \brief Compute the overlap between the two functions
	 *
	 * Compute the overlap between the two functions in this pair. The
	 * default implementation in this class calls oneElectron(), and returns
	 * only the overlap. This function can be overridden by potentially more
	 * efficient methods in classes deriving from AbstractBFPair.
	 */
	virtual double overlap() const
	{
		double S, T;
		oneElectron(S, T);
		return S;
	}
	/*!
	 * \brief Compute the kinetic energy between the two functions
	 *
	 * Compute the kinetic energy integral between the two functions in this
	 * pair. The default implementation in this class calls oneElectron(),
	 * and returns only the kinetic energy. This function can be overridden
	 * by potentially more efficient methods in classes deriving from
	 * AbstractBFPair.
	 */
	virtual double kineticEnergy() const
	{
		double S, T;
		oneElectron(S, T);
		return T;
	}
	/*!
	 * \brief Compute the one-electron integrals
	 *
	 * Compute the overlap and kinetic energy integral between the two
	 * functions in this pair. Child classes should implement this function.
	 * \param S Place to store the overlap
	 * \param T Place to store the kinetic energy
	 */
	virtual void oneElectron(double& S, double& T) const = 0;

	/*!
	 * Compute the nuclear attraction integrals, due to the nuclei with
	 * positions \a nuc_pos and charges \a nuc_charge, between the functions
	 * in this pair. Child classes should implement this function.
	 * \param nuc_pos    The positions of the nuclei
	 * \param nuc_charge The nuclear charges
	 */
	virtual double nuclearAttraction(const Eigen::MatrixXd& nuc_pos,
		const Eigen::VectorXd& nuc_charge) const = 0;
		
	//virtual double electronRepulsion(const AbstractBFPair& pair) const = 0;

private:
	//! The first function in the pair
	const AbstractBF& _f;
	//! The second function in the pair
	const AbstractBF& _g;
};

#endif // ABSTRACTBFPAIR_HH
