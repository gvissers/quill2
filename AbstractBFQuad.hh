#ifndef ABSTRACTBFQUAD_HH
#define ABSTRACTBFQUAD_HH

/*!
 * \file AbstractBFQuad.hh
 * \brief Definition of the AbstractBFQuad class
 */

#include "AbstractBFPair.hh"

/*!
 * \brief Class for a pair of basis functions
 *
 * Class AbstractBFQuad defines the interface for a tuple of four basis
 * functions, used in computing electron repulsion integrals.
 */
class AbstractBFQuad
{
public:
	//! Return the first pair in the quartet
	const AbstractBFPair& p() const { return _p; }
	//! Return the first pair in the quartet
	const AbstractBFPair& q() const { return _q; }
	
	/*!
	 * Compute the electron repulsion integral
	 * \f$\langle ab|\frac{1}{r_{12}}|cd\rangle\f$, where orbitals \f$a\f$
	 * and \f$b\f$ are stored in the first orbital pair of this quartet, and 
	 * \f$c\f$ and \f$d\f$ are stored in the second pair. Child classes
	 * should implement this function.
	 */
	virtual double electronRepulsion() const = 0;

protected:
	/*!
	 * \brief Constructor
	 * 
	 * Create a new AbstractBFQuad. For computing the electron repulsion
	 * integral \f$\langle ab|\frac{1}{r_{12}}|cd\rangle\f$, \a p should
	 * contain orbitals \f$a\f$ and \f$b\f$, and \a q should consist of
	 * \f$c\f$ and \f$d\f$.
	 * \param p The first orbital pair in the quartet
	 * \param q The second orbital pair in the quartet
	 */
	AbstractBFQuad(const AbstractBFPair& p,
		const AbstractBFPair& q): _p(p), _q(q) {}

private:
	//! The first orbital pair in the quartet
	const AbstractBFPair& _p;
	//! The second orbital pair in the quartet
	const AbstractBFPair& _q;
};

#endif // ABSTRACTBFQUAD_HH