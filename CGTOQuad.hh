#ifndef CGTOQUAD_HH
#define CGTOQUAD_HH

#include "AbstractBFQuad.hh"
#include "CGTOPair.hh"

// Forward declarations
class BFQuadPool;

class CGTOQuad: public AbstractBFQuad
{
public:
	/*!
	 * \brief Constructor
	 *
	 * Create a new CGTOQuad from pairs \a p and \a q.
	 * \param p       The first orbital pair in the quartet.
	 * \param q       The second orbital pair in the quartet.
	 */
	CGTOQuad(const CGTOPair& pp, const CGTOPair& qq);
	//! Destructor
	virtual ~CGTOQuad() {}

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
	
	/*!
	 * Compute the electron repulsion integral
	 * \f$\langle ac|\frac{1}{r_{12}}|bd\rangle\f$, where orbitals \f$a\f$
	 * and \f$b\f$ are stored in the first orbital pair of this quartet, and 
	 * \f$c\f$ and \f$d\f$ are stored in the second pair.
	 */
	double electronRepulsion() const
	{
		return _shell_quad.getEri(lA(0), lA(1), lA(2),
			lB(0), lB(1), lB(2),
			lC(0), lC(1), lC(2),
			lD(0), lD(1), lD(2)) * _norm;
	}

	//! Return the angular momentum in the \a i direction for the first orbital
	int lA(int i) const
	{
		return p().lA(i);
	}
	//! Return the angular momentum in the \a i direction for the second orbital
	int lB(int i) const
	{
		return p().lB(i);
	}
	//! Return the angular momentum in the \a i direction for the third orbital
	int lC(int i) const
	{
		return q().lA(i);
	}
	//! Return the angular momentum in the \a i direction for the fourth orbital
	int lD(int i) const
	{
		return q().lB(i);
	}

	/*!
	 * \brief Create a new CGTOQuad
	 *
	 * Create a new quartet of CGTOs. This pseudo-constructor provides a
	 * common interface that allows the Dispatcher to add quartet creation
	 * functions for arbitrary quartets of basis function types.
	 * \param p The first orbital pair in the quartet. Should be a CGTOPair.
	 * \param q The second orbital pair in the quartet. Should be a CGTOPair.
	 */
	static AbstractBFQuad *create(const AbstractBFPair& p,
		const AbstractBFPair& q, BFQuadPool& pool);

protected:
	/*!
	 * Multiply the contributions to an integral \a C of each primitive pair
	 * with the weights of the primitives, and sum the results to compute
	 * the integral over the contraction.
	 * \param C the primitve integrals
	 */
	template <typename Derived>
	double mulWeights(const Eigen::ArrayBase<Derived>& C) const
	{
		return _shell_quad.mulWeights(C) * _norm;
	}

private:
	//! Normalization factor for two-electron integrals
	double _norm;
	//! The combined index of this quartet's shells in the CGTOShellList
	int _ishell_quad;
	//! The shells for this quartet of orbitals
	const CGTOShellQuad& _shell_quad;
};

#endif // CGTOQUAD_HH
