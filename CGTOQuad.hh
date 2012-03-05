#ifndef CGTOQUAD_HH
#define CGTOQUAD_HH

#include "AbstractBFQuad.hh"
#include "CGTOPair.hh"

// Forward declarations
class FmCoefs;
class BFQuadPool;

class CGTOQuad: public AbstractBFQuad
{
public:
	typedef CGTOShellQuad::ColArray ColArray;
	typedef CGTOShellQuad::RowArray RowArray;
	typedef CGTOShellQuad::ColArrayMap ColArrayMap;
	typedef CGTOShellQuad::RowArrayMap RowArrayMap;
	typedef CGTOShellQuad::WidthsReducedExpression WidthsReducedExpression;
	typedef CGTOShellQuad::DPQExpression DPQExpression;
	typedef CGTOShellQuad::DXPExpression DXPExpression;
	typedef CGTOShellQuad::DXQExpression DXQExpression;
	typedef Eigen::VectorXd::ConstAlignedMapType VectorMap;
	typedef CGTOShellQuad::DPWExpression DPWExpression;
	typedef CGTOShellQuad::DQWExpression DQWExpression;
	typedef CGTOShellQuad::Rho1Expression Rho1Expression;
	typedef CGTOShellQuad::Rho2Expression Rho2Expression;
	typedef ColArrayMap HInvWidthsABExpression;
	typedef RowArrayMap HInvWidthsCDExpression;
	typedef CGTOShellQuad::KKWExpression KKWExpression;

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
	virtual double electronRepulsion() const;

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
	//! Return the sum of angular momenta for the orbitals in the first pair
	int lAB(int i) const
	{
		return p().lsum(i);
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
	//! Return the sum of angular momenta for the orbitals in the first pair
	int lCD(int i) const
	{
		return q().lsum(i);
	}
	//! Return the sum of all for angular momenta in direction \a i.
	int lsum(int i) const
	{
		return lAB(i) + lCD(i);
	}

	//! Sums of primitive widths for the first pair
	ColArrayMap widthsAB() const
	{
		return _shell_quad.widthsAB();
	}
	//! Sums of primitive widths for the second pair
	RowArrayMap widthsCD() const
	{
		return _shell_quad.widthsCD();
	}
	//! Sums of primitive widths, for all combinations of primitives GTOs
	const Eigen::ArrayXXd& invWidthsSum() const
	{
		return _shell_quad.invWidthsSum();
	}
	//! Reduced widths \f$\rho = \frac{\zeta\eta}{\zeta+\eta}\f$
	WidthsReducedExpression widthsReduced() const
	{
		return _shell_quad.widthsReduced();
	}

	//! Return the \a i coordinate of the first orbital center
	double centerA(int i) const
	{
		return p().centerA(i);
	}
	//! Return the \a i coordinate of the second orbital center
	double centerB(int i) const
	{
		return p().centerB(i);
	}
	//! Return the \a i coordinate of the third orbital center
	double centerC(int i) const
	{
		return q().centerA(i);
	}
	//! Return the \a i coordinate of the fourth orbital center
	double centerD(int i) const
	{
		return q().centerB(i);
	}
	//! Return the distance in the \a i direction between the first two orbital centers
	double rAB(int i) const
	{
		return p().r(i);
	}
	//! Return the distance in the \a i direction between the second two orbital centers
	double rCD(int i) const
	{
		return q().r(i);
	}

	ColArrayMap P(int i) const
	{
		return _shell_quad.P(i);
	}
	DXPExpression dxP(int i, double x) const
	{
		return _shell_quad.dxP(i, x);
	}
	DPWExpression dPW(int i) const
	{
		return _shell_quad.dPW(i);
	}
	RowArrayMap Q(int i) const
	{
		return _shell_quad.Q(i);
	}
	DXQExpression dxQ(int i, double x) const
	{
		return _shell_quad.dxQ(i, x);
	}
	DQWExpression dQW(int i) const
	{
		return _shell_quad.dQW(i);
	}
	DPQExpression dPQ(int i) const
	{
		return _shell_quad.dPQ(i);
	}
	HInvWidthsABExpression hInvWidthsAB() const
	{
		return _shell_quad.hInvWidthsAB();
	}
	HInvWidthsCDExpression hInvWidthsCD() const
	{
		return _shell_quad.hInvWidthsCD();
	}
	Rho1Expression rho1() const
	{
		return _shell_quad.rho1();
	}
	Rho2Expression rho2() const
	{
		return _shell_quad.rho2();
	}

	KKWExpression KKW() const
	{
		return _shell_quad.KKW();
	}

	const Eigen::ArrayXXd& T() const
	{
		return _shell_quad.T();
	}
	const Eigen::ArrayXXd& expmT() const
	{
		return _shell_quad.expmT();
	}
	Eigen::ArrayXXd Fm(int m) const
	{
		return _shell_quad.Fm(m);
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
	template <typename Derived>
	double mulWeights(const Eigen::ArrayBase<Derived>& C) const
	{
		return ((C.colwise() * ColArray::MapAligned(p().weights().data(), p().size())).colwise().sum()
			* RowArray::MapAligned(q().weights().data(), q().size())).sum() * p().norm() * q().norm();
	}

private:
	int _ishell_quad;
	const CGTOShellQuad& _shell_quad;

	void elecRepPrim1d_aacc_psss(int i, FmCoefs& Cm) const;
	void elecRepPrim1d_abcc_psss(int i, FmCoefs& Cm) const;
	void elecRepPrim1d_aacd_psss(int i, FmCoefs& Cm) const;
	void elecRepPrim1d_abcd_psss(int i, FmCoefs& Cm) const;

	void elecRepPrim1d_abcc_ppss(int i, FmCoefs& Cm) const;
	void elecRepPrim1d_aacd_ppss(int i, FmCoefs& Cm) const;
	void elecRepPrim1d_abcd_ppss(int i, FmCoefs& Cm) const;

	void elecRepPrim1d_abcc_psps(int i, FmCoefs& Cm) const;
	void elecRepPrim1d_aacd_psps(int i, FmCoefs& Cm) const;
	void elecRepPrim1d_abcd_psps(int i, FmCoefs& Cm) const;

	void elecRepPrim1d_abcc_dsss(int i, FmCoefs& Cm) const;
	void elecRepPrim1d_aacd_dsss(int i, FmCoefs& Cm) const;
	void elecRepPrim1d_abcd_dsss(int i, FmCoefs& Cm) const;

	void elecRepPrim1d_aaaa(int i, FmCoefs& Cm) const;
	void elecRepPrim1d_aacc(int i, FmCoefs& Cm) const;
	void elecRepPrim1d_abcc(int i, FmCoefs& Cm) const;
	void elecRepPrim1d_aacd(int i, FmCoefs& Cm) const;
	void elecRepPrim1d_abcd(int i, FmCoefs& Cm) const;

	double electronRepulsion_aaaa() const;
	double electronRepulsion_aacc() const;
	double electronRepulsion_abcc() const;
	double electronRepulsion_aacd() const;
	double electronRepulsion_abcd() const;
};

#endif // CGTOQUAD_HH
