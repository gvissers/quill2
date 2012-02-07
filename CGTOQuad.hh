#ifndef CGTOQUAD_HH
#define CGTOQUAD_HH

#include "AbstractBFQuad.hh"
#include "CGTOPair.hh"

class FmCoefs;

class CGTOQuad: public AbstractBFQuad
{
public:
	//! Enumeration for the symmetry in center positions of the four orbitals
	enum PositionSymmetry
	{
		//! All orbitals on different centers
		POS_SYM_ABCD,
		//! First pair is on the same center
		POS_SYM_AACD,
		//! Second pair is on the same center
		POS_SYM_ABCC,
		//! First pair is on one center, second pair on another
		POS_SYM_AACC,
		//! All orbitals are on the same center
		POS_SYM_AAAA
	};
	
	/*!
	 * \brief Constructor
	 * 
	 * Create a new CGTOQuad from pairs \a p and \a q.
	 * \param p The first orbital pair in the quartet.
	 * \param q The second orbital pair in the quartet.
	 */
	CGTOQuad(const CGTOPair& p, const CGTOPair& q);
	
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
	 * \f$\langle ab|\frac{1}{r_{12}}|cd\rangle\f$, where orbitals \f$a\f$
	 * and \f$b\f$ are stored in the first orbital pair of this quartet, and 
	 * \f$c\f$ and \f$d\f$ are stored in the second pair.
	 */
	double electronRepulsion() const;

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

	/*!
	 * \brief Sums of primitive widths for the first pair, for all
	 *     primitives in the second pair
	 */
	Eigen::ArrayXXd widthsAB() const
	{
		return Eigen::ArrayXd::MapAligned(p().widthsSum().data(), p().size())
			.replicate(1, q().size());
	}
	/*!
	 * \brief Sums of primitive widths for the second pair, for all
	 *    primitives in the first pair
	 */
	Eigen::ArrayXXd widthsCD() const
	{
		return Eigen::ArrayXd::MapAligned(q().widthsSum().data(), q().size())
			.transpose().replicate(p().size(), 1);
	}
	//! Sums of primitive widths, for all combinations of primitives GTOs
	Eigen::ArrayXXd widthsSum() const
	{
		//return widthsAB() + widthsCD();
		return widthsAB().rowwise()
			+ Eigen::ArrayXd::MapAligned(q().widthsSum().data(), q().size()).transpose();
	}
	//! Products of widths sums of the two pairs, for all combinations of primitives GTOs
	Eigen::ArrayXXd widthsProduct() const
	{
		return Eigen::VectorXd::MapAligned(p().widthsSum().data(), p().size())
			* Eigen::VectorXd::MapAligned(q().widthsSum().data(), q().size()).transpose();
	}
	//! Reduced widths \f$\rho = \frac{\zeta\eta}{\zeta+\eta}\f$
	Eigen::ArrayXXd widthsReduced() const
	{
		return widthsProduct() / widthsSum();
	}
	
	Eigen::VectorXd weightsAB() const
	{
		return Eigen::VectorXd::MapAligned(p().weights().data(), p().size());
	}
	Eigen::VectorXd weightsCD() const
	{
		return Eigen::VectorXd::MapAligned(q().weights().data(), q().size());
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

	Eigen::ArrayXXd P(int i) const
	{
		return Eigen::ArrayXd::MapAligned(p().P(i).data(), p().size())
			.replicate(1, q().size());
	}
	Eigen::ArrayXXd Q(int i) const
	{
		return Eigen::ArrayXd::MapAligned(q().P(i).data(), q().size())
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
		return Eigen::VectorXd::MapAligned(p().K().data(), p().size())
			* Eigen::VectorXd::MapAligned(q().K().data(), q().size()).transpose();
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

private:
	PositionSymmetry _pos_sym;
	
	void elecRepPrim1d_aacc_psss(int i, double Pi, double Qi,
		FmCoefs& Cm) const;
	void elecRepPrim1d_abcd_psss(int i,
		const Eigen::ArrayXXd& Pi, const Eigen::ArrayXXd& Qi,
		FmCoefs& Cm) const;

	void elecRepPrim1d_abcd_ppss(int i,
		const Eigen::ArrayXXd& Pi, const Eigen::ArrayXXd& Qi,
		FmCoefs& Cm) const;
	void elecRepPrim1d_abcd_psps(int i,
		const Eigen::ArrayXXd& Pi, const Eigen::ArrayXXd& Qi,
		FmCoefs& Cm) const;
	void elecRepPrim1d_abcd_dsss(int i,
		const Eigen::ArrayXXd& Pi, const Eigen::ArrayXXd& Qi,
		FmCoefs& Cm) const;

	void elecRepPrim1d_aaaa(int i, FmCoefs& Cm) const;
	void elecRepPrim1d_aacc(int i, double Pi, double Qi, FmCoefs& Cm) const;
	void elecRepPrim1d_abcd(int i,
		const Eigen::ArrayXXd& Pi, const Eigen::ArrayXXd& Qi,
		FmCoefs& Cm) const;

	double electronRepulsion_aaaa() const;
	double electronRepulsion_aacc() const;
	double electronRepulsion_abcd() const;

	double electronRepulsion_aaaa_ssss() const;
	double electronRepulsion_aacc_ssss() const;
	double electronRepulsion_abcd_ssss() const;
	double electronRepulsion_aacc_psss() const;
	double electronRepulsion_abcd_psss() const;
};

#endif // CGTOQUAD_HH