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

	typedef Eigen::ArrayXd ColArray;
	typedef Eigen::Array<double, 1, Eigen::Dynamic> RowArray;
	typedef ColArray::ConstAlignedMapType ColArrayMap;
	typedef RowArray::ConstAlignedMapType RowArrayMap;
	typedef Eigen::VectorXd::ConstAlignedMapType VectorMap;
	typedef Eigen::CwiseUnaryOp<
		Eigen::internal::scalar_add_op<double>,
		ColArrayMap > DXPExpression;
	typedef Eigen::CwiseUnaryOp<
		Eigen::internal::scalar_add_op<double>,
		RowArrayMap > DXQExpression;
	typedef Eigen::CwiseBinaryOp<
		Eigen::internal::scalar_product_op<double, double>,
		const Eigen::CwiseBinaryOp<
			Eigen::internal::scalar_product_op<double, double>,
			const Eigen::ArrayXXd,
			const Eigen::ArrayXXd>,
		const Eigen::Replicate<
			// RowArrayMap doesn't work???
			Eigen::Map<
				const RowArray,
				Eigen::Aligned,
				Eigen::Stride<0, 0> >,	
			Eigen::Dynamic,
			1> > DPWExpression;
	typedef Eigen::CwiseBinaryOp<
		Eigen::internal::scalar_product_op<double, double>,
		const Eigen::CwiseBinaryOp<
			Eigen::internal::scalar_product_op<double, double>,
			const Eigen::ArrayXXd,
			const Eigen::ArrayXXd >,
		const Eigen::Replicate<
			Eigen::CwiseUnaryOp<
				Eigen::internal::scalar_opposite_op<double>,
				ColArrayMap>,
			1,
			Eigen::Dynamic> > DQWExpression;
	typedef Eigen::CwiseUnaryOp<
		Eigen::internal::scalar_multiple_op<double>,
		const Eigen::CwiseUnaryOp<
			Eigen::internal::scalar_inverse_op<double>,
			ColArrayMap> > HInvWidthsABExpression;
	typedef Eigen::CwiseUnaryOp<
		Eigen::internal::scalar_multiple_op<double>,
		const Eigen::CwiseUnaryOp<
			Eigen::internal::scalar_inverse_op<double>,
			RowArrayMap> > HInvWidthsCDExpression;
	typedef Eigen::CwiseBinaryOp<
		Eigen::internal::scalar_sum_op<double>,
		const Eigen::Replicate<
			HInvWidthsABExpression,
			Eigen::Dynamic,
			Eigen::Dynamic>,
		const Eigen::Replicate<
			HInvWidthsCDExpression,
			Eigen::Dynamic,
			1> > DWidthsReducedExpression;
	typedef Eigen::CwiseBinaryOp<
		Eigen::internal::scalar_product_op<double, double>,
		const Eigen::CwiseBinaryOp<
			Eigen::internal::scalar_product_op<double, double>,
				const Eigen::ArrayXXd,
				const Eigen::Replicate<
					HInvWidthsABExpression,
					1,
					Eigen::Dynamic> >,
		const Eigen::Replicate<
			// RowArrayMap doesn't work ???
			Eigen::Map<
				const RowArray,
				1,
				Eigen::Stride<0, 0> >,
			Eigen::Dynamic,
			1> > Rho1Expression;
	typedef Eigen::CwiseBinaryOp<
		Eigen::internal::scalar_product_op<double, double>,
		const Eigen::CwiseBinaryOp<
			Eigen::internal::scalar_product_op<double, double>,
			const Eigen::ArrayXXd,
			const Eigen::Replicate<
				// ColArrayMap doesn't work???
				Eigen::Map<
					const ColArray,
					1,
					Eigen::Stride<0, 0> >,
				1,
				Eigen::Dynamic> >,
		const Eigen::Replicate<
			HInvWidthsCDExpression,
			Eigen::Dynamic,
			1> > Rho2Expression;
	typedef Eigen::CwiseUnaryOp<
		Eigen::internal::scalar_multiple_op<double>,
		const Eigen::CwiseUnaryOp<
			Eigen::internal::scalar_multiple_op<double>,
			const Eigen::CwiseBinaryOp<
				Eigen::internal::scalar_product_op<double, double>,
				const Eigen::CwiseBinaryOp<
					Eigen::internal::scalar_product_op<double, double>,
					const Eigen::CwiseUnaryOp<
						Eigen::internal::scalar_sqrt_op<double>,
						const Eigen::ArrayXXd >,
					const Eigen::Replicate<
						// ColArrayMap doesn't work???
						Eigen::Map<
							const ColArray,
							1,
							Eigen::Stride<0, 0> >,
						1,
						Eigen::Dynamic > >,
				const Eigen::Replicate<
					Eigen::Map<
						const RowArray,
						1,
						Eigen::Stride<0, 0> >,
					Eigen::Dynamic,
					1 > > > > KKWExpression;
	typedef Eigen::CwiseBinaryOp<
		Eigen::internal::scalar_product_op<double, double>,
		const Eigen::CwiseBinaryOp<
			Eigen::internal::scalar_product_op<double, double>,
			const Eigen::ArrayXXd,
			const Eigen::Replicate<
				Eigen::Map<
					const ColArray,
					1,
					Eigen::Stride<0, 0> >,
				1,
				Eigen::Dynamic> >,
		const Eigen::Replicate<
			Eigen::Map<
				const RowArray,
				1,
				Eigen::Stride<0, 0> >,
			Eigen::Dynamic,
			1> > WidthsReducedExpression;
	
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

	//! Sums of primitive widths for the first pair
	ColArrayMap widthsAB() const
	{
		return ColArray::MapAligned(p().widthsSum().data(), p().size());
	}
	//! Sums of primitive widths for the second pair
	RowArrayMap widthsCD() const
	{
		return RowArray::MapAligned(q().widthsSum().data(), q().size());
	}
	//! Sums of primitive widths, for all combinations of primitives GTOs
	const Eigen::ArrayXXd& invWidthsSum() const
	{
		return _inv_widths_sum;
	}
	//! Reduced widths \f$\rho = \frac{\zeta\eta}{\zeta+\eta}\f$
	WidthsReducedExpression widthsReduced() const
	{
		return (invWidthsSum().colwise() * widthsAB()).rowwise()
			* widthsCD();
	}
	
	VectorMap weightsAB() const
	{
		return Eigen::VectorXd::MapAligned(p().weights().data(), p().size());
	}
	VectorMap weightsCD() const
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

	ColArrayMap P(int i) const
	{
		return ColArray::MapAligned(p().P(i).data(), p().size());
	}
	DXPExpression dxP(int i, double x) const
	{
		return P(i) - x;
	}
	DPWExpression dPW(int i) const
	{
		return (dPQ(i) * invWidthsSum()).rowwise() * widthsCD();
	}
	RowArrayMap Q(int i) const
	{
		return RowArray::MapAligned(q().P(i).data(), q().size());
	}
	DXQExpression dxQ(int i, double x) const
	{
		return Q(i) - x;
	}
	DQWExpression dQW(int i) const
	{
		return (dPQ(i) * invWidthsSum()).colwise() * (-widthsAB());
	}
	const Eigen::ArrayXXd& dPQ(int i) const
	{
		if (!_dPQ) setDPQ();
		return _dPQ[i];
	}
	HInvWidthsABExpression hInvWidthsAB() const
	{
		return 0.5 * widthsAB().inverse();
	}
	HInvWidthsCDExpression hInvWidthsCD() const
	{
		return 0.5 * widthsCD().inverse();
	}
	DWidthsReducedExpression dWidthsReduced() const
	{
		return hInvWidthsAB().replicate(1, q().size()).rowwise()
			+ hInvWidthsCD();
	}
	Rho1Expression rho1() const
	{
		return (invWidthsSum().colwise() * hInvWidthsAB()).rowwise()
			* widthsCD();
	}
	Rho2Expression rho2() const
	{
		return (invWidthsSum().colwise() * widthsAB()).rowwise()
			* hInvWidthsCD();
	}
	
	KKWExpression KKW() const
	{
		return (invWidthsSum().sqrt().colwise()
			* ColArray::MapAligned(p().gaussReduced().data(), p().size())).rowwise()
			* RowArray::MapAligned(q().gaussReduced().data(), q().size())
			* Constants::sqrt_2_pi_5_4 * Constants::sqrt_2_pi_5_4;
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
	Eigen::ArrayXXd _inv_widths_sum;
	mutable Eigen::ArrayXXd* _dPQ;
	
	void elecRepPrim1d_aacc_psss(int i, FmCoefs& Cm) const;
	void elecRepPrim1d_abcc_psss(int i, FmCoefs& Cm) const;
	void elecRepPrim1d_abcd_psss(int i, FmCoefs& Cm) const;

	void elecRepPrim1d_abcc_ppss(int i, FmCoefs& Cm) const;
	void elecRepPrim1d_abcd_ppss(int i, FmCoefs& Cm) const;

	void elecRepPrim1d_abcc_psps(int i, FmCoefs& Cm) const;
	void elecRepPrim1d_abcd_psps(int i, FmCoefs& Cm) const;

	void elecRepPrim1d_abcc_dsss(int i, FmCoefs& Cm) const;
	void elecRepPrim1d_abcd_dsss(int i, FmCoefs& Cm) const;

	void elecRepPrim1d_aaaa(int i, FmCoefs& Cm) const;
	void elecRepPrim1d_aacc(int i, FmCoefs& Cm) const;
	void elecRepPrim1d_abcc(int i, FmCoefs& Cm) const;
	void elecRepPrim1d_abcd(int i, FmCoefs& Cm) const;

	double electronRepulsion_aaaa_ssss() const;
	double electronRepulsion_aacc_ssss() const;
	double electronRepulsion_abcd_ssss() const;
	double electronRepulsion_abcc_ssss() const;
	double electronRepulsion_aacc_psss() const;
	double electronRepulsion_abcc_psss() const;
	double electronRepulsion_abcd_psss() const;

	double electronRepulsion_aaaa() const;
	double electronRepulsion_aacc() const;
	double electronRepulsion_abcc() const;
	double electronRepulsion_abcd() const;

	void setDPQ() const;
	void freeDPQ() const { delete[] _dPQ; _dPQ = 0; }
};

#endif // CGTOQUAD_HH