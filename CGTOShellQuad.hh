#ifndef CGTOSHELLQUAD_HH
#define CGTOSHELLQUAD_HH

/*!
 * \file CGTOShellQuad.hh
 * \brief Definition of the CGTOShellQuad class
 */

#include <memory>
#include "CGTOShellPair.hh"
#include "EriCoefs.hh"
#include "Fms.hh"
#include "AngMomPairIterator.hh"
#include "constants.hh"

class CGTOShellQuad
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
		POS_SYM_AAAA,
		//! The number of different symmetries
		POS_SYM_COUNT
	};

	typedef Eigen::ArrayXd ColArray;
	typedef Eigen::Array<double, 1, Eigen::Dynamic> RowArray;
	typedef ColArray::ConstAlignedMapType ColArrayMap;
	typedef RowArray::ConstAlignedMapType RowArrayMap;
	typedef MultiArray::ConstBlock InvWidthsSumExpression;
	typedef Eigen::CwiseBinaryOp<
		Eigen::internal::scalar_product_op<double, double>,
		const Eigen::CwiseBinaryOp<
			Eigen::internal::scalar_product_op<double, double>,
			InvWidthsSumExpression,
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
	typedef MultiArray::ConstBlock DPQExpression;
	typedef Eigen::CwiseUnaryOp<
		Eigen::internal::scalar_add_op<double>,
		ColArrayMap > DXPExpression;
	typedef Eigen::CwiseUnaryOp<
		Eigen::internal::scalar_add_op<double>,
		RowArrayMap > DXQExpression;
	typedef MultiArray::ConstBlock DPWExpression;
	typedef MultiArray::ConstBlock DQWExpression;
	typedef Eigen::CwiseBinaryOp<
		Eigen::internal::scalar_product_op<double, double>,
		const Eigen::CwiseBinaryOp<
			Eigen::internal::scalar_product_op<double, double>,
				InvWidthsSumExpression,
				const Eigen::Replicate<
					// HInvWidthsABExpression doesn't work?
					Eigen::Map<
						const ColArray,
						1,
						Eigen::Stride<0, 0> >,
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
			InvWidthsSumExpression,
			const Eigen::Replicate<
				// ColArrayMap doesn't work???
				Eigen::Map<
					const ColArray,
					1,
					Eigen::Stride<0, 0> >,
				1,
				Eigen::Dynamic> >,
		const Eigen::Replicate<
			// HInvWidthsCDExpression doesn't work?
			Eigen::Map<
				const RowArray,
				1,
				Eigen::Stride<0, 0> >,
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
						InvWidthsSumExpression>,
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

	//! Constructor
	CGTOShellQuad(const CGTOShellPair& pAB, const CGTOShellPair& pCD);

	//! Return the total number of primitive combinations in this quartet
	int size() const { return _pAB.size() * _pCD.size(); }
	//! Return the symmetry in positions of the four centers
	PositionSymmetry positionSymmetry() const { return _pos_sym; }

	//! Sums of primitive widths for the first pair
	ColArrayMap widthsAB() const
	{
		return ColArray::MapAligned(_pAB.widthsSum().data(), _pAB.size());
	}
	//! Sums of primitive widths for the second pair
	RowArrayMap widthsCD() const
	{
		return RowArray::MapAligned(_pCD.widthsSum().data(), _pCD.size());
	}
	//! Return half the inverse width sum for the first pair
	ColArrayMap hInvWidthsAB() const
	{
		return ColArray::MapAligned(_pAB.hInvWidths().data(), _pAB.size());
	}
	//! Return half the inverse width sum for the second pair
	RowArrayMap hInvWidthsCD() const
	{
		return RowArray::MapAligned(_pCD.hInvWidths().data(), _pCD.size());
	}
	//! Sums of primitive widths, for all combinations of primitives GTOs
	InvWidthsSumExpression invWidthsSum() const
	{
		return _data[0];
	}
	//! Reduced widths \f$\rho = \frac{\zeta\eta}{\zeta+\eta}\f$
	WidthsReducedExpression widthsReduced() const
	{
		return (invWidthsSum().colwise() * widthsAB()).rowwise()
			* widthsCD();
	}
	//! Return the weighted \a i coordinate of the first pair
	ColArrayMap P(int i) const
	{
		return ColArray::MapAligned(_pAB.P(i).data(), _pAB.size());
	}
	/*!
	 * \brief Return the distance from the first shell to the weighted
	 *    center in the first pair.
	 */
	DXPExpression dAP(int i) const
	{
		return P(i) - _pAB.centerA(i);
	}
	/*!
	 * \brief Return the distance in dimension \a i from the weighted center
	 *    of the first pair to that of the quartet
	 */
	DPWExpression dPW(int i) const
	{
		return _data[4+i];
	}
	//! Return the weighted \a i coordinate of the second pair
	RowArrayMap Q(int i) const
	{
		return RowArray::MapAligned(_pCD.P(i).data(), _pCD.size());
	}
	/*!
	 * \brief Return the distance from the first shell to the weighted
	 *    center in the second pair.
	 */
	DXQExpression dCQ(int i) const
	{
		return Q(i) - _pCD.centerA(i);
	}
	/*!
	 * \brief Return the distance in dimension \a i from the weighted center
	 *    of the second pair to that of the quartet
	 */
	DQWExpression dQW(int i) const
	{
		return _data[7+i];
	}
	/*!
	 * \brief Return the distance in the \a i dimension between the weighted
	 *    centers of the two pairs.
	 */
	DPQExpression dPQ(int i) const
	{
		return _data[1+i];
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
			* ColArray::MapAligned(_pAB.gaussReduced().data(), _pAB.size())).rowwise()
			* RowArray::MapAligned(_pCD.gaussReduced().data(), _pCD.size())
			* Constants::sqrt_2_pi_5_4 * Constants::sqrt_2_pi_5_4;
	}

	double getEri(int lx1, int ly1, int lz1, int lx2, int ly2, int lz2) const
	{
		if (!_have_eri)
			setEri();
		return eri(lx1, ly1, lz1, lx2, ly2, lz2);
	}
	double getEri(int lxA, int lyA, int lzA,
		int lxB, int lyB, int lzB,
		int lxC, int lyC, int lzC,
		int lxD, int lyD, int lzD) const
	{
		if (!_have_eri)
			setEri();
		return eri(lxA+lxB, lyA+lyB, lzA+lzB, lxB, lyB, lzB,
			lxC+lxD, lyC+lyD, lzC+lzD, lxD, lyD, lzD);
	}

	/*!
	 * Multiply the contributions to an integral \a C of each primitive pair
	 * with the weights of the primitives, and sum the results to compute
	 * the integral over the contraction.
	 * \param C the primitive integrals
	 */
	template <typename Derived>
	double mulWeights(const Eigen::ArrayBase<Derived>& C) const
	{
		return ((C.colwise() * ColArray::MapAligned(_pAB.weights().data(), _pAB.size())).colwise().sum()
			* RowArray::MapAligned(_pCD.weights().data(), _pCD.size())).sum();
	}
//	double mulWeights(const Eigen::ArrayXXd& C) const;

	/*!
	 * \brief Determine the positional symmetry of the shell quartet
	 *    consisting of the pairs \a pAB and \a pCD.
	 */
	static PositionSymmetry symmetry(const CGTOShellPair& pAB,
		const CGTOShellPair& pCD);

private:
	//! The first pair of shells in this quartet
	const CGTOShellPair& _pAB;
	//! The second pair of shells in this quartet
	const CGTOShellPair& _pCD;
	//! The symmetry in positions of the four orbitals
	PositionSymmetry _pos_sym;
	//! The total angular momentum of the first pair
	int _lAB;
	//! The total angular momentum of the second pair
	int _lCD;
	//! The total angular momentum of the four shells combined
	int _lsum;
	/*!
	 * \brief Data for this shell quad
	 *
	 * Data for this shell quartet. In order
	 * 0 :  Inverse sums of widths
	 * 1-3: Distance between weighted centers of the pairs
	 * 4-6: dPW
	 * 7-9: dQW
	 */
	MultiArray _data;
	//! The electron repulsion integrals themselves
	mutable Eigen::ArrayXXd _ints;
	//! Whether the electron repulsion integrals have been computed
	mutable bool _have_eri;
	
	double eri_xx(int lx1, int lx2, const EriCoefs& Cx, const Fms& fms,
		Eigen::ArrayXXd& Ctot) const;
	double eri_xx(int lx1, int ly1, int lx2, int ly2,
		const EriCoefs& Cx, const EriCoefs& Cy,
		const Fms& fms, Eigen::ArrayXXd& Cm, Eigen::ArrayXXd& Ctot) const;
	double eri_xx(const AngMomPairIterator& iter,
		const EriCoefs& Cx, const EriCoefs& Cy, const EriCoefs& Cz,
		const Fms& fms, Eigen::ArrayXXd& Cm, Eigen::ArrayXXd& Ctot) const;
	double eri_10(const EriCoefs::AllMBlock& Cxl, const Fms& fms) const;
	double eri_20(const EriCoefs::AllMBlock& Cxl, const Fms& fms) const;
	double eri_11(const EriCoefs::AllMBlock& Cxl, const EriCoefs::AllMBlock& Cyl,
		const Fms& fms, Eigen::ArrayXXd& Ctot) const;

	void elecRepPrim1d_abcd(int i, EriCoefs& coefs) const;

	void setEri() const;

	double& eri(int lx1, int ly1, int lz1, int lx2, int ly2, int lz2) const
	{
		return _ints((lx1*(_lAB+1) + ly1)*(_lAB+1) + lz1,
			(lx2*(_lCD+1) + ly2)*(_lCD+1) + lz2);
	}
	double& eri(const AngMomPairIterator& iter) const
	{
		return eri(iter.lx1(), iter.ly1(), iter.lz1(),
			iter.lx2(), iter.ly2(), iter.lz2());
	}

	double eri(int lx1, int ly1, int lz1,
		int lx2, int ly2, int lz2,
		int lzD) const;
	double eri(int lx1, int ly1, int lz1,
		int lx2, int ly2, int lz2,
		int lyD, int lzD) const;
	double eri(int lx1, int ly1, int lz1,
		int lx2, int ly2, int lz2,
		int lxD, int lyD, int lzD) const;
	double eri(int lx1, int ly1, int lz1,
		int lzB,
		int lx2, int ly2, int lz2,
		int lxD, int lyD, int lzD) const;
	double eri(int lx1, int ly1, int lz1,
		int lyB, int lzB,
		int lx2, int ly2, int lz2,
		int lxD, int lyD, int lzD) const;
	double eri(int lx1, int ly1, int lz1,
		int lxB, int lyB, int lzB,
		int lx2, int ly2, int lz2,
		int lxD, int lyD, int lzD) const;
};

#endif // CGTOSHELLQUAD_HH