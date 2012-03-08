#ifndef CGTOSHELLQUAD_HH
#define CGTOSHELLQUAD_HH

/*!
 * \file CGTOShellQuad.hh
 * \brief Definition of the CGTOShellQuad class
 */

#include <memory>
#include "CGTOShellPair.hh"
#include "EriCoefs.hh"
#include "constants.hh"

// Forward declaration
class Fms;

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
	typedef const Eigen::Block<const Eigen::ArrayXXd> DPQExpression;
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
			DPQExpression,
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
			DPQExpression,
			const Eigen::ArrayXXd >,
		const Eigen::Replicate<
			Eigen::CwiseUnaryOp<
				Eigen::internal::scalar_opposite_op<double>,
				ColArrayMap>,
			1,
			Eigen::Dynamic> > DQWExpression;
		typedef Eigen::CwiseBinaryOp<
		Eigen::internal::scalar_product_op<double, double>,
		const Eigen::CwiseBinaryOp<
			Eigen::internal::scalar_product_op<double, double>,
				const Eigen::ArrayXXd,
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

	CGTOShellQuad(const CGTOShellPair& pAB, const CGTOShellPair& pCD);

	//! Return the total number of primitive combinations in this qu
	int size() const { return _inv_widths_sum.size(); }
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
	//! Return the weighted \a i coordinate of the first pair
	ColArrayMap P(int i) const
	{
		return ColArray::MapAligned(_pAB.P(i).data(), _pAB.size());
	}
	//! Return the distance from \a x to the first weighted center
	DXPExpression dxP(int i, double x) const
	{
		return P(i) - x;
	}
	/*!
	 * \brief Return the distance in dimension \a i from the weighted center
	 *    of the first pair to that of the quartet
	 */
	DPWExpression dPW(int i) const
	{
		return (dPQ(i) * invWidthsSum()).rowwise() * widthsCD();
	}
	//! Return the weighted \a i coordinate of the second pair
	RowArrayMap Q(int i) const
	{
		return RowArray::MapAligned(_pCD.P(i).data(), _pCD.size());
	}
	//! Return the distance from \a x to the second weighted center
	DXQExpression dxQ(int i, double x) const
	{
		return Q(i) - x;
	}
	/*!
	 * \brief Return the distance in dimension \a i from the weighted center
	 *    of the second pair to that of the quartet
	 */
	DQWExpression dQW(int i) const
	{
		return (dPQ(i) * invWidthsSum()).colwise() * (-widthsAB());
	}
	/*!
	 * \brief Return the distance in the \a i dimension between the weighted
	 *    centers of the two pairs.
	 */
	DPQExpression dPQ(int i) const
	{
		return _dPQ.block(0, i*_pCD.size(), _pAB.size(), _pCD.size());
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

	const Eigen::ArrayXXd& T() const { return _T; }
	const Eigen::ArrayXXd& expmT() const { return _expmT; }
	Eigen::ArrayXXd Fm(int m) const;
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
		return eri(lxA, lyA, lzA, lxB, lyB, lzB, lxC, lyC, lzC, lxD, lyD, lzD);
	}

	/*!
	 * Multiply the contributions to an integral \a C of each primitive pair
	 * with the weights of the primitives, and sum the results to compute
	 * the integral over the contraction.
	 * \param C the primitve integrals
	 */
	template <typename Derived>
	double mulWeights(const Eigen::ArrayBase<Derived>& C) const
	{
		return ((C.colwise() * ColArray::MapAligned(_pAB.weights().data(), _pAB.size())).colwise().sum()
			* RowArray::MapAligned(_pCD.weights().data(), _pCD.size())).sum();
	}

private:
	//! The first pair of shells in this quartet
	const CGTOShellPair& _pAB;
	//! The second pair of shells in this quartet
	const CGTOShellPair& _pCD;
	//! Inverse of the sums of widths, for all primitive combinations
	Eigen::ArrayXXd _inv_widths_sum;
	//! Distance between weighted centers of the two pairs
	Eigen::ArrayXXd _dPQ;
	//! The symmetry in positions of the four orbitals
	PositionSymmetry _pos_sym;
	//! The total angular momentum of the four shells combined
	int _lsum;

	Eigen::ArrayXXd _T;
	Eigen::ArrayXXd _expmT;
	mutable int _m;
	mutable Eigen::ArrayXXd _Fm;
	mutable Eigen::ArrayXXd _ints;
	//! Whether the electron repulsion integrals have been computed
	mutable bool _have_eri;

	void elecRepPrim1d_abcd(int i, EriCoefs& coefs) const;
	
	double eri_xx(int lx1, int ly1, int lz1, int lx2, int ly2, int lz2,
		const EriCoefs& Cx, const EriCoefs& Cy, const EriCoefs& Cz,
		const Fms& fms, Eigen::ArrayXXd& Cm, Eigen::ArrayXXd& Ctot) const;
	double eri_10(const EriCoefs::AllMBlock& Cxl, const Fms& fms) const;
	double eri_20(const EriCoefs::AllMBlock& Cxl, const Fms& fms) const;
	double eri_11(const EriCoefs::AllMBlock& Cxl,
		const EriCoefs::AllMBlock& Cyl, const Fms& fms) const;
	void setEri() const;

	double& eri(int lx1, int ly1, int lz1, int lx2, int ly2, int lz2) const
	{
		return _ints(lx1*(_lsum+1)*(_lsum+1)+ly1*(_lsum+1)+lz1,
			lx2*(_lsum+1)*(_lsum+1)+ly2*(_lsum+1)+lz2);
	}
	double eri(int lx1, int ly1, int lz1,
		int lx2, int ly2, int lzC,
		int lzD) const;
	double eri(int lx1, int ly1, int lz1,
		int lx2, int lyC, int lzC,
		int lyD, int lzD) const;
	double eri(int lx1, int ly1, int lz1,
		int lxC, int lyC, int lzC,
		int lxD, int lyD, int lzD) const;
	double eri(int lx1, int ly1, int lzA,
		int lzB,
		int lxC, int lyC, int lzC,
		int lxD, int lyD, int lzD) const;
	double eri(int lx1, int lyA, int lzA,
		int lyB, int lzB,
		int lxC, int lyC, int lzC,
		int lxD, int lyD, int lzD) const;
	double eri(int lxA, int lyA, int lzA,
		int lxB, int lyB, int lzB,
		int lxC, int lyC, int lzC,
		int lxD, int lyD, int lzD) const;
};

#endif // CGTOSHELLQUAD_HH