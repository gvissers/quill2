#ifndef CGTOSHELLQUAD_HH
#define CGTOSHELLQUAD_HH

/*!
 * \file CGTOShellQuad.hh
 * \brief Definition of the CGTOShellQuad class
 */

#include <memory>
#include "CGTOShellPair.hh"
#include "constants.hh"

class CGTOShellQuad
{
public:
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
	typedef const Eigen::Block<Eigen::ArrayXXd> DPQExpression;
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
			const Eigen::Block<Eigen::ArrayXXd>,
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
			const Eigen::Block<Eigen::ArrayXXd>,
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

	CGTOShellQuad(const CGTOShellPair& pAB, const CGTOShellPair& pCD):
		_pAB(pAB), _pCD(pCD),
		_inv_widths_sum((widthsAB().replicate(1, pCD.size()).rowwise()
			+ widthsCD()).inverse()) {}

	//! Return the total number of primitive combinations in this qu
	int size() const { return _inv_widths_sum.size(); }

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
		if (!_dPQ) setDPQ();
		return _dPQ->block(0, i*_pCD.size(), _pAB.size(), _pCD.size());
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

private:
	//! The first pair of shells in this quartet
	const CGTOShellPair& _pAB;
	//! The second pair of shells in this quartet
	const CGTOShellPair& _pCD;
	//! Inverse of the sums of widths, for all primitive combinations
	Eigen::ArrayXXd _inv_widths_sum;
	//! Distance between weighted centers of the two pairs
	mutable std::unique_ptr<Eigen::ArrayXXd> _dPQ;

	void setDPQ() const;
};

#endif // CGTOSHELLQUAD_HH