#ifndef CGTOSHELLPAIR_HH
#define CGTOSHELLPAIR_HH

/*!
 * \file CGTOShellPair.hh
 * \brief Definition of the CGTOShellPair class
 */

#include <Eigen/Core>
#include "CGTOShell.hh"

class CGTOShellPair
{
public:
	CGTOShellPair(const CGTOShell& shA, const CGTOShell& shB):
		_shA(shA), _shB(shB),
		_widthsA(shA.widths().replicate(1, shB.size())),
		_widthsB(shB.widths().transpose().replicate(shA.size(), 1)),
		_widths_sum(_widthsA + _widthsB),
		_widths_red(_widthsA * _widthsB / _widths_sum),
		_gauss_red((-(shA.center()-shB.center()).squaredNorm() * _widths_red).exp()
			/ _widths_sum),
		_hinv_widths(0.5*_widths_sum.inverse()),
		_lsum(shA.lsum() + shB.lsum())
	{
		for (int i = 0; i < 3; ++i)
		{
			_P[i] = (widthsA()*shA.center(i) + widthsB()*shB.center(i))
				/ widthsSum();
		}
	}

	//! Return the number of primitive combinations in this pair
	int size() const { return widthsSum().size(); }
	
	/*!
	 * \brief Return the widths of the primitives in the first shell,
	 *    for all primitives in the second shell.
	 */
	const Eigen::ArrayXXd& widthsA() const { return _widthsA; }
	/*!
	 * \brief Return the widths of the primitives in the second shell,
	 *    for all primitives in the first shell.
	 */
	const Eigen::ArrayXXd& widthsB() const { return _widthsB; }
	//! The sums of primitive widths, equivalent to widths() + widthB()
	const Eigen::ArrayXXd& widthsSum() const { return _widths_sum; }
	//! The "reduced" primitive widths \f$\xi = \alpha\beta / (\alpha+\beta)\f$
	const Eigen::ArrayXXd& widthsReduced() const { return _widths_red; }
	/*!
	 * \brief \f\frac{$\exp(-\xi r^2)}{\alpha+\beta}\f$ with \f$r\f$ the
	 *    distance between the centers of the two orbitals.
	 */
	const Eigen::ArrayXXd& gaussReduced() const { return _gauss_red; }
	//! Return half ithe inverse widfths summs, \f$\frac{1}{2(\alpha+\beta)}\f$
	const Eigen::ArrayXXd& hInvWidths() const { return _hinv_widths; }
	/*!
	 * \brief Return the weighted average \a i coordinate, for each
	 *    combination of primitive weights.
	 */
	const Eigen::ArrayXXd& P(int i) const
	{
		return _P[i];
	}

	//! Return the position index of the first shell
	int positionIdA() const
	{
		return _shA.positionId();
	}
	//! Return the position index of the first shell
	int positionIdB() const
	{
		return _shB.positionId();
	}
	//! Return \c true if both shells are located on the same center
	bool samePositionId() const
	{
		return positionIdA() == positionIdB();
	}

private:
	//! First shell in this pair
	const CGTOShell& _shA;
	//! Second shell in this pair
	const CGTOShell& _shB;
	//! Primitive widths of the first orbital
	Eigen::ArrayXXd _widthsA;
	//! Primitive widths of the second orbital
	Eigen::ArrayXXd _widthsB;
	//! Sum of primitive widths, for all combinations of primitives
	Eigen::ArrayXXd _widths_sum;
	//! Reduced primitive widths, for all combinations of primitives
	Eigen::ArrayXXd _widths_red;
	//! \f$\frac{\exp(-\xi r^2)}{\alpha+\beta}\f$
	Eigen::ArrayXXd _gauss_red;
	//! Half of inverse width sum
	Eigen::ArrayXXd _hinv_widths;
	//! Weighted average coordinates
	Eigen::ArrayXXd _P[3];
	//! The total angular momentum of the two shells
	int _lsum;
};

#endif // CGTOSHELLPAIR_HH