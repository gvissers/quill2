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
	/*!
	 * \brief Constructor
	 *
	 * Create a new pair of shells \a shA and \a shB
	 */ 
	CGTOShellPair(int index, const CGTOShell& shA, const CGTOShell& shB);

	//! Return the number of primitive combinations in this pair
	int size() const { return widthsSum().size(); }

	//! Return the index of this shell pair in the CGTOShellList
	int index() const { return _index; }
	//! Return the total angular momentum of the swo shells combined
	int lsum() const { return _lsum; }
	//! Return the unnormalized weights of the primitive combinations
	const Eigen::ArrayXXd& weights() const { return _weights; }
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
	//! Return half the inverse widths sums, \f$\frac{1}{2(\alpha+\beta)}\f$
	const Eigen::ArrayXXd& hInvWidths() const { return _hinv_widths; }

	//! Return the \a i coordinate of the first shell in this pair
	double centerA(int i) const { return _shA.center(i); }
	//! Return the \a i coordinate of the second shell in this pair
	double centerB(int i) const { return _shB.center(i); }
	//! Return the distance between the two centers in the \a i dimension
	double dAB(int i) const { return _shB.center(i) - _shA.center(i); }
	
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

	/*!
	 * Multiply the contributions to an integral \a C of each primitive pair
	 * with the weights of the primitives, and sum the results to compute
	 * the integral over the contraction.
	 * \param C the primitve integrals
	 */
	template <typename Derived>
	double mulWeights(const Eigen::ArrayBase<Derived>& C) const
	{
		return ((C.colwise() * _shA.weights()).colwise().sum().transpose()
			* _shB.weights()).sum();
	}

private:
	//! Index of this pair in the CGTOShellList
	int _index;
	//! First shell in this pair
	const CGTOShell& _shA;
	//! Second shell in this pair
	const CGTOShell& _shB;
	//! Unnormalized weights of the primitive combinations
	Eigen::ArrayXXd _weights;
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