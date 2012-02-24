#ifndef ERICOEFS_HH
#define ERICOEFS_HH

/*!
 * \file EriCoefs.hh
 * \brief Definition of the EriCoefs class
 */

#include <Eigen/Core>

/*!
 * \brief Class for storing ERI coefficients
 *
 * Class EriCoefs stores coefficients used in the calculation of electron
 * repulsion integrals using contracted Gaussian type orbitals in a single
 * matrix, and provides access methods to set and retrieve those values.
 */
class EriCoefs
{
public:
	/*!
	 * \brief Constructor
	 *
	 * Create a new EriCoefs object that can hold all coefficients needed
	 * for computing an ERI with total angular momentum \a l1 on the first
	 * shell pair, and \a l2 on the second pair. Parameters \a p_size and
	 * \a q_size are the total number of primitive gaussian combinations in
	 * the two orbitals pairs.
	 * \param l1     The total angular momentum of the first shell pair
	 * \param l2     The total angular momentum of the second shell pair
	 * \param p_size The number of primitive combinations in the first pair
	 * \param p_size The number of primitive combinations in the second pair
	 */
	EriCoefs(int l1, int l2, int p_size, int q_size):
		_nrl1(l1+1), _nrl2(l2+1), _p_size(p_size), _q_size(q_size),
		_C(p_size, (_nrl1*_nrl2*(_nrl1+_nrl2)*q_size)/2) {}

	/*!
	 * \brief Return the index where the coefficients for auxiliary integral
	 *    \f$[i0,j0]^m\f$ are located.
	 */
	int index(int i, int j, int m) const
	{
		return (i*_nrl2*(i+_nrl2) + j*(2*i+j+1)) / 2 + m;
	}

	//! Read-write access to the coefficients for auxiliary integral \f$[i0,j0]^m\f$.
	Eigen::Block<Eigen::ArrayXXd> operator()(int i, int j, int m)
	{
		return (*this)(index(i,j,m));
	}
	//! Read-only access to the coefficients for auxiliary integral \f$[i0,j0]^m\f$.
	const Eigen::Block<const Eigen::ArrayXXd> operator()(int i, int j,
		int m) const
	{
		return (*this)(index(i,j,m));
	}
	//! Read-write access to the coefficients at block \a idx.
	Eigen::Block<Eigen::ArrayXXd> operator()(int idx)
	{
		return _C.block(0, idx*_q_size, _p_size, _q_size);
	}
	//! Read-only access to the coefficients at block \a idx.
	const Eigen::Block<const Eigen::ArrayXXd> operator()(int idx) const
	{
		return _C.block(0, idx*_q_size, _p_size, _q_size);
	}

	/*!
	 * \brief Read-write access to the coefficients for auxiliary integrals
	 *    \f$[i0,j0]^m\f$ up to \f$[i0,j0]^{m+n}\f$.
	 */
	Eigen::Block<Eigen::ArrayXXd> operator()(int i, int j, int m, int n)
	{
		return _C.block(0, (index(i,j,m))*_q_size, _p_size,
			n*_q_size);
	}
	/*!
	 * \brief Read-only access to the coefficients for auxiliary integrals
	 *    \f$[i0,j0]^m\f$ up to \f$[i0,j0]^{m+n}\f$.
	 */
	const Eigen::Block<const Eigen::ArrayXXd> operator()(int i, int j,
		int m, int count) const
	{
		return _C.block(0, (index(i,j,m))*_q_size, _p_size,
			count*_q_size);
	}

private:
	//! The number of possible different values for i
	int _nrl1;
	//! The number of possible different values for j
	int _nrl2;
	//! The number of primitive combinations in the first shell pair
	int _p_size;
	//! The number of primitive combinations in the second shell pair
	int _q_size;
	//! The coefficient matrix
	Eigen::ArrayXXd _C;
};

#endif // ERICOEFS_HH