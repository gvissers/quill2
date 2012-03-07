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
	//! Local definition of the number of elements in a (SSE) packet
	static const int pkt_size = Eigen::internal::packet_traits<double>::size;

public:
	//! Type for the set of coefficients for a single \f$[i0,j0]^m\f$ integral
	typedef Eigen::ArrayXXd::AlignedMapType Block;
	//! Constant type for the set of coefficients for a single \f$[i0,j0]^m\f$ integral
	typedef Eigen::ArrayXXd::ConstAlignedMapType ConstBlock;
	//! Type for the set of coefficients for \f$[i0,j0]^m\f$ integrals, for all values of \f$m\f$
	class AllMBlock
	{
	public:
		//! Constructor
		AllMBlock(const EriCoefs& coefs, int i, int j):
			_coefs(coefs), _off(coefs.index(i, j, 0)) {}

		//! Retrieve the coefficients for the \f$[i0,j0]^m\f$ integral
		ConstBlock operator()(int m) const { return _coefs(_off+m); }

	private:
		//! Coefficients we refer to
		const EriCoefs& _coefs;
		//! Column offset in the coefficients matrix
		int _off;
	};

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
		_block_size(pkt_size * ((_p_size*_q_size+pkt_size-1) / pkt_size)),
		_C(_block_size, (_nrl1*_nrl2*(_nrl1+_nrl2))/2) {}

	/*!
	 * \brief Return the index where the coefficients for auxiliary integral
	 *    \f$[i0,j0]^m\f$ are located.
	 */
	int index(int i, int j, int m) const
	{
		return (i*_nrl2*(i+_nrl2) + j*(2*i+j+1)) / 2 + m;
	}

	//! Read-write access to the coefficients for auxiliary integral \f$[i0,j0]^m\f$.
	Block operator()(int i, int j, int m)
	{
		return (*this)(index(i,j,m));
	}
	//! Read-only access to the coefficients for auxiliary integral \f$[i0,j0]^m\f$.
	ConstBlock operator()(int i, int j, int m) const
	{
		return (*this)(index(i,j,m));
	}
	//! Read-write access to the coefficients at block \a idx.
	Block operator()(int idx)
	{
		return Eigen::ArrayXXd::MapAligned(_C.data()+idx*_block_size,
			_p_size, _q_size);
	}
	//! Read-only access to the coefficients at block \a idx.
	ConstBlock operator()(int idx) const
	{
		return Eigen::ArrayXXd::MapAligned(_C.data()+idx*_block_size,
			_p_size, _q_size);
	}
	//! Read-only access to coefficients for \f$[i0,j0]^m\f$, for all \f$m\f$.
	AllMBlock allM(int i, int j) const
	{
		return AllMBlock(*this, i, j);
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
	//! The size of a single coefficients block
	int _block_size;
	//! The coefficient matrix
	Eigen::ArrayXXd _C;
};

#endif // ERICOEFS_HH