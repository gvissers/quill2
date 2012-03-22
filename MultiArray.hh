#ifndef MULTIARRAY_HH
#define MULTIARRAY_HH

#include <Eigen/Core>

/*!
 * \brief Class to store multiple matrices
 *
 * Class MultiArray stores multiple two-dimenional arrays in a single, bigger
 * array, making sure that each subarray is aligned on an appropriate memory
 * boundary.
 */
class MultiArray
{
	//! The number of elements in a single (SSE) packet
	static const int pkt_size = Eigen::internal::packet_traits<double>::size;

public:
	//! Return type of a sub-array
	typedef Eigen::ArrayXXd::AlignedMapType Block;
	//! Read-only type of a sub-array
	typedef Eigen::ArrayXXd::ConstAlignedMapType ConstBlock;

	/*!
	 * \brief Constructor
	 *
	 * Create a new MultiArray object that can hold \a m matrices of size
	 * \a n1 x \a n2.
	 * \param n1 The number of rows in each sub-array
	 * \param n2 The number of columns in each sub-array
	 * \param m  The number of sub-arrays in this object
	 */
	MultiArray(int n1, int n2, int m):
		_n1(n1), _n2(n2),
		_nn(pkt_size * ((_n1*_n2+pkt_size-1) / pkt_size)),
		_A(_nn, m) {}

	//! Return the number of sub-arrays in this MultiArray
	int blocks() const { return _A.cols(); }

	//! Return the \a m'th sub-array, writable
	Block operator[](int m)
	{
		return Eigen::ArrayXXd::MapAligned(_A.data()+m*_nn, _n1, _n2);
	}
	//! Return the \a m'th sub-array, read-only
	ConstBlock operator[](int m) const
	{
		return Eigen::ArrayXXd::MapAligned(_A.data()+m*_nn, _n1, _n2);
	}

	//! Copy \a count sub-arrays from \a arr into this array
	void copy(int m, const MultiArray& arr, int am, int count)
	{
		_A.block(0, m, _nn, count) = arr._A.block(0, am, _nn, count);
	}

private:
	MultiArray(int n1, int n2, int nn, Eigen::ArrayXXd A):
		_n1(n1), _n2(n2), _nn(nn), _A(A) {}

	//! The number of rows in each sub-array
	int _n1;
	//! The number of columns in each sub-array
	int _n2;
	//! The number of elements in each sub-array, rounded up to an entire packet
	int _nn;
	//! The actual data matrix
	Eigen::ArrayXXd _A;
};

#endif // MULTIARRAY_HH