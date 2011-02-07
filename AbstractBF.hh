#ifndef ABSTRACTBF_HH
#define ABSTRACTBF_HH

/*!
 * \file AbstractBF.hh
 * Definition of the AbstractBF class
 */

#include <Eigen/Core>

/*!
 * \brief Base class for basis functions
 *
 * Class AbstractBF acts as a base class for basis functions, and defines the
 * interface that basis functions should adhere to. A basis function is
 * typically obtained by applying a basis function definition on an atom
 * position.
 * \tparam ResType The result type of the function, typically \c double or
 *    \c std::complex<double>.
 */
struct AbstractBF
{
	//! Destructor
	virtual ~AbstractBF() {}

	//! Evaluate this basis function in point \a pos.
	virtual double eval(const Eigen::Vector3d& pos) const = 0;

	/*!
	 * \brief Print this basis function
	 *
	 * Write a textual representation of this basis function to output
	 * stream \a os.
	 * \param os The output stream to write to
	 * \return The updated output stream
	 */
	virtual std::ostream& print(std::ostream& os) const = 0;
};

namespace {

/*!
 * \brief Print this basis function
 *
 * Write a textual representation of basis function \a bf to output stream
 * \a os.
 * \param os The output stream to write to
 * \param bf The basis function to print
 * \return The updated output stream
 */
inline std::ostream& operator<<(std::ostream& os, const AbstractBF& bf)
{
	return bf.print(os);
}

} // namespace

#endif // ABSTRACTBF_HH
