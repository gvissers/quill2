#ifndef ABSTRACTBF_HH
#define ABSTRACTBF_HH

/*!
 * \file AbstractBF.hh
 * \brief Definition of the AbstractBF class
 */

#include <iostream>

/*!
 * \brief Abstract basis function
 *
 * Class AbstractBF acts as a parent class for basis function, specifying the
 * the interface that concrete basis function classes should conform to.
 */
struct AbstractBF
{
	//! Destructor
	virtual ~AbstractBF() {}
	/*!
	 * \brief Print a basis function
	 *
	 * Print a textual representation of this basis function on output
	 * stream \a os.
	 * \param os The output stream to write to
	 * \return The updated output stream
	 */
	virtual std::ostream& print(std::ostream& os) const = 0;
};

namespace
{

/*!
 * \brief Print a basis function
 *
 * Print a textual representation of basis function \a bf on output stream
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
