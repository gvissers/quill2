#ifndef ABSTRACTBFDEF_HH
#define ABSTRACTBFDEF_HH

/*!
 * \file AbstractBFDef.hh
 * \brief Definition of the AbstractBFDef class
 */

#include <iostream>
#include "Basis.hh"

/*!
 * \brief Abstract basis function definition
 *
 * Class AbstractBF acts as a parent class for basis function definitions,
 * specifying the interface that concrete basis function definition classes
 * should conform to.
 */
struct AbstractBFDef
{
	//! Destructor
	virtual ~AbstractBFDef() {}

	/*!
	 * \brief Expand this basis function definition
	 *
	 * Expand this definition, adding basis function for position \a pos
	 * to basis \a basis.
	 * \param ipos  Position identifier
	 * \param pos   Position on which the functions are centered
	 * \param basis The basis to which the functions are added
	 */
	virtual void expand(int ipos, const Eigen::Vector3d& pos,
		Basis *basis) const = 0;

	/*!
	 * \brief Print a basis function definition
	 *
	 * Print a textual representation of this basis function definition
	 * on output stream \a os.
	 * \param os The output stream to write to
	 * \return The updated output stream
	 */
	virtual std::ostream& print(std::ostream& os) const = 0;
};

namespace
{

/*!
 * \brief Print a basis function definition
 *
 * Print a textual representation of basis function definition \a bf on
 * output stream \a os.
 * \param os The output stream to write to
 * \param bf The basis function to print
 * \return The updated output stream
 */
inline std::ostream& operator<<(std::ostream& os, const AbstractBFDef& bf)
{
	return bf.print(os);
}

} // namespace

#endif // ABSTRACTBFDEF_HH
