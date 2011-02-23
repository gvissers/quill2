#ifndef ABSTRACTBF_HH
#define ABSTRACTBF_HH

/*!
 * \file AbstractBF.hh
 * \brief Definition of the AbstractBF class
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
	//! Unique class ID, used in looking up integral calculation functions
	size_t cid;

	//! Constructor
	AbstractBF(size_t cid): cid(cid) {}
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

/*!
 * \brief The overlap between basis functions
 *
 * Compute the overlap between basis functions \a f and \a g. This looks up
 * the overlap function to use for the concrete types of \a f and \a g in
 * the Dispatcher, and executes the function found.
 * \param f The first basis function
 * \param g The second basis function
 * \return The overlap between \a f and \a g
 */
double overlap(const AbstractBF& f, const AbstractBF& g);
/*!
 * \brief The kinetic energy integral between basis functions
 *
 * Compute the kinetic energy matrix element over basis functions \a f and
 * \a g. This looks up the function to use for the concrete types of \a f
 * and \a g in the Dispatcher, and executes the function found.
 * \param f The first basis function
 * \param g The second basis function
 * \return The kinetic energy integral between \a f and \a g
 */
double kineticEnergy(const AbstractBF& f, const AbstractBF& g);
/*!
 * \brief The one electron integrals between basis functions
 *
 * Compute the overlap and kinetic energy matrix elements over basis functions
 * \a f and \a g. This looks up the function to use for the concrete types of
 * \a f and \a g in the Dispatcher, and executes the function found.
 * \param f The first basis function
 * \param g The second basis function
 * \param S Place to store the overlap between \a f and \a g
 * \param T Place to store the kinetic energy integral between \a f and \a g
 */
void oneElectron(const AbstractBF& f, const AbstractBF& g, double *S, double *T);

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
